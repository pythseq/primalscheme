import logging
import primer3
from primal import settings
from .exceptions import MaxGapReached, NoSuitablePrimers
from .models import _primer, _primerPair, _candidatePrimer, _candidatePrimerPair, _region
from Bio.Align import MultipleSeqAlignment
import sys
from itertools import product
from pprint import pprint

logger = logging.getLogger('Primal Log')

class poaMultiplexScheme(object):
    """A complete multiplex primer scheme."""

    def __init__(self, references, amplicon_length, min_overlap, max_gap, max_alts, max_candidates,
                 step_size, max_variation, prefix='PRIMAL_SCHEME'):
        self.references = references
        self.amplicon_length = amplicon_length
        self.min_overlap = min_overlap
        self.max_gap = max_gap
        self.max_alts = max_alts
        self.max_candidates = max_candidates
        self.step_size = step_size
        self.max_variation = max_variation
        self.prefix = prefix
        self.regions = []

        self.run()

    @property
    def primaryReference(self):
        return self.references[0]

    def run(self):
        regions = []
        region_num = 0
        is_last_region = False

        while True:
            region_num += 1
            # Get the previous region in each pool
            prev_pair = regions[-1].candidatePairs[0] if len(regions) >= 1 else None
            prev_pair_same_pool = regions[-2].candidatePairs[0] if len(regions) > 2 else None

            # If there are two regions or more
            if prev_pair_same_pool:
                # Gap opened between -1 and -2 regions
                if prev_pair.left.startEnd('fwd')[0] > prev_pair_same_pool.right.startEnd('rev')[0]:
                    # If there was a gap left primer cannot overlap with -1 region
                    left_primer_left_limit = prev_pair.left.startEnd('fwd')[1] + 1
                else:
                    # Left primer cannot overlap -2 region
                    left_primer_left_limit = prev_pair_same_pool.right.startEnd('rev')[0] + 1
            # If there is more than one region
            elif prev_pair:
                # Left primer cannot overlap with -1 region or you don't move
                left_primer_left_limit = prev_pair.left.startEnd('fwd')[1] + 1
            else:
                # Region one only limit is 0
                # TODO: start when 90% of the genomes have sequence
                left_primer_left_limit = 0

            # Right start limit maintains the minimum_overlap
            left_primer_right_limit = prev_pair.right.startEnd('rev')[1] - self.min_overlap - 1 if prev_pair else self.max_gap

            # Last region if less than one amplicon length remaining
            if prev_pair:
                if (len(self.primaryReference) - prev_pair.right.startEnd('rev')[1]) < self.amplicon_length:
                    is_last_region = True
                    logger.debug('Region {}: is last region'.format(region_num))

            # Log limits
            logger.debug('Region {}: forward primer limits {}:{}'.format(region_num, left_primer_left_limit, left_primer_right_limit))

            # Find primers or handle no suitable error
            try:
                region = self.poaFindPrimers(region_num, left_primer_left_limit, left_primer_right_limit, is_last_region)
                regions.append(region)
            except NoSuitablePrimers:
                logger.debug('Region {}: no suitable primer error'.format(region_num))
                break

            # Handle the end
            if is_last_region:
                logger.debug('Region {}: ending normally'.format(region_num))
                break

            """
            # Report scores and alignments
            for i in range(0, len(self.references)):
                # Don't display alignment to reference
                logger.debug(regions[-1].candidate_pairs[0].left.alignments[i].formatted_alignment)
            logger.debug('Identities for sorted left candidates: ' + ','.join(['%.2f' %each.left.percent_identity for each in regions[-1].candidate_pairs]))
            logger.debug('Left start for sorted candidates: ' + ','.join(['%i' %each.left.start for each in regions[-1].candidate_pairs]))
            logger.debug('Left end for sorted candidates: ' + ','.join(['%i' %each.left.end for each in regions[-1].candidate_pairs]))
            logger.debug('Left length for sorted candidates: ' + ','.join(['%i' %each.left.length for each in regions[-1].candidate_pairs]))

            for i in range(0, len(self.references)):
                logger.debug(regions[-1].candidate_pairs[0].right.alignments[i].formatted_alignment)
            logger.debug('Identities for sorted right candidates: ' + ','.join(['%.2f' %each.right.percent_identity for each in regions[-1].candidate_pairs]))
            logger.debug('Right start for sorted candidates: ' + ','.join(['%i' %each.right.start for each in regions[-1].candidate_pairs]))
            logger.debug('Right end for sorted candidates: ' + ','.join(['%i' %each.right.end for each in regions[-1].candidate_pairs]))
            logger.debug('Right length for sorted candidates: ' + ','.join(['%i' %each.right.length for each in regions[-1].candidate_pairs]))

            logger.debug('Totals for sorted pairs: ' + ','.join(['%.2f' %each.mean_percent_identity for each in regions[-1].candidate_pairs]))
            """
            if len(regions) > 1:
            # Remember, results now include this one, so -2 is the other pool
                trimmed_overlap = regions[-2].candidatePairs[0].right.startEnd('rev')[1] - regions[-1].candidatePairs[0].left.startEnd('fwd')[1] - 1
                logger.info("Region %i: highest scoring product %i:%i, length %i, trimmed overlap %i" % (region_num, regions[-1].candidatePairs[0].left.startEnd('fwd')[0], regions[-1].candidatePairs[0].right.startEnd('rev')[0], regions[-1].candidatePairs[0].productLength, trimmed_overlap))
            else:
                logger.info("Region %i: highest scoring product %i:%i, length %i" % (region_num, regions[-1].candidatePairs[0].left.startEnd('fwd')[0], regions[-1].candidatePairs[0].right.startEnd('rev')[0], regions[-1].candidatePairs[0].productLength))


        # Return regions
        self.regions = regions

    def digestSeq(self, k, start, seq):
            return [(int(start)+i, seq[i:i+k]) for i in range((len(seq)-k)+1)]

    def poaFindPrimers(self, region_num, left_primer_left_limit, left_primer_right_limit, is_last_region):
        """
        Find primers for a given region.

        Return a list of Region objects containing candidate
        primer pairs sorted by mean percent identity against all references.
        """

        # Calculate where to slice the reference
        if region_num == 1:
            chunk_start = 0
            chunk_end = int((1 + self.max_variation / 2) * self.amplicon_length)
        elif is_last_region:
            # Last region work backwards
            chunk_start = int(len(self.primaryReference) - ((1 + self.max_variation / 2) * self.amplicon_length))
            chunk_end = len(self.primaryReference)
        else:
            # Start is right limit - delta min/max product length - max primer length
            chunk_start = int(left_primer_right_limit - (self.max_variation * self.amplicon_length) - settings.global_args['PRIMER_MAX_SIZE'])
            chunk_end = int(chunk_start + ((1 + self.max_variation/2) * self.amplicon_length))
        initial_chunk_start = chunk_start
        initial_chunk_end = chunk_end

        # Primer3 setup
        p3_global_args = settings.global_args
        p3_seq_args = settings.seq_args
        p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [
            [int(self.amplicon_length * (1 - self.max_variation / 2)), int(self.amplicon_length * (1 + self.max_variation / 2))]]
        p3_global_args['PRIMER_NUM_RETURN'] = self.max_candidates

        # Run digestSeq until unique primers are found
        hit_left_limit = False
        while True:
            # Slice primary reference
            #seq = str(self.primaryReference.seq[chunk_start:chunk_end])
            #p3_seq_args['SEQUENCE_TEMPLATE'] = seq
            #p3_seq_args['SEQUENCE_INCLUDED_REGION'] = [0, len(seq) - 1]
            #p3_seq_args['SEQUENCE_PRIMER'] = 'TCTTTTGTGTGCGAATAACTATGAGGA'
            #print(p3_seq_args['SEQUENCE_TEMPLATE'], p3_seq_args['SEQUENCE_INCLUDED_REGION'])
            logger.info("Region %i: reference chunk %i:%i, length %i" %(region_num, chunk_start, chunk_end, chunk_end-chunk_start))
            #primer3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
            #pprint(primer3_output)

            #Digest the references (except the consensus) into candidatePrimers
            allKmers = set()
            for ref in self.references[1:]:
                seq = str(ref.seq[chunk_start:chunk_end])
                for k in range(22, 30+1):
                    allKmers.update(self.digestSeq(k, chunk_start, seq))

            #Filter for valid start and end position
            #fwdPos = [p for p in allKmers if p.startEnd('fwd')[0] >= chunk_start and p.startEnd('fwd')[1] <= chunk_start+40]
            #revPos = [p for p in allKmers if p.startEnd('rev')[0] <= chunk_end and p.startEnd('rev')[1] >= chunk_end-40]
            fwdPos = [_candidatePrimer(k[0], k[1], 'fwd') for k in allKmers if chunk_start <= _candidatePrimer(k[0], k[1], 'fwd').start < chunk_start+40]
            revPos = [_candidatePrimer(k[0], k[1], 'rev') for k in allKmers if chunk_end-40 <= _candidatePrimer(k[0], k[1], 'rev').start < chunk_end]
            #print(vars([p for p in fwdPos if p.seq == 'AGAATTTTTAGGATCTTTTGTGTGCGA'][0]))
            #print(vars([p for p in revPos if p.seq == 'GGGTTGTTGAATCTCCAATCCTCT'][0]))

            #Perform the hard filtering
            fwdThermo = [p for p in fwdPos if (30 <= p.gc <= 55) and (60 <= p.tm <= 63) and (p.hairpin <= 50.0) and (p.maxPoly <= 5)]
            revThermo = [p for p in revPos if (30 <= p.gc <= 55) and (60 <= p.tm <= 63) and (p.hairpin <= 50.0) and (p.maxPoly <= 5)]

            #Filter for valid length
            pairs = [_primerPair(f,r) for f in fwdThermo for r in revThermo if 380 <= _primerPair(f,r).productLength <= 420]
            logger.info("Region %i: current position returned %i left and %i right candidate primers" %(region_num, len(fwdThermo), len(revThermo)))
            logger.info("Region %i: current position returned %i candidate primer pairs" %(region_num, len(pairs)))

            if pairs:
                #Sort pairs on left ref coverage, right ref coverage
                #Sort pairs on pairPenalty
                scoredPairs = [_candidatePrimerPair(p.left, p.right) for p in pairs]
                sortPairs = sorted(scoredPairs, key=lambda x: x.pairPenalty)
                #sortPairs = [p for p in scoredPairs if p.left.seq == 'AGAATTTTTAGGATCTTTTGTGTGCGA' and p.right.seq == 'GGGTTGTTGAATCTCCAATCCTCT']
                #sortPairs = sorted(scoredPairs, key=lambda x: (len(x.left.queryMatch(self.references[1:])), len(x.right.queryMatch(self.references[1:]))), reverse=True)
                #print('top pair', len(sortPairs[0].left.queryMatch(self.references[1:])), len(sortPairs[0].right.queryMatch(self.references[1:])))
                print(sorted([p.pairPenalty for p in sortPairs][:20]))

                pprint(vars(sortPairs[0]))
                print(sortPairs[0].productLength)
                pprint(vars(sortPairs[0].left))
                pprint(vars(sortPairs[0].right))

                print('gc', sortPairs[0].left.gc, sortPairs[0].right.gc)
                print('hairpin', sortPairs[0].left.hairpin, sortPairs[0].right.hairpin)
                print('homodimer', sortPairs[0].left.homodimer, sortPairs[0].right.homodimer)
                print('tm', sortPairs[0].left.tm, sortPairs[0].right.tm)

                #Get list of alts or None if failed to cover all references
                leftAlts = sortPairs[0].fwdAlts(self.references, pairs, sortPairs)
                rightAlts = sortPairs[0].revAlts(self.references, pairs, sortPairs)

                #If there is a list of alts return (even if empty)
                if not leftAlts == None:
                    print('left alts', len(leftAlts), [len(p.queryMatch(self.references[1:])) for p in leftAlts])
                    if not rightAlts == None:
                        print('right alts', len(rightAlts), [len(p.queryMatch(self.references[1:])) for p in rightAlts])
                        return _region(region_num, chunk_start, sortPairs, leftAlts, rightAlts, self.references, self.prefix, self.max_alts)

            # Move right if first region or to open gap
            if region_num == 1 or hit_left_limit:
                logger.debug("Region %i: stepping right, position %i" %(region_num, chunk_start))
                chunk_start += self.step_size
                chunk_end += self.step_size
                # Hit end of regerence
                if chunk_end > len(self.primaryReference):
                    logger.debug("Region %i: hit right limit %i" %(region_num, len(self.primaryReference)))
                    raise NoSuitablePrimers("No suitable primers in region")
            else:
                # Move left for all other regions
                logger.debug("Region %i: stepping left, position %i, limit %s" %(region_num, chunk_start, left_primer_left_limit))
                chunk_start -= self.step_size
                chunk_end -= self.step_size
                if chunk_start <= left_primer_left_limit:
                    # Switch direction to open gap
                    logger.debug("Region %i: hit left limit" %(region_num))
                    chunk_start = initial_chunk_start
                    chunk_end = initial_chunk_end
                    hit_left_limit = True
