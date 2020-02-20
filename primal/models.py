import os
import sys
import logging
from primal import settings
from itertools import groupby

from Bio.Seq import Seq
from primer3 import calcTm, calcHairpin, calcHomodimer
from Porechop.porechop.cpp_function_wrappers import adapter_alignment

logger = logging.getLogger('Primal Log')

class _primer(object):
    """A simple primer."""

    def __init__(self, position, seq):
        self.position = position
        self.seq = seq
        self.name = None
        #self.penalty = None

    def startEnd(self, direction):
        if direction == 'fwd':
            return (self.position, self.position + self.length)
        if direction == 'rev':
            return (self.position + self.length, self.position)

    #Stability of homodimer
    def homodimer(self, direction):
        if direction == 'fwd':
            return calcHomodimer(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm
        if direction == 'rev':
            return calcHomodimer(Seq(self.seq).reverse_complement()._data, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

    #Stability of hairpin
    def hairpin(self, direction):
        if direction == 'fwd':
            return calcHairpin(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm
        if direction == 'rev':
            return calcHairpin(Seq(self.seq).reverse_complement()._data, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

    #Check the thermo calculations are the same for fwd and rev
    def revComp(self):
        return Seq(self.seq).reverse_complement()

    @property
    def tm(self):
        return calcTm(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6)

    #Stability of last 5 3' bases
    @property
    def endStability(self):
        #Only works for fwd primers
        return calcTm(self.seq[-5:], mv_conc=50, dv_conc=1.5, dntp_conc=0.6)

    #GC content
    @property
    def gc(self):
        return 100.0 * (self.seq.count('G') + self.seq.count('C')) / len(self.seq)

    #Max homopolymer length using itertools
    @property
    def maxPoly(self):
        return sorted([(len(list(g))) for k,g in groupby(self.seq)], reverse=True)[0]

    #Length of primer
    @property
    def length(self):
        return len(self.seq)


class _candidatePrimer(_primer):
    """A candidate primer for a region."""

    def __init__(self, position, seq):
        super(_candidatePrimer, self).__init__(position, seq)
        #primerCov = None

    def __eq__(self, other):
        return self.seq == other.seq

    def __hash__(self):
        return hash(self.seq)

    def queryAlign(self, references):
        return [ref.id for ref in references if ref[self.startEnd('fwd')[0]:self.startEnd('fwd')[1]].seq == self.seq]

class _primerPair(object):
    """A pair of primers for a region."""

    def __init__(self, left, right):
        self.left = left
        self.right = right
        #self.penalty = None

class _candidatePrimerPair(object):
    """A pair of candidate primers for a region."""

    def __init__(self, left, right):
        self.left = left
        self.right = right
        #self.leftAlts = []
        #self.rightAlts = []

    def fwdAlts(self, references, pairs, sortPairs, max_alts=5):
        #Update set of refs covered
        fwdCov = set(self.left.queryAlign(references[1:]))
        leftAlts = []
        #breakOut = False
        #Generate left alts
        while len(fwdCov) < len(references[1:]):
            #Refs not covered
            fwdReq = [r for r in references[1:] if r.id not in fwdCov]
            #Alts with valid length product
            fwdAlts = [p.left for p in pairs if p.right == sortPairs[0].right]
            #Sort alts on required refs then all refs
            sortFwdAlts = sorted(fwdAlts, key=lambda x: (len(x.queryAlign(fwdReq)), len(x.queryAlign(references[1:]))), reverse=True)
            #Store alts
            leftAlts.append(sortFwdAlts[0])
            #Update refs covered
            lastCov = len(fwdCov)
            fwdCov.update(sortFwdAlts[0].queryAlign(references[1:]))
            if lastCov == len(fwdCov):
                return None
        return leftAlts

    def revAlts(self, references, pairs, sortPairs):
        revCov = set(sortPairs[0].right.queryAlign(references[1:]))
        rightAlts = []
        #breakOut = False
        while len(revCov) < len(references[1:]):
            revReq = [r for r in references[1:] if r.id not in revCov]
            revAlts = [p.right for p in pairs if p.left == sortPairs[0].left]
            sortRevAlts = sorted(revAlts, key=lambda x: (len(x.queryAlign(revReq)), len(x.queryAlign(references[1:]))), reverse=True)
            rightAlts.append(sortRevAlts[0])
            lastCov = len(revCov)
            revCov.update(sortRevAlts[0].queryAlign(references[1:]))
            if lastCov == len(revCov):
                return None
        return rightAlts

    @property
    def productLength(self):
        return self.right.startEnd('rev')[0] - self.left.startEnd('fwd')[0] + 1


class _region(object):
    """A region that forms part of a scheme."""
    def __init__(self, region_num, chunk_start, candidatePairs, fwdAlternates, revAlternates, references, prefix, max_alts=0):
        self.region_num = region_num
        self.prefix = prefix
        self.pool = '%s_2' %(self.prefix) if self.region_num % 2 == 0 else '%s_1' %(self.prefix)
        self.candidatePairs = candidatePairs
        self.fwdAlternates = fwdAlternates
        self.revAlternates = revAlternates

    @property
    def topPair(self):
        return self.candidatePairs[0]



class Primer(object):
    """A simple primer."""

    def __init__(self, direction, name, seq):
        self.direction = direction
        self.name = name
        self.seq = seq
        self.tm = calcTm(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6)
        self.homodimer = calcHomodimer(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm
        self.hairpin = calcHairpin(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm
        self.gc = 100.0 * (seq.count('G') + seq.count('C')) / len(seq)

    @property
    def length(self):
        return len(self.seq)


class CandidatePrimer(Primer):
    """A candidate primer for a region."""

    def __init__(self, direction, name, seq, start):
        super(CandidatePrimer, self).__init__(direction, name, seq)
        self.start = start
        self.percent_identity = 0
        self.alignments = []

    def align(self, references):
        #percent_identity = 0
        for ref in references:
            alignment = CAlignment(self, ref)
            self.alignments.append(alignment)
            #percent_identity += alignment.percent_identity
        # Calculate average percent identity
        #print(self.name, self.seq)
        idents = [i.percent_identity for i in self.alignments if i.percent_identity]
        if idents:
            self.percent_identity = sum(idents) / len(idents)
        #self.percent_identity = percent_identity / len(self.alignments)
        return(self)

    @property
    def end(self):
        if self.direction == 'LEFT':
            return self.start + self.length
        else:
            return self.start - self.length

class CandidatePrimerPair(object):
    """A pair of candidate primers for a region."""

    def __init__(self, left, right):
        self.left = left
        self.right = right
        # Calculate mean percent identity
        self.mean_percent_identity = (left.percent_identity + right.percent_identity) / 2

    @property
    def product_length(self):
        return self.right.start - self.left.start + 1


class Region(object):
    """A region that forms part of a scheme."""
    def __init__(self, region_num, chunk_start, candidate_pairs, references, prefix, max_alts=0):
        self.region_num = region_num
        self.prefix = prefix
        self.pool = '%s_2' %(self.prefix) if self.region_num % 2 == 0 else '%s_1' %(self.prefix)
        self.candidate_pairs = candidate_pairs
        self.alternates = []

        # Sort by highest scoring pair with the rightmost position
        self.candidate_pairs.sort(key=lambda x: (x.mean_percent_identity, x.right.end), reverse=True)

        # Align candidate pairs
        for pair in self.candidate_pairs:
            pair.left.align(references)
            pair.right.align(references)

        # Get a list of alts based on the alignments
        #for each in self.candidate_pairs[0].left.alignments:
        #    print(vars(each))
        left_alts = [each.aln_ref for each in self.candidate_pairs[0].left.alignments if each.aln_ref != self.candidate_pairs[0].left.seq]
        right_alts = [each.aln_ref for each in self.candidate_pairs[0].right.alignments if each.aln_ref != self.candidate_pairs[0].right.seq]

        left_alts_filt = [str(alt) for alt in left_alts if alt is not None]
        right_alts_filt = [str(alt) for alt in right_alts if alt is not None]

        # Get the counts for the alts to prioritise
        left_alts_counts = [(alt, left_alts_filt.count(alt)) for alt in set(left_alts_filt)]
        left_alts_counts.sort(key=lambda x: x[1], reverse=True)

        # Make tuples of unique primers and frequency and sort
        right_alts_counts = [(alt, right_alts_filt.count(alt)) for alt in set(right_alts_filt)]
        right_alts_counts.sort(key=lambda x: x[1], reverse=True)

        # For up to max_alts and if it occurs more than once generate a CandidatePrimer and add it to alternates list in Region
        for n, left_alt in enumerate(left_alts_counts):
            # max_alts is 1-indexed
            if n <= max_alts - 1 and left_alt[1] > 1:
                logger.debug('Found an alternate primer {} which covers {} reference sequences'.format(n, left_alt[1]))
                left = CandidatePrimer('LEFT', self.candidate_pairs[0].left.name + '_alt%i' %(n+1), left_alt[0], self.candidate_pairs[0].left.start).align(references)
                self.alternates.append(left)

        for n, right_alt in enumerate(right_alts_counts):
            if n <= max_alts - 1 and right_alt[1] > 1:
                logger.debug('Found an alternate primer {} which covers {} reference sequences'.format(n, right_alt[1]))
                right = CandidatePrimer('RIGHT', self.candidate_pairs[0].right.name + '_alt%i' %(n+1), right_alt[0], self.candidate_pairs[0].right.start).align(references)
                self.alternates.append(right)

    @property
    def top_pair(self):
        return self.candidate_pairs[0]

    @property
    def unique_candidates(self):
        unique_left = len(set(pair.left.seq for pair in self.candidate_pairs))
        unique_right = len(set(pair.right.seq for pair in self.candidate_pairs))
        return (unique_left, unique_right)


class CAlignment(object):
    """An seqan alignment of a primer against a reference."""

    def __init__(self, primer, ref):
        self.start = None
        self.end = None
        self.length = None
        self.percent_identity = None
        self.aln_query = None
        self.aln_ref = None
        self.aln_ref_comp = None
        self.ref_id = None
        self.mm_3prime = None
        self.cigar = None
        self.formatted_alignment = None

        if primer.direction == 'LEFT':
            alignment_result = adapter_alignment(str(ref.seq), str(primer.seq), [2, -1, -2, -1])
        elif primer.direction == 'RIGHT':
            alignment_result = adapter_alignment(str(ref.seq.reverse_complement()), str(primer.seq), [2, -1, -2, -1])
        result_parts = alignment_result.split(',')
        ref_start = int(result_parts[0])
        full_primer_percent_identity = float(result_parts[6])

        # If the read start is -1, that indicates that the alignment failed completely.
        if ref_start == -1 or full_primer_percent_identity < 70:
            return
            #self.percent_identity = 0.0
            #self.formatted_alignment = 'None'
        else:
            ref_end = int(result_parts[1]) + 1
            #aligned_region_percent_identity = float(result_parts[5])
            #full_primer_percent_identity = float(result_parts[6])

            if primer.direction == 'LEFT':
                self.start = ref_start
                self.end = ref_end
                self.length = self.end - self.start
            else:
                self.start = len(ref) - ref_start
                self.end = len(ref) - (int(result_parts[1]) + 1)
                self.length = self.start - self.end

            # Percentage identity for glocal alignment
            self.percent_identity = full_primer_percent_identity

            # Get alignment strings
            self.aln_query = result_parts[8][ref_start:ref_end]
            self.aln_ref = result_parts[7][ref_start:ref_end]
            self.aln_ref_comp = Seq(str(self.aln_ref)).complement()
            self.ref_id = ref.id
            self.mm_3prime = False

            # Make cigar
            self.cigar = ''
            for a, b in zip(self.aln_query, self.aln_ref):
                if a == '-' or b == '-':
                    self.cigar += ' '
                    continue
                if a != b:
                    self.cigar += '*'
                    continue
                else:
                    self.cigar += '|'

            # Format alignment
            short_primer = primer.name[:30] if len(primer.name) > 30 else primer.name
            short_ref = ref.id[:30] if len(ref.id) > 30 else ref.id
            self.formatted_alignment = "\n{: <30}5\'-{}-3\'\n{: <33}{}\n{: <30}3\'-{}-5\'".format(short_primer, self.aln_query, '', self.cigar, short_ref, self.aln_ref_comp)

            # Check 3' mismatches
            if set([self.aln_query[-1], self.aln_ref_comp[-1]]) in settings.MISMATCHES:
                self.mm_3prime = True
                self.percent_identity = 0
