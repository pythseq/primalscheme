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

    #Stability of 3' ends
    def endStability(self, seq2, direction):
        if direction == 'fwd':
            return calcEndStability(self.seq, seq2, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm
        if direction == 'rev':
            return calcEndStability(Seq(self.seq).reverse_complement()._data, seq2, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

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
        penalty = 0

    def __eq__(self, other):
        return self.seq == other.seq

    def __hash__(self):
        return hash(self.seq)

    def calcPenalty(self):
        #Tm high
        if self.tm > settings.global_args['PRIMER_OPT_TM']:
            self.penalty += settings.global_args['PRIMER_WT_TM_GT'] * (self.tm - settings.global_args['PRIMER_OPT_TM'])
        #Tm low
        if self.tm < settings.global_args['PRIMER_OPT_TM']:
            self.penalty += settings.global_args['PRIMER_WT_TM_LT'] * (settings.global_args['PRIMER_OPT_TM'] - self.tm)
        #High GC
        if self.gc > settings.global_args['PRIMER_OPT_GC_PERCENT']:
            self.penalty += settings.global_args['PRIMER_WT_GC_PERCENT_GT'] * (self.gc - settings.global_args['PRIMER_OPT_GC_PERCENT'])
        #Low GC
        if self.gc < settings.global_args['PRIMER_OPT_GC_PERCENT']:
            self.penalty += settings.global_args['PRIMER_WT_GC_PERCENT_LT'] * (settings.global_args['PRIMER_OPT_GC_PERCENT'] - self.gc)
        #Length high
        if self.length > settings.global_args['PRIMER_OPT_SIZE']:
            self.penalty += settings.global_args['PRIMER_WT_SIZE_GT'] * (self.length - settings.global_args['PRIMER_OPT_SIZE'])
        #Length low
        if self.length < settings.global_args['PRIMER_OPT_SIZE']:
            self.penalty += settings.global_args['PRIMER_WT_SIZE_LT'] * (settings.global_args['PRIMER_OPT_SIZE'] - self.length)
        #Self any
        if (self.tm - 5) <= self.selfAny:
            self.penalty += settings.global_args['PRIMER_WT_SELF_ANY_TH'] * (self.selfAny - (self.tm - 5 - 1))
        elif (self.tm - 5) > self.selfAny:
            self.penalty += settings.global_args['PRIMER_WT_SELF_ANY_TH'] * (1/(self.tm - 5 + 1 - self.selfAny))
        #Self end
        if (self.tm - 5) <= self.selfEnd:
            self.penalty += settings.global_args['PRIMER_WT_SELF_END_TH'] * (self.selfEnd - (self.tm - 5 - 1))
        elif (self.tm - 5) > self.selfAny:
            self.penalty += settings.global_args['PRIMER_WT_SELF_END_TH'] * (1/(self.tm - 5 + 1 - self.selfEnd))
        #Hairpin
        if (self.tm - 5) <= self.hairpin:
            self.penalty += settings.global_args['PRIMER_WT_HAIRPIN_TH'] * (self.hairpin - (self.tm - 5 - 1))
        elif (self.tm - 5) > self.selfAny:
            self.penalty += settings.global_args['PRIMER_WT_HAIRPIN_TH'] * (1/(self.tm - 5 + 1 - self.hairpin))

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
