import os
import sys
import logging
from primal import settings
from itertools import groupby

from Bio.Seq import Seq
from primer3 import calcTm, calcHairpin, calcHomodimer, calcHeterodimer
from primer3.bindings import calcEndStability
from Porechop.porechop.cpp_function_wrappers import adapter_alignment

logger = logging.getLogger('Primal Log')

class _primer(object):
    """A simple primer class"""

    def __init__(self, position, seq, direction):
        self.direction = direction
        self.name = None
        if self.direction == 'fwd':
            self.seq = seq
            self.start = position
            self.end = position + len(self.seq)
        elif self.direction == 'rev':
            self.seq = self.revComp(seq)
            self.start = position + len(self.seq)
            self.end = position

    #Rev comp reverse primers
    def revComp(self, seq):
        return Seq(seq).reverse_complement()._data

class _candidatePrimer(_primer):
    """A candidate primer super class"""

    def __init__(self, position, seq, direction):
        super(_candidatePrimer, self).__init__(position, seq, direction)
        self.penalty = 0
        self.calcPenalty()

    def __eq__(self, other):
        return self.seq == other.seq

    def __hash__(self):
        return hash(self.seq)

    def calcPenalty(self):
        #As per penalty routine described in http://primer3.ut.ee/primer3web_help.htm
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
        """
        #All of the default weights are 0 so ignore
        #Homodimer
        if (self.tm - 5) <= self.homodimer('fwd'):
            self.penalty += settings.global_args['PRIMER_WT_SELF_ANY_TH'] * (self.homodimer('fwd') - (self.tm - 5 - 1))
        elif (self.tm - 5) > self.homodimer('fwd'):
            self.penalty += settings.global_args['PRIMER_WT_SELF_ANY_TH'] * (1/(self.tm - 5 + 1 - self.homodimer('fwd')))
        #End stability
        if (self.tm - 5) <= self.endStability('fwd'):
            self.penalty += settings.global_args['PRIMER_WT_SELF_END_TH'] * (self.endStability('fwd') - (self.tm - 5 - 1))
        elif (self.tm - 5) > self.endStability('fwd'):
            self.penalty += settings.global_args['PRIMER_WT_SELF_END_TH'] * (1/(self.tm - 5 + 1 - self.endStability('fwd')))
        #Hairpin
        if (self.tm - 5) <= self.hairpin('fwd'):
            self.penalty += settings.global_args['PRIMER_WT_HAIRPIN_TH'] * (self.hairpin('fwd') - (self.tm - 5 - 1))
        elif (self.tm - 5) > self.hairpin('fwd'):
            self.penalty += settings.global_args['PRIMER_WT_HAIRPIN_TH'] * (1/(self.tm - 5 + 1 - self.hairpin('fwd')))
        """

    #Get reference coverage
    def queryMatch(self, references):
        if self.direction == 'fwd':
            return [ref.id for ref in references if ref[self.start:self.end].seq == self.seq]
        elif self.direction == 'rev':
            return [ref.id for ref in references if ref[self.end:self.start].seq == self.revComp(self.seq)]

    #Stability of homodimer
    @property
    def homodimer(self):
        return calcHomodimer(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

    #Stability of hairpin
    @property
    def hairpin(self):
        return calcHairpin(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

    #Stability of 3' ends
    @property
    def endStability(self):
        return calcEndStability(self.seq, self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

    #Tm
    @property
    def tm(self):
        return calcTm(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6)

    #Stability of last 5 3' bases
    @property
    def endStability(self):
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

class _primerPair(object):
    """A simple primer pair class"""

    def __init__(self, left, right):
        self.left = left
        self.right = right

    @property
    def productLength(self):
        return self.right.start - self.left.start + 1

class _candidatePrimerPair(_primerPair):
    """A primer pair super class"""

    def __init__(self, left, right):
        super(_candidatePrimerPair, self).__init__(left, right)
        #self.left.calcPenalty()
        #self.right.calcPenalty()
        self.pairPenalty = self.left.penalty + self.right.penalty

    @property
    def heterodimer(self):
        return calcHeterodimer(self.left.seq, self.right.revComp._data, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

    @property
    def endStability(self):
        return calcEndStability(self.left.seq, self.right.revComp._data, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm

    def fwdAlts(self, references, pairs, sortPairs, max_alts=5):
        #Update set of refs covered
        fwdCov = set(self.left.queryMatch(references[1:]))
        leftAlts = []
        #Generate left alts
        while len(fwdCov) < len(references[1:]):
            #Refs not covered
            fwdReq = [r for r in references[1:] if r.id not in fwdCov]
            #Alts with valid length product
            fwdAlts = [p.left for p in pairs if p.right == sortPairs[0].right]
            #Sort alts on required refs then all refs
            sortFwdAlts = sorted(fwdAlts, key=lambda x: (len(x.queryMatch(fwdReq)), len(x.queryMatch(references[1:]))), reverse=True)
            #Store alts
            leftAlts.append(sortFwdAlts[0])
            #Update refs covered
            lastCov = len(fwdCov)
            fwdCov.update(sortFwdAlts[0].queryMatch(references[1:]))
            if lastCov == len(fwdCov):
                return None
        return leftAlts

    def revAlts(self, references, pairs, sortPairs):
        revCov = set(sortPairs[0].right.queryMatch(references[1:]))
        rightAlts = []
        print(len(revCov), len(references[1:]))
        while len(revCov) < len(references[1:]):
            revReq = [r for r in references[1:] if r.id not in revCov]
            revAlts = [p.right for p in pairs if p.left == sortPairs[0].left]
            sortRevAlts = sorted(revAlts, key=lambda x: (len(x.queryMatch(revReq)), len(x.queryMatch(references[1:]))), reverse=True)
            rightAlts.append(sortRevAlts[0])
            lastCov = len(revCov)
            revCov.update(sortRevAlts[0].queryMatch(references[1:]))
            if lastCov == len(revCov):
                return None
        return rightAlts

    @property
    def productLength(self):
        return self.right.start - self.left.start + 1


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
