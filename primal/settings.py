global_args = {
    'PRIMER_OPT_SIZE': 22,
    'PRIMER_MIN_SIZE': 22,
    'PRIMER_MAX_SIZE': 30,
    'PRIMER_WT_SIZE_GT': 1.0,
    'PRIMER_WT_SIZE_LT': 1.0,
    'PRIMER_OPT_TM': 61.5,
    'PRIMER_MIN_TM': 60.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_WT_TM_GT': 1.0,
    'PRIMER_WT_TM_LT': 1.0,
    'PRIMER_MIN_GC': 30.0,
    'PRIMER_MAX_GC': 55.0,
    'PRIMER_OPT_GC_PERCENT': 50.0,
    'PRIMER_WT_GC_PERCENT_GT': 0.0,
    'PRIMER_WT_GC_PERCENT_LT': 0.0,
    'PRIMER_MAX_SELF_ANY_TH': 47.0,
    'PRIMER_WT_SELF_ANY_TH': 0.0,
    'PRIMER_MAX_SELF_END_TH': 47.0,
    'PRIMER_WT_SELF_END_TH': 0.0,
    'PRIMER_MAX_HAIRPIN_TH': 47.0,
    'PRIMER_WT_HAIRPIN_TH': 0.0,
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': 47.0,
    'PRIMER_PAIR_WT_COMPL_ANY_TH',: 0.0,
    'PRIMER_PAIR_MAX_COMPL_END_TH': 47.0,
    'PRIMER_PAIR_WT_COMPL_END_TH': 0.0,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_MAX_POLY_X': 5,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
}

seq_args = {
    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [-1, -1, -1, -1],
    }

MATCHES = [
    set(['A', 'T']),
    set(['C', 'G']),
    set(['G', 'T']),
    set(['C', 'T']),
    set(['T', 'T'])
],

MISMATCHES = [
    set(['A', 'A']),
    set(['A', 'C']),
    set(['C', 'C']),
    set(['G', 'A']),
    set(['G', 'G']),
]

NATIVE_DICT = {
    'NB01': 'AAGAAAGTTGTCGGTGTCTTTGTG',
    'NB02': 'TCGATTCCGTTTGTAGTCGTCTGT',
    'NB03': 'GAGTCTTGTGTCCCAGTTACCAGG',
    'NB04': 'TTCGGATTCTATCGTGTTTCCCTA',
    'NB05': 'CTTGTCCAGGGTTTGTGTAACCTT',
    'NB06': 'TTCTCGCAAAGGCAGAAAGTAGTC',
    'NB07': 'GTGTTACCGTGGGAATGAATCCTT',
    'NB08': 'TTCAGGGAACAAACCAAGTTACGT',
    'NB09': 'AACTAGGCACAGCGAGTCTTGGTT',
    'NB10': 'AAGCGTTGAAACCTTTGTCCTCTC',
    'NB11': 'GTTTCATCTATCGGAGGGAATGGA',
    'NB12': 'CAGGTAGAAAGAAGCAGAATCGGA',
}

SISPA_PRIMER = 'TTAACCGGCAGTGACACTGC'

RLBseq = 'TTTTTCGTGCGCCGCTTCAAC'

# Length 22, 30, 24.43
# Tm 59.96, 62.72, 61.33
# GC 33.33, 54.55, 45.71
# Longest homopolymer 5
# Ambiguous 0
