"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Dr Josh Quick
Contributions from Andrew Smith
www.github.com/aresti/primalscheme

This module contains the main script for PrimalScheme.
It is executed when the user runs 'primalscheme' after installation,
or 'primal.py' (directly from the source directory).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Alphabet import AlphabetEncoder, _verify_alphabet, IUPAC, generic_dna
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from .multiplex_reporting import MultiplexReporter

logger = logging.getLogger('Primal Log')


def main():
    """
    Entry point for primalscheme. Proccesses args, sets up logging and calls
    the appropriate scheme function.
    """

    args = get_arguments()
    output_path = args.output_path
    setup_logging(output_path, debug=args.debug, prefix=args.prefix)

    try:
        check_output_dir(output_path, force=args.force)

        logger.info('PrimalScheme started...')
        for arg in vars(args):
            logger.debug('{}: {}'.format(arg, str(vars(args)[arg])))
        # Run
        args.func(args)
    except Exception as e:
        logger.error('Error: {}'.format(e))
        sys.exit()


def multiplex(args):
    """
    Multipex scheme runner.
    """

    user_references = parse_user_input_fasta(args.fasta)

    for record in user_references:
        logger.info('Reference: {}'.format(record.id))

    consensus_sequences = find_consensus_sequences(
        user_references, args.output_path, args.prefix)

    scheme = MultiplexReporter(
        consensus_sequences, args.amplicon_length, min_overlap=args.min_overlap,
        max_gap=args.max_gap, max_alts=args.max_alts,
        max_candidates=args.max_candidates, step_size=args.step_size,
        max_variation=args.max_variation, prefix=args.prefix)

    scheme.write_all(args.output_path)


def find_consensus_sequences(user_references, output_path, prefix):
    """
    remove gaps and make uppercase
    """

    cleaned_poa_input = os.path.join(output_path, '{}.cleaned_poa_input.fasta'.format(prefix))
    SeqIO.write(map(SeqRecord.lower, user_references), cleaned_poa_input, format='fasta')

    basedir = '/Users/josh/bioinfo'
    poa_cmd = [basedir + '/bin/poaV2/poa', '-read_fasta', cleaned_poa_input \
    , basedir + '/bin/poaV2/blosum80.mat', '-do_global', \
    '-hb', '-best', '-pir', '/dev/stdout']

    poa_output = subprocess.run(poa_cmd, capture_output=True, check=True)

    with tempfile.NamedTemporaryFile() as temp_file:
        temp_file.write(poa_output.stdout)
        records = list(
            r.upper() for r in SeqIO.parse(temp_file.name, 'fasta'))

    return records[len(user_references):]


def parse_user_input_fasta(file_path):
    """
    Parse and validate the user multi-fasta file.
    """

    alphabet = AlphabetEncoder(IUPAC.unambiguous_dna, 'N')
    references = list(SeqIO.parse(file_path, 'fasta', alphabet=alphabet))  # may raise
    references.sort(key=len)

    # Check for too few or too many references
    if not (1 <= len(references) <= 100):
        raise ValueError('Between 1 and 100 reference genomes are required.')

    # Check that the shortest length is at least 0.9x the longest length
    shortest = references[0]
    longest = references[-1]
    if len(shortest) < 0.9 * len(longest):
        raise ValueError(
            'Your shortest reference must be at least 90% as long as your'
            'longest reference.'
        )

    # Check for a valid alphabet
    if any(not _verify_alphabet(r.seq) for r in references):
        raise ValueError(
            'One or more of your fasta sequences contain invalid '
            "nucleotide codes. The supported alphabet is '{}'. "
            'Ambiguity codes and gaps are not currently supported.'
            .format(alphabet.letters)
        )

    return references


def setup_logging(output_path, debug=False, prefix='primalscheme'):
    """
    Setup logging output and verbosity.
    """

    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    fh = logging.FileHandler(
        os.path.join(output_path, '{}.log'.format(prefix)))
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh_formatter = logging.Formatter('%(message)s')
    sh.setFormatter(sh_formatter)
    logger.addHandler(sh)


def check_output_dir(output_path, force=False):
    """
    Check for an existing output dir, require --force to overwrite.
    Otherwise, create a new one.
    """

    if os.path.isdir(output_path) and not force:
        raise IOError('Directory exists add --force to overwrite')
    if not os.path.isdir(output_path):
        os.mkdir(output_path)


def get_arguments():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(
        prog='primalscheme',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    # Standard multiplex scheme
    parser_scheme = subparsers.add_parser(
        'multiplex', help='Multiplex PCR scheme')
    parser_scheme.add_argument(
        'fasta', help='FASTA file')
    parser_scheme.add_argument(
        'prefix', help='Prefix')
    parser_scheme.add_argument(
        '--amplicon-length', type=int, default=400,
        help='Amplicon length (default: %(default)i)')
    parser_scheme.add_argument(
        '--min-overlap', type=int, default=0,
        help='Minimum overlap length (default: %(default)i)')
    parser_scheme.add_argument(
        '--max-gap', type=int, default=200,
        help='Maximum gap to introduce before failing (default: %(default)i)')
    parser_scheme.add_argument(
        '--max-alts', type=int, default=2,
        help='Maximum number of alternate primers to output '
        '(default: %(default)i)')
    parser_scheme.add_argument(
        '--max-candidates', type=int, default=10,
        help='Maximum candidate primers (default: %(default)i)')
    parser_scheme.add_argument(
        '--step-size', type=int, default=11,
        help='Step size when moving left or right (default: %(default)i)')
    parser_scheme.add_argument(
        '--max-variation', type=float, default=0.1,
        help='Variation in allowed product length (default: %(default)i)')
    parser_scheme.add_argument(
        '--output-path', default='./',
        help='Output directory to save files (default: %(default)s)')
    parser_scheme.add_argument(
        '--force', action='store_true', help='Force overwrite')
    parser_scheme.add_argument(
        '--debug', action='store_true', help='Verbose logging')
    parser_scheme.set_defaults(func=multiplex)

    # Generate args
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
