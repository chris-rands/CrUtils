import sys
import time
import logging
import argparse
from Bio import SeqIO

# Constants
STOP_SYMBOL = '*'  # the default for biopython


# Functions
def run_argparse():
    """User arguments"""
    parser = argparse.ArgumentParser(description='Translate fasta nucleotide sequence to fasta protein sequence, getting single longest sequence from start to stop codon. Default parameters for bacteria. @Chris Rands 2018')
    parser.add_argument('-i', '--input', required=True, help='fasta input file')
    parser.add_argument('-o', '--output', required=True, help='fasta output file')
    parser.add_argument('-sc', '--start_codons', nargs='*', default=['M', 'V', 'L'], help='start codons to use, default is bacterial start codonds')
    parser.add_argument('-tn', '--table_number', nargs='?', type=int, default=11, help='NCBI translation table number, default is for bacteria')
    parser.add_argument('-lf', '--log_file', nargs='?', default='faTranslate_log.txt', help='Log file name')
    return parser.parse_args()


def pad_seq(seq):
    """Pad sequence to multiple of 3 with Ns"""
    return {0: seq, 1: seq+'NN', 2: seq+'N'}[len(seq) % 3]


def translate_records(in_file, start_codons, table_number):
    """Translate records by taking single longest protein from start to stop codon for each record"""
    for record in SeqIO.parse(in_file, 'fasta'):
        # Translate into all six frames
        d = {}

        # Forward
        for i in range(3):
            d['fwd_frame_{}'.format(i+1)] = pad_seq(record.seq[i:]).translate(table=table_number, stop_symbol=STOP_SYMBOL)

        # Reverse
        rev_complement = record.seq.reverse_complement()
        for i in range(3):
            d['rvs_frame_{}'.format(i+4)] = pad_seq(rev_complement[i:]).translate(table=table_number, stop_symbol=STOP_SYMBOL)

        # Get all longest proteins for each frame from start to end codons
        d2 = {}
        for frame, each_seq in d.items():
            # For each of 6 frames
            for start_codon in start_codons:
                # For each of the codons
                prev_start_idx = -1
                # TODO: use prev_end_idx too to improve algorithm
                while True:
                    start_idx = each_seq.find(start_codon, prev_start_idx + 1)
                    if start_idx == -1:
                        break
                    end_idx = each_seq.find(STOP_SYMBOL, start_idx + 1) + 1
                    assert end_idx != -1
                    if end_idx == 0:
                        end_idx = None
                    tmp_seq = each_seq[start_idx:end_idx]
                    tmp_len = len(tmp_seq)
                    if frame in d2:
                        if tmp_len > d2[frame][0]:
                            d2[frame] = tmp_len, tmp_seq
                    else:
                        d2[frame] = tmp_len, tmp_seq
                    if end_idx is None:  # any start codon after this one to the end will be a shorter protein
                        break
                    prev_start_idx = start_idx

        if not d2:
            logging.error('No start codons found in record {}, exiting'.format(record.id))
            raise SystemExit

        # Get longest protein
        frame, (len_seq, seq) = max(d2.items(), key=lambda x: (x[1][0], x[1][1].startswith('M')))  # for a tie, take the sequence starting with M
        seq = seq.rstrip(STOP_SYMBOL)
        len_seq = len(seq)
        assert len_seq > 0
        if len_seq < 30:
            logging.warning('Fewer than 30 aa ({}), short protein for record {}'.format(len_seq, record.id))
        perc_of_original_len = round(float(len_seq) / (len(record.seq) // 3) * 100, 1)
        if perc_of_original_len < 50:
            logging.warning('Less than half ({}%) of nucleotide sequence translated for record {}'.format(perc_of_original_len, record.id))
        logging.info('Translating record: {} ; start codon: {} ;  frame: {} ; % of seq left: {}'.format(record.id, seq[0], frame, perc_of_original_len))
        record.seq = seq
        yield record


def main(in_file, start_codons, table_number, out_file):
    """Main work"""
    SeqIO.write(translate_records(in_file, start_codons, table_number), out_file, 'fasta')


# Run
if __name__ == '__main__':
    # Start time
    start_time = time.time()

    # Arguments setup
    args = run_argparse()

    # Logging setup
    logging.basicConfig(filename=args.log_file, filemode='a', level=logging.DEBUG)
    print('Logging info to: {}'.format(args.log_file))
    logging.info('Starting, argparse arguments: {}'.format(args))

    # Check python verion
    if sys.version_info[0] < 3:
        logging.warning('Only tested on Python 3, you are running an older Python version')

    # Check start codons
    args.start_codons = [codon.upper() for codon in args.start_codons]
    if 'M' not in args.start_codons:
        logging.error('"M" should be specified in start codons {}, exiting'.format(args.start_codons))
        raise SystemExit

    # Execute main
    main(args.input, args.start_codons, args.table_number, args.output)

    # Log time
    end_time = time.time()
    logging.info('Complete! Time to run {} seconds'.format(round(end_time-start_time)))
