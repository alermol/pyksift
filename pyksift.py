from Bio import SeqIO
import sys
from statistics import median
from pathlib import Path
from multiprocessing import Pool
from itertools import repeat


def parser_resolve_path(path):
    return Path(path).resolve()


def calc_AT_frac(seq):
    res = (seq.count('A') + seq.count('T')) / len(seq)
    return res


def count_kmers(seq, k, max_at):
    assert k < len(seq), "k-mer must be less than sequence length"
    assert k > 1, "k-mer must be greater than 1"
    kmers = {}
    for i in range(k, len(seq)):
        if calc_AT_frac(seq[(i - k):i]) > max_at:
            continue
        else:
            if seq[(i - k):i] not in kmers.keys():
                kmers[seq[(i - k):i]] = 1
            else:
                kmers[seq[(i - k):i]] += 1
    if len(kmers) == 0:
        return None
    res = tuple((kmer, i) for kmer, i in kmers.items())
    return sorted(res, key=lambda x: x[1], reverse=True)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(usage='python pyksift.py -o out.fa sequences.fa')

    parser.add_argument('infasta', help='input fasta file', type=parser_resolve_path)
    parser.add_argument('-o', help='name of output fasta', type=parser_resolve_path, metavar='outFasta', required=True)
    parser.add_argument('-t', help='number of CPUs threads to use (default: 1)', type=int, metavar='numthreads', default=1)
    parser.add_argument('-k', help='length of kmer to use (default: 10)', type=int, metavar='klegnth', default=10)
    parser.add_argument('-a', help='max AT proportion allowed in a kmer for use in calculating median kmer abundance (default: 0.7)',
                        type=float, metavar='maxAT', default=0.7)
    parser.add_argument('-j', help='number of top kmers accounted for median calculation (default: 10)', type=int, metavar='minTopKmers',
                        default=10)
    parser.add_argument('-m', help='minimum median of the top J kmers (default: 5)', type=int, metavar='minMedian', default=5)
    parser.add_argument('-l', help='minumum length of the sequence to consider (default: 1000)', type=int, metavar='minLength', default=1000)
    parser.add_argument('-n', help='include reads containing N', action='store_true')

    args = parser.parse_args()


    def worker(record, args):
        if (len(record.seq) < args.l) or (args.n and ('N' in str(record.seq))):
            return ''
        else:
            kmers = count_kmers(str(record.seq), args.k, args.a)[:(args.j + 1)]
            if kmers is None:
                return ''
            elif median([i[1] for i in kmers]) >= args.m:
                return f'>{record.id}\n{record.seq}\n'
            else:
                return ''

    print('Reading fasta file...', file=sys.stderr)
    records = list(SeqIO.parse(args.infasta, 'fasta'))

    print(f'Processing sequences in {args.t} threads...', file=sys.stderr)
    with Pool(processes=args.t) as pool:
        result = pool.starmap(worker, zip(records, repeat(args)))

    print('Writing selected sequences in file...', file=sys.stderr)
    with open(args.o, 'w') as outfile:
        outfile.write(''.join(result))

    print('Done', file=sys.stderr)
