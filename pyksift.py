import os
import shutil
import sys
import tempfile
from itertools import repeat
from multiprocessing import Pool
from pathlib import Path
from statistics import median
from subprocess import run

from Bio import SeqIO


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
    parser.add_argument('-y', help='additional parameters for yass, passing as a string (default: "-E 1.0E-5 -w 0")', type=str, default='-E 1.0E-5 -w 0',
                        metavar='yassParams')

    args = parser.parse_args()

    assert args.t > 0, 'Numbers of threads (-t) must be >0'
    assert args.k > 0, 'Kmer length (-k) must be >0'
    assert (args.a > 0) and (args.a <= 1), 'Max AT proportion (-a) must be in (0, 1]'
    assert args.j > 0, 'Number of top kmers (-j) must be >0'
    assert args.m > 0, 'minimum median of the top kmers (-m) must be >0'
    assert args.l > 0, 'minumum length of the sequence (-l) must be >0'
    assert args.k <= args.l, 'kmer length (-k) must be greater than minumum length of the sequence (-l)'


    if shutil.which('yass') is None:
        sys.exit('yass executable is not found in PATH variable. Exit.')
    else:
        yass_ex = shutil.which('yass')
        print('yass executable is found in PATH variable. Continue.', file=sys.stderr)

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
    result = [i for i in result if i != '']

    print(f'Writing {len(result)} selected sequences in file...', file=sys.stderr)
    with open(args.o, 'w') as outfile:
        outfile.write(''.join(result))

    print(f'Self-alignment of selected reads in {args.t} threads...', file=sys.stderr)


    def worker(seq, args, yass_ex_path):
        seq_name = seq.strip('>').split('\n')[0]
        with tempfile.NamedTemporaryFile(mode='w', dir=args.o.parent, suffix='.fa') as tfasta:
            tfasta.write(seq)
            yass_out = run(f'{yass_ex} -d 3 {args.y} {tfasta.name} {tfasta.name}', shell=True, capture_output=True)
        return [f'{seq_name}\t{i}' for i in yass_out.stdout.decode().split('\n')[1:] if i != '']


    with Pool(processes=args.t) as pool:
        yass_result = pool.starmap(worker, zip(result, repeat(args), repeat(yass_ex)))

    print(f'Writing yass output in file...', file=sys.stderr)
    os.remove(f'{args.o}.yass.tsv') if Path(f'{args.o}.yass.tsv').exists() else ''
    add_header = True
    with open(f'{args.o}.yass.tsv', 'a') as outfile:
        for output in yass_result:
            for line in output:
                outfile.write(f'{line}\n')

    print('Done', file=sys.stderr)
