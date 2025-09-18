#!/usr/bin/env python3

import argparse
import math
import zlib

def clean_seq(seq):
    seq = seq.upper()
    return ''.join([c if c in "ACGT" else "N" for c in seq])

def read_fasta(path):
    with open(path) as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    return clean_seq(seq)

def get_modimizers(seq, k, m):
    mods = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        h = zlib.crc32(kmer.encode('utf-8')) & 0xffffffff
        if h % m == 0:
            mods.add(kmer)
    return mods

def compute_jaccard(A, B):
    intersection = len(A & B)
    union = len(A | B)
    return intersection / union if union > 0 else 0.0

def compute_ani(jaccard, k):
    if jaccard <= 0.0:
        return 0.0
    return jaccard ** (1 / k)

def compute_ani_approx(jaccard, k):
    if jaccard <= 0.0:
        return 0.0
    return 1 + (math.log(jaccard) / k)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', required=True, help='Reference FASTA')
    parser.add_argument('-b', required=True, help='Mutated FASTA')
    parser.add_argument('-k', required=True, type=int, help='k-mer size')
    parser.add_argument('-m', required=True, type=int, help='mod value for modimizer sampling')
    args = parser.parse_args()

    seq_a = read_fasta(args.a)
    seq_b = read_fasta(args.b)
    k = args.k
    m = args.m

    mods_a = get_modimizers(seq_a, k, m)
    mods_b = get_modimizers(seq_b, k, m)

    jaccard = compute_jaccard(mods_a, mods_b)
    ani = compute_ani(jaccard, k)
    ani_approx = compute_ani_approx(jaccard, k)

    print(f'filename_a\tfilename_b\tmod_value\tmodimizers_a\tmodimizers_b\tjaccard\tani_exact\tani_approx')
    print(f'{args.a}\t{args.b}\t{m}\t{len(mods_a)}\t{len(mods_b)}\t{jaccard:.6f}\t{ani:.6f}\t{ani_approx:.6f}')

if __name__ == '__main__':
    main()