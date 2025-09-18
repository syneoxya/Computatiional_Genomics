import argparse
import math

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

def get_kmers(seq, k):
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))

def compute_jaccard(A, B):
    intersection = len(A & B)
    union = len(A | B)
    return intersection / union if union else 0.0

def compute_ani(jaccard, k):
    # Exact formula
    if jaccard <= 0.0:
        return 0.0
    return jaccard ** (1 / k)

def compute_ani_approx(jaccard, k):
    # Approximation
    if jaccard <= 0.0:
        return 0.0
    return 1 + (math.log(jaccard) / k)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', required=True, help='Reference FASTA')
    parser.add_argument('-b', required=True, help='Mutated FASTA')
    parser.add_argument('-k', required=True, type=int, help='k-mer size')
    args = parser.parse_args()

    seq_a = read_fasta(args.a)
    seq_b = read_fasta(args.b)
    k = args.k

    kmers_a = get_kmers(seq_a, k)
    kmers_b = get_kmers(seq_b, k)
    jaccard = compute_jaccard(kmers_a, kmers_b)
    ani = compute_ani(jaccard, k)
    ani_approx = compute_ani_approx(jaccard, k)

    print(f'filename_a\tfilename_b\tjaccard\tani_exact\tani_approx')
    print(f'{args.a}\t{args.b}\t{jaccard:.6f}\t{ani:.6f}\t{ani_approx:.6f}')

if __name__ == '__main__':
    main()