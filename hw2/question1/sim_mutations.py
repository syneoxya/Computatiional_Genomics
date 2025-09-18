import argparse
import random
from Bio import SeqIO

def mutate_sequence(seq, mutation_rate, seed):
    random.seed(seed)
    seq = list(seq)  # Convert to list for mutation
    num_mutations = int(len(seq) * mutation_rate)

    print(f"Total bases: {len(seq)}, Mutating {num_mutations} bases ({mutation_rate*100}%)")

    mutation_positions = random.sample(range(len(seq)), num_mutations)
    bases = ['A', 'C', 'G', 'T']

    for pos in mutation_positions:
        original_base = seq[pos].upper()
        if original_base not in bases:
            continue  # skip N or other ambiguous bases
        new_base = random.choice([b for b in bases if b != original_base])
        seq[pos] = new_base

    return ''.join(seq)

def main():
    parser = argparse.ArgumentParser(description="Introduce random substitutions into a FASTA sequence.")
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file with mutations')
    parser.add_argument('-m', '--mutation_rate', type=float, required=True, help='Mutation rate (e.g., 0.015 for 1.5%)')
    parser.add_argument('-s', '--seed', type=int, required=True, help='Random seed for reproducibility')

    args = parser.parse_args()

    with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            mutated_seq = mutate_sequence(str(record.seq), args.mutation_rate, args.seed)
            record.seq = mutated_seq
            SeqIO.write(record, outfile, 'fasta')

if __name__ == '__main__':
    main()