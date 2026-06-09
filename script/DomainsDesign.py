from itertools import product

def hamming_distance(seq1, seq2):
    """Compute the Hamming distance between two sequences"""
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def is_valid_sequence(candidate, sequences, min_distance):
    """Check if a new sequence meets the Hamming distance requirement within the set"""
    for seq in sequences:
        if hamming_distance(candidate, seq) < min_distance:
            return False
    return True

def gc_content(sequence):
    """Calculate the total number of G and C in the sequence"""
    return sequence.count('G') + sequence.count('C')

def find_gc_filtered_sequences(sequence_length, symbols, min_distance, gc_min, gc_max):
    """Find sequences satisfying GC content constraints and Hamming distance requirements"""
    all_sequences = list(product(symbols, repeat=sequence_length))
    valid_sequences = []

    for seq in all_sequences:
        if is_valid_sequence(seq, valid_sequences, min_distance):
            valid_sequences.append(seq)

    gc_filtered_sequences = [
        seq for seq in valid_sequences
        if gc_min <= gc_content(seq) <= gc_max
    ]

    return len(gc_filtered_sequences), ["".join(seq) for seq in gc_filtered_sequences]

sequence_length = 6
symbols = ['A', 'G', 'C', 'T']
min_distance = 4
gc_min = 2
gc_max = 4

count, sequences = find_gc_filtered_sequences(
    sequence_length,
    symbols,
    min_distance,
    gc_min,
    gc_max
)

print(f"Number of sequences satisfying the conditions: {count}")
print("The sequences are:")
print(sequences)
