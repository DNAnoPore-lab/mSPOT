# mSPOT

Given that pod5 file is too large, we uploaded a short fasta file used for testing after step 1-3. Reviewers can run Part B step 4-12 to verify data processing. “calls.fa” is the initial file.

## Part A. Design of domains of multiplexed DNAzyme and barcode sequence of Barcoded reporters

Python:

## 1. Design of Domains of Multiplexed DNAzyme

Generation of DNAzyme domains with constrained sequence length, GC content, and minimum Hamming distance to minimize cross-reactivity among multiplexed DNAzymes.

```python
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
```

---

## 2. Design of Barcode Sequences for Nanopore Sequencing

Generation of 24-nt barcode sequences for nanopore sequencing with constraints on Levenshtein distance, GC content, and homopolymer length to ensure robust barcode discrimination and sequencing accuracy.

```python
import random
import csv

random.seed(42)

num_sequences = 112
seq_length = 24
min_levenshtein_distance = 10
max_homopolymer = 2
gc_content_range = (0.4, 0.6)

nucleotides = ['A', 'T', 'C', 'G']

def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)

    for i, c1 in enumerate(s1):
        current_row = [i + 1]

        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)

            current_row.append(
                min(insertions, deletions, substitutions)
            )

        previous_row = current_row

    return previous_row[-1]

def has_long_homopolymer(seq, max_run):
    count = 1
    last = seq[0]

    for nt in seq[1:]:
        if nt == last:
            count += 1

            if count > max_run:
                return True
        else:
            count = 1
            last = nt

    return False

def gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq)

def generate_random_sequence():
    return ''.join(
        random.choices(
            nucleotides,
            k=seq_length
        )
    )

sequences = []
attempts = 0
max_attempts = 100000

while len(sequences) < num_sequences and attempts < max_attempts:

    attempts += 1

    candidate = generate_random_sequence()

    if has_long_homopolymer(candidate, max_homopolymer):
        continue

    gc_frac = gc_content(candidate)

    if not (
        gc_content_range[0]
        <= gc_frac
        <= gc_content_range[1]
    ):
        continue

    if all(
        levenshtein(candidate, existing)
        >= min_levenshtein_distance
        for existing in sequences
    ):
        sequences.append(candidate)

output_file = 'dna_sequences_levenshtein.csv'

with open(output_file, 'w', newline='') as csvfile:

    writer = csv.writer(csvfile)

    writer.writerow(['Index', 'Sequence'])

    for idx, seq in enumerate(sequences, 1):
        writer.writerow([idx, seq])

print(
    f"Generation complete: "
    f"{len(sequences)} sequences saved to {output_file}"
)
```
## Part B. Data processing of nanopore sequencing

Linux

Basecalling and read length filter:

1. dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/ --kit-name SQK-NBD114-24 --device cuda:all --no-trim > calls_notrim.bam

2. dorado demux -t 200 --emit-summary  --emit-fastq --kit-name SQK-NBD114-24 -v --output-dir calls.fq calls_notrim.bam

3. seqkit fq2fa calls.fq > calls.fa

4. seqkit seq -m55 -M77 calls.fa > calls_m55M77.fa

Mapping and counting:

5. bwa mem -k10 -x ont2d ref/ref.fa calls_m55M77.fa -t 200 > calls_m55M77.sam

6. samtools sort -@ 200 -O bam -o calls_m55M77.sorted.bam calls_m55M77.sam

7. samtools index calls_m55M77.sorted.bam

8. samtools idxstats calls_m55M77.sorted.bam > calls_m55M77.txt


Read length distribution:

9. seqkit seq -M150 calls.fa > calls_M150.fa

10. seqkit fx2tab -j 30 -l  -n -i -H  calls_M150.fa  |cut -f 2 > Length_calls_M150.txt
 
11. awk 'NR > 1 { print $1 }' Length_calls_M150.txt | sort | uniq -c | awk '{ print $2, $1 }' > Length_distribution_calls_M150.txt

12. Process Length_distribution_calls_M150.txt in excel to calculate percentage for each read length.

Version:
seqkit: 2.6.1
https://github.com/shenwei356/seqkit
bwa：0.7.17-r1188
https://github.com/lh3/bwa
samtools: 1.13
https://github.com/samtools/samtools

Installation guide can be obtained via their github. Their installation is very fast and can be used directly after installation.
