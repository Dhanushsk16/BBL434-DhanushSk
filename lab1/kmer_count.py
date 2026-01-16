from Bio import SeqIO

filename = "seq.fa"
k = 3  # Change this number to count different lengths (e.g., 4, 5)
kmer_counts = {}

try:
    for record in SeqIO.parse(filename, "fasta"):
        sequence = record.seq
        # Iterate through the sequence to find k-mers
        for i in range(len(sequence) - k + 1):
            kmer = str(sequence[i : i + k])
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1

    print(f"Dictionary of unique {k}-mers and their counts:")
    print(kmer_counts)

except FileNotFoundError:
    print(f"Error: The file '{filename}' was not found.")
