import sys
from Bio import SeqIO

# Define the filename
filename = "seq.fa"

try:
    # Parse the FASTA file
    for record in SeqIO.parse(filename, "fasta"):
        print(f"Sequence ID: {record.id}")
        print(f"Sequence Length: {len(record.seq)}")

except FileNotFoundError:
    print(f"Error: The file '{filename}' was not found.")
