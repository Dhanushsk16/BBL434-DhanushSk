from Bio import SeqIO

# Define the filename
filename = "seq.mfa"

try:
    # SeqIO.parse returns an iterator. 
    # Converting it to a list loads all records, allowing us to count them easily.
    records = list(SeqIO.parse(filename, "fasta"))
    
    count = len(records)
    print(f"Number of FASTA records in '{filename}': {count}")

except FileNotFoundError:
    print(f"Error: The file '{filename}' was not found.")
