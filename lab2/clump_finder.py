import sys
from Bio import SeqIO
from collections import Counter

# Parameters from the PDF Step 8
L = 1000  # Window length
K = 8     # k-mer length
T = 3     # Threshold (min occurrences)

def find_clumps(sequence, k, L, t):
    """
    Finds k-mers that appear at least t times in a window of length L.
    Uses a sliding window approach for efficiency.
    """
    clumps = set()
    n = len(sequence)
    
    # 1. Initialize the first window
    window_seq = sequence[0:L]
    kmer_counts = Counter()
    
    # Count k-mers in the first window
    for i in range(len(window_seq) - k + 1):
        kmer = window_seq[i : i+k]
        kmer_counts[kmer] += 1
    
    # Check for clumps in first window
    for kmer, count in kmer_counts.items():
        if count >= t:
            clumps.add(kmer)
            
    # 2. Slide the window across the genome
    # When we slide by 1 bp:
    # - We remove the k-mer at the very start of the OLD window
    # - We add the k-mer at the very end of the NEW window
    
    for i in range(1, n - L + 1):
        # Remove the k-mer that is sliding out (from the left)
        # previous window started at i-1
        removed_kmer = sequence[i-1 : i-1+k]
        kmer_counts[removed_kmer] -= 1
        
        # Add the k-mer that is sliding in (from the right)
        # new window ends at i+L
        added_kmer = sequence[i+L-k : i+L]
        kmer_counts[added_kmer] += 1
        
        # Check if the new k-mer makes a clump
        if kmer_counts[added_kmer] >= t:
            clumps.add(added_kmer)
            
    return clumps

def main():
    filename = "genomic.fa"
    print(f"Loading {filename}...")
    
    try:
        record = next(SeqIO.parse(filename, "fasta"))
        seq = str(record.seq).upper()
    except FileNotFoundError:
        print(f"Error: {filename} not found. Did you run 'gunzip'?")
        return

    print(f"Scanning genome ({len(seq)} bp) for clumps...")
    print(f"Parameters: L={L}, k={K}, t={T}")
    
    found_clumps = find_clumps(seq, K, L, T)
    
    print(f"\nTotal unique k-mers forming clumps: {len(found_clumps)}")
    print(f"First 10 clumps found: {list(found_clumps)[:10]}")

if __name__ == "__main__":
    main()
