import sys
from Bio import SeqIO
from collections import Counter

# Parameters from PDF Step 8
K = 8
WINDOW_SIZE = 5000
STEP = 500

def get_gc_skew(sequence):
    """Calculates GC skew for the entire sequence using sliding windows."""
    skews = []
    indices = []
    
    # Calculate skew across the genome
    for i in range(0, len(sequence) - WINDOW_SIZE, STEP):
        window = sequence[i : i + WINDOW_SIZE]
        g = window.count('G')
        c = window.count('C')
        
        if g + c == 0:
            skew = 0
        else:
            skew = (g - c) / (g + c)
            
        skews.append(skew)
        indices.append(i)
        
    return skews, indices

def main():
    filename = "genomic.fa"
    print(f"Loading {filename}...")
    
    try:
        record = next(SeqIO.parse(filename, "fasta"))
        seq = str(record.seq).upper()
    except FileNotFoundError:
        print("Error: genomic.fa not found.")
        return

    # 1. Find the location of Minimum GC Skew
    print("Calculating GC Skew to find minimum...")
    skews, indices = get_gc_skew(seq)
    
    min_skew = min(skews)
    min_skew_index = skews.index(min_skew)
    # The actual genomic position is stored in indices[] at that index
    ori_candidate_pos = indices[min_skew_index]
    
    print(f"Minimum GC Skew ({min_skew:.4f}) found at position: {ori_candidate_pos} bp")
    print("This is the likely Origin of Replication (ORI).")

    # 2. Check for K-mer enrichment in that specific region
    # We look at a slightly larger window around the candidate position
    region_start = max(0, ori_candidate_pos - 1000)
    region_end = min(len(seq), ori_candidate_pos + WINDOW_SIZE + 1000)
    ori_region = seq[region_start : region_end]
    
    print(f"\nScanning candidate ORI region ({region_start}-{region_end} bp) for frequent {K}-mers...")
    
    # Count k-mers in this specific region
    kmer_counts = Counter()
    for i in range(len(ori_region) - K + 1):
        kmer = ori_region[i : i+K]
        kmer_counts[kmer] += 1
        
    # Display results
    top_kmers = kmer_counts.most_common(5)
    print("\nTop 5 k-mers in the ORI region:")
    for kmer, count in top_kmers:
        print(f"  {kmer}: {count} times")

if __name__ == "__main__":
    main()
