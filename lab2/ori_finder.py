import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import collections

# Parameters from PDF 
K = 8
WINDOW_SIZE = 5000
STEP = 500

def get_kmer_counts(sequence, k):
    """Counts frequency of each k-mer in a sequence."""
    counts = collections.Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        counts[kmer] += 1
    return counts

def main():
    # Load the sequence (assuming genomic.fa is in the same folder)
    filename = "genomic.fa"
    
    try:
        record = next(SeqIO.parse(filename, "fasta"))
        seq = str(record.seq).upper()
    except FileNotFoundError:
        print(f"Error: Could not find {filename}. Make sure you ran 'gunzip'.")
        return

    print(f"Analyzing sequence length: {len(seq)} bp")
    
    # Simple logic to find most frequent k-mer in the whole genome first
    # (To simplify the plot, we often track specific 'signal' k-mers like DnaA boxes)
    print("Finding most frequent k-mers (this might take a moment)...")
    global_counts = get_kmer_counts(seq, K)
    top_kmers = global_counts.most_common(3)
    print(f"Top 3 k-mers: {top_kmers}")
    
    target_kmer = top_kmers[0][0] # Let's plot the most frequent one
    print(f"Plotting density for: {target_kmer}")

    # Sliding Window Analysis
    x_points = []
    y_points = []
    
    for i in range(0, len(seq) - WINDOW_SIZE, STEP):
        window = seq[i : i + WINDOW_SIZE]
        # Count occurrences of the target k-mer in this window
        count = window.count(target_kmer)
        
        x_points.append(i)
        y_points.append(count)

    # Plotting
    plt.figure(figsize=(12, 6))
    plt.plot(x_points, y_points, label=f"Count of {target_kmer}")
    plt.xlabel("Genome Position (bp)")
    plt.ylabel(f"Count in {WINDOW_SIZE}bp Window")
    plt.title(f"K-mer Enrichment: {target_kmer}")
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    output_file = "ori_plot.png"
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    main()
