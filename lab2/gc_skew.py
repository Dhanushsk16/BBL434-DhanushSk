import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

# Parameters from PDF Step 8 (GC Skew section)
WINDOW_SIZE = 5000
STEP = 500

def calculate_gc_skew(sequence):
    """Calculates (G-C)/(G+C) for a given sequence string."""
    g = sequence.count('G')
    c = sequence.count('C')
    
    if g + c == 0:
        return 0.0
        
    return (g - c) / (g + c)

def main():
    filename = "genomic.fa"
    print(f"Loading {filename}...")
    
    try:
        record = next(SeqIO.parse(filename, "fasta"))
        seq = str(record.seq).upper()
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return

    print(f"Calculating GC Skew (Window={WINDOW_SIZE}, Step={STEP})...")
    
    skew_values = []
    positions = []
    
    # Sliding window loop
    for i in range(0, len(seq) - WINDOW_SIZE, STEP):
        window = seq[i : i + WINDOW_SIZE]
        skew = calculate_gc_skew(window)
        
        skew_values.append(skew)
        positions.append(i)

    # Plotting
    plt.figure(figsize=(12, 6))
    plt.plot(positions, skew_values, color='purple', label='GC Skew')
    
    # Adding a horizontal line at 0 for reference
    plt.axhline(0, color='black', linestyle='--', linewidth=0.5)
    
    plt.title("GC Skew Analysis")
    plt.xlabel("Genome Position (bp)")
    plt.ylabel("GC Skew (G-C)/(G+C)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    output_file = "gc_skew_plot.png"
    plt.savefig(output_file)
    print(f"Done! Plot saved to {output_file}")

if __name__ == "__main__":
    main()
