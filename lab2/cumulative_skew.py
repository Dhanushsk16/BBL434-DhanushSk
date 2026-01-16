import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

# Parameters from PDF Step 8
WINDOW_SIZE = 5000  # [cite: 31]
STEP = 500          # [cite: 31]

def main():
    filename = "genomic.fa"
    print(f"Loading {filename}...")
    
    try:
        record = next(SeqIO.parse(filename, "fasta"))
        seq = str(record.seq).upper()
    except FileNotFoundError:
        print("Error: genomic.fa not found. Please run 'gunzip genomic.fa.gz' first.")
        return

    print("Calculating Cumulative GC Skew...")
    
    skew_values = []
    positions = []
    current_cumulative_skew = 0  # This variable tracks the running total
    
    # We slide across the genome
    for i in range(0, len(seq) - WINDOW_SIZE, STEP):
        window = seq[i : i + WINDOW_SIZE]
        
        g = window.count('G')
        c = window.count('C')
        
        # Calculate simple skew for this window 
        if g + c == 0:
            local_skew = 0
        else:
            local_skew = (g - c) / (g + c)
            
        # ADD to the cumulative total 
        current_cumulative_skew += local_skew
        
        skew_values.append(current_cumulative_skew)
        positions.append(i)

    # Find the minimum point (The predicted ORI)
    min_skew = min(skew_values)
    min_index = skew_values.index(min_skew)
    min_pos = positions[min_index]
    
    print(f"Minimum Cumulative Skew found at approx: {min_pos} bp")

    # Plotting
    plt.figure(figsize=(12, 6))
    
    # Plot the cumulative line
    plt.plot(positions, skew_values, color='green', label='Cumulative GC Skew')
    
    # Mark the minimum point (Predicted ORI)
    plt.axvline(min_pos, color='red', linestyle='--', label=f'Predicted ORI ({min_pos})')
    
    plt.title(f"Cumulative GC Skew Plot\n(Minimum usually indicates Origin of Replication)")
    plt.xlabel("Genome Position (bp)")
    plt.ylabel("Cumulative Skew")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    output_file = "cumulative_skew_plot.png"
    plt.savefig(output_file)
    print(f"Proper plot saved to {output_file}")

if __name__ == "__main__":
    main()
