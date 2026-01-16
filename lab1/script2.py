from collections import Counter

seq = "ATGCGATCGATCGATCG"
k = 3

# Creates the dictionary in one go
counts = Counter([seq[i:i+k] for i in range(len(seq) - k + 1)])

print(dict(counts))
