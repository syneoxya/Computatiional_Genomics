reads = [
    "ATTCA",
    "ATTGA",
    "CATTG",
    "CTTAT",
    "GATTG",
    "TATTT",
    "TCATT",
    "TCTTA",
    "TGATT",
    "TTATT",
    "TTCAT",
    "TTCTT",
    "TTGAT"
]

k = 3
  # length of kmer
edges = set()
nodes = set()

for read in reads:
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        src = kmer[:k-1]
        dst = kmer[1:]
        nodes.add(src)
        nodes.add(dst)
        edges.add((src, dst))

# Output in DOT format
print('digraph debruijn {')
for src, dst in edges:
    print(f'    "{src}" -> "{dst}";')
print('}')
