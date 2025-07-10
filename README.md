SeqGraph for Bioinformatics Applications
- Utilizes Graphviz library
- Programmed in Python, this program generates an overlap graph or de Bruijn graph by applying prefix-suffix matching to the given FASTA file.
- To run,

Input: Small FASTA file of sequenced reads, where the flags are the following options:
      -O for overlap graphs
      -D for deBruijn graphs
Output: Overlap graph or deBruijn graph written to new files

Examples:
$ python SeqGraph.py testOData.fa -O 
$ ls overlapG.png

$ python SeqGraph.py testDBData.fa -D
$ ls deBruijnG.png

