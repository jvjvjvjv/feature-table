from Bio import SeqIO
import sys

gene = input("Gene: ")
pd = input("Product: ")
table = input("Translation Table: ")
fasta = sys.argv[1]
sequence_record = list(SeqIO.parse(fasta, "fasta"))

for n,i in enumerate(sequence_record):
    frame_1 = i.seq.translate(5)
    frame_2 = i.seq[1:].translate(5)
    frame_3 = i.seq[2:].translate(5)
    for m,j in enumerate([frame_1, frame_2, frame_3]):
        if "*" not in j:
            print(">Feature ", i.name)
            print("<1", ">"+ str(len(i.seq)), "gene", sep = "\t")
            print("","","","gene", gene, sep = "\t")
            print("<1", ">"+ str(len(i.seq)), "CDS", sep = "\t")
            print("","","", "product", pd, sep = "\t")
            print("","","", "codon_start", m+1, sep = "\t")
            print("","","", "transl_table", table, sep = "\t")
