from Bio import SeqIO
import sys

gene = input("Gene: ")
pd = input("Product: ")
table = input("Translation Table: ")
fasta = sys.argv[1]
sequence_record = list(SeqIO.parse(fasta, "fasta"))

for i in sequence_record:
    frame_1 = i.seq.translate(table)
    frame_2 = i.seq[1:].translate(table)
    frame_3 = i.seq[2:].translate(table)
    for m,j in enumerate([frame_1, frame_2, frame_3]):
        if "*" not in j[:-1]:
            if j[0] is "M":   #if first codon is start
                start = "1"
            else:
                start = "<1"
            if j[-1] is "*":   #if last codon is stop
                end = str(len(i.seq))
            else:
                end = ">" + str(len(i.seq))
            print(">Feature ", i.name)
            print(start, end, "gene", sep = "\t")
            print("","","","gene", gene, sep = "\t")
            print(start, end, "CDS", sep = "\t")
            print("","","", "product", pd, sep = "\t")
            print("","","", "codon_start", m+1, sep = "\t")
            print("","","", "transl_table", table, sep = "\t")
