from Bio import SeqIO
import sys

gene = input("Gene: ")
pd = input("Product: ")
table = input("Translation Table: ")
fasta = sys.argv[1]
sequence_record = list(SeqIO.parse(fasta, "fasta"))

if len(sys.argv) == 3:
    sys.stdout = open(sys.argv[2], "w")

for i in sequence_record:
    frames = [i.seq, i.seq[1:], i.seq[2:]]
    reading_frames_trimmed_translated = [frame[:(len(frame)-len(frame)%3)].translate(table) for frame in frames]

    for n,j in enumerate(reading_frames_trimmed_translated):
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
            print("","","", "codon_start", n+1, sep = "\t")
            print("","","", "transl_table", table, sep = "\t")
