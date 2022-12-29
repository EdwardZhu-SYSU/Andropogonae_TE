import csv
from Bio import SeqIO


class GFF:
    def __init__(self, gff_entry):
        self.seqid = gff_entry[0]
        self.source = gff_entry[1]
        self.type = gff_entry[2]
        self.start = int(gff_entry[3])
        self.end = int(gff_entry[4])
        self.score = gff_entry[5]
        self.strand = gff_entry[6]
        self.phase = gff_entry[7]
        self.attributes = gff_entry[8]


def GFF_and_fai_parse(TE_gff, fasta_index):
    gff_list = []
    seq_len_dict = {}  # Seq1 -> 500
    with open(TE_gff, mode='r') as OF:
        reader = csv.reader(OF, delimiter='\t')
        for row in reader:
            if (row[0][0] != '#') and (len(row) > 1):
                gff_entry = GFF(row)
                gff_list.append(gff_entry)

    with open(fasta_index, mode='r') as fai_fh:
        for indv_seq in fai_fh:
            seq_len_dict[indv_seq.split("\t")[0]] = indv_seq.split("\t")[1]

    return gff_list, seq_len_dict


def calc_feature_density(_seq_len_dict, _gff_list, windowsize, windowstep):
    density_list = []
    sequence_counter = 1
    for seq_name, seq_length in _seq_len_dict.items():
        print("Validating sequence " + str(sequence_counter) + ":" + seq_name + ", length: " + str(seq_length))
        if int(seq_length) >= windowsize:
            start = 0
            end = windowsize
            while end < int(seq_length):
                reduced_gff = [gene for gene in gff_list if
                               gene.seqid == seq_name and start <= gene.start <= end]
                num_genes = len(reduced_gff)
                entry = [seq_name, start, end, num_genes]
                density_list.append(entry)
                start = start + windowstep
                end = end + windowstep

        elif int(seq_length) < windowsize:
            reduced_gff = [gene for gene in gff_list if gene.seqid == seq_length]
            num_genes = len(reduced_gff)
            start = 0
            end = int(seq_length)
            entry = [seq_length, start, end, num_genes]
            print(entry)
            density_list.append(entry)

        sequence_counter += 1

    return density_list


#########################################################################################
if __name__ == "__main__":
    windowsize = 100000
    windowstep = 10000

    # Calalog GFF and fasta index
    gff_list, seq_len_dict = GFF_and_fai_parse(
        "E:/Bioinformatics_General_Data/TE_Module_B73v4/GCF_000005005.2_Sequence_Converted/GCF_000005005.2_TEanno_namefinal.gff3",
        "E:/Bioinformatics_General_Data/TE_Module_B73v4/GCF_000005005.2_Sequence_Converted/GCF_000005005.2_genome.fa.gz.fai"
    )
    # Specify Feature Name
    helitron_list = [entry for entry in gff_list if entry.type == 'helitron']
    LTR_list = [entry for entry in gff_list if entry.type in ['', '', '']]
    print("Data preparation completed, %d sequences cataloged." % len(seq_len_dict.keys()))


    gene_density = calc_feature_density(seq_len_dict, helitron_list, windowsize, windowstep)

    with open("./density.txt", 'w') as outFile:
        outfile_writer = csv.writer(outFile, delimiter='\t')
        outfile_writer.writerow(['scaffold', 'start', 'end', 'feature_name'])
        for entry in gene_density:
            outfile_writer.writerow(entry)
