import sys


def extract_lLTR(EDTA_annot_GFF3, output_lLTR_gff3, output_rLTR_gff3):
    lLTR_to_output = []
    rLTR_to_output = []

    with open(EDTA_annot_GFF3, mode='r') as lines:
        for line in lines:
            line_data = line.split("\t")
            if len(line_data) == 9 and line_data[2] == "long_terminal_repeat" and \
                    line_data[8].split(";")[0].startswith("ID=lLTR_"):
                lLTR_to_output.append(line)
            if len(line_data) == 9 and line_data[2] == "long_terminal_repeat" and \
                    line_data[8].split(";")[0].startswith("ID=rLTR_"):
                rLTR_to_output.append(line)

    with open(output_lLTR_gff3, mode='w') as output:
        for lLTR_line in lLTR_to_output:
            output.write(lLTR_line)
    with open(output_rLTR_gff3, mode='w') as output:
        for rLTR_line in rLTR_to_output:
            output.write(rLTR_line)


if __name__ == "__main__":
    # extract_lLTR(sys.argv[1], sys.argv[2])
    extract_lLTR(
        "E:/2022_Python_Projects\LTR Andropogoneae Scripts\S1_TE_Annotation_EDTA\Chrysopogon_serrulatus_EDTA_0426/GCA_015844335.1.fa.mod.EDTA.TEanno.gff3",
        "lLTR_0506.gff3", "rLTR_0506.gff3")
