import os
import sys
from Split_EDTATElib_into_LTRRT_and_NONLTRRT_FA import Extract_LTR_RT_fa


def create_TEs_Linc_task(ori_directory):
    current_dir = os.getcwd() + "/"

    if ori_directory.endswith("/"):
        ori_directory = ori_directory.rstrip("/")

    # Link TE annotation files generated from the EDTA pipeline.
    accession_name = ori_directory.split("/")[-1]
    TElib = ori_directory + "/" + accession_name + ".fa.mod.EDTA.TElib.fa"
    TEannot_gff3 = ori_directory + "/" + accession_name + ".fa.mod.EDTA.TEanno.gff3"
    TEannot_sum = ori_directory + "/" + accession_name + ".fa.mod.EDTA.TEanno.sum"
    sequence_translation_tab = ori_directory + "/" + accession_name + "_FA_sequence_id_translation/" + accession_name + ".fa.translation.tab"

    TElib_destination = current_dir + accession_name + "_TElib.fa"
    TEannot_gff3_destination = current_dir + accession_name + "_TEanno.gff3"
    TEannot_sum_destination = current_dir + accession_name + "_TEanno.sum"
    sequence_translation_tab_destination = current_dir + accession_name + "_Sequence_Translation.tab"

    os.system("ln -s %s %s" % (TElib, TElib_destination))
    os.system("ln -s %s %s" % (TEannot_gff3, TEannot_gff3_destination))
    os.system("ln -s %s %s" % (TEannot_sum, TEannot_sum_destination))
    os.system("ln -s %s %s" % (sequence_translation_tab, sequence_translation_tab_destination))

    Extract_LTR_RT_fa(TElib_destination)
    os.system("TEsorter ./EDTA_TElib_LTRRT.fa -db rexdb-plant -p 20")


if __name__ == "__main__":
    create_TEs_Linc_task(sys.argv[1])
    # create_TEs_Linc_task('/data5/Andropogoneae_LTR_Project/Step0_Genome_Download/GCA_019095995.1')
    # The original directory for EDTA results

