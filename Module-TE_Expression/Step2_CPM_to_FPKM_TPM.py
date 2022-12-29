import os
import re

import rpy2.robjects as robjects


# Env req: /data5/condaenvs/rpy3dev

class IndividualAccession(object):
    def __init__(self):
        # Root folder for the project, containing all links.
        self._accession_name = os.getcwd().split("/")[-1]
        self._current_dir = os.getcwd() + "/"
        # bam-style alignment file for a certain species
        self._alignment_bam_dir_list = []

    def CPM_Normalizar(self):
        # ============================================================================
        # Split and format raw TEcount files
        def cntTable_simplifier(TEcount_cntTable, TEout_cntTable, Geneout_cntTable):
            with open(TEcount_cntTable, mode='r') as in_fh:
                with open(TEout_cntTable, mode='w') as out_fh1:
                    with open(Geneout_cntTable, mode='w') as out_fh2:

                        te_name_dict = {}

                        for entry in in_fh:
                            if re.match(r'^".*"$', entry.split("\t")[0]):  # Gene
                                out_fh2.write("\t".join([entry.split("\t")[0].rstrip("\"").lstrip("\""),
                                                         entry.split("\t")[1]]))

                            elif re.match(r'^.*:.*:.*:.*:.*$', entry.split("\t")[0]):  # TE_longver
                                out_fh1.write("\t".join(
                                    [":".join([entry.split("\t")[0].split(":")[0], entry.split("\t")[0].split(":")[1]]),
                                     entry.split("\t")[1].rstrip("\n"),
                                     entry.split("\t")[0].split(":")[-1]]) + "\n")

                            elif re.match(r'^.*:.*:.*$', entry.split("\t")[0]):  # TE_shortver
                                out_fh1.write(
                                    "\t".join([entry.split("\t")[0].split(":")[0], entry.split("\t")[1].rstrip("\n"),
                                               entry.split("\t")[0].split(":")[2]]) + "\n")
                            else:
                                print("Skipped line, ", entry)
            print("Done processing cntTable")

        # ============================================================================

        # ============================================================================
        # Calculation of FPKM and TPM values
        def fpkm_tpm_calculation(_gene_count, _TE_count, _gene_lentable, _TE_lentable, _outdir):
            """
            All inputs are in dataframes.
            :param sample_id:
            :param _gene_count:
            :param _TE_count:
            :param _gene_lentable:
            :param _TE_lentable:
            :return:
            """

            robjects.globalenv['gene_count'] = _gene_count
            robjects.globalenv['TE_count'] = _TE_count
            robjects.globalenv['gene_lentable'] = _gene_lentable
            robjects.globalenv['TE_lentable'] = _TE_lentable
            robjects.globalenv['outdir'] = _outdir

            r_norm_script = '''
            fpkm_tpm_calculator = function (gene_count, TE_count, gene_lentable, TE_lentable, outdir) {
              #rm(list=ls())
              library(dplyr)

              #gene_CPM = read.table(gene_count,header=F,row.names=1)
              
              te_CPM_withdup = read.table(TE_count,header=F)
              #gene_length = read.table(gene_lentable,header=F,row.names=1)
              te_length = read.table(TE_lentable,header=F,row.names=1)
              out_dir = outdir

              # Remove duplicates for te_CPM entries, keeping the max CPM for overlapping samples
              #head(te_CPM_withdup)
              colnames(te_CPM_withdup) = c("V1", "V2", "V3")
              te_CPM = te_CPM_withdup %>%
                group_by(V1) %>%
                #summarize_all(max) %>%
                # Keep the max count number when encountering duplicates
                distinct(V1, .keep_all = T)
              colnames(te_CPM) = c("TE_id", "CPM", "TE_class")
              rownames(te_CPM) = te_CPM$TE_id
              rm(te_CPM_withdup)
              
              print("te_CPM_non-overlapping frame: ")
              #print(te_CPM)
              
              
              


              # Define functions for CPM normalization
              countToFpkm = function(effcount, efflen){
                N = sum(effcount)
                exp(log(effcount)+log(1e9)-log(efflen)-log(N))
              }

              fpkmToTpm = function(fpkm){
                exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
              }

              # ================================
              # Process TE expression data
              # ================================

              # Intersect, re-order entries
              colnames(te_length) = c("fam_length")
              te_length$TE_id = rownames(te_length)

              te_align_order = intersect(rownames(te_CPM), rownames(te_length))
              te_CPM = te_CPM[te_align_order,]
              te_length = te_length[te_align_order,]

              rownames(te_length) = 1:nrow(te_length)

              te_CPM$CPM = as.numeric(te_CPM$CPM)
              te_length$fam_length = as.numeric((te_length$fam_length))

              rm(te_align_order)


              # Calculation of FPKM and TPM values
              te_total = te_CPM
              te_total$fam_length = te_length$fam_length
              te_total$fpkm = countToFpkm(te_CPM$CPM, te_length$fam_length)
              te_total$tpm = fpkmToTpm(te_total$fpkm)
              #print(te_total)


              # Cataloging TE types
              te_total$TE_class_group = "TE"

              DNA_TE_list = c("DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron")
              LTR_TE_list = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown")
              MITE_TE_list = c("MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", "MITE/DTT")

              te_total[which(te_total$TE_class %in% DNA_TE_list),7] = "DNA"
              te_total[which(te_total$TE_class %in% LTR_TE_list),7] = "LTR"
              te_total[which(te_total$TE_class %in% MITE_TE_list),7] = "MITE"

              te_total$TE_class_group = as.factor(te_total$TE_class_group)
              #print(te_total)

              # Delete CPM column
              te_total = te_total[,-2]

              # Filter expressed TE FAMILIES (TPM>1)
              te_total_expressed = te_total[which(te_total$tpm > 1),]
              
              write.table(te_total, file = paste0(out_dir, "TE.csv"), row.names=FALSE, col.names=TRUE, sep=",")
            }
            '''
            # Activate function
            robjects.r(r_norm_script)
            # Execute function
            robjects.r['fpkm_tpm_calculator'](_gene_count, _TE_count, _gene_lentable,
                                              _TE_lentable, _outdir)
            # Output1: Normalized FPKM & TPM value - genes
            # Module6_TE_expression/TEcount_final/SRR111111_gene.tab
            #
            # Output2: Normalized FPKM & TPM value - TEs
            # Module6_TE_expression/TEcount_final/SRR111111_te.tab

        # ============================================================================

        # Search for TEcount-generated cntTable files
        print("Filtering TEcount tables and normalizing CPM to FPKM & TPM.")
        for indv_file in os.listdir(self._current_dir + "Module6_TE_expression/TEcount_raw"):  # ex. 1.cntTable
            if indv_file.endswith(".cntTable"):
                print("Splitting: ", indv_file)
                cntTable_simplifier(self._current_dir + "Module6_TE_expression/TEcount_raw/" + indv_file,
                                    self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                                        ".cntTable") + "_TE.count",
                                    self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                                        ".cntTable") + "_gene.count",
                                    )
                print("Normalizing CPM: ", indv_file)
                fpkm_tpm_calculation(self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                    ".cntTable") + "_gene.count",
                                     self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                                         ".cntTable") + "_TE.count",
                                     self._current_dir + "Module6_TE_expression/genome_gene.length",
                                     self._current_dir + "Module6_TE_expression/TE_family.length",
                                     self._current_dir + "Module6_TE_expression/TEcount_final/" + indv_file.rstrip(
                                         ".cntTable") + "_"
                                     )
                print("Done:", indv_file)
                # break  # For diagnosis
        print("All done.")


mainAccession = IndividualAccession()
mainAccession.CPM_Normalizar()
