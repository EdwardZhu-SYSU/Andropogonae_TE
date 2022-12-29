
fpkm_tpm_calculator = function (gene_count, TE_count, gene_lentable, TE_lentable, outdir) {
  rm(list=ls())
  library(dplyr)

  gene_CPM = gene_count
  te_CPM_withdup = TE_count
  gene_length = gene_lentable
  te_length = TE_lentable


  # Remove duplicates for te_CPM entries, keeping the max CPM for overlapping samples
  colnames(te_CPM_withdup) = c("V1", "V2", "V3")
  te_CPM = te_CPM_withdup %>%
    group_by(V1) %>%
    summarize_all(max) %>%
    # Keep the max count number when encountering duplicates
    distinct(V1, .keep_all = T)
  colnames(te_CPM) = c("TE_id", "CPM", "TE_class")
  rownames(te_CPM) = te_CPM$TE_id
  #rm(te_CPM_withdup)


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


  # Cataloging TE types
  te_total$TE_class_group = "TE"

  DNA_TE_list = c("DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron")
  LTR_TE_list = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown")
  MITE_TE_list = c("MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", "MITE/DTT")

  te_total[which(te_total$TE_class %in% DNA_TE_list),7] = "DNA"
  te_total[which(te_total$TE_class %in% LTR_TE_list),7] = "LTR"
  te_total[which(te_total$TE_class %in% MITE_TE_list),7] = "MITE"

  te_total$TE_class_group = as.factor(te_total$TE_class_group)

  # Delete CPM column
  te_total = te_total[,-2]

  # Filter expressed TE FAMILIES (TPM>1)
  te_total_expressed = te_total[which(te_total$tpm > 1),]

  # ================================
  # Process gene expression data
  # ================================

  # Intersect, re-order entries
  colnames(gene_length) = c("gene_length")
  gene_length$gene_id = rownames(gene_length)

  colnames(gene_CPM) = c("CPM")
  gene_CPM$gene_id = rownames(gene_CPM)

  gene_align_order = intersect(rownames(gene_CPM), rownames(gene_length))
  gene_CPM = gene_CPM[gene_align_order,]
  gene_length = gene_length[gene_align_order,]

  rownames(gene_length) = 1:nrow(gene_length)
  rownames(gene_CPM) = 1:nrow(gene_CPM)

  rm(gene_align_order)


  # Calculation of FPKM and TPM values
  gene_total = gene_CPM
  gene_total$gene_length = gene_length$gene_length
  gene_total$fpkm = countToFpkm(gene_CPM$CPM, gene_length$gene_length)
  gene_total$tpm = fpkmToTpm(gene_total$fpkm)

  # Delete CPM column
  gene_total = gene_total[,-1]

  # Filter expressed genes (TPM>1)
  gene_total_expressed = gene_total[which(gene_total$tpm > 1),]

  write.table(gene_total, file = paste0(outdir, "gene_total_converted.tab"))  # outdir need to be ended with "/"
  write.table(te_total, file = paste0(outdir, "te_total_converted.tab"))

}