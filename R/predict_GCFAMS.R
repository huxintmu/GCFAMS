#' Calculate Gastric Cancer Fatty Acid Metabolism Score (GCFAMS)
#'
#' This function calculates the GCFAMS based on a given gene expression matrix.
#' The expression matrix should have rows as genes and columns as samples.
#' If the matrix is missing required genes, an error will be thrown.
#'
#' @param my_expr A numeric matrix with rows as genes and columns as samples.
#' @return A data frame with sample names and their Fatty Acid Metabolism Score.
#' @export
predict_GCFAMS <- function(my_expr) {
  required_genes <- c("CCDC80", "ZDHHC15", "AKT3", "CNRIP1", "GDF6", "LRCH2", "SVEP1",
                      "ZCCHC24", "MCC", "BCHE", "RERG", "RGAG4", "PDE2A", "CHRDL1",
                      "RUNX1T1", "CRTAC1", "NUDT10", "NAP1L2", "CYP1B1", "SFRP1", "ANTXR1",
                      "DAAM2", "FZD4", "SPON1", "FBN1", "OGN", "LGI4", "ZNF423",
                      "SOCS2", "HEG1", "SERPINF1", "KCNJ8", "C1QTNF7", "GALNT16", "MFAP5",
                      "TLL1", "MMRN2", "VTN", "MMRN1")
  missing_genes <- required_genes[!required_genes %in% rownames(my_expr)]
  if (length(missing_genes) > 0) {
    stop("Error: The following genes are missing in the expression matrix and are required for GCFAMS calculation: ",
         paste(missing_genes, collapse = ", "), ". Please provide a complete expression matrix.")
  }
  load(system.file("data/my_coef.rdata", package = "GCFAMS"))
  GCFAMS <- cbind(colnames(my_expr), my_GCFAMS=colSums(my_expr[required_genes,]*my_coef$coef))
  return(GCFAMS)
}


