#' Classify cells based on reference profiles
#' 
#' Supervised classification of cells. Each cell is assigned to the cell type 
#'  under which its observed expression profile is most likely. 
#' @param counts Counts matrix (or dgCMatrix), cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param bg Expected background
#' @param cohort Vector of cells' cohort memberships
#' @param reference_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param align_genes Logical, for whether to align the counts matrix and the reference_profiles by gene ID.
#' @return A list, with the following elements:
#' \enumerate{
#' \item clust: a vector given cells' cluster assignments
#' \item prob: a vector giving the confidence in each cell's cluster
#' \item profiles: Matrix of clusters' mean background-subtracted profiles
#' \item logliks: Matrix of cells' log-likelihoods under each cluster. Cells in rows, clusters in columns.
#' }
#' @export
insitutypeML <- function(counts, neg = NULL, bg = NULL, cohort = NULL, reference_profiles, reference_sds, nb_size = 10, assay_type, align_genes = TRUE) {
  
  if (any(rowSums(counts) == 0)) {
    stop("Cells with 0 counts were found. Please remove.")
  }
  
  # get vector of expected background:
  bg <- estimateBackground(counts = x, neg = neg, bg = bg)
  
  # align genes:
  if (align_genes) {
    x <- alignGenes(counts = x, profiles = reference_profiles)
    reference_profiles <- reference_profiles[colnames(x), ]

  }
  
  # prep cohort vector:
  if (is.null(cohort)) {
    cohort <- rep("all", length(bg))
  }
  
  # get logliks
  if(assay_type=="RNA"){
    logliks <- parallel::mclapply(asplit(reference_profiles, 2),
                                  lldist,
                                  xsd=NULL,
                                  mat = counts,
                                  bg = bg,
                                  size = size,
                                  assay_type=assay_type,
                                  mc.cores = numCores())
    
    logliks <- do.call(cbind, logliks)
  }else{
    
    logliks <- parallel::mcmapply(function(x, y){lldist(x=x, xsd=y, mat=counts, bg = bg, size = size, assay_type=assay_type)},
                                  asplit(reference_profiles, 2), asplit(reference_sds, 2), mc.cores = numCores())

  }
  
  # update logliks based on frequencies within cohorts:
  logliks <- update_logliks_with_cohort_freqs(logliks = logliks, 
                                              cohort = cohort, 
                                              minfreq = 1e-4, 
                                              nbaselinecells = 100) 
  
  # get remaining outputs
  clust <- colnames(logliks)[apply(logliks, 1, which.max)]
  names(clust) <- rownames(logliks)
  
  probs <- logliks2probs(logliks)
  prob <- apply(probs, 1, max)
  names(prob) <- names(clust)
  profiles_info <- Estep(counts, 
                    clust = clust,
                    neg = neg, 
                    assay_type=assay_type)
  
  profiles <- profiles_info$profiles
  sds <- profiles_info$sds
  
  # aligns profiles and logliks, removing lost clusters:
  logliks_from_lost_celltypes <- logliks[, !is.element(colnames(logliks), unique(clust)), drop = FALSE]
  logliks <- logliks[, is.element(colnames(logliks), clust), drop = FALSE]
  profiles <- profiles[, colnames(logliks), drop = FALSE]
  
  out <- list(clust = clust,
             prob = prob,
             profiles = profiles,
             sds = sds,
             logliks = round(logliks, 4),
             logliks_from_lost_celltypes = round(logliks_from_lost_celltypes, 4))
  return(out)    
}
