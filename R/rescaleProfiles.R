
#' Update reference profiles
#'
#' Update reference profiles using pre-specified anchor cells, or if no anchors
#' are specified, by first choosing anchor cells. Option to return reference 
#' profiles rescaled for platform effect and/or to return further refitted profiles 
#' based on the observed profiles of anchor cells.
#' @param reference_profiles Matrix of reference profiles, genes * cell types
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param anchors named vector giving "anchor" cell types with cell_id in names, 
#' for use in semi-supervised clustering. Vector elements will be mainly NA's 
#' (for non-anchored cells) and cell type names for cells to be held constant 
#' throughout iterations.
#' @param n_anchor_cells For semi-supervised learning. Maximum number of anchor
#'   cells to use for each cell type.
#' @param min_anchor_cosine For semi-supervised learning. Cells must have at
#'   least this much cosine similarity to a fixed profile to be used as an
#'   anchor.
#' @param min_anchor_llr For semi-supervised learning. Cells must have
#'   (log-likelihood ratio / totalcounts) above this threshold to be used as an
#'   anchor
#' @param refinement flag to further refine the anchors via UMAP projection (default = TRUE)
#' @param blacklist vector of genes to be excluded for cell typing 
#' @param rescale flag to rescale reference profiles for platform effect 
#' @param refit flag to refit reference profiles based on the anchor cells, run after rescale if rescale = TRUE
#' @return a list 
#' \describe{
#'     \item{updated_profiles}{a genes * cell types matrix for final updated reference profiles}
#'     \item{blacklist} {a vector of genes to be excluded from cell typing}
#'     \item{anchors} {a named vector for final anchors used for reference profile update}
#'     \item{rescale_res}{a list of 3 elements, `profiles`, `anchors` and `platform_effect`, for platform effect correction outputs, return when rescale = TRUE}
#'     \item{refit_res}{a list of 2 elements, `profiles` and `anchors`, for anchor-based profile refitting outputs, return when refit = TRUE}
#' }
#' @export
#' @examples
#' data("mini_nsclc")
#' data("ioprofiles")
#' counts <- mini_nsclc$counts
#' astats <- get_anchor_stats(counts = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  profiles = ioprofiles)
#'
#' # estimate per-cell bg as a fraction of total counts:
#' negmean.per.totcount <- mean(rowMeans(mini_nsclc$neg)) / mean(rowSums(counts))
#' per.cell.bg <- rowSums(counts) * negmean.per.totcount
#'
#' # now choose anchors:
#' anchors <- choose_anchors_from_stats(counts = counts, 
#'                                     neg = mini_nsclc$negmean, 
#'                                     bg = per.cell.bg,
#'                                     anchorstats = astats, 
#'                                     # a very low value chosen for the mini
#'                                     # dataset. Typically hundreds of cells
#'                                     # would be better.
#'                                     n_cells = 50, 
#'                                     min_cosine = 0.4, 
#'                                     min_scaled_llr = 0.03, 
#'                                     insufficient_anchors_thresh = 5)
#'
#' # The next step is to use the anchors to update the reference profiles:
#'
#' updateReferenceProfiles(reference_profiles = ioprofiles, 
#'                        counts = mini_nsclc$counts, 
#'                        neg = mini_nsclc$neg, 
#'                        bg = per.cell.bg,
#'                        anchors = anchors) 
updateReferenceProfiles <-
  function(reference_profiles,
           counts,
           neg,
           bg = NULL,
           nb_size = 10,
           anchors = NULL,
           n_anchor_cells = 2000,
           min_anchor_cosine = 0.3,
           min_anchor_llr = 0.01,
           refinement = TRUE, 
           blacklist = NULL,
           rescale = TRUE, 
           refit = TRUE) {
    
    if(!any(rescale, refit)){
      stop("At least one of `rescale` or `refit` must be TRUE to update the reference profiles.")
    }
    
    # align genes:
    sharedgenes <- intersect(colnames(counts), rownames(reference_profiles))
    if(!is.null(blacklist)){
      sharedgenes <- setdiff(sharedgenes, blacklist)
    }
    
    ## step 1: if no initial anchors are provided, select them automatically:
    if (is.null(anchors)) {
      message("Automatically selecting anchor cells with the best fits to fixed profiles.")
      
      anchors <- find_anchor_cells(counts = counts[, sharedgenes], 
                                   neg = neg, 
                                   bg = bg, 
                                   profiles = reference_profiles[sharedgenes, ], 
                                   size = nb_size, 
                                   n_cells = n_anchor_cells, 
                                   min_cosine = min_anchor_cosine, 
                                   min_scaled_llr = min_anchor_llr,
                                   insufficient_anchors_thresh = 0, 
                                   align_genes = FALSE, 
                                   refinement = refinement) 
    } else {
      # test anchors are valid:
      if (length(anchors) != nrow(counts)) {
        stop("anchors must have length equal to the number of cells (row) in counts")
      }
      if(is.null(names(anchors))){
        names(anchors) <- rownames(counts)
      } else {
        anchors <- anchors[rownames(counts)]
      }
      
      
      # refine existing anchors 
      if (refinement){
        anchors <- refineAnchors(counts = counts[, sharedgenes], 
                                 neg = neg, 
                                 bg = bg, 
                                 align_genes = FALSE,
                                 profiles = reference_profiles[sharedgenes, ],  
                                 anchor_candidates = anchors, 
                                 nn_cells = n_cells)
      }
    }
    
    if (is.null(anchors))  {
      stop("No anchors were selected. The algorithm can't run under these conditions. 
         Solutions include: 1. make anchor selection more generous, try `refinement = FALSE`. 2. select anchors by hand.")
    }
    
    outs <- list()
    
    ## step 2: rescale profiles for platform effect based on high confidence of anchors
    if(rescale){
      message("Rescale reference profiles for platform effect.")
      platformEff_res <- estimatePlatformEffects(counts[, sharedgenes], 
                                                 neg = neg, 
                                                 bg = bg, 
                                                 profiles = reference_profiles[sharedgenes, ], 
                                                 anchors = anchors)
      
      outs[['rescale_res']] <- list(profiles = platformEff_res[['adj_profiles']],
                                    anchors = anchors, 
                                    platform_effect = platformEff_res[['platform_effect']])  
      
      # update variables for the platform effect corrected outcomes 
      reference_profiles <- platformEff_res[['adj_profiles']]
      blacklist <- unique(c(blacklist, platformEff_res[['blacklist']]))
      if(!is.null(blacklist)){
        sharedgenes <- setdiff(sharedgenes, blacklist)
      }
    }
    
    
    ## step 3: second around of anchor selection based on corrected profiles
    if (rescale & refit){
      message("Second round of anchor selection given the rescaled reference profiles: ")
      anchors_second <- find_anchor_cells(counts = counts[, sharedgenes], 
                                          neg = neg, 
                                          bg = bg, 
                                          profiles = reference_profiles[sharedgenes, ], 
                                          size = nb_size, 
                                          n_cells = n_anchor_cells, 
                                          min_cosine = min_anchor_cosine, 
                                          min_scaled_llr = min_anchor_llr,
                                          insufficient_anchors_thresh = 0, 
                                          align_genes = FALSE, 
                                          refinement = refinement) 
      if (is.null(anchors_second)){
        message("No anchors were selected in the second around. Fall back to use first round of anchors for refitting. 
                Consider to 1. make anchor selection more generous, try `refinement = FALSE`. 2. select anchors by hand. 3. do `rescale` and `refit` in separate step. ")
      } else {
        # combine both rounds of anchors, which are in same cell_id orders as counts 
        use <- which(!(is.na(anchors) & is.na(anchors_second)))
        combined_anchors <- sapply(
          seq_along(anchors), 
          function(idx){
            ct <- setdiff(c(anchors[idx], anchors_second[idx]), NA)
            if(length(ct)!=1){
              # not to use if conflicts btw 2 rounds of anchor selection
              return(NA)
            }else {
              return(ct)
            }
          })
        names(combined_anchors) <- names(anchors)
        anchors <- combined_anchors
      }
    }
    

    ## step 4: refit the reference profiles using the second around of anchors 
    if(refit){
      reference_profiles <- updateProfilesFromAnchors(counts = counts[, sharedgenes],  
                                                    neg = neg, 
                                                    bg = bg,
                                                    anchors = anchors, 
                                                    reference_profiles = reference_profiles[sharedgenes, ],
                                                    align_genes = FALSE, 
                                                    nb_size = nb_size)
      outs[['refit_res']] <- list(profiles = reference_profiles, 
                                  anchors = anchors)
    }
    
    outs[['updated_profiles']] <- reference_profiles
    outs[['anchors']] <- anchors
    outs[['blacklist']] <- blacklist

    return(outs)
  }



#' Use anchor cells to update reference profiles, simply by taking the mean
#' profile of the anchors.
#'
#' Uses anchor cells to estimate platform effects / scaling factors to be
#' applied to the genes/rows of the reference profile matrix. Then uses Bayesian
#' math to update the individual elements on X.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided
#' @param bg Expected background
#' @param anchors Vector of anchor assignments
#' @param reference_profiles Matrix of expression profiles of pre-defined
#'   clusters, e.g. from previous scRNA-seq. These profiles will not be updated
#'   by the EM algorithm. Colnames must all be included in the init_clust
#'   variable.
#' @param align_genes Logical, for whether to align the counts matrix and the
#'   reference_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param max_rescaling Scaling factors will be truncated above by this value
#'   and below by its inverse (at 1/value and value)
#' @return \enumerate{ \item profiles: A profiles matrix with the rows rescaled
#' according to platform effects and individual elements updated further \item
#' scaling_factors: A vector of genes' scaling factors (what they were
#' multiplied by when updating the reference profiles). }
#' @export
updateProfilesFromAnchors <-
  function(counts,
           neg,
           bg = NULL, 
           anchors,
           reference_profiles,
           align_genes = TRUE,
           nb_size = 10,
           max_rescaling = 5) {
    bg <- estimateBackground(counts, neg, bg)
    use <- !is.na(anchors)
    updated_profiles <- Estep(counts = counts[use, ],
                              clust = anchors[use],
                              neg = bg[use])
    return(updated_profiles)
  }


#' Platform effect adjustment on reference profiles based on the expression profiles of high confident anchors 
#' 
#' Calculates gene-wise scaling factor between reference profiles and the provided 
#' cluster mean profiles of high confidence anchors for high expressor genes in 
#' either profiles, and then adjusts the reference profiles accordingly to get 
#' the platform effect corrected profiles
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised
#'   clustering. Vector elements will be mainly NA's (for non-anchored cells)
#'   and cell type names for cells to be held constant throughout iterations.
#' @param min_ref Minimal expression values in reference profiles to define high expressors in reference, default = 1
#' @param cutoff_ref Cutoff in quantile to define high expressors of each cell type in reference profiles, default = 0.6
#' @param cutoff_obs Cutoff in quantile to define high expressors of each cell type in observed profiles, default = 0.9 
#' @param scaleFactor_limits Vector of 2 elements to define the range of per gene scale factors to be consider in the adjusted reference profiles, set to NULL if not to exclude based on scale factor value, default = c(0.01, 100)
#' @return A list with three elements: 
#' \describe{
#'     \item{adj_profiles}{Matrix of adjusted reference profiles in genes x cell types format}
#'     \item{platform_effect}{a named vector of gene-specific scaling factor for genes in `profiles`, with names in gene name}
#'     \item{blacklist}{a named vector with genes in name, scale factor in values, for genes with extreme scale factor beyond the range of `scaleFactor_limits` and thus not included in the returned `adj_profiles`, return NULL if no blacklist gene.}
#' }
#' @importFrom reshape2 melt
#' @export
estimatePlatformEffects <- function(profiles, counts, neg, bg = NULL, anchors = NULL,
                                    min_ref = 1, cutoff_ref = 0.6, 
                                    cutoff_obs = 0.9, 
                                    scaleFactor_limits = c(0.01, 100)){
  # get observed profiles given the provided anchors
  cells_to_use <- intersect(names(anchors)[!is.na(anchors)], rownames(counts))
  bg <- estimateBackground(counts, neg, bg)
  obs_profiles <- Estep(counts = counts[cells_to_use, ], 
                        clust = anchors[cells_to_use],
                        neg = bg[cells_to_use])
  
  
  cts_to_use <- intersect(colnames(profiles), colnames(obs_profiles))
  if(length(cts_to_use)<3){
    stop(sprintf('%d cell types shared between `profiles` and ` obs_profiles`. Must have at least 3 shared cell types for platform effect adjustment.', length(cts_to_use)))
  }
  
  scaleFactor_DF <- data.frame(GeneName = rownames(profiles), 
                               scale_factor = rep(1, nrow(profiles)), 
                               type = rep(NA, nrow(profiles)))
  rownames(scaleFactor_DF) <- scaleFactor_DF$GeneName
  
  sharedgenes <- intersect(rownames(profiles), rownames(obs_profiles))
  scaleFactor_DF[sharedgenes, 'type'] <- 'shared'
  
  
  # get observed / reference efficiency ratio for each gene under each cell type ----
  # cell types x genes ratio matrix
  ratioMat <- mapply(function(x, y) x/y, 
                     as.data.frame(t(obs_profiles[sharedgenes, cts_to_use])),
                     as.data.frame(t(profiles[sharedgenes, cts_to_use])))
  rownames(ratioMat) <- cts_to_use
  
  # flag genes with < min_ref or cutoff_ref quantile expression in given cell type for profiles
  cutoff_perCT <- pmax(apply(profiles[sharedgenes, cts_to_use], 2, quantile, cutoff_ref), min_ref)
  mat_low_ref <- (sweep(profiles[sharedgenes, cts_to_use], 2, 
                        cutoff_perCT, '-') < 0)
  # flag high expressors in df
  high_genes <- sharedgenes[rowSums(mat_low_ref) < ncol(mat_low_ref)]
  scaleFactor_DF[high_genes, 'type'] <- 'high_ref'
  
  # flag genes with < cutoff_obs quantile expression in given cell type for obs_profiles
  cutoff_perCT <- apply(obs_profiles[sharedgenes, cts_to_use], 2, quantile, cutoff_obs)
  mat_low_obs <- (sweep(obs_profiles[sharedgenes, cts_to_use], 2, 
                        cutoff_perCT, '-') < 0)
  # flag high expressors in df
  high_genes2 <- sharedgenes[rowSums(mat_low_obs) < ncol(mat_low_obs)]
  scaleFactor_DF[high_genes2, 'type'] <- 'high_obs'
  scaleFactor_DF[intersect(high_genes, high_genes2), 'type'] <- 'high_both'
  
  # genes x cell types matrix to flag genes with low expression in both profiles
  mat_low_both <- mat_low_ref * mat_low_obs
  
  
  # genes x cell types matrix for obs vs. ref ratio of high expressors in either profiles
  cleanMat <- t(ratioMat)
  cleanMat[which(mat_low_both >0)] <- NA
  cleanMat <- cleanMat[rowSums(cleanMat, na.rm = T) >0, ]
  
  
  
  # get baseline observed / reference efficiency ratio for genes low in both profiles ----
  genes_to_use <- setdiff(sharedgenes, rownames(cleanMat))
  
  # get linear regression of those genes for obs vs. ref
  my_data <- reshape2::melt(obs_profiles[genes_to_use, cts_to_use])
  
  my_data <- merge(my_data , 
                   reshape2::melt(profiles[genes_to_use, cts_to_use]), 
                   by = c("Var1", "Var2"))
  colnames(my_data) <- c("GeneName", "CellType", "obs", "ref")
  
  # remove genes with zero in either ones
  my_data <- my_data[apply(my_data[, c("obs", "ref")], 1, min) >0, ]
  
  my.lm <- lm(obs~ref, my_data)
  baseline_ratio <- my.lm$coefficients[2]
  
  
  # per gene scale factor for high expressors across cell types ---
  scaleFactor <- rowMeans(cleanMat, na.rm = T)/ baseline_ratio
  scaleFactor_DF[names(scaleFactor), 'scale_factor'] <- unname(scaleFactor)
  
  # platform effect corrected reference
  adj_profiles <- profiles
  adj_profiles[names(scaleFactor), ] <- adj_profiles[names(scaleFactor), ] * scaleFactor
  
  # remove blacklist genes with extreme platform effect
  if(!is.null(scaleFactor_limits)){
    blacklist <- scaleFactor[scaleFactor < min(scaleFactor_limits) | scaleFactor > max(scaleFactor_limits)]
    
    scaleFactor_DF[names(blacklist), 'type'] <- 'beyond_limits'
    adj_profiles <- adj_profiles[!(rownames(adj_profiles) %in% names(blacklist)), ]
  } else {
    blacklist <- NULL
  } 
  
  if (length(blacklist)>0){
    message(sprintf('%d genes with extreme scale factor values are excluded in the `adj_profiles`: %s', 
                    length(blacklist), list(round(blacklist, 4))))
  }  else {
    blacklist <- NULL
  }
  
  genes_to_use <- scaleFactor_DF$GeneName[! scaleFactor_DF$type %in% c(NA, 'beyond_limits')]
  
  outs <- list(adj_profiles = adj_profiles, 
               platform_effect = setNames(scaleFactor_DF[genes_to_use, 'scale_factor'], 
                                       genes_to_use), 
               blacklist = blacklist)
  
  return(outs)
  
}

