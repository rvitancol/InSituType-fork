#' Get anchor stats
#' 
#' Compute the statistics used in finding anchor cells.
#' Often the anchor cell selection process will involve some trial-and-error. 
#' This function performs the computationally-expensive steps that only need to 
#' happen once.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param align_genes Logical, for whether to align the columns of the counts matrix and the rows of
#'  the profiles matrix based on their names. 
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param size Negative binomial size parameter to be used in likelihood calculation.
#' @param min_cosine Cells must have at least this much cosine similarity to a fixed profile to be used as an anchor.
#' @return A list with two elements: cos, the matrix of cosine distances;
#'  and llr, the matrix of log likelihood ratios of each cell under each cell type vs. the 2nd best cell type.
#' @importFrom lsa cosine
#' @export
get_anchor_stats <- function(counts, neg = NULL, bg = NULL, align_genes = TRUE,
                             profiles, size = 10, 
                             min_cosine = 0.3) {
  
  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg) && is.null(neg)) {
    stop("Must provide either bg or neg")
  }
  if (is.null(bg)) {
    s <- Matrix::rowMeans(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  ### align genes in counts and fixed_profiles
  if (align_genes && !is.null(profiles)) {
    sharedgenes <- intersect(rownames(profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    profiles <- profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) && length(lostgenes < 50)) {
      message(
        paste0(
          "The following genes in the count data are missing from fixed_profiles and will be omitted from anchor selection: ",
          paste0(lostgenes, collapse = ",")
        )
      )
    }
    if (length(lostgenes) > 50) {
      message(
        paste0(
          length(lostgenes),
          " genes in the count data are missing from fixed_profiles and will be omitted from anchor selection"
        )
      )
    }
  }
  
  # get cosine distances:
  cos <- matrix(NA, nrow(counts), ncol(profiles),
                dimnames = list(rownames(counts), colnames(profiles)))
  cos <- sapply(colnames(profiles), function(cell) {
    cos[, cell] <- apply(counts, 1, cosine, profiles[, cell])
  })
  
  # stats for which cells to get loglik on: 
  # get 3rd hightest cosines of each cell:
  cos3 <- apply(cos, 1, function(x) {
    return(x[order(x, decreasing = TRUE)[3]])
  })
  # get cells with sufficient cosine:
  cells_with_high_cos <- apply(cos, 1, max) > min_cosine
  # get logliks (only when the cosine similarity is high enough to be worth considering):
  logliks <- sapply(colnames(profiles), function(cell) {
    templl <- cos[, cell] * NA
    usecells <- which((cos[, cell] >= pmin(0.75 * min_cosine, cos3)) & cells_with_high_cos)
    if (length(usecells) > 0) {
      templl[usecells] <- Mstep(counts = counts[usecells, ], 
                                means = profiles[, cell, drop = FALSE],
                                cohort = rep("all", length(usecells)),
                                bg = bg[usecells], 
                                size = size, 
                                digits = 3, return_loglik = TRUE) 
      # scale the logliks by total counts:
      templl[usecells] <- templl[usecells] / rowSums(counts[usecells,, drop = FALSE])
    }
    return(templl)
  })
  # convert loglik to LLR:
  getsecondbest <- function(x) {
    return(x[order(x, decreasing = TRUE)][2])
  }
  secondbest <- apply(logliks, 1, getsecondbest)
  llr <- suppressWarnings(sweep(logliks, 1, secondbest, "-"))
  
  out <- list(cos = cos, llr = llr)
  return(out)
}


#' Choose anchor cells given anchor stats
#'
#' Starting with cosine distances and log likelihood ratios, choose anchor
#' cells.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param anchorstats Output from get_anchor_stats. Must provide either this or
#'   both cos and llr matrices.
#' @param cos Matrix of cosine distances from reference profiles. Cells in rows,
#'   cell types in columns.
#' @param llr Matrix of log likelihood ratios from reference profiles. Cells in
#'   rows, cell types in columns.
#' @param n_cells Up to this many cells will be taken as anchor points
#' @param min_cosine Cells must have at least this much cosine similarity to a
#'   fixed profile to be used as an anchor
#' @param min_scaled_llr Cells must have (log-likelihood ratio / totalcounts)
#'   above this threshold to be used as an anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than
#'   this many anchors will be discarded.
#' @return A vector holding anchor cell assignments (or NA) for each cell in the
#'   counts matrix
#' @export
choose_anchors_from_stats <-
  function(counts,
           neg = NULL,
           bg,
           anchorstats = NULL,
           cos = NULL,
           llr = NULL,
           n_cells = 500,
           min_cosine = 0.3,
           min_scaled_llr = 0.01,
           insufficient_anchors_thresh = 20) {
    
  
  if (is.null(anchorstats) && (is.null(cos) || is.null(llr))) {
    stop("Must provide either anchorstats or both cos and llr matrices.")
  }
  
  # get input:
  if (!is.null(anchorstats)) {
    cos <- anchorstats$cos
    llr <- anchorstats$llr
  }
  
  # apply thresholds:
  cos <- cos * (cos > min_cosine)
  llr <- llr * (llr > min_scaled_llr)
  
  # choose anchors for each cell type:
  anchors <- rep(NA, nrow(cos))
  names(anchors) <- rownames(cos)
  for (cell in colnames(cos)) {
    score <- llr[, cell] * cos[, cell] * (llr[, cell] > min_scaled_llr) * (cos[, cell] > min_cosine)
    topn <- order(score, decreasing = TRUE)[seq_len(min(n_cells, sum(score > 0, na.rm = TRUE)))]
    rm(score)
    anchors[topn] <- cell
  }
  
  # anchor consolidation: identify and remove anchor cells that are poor fits to
  # the mean anchor profile for the cell type:
  for (cell in setdiff(unique(anchors), NA)) {
    
    use <- (anchors == cell) & !is.na(anchors)
    # get centroid:
    if (!is.null(neg)) {
      mean_anchor_profile <- Estep(counts = counts[use, , drop = FALSE], clust = cell, neg = neg[use])
    } else {
      mean_anchor_profile <- Estep(counts = counts[use, , drop = FALSE], clust = cell, neg = bg[use])
    }
    
    # get anchors' cosine distances from centroid:
    newcos <- apply(counts[use, , drop = FALSE], 1, cosine, mean_anchor_profile)
    updated_anchors <- replace(anchors[use], (newcos < min_cosine), NA)
    anchors[names(updated_anchors)] <- updated_anchors
  }
  
  # remove all anchors from cells with too few total anchors:
  anchornums <- table(anchors)
  too_few_anchors <- names(anchornums)[anchornums <= insufficient_anchors_thresh]
  anchors[is.element(anchors, too_few_anchors)] <- NA
  
  if (length(setdiff(colnames(cos), unique(anchors))) > 0) {
    message(paste0("The following cell types had too few anchors and so are being removed from consideration: ",
                   paste0(setdiff(colnames(cos), unique(anchors)), collapse = ", ")))
  }
  
  if (all(is.na(anchors))) {
    warning("No anchor cells were selected - not enough cells met the selection criteria.")
    anchors <- NULL
  }
  return(anchors)  
}

  


  
#' Choose anchor cells
#'
#' Finds cells with very good fits to the reference profiles, and saves these
#' cells for use as "anchors" in the semi-supervised learning version of
#' nbclust.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param align_genes Logical, for whether to align the columns of the counts
#'   matrix and the rows of the profiles matrix based on their names.
#' @param profiles Matrix of reference profiles holding mean expression of genes
#'   x cell types. Input linear-scale expression, with genes in rows and cell
#'   types in columns.
#' @param size Negative binomial size parameter to be used in loglikelihood
#'   calculatoin
#' @param n_cells Up to this many cells will be taken as anchor points
#' @param min_cosine Cells must have at least this much cosine similarity to a
#'   fixed profile to be used as an anchor
#' @param min_scaled_llr Cells must have (log-likelihood ratio / totalcounts)
#'   above this threshold to be used as an anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than
#'   this many anchors will be discarded.
#' @return A vector holding anchor cell assignments (or NA) for each cell in the
#'   counts matrix
#' @importFrom lsa cosine
#' @export
find_anchor_cells <- function(counts, neg = NULL, bg = NULL, align_genes = TRUE,
                              profiles, size = 10, n_cells = 500, 
                              min_cosine = 0.3, min_scaled_llr = 0.01, 
                              insufficient_anchors_thresh = 20) {
  
  # get cos and llr stats:
  anchorstats <- get_anchor_stats(counts = counts,
                                   neg = neg,
                                   bg = bg, 
                                   align_genes = TRUE,
                                   profiles = profiles, 
                                   size = size, 
                                   min_cosine = min_cosine)  
  
  # select anchors based on stats:
  anchors <- choose_anchors_from_stats(counts = counts, 
                                       neg = neg,
                                       bg = bg,
                                       anchorstats = anchorstats, 
                                       cos = NULL, 
                                       llr = NULL, 
                                       n_cells = n_cells, 
                                       min_cosine = min_cosine, 
                                       min_scaled_llr = min_scaled_llr, 
                                       insufficient_anchors_thresh = insufficient_anchors_thresh) 
  
  return(anchors)
  
}



#' Filter anchor candidates via projection of reference profiles to anchor-derived UMAP
#' 
#' Calculates expression UMAP model for anchor candidates, then projects reference 
#' profiles to the anchor-derived UMAP and select anchor candidates within top 
#' nearest neighbors of the projected reference profiles of same cell type in the 
#' UMAP as the final anchor cells. 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param align_genes Logical, for whether to align the columns of the counts matrix and the rows of
#'  the profiles matrix based on their names. 
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param anchor_candidates Named vector of anchor candidates with cell_ID in name and corresponding cell type in values. 
#' @param nn_cells Number of top nearest neighbors to the projected reference profiles to be selected as final anchor cells. 
#' @return A list with two elements: anchors, a named vector for the final anchor cells;
#'  and anc_umap, the umap model derived from the expression profiles of anchor candidates.
#' @importFrom uwot umap umap_transform
#' @importFrom spatstat.geom ppp nncross
#' @export
filter_anchors_by_umap_projection <- function(counts, 
                                              neg = NULL, bg = NULL, 
                                              align_genes = TRUE,
                                              profiles, anchor_candidates, 
                                              nn_cells = 500) {
  # anchor candidates 
  cells_to_use <- names(anchor_candidates)[!is.na(anchor_candidates)]
  cells_to_use <- intersect(cells_to_use, rownames(counts))
  cts <- intersect(unique(anchor_candidates[cells_to_use]), colnames(profiles))
  
  if(length(cells_to_use)<20 | length(cts)<3){
    stop(sprintf("Only %d anchor candidates for %d cell types shared among `anchor_candidates`, `counts` and `profiles`. Must have at least 20 candidates for 3 cell types to use the projection based filtering.", 
         length(cells_to_use), length(cts)))
  } else {
    message(sprintf("Start filtering on %d anchor candidates for %d cell types.", 
                 length(cells_to_use), length(cts)))
  }
  
  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg) && is.null(neg)) {
    stop("Must provide either bg or neg")
  }
  if (is.null(bg)) {
    s <- Matrix::rowMeans(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  
  ### align genes in counts and fixed_profiles
  if (align_genes && !is.null(profiles)) {
    sharedgenes <- intersect(rownames(profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    profiles <- profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) && length(lostgenes < 50)) {
      message(
        paste0(
          "The following genes in the count data are missing from fixed_profiles and will be omitted from anchor selection: ",
          paste0(lostgenes, collapse = ",")
        )
      )
    }
    if (length(lostgenes) > 50) {
      message(
        paste0(
          length(lostgenes),
          " genes in the count data are missing from fixed_profiles and will be omitted from anchor selection"
        )
      )
    }
  }
  
  # net expression profiles of anchor candidates, cell x gene 
  netExpr <- pmax(sweep(counts[cells_to_use, ], 1, bg[cells_to_use], "-"), 0)
  netExpr <- as.matrix(netExpr)[, Matrix::colSums(netExpr)>0]
  
  sharedgenes <- intersect(rownames(profiles), colnames(netExpr))

  # get umap model for net expression of anchor candidates
  anc_umap <- uwot::umap(netExpr[, sharedgenes], metric = "cosine", ret_model = T)
  
  # projected ref in anc_umap
  projRef_umapcoord <- uwot::umap_transform(t(profiles)[, sharedgenes], anc_umap)
  
  # get top nearest neighbors for each projected reference within the anchors
  topNN <- as.matrix(spatstat.geom::nncross(
    # query point pattern for projected Ref
    X = spatstat.geom::ppp(x = projRef_umapcoord[, 1], 
                           y = projRef_umapcoord[, 2], 
                           range(projRef_umapcoord[, 1]), 
                           range(projRef_umapcoord[, 2]), 
                           marks = factor(rownames(projRef_umapcoord))), 
    # nearest neighbors in anchor data
    Y = spatstat.geom::ppp(x = anc_umap$embedding[, 1], 
                           y = anc_umap$embedding[, 2], 
                           range(anc_umap$embedding[, 1]), 
                           range(anc_umap$embedding[, 2]), 
                           marks = factor(anchor_candidates[rownames(anc_umap$embedding)])), 
    what = "which", k = 1:nn_cells))
  rownames(topNN) <- rownames(projRef_umapcoord)
  
  # get anchors within top Nearest Neighbors of consistent cell types
  anchors <- lapply(
    rownames(topNN),  
    function(ct){
      n500_cells <- rownames(anc_umap$embedding)[topNN[ct, ]]
      n500_cells <- intersect(n500_cells, 
                              names(anchor_candidates)[anchor_candidates == ct])
      ct_anchors <- rep(ct, length(n500_cells))
      names(ct_anchors) <- n500_cells
      return(ct_anchors)
    }
  )
  anchors <- do.call(c, anchors)
  
  message(sprintf("%d out of %d anchors are within top %d nearest neighbors of projected refProfiles of same cell types for %d cell types. ", 
                  length(anchors), length(cells_to_use), nn_cells, length(setdiff(unique(anchors), NA))))
  
  # full vector for final anchors including all cells in anchor_candidates
  anchors <- setNames(anchors[names(anchor_candidates)], 
                      names(anchor_candidates))
  
  out <- list(anchors = anchors,
              anc_umap = anc_umap)

  return(out)
}


#' Platform effect adjustment on reference profiles based on the expression profiles of high confident anchors 
#' 
#' Calculates gene-wise scaling factor between reference profiles and the provided 
#' cluster mean profiles of high confidence anchors for high expressor genes in 
#' either profiles, and then adjusts the reference profiles accordingly to get 
#' the platform effect corrected profiles
#' @param ref_profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param obs_profiles Matrix of observed cluster mean profiles derived from high confident anchors in genes x cell types format. 
#' @param min_ref Minimal expression values in reference profiles to define high expressors in reference, default = 1
#' @param cutoff_ref Cutoff in quantile to define high expressors of each cell type in reference profiles, default = 0.6
#' @param cutoff_obs Cutoff in quantile to define high expressors of each cell type in observed profiles, default = 0.9 
#' @param scaleFactor_limits Vector of 2 elements to define the range of per gene scale factors to be consider in the adjusted reference profiles, set to NULL if not to exclude based on scale factor value, default = c(0.01, 100)
#' @return A list with three elements: adj_profiles, Matrix of adjusted reference profiles in genes x cell types format; 
#' scaleFactor_DF, a data.frame with all genes of `ref_profiles` in rows and columns for `GeneName`, `scale_factor`, `type`; 
#'  and blacklist, a named vector with genes in name, scale factor in values, for genes with extreme scale factor beyond the range of `scaleFactor_limits` and thus not included in the returned `adj_profiles`, return NULL if no blacklist gene.
#' @importFrom reshape2 melt
#' @export
adjust_reference_for_platformEff <- function(ref_profiles, obs_profiles, 
                                             min_ref = 1, cutoff_ref = 0.6, 
                                             cutoff_obs = 0.9, 
                                             scaleFactor_limits = c(0.01, 100)){
  cts_to_use <- intersect(colnames(ref_profiles), colnames(obs_profiles))
  if(length(cts_to_use)<3){
    stop(sprintf('%d cell types shared between `ref_profiles` and ` obs_profiles`. Must have at least 3 shared cell types for platform effect adjustment.', length(cts_to_use)))
  }
  
  scaleFactor_DF <- data.frame(GeneName = rownames(ref_profiles), 
                               scale_factor = rep(1, nrow(ref_profiles)), 
                               type = rep(NA, nrow(ref_profiles)))
  rownames(scaleFactor_DF) <- scaleFactor_DF$GeneName
  
  sharedgenes <- intersect(rownames(ref_profiles), rownames(obs_profiles))
  scaleFactor_DF[sharedgenes, 'type'] <- 'shared'
  
  
  # get observed / reference efficiency ratio for each gene under each cell type ----
  # cell types x genes ratio matrix
  ratioMat <- mapply(function(x, y) x/y, 
                     as.data.frame(t(obs_profiles[sharedgenes, cts_to_use])),
                     as.data.frame(t(ref_profiles[sharedgenes, cts_to_use])))
  rownames(ratioMat) <- cts_to_use
  
  # flag genes with < min_ref or cutoff_ref quantile expression in given cell type for ref_profiles
  cutoff_perCT <- pmax(apply(ref_profiles[sharedgenes, cts_to_use], 2, quantile, cutoff_ref), min_ref)
  mat_low_ref <- (sweep(ref_profiles[sharedgenes, cts_to_use], 2, 
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
                     reshape2::melt(ref_profiles[genes_to_use, cts_to_use]), 
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
  adj_profiles <- ref_profiles
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
 
  
  outs <- list(adj_profiles = adj_profiles, 
               scaleFactor_DF = scaleFactor_DF, 
               blacklist = blacklist)
  
  return(outs)
  
}