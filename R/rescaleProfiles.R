
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
#' @param insufficient_anchors_thresh Cell types that end up with fewer than
#'   this many anchors will be discarded.
#' @param refinement Logical, flag for further anchor refinement via UMAP projection (default = FALSE)
#' @param blacklist vector of genes to be excluded for cell typing (default = NULL)
#' @param rescale Logical, flag for platform effect correction (default = FALSE)
#' @param refit Logical, flag for fitting reference profiles to anchors, run after rescale if rescale = TRUE (default = TRUE)
#' @return a list 
#' \describe{
#'     \item{updated_profiles}{a genes * cell types matrix for final updated reference profiles}
#'     \item{blacklist}{a vector of genes excluded from the final updated reference profiles}
#'     \item{anchors}{a named vector for final anchors used for reference profile update}
#'     \item{rescale_res}{a list of 5 elements, `rescaled_profiles`, `platformEff_statsDF`, `anchors`, `blacklist` and `lostgenes`, for platform effect correction outputs, return when rescale = TRUE}
#'     \item{refit_res}{a list of 2 elements, `refitted_profiles` and `anchors`, for anchor-based profile refitting outputs, return when refit = TRUE}
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
           insufficient_anchors_thresh = 20,
           refinement = FALSE, 
           blacklist = NULL,
           rescale = FALSE, 
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
                                   insufficient_anchors_thresh = insufficient_anchors_thresh, 
                                   align_genes = FALSE, 
                                   refinement = refinement) 
    } else {
      # test anchors are valid:
      if(is.null(names(anchors))){
        if (length(anchors) != nrow(counts)) {
          stop("anchors must have length equal to the number of cells (row) in counts")
        }
        names(anchors) <- rownames(counts)
      } else {
        anchors <- anchors[rownames(counts)]
        
        if (all(is.na(anchors))) {
          stop("The provided `anchors` have no shared cell_id in names as the rownames in `counts`.")
        } 
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
      outs[['rescale_res']] <- estimatePlatformEffects(counts[, sharedgenes], 
                                                       neg = neg, 
                                                       bg = bg, 
                                                       profiles = reference_profiles[sharedgenes, ], 
                                                       anchors = anchors)
      # add outliers from platform effect estimation, but not include lostgenes
      blacklist <- unique(c(blacklist, outs[['rescale_res']][['blacklist']]))
      if(!is.null(blacklist)){
        sharedgenes <- setdiff(sharedgenes, blacklist)
      }
      
      # if just rescale but no refit, blacklist contains lostgenes
      outs[['updated_profiles']] <- outs[['rescale_res']][['rescaled_profiles']]
      outs[['blacklist']] <- unique(c(blacklist, outs[['rescale_res']][['lostgenes']]))
    }
    
    
    ## step 3: second around of anchor selection based on corrected profiles
    if (rescale & refit){
      message("Second round of anchor selection given the rescaled reference profiles: ")
      genes_for_anchors <- rownames(outs[['rescale_res']][['rescaled_profiles']])
      anchors_second <- find_anchor_cells(counts = counts[, genes_for_anchors], 
                                          neg = neg, 
                                          bg = bg, 
                                          profiles = outs[['rescale_res']][['rescaled_profiles']], 
                                          size = nb_size, 
                                          n_cells = n_anchor_cells, 
                                          min_cosine = min_anchor_cosine, 
                                          min_scaled_llr = min_anchor_llr,
                                          insufficient_anchors_thresh = insufficient_anchors_thresh, 
                                          align_genes = FALSE, 
                                          refinement = refinement) 
      if (is.null(anchors_second)){
        message("No anchors were selected in the second around. Fall back to use first round of anchors for refitting. 
                Consider to 1. make anchor selection more generous, try `refinement = FALSE`. 2. select anchors by hand. 3. do `rescale` and `refit` in separate step. ")
      } else {
        # combine both rounds of anchors, which are in same cell_id orders as counts 
        anchors_second <- anchors_second[names(anchors)]
        combined_anchors <- sapply(
          seq_along(anchors), 
          function(idx){
            ct <- setdiff(unique(c(anchors[idx], anchors_second[idx])), NA)
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
      # refit original reference profiles given the anchors, will include lostgenes from platform effect estimation  
      refitted_profiles <- updateProfilesFromAnchors(counts = counts[, sharedgenes],  
                                                     neg = neg, 
                                                     bg = bg,
                                                     anchors = anchors, 
                                                     reference_profiles = reference_profiles[sharedgenes, ],
                                                     align_genes = FALSE, 
                                                     nb_size = nb_size)
      outs[['refit_res']] <- list(refitted_profiles = refitted_profiles, 
                                  anchors = anchors)
      outs[['updated_profiles']] <- refitted_profiles
      outs[['blacklist']] <- blacklist
    }
    
    outs[['anchors']] <- anchors
    
    
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


#' Platform effect adjustment on reference profiles based on the expression profiles of anchors 
#' 
#' Calculates gene-wise scaling factor between reference profiles and the observed profiles of the provided anchors.
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised
#'   clustering. Vector elements will be mainly NA's (for non-anchored cells)
#'   and cell type names for cells to be held constant throughout iterations.
#' @param blacklist vector of user-defined genes to be excluded for cell typing (default = NULL)
#' @return A list with three elements: 
#' \describe{
#'     \item{rescaled_profiles}{genes * cell types Matrix of rescaled reference profiles with platform effect corrected }
#'     \item{platformEff_statsDF}{a data.frame for statistics on platform effect estimation with genes in rows and columns for `Gene`, `Beta`, `beta_SE`.}
#'     \item{anchors}{a named vector of anchors used for platform effect estimation}
#'     \item{blacklist}{a vector of genes excluded from cell typing, including both outliers identified in platform effect estimation and the user-defined genes}
#'     \item{lostgenes}{a vector of genes excluded from platform effect estiamtion and the returned `rescaled_profiles`}
#' }
#' @description The general workflow would be: (1) extract the anchor cells from input; 
#' (2) Run poisson regression with anchor cells; (3) Filter user defined genes(if any) 
#' and genes with negative betas; (4) Re-scale Profile with Beta estimates. 
#' @importFrom data.table data.table melt
#' @importFrom Matrix colSums t rowSums
#' @importFrom parallel mclapply
#' @importFrom stats glm poisson
#' @export
estimatePlatformEffects <- 
  function(profiles,
           counts,
           neg, 
           bg=NULL, 
           anchors, 
           blacklist=NULL){
    
    #### Step1: clean up and prepare inputs, foucs on anchors only
    bg <- estimateBackground(counts = counts, neg = neg, bg = bg)
    # non-NA anchors 
    anchors <- anchors[!is.na(anchors)]
    anchors <- anchors[(anchors %in% colnames(profiles)) & (names(anchors) %in% rownames(counts))]
    if(length(anchors)<1){
      stop("The provided non-NA `anchors` are either not present in `counts` or have no shared cell types with `profiles`, check if missing names for cell_id in `anchors`.") 
    }
    
    # focus on shared genes and non-zero genes
    count_data <- alignGenes(counts = counts[names(anchors), ], 
                             profiles = profiles)
    # expanded cell types * genes matrix 
    reference_data <- Matrix::t(profiles[colnames(count_data), unname(anchors)])
    
    bad_genes <- (Matrix::colSums(count_data)==0 | Matrix::colSums(reference_data)==0)
    if(sum(bad_genes)>0){
      lostgenes <- colnames(count_data)[bad_genes]
      message(paste0("The following ", sum(bad_genes), 
                     " genes with zero counts in anchor cells or reference of corresponding cell types, are excluded from platform effect estimation: ", 
                     paste0(lostgenes, collapse = ",")))
      
      count_data <- count_data[,!bad_genes, drop = F]
      reference_data <- reference_data[,!bad_genes, drop = F]
    } else {
      lostgenes <- NULL
    }
    
    if(ncol(count_data)<1){
      stop("No shared genes in `anchors` with `counts` above zero.")
    }
    
    # fast conversion to data.frame for glm
    count_vector <- { dt = data.table::melt(data.table::data.table(count_data, keep.rownames = TRUE) , id.vars = c("rn"))} 
    colnames(count_vector) <- c("Cell_ID","Gene","Counts")
    
    reference_vector <- data.table::melt(data.table::data.table(reference_data, keep.rownames = TRUE) , id.vars = c("rn"))
    colnames(reference_vector) <- c("Cell_Type","Gene","Reference")
    
    ### Step2: Run poisson regression with link being Identity to estimate gene-level platform effect based on selected anchor cells
    query_DF <- as.data.frame(cbind(count_vector, reference_vector[, Gene:=NULL]))
    query_DF[['BG']] <- bg[count_vector$Cell_ID]
    # cell level scaling factor between obs vs. reference 
    query_DF[['Cell_SF']] <- Matrix::rowSums(count_data)/Matrix::rowSums(reference_data)
    
    
    PlatformEstimator <- function(gene_ID){
      GLM_Fit<- stats::glm(Counts ~  Reference :  Cell_SF + offset(BG) -1 , family= stats::poisson(link = "identity"),
                           data=query_DF[query_DF$Gene==gene_ID,],start=c(0))
      
      return(data.frame(Gene=gene_ID,
                        Beta=GLM_Fit$coefficients["Reference:Cell_SF"],
                        beta_SE=summary(GLM_Fit)$coefficients["Reference:Cell_SF","Std. Error"]))
    }
    
    PlatformEff <- parallel::mclapply(colnames(count_data), PlatformEstimator, mc.cores = numCores())
    PlatformEff <- as.data.frame(do.call(rbind, PlatformEff))
    rownames(PlatformEff) <- PlatformEff$Gene
    
    
    ### Step3: Filtering genes with negative betas and pre-specified by the user
    blacklist_addon <- as.vector(PlatformEff$Gene[PlatformEff$Beta<0])
    blacklist <- unique(c(blacklist,blacklist_addon))
    genes_to_keep <- as.vector(PlatformEff$Gene[!(PlatformEff$Gene %in%blacklist)])
    
    ### Step4: Rescale the raw reference profile
    rownames(PlatformEff) <- PlatformEff$Gene
    rescaled_profiles <- diag(PlatformEff[genes_to_keep,]$Beta) %*% profiles[genes_to_keep,]
    rownames(rescaled_profiles) <- genes_to_keep
    
    return(list(rescaled_profiles = rescaled_profiles,
                platformEff_statsDF = PlatformEff,
                anchors = anchors,
                blacklist = blacklist,
                lostgenes = lostgenes))
    
  }

