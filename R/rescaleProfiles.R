
#' Update reference profiles
#'
#' Update reference profiles using pre-specified anchor cells, or if no anchors
#' are specified, by first choosing anchor cells
#' @param reference_profiles Matrix of reference profiles, genes * cell types
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised
#'   clustering. Vector elements will be mainly NA's (for non-anchored cells)
#'   and cell type names for cells to be held constant throughout iterations.
#' @param n_anchor_cells For semi-supervised learning. Maximum number of anchor
#'   cells to use for each cell type.
#' @param min_anchor_cosine For semi-supervised learning. Cells must have at
#'   least this much cosine similarity to a fixed profile to be used as an
#'   anchor.
#' @param min_anchor_llr For semi-supervised learning. Cells must have
#'   (log-likelihood ratio / totalcounts) above this threshold to be used as an
#'   anchor
#' @return updated reference profiles
#' @export
updateReferenceProfiles <-
  function(reference_profiles,
           counts,
           neg,
           bg = NULL,
           nb_size = 10,
           anchors = NULL,
           n_anchor_cells = 2000,
           min_anchor_cosine = 0.3,
           min_anchor_llr = 0.01) {
    
  
  
  ## step 1: if no anchors are provided, select them automatically:
  if (is.null(anchors)) {
    message("automatically selecting anchor cells with the best fits to fixed profiles")
    # align genes:
    sharedgenes <- intersect(colnames(counts), rownames(reference_profiles))
    anchors <- find_anchor_cells(counts = counts[, sharedgenes], 
                                 neg = neg, 
                                 bg = bg, 
                                 profiles = reference_profiles[sharedgenes, ], 
                                 size = nb_size, 
                                 n_cells = n_anchor_cells, 
                                 min_cosine = min_anchor_cosine, 
                                 min_scaled_llr = min_anchor_llr,
                                 insufficient_anchors_thresh = 0) 
  }
  
  if (is.null(anchors))  {
    stop("No anchors were selected. The algorithm can't run under these conditions. 
         Solutions include: 1. make anchor selection more generous. 2. select anchors by hand.")
  }
  
  # test anchors are valid:
  if (!is.null(anchors) && (length(anchors) != nrow(counts))) {
    stop("anchors must have length equal to the number of cells (row) in counts")
  }
  names(anchors) <- rownames(counts)
  
  ## step 2: use the anchors to update the reference profiles
  updated_profiles <- updateProfilesFromAnchors(counts = counts, 
                                                neg = neg, 
                                                anchors = anchors, 
                                                reference_profiles = reference_profiles, 
                                                align_genes = TRUE, nb_size = 10, max_rescaling = 5)
  out <- list(updated_profiles = updated_profiles,
              anchors = anchors)
  return(out)
}



#' Use anchor cells to update reference profiles, simply by taking the mean
#' profile of the anchors.
#'
#' Uses anchor cells to estimate platform effects / scaling factors to be
#' applied to the genes/rows of the reference profile matrix. Then uses Bayesian
#' math to update the individual elements on X.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided
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
           anchors,
           reference_profiles,
           align_genes = TRUE,
           nb_size = 10,
           max_rescaling = 5) {
  use <- !is.na(anchors)
  updated_profiles <- Estep(counts = counts[use, ],
                            clust = anchors[use],
                            neg = neg[use])
  return(updated_profiles)
  }



#' Rescale_Reference is the function for reference profile re-scaling
#' Inputs@@
#' @param ioprofiles_Raw = raw reference profile
#' @param counts_input= cosmx counts data with genes on column and cells on row
#' @param negmeans_counts =mean expression for the negative probes for all cells
#' @param anchors_io_raw = a vector with name being cell ID and value being either cell type or NA
#' @param Gene_Out =user specified list of genes to be excluded for the cell typing, for example, 103 problematic probes identified during 6K panel development 
#' Output@@
#' @return An list consist of following things
#' 1: rescaled reference profile, which can be ready used as input reference for regular InStituType cell typing or you can do another round of anchor selection and profile updation with rescaled profile
#' 2: gene-level platform effect
#' 3: a vector of genes that are filtered out
#' General workflow@@
#' extract the anchor cells from input
#' Run poisson regression with anchor cells
#' Filter user defined genes(if any) and genes with negative betas 
#' Re-scale Profile with Beta estimates


updateReferenceProfiles<- function(ioprofiles_Raw=ioprofiles_Raw,counts_input=counts_input,negmeans_counts =negmeans_counts, 
                                   anchors_io_raw=anchors_io_raw, Gene_Out=Gene_Out, ncores=20 ){
  
  
  #### Step1: extract anchors based on input anchor information, and only keep anchor cells in CosMx data
  anchors_io_raw =anchors_io_raw[names(anchors_io_raw) %in% rownames(counts_input)]
  Common_Genes<-intersect(colnames(counts_input), rownames(ioprofiles_Raw))
  
  Anchor_index<-names(anchors_io_raw )[!is.na(anchors_io_raw)]
  CellType<-anchors_io_raw[Anchor_index]
  
  Count_Data<- counts_input[Anchor_index,Common_Genes]
  Reference_Data<- t(ioprofiles_Raw[Common_Genes,CellType])
  
  #get ride of genes that are purely 0s in reference or counts data
  Bad_Genes<-(colSums(Count_Data)==0 | colSums(Reference_Data)==0)
  Count_Data<- Count_Data[,!Bad_Genes]
  Reference_Data<- Reference_Data[,!Bad_Genes]
  
  Count_vector<-matrix(Count_Data, dimnames=list(t(outer(colnames(Count_Data), rownames(Count_Data), FUN=paste)), NULL))
  Count_vector= { dt = melt(data.table(Count_Data, keep.rownames = TRUE) , id.vars = c("rn"))} 
  colnames(Count_vector)<-c("Cell_ID","Gene","Count")
  
  Reference_vector= { dt = melt(data.table(Reference_Data, keep.rownames = TRUE) , id.vars = c("rn"))} 
  colnames(Reference_vector)<-c("Cell_Type","Gene","Count")
  
  #for cell level scaling factor, always make sure the orders are correct
  TotalCounts= rowSums(Count_Data)
  TotalExpression=rowSums(Reference_Data)
  
  #Determine background calculation:s =Vector of mean target counts per cell  
  s <- Matrix::rowMeans(counts_input)
  bgmod <- stats::lm(negmeans_counts ~ s - 1)
  bg <- bgmod$fitted
  names(bg) <- names(s)
  
  
  
  ### Step2: Run poisson regression with link being Identity to estimate gene-level platform effect based on selected anchor cells
  PlatformEffects_Data<-data.frame(Cell_ID=Count_vector$Cell_ID,
                                   Gene=Count_vector$Gene,
                                   Counts=Count_vector$Count,
                                   Cell_Type=anchors_io_raw[Count_vector$Cell_ID],
                                   Reference=Reference_vector$Count,
                                   Cell_SF=TotalCounts/TotalExpression,
                                   BG=bg[Count_vector$Cell_ID])
  
  
  Gene_list=unique(PlatformEffects_Data$Gene)
  
  estimatePlatformEffects<-function(gene_ID){
    GLM_Fit<- glm(Counts ~  Reference :  Cell_SF + offset(BG) -1 , family=poisson(link = identity),
                  data=PlatformEffects_Data[PlatformEffects_Data$Gene==gene_ID,],start=c(0))
    
    return(data.frame(Gene=gene_ID,
                      Beta=GLM_Fit$coefficients["Reference:Cell_SF"],
                      beta_SE=summary(GLM_Fit)$coefficients["Reference:Cell_SF","Std. Error"]))
  }
  
  PlatformEff<-mclapply(Gene_list, estimatePlatformEffects, mc.cores = ncores)
  PlatformEff<- as.data.frame(do.call(rbind, PlatformEff))
  rownames( PlatformEff)= PlatformEff$Gene
  
  
  ### Step3: Filtering genes with negative betas and pre-specified by the user
  Gene_Out_addon<-as.vector(PlatformEff$Gene[PlatformEff$Beta<0])
  Gene_Out<-unique(c(Gene_Out,Gene_Out_addon))
  Shared_Gene<-as.vector(PlatformEff$Gene[!(PlatformEff$Gene %in%Gene_Out)])
  
  ### Step4: Rescale the raw reference profile
  rownames(PlatformEff)=PlatformEff$Gene
  Rescaled_ioprofiles =diag(PlatformEff[Shared_Gene,]$Beta) %*% ioprofiles_Raw[Shared_Gene,]
  rownames(Rescaled_ioprofiles)=Shared_Gene
  
  return(list(Rescaled_ioprofiles=Rescaled_ioprofiles,
              PlatformEffect=PlatformEff[Shared_Gene,]$Beta,
              Gene_Out=Gene_Out))
  
}


