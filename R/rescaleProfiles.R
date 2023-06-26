
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





#' estimatePlatformEffects is the function for gene-level platform effect estimation and reference profile re-scaling
#' Inputs@@
#' @param profiles = raw reference profile
#' @param counts= cosmx counts data with genes on column and cells on row
#' @param neg =mean expression for the negative probes for all cells
#' @param bg= expected bg nosice for all cells
#' @param anchors = a vector with name being cell ID and value being either cell type or NA
#' @param blacklist =user specified list of genes to be excluded for the cell typing, for example, 103 problematic probes identified during 6K panel development 
#' Output@@
#' @return An list consist of following things
#' @return 1: rescaled reference profile, which can be ready used as input reference for regular InStituType cell typing or you can do another round of anchor selection and profile updation with rescaled profile
#' @return 2: gene-level platform effect
#' @return 3: a vector of genes that are filtered out
#' General workflow@@
#' extract the anchor cells from input
#' Run poisson regression with anchor cells
#' Filter user defined genes(if any) and genes with negative betas 
#' Re-scale Profile with Beta estimates


estimatePlatformEffects<- function(profiles=profiles,counts =counts,neg =neg, bg=bg, 
                                   anchors=anchors, blacklist=blacklist, ncores=20 ){
  
  
  #### Step1: extract anchors based on input anchor information, and only keep anchor cells in CosMx data
  anchors =anchors[names(anchors) %in% rownames(counts)]
  Common_Genes<-intersect(colnames(counts), rownames(profiles))
  
  Anchor_index<-names(anchors )[!is.na(anchors)]
  CellType<-anchors[Anchor_index]
  
  Count_Data<- counts[Anchor_index,Common_Genes]
  Reference_Data<- t(profiles[Common_Genes,CellType])
  
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
  
  ### Step2: Run poisson regression with link being Identity to estimate gene-level platform effect based on selected anchor cells
  PlatformEffects_Data<-data.frame(Cell_ID=Count_vector$Cell_ID,
                                   Gene=Count_vector$Gene,
                                   Counts=Count_vector$Count,
                                   Cell_Type=anchors[Count_vector$Cell_ID],
                                   Reference=Reference_vector$Count,
                                   Cell_SF=TotalCounts/TotalExpression,
                                   BG=bg[Count_vector$Cell_ID])
  
  
  Gene_list=unique(PlatformEffects_Data$Gene)
  
  PlatformEstimator<-function(gene_ID){
    GLM_Fit<- glm(Counts ~  Reference :  Cell_SF + offset(BG) -1 , family=poisson(link = identity),
                  data=PlatformEffects_Data[PlatformEffects_Data$Gene==gene_ID,],start=c(0))
    
    return(data.frame(Gene=gene_ID,
                      Beta=GLM_Fit$coefficients["Reference:Cell_SF"],
                      beta_SE=summary(GLM_Fit)$coefficients["Reference:Cell_SF","Std. Error"]))
  }
  
  PlatformEff<-mclapply(Gene_list, PlatformEstimator, mc.cores = ncores)
  PlatformEff<- as.data.frame(do.call(rbind, PlatformEff))
  rownames( PlatformEff)= PlatformEff$Gene
  
  
  ### Step3: Filtering genes with negative betas and pre-specified by the user
  blacklist_addon<-as.vector(PlatformEff$Gene[PlatformEff$Beta<0])
  blacklist<-unique(c(blacklist,blacklist_addon))
  Shared_Gene<-as.vector(PlatformEff$Gene[!(PlatformEff$Gene %in%blacklist)])
  
  ### Step4: Rescale the raw reference profile
  rownames(PlatformEff) <- PlatformEff$Gene
  Rescaled_ioprofiles <- diag(PlatformEff[Shared_Gene,]$Beta) %*% profiles[Shared_Gene,]
  rownames(Rescaled_ioprofiles)<- Shared_Gene
  
  PlatformEffect<-PlatformEff[Shared_Gene,]$Beta
  names(PlatformEffect)<- Shared_Gene
  
  return(list(Rescaled_ioprofiles=Rescaled_ioprofiles,
              PlatformEffect=PlatformEffect,
              blacklist=blacklist))
  
}





