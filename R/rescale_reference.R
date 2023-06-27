library(parallel)
library(InSituType)

#' rescale_reference is the function for reference profile re-scaling.
#' @description
#' \enumerate{\item anchor selection with the original reference profile
#' \item with the anchor cells, run poisson regression for RNA data and gaussian regression for protein data
#' \item filter user defined genes or protein (if any) and genes with negative betas
#' \item re-scale profile with beta estimates }
#' 
#' @param counts Counts matrix (or dgCMatrix), cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param reference_profiles Matrix of mean expression profiles of pre-defined
#'   clusters, e.g. from previous scRNA-seq. These profiles will not be updated
#'   by the EM algorithm. Columns must all be included in the init_clust
#'   variable.
#' @param reference_sds Matrix of standard deviation profiles of pre-defined
#'   clusters. These SD profiles also will not be updated by the EM algorithm. 
#'   Columns must all be included in the init_clust variable. This parameter should
#'   be defined if assay_type is protein. Default is NULL. 
#' @param exclude_genes A vector of genes to be excluded for the cell typing, for example, 103 problematic probes identified during 6K panel development 
#' @param assay_type Assay type of RNA, protein 
#' @export
#' @return A list, with the following elements: \enumerate{ \item reference_profile: a
#'   matrix of rescaled mean reference profile \item reference_sds: a
#'   matrix of rescaled SD reference profile 

rescale_reference <- function(counts, neg, reference_profiles, reference_sds, exclude_genes, ncores=20 ){
  
  #Determine background calculation:s =Vector of mean target counts per cell  
  s <- Matrix::rowMeans(counts)
  sum(names(neg )==names(s))==length(s)
  bgmod <- stats::lm(neg ~ s - 1)
  bg <- bgmod$fitted
  names(bg) <- names(s)
  
  
  
  #### Step1: find anchors based on raw reference profile; 
  #### there is no need to use less stringent cutoff for anchor selections as those selected anchors are just for gene-level platform effect estimation
  astats_io_raw <- get_anchor_stats(counts = counts,
                                    neg = neg,
                                    profiles = reference_profiles,
                                    min_cosine = 0.05)
  
  anchors_io_raw <- choose_anchors_from_stats(counts=counts,
                                              neg=neg,
                                              bg = bg,
                                              anchorstats=astats_io_raw ,
                                              n_cells = 400,
                                              min_cosine=0.15,
                                              min_scaled_llr = 0.03,
                                              insufficient_anchors_thresh=20)
  
  ### Step2: Run poisson regression with link being Indentity to estimate gene-level platform effect based on selected anchor cells
  Anchor_index<-names(anchors_io_raw )[!is.na(anchors_io_raw)]
  CellType<-anchors_io_raw[Anchor_index]
  
  Count_Data<- counts[Anchor_index,Common_Genes]
  Reference_Data<- t(reference_profiles[Common_Genes,CellType])
  
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

  
  GLMNB_Data<-data.frame(Cell_ID=Count_vector$Cell_ID,
                         Gene=Count_vector$Gene,
                         Counts=Count_vector$Count,
                         Cell_Type=anchors_io_raw[Count_vector$Cell_ID],
                         Reference=Reference_vector$Count,
                         Cell_SF=TotalCounts/TotalExpression,
                         BG=bg[Count_vector$Cell_ID],
                         Per_cell_bg=per.cell.bg[Count_vector$Cell_ID])
  
  
  Gene_list=unique(GLMNB_Data$Gene)
  
  GlM_Run<-function(gene_ID){
    GLM_Fit<- glm(Counts ~  Reference :  Cell_SF + offset(BG) -1 , family=poisson(link = identity),data=GLMNB_Data[GLMNB_Data$Gene==gene_ID,],start=c(0))
    return(data.frame(Gene=gene_ID,
                      Beta=GLM_Fit$coefficients["Reference:Cell_SF"],
                      beta_SE=summary(GLM_Fit)$coefficients["Reference:Cell_SF","Std. Error"]))
  }
  
  GLM_Results<-mclapply(Gene_list, GlM_Run, mc.cores = ncores)
  GLM_Results<- as.data.frame(do.call(rbind, GLM_Results))
  rownames( GLM_Results)= GLM_Results$Gene
  
  
  ### Step3: Filtering genes with negative betas and pre-specified by the user
  Gene_Out_addon<-as.vector(GLM_Results$Gene[GLM_Results$Beta<0])
  exclude_genes<-unique(c(exclude_genes,Gene_Out_addon))
  Shared_Gene<-as.vector(GLM_Results$Gene[!(GLM_Results$Gene %in%exclude_genes)])
  
  ### Step4: Rescale the raw reference profile
  rownames(GLM_Results)=GLM_Results$Gene
  Rescaled_ioprofiles =diag(GLM_Results[Shared_Gene,]$Beta) %*% reference_profiles[Shared_Gene,]
  rownames(Rescaled_ioprofiles)=Shared_Gene
  
  return(Rescaled_ioprofiles)
  
}



####################Toy Example ##################
#SeuratObj is just any cosmx seurat object
counts <-t(as.matrix(SeuratObj@assays$RNA@counts))
neg <- rowMeans(t(SeuratObj@assays$negprobes@counts))


dim(counts)
system.time( 
ioprofiles_Rescaled<- Rescale_Reference(reference_profiles=reference_profiles, counts=counts, neg=neg,
                              exclude_genes=exclude_genes, ncores=20 ) 
)



  