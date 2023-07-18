#' Generate the mean reference profile and its SD reference profile from an annotation file
#' This function is only for protein data set with known anchor cells and their cell types
#'
#' @param exp.mat a matrix of raw protein expression data. cells are in rows and proteins are in columns
#' @param anno a data frame or matrix of cell types for anchor cells or manually annotated cell typing information for some cells. Should include cell_ID and celltype at least. 
#' 
#' @return A list, with the following elements:
#' \enumerate{
#' \item mean.ref.profile: a matrix of cluster-specific expression profiles. proteins * cell types
#' \item SDs.ref.profile: a matrix of standard deviation profiles of pre-defined clusters. proteins * cell types
#' \item anchors: a vector giving "anchor" cell types. Vector elements will be mainly NA's (for non-anchored cells)
#' }

gen_profiles_protein_annotation <- function(exp.mat, anno) {

  library(dplyr)
  library(tibble)
  
  anno_ref_mat <- merge(exp.mat %>% as.data.frame() %>% rownames_to_column(var="cell_ID"), anno %>% dplyr::select(c(cell_ID, cellType)), by="cell_ID") %>% column_to_rownames(var="cell_ID")
  
  mean.ref.profile <- anno_ref_mat %>% group_by(cellType) %>% summarise_all(mean) %>% column_to_rownames(var="cellType") %>% t()
  SDs.ref.profile <- anno_ref_mat %>% group_by(cellType) %>% summarise_all(sd) %>% column_to_rownames(var="cellType") %>% t()
  
  ## Set NAs for non-anchor cells' cell types
  anchors <- rbind(anno %>% dplyr::select(c(cell_ID, cellType)), data.frame(cell_ID = setdiff(rownames(exp.mat), rownames(anno_ref_mat)), cellType=NA)) 
  rownames(anchors) <- anchors$cell_ID
  anchors$cell_ID <- NULL
  anchors <- anchors %>% t()
  anchors <- anchors[1,]
  
  out <- list(mean.ref.profile=mean.ref.profile,
              SDs.ref.profile=SDs.ref.profile,
              anchors=anchors[rownames(exp.mat)])
  return(out)
}


#' Generate the mean reference profile and its SD reference profile based on the data itself
#' This function is based on signature matrix included in CELESTA package 
#' First, we rebuild a nested cell typing lists based on the 2-D signature matrix
#' Second, we identify anchor cells ranked by their expression level for each cell type's protein marker
#' Third, we estimate averaged expression level and SDs for proteins and cell types using the anchors
#'
#' @param exp.mat a matrix of raw protein expression data. cells are in rows and proteins are in columns
#' @param sig_mat a signature matrix of cell types. cell types x protein markers 
#' @param cutoff a cutoff of quantile. e.g) cutoff=0.9 means that top 90 percentiles of cells are called anchors for the protein expression
#' @param min.num.cells a minimum number of cells each cell type to estimate its mean or SDs. default value is 30.
#' 
#' @return A list, with the following elements:
#' \enumerate{
#' \item mean.ref.profile: a matrix of cluster-specific expression profiles. proteins x cell types
#' \item SDs.ref.profile: a matrix of standard deviation profiles of pre-defined clusters. proteins x cell types
#' \item anchors: a vector giving "anchor" cell types. Vector elements will be mainly NA's (for non-anchored cells)
#' }

gen_profiles_protein_expression_old <- function(exp.mat, sig_mat=NULL, cutoff=0.9, min.num.cells=30){
  library(dplyr)
  library(data.table)
  if(is.null(sig_mat)){
    sig_mat = fread(paste0(system.file("extdata", package="smiProtein"), "/default_signature_matrix.csv"))
  }
  
  ## Split Lineage levels into columns
  sig_mat[is.na(sig_mat)] <- 0
  sig_mat$level1 <- lapply(strsplit(sig_mat$Lineage_level, "_"), function(x){x[1]}) %>% unlist()
  sig_mat$level2 <- lapply(strsplit(sig_mat$Lineage_level, "_"), function(x){x[2]}) %>% unlist()
  sig_mat$level3 <- lapply(strsplit(sig_mat$Lineage_level, "_"), function(x){x[3]}) %>% unlist()
  
  ## data.frame including cellType, and its protein marker and its upper cell type where this cell type belong to its upper cell type (level 1 is parental cell types)
  markerProtein_celltype_level1 <- data.frame(celltype = sig_mat[sig_mat$level1==1,]$celltype, 
                                              marker_protein=apply(sig_mat[sig_mat$level1==1,], 1, function(x){colnames(sig_mat[sig_mat$level1==1,])[which(x==1)[1]]}),
                                              upper_celltype = "Parent")
  
  ## level 2 is subset of a parental cell type (e.g. t/b-cell, NK, dendritic cell are subgroups of Immune cell type (parent celltype))
  markerProtein_celltype_level2 <- data.frame(celltype = sig_mat[sig_mat$level1==2,]$celltype, 
                                              marker_protein=apply(sig_mat[sig_mat$level1==2,], 1, function(x){colnames(sig_mat[sig_mat$level1==2,])[which(x==1)[1]]}),
                                              upper_celltype = sig_mat$celltype[which(sig_mat$level3==unique(sig_mat[sig_mat$level1==2,]$level2))])
  
  ## level 3 is subset of a level 1 cell type (e.g. CD4-t-cell, CD8-t-cell)
  markerProtein_celltype_level3 <- data.frame(celltype = sig_mat[sig_mat$level1==3,]$celltype, 
                                              marker_protein=apply(sig_mat[sig_mat$level1==3,], 1, function(x){colnames(sig_mat[sig_mat$level1==3,])[which(x==1)[1]]}),
                                              upper_celltype = sig_mat$celltype[which(sig_mat$level3==unique(sig_mat[sig_mat$level1==3,]$level2))])
  
  
  dat_mat_level1 <- lapply(markerProtein_celltype_level1$marker_protein, function(x){
    if(max(exp.mat)<=1){
      cutoff <- 0.9
    }else{
      cutoff <- quantile(exp.mat[, x], prob=0.9)
    }
    
    rownames(exp.mat)[which(exp.mat[, x] > cutoff)]
  })
  names(dat_mat_level1) <- markerProtein_celltype_level1$marker_protein
  
  
  dat_mat_level2 <- vector("list", nrow(markerProtein_celltype_level2))
  names(dat_mat_level2) <- markerProtein_celltype_level2$marker_protein
  for(i in markerProtein_celltype_level2$marker_protein){
    tempMar <- markerProtein_celltype_level1 %>% filter(celltype==markerProtein_celltype_level2$upper_celltype[1])
    tempD <- exp.mat[rownames(exp.mat) %in% dat_mat_level1[[tempMar$marker_protein]], ]
    
    if(max(exp.mat)<=1){
      cutoff <- 0.9
    }else{
      cutoff <- quantile(tempD[, i], prob=0.9)
    }
    
    tempID <- rownames(tempD)[which(tempD[, i] > cutoff)]
    
    dat_mat_level1[[tempMar$marker_protein]] <- setdiff(dat_mat_level1[[tempMar$marker_protein]], tempID)
    dat_mat_level2[[i]] <- tempID
  }
  
  
  dat_mat_level3 <- vector("list", nrow(markerProtein_celltype_level3))
  names(dat_mat_level3) <- markerProtein_celltype_level3$marker_protein
  for(i in markerProtein_celltype_level3$marker_protein){
    tempMar <- markerProtein_celltype_level2 %>% filter(celltype==markerProtein_celltype_level3$upper_celltype[1])
    tempD <- exp.mat[rownames(exp.mat) %in% dat_mat_level2[[tempMar$marker_protein]], ]
    
    if(max(exp.mat)<=1){
      cutoff <- 0.9
    }else{
      cutoff <- quantile(tempD[, i], prob=0.9)
    }
    
    tempID <- rownames(tempD)[which(tempD[, i] > cutoff)]
    
    dat_mat_level2[[tempMar$marker_protein]] <- setdiff(dat_mat_level2[[tempMar$marker_protein]], tempID)
    
    dat_mat_level3[[i]] <- tempID
  }
  
  markerProtein_celltype_all <- rbind(markerProtein_celltype_level1, markerProtein_celltype_level2, markerProtein_celltype_level3)
  marker_id_cell_type <- c(dat_mat_level1, dat_mat_level2, dat_mat_level3)
  names(marker_id_cell_type) <- markerProtein_celltype_all$celltype
  
  
  marker_id_cell_type_insitu <- lapply(1:length(marker_id_cell_type), function(x){data.frame(cell_ID=marker_id_cell_type[[x]],
                                                                                             celltype=rep(names(marker_id_cell_type[x]), length(marker_id_cell_type[[x]])))})
  anchors <- do.call("rbind", marker_id_cell_type_insitu) %>% as.data.frame()
  anchors_duplicate <- anchors[which(duplicated(anchors$cell_ID)==TRUE),]$cell_ID
  
  marker_id_cell_type_unique <- lapply(marker_id_cell_type, 
                                       function(x) {
                                         tempV <- setdiff(x, anchors_duplicate)
                                         if(length(tempV) > 20){
                                           tempV <- tempV
                                         }else{
                                           tempV <- NULL
                                         }
                                         return(tempV)})

  marker_id_cell_type_unique <- Filter(Negate(is.null), marker_id_cell_type_unique)
  
  anchors <- anchors[which(duplicated(anchors$cell_ID)==FALSE),] 
  anchors <- anchors %>% filter(celltype %in% names(marker_id_cell_type_unique))
  
  anchors <- rbind(anchors, data.frame(cell_ID = setdiff(rownames(exp.mat), anchors$cell_ID), celltype=NA))
  rownames(anchors) <- anchors$cell_ID
  anchors$cell_ID <- NULL
  anchors <- t(anchors)[1,]
  
  ############################ Estimate averaged protein expression each cell type with its anchor cells ######################################
  protein_exp_means_list <- lapply(marker_id_cell_type_unique, function(x){
    
      mean.exp <- exp.mat[rownames(exp.mat) %in% x, ] %>% colMeans()

  })

  mean.ref.profile <- do.call("rbind", protein_exp_means_list) %>% t() %>% as.data.frame()
  
  protein_exp_SDs_list <- lapply(marker_id_cell_type_unique, function(x){
    apply(exp.mat[rownames(exp.mat) %in% x, ], 2, sd )
  })
  names(protein_exp_SDs_list) <- names(marker_id_cell_type_unique)
  
  SDs.ref.profile <- do.call("rbind", protein_exp_SDs_list) %>% t() %>% as.data.frame()
  
  out <- list(mean.ref.profile=mean.ref.profile, SDs.ref.profile=SDs.ref.profile, anchors=anchors[rownames(exp.mat)])
  return(out)
}
