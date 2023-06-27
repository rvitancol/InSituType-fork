#' Prepare bg data for other functions 
#' 
#' Process neg data or bg to get background for each cell
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @return A named vector for the estimated background of each cell

estimateBackground <- function(counts, neg, bg = NULL){
  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg) && is.null(neg)) {
    stop("Must provide either bg or neg")
  }
  
  if (is.null(bg)) {
    ## get neg in condition 
    if (is.null(names(neg))) {
      names(neg) <- rownames(counts)
    }
    if (length(neg) != nrow(counts)) {
      stop("length of neg should equal nrows of counts.")
    }
    
    s <- Matrix::rowMeans(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  return(bg)
  
}


#' align genes in counts to profiles for other functions 
#' 
#' Process counts to have genes shared with profiles
#' @param counts Counts matrix, cells * genes.
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @return a cells * genes count matrix for shared genes only
alignGenes <- function(counts, profiles){
  sharedgenes <- intersect(rownames(profiles), colnames(counts))
  lostgenes <- setdiff(colnames(counts), rownames(profiles))
  
  # subset:
  counts <- counts[, sharedgenes]
  
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

  return(counts)
}


#' Get number of cores for parallelized operations
#'
#' @return number of cores to use for mclapply
#' @export
numCores <- function() {
  num_cores <- 1
  if (.Platform$OS.type == "unix") {
    if (is.null(getOption("mc.cores"))) {
      num_cores <- parallel::detectCores() - 2
    } else {
      num_cores <- getOption("mc.cores")
    }
    
  }
  return(num_cores)
}