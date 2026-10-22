# Input contracts: all matrices are feature-by-sample, metadata columns sample_id, group.
validate_metadata <- function(meta, samples) {
  stopifnot(is.data.frame(meta))
  if (!all(c('sample_id','group') %in% names(meta))) stop('metadata needs sample_id and group columns')
  meta$sample_id <- as.character(meta$sample_id)
  meta$group <- as.character(meta$group)
  if (anyDuplicated(meta$sample_id) || anyNA(meta$group)) stop('Duplicate sample_id or missing groups')
  if (!all(samples %in% meta$sample_id)) stop('Metadata is missing matrix sample IDs')
  meta[match(samples,meta$sample_id),,drop=FALSE]
}
read_feature_matrix <- function(path, filename=path) {
  ext <- tolower(tools::file_ext(filename))
  if (ext=='rds') {
    obj <- readRDS(path)
    if (inherits(obj,'Seurat')) stop('Export ADT assay counts as a feature x sample matrix .rds before upload')
    m <- as.matrix(obj)
  } else {
    dat <- utils::read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
    if(ncol(dat)<3) stop('Matrix CSV needs a feature column plus at least two samples')
    ids <- as.character(dat[[1]]); if(anyDuplicated(ids) || any(!nzchar(ids))) stop('Feature IDs must be unique and nonempty')
    m <- as.matrix(dat[,-1,drop=FALSE]); rownames(m)<-ids
  }
  suppressWarnings(storage.mode(m)<-'double')
  if(is.null(rownames(m)) || is.null(colnames(m)) || anyDuplicated(colnames(m))) stop('Matrix must have feature row names and unique sample column names')
  if(any(!is.finite(m[!is.na(m)]))) stop('Matrix contains infinite entries')
  if(nrow(m)<2 || ncol(m)<4) stop('At least 2 features and 4 samples required')
  m
}
normalize_adt <- function(mat, method='CLR') {
  if(anyNA(mat) || any(mat<0) || any(abs(mat-round(mat))>1e-7)) stop('ADT counts must be nonnegative whole numbers without missing values')
  if(method=='CLR') {
    # Within-sample CLR on log2(count+1) so logFC is on the same log2 scale as MS.
    logmat<-log2(mat+1); sweep(logmat,2,colMeans(logmat),'-')
  } else if(method=='logCPM') {
    totals<-colSums(mat); if(any(totals==0)) stop('Zero-library samples cannot be normalized')
    log2(t(t(mat)/totals)*1e6+1)
  } else stop('Unsupported ADT normalization')
}
normalize_ms <- function(mat, logged=FALSE) {
  if(any(mat[!is.na(mat)]<0) && !logged) stop('Raw MS intensities must be nonnegative')
  if(!logged) mat<-log2(mat+1)
  med<-apply(mat,2,stats::median,na.rm=TRUE)
  sweep(mat,2,med,'-') + stats::median(med,na.rm=TRUE)
}
fit_differential <- function(mat, meta, case='Tumor', control='Normal', paired=FALSE) {
  meta<-validate_metadata(meta,colnames(mat))
  keep<-meta$group %in% c(case,control)
  x<-mat[,keep,drop=FALSE]; md<-meta[keep,,drop=FALSE]
  if(sum(md$group==case)<2 || sum(md$group==control)<2) stop('Need at least 2 case and 2 control samples')
  if(paired) stop('Paired design requires explicit patient/subject blocking; not available in this workflow')
  if(!requireNamespace('limma',quietly=TRUE)) stop('Install Bioconductor package limma for differential expression')
  design<-stats::model.matrix(~ factor(md$group,levels=c(control,case)))
  colnames(design)<-c('Intercept','CaseVsControl')
  if(anyNA(x)) { # limma lmFit handles missing measurements, but cannot fit all-missing rows
    valid<-rowSums(is.finite(x[,md$group==case,drop=FALSE]))>=2 & rowSums(is.finite(x[,md$group==control,drop=FALSE]))>=2
  } else valid<-rep(TRUE,nrow(x))
  if(!any(valid)) stop('No features with sufficient measurements in both groups')
  fit<-limma::lmFit(x[valid,,drop=FALSE],design)
  fit<-limma::eBayes(fit)
  tt<-limma::topTable(fit,coef='CaseVsControl',number=Inf,sort.by='none')
  tt$feature<-rownames(tt)
  ans<-data.frame(feature=tt$feature,logFC=tt$logFC,P.Value=tt$P.Value,adj.P.Val=tt$adj.P.Val,t=tt$t,stringsAsFactors=FALSE)
  ans[order(ans$adj.P.Val),,drop=FALSE]
}
read_annotations <- function(path) {
  x<-utils::read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  needed<-c('feature','surface_evidence','normal_tissue_risk')
  if(!all(needed %in% names(x))) stop('Annotation CSV needs: feature, surface_evidence, normal_tissue_risk')
  if(anyDuplicated(x$feature)) stop('Annotation IDs must be unique')
  for(nm in c('surface_evidence','normal_tissue_risk')) {
    x[[nm]]<-suppressWarnings(as.numeric(x[[nm]]))
    if(anyNA(x[[nm]]) || any(x[[nm]]<0 | x[[nm]]>1)) stop(paste(nm,'must be between 0 and 1'))
  }
  x
}
