plot_roc <- function(res) {
  d<-data.frame(FPR=1-res$roc$specificities,TPR=res$roc$sensitivities)
  ggplot2::ggplot(d,ggplot2::aes(FPR,TPR))+ggplot2::geom_line(color='#2563eb',linewidth=1.1)+ggplot2::geom_abline(slope=1,intercept=0,linetype=2,color='grey60')+ggplot2::coord_equal()+ggplot2::labs(title=sprintf('ROC · AUROC %.3f',res$auc),x='False positive rate',y='True positive rate')+ggplot2::theme_minimal(base_size=12)
}
plot_pr <- function(res) {
  d<-data.frame(Recall=res$pr$curve[,1],Precision=res$pr$curve[,2])
  ggplot2::ggplot(d,ggplot2::aes(Recall,Precision))+ggplot2::geom_line(color='#f97316',linewidth=1.1)+ggplot2::geom_hline(yintercept=res$prevalence,linetype=2,color='grey60')+ggplot2::coord_cartesian(ylim=c(0,1))+ggplot2::labs(title=sprintf('PR · AUPRC %.3f; prevalence %.2f',res$auprc,res$prevalence))+ggplot2::theme_minimal(base_size=12)
}
plot_calibration <- function(res,bins=10) {
  df<-res$df; br<-unique(stats::quantile(df$scoreA,probs=seq(0,1,length.out=bins+1)))
  if(length(br)<3) return(ggplot2::ggplot()+ggplot2::annotate('text',x=0,y=0,label='Not enough unique probabilities for calibration')+ggplot2::theme_void())
  grp<-cut(df$scoreA,breaks=br,include.lowest=TRUE)
  d<-stats::aggregate(cbind(pred=df$scoreA,truth=df$truth),by=list(bin=grp),FUN=mean)
  ggplot2::ggplot(d,ggplot2::aes(pred,truth))+ggplot2::geom_abline(linetype=2,color='grey60')+ggplot2::geom_point(color='#16a34a',size=2.5)+ggplot2::geom_line(color='#16a34a')+ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+ggplot2::labs(title='Calibration · quantile bins',x='Mean predicted probability',y='Observed positive fraction')+ggplot2::theme_minimal(base_size=12)
}
plot_volcano <- function(d,title='Differential expression') {
  d$neglogp <- -log10(pmax(d$adj.P.Val,.Machine$double.xmin))
  d$hit<-!is.na(d$adj.P.Val) & d$adj.P.Val<.05 & abs(d$logFC)>=1
  ggplot2::ggplot(d,ggplot2::aes(logFC,neglogp,color=hit))+ggplot2::geom_point(alpha=.7)+ggplot2::scale_color_manual(values=c('FALSE'='grey65','TRUE'='#f97316'))+ggplot2::labs(title=title,x='log2 fold-change',y='-log10 adjusted p-value',color='FDR<0.05 & |logFC|>=1')+ggplot2::theme_minimal(base_size=12)
}
plot_heatmap <- function(mat,top=25) {
  v<-apply(mat,1,stats::var,na.rm=TRUE); idx<-head(order(v,decreasing=TRUE,na.last=NA),top)
  x<-mat[idx,,drop=FALSE]; if(!length(x)) stop('No features to display')
  z<-t(scale(t(x))); z[!is.finite(z)]<-0
  d<-data.frame(feature=rep(rownames(z),times=ncol(z)),sample=rep(colnames(z),each=nrow(z)),z=as.vector(z))
  ggplot2::ggplot(d,ggplot2::aes(sample,feature,fill=z))+ggplot2::geom_tile()+ggplot2::scale_fill_gradient2(low='#2563eb',mid='white',high='#f97316')+ggplot2::labs(title='Top variable features · row z-scores',x=NULL,y=NULL)+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
}
plot_distribution <- function(mat,group='All') {
  d<-data.frame(value=as.vector(mat),sample=rep(colnames(mat),each=nrow(mat)))
  ggplot2::ggplot(d,ggplot2::aes(sample,value))+ggplot2::geom_boxplot(fill='#bfdbfe',outlier.alpha=.2)+ggplot2::labs(title=paste('Expression distribution ·',group),x=NULL,y='Normalized abundance')+ggplot2::theme_minimal()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1))
}
plot_concordance <- function(d) {
  ggplot2::ggplot(d,ggplot2::aes(logFC_ADT,logFC_MS))+ggplot2::geom_hline(yintercept=0,linetype=3)+ggplot2::geom_vline(xintercept=0,linetype=3)+ggplot2::geom_point(color='#2563eb',alpha=.65)+ggplot2::geom_smooth(method='lm',se=FALSE,color='#f97316')+ggplot2::labs(title='ADT vs proteomics differential concordance',x='ADT logFC',y='MS logFC')+ggplot2::theme_minimal()
}
