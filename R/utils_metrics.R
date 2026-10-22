check_predictions <- function(df) {
  if(!is.data.frame(df) || !all(c('truth','scoreA') %in% names(df))) stop('Predictions require truth and scoreA columns')
  if(anyNA(df$truth) || !all(as.character(df$truth) %in% c('0','1'))) stop('truth must contain only 0 and 1')
  df$truth<-as.integer(as.character(df$truth))
  if(length(unique(df$truth))!=2L) stop('Both outcome classes required')
  for(nm in intersect(c('scoreA','scoreB'),names(df))) {
    df[[nm]]<-suppressWarnings(as.numeric(as.character(df[[nm]])))
    if(any(!is.finite(df[[nm]])) || any(df[[nm]]<0 | df[[nm]]>1)) stop(paste(nm,'must contain finite probabilities from 0 to 1'))
  }
  df
}
metric_auc <- function(y,p) {
  r<-pROC::roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)
  as.numeric(pROC::auc(r))
}
metric_pr <- function(y,p) {
  as.numeric(PRROC::pr.curve(scores.class0=p[y==1],scores.class1=p[y==0],curve=FALSE)$auc.integral)
}
metric_brier <- function(y,p) mean((y-p)^2)
confusion_at <- function(y,p,threshold=.5) {
  pred<-as.integer(p>=threshold)
  tp<-sum(y==1 & pred==1); tn<-sum(y==0 & pred==0)
  fp<-sum(y==0 & pred==1); fn<-sum(y==1 & pred==0)
  data.frame(TP=tp,FP=fp,TN=tn,FN=fn,Sensitivity=if((tp+fn)>0) tp/(tp+fn) else NA_real_,Specificity=if((tn+fp)>0) tn/(tn+fp) else NA_real_,Precision=if((tp+fp)>0) tp/(tp+fp) else NA_real_)
}
bootstrap_metrics <- function(y,p,B=200,seed=42) {
  B<-as.integer(B); if(B<20 || B>5000) stop('Bootstrap iterations must be 20–5000')
  set.seed(seed); pos<-which(y==1); neg<-which(y==0)
  vals<-replicate(B, {
    ids<-c(sample(pos,length(pos),replace=TRUE),sample(neg,length(neg),replace=TRUE))
    c(AUROC=metric_auc(y[ids],p[ids]),AUPRC=metric_pr(y[ids],p[ids]),Brier=metric_brier(y[ids],p[ids]))
  })
  t(apply(vals,1,stats::quantile,probs=c(.025,.975),na.rm=TRUE))
}
evaluate_scores <- function(raw, threshold=.5, B=200) {
  df<-check_predictions(raw)
  y<-df$truth;p<-df$scoreA
  rocA<-pROC::roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)
  prA<-PRROC::pr.curve(scores.class0=p[y==1],scores.class1=p[y==0],curve=TRUE)
  comp<-NULL
  if('scoreB' %in% names(df)) {
    rocB<-pROC::roc(y,df$scoreB,levels=c(0,1),direction='<',quiet=TRUE)
    test<-pROC::roc.test(rocA,rocB,method='delong',paired=TRUE)
    comp<-list(aucB=as.numeric(pROC::auc(rocB)),delta=as.numeric(pROC::auc(rocA)-pROC::auc(rocB)),pvalue=as.numeric(test$p.value),rocB=rocB)
  }
  list(df=df,auc=as.numeric(pROC::auc(rocA)),auprc=prA$auc.integral,brier=metric_brier(y,p),roc=rocA,pr=prA,prevalence=mean(y),ci=bootstrap_metrics(y,p,B=B),confusion=confusion_at(y,p,threshold),threshold=threshold,comparison=comp,iterations=B)
}
