# Explicit held-out and leave-one-cohort-out predictions. Fit only on training rows.
run_external_model <- function(df, feature_cols, validation='holdout', frac=.75, seed=2026) {
 if(!all(c('truth','cohort') %in% names(df))) stop('ML data needs truth and cohort columns')
 if(!all(as.character(df$truth) %in% c('0','1'))) stop('Binary truth 0/1 required')
 df$truth<-as.integer(as.character(df$truth)); df$cohort<-as.character(df$cohort)
 if(length(unique(df$truth))!=2) stop('Both classes required')
 if(length(feature_cols)<1 || any(!feature_cols %in% names(df))) stop('Select numeric features')
 if(anyNA(df[,c('truth','cohort',feature_cols),drop=FALSE])) stop('Missing values: preprocess upstream using training-only statistics')
 if(any(!vapply(df[feature_cols],is.numeric,logical(1)))) stop('All selected predictors must be numeric')
 if(validation=='LOCO') {
  folds<-lapply(unique(df$cohort),function(co) which(df$cohort==co)); names(folds)<-unique(df$cohort)
  if(length(folds)<2) stop('LOCO needs at least two cohorts')
 } else {
  set.seed(seed)
  # Stratified split for reproducibility
  train<-unlist(lapply(split(seq_len(nrow(df)),df$truth),function(ids) sample(ids,max(1,min(length(ids)-1,floor(length(ids)*frac))))))
  folds<-list(Heldout=setdiff(seq_len(nrow(df)),train))
 }
 predictions<-list()
 for(fold in names(folds)) {
  test<-folds[[fold]]; train<-setdiff(seq_len(nrow(df)),test)
  if(length(unique(df$truth[train]))<2 || length(unique(df$truth[test]))<2) stop(paste('Fold',fold,'does not contain both classes in train and test'))
  xtrain<-df[train,feature_cols,drop=FALSE]; xtest<-df[test,feature_cols,drop=FALSE]
  mu<-vapply(xtrain,mean,numeric(1)); sdv<-vapply(xtrain,stats::sd,numeric(1));sdv[sdv==0]<-1
  tx<-as.data.frame(sweep(sweep(as.matrix(xtrain),2,mu,'-'),2,sdv,'/'))
  vx<-as.data.frame(sweep(sweep(as.matrix(xtest),2,mu,'-'),2,sdv,'/'))
  names(tx)<-names(vx)<-make.names(feature_cols,unique=TRUE)
  tx$truth<-df$truth[train]
  fit<-suppressWarnings(stats::glm(truth~.,family=stats::binomial(),data=tx))
  p<-suppressWarnings(stats::predict(fit,newdata=vx,type='response'))
  if(any(!is.finite(p))) stop(paste('Model produced non-finite probabilities in',fold,'— possible complete separation'))
  # Reference prevalence model fitted ONLY on training data
  baseline<-rep(mean(df$truth[train]),length(test))
  predictions[[fold]]<-data.frame(sample_id=if('sample_id' %in% names(df)) df$sample_id[test] else as.character(test),cohort=df$cohort[test],fold=fold,truth=df$truth[test],scoreA=p,scoreB=baseline)
 }
 do.call(rbind,predictions)
}
mod_ml_ui <- function(id){ns<-shiny::NS(id);bslib::card(bslib::card_header('External validation · reproducible logistic regression'),bslib::card_body(
 shiny::helpText('CSV columns truth (0/1), cohort, optional sample_id, and numerical feature columns. This module trains GLM on training partitions only. Baseline model B predicts training prevalence.'),
 shiny::fileInput(ns('file'),'Model feature table CSV',accept='.csv'),
 shiny::selectInput(ns('strategy'),'Validation strategy',c('Stratified hold-out'='holdout','Leave-one-cohort-out'='LOCO')),
 shiny::sliderInput(ns('frac'),'Hold-out training fraction',.5,.9,.75,step=.05),
 shiny::uiOutput(ns('features_ui')),shiny::actionButton(ns('fit'),'Fit and validate',class='btn-primary'),
 shiny::verbatimTextOutput(ns('status')),DT::DTOutput(ns('pred_table'))))}
mod_ml_server <- function(id) shiny::moduleServer(id,function(input,output,session){
 uploaded<-shiny::reactiveVal(NULL); preds<-shiny::reactiveVal(NULL);msg<-shiny::reactiveVal('Upload training features')
 shiny::observeEvent(input$file,{tryCatch({uploaded(utils::read.csv(input$file$datapath,check.names=FALSE));preds(NULL);msg('Loaded feature table. Choose predictors and validate.')},error=function(e)msg(conditionMessage(e)))})
 output$features_ui<-shiny::renderUI({d<-shiny::req(uploaded());cols<-setdiff(names(d)[vapply(d,is.numeric,logical(1))],c('truth','cohort','sample_id','fold'));shiny::selectInput(session$ns('features'),'Predictors',choices=cols,selected=head(cols,5),multiple=TRUE)})
 shiny::observeEvent(input$fit,{tryCatch({shiny::req(uploaded(),input$features);preds(run_external_model(uploaded(),input$features,input$strategy,input$frac));msg(paste('OUT-OF-SAMPLE predictions:',nrow(preds()),'samples across',length(unique(preds()$fold)),'validation fold(s). Select Evaluation to inspect metrics.'))},error=function(e){preds(NULL);msg(paste('ERROR:',conditionMessage(e)))})})
 output$status<-shiny::renderText(msg())
 output$pred_table<-DT::renderDT({DT::datatable(shiny::req(preds()),options=list(pageLength=8,scrollX=TRUE))})
 shiny::reactive(preds())
})
