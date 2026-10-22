rank_targets <- function(integrated,annotation,weights=c(specificity=.35,surface=.30,safety=.25,concordance=.10),fdr=.05) {
 if(is.null(integrated) || !nrow(integrated)) stop('Integrate ADT and MS first')
 if(length(weights)!=4 || any(weights<0) || sum(weights)<=0) stop('Weights must be nonnegative and have positive sum')
 weights<-weights/sum(weights)
 d<-merge(integrated,annotation,by='feature')
 if(!nrow(d)) stop('No overlap with annotation feature IDs')
 # Significance gate: keep features significant (BH FDR) in at least one modality.
 if(all(c('adj.P.Val_ADT','adj.P.Val_MS') %in% names(d))) {
  sig<-(!is.na(d$adj.P.Val_ADT) & d$adj.P.Val_ADT<fdr) | (!is.na(d$adj.P.Val_MS) & d$adj.P.Val_MS<fdr)
  d<-d[sig,,drop=FALSE]
  if(!nrow(d)) stop(paste0('No features pass FDR < ',fdr,' in either modality'))
 }
 # Tumor specificity proxy is positive tumor-normal ADT and MS effect sizes.
 # This is not a substitute for external normal-tissue panel.
 d$specificity<-pmin(1,pmax(0,(d$logFC_ADT+d$logFC_MS)/4))
 d$safety<-1-d$normal_tissue_risk
 d$concordance<-as.numeric(d$logFC_ADT>0 & d$logFC_MS>0) # concordant UP-regulation only
 d$rank_score<-weights['specificity']*d$specificity+weights['surface']*d$surface_evidence+weights['safety']*d$safety+weights['concordance']*d$concordance
 d<-d[order(d$rank_score,decreasing=TRUE),,drop=FALSE]
 d$rank<-seq_len(nrow(d)); d
}
mod_prioritize_ui <- function(id){ns<-shiny::NS(id); bslib::card(bslib::card_header('Therapeutic target prioritization'),bslib::card_body(
 shiny::fileInput(ns('ann'),'Reference annotation CSV',accept='.csv'),
 shiny::helpText('Must include feature, surface_evidence (0–1), normal_tissue_risk (0–1). Values must come from documented sources; the tool does not fetch clinical expression references.'),
 shiny::actionButton(ns('demo'),'Load SIMULATED annotations'),
 shiny::sliderInput(ns('ws'),'Tumor-specificity weight',0,1,.35,step=.05),
 shiny::sliderInput(ns('we'),'Surface evidence weight',0,1,.30,step=.05),
 shiny::sliderInput(ns('wr'),'Normal-tissue safety weight',0,1,.25,step=.05),
 shiny::sliderInput(ns('wc'),'Concordance weight',0,1,.10,step=.05),
 shiny::numericInput(ns('fdr'),'FDR gate (significant in ≥1 modality)',.05,min=.001,max=1,step=.01),
 shiny::verbatimTextOutput(ns('status')),DT::DTOutput(ns('ranking'))))}
mod_prioritize_server <- function(id,integrated) shiny::moduleServer(id,function(input,output,session){
 annot<-shiny::reactiveVal(NULL); msg<-shiny::reactiveVal('Upload an annotation table to calculate rankings')
 shiny::observeEvent(input$ann,{tryCatch({annot(read_annotations(input$ann$datapath));msg('Loaded USER reference annotations')},error=function(e)msg(conditionMessage(e)))})
 shiny::observeEvent(input$demo,{tryCatch({d<-shiny::req(integrated());set.seed(44);annot(data.frame(feature=d$feature,surface_evidence=stats::runif(nrow(d)),normal_tissue_risk=stats::runif(nrow(d))));msg('SIMULATED annotation evidence — NOT validated therapeutic safety')},error=function(e)msg(conditionMessage(e)))})
 ranked<-shiny::reactive({shiny::req(annot(),integrated());rank_targets(integrated(),annot(),weights=c(specificity=input$ws,surface=input$we,safety=input$wr,concordance=input$wc),fdr=input$fdr)})
 output$status<-shiny::renderText(msg())
 output$ranking<-DT::renderDT({DT::datatable(ranked(),options=list(pageLength=10,scrollX=TRUE))})
 ranked
})
