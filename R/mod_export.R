mod_export_ui <- function(id){ns<-shiny::NS(id);bslib::card(bslib::card_header('Publication-ready export & provenance'),bslib::card_body(
 shiny::helpText('PDF contains only analyses actually run. Generated PDF is vector-based using the R PDF device.'),
 shiny::downloadButton(ns('pdf'),'All available figures · vector PDF',class='btn-warning w-100 mb-2'),
 shiny::downloadButton(ns('bundle'),'Analysis objects · RDS',class='btn-outline-primary w-100 mb-2'),
 shiny::downloadButton(ns('targets'),'Target ranking · CSV',class='btn-outline-success w-100 mb-2'),
 shiny::downloadButton(ns('methods'),'Methods and parameters · TXT',class='btn-outline-secondary w-100 mb-2'),
 shiny::verbatimTextOutput(ns('provenance'))))}
mod_export_server <- function(id,adt,ms,integrated,ranking,evaluation,predictions) shiny::moduleServer(id,function(input,output,session){
 safe<-function(rx) tryCatch(shiny::isolate(rx()),error=function(e) NULL)
 payload<-shiny::reactive(list(adt=safe(adt),ms=safe(ms),integration=safe(integrated),ranking=safe(ranking),evaluation=safe(evaluation),external_predictions=safe(predictions),generated=as.character(Sys.time()),session=utils::capture.output(sessionInfo())))
 methods_text<-function(x) c('SURFACEOME PIPELINE V3 — ANALYSIS PROVENANCE',paste('Timestamp:',x$generated),'All simulated datasets and annotations explicitly labelled in UI.','ADT: within-sample CLR = log2(count+1) minus sample mean log2(count+1); or log2(CPM+1). All logFC on log2 scale.','Proteomics: log2(x+1) unless prelogged, then per-sample median centering.','Differential: limma eBayes, unpaired two-group design Tumor vs Normal, BH FDR.','Integration: exact feature IDs; no automatic UniProt/HGNC mapping.','Prioritization: BH FDR gate (significant in >=1 modality), concordance = up in both ADT and MS; normalized weighted score; logFC proxy for specificity; annotation evidence user-provided or simulated.','Evaluation: pROC AUROC, PRROC AUPRC, Brier mean squared probability error.','Stratified bootstrap 95% percentile CI; paired DeLong when scoreB present.','ML: GLM with train-only scaling, stratified hold-out or leave-one-cohort-out.','Baseline: training outcome prevalence. No hyperparameter optimization.',paste('ADT source:',if(is.null(x$adt)) 'not run' else x$adt$source),paste('MS source:',if(is.null(x$ms)) 'not run' else x$ms$source),paste('Ranked targets:',if(is.null(x$ranking)) 'none' else nrow(x$ranking)),'','R SESSION:',x$session)
 output$pdf<-shiny::downloadHandler(filename=function()paste0('surfaceome_figures_',Sys.Date(),'.pdf'),content=function(file){x<-payload();grDevices::pdf(file,width=10,height=7,useDingbats=FALSE);on.exit(grDevices::dev.off(),add=TRUE);count<-0L
  emit<-function(p){print(p);count<<-count+1L}
  if(!is.null(x$adt)){emit(plot_distribution(x$adt$norm,'ADT'));emit(plot_heatmap(x$adt$norm));emit(plot_volcano(x$adt$diff,'ADT differential'))}
  if(!is.null(x$ms)){emit(plot_distribution(x$ms$norm,'MS'));emit(plot_heatmap(x$ms$norm));emit(plot_volcano(x$ms$diff,'MS differential'))}
  if(!is.null(x$integration))emit(plot_concordance(x$integration))
  if(!is.null(x$evaluation)){emit(plot_roc(x$evaluation));emit(plot_pr(x$evaluation));emit(plot_calibration(x$evaluation))}
  if(count==0) {graphics::plot.new();graphics::text(.5,.5,'No completed analyses to export')}
 })
 output$bundle<-shiny::downloadHandler(filename=function()paste0('surfaceome_analysis_',Sys.Date(),'.rds'),content=function(file)saveRDS(payload(),file))
 output$targets<-shiny::downloadHandler(filename=function()paste0('surfaceome_targets_',Sys.Date(),'.csv'),content=function(file){shiny::req(ranking());utils::write.csv(ranking(),file,row.names=FALSE)})
 output$methods<-shiny::downloadHandler(filename=function()paste0('surfaceome_methods_',Sys.Date(),'.txt'),content=function(file)writeLines(methods_text(payload()),file))
 output$provenance<-shiny::renderText({paste('R version:',getRversion(),'| Export includes sessionInfo, actual objects and simulated-data flags.')})
})
