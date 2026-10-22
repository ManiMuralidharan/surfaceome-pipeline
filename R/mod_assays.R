mod_assay_ui <- function(id,label,is_adt=TRUE) {
 ns<-shiny::NS(id)
 bslib::card(bslib::card_header(label),bslib::card_body(
  shiny::fileInput(ns('matrix'),'Feature × sample CSV / matrix RDS',accept=c('.csv','.rds')),
  shiny::fileInput(ns('meta'),'Sample metadata CSV (sample_id, group)',accept='.csv'),
  if(!is_adt) shiny::checkboxInput(ns('logged'),'Input intensities already log2 transformed',FALSE),
  shiny::selectInput(ns('method'),'Normalization',if(is_adt)c('CLR','logCPM') else 'Median-centered log2'),
  shiny::textInput(ns('case'),'Case group','Tumor'),shiny::textInput(ns('control'),'Control group','Normal'),
  shiny::actionButton(ns('run'),'Process & fit differential model',class='btn-primary'),
  shiny::actionButton(ns('demo'),'Load clearly labelled DEMO data',class='btn-outline-secondary'),
  shiny::verbatimTextOutput(ns('status')),shiny::tableOutput(ns('qc')),
  bslib::layout_columns(shiny::plotOutput(ns('dist'),height=270),shiny::plotOutput(ns('volcano'),height=270)),
  shiny::plotOutput(ns('heatmap'),height=360),DT::DTOutput(ns('table'))))
}
mod_assay_server <- function(id,is_adt=TRUE) shiny::moduleServer(id,function(input,output,session){
 state<-shiny::reactiveVal(NULL); err<-shiny::reactiveVal('Upload assay and metadata, then run.');
 process<-function(m,md,origin) {
  md<-validate_metadata(md,colnames(m)); norm<-if(is_adt) normalize_adt(m,input$method) else normalize_ms(m,isTRUE(input$logged))
  diff<-fit_differential(norm,md,input$case,input$control)
  state(list(raw=m,norm=norm,meta=md,diff=diff,source=origin,type=if(is_adt)'ADT' else 'MS',method=if(is_adt)input$method else 'log2 + median center'))
  err(paste('Processed',nrow(m),'features across',ncol(m),'samples. Source:',origin))
 }
 shiny::observeEvent(input$run,{
  tryCatch({shiny::req(input$matrix,input$meta);m<-read_feature_matrix(input$matrix$datapath,input$matrix$name);md<-utils::read.csv(input$meta$datapath);process(m,md,'USER DATA')},error=function(e){state(NULL);err(paste('ERROR:',conditionMessage(e)))})
 })
 shiny::observeEvent(input$demo,{
  tryCatch({set.seed(if(is_adt) 101 else 102); n<-120;p<-12; groups<-rep(c('Normal','Tumor'),each=6);m<-matrix(if(is_adt)stats::rpois(n*p,20) else stats::rlnorm(n*p,8,.6),nrow=n,dimnames=list(paste0('P',seq_len(n)),paste0('S',seq_len(p))));m[1:15,groups=='Tumor']<-m[1:15,groups=='Tumor']*if(is_adt)3 else 2;md<-data.frame(sample_id=colnames(m),group=groups);process(m,md,'SIMULATED DEMO — NOT BIOLOGICAL RESULTS')},error=function(e){state(NULL);err(paste('ERROR:',conditionMessage(e)))})
 })
 output$status<-shiny::renderText(err())
 output$qc<-shiny::renderTable({x<-shiny::req(state());data.frame(Features=nrow(x$raw),Samples=ncol(x$raw),Missing=round(mean(is.na(x$raw))*100,2),Significant=sum(x$diff$adj.P.Val<.05,na.rm=TRUE))})
 output$dist<-shiny::renderPlot({x<-shiny::req(state());plot_distribution(x$norm,x$type)})
 output$volcano<-shiny::renderPlot({x<-shiny::req(state());plot_volcano(x$diff,paste(x$type,'differential'))})
 output$heatmap<-shiny::renderPlot({plot_heatmap(shiny::req(state())$norm)})
 output$table<-DT::renderDT({DT::datatable(shiny::req(state())$diff,options=list(pageLength=8,scrollX=TRUE))})
 state
})
