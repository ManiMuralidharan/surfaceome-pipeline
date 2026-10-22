mod_eval_ui <- function(id) {ns<-shiny::NS(id);shiny::tagList(
 bslib::card(bslib::card_header('External predictions & model comparison'),bslib::card_body(
 shiny::fileInput(ns('file'),'Prediction CSV: truth, scoreA, optional scoreB',accept='.csv'),
 shiny::actionButton(ns('demo'),'Use SIMULATED predictions'),
 shiny::sliderInput(ns('threshold'),'Probability decision threshold',0,1,.5,step=.01),
 shiny::numericInput(ns('boot'),'Stratified bootstrap replicates',200,min=20,max=5000,step=20),
 shiny::verbatimTextOutput(ns('status')))),
 bslib::layout_columns(
 bslib::value_box('AUROC',shiny::textOutput(ns('auc')),theme='primary'),
 bslib::value_box('AUPRC',shiny::textOutput(ns('pr')),theme='warning'),
 bslib::value_box('Brier',shiny::textOutput(ns('brier')),theme='primary'),
 bslib::value_box('Paired DeLong',shiny::textOutput(ns('delong')),theme='success')),
 bslib::layout_columns(bslib::card(shiny::plotOutput(ns('roc_plot'))),bslib::card(shiny::plotOutput(ns('pr_plot'))),bslib::card(shiny::plotOutput(ns('cal_plot')))),
 bslib::card(bslib::card_header('95% stratified bootstrap confidence intervals'),shiny::tableOutput(ns('ci'))),
 bslib::card(bslib::card_header('Confusion matrix & threshold metrics'),shiny::tableOutput(ns('conf'))))}
mod_eval_server <- function(id,predictions=NULL) shiny::moduleServer(id,function(input,output,session){
 data<-shiny::reactiveVal(NULL); msg<-shiny::reactiveVal('Upload predictions or use simulated data')
 shiny::observeEvent(input$file,{tryCatch({z<-check_predictions(utils::read.csv(input$file$datapath));data(z);msg(paste('USER PREDICTIONS ·',nrow(z),'records'))},error=function(e){data(NULL);msg(paste('ERROR:',conditionMessage(e)))})})
 shiny::observeEvent(input$demo,{set.seed(22);n<-400;y<-stats::rbinom(n,1,.1);p1<-stats::plogis(-2.4+1.7*y+stats::rnorm(n,sd=1));p2<-stats::plogis(-2.4+.9*y+stats::rnorm(n,sd=1));data(data.frame(truth=y,scoreA=p1,scoreB=p2));msg('SIMULATED PREDICTIONS — illustrative only')})
 if(!is.null(predictions)) shiny::observeEvent(predictions(),{if(!is.null(predictions())) tryCatch({data(check_predictions(predictions()));msg('Predictions passed from model validation module')},error=function(e){data(NULL);msg(paste('ERROR:',conditionMessage(e)))})})
 result<-shiny::reactive({shiny::req(data());tryCatch(evaluate_scores(data(),input$threshold,input$boot),error=function(e){msg(paste('ERROR:',conditionMessage(e)));NULL})})
 output$status<-shiny::renderText(msg())
 output$auc<-shiny::renderText({sprintf('%.3f',shiny::req(result())$auc)})
 output$pr<-shiny::renderText({sprintf('%.3f',shiny::req(result())$auprc)})
 output$brier<-shiny::renderText({sprintf('%.3f',shiny::req(result())$brier)})
 output$delong<-shiny::renderText({r<-shiny::req(result());if(is.null(r$comparison)) 'Need scoreB' else sprintf('Δ %.3f; p=%.3g',r$comparison$delta,r$comparison$pvalue)})
 output$roc_plot<-shiny::renderPlot(plot_roc(shiny::req(result())))
 output$pr_plot<-shiny::renderPlot(plot_pr(shiny::req(result())))
 output$cal_plot<-shiny::renderPlot(plot_calibration(shiny::req(result())))
 output$ci<-shiny::renderTable({r<-shiny::req(result());data.frame(metric=c('AUROC','AUPRC','Brier'),estimate=c(r$auc,r$auprc,r$brier),low=r$ci[,1],high=r$ci[,2])},digits=4)
 output$conf<-shiny::renderTable(shiny::req(result())$confusion,digits=3)
 result
})
