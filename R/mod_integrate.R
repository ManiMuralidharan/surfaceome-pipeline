mod_integrate_ui <- function(id) {
 ns<-shiny::NS(id)
 bslib::card(bslib::card_header('Cross-modality integration'),bslib::card_body(
  shiny::helpText('Exact feature-ID overlap. Map UniProt/HGNC identifiers before integration; unmatched features are excluded.'),
  shiny::verbatimTextOutput(ns('summary')),shiny::plotOutput(ns('corr'),height=340),DT::DTOutput(ns('table'))))
}
mod_integrate_server <- function(id,adt,ms) shiny::moduleServer(id,function(input,output,session){
 out<-shiny::reactive({a<-shiny::req(adt());b<-shiny::req(ms());d<-merge(a$diff,b$diff,by='feature',suffixes=c('_ADT','_MS'));if(!nrow(d)) stop('No matching feature IDs across modalities');d})
 output$summary<-shiny::renderText({d<-out();paste(nrow(d),'shared feature IDs; Spearman rho =',round(if(nrow(d)>2) stats::cor(d$logFC_ADT,d$logFC_MS,method='spearman',use='complete.obs') else NA_real_,3),'(effect-size concordance, not proof of surface localization)')})
 output$corr<-shiny::renderPlot({plot_concordance(out())})
 output$table<-DT::renderDT({DT::datatable(out(),options=list(pageLength=8,scrollX=TRUE))})
 out
})
