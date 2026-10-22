# Surfaceome Pipeline v3 · reproducible R Shiny application
required <- c('shiny','bslib','ggplot2','DT','pROC','PRROC','limma')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install missing packages: ',paste(missing,collapse=', '),'\nlimma is installed through BiocManager.')
library(shiny);library(bslib);library(ggplot2)
for(f in c('utils_data','utils_metrics','utils_plotting','mod_assays','mod_integrate','mod_prioritize','mod_ml','mod_eval','mod_export')) source(file.path('R',paste0(f,'.R')),local=FALSE)
app_theme<-bs_theme(version=5,primary='#2563eb',secondary='#f97316',success='#16a34a')
css<-'body{background:#f8fafc}.metric-box{border-radius:12px;padding:14px;background:#fff;border:1px solid #e2e8f0}.sidebar-step{display:block;width:100%;text-align:left;margin-bottom:6px}.app-header{background:white;border-bottom:1px solid #e2e8f0;padding:12px 20px}.app-header strong{color:#2563eb}'
ui<-page_fluid(theme=app_theme,tags$head(tags$style(HTML(css))),
 div(class='app-header',strong('Surfaceome Pipeline'),span(class='badge bg-primary ms-2','v3.0'),span(class='badge bg-warning text-dark ms-2','Research use · not clinically validated')),
 layout_sidebar(sidebar=sidebar(width=260,open='desktop',h6('PIPELINE STEPS',class='text-muted'),
 actionButton('nav_overview','0. Overview',class='sidebar-step'),actionButton('nav_adt','1. CITE-seq / ADT',class='sidebar-step'),
 actionButton('nav_ms','2. Mass Spec',class='sidebar-step'),actionButton('nav_integrate','3. Integrate',class='sidebar-step'),
 actionButton('nav_targets','4. Target Ranking',class='sidebar-step'),actionButton('nav_ml','5. External Validation',class='sidebar-step'),
 actionButton('nav_eval','6. Evaluation',class='sidebar-step'),actionButton('nav_export','7. Export',class='sidebar-step')),
 navset_hidden(id='pages',selected='overview',
 nav_panel_hidden('overview',div(class='p-3',h3('Surfaceome Pipeline · Seven analytical layers'),p('Reproducible CITE-seq and proteomics integration, model evaluation, and therapeutic target prioritization.'),
 card(card_header('Input schemas'),card_body(tags$ul(tags$li('Feature matrix CSV: first column feature ID; remaining columns sample IDs. Alternatively a matrix RDS.'),tags$li('Metadata CSV: sample_id, group (Tumor or Normal).'),tags$li('Annotations CSV: feature, surface_evidence (0–1), normal_tissue_risk (0–1).'),tags$li('Predictions CSV: truth, scoreA, optional scoreB (probabilities).'),tags$li('ML table: truth, cohort, numerical feature columns.')))),
 card(card_header('Important methodological limits'),card_body('Exact-ID integration only. Feature effects are proxy measurements; tissue-level evidence must be supplied. No patient data or external reference is uploaded by this app. Demonstration samples are simulated. No diagnostic or clinical decisions should be based on scores.')))),
 nav_panel_hidden('adt',div(class='p-3',h4('CITE-seq QC / Normalization'),mod_assay_ui('adt_mod','ADT matrix & differential testing',TRUE))),
 nav_panel_hidden('ms',div(class='p-3',h4('Mass Spectrometry'),mod_assay_ui('ms_mod','MS abundance & differential testing',FALSE))),
 nav_panel_hidden('integrate',div(class='p-3',h4('Multi-omics Integration'),mod_integrate_ui('integrate_mod'))),
 nav_panel_hidden('targets',div(class='p-3',h4('Target Prioritization'),mod_prioritize_ui('target_mod'))),
 nav_panel_hidden('ml',div(class='p-3',h4('Independent / Cohort-held-out Validation'),mod_ml_ui('ml_mod'))),
 nav_panel_hidden('eval',div(class='p-3',h4('Model Evaluation'),mod_eval_ui('eval_mod'))),
 nav_panel_hidden('export',div(class='p-3',h4('Methods & Export'),mod_export_ui('export_mod'))))))
server<-function(input,output,session){
 for(key in c('overview','adt','ms','integrate','targets','ml','eval','export')) local({k<-key;observeEvent(input[[paste0('nav_',k)]],nav_select('pages',k))})
 adt<-mod_assay_server('adt_mod',TRUE)
 ms<-mod_assay_server('ms_mod',FALSE)
 integrated<-mod_integrate_server('integrate_mod',adt,ms)
 targets<-mod_prioritize_server('target_mod',integrated)
 predicted<-mod_ml_server('ml_mod')
 eval<-mod_eval_server('eval_mod',predictions=predicted)
 mod_export_server('export_mod',adt,ms,integrated,targets,eval,predicted)
}
shinyApp(ui,server)
