server <- function(input, output, session) {
  
  #Load functions
  casecontrol_cache <<- data.frame("Coverage"=numeric(), "nr_cases"=numeric(),  "r"=numeric(),"VE"=numeric(), "VE_lower"=numeric())
  cohort_cache <<- data.frame("Coverage"=numeric(), "SS"=numeric(),  "AR"=numeric(),"VE"=numeric(), "VE_lower"=numeric())
  MAX_CORES = 6
  SPLINE_DF = 6
  
  current_data = NULL
  study = results = analysis = VE = NULL
  
  source("src/aux_functions.R", local = T)
  for(file_name in list.files(path="src/samp_size/", pattern = "*.R", full.names = T)){
    source(file_name, local = T)  
  }
  
  options(shiny.maxRequestSize=100*1024^2) 
  
  validateValue = function(curValue, min, max, name){
    res = !is.na(curValue) & is.numeric(curValue) & curValue >= min & curValue <= max
    if(res==F){
      shinyalert(text = paste0("'",name, "' should be a number between ", min, " and ", max), 
                 type = "error", 
                 title = "Invalid parameter value")
    }
    return(res)
  }
  
  validateMultipleValues = function(curValues, min, max, name){
    res = !sapply(curValues, is.na) & sapply(curValues, is.numeric) & curValues >= min & curValues <= max
    
    if(any(res==F)){
      shinyalert(text = paste0("All '",name, "' should be numbers between ", min, " and ", max), 
                 type = "error",
                 title = 'Invalid parameter values')
    }
    return(all(res))
  }
  
  validateSequence = function(min, max, step, name){
    res = validateValue(min, 1, 2e6, paste("Minimum", name)) &
      validateValue(max, 1, 2e6, paste("Maximum", name)) & 
      validateValue(step, 0, 2e6, paste("Step", name))
    
    if(res==T){
      res = validateValue(min, 1, max, paste("Minimum", name))
    }
    
    return(res)
  }
  
  #When you press calculate this is what happens
  #probably the worst if-else I have written
  #using switch is a better strategy than if-else
  observeEvent(input$calcButton, {
    updateTabsetPanel(session, "myTabs", selected = "Calculation")
    
    #Global
    study <<- isolate(input$study)
    results <<- isolate(input$results)
    analysis <<- isolate(input$analysis)
    VE <<- isolate(input$VE)
    current_data <<- NULL
    
    #Local
    power = input$power
    alpha = input$alpha
    
    if(validateValue(alpha, 0 , 1, "Significance level")){
      
      AR=input$AR
      propbrand_cohort = as.numeric(unlist(strsplit(input$propbrand_cohort,",")))
      COV_cohort = input$COV_cohort
      COVrange_cohort=as.numeric(unlist(strsplit(input$COVrange_cohort,",")))
      
      if(study=="Cohort"){
        
        if(validateValue(AR, 0, 100, "Attack rate") & 
           validateSequence(input$minSizeCohort, input$maxSizeCohort, input$bySizeCohort, "Total sample size")){
          
          N_cohort = seq(from = input$minSizeCohort, to = input$maxSizeCohort, by = input$bySizeCohort)
          
          if(analysis=="MindetVE"){
            if(validateValue(power, 0 , 1, "Power")){
              if (results=="Overall"){
                if(validateMultipleValues(COVrange_cohort, 0, 100, "Vaccination coverage in population")){
                  current_data <<- Cohort_overall_mindetVE(power=power, alpha=alpha, COV=COVrange_cohort, AR=AR,N=N_cohort) 
                  current_data$minVE <<- round(current_data$minVE, 2)
                }
              }else if (results=="Brand-specific"){
                if(validateValue(COV_cohort,0, 100, "Vaccination coverage in population") &
                   validateMultipleValues(propbrand_cohort, 0, 100, "Brand proportions")){
                  current_data <<- Cohort_brandspecific_mindetVE(power=power, alpha=alpha, COV=COV_cohort, AR=AR, propbrand = propbrand_cohort, N=N_cohort)
                  current_data$minVE <<- round(current_data$minVE,2)
                }
              }
            }
          }else if(analysis=="Precision"){
            if(validateValue(VE, 0, 100, "Anticipated VE")){
              if (results=="Overall"){
                current_data <<- Cohort_overall_precision(alpha = alpha, VE = VE, COV = COVrange_cohort, AR = AR, N=N_cohort)
              }else if (results=="Brand-specific"){
                if(validateValue(COV_cohort,0, 100, "Vaccination coverage in population") &
                   validateMultipleValues(propbrand_cohort, 0, 100, "Brand proportions")){
                  current_data <<- Cohort_brandspecific_precision(alpha = alpha, VE = VE, COV=COV_cohort, AR=AR, propbrand = propbrand_cohort, N=N_cohort)
                }
              }
            }
          }
        }
        
      }else if(study=="Casecontrol"){
        
        r=input$r
        propbrand_casecontrol = as.numeric(unlist(strsplit(input$propbrand_casecontrol,",")))
        COV_casecontrol = input$COV_casecontrol
        COVrange_casecontrol = as.numeric(unlist(strsplit(input$COVrange_casecontrol,",")))
        
        if(validateValue(r, 1, 100, "Number of controls per Case") &
           validateSequence(input$minSizeCaseControl, input$maxSizeCaseControl, input$bySizeCaseControl, "Number of cases")){
          
          Ncases_casecontrol = seq(from = input$minSizeCaseControl, to = input$maxSizeCaseControl, by = input$bySizeCaseControl)
          
          if(analysis=="MindetVE"){
            if(validateValue(power, 0 , 1, "Power")){
              if (results=="Overall"){
                if(validateMultipleValues(COVrange_casecontrol, 0, 100, "Vaccination coverage in controls")){
                  current_data <<- Casecontrol_overall_mindetVE(power=power, alpha=alpha, COV = COVrange_casecontrol, r = r, ncases = Ncases_casecontrol)
                  current_data$minVE <<- round(current_data$minVE, 2)
                }
              }else if (results=="Brand-specific"){
                if(validateValue(COV_casecontrol,0, 100, "Vaccination coverage in controls") &
                   validateMultipleValues(propbrand_casecontrol, 0, 100, "Brand proportions")){
                  current_data <<- Casecontrol_brandspecific_mindetVE(power=power, alpha=alpha, COV = COV_casecontrol, r = r, propbrand = propbrand_casecontrol, ncases=Ncases_casecontrol)
                  current_data$minVE <<- round(current_data$minVE, 2)
                }
              }
            }
          }else if(analysis=="Precision"){
            if(validateValue(VE, 0, 100, "Anticipated VE")){
              if (input$results=="Overall"){
                current_data <<- Casecontrol_overall_precision(alpha = alpha, VE = VE, COV = COVrange_casecontrol, r = r, ncases = Ncases_casecontrol)
              } else if (results=="Brand-specific"){
                if(validateValue(COV_casecontrol,0, 100, "Vaccination coverage in controls") &
                   validateMultipleValues(propbrand_casecontrol, 0, 100, "Brand proportions")){
                  current_data <<- Casecontrol_brandspecific_precision(alpha = alpha, VE = VE, COV = COV_casecontrol, r = r, propbrand = propbrand_casecontrol, ncases=Ncases_casecontrol)
                }
              }
            }
          }
          
        }
        
      }
    }
    
    updatePlot()
    updateTable()
  })
  
  updatePlot = function(){
    output$plots<-renderUI({
      if(is.null(current_data)){
        return(NULL)
      }
      output$plot <- renderHighchart({ 
        export <- list(
          list(
            text = "PNG",
            onclick = JS("function () {
                       this.exportChart({ type: 'image/png' }); }")
          ),
          list(
            text = "JPEG",
            onclick = JS("function () {
                       this.exportChart({ type: 'image/jpeg' }); }")
          ),
          list(
            text = "TIFF",
            onclick = JS("function () {
                       this.exportChart({ type: 'image/tiff' }); }")
          ),
          list(
            text = "PDF",
            onclick = JS("function () {
                       this.exportChart({ type: 'application/pdf' }); }")
          )
        )
        
        if (study=="Cohort" & results=="Overall" & analysis=="MindetVE"){ #1 Cohort overall min detectable IVE
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "SS", y = "minVE", group= "Coverage"), tooltip = list(headerFormat = "", pointFormat = "Total Sample Size = {point.x}<br>MD VE = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Total Sample Size"))
          hc <- hc_yAxis(hc = hc, title = list(text = "Minimum Detectable VE [%]"),min=0, max=100)
          hc <- hc_legend(hc = hc, title = list(text = "Vaccination Coverage [%]"))}
        else if (study=="Cohort" & results=="Brand-specific" & analysis=="MindetVE"){ #2 Cohort brand-specific min detectable IVE
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "SS", y = "minVE", group= "prop_brand"), tooltip = list(headerFormat = "", pointFormat = "Total Sample Size = {point.x}<br>MD BS VE = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Total Sample Size"))
          hc <- hc_yAxis(hc = hc, title = list(text = "Minimum Detectable Brand-Specific VE [%]"),min=0, max=100)
          hc <- hc_legend(hc = hc, title = list(text = "Brand Proportions [%]"))}
        else if (study=="Cohort" & results=="Overall" & analysis=="Precision"){ #3 Cohort overall precision
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "SS", y = "VE_lower", group= "Coverage"), tooltip = list(headerFormat = "", pointFormat = "Total Sample Size = {point.x}<br>VE 95% CI lower = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Total Sample Size"))
          hc <- hc_yAxis(hc = hc, title = list(text = "VE Lower 95% CI [%]"), plotLines = list(list(label = NULL, color = "#000000", width = 3,value = VE)))
          hc <- hc_legend(hc = hc, title = list(text = "Vaccination Coverage [%]"))}
        else if (study=="Cohort" & results=="Brand-specific" & analysis=="Precision"){ #4 Cohort brand-specific precision
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "SS", y = "VE_lower", group= "prop_brand"), tooltip = list(headerFormat = "", pointFormat = "Total Sample Size = {point.x}<br>BS VE 95% CI lower = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Total Sample Size"))
          hc <- hc_yAxis(hc = hc, title = list(text = "Brand-Specific VE Lower 95% CI [%]"), plotLines = list(list(label = NULL, color = "#000000", width = 3,value = VE)))
          hc <- hc_legend(hc = hc, title = list(text = "Brand Proportions [%]"))}
        else if (study=="Casecontrol" & results=="Overall" & analysis=="MindetVE"){ #5 TND overall min detectable IVE
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "nr_cases", y = "minVE", group= "Coverage"), tooltip = list(headerFormat = "", pointFormat = "Nr. Cases = {point.x}<br>MD VE = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Number of Cases"))
          hc <- hc_yAxis(hc = hc, title = list(text = "Minimum Detectable VE [%]"),min=0, max=100)
          hc <- hc_legend(hc = hc, title = list(text = "Vaccination Coverage [%]"))}
        else if (study=="Casecontrol" & results=="Brand-specific" & analysis=="MindetVE"){ #6 TND brand-specific min detectable IVE
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "nr_cases", y = "minVE", group= "prop_brand"), tooltip = list(headerFormat = "", pointFormat = "Nr. Cases = {point.x}<br>MD BS VE = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Number of Cases"))
          hc <- hc_yAxis(hc = hc, title = list(text = "Minimum Detectable Brand-Specific VE [%]"),min=0, max=100)
          hc <- hc_legend(hc = hc, title = list(text = "Brand Proportions [%]"))}
        else if (study=="Casecontrol" & input$results=="Overall" & analysis=="Precision"){ #7 TND overall precision
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "nr_cases", y = "VE_lower", group= "Coverage"), tooltip = list(headerFormat = "", pointFormat = "Nr. Cases = {point.x}<br>VE 95% CI lower = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Number of Cases"))
          hc <- hc_yAxis(hc = hc, title = list(text = "VE Lower 95% CI [%]"), plotLines = list(list(label = NULL, color = "#000000", width = 3,value = VE)))
          hc <- hc_legend(hc = hc, title = list(text = "Vaccination Coverage [%]"))
        } else{ #8 TND brand-specific precision
          hc <- highchart()
          hc <- hc_add_series(hc = hc, data = current_data, type = "line", hcaes(x= "nr_cases", y = "VE_lower", group= "prop_brand"), tooltip = list(headerFormat = "", pointFormat = "Nr. Cases = {point.x}<br>BS VE 95% CI lower = {point.y}%"))
          hc <- hc_xAxis(hc = hc, title = list(text = "Number of Cases"))
          hc <- hc_yAxis(hc = hc, title = list(text = "Brand-Specific VE Lower 95% CI [%]"), plotLines = list(list(label = NULL, color = "#000000", width = 3,value = VE)))
          hc <- hc_legend(hc = hc, title = list(text = "Brand Proportions [%]"))
        }
        
        hc <- hc_exporting(hc = hc, enabled = T)
        return(hc)
      })
      return(highchartOutput("plot"))
    }) 
  }
  
  tableData = function(){
    
    #Rename the columns
    if (study=="Cohort" & results=="Overall" & analysis=="MindetVE"){
      dataPlot = current_data
      colnames(dataPlot)<-c("Attack rate in the unvaccinated","Coverage", "Sample size", "Minimum Detectable VE")
    }
    else if (study=="Cohort" & results=="Brand-specific" & analysis=="MindetVE"){ 
      dataPlot = current_data
      colnames(dataPlot)<-c("Attack rate in the unvaccinated","Coverage", "Sample size", "Minimum Detectable VE", "Proportion of Brand")
    }
    else if (study=="Cohort" & results=="Overall" & analysis=="Precision"){
      dataPlot = current_data
      colnames(dataPlot)<-c("Coverage", "Sample size", "Attack rate in the unvaccinated","Anticipated VE", "VE 95% CI lower")
    }
    else if (study=="Cohort" & results=="Brand-specific" & analysis=="Precision"){
      dataPlot = current_data
      colnames(dataPlot)<-c("Coverage", "Sample size", "Attack rate in the unvaccinated", "Proportion of Brand", "Anticipated VE", "VE 95% CI lower")
    }
    else if (study=="Casecontrol" & results=="Overall" & analysis=="MindetVE"){
      dataPlot = current_data
      colnames(dataPlot)<-c("Coverage", "Number of cases", "Controls per case", "Minimum Detectable VE")
    }
    else if (study=="Casecontrol" & results=="Brand-specific" & analysis=="MindetVE"){
      dataPlot = current_data
      colnames(dataPlot)<-c("Coverage", "Proportion of Brand", "Number of cases", "Controls per case", "Minimum Detectable VE")
    }
    else if (study=="Casecontrol" & results=="Overall" & analysis=="Precision"){
      dataPlot = current_data
      colnames(dataPlot)<-c("Coverage", "Number of cases", "Controls per case", "Anticipated VE", "VE 95% CI lower")
    }
    else{
      dataPlot = current_data
      colnames(dataPlot)<-c("Coverage", "Number of cases", "Controls per case", "Proportion of Brand", "Anticipated VE", "VE 95% CI lower")
    }
    
    return(dataPlot)
  }
  
  updateTable = function(){
    output$tableBox <- renderUI({
      if(is.null(current_data)){
        return(NULL)
      }
      output$DTTable <- DT::renderDataTable({
        # if(is.null(dataUpload1())) {
        #   return(NULL)
        # }
        DT::datatable(tableData(),rownames = FALSE,options = list(paging = TRUE, searching = TRUE, info = FALSE, ordering = TRUE, scrollX = TRUE, 
                                                                  columnDefs = list(list(className = 'dt-head-center dt-center', targets = "_all"))))
      })
      return(box(width = 12, DT::dataTableOutput("DTTable"), downloadLink('downloadData', 'Download Table')))
    })
  }
  
  #this is still reactive. won't change it
  output$mainPanelTopMessage <- renderUI({
    if(is.null(current_data)){
      return(HTML('<h4 style= "margin-left: 15px;">Press the <b>Calculate</b> button in the lower left corner to perform the sample size calculations</h4>'))
    }else{
      return(
        HTML(
          paste0(
            '<h4 style= "margin-left: 15px;">',
            '<b>Sample size calculations for: </b>', 
            ifelse(analysis=="MindetVE", 
                   "Minimum detectable vaccine effectiveness (MD VE)", 
                   "Expected lower limit of the VE 95% confidence interval (VE 95% CI lower)"
            ),
            '<br/>',
            '<b> Study design: </b>', 
            ifelse(
              study=="Cohort", 
              "Cohort study", 
              "Test negative design study"
            ),
            '<br/>',
            '<b>Outcome of interest: </b>', 
            ifelse(
              results=="Overall", 
              "Overall VE", 
              "Brand-specific VE"
            ), 
            '</h4>'
          )
        )
      )
    }
  })
  
  # output$downLink <- renderUI({
  #   if(is.null(dataPlot())) {
  #     return(NULL)
  #   }
  #   downloadLink('downloadData', 'Download Table')
  # })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('DRIVE_SS_', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      df <- tableData()
      df$Design <- study
      df$Outcome <- results
      df$'Significance Level' <- input$alpha
      df$Power <- input$power
      write.csv(df, con)
    }
  )
  
  output$contact <- renderUI({
    html <- "<div style='margin-top: 50px; margin-left: 15px;'>"
    html <- paste0(html, "<h4>Questions or Suggestions?</h4>")
    html <- paste0(html, "<p>Please contact: <a href='mailto:nick.desmedt@p-95.com'>nick.desmedt@p-95.com</a></p>")
    html <- paste0(html, "</div>")
    HTML(html)
  })
  
  output$info <- renderUI({
    html <- "<div style='margin: 15px;'>"
    html <- paste0(html, "<h3>Background & Citation</h4>")
    html <- paste0(html, "<p>This web application allows to perform sample size calculations for cohort and test negative design studies on brand-specific and overall vaccine effectiveness. Sample size requirements are expressed in terms of both minimal detectable vaccine effectiveness and precision (here defined as the expected lower limit of the 95% confidence interval of the anticipated Vaccine effectiveness). The application allows uses to input parameters (e.g. study design, case-control ratio, disease attack rate, overall and brand-specific vaccination coverage) and to obtain graphical and tabular sample size outputs.</p>")
    html <- paste0(html, "<p>This web application was developed as part of the DRIVE (Development of Robust and Innovative Vaccine Effectiveness) project, funded by the Innovative Medicines Initiative (IMI). DRIVE is a public-private partnership aiming to build capacity for estimating brand-specific influenza vaccine effectiveness (IVE) in Europe. This work is free to use upon acknowledgement by citing the following paper:</p>")
    html <- paste0(html, "paper: link coming soon<br>") 
    html <- paste0(html, "user manual:  link coming soon<br>")
    html <- paste0(html, "DRIVE website: <a href='https://www.drive-eu.org/'>https://www.drive-eu.org/</a>")
    html <- paste0(html, "</div>")
    HTML(html)
  })
  
  
}#End of server