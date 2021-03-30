dashboardHeader <- function(..., title = NULL, disable = FALSE,title.navbar=NULL, .list = NULL) {
  items <- c(list(...), .list)
  #lapply(items, tagAssert, type = "li", class = "dropdown")
  tags$header(
    class = "main-header",
    
    style = if (disable) "display: none;",
    
    span(class = "logo", title),
    
    tags$nav(
      class = "navbar navbar-static-top", role = "navigation",
      
      # Embed hidden icon so that we get the font-awesome dependency
      span(shiny::icon("bars"), style = "display:none;"),
      
      # Sidebar toggle button
      a(
        href="#", class="sidebar-toggle",
        `data-toggle`="offcanvas",
        role="button",
        span(class="sr-only", "Toggle navigation")
      ),
      
      title.navbar,
      
      div(
        class = "navbar-custom-menu",
        tags$ul(
          class = "nav navbar-nav",
          items
        )
      )
    )
  )
}

title.navbar.html <- tags$div(
  style="margin-left: 20%; color:#29405B; text-align:center; vertical-align:middle; display:inline-block; font-size: 200%; font-weight: 400; padding-left: 5px; padding-top: 5px;",
  span(style="color: #FFFFFF; font: bold", "Sample Size Calculation for Vaccine Effectiveness Studies")
)

header <- dashboardHeader(
  title = img(src="DRIVE_LOGO_LAND_WC.png", height =  40),
  title.navbar = title.navbar.html
)

# Define UI for data upload app ----
ui <- fluidPage(
  useShinyalert(),
  dashboardPage(
    # App title ----
    header,
    
    # Sidebar panel for inputs ----
    dashboardSidebar(
      tags$head(
        tags$link(rel = "stylesheet", 
                  type = "text/css", 
                  href = "BI1.css"),
        tags$script(type="text/javascript", src = "busy.js")
      ),
      
      width = 400,
      tags$hr(),
      h4("Options"),
      #Initial choices for the client
      radioButtons("study", "Study Design",
                   choices = c("Cohort study" = "Cohort",
                               "Test negative design study" = "Casecontrol"),
                   #"Case-Control Study" = "Casecontrol"), 
                   selected = "Cohort"),
      radioButtons("analysis", "Sample size calculations for",
                   choices = c("Minimum detectable vaccine effectiveness (MD VE)" = "MindetVE",
                               "Expected lower limit of the VE 95% confidence interval (VE 95% CI lower)" = "Precision"),
                   selected = "MindetVE", width=350),
      radioButtons("results", "Outcome of interest",
                   choices = c("Overall VE" = "Overall",
                               "Brand-specific VE (BS VE)" = "Brand-specific"),
                   selected = "Overall"),
      
      tags$hr(),
      h4("Parameters"),
      
      fluidRow(
        column(
          width = 6, 
          numericInput("alpha","Significance level [0-1]",value = 0.05, min = 0.01, max = 1)
        ),
        column(
          width = 6, 
          conditionalPanel(
            condition = "input.analysis == 'MindetVE'", 
            numericInput("power","Power [0-1]",value = 0.80, min = 0.01, max = 1)
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.study == 'Cohort'", 
        numericInput("AR","Attack rate in the unvaccinated [%]",value = 5, min = 1, max = 100)
      ),
      
      conditionalPanel(
        condition = "input.study == 'Casecontrol'", 
        numericInput("r","Number of controls per case", value = 1, min = 1, max = 100)
      ),
      conditionalPanel(
        condition = "input.analysis == 'Precision'",
        numericInput("VE","Anticipated VE [%]", value = 50, min = 1, max = 100)
      ),
      
      conditionalPanel(
        condition = "input.results == 'Overall' & input.study == 'Casecontrol'",
        textInput("COVrange_casecontrol","Vaccination coverage in controls [%]", width='88%', value="5, 20, 50, 70")
      ),
      
      conditionalPanel(
        condition = "input.results == 'Overall' & input.study == 'Cohort'",
        textInput("COVrange_cohort","Vaccination coverage in population [%]", width='88%', value="5, 20, 50, 70")
      ),
      
      conditionalPanel(
        condition = "input.results == 'Brand-specific' & input.study =='Casecontrol'", 
        numericInput ("COV_casecontrol","Vaccination coverage in controls [%]", value =50, min = 1, max= 100),
        textInput("propbrand_casecontrol","Brand proportions [%]", width='88%', value="10,20,40,60,80,100")
      ),
      
      conditionalPanel(
        condition = "input.results == 'Brand-specific' & input.study == 'Cohort'", 
        numericInput ("COV_cohort","Vaccination coverage in population [%]", value =50, min = 1, max= 100),
        textInput("propbrand_cohort","Brand proportions [%]", width='88%', value="10,20,40,60,80,100")
      ),
      
      conditionalPanel(
        condition = "input.study=='Cohort'", 
        tags$label("Total sample size"), 
        tags$br(),
        
        fluidRow(
          column(
            width = 4,
            numericInput(inputId = "minSizeCohort", label = "Min.", value = 10000, min = 0, max = 2e6, step = 1)
          ),
          column(
            width = 4,
            numericInput(inputId = "maxSizeCohort", label = "Max.", value = 50000, min = 0, max = 2e6, step = 1)
          ),
          column(
            width = 4,
            numericInput(inputId = "bySizeCohort", label = "Step", value = 500, min = 0, max = 10000, step = 1)
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.study=='Casecontrol'", 
        tags$label("Number of cases"), 
        tags$br(),
        fluidRow(
          column(
            width = 4,
            numericInput(inputId = "minSizeCaseControl", label = "Min.", value = 100, min = 0, max = 2e6, step = 1)
          ),
          column(
            width = 4,
            numericInput(inputId = "maxSizeCaseControl", label = "Max.", value = 500, min = 0, max = 2e6, step = 1)
          ),
          column(
            width = 4,
            numericInput(inputId = "bySizeCaseControl", label = "Step", value = 25, min = 0, max = 10000, step = 1)
          )
        )
      ),
      
      actionButton(
        inputId = "calcButton", label = "Calculate", 
        class = "btn btn-primary", 
        style="color: white; margin-left: 0px",
      ),
      
      uiOutput("contact"),
      img(src="logo_P95.png", height =  80, style = "margin-left: 15px")
      
      #, submitButton("Submit")
    ), #End of dashboardSidebar
    dashboardBody( 
      tabsetPanel(
        id = "myTabs",
        tabPanel(
          "Calculation", 
          div(class = "busy", style = "font-size:150%;background:rgba(0,0,0,.25);color:#222222;position:fixed;bottom: 1%;right:1%;overflow: visible;z-index: 999;",
              HTML("Calculation in progress..."),
              img(src="ajaxloaderq.gif")
          ),
          uiOutput("mainPanelTopMessage"), 
          box(width=12, uiOutput("plots")), 
          uiOutput("tableBox")
        ),
        tabPanel("Background Info", uiOutput("info"))
      )
    )
  )# End of dashboardPage
)#End of ShinyUI
