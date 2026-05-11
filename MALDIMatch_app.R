options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Check and install required packages
required_packages <- c("shiny", "shinyFiles", "shinyjs", "data.table", "ggplot2", "VennDiagram", "grid", "DT")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(shiny)
library(shinyFiles)
library(shinyjs)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(grid)
library(DT)

# Source the functions from the same directory
source("ComputeIon.R")
source("MALDIMatches.R")
source("PlotVenn.R")

custom_css <- "
body { background-color: #f5f5f5; }
.navbar { background-color: #1a1a1a !important; border-color: #ff8c00; }
.navbar-brand { color: #ff8c00 !important; font-weight: bold; font-size: 18px; }
.btn-primary { background-color: #ff8c00; border-color: #ff8c00; }
.btn-primary:hover { background-color: #e67e00; border-color: #e67e00; }
h1, h2, h3, h4 { color: #1a1a1a; }
.stats-box { background-color: white; border-left: 4px solid #ff8c00; padding: 15px; margin: 10px 0; border-radius: 3px; }
.stats-number { font-size: 24px; font-weight: bold; color: #ff8c00; }
#progress_container { display: none; }
#progress_container.show { display: block; }
.task-item { padding: 8px 0; font-size: 14px; color: #666; display: flex; align-items: center; }
.task-item.completed { color: #ff8c00; font-weight: bold; }
.task-checkmark { width: 20px; height: 20px; margin-right: 10px; background-color: #e0e0e0; border-radius: 50%; display: flex; align-items: center; justify-content: center; }
.task-item.completed .task-checkmark { background-color: #ff8c00; color: white; font-weight: bold; }
@keyframes spinner_rotate {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}
.spinner {
  display: inline-block;
  width: 16px;
  height: 16px;
  border: 3px solid #e0e0e0;
  border-top: 3px solid #ff8c00;
  border-radius: 50%;
  animation: spinner_rotate 0.8s linear infinite;
  margin-right: 8px;
}

"

ui <- navbarPage(
  "MALDImatch Explorer",
  header = tags$head(tags$style(HTML(custom_css))),
  shinyjs::useShinyjs(),
  tabPanel("Upload & Process",
    sidebarLayout(
      sidebarPanel(
        h3("File Upload", style = "color: #ff8c00;"),
        br(),
        fileInput("file_maldi", "Upload MALDI CSV (peaks file)", accept = c(".csv", "text/csv")),
        fileInput("file_lcms", "Upload msFragger File (LC-MS)", accept = c(".tsv", ".csv", ".txt", "text/plain", "text/csv")),
        hr(),
        h4("Matching Parameters", style = "color: #ff8c00;"),
        numericInput("tolerance_ppm", "Mass Tolerance (ppm):", value = 10, min = 1, max = 100, step = 1),
        checkboxInput("calc_isotopes", "Calculate Isotopic Distribution", value = TRUE),
        br(),
        actionButton("btn_process", "Process Files", class = "btn-lg btn-primary", width = "100%")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Summary",
            br(),
            h3("Processing Status", style = "color: #1a1a1a;"),
            uiOutput("status_ui"),
            br(), br(),
            h3("Statistics", style = "color: #1a1a1a;"),
            div(class = "row",
              div(class = "col-md-6", div(class = "stats-box",
                p("Total MALDI m/z values", style = "margin: 0; font-weight: bold; color: #666;"),
                p(textOutput("stat_maldi"), class = "stats-number"))),
              div(class = "col-md-6", div(class = "stats-box",
                p("Total LC-MS ions", style = "margin: 0; font-weight: bold; color: #666;"),
                p(textOutput("stat_lcms"), class = "stats-number")))),
            div(class = "row",
              div(class = "col-md-6", div(class = "stats-box",
                p("Total Matches Found", style = "margin: 0; font-weight: bold; color: #666;"),
                p(textOutput("stat_matches"), class = "stats-number"))),
              div(class = "col-md-6", div(class = "stats-box",
                p("Match Percentage", style = "margin: 0; font-weight: bold; color: #666;"),
                p(textOutput("stat_percentage"), class = "stats-number"))))),
          
          tabPanel("Matching Table",
            br(),
            downloadButton("btn_download_matches", "Download Matches CSV", class = "btn-primary"),
            br(), br(),
            DTOutput("matching_table")),
          
          tabPanel("Venn Diagram",
            br(),
            downloadButton("btn_download_venn", "Download Venn Diagram PNG", class = "btn-primary"),
            br(), br(),
            plotOutput("venn_plot", height = "600px", width = "600px"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  shinyjs::useShinyjs()
  
  reactive_data <- reactiveValues(
    maldi = NULL, lcms = NULL, lcms_annotated = NULL, matches = NULL, venn_plot = NULL, processing = FALSE
  )
  
  output$status_ui <- renderUI({
    div(textOutput("status_message"), br(),
      div(id = "progress_container",
        div(style = "background-color: #f0f0f0; padding: 15px; border-radius: 5px; border: 2px solid #ff8c00;",
          div(style = "display: flex; justify-content: space-between; margin-bottom: 8px;",
            span(style = "font-weight: bold; color: #1a1a1a;", "Progress"),
            span(id = "progress_percent", style = "color: #ff8c00; font-weight: bold;", "0%")),
          div(style = "background-color: #e0e0e0; height: 25px; border-radius: 3px; overflow: hidden;",
            div(id = "progress_bar_fill", style = "background-color: #ff8c00; height: 100%; width: 0%; transition: width 0.3s; display: flex; align-items: center; justify-content: center; color: white; font-weight: bold; font-size: 12px;", "")),
          br(),
          div(id = "progress_detail", style = "color: #666; font-size: 13px; margin-top: 10px;",
            span(id = "spinner", class = "spinner", style = "display: none;"),
            span(id = "detail_text", "Initializing...")),
          br(),
          div(style = "border-top: 1px solid #ddd; padding-top: 10px; margin-top: 10px;",
            div(class = "task-item", id = "task_1", div(class = "task-checkmark", "✓"), "Reading MALDI file"),
            div(class = "task-item", id = "task_2", div(class = "task-checkmark", "✓"), "Reading LC-MS file"),
            div(class = "task-item", id = "task_3", div(class = "task-checkmark", "✓"), "Computing ion masses & isotopes"),
            div(class = "task-item", id = "task_4", div(class = "task-checkmark", "✓"), "Matching MALDI peaks"),
            div(class = "task-item", id = "task_5", div(class = "task-checkmark", "✓"), "Creating Venn diagram"))))
    )
  })
  
  observeEvent(input$btn_process, {
    output$status_message <- renderText({
      if (is.null(input$file_maldi)) return("⚠️ Please upload a MALDI CSV file")
      if (is.null(input$file_lcms)) return("⚠️ Please upload an LC-MS file (TSV or CSV)")
      "Processing..."
    })
    
    reactive_data$processing <- TRUE
    shinyjs::addClass("progress_container", "show")
    
    tryCatch({
      # Read MALDI
      reactive_data$progress <- 10
      shinyjs::runjs("document.getElementById('progress_bar_fill').style.width = '10%'; document.getElementById('progress_percent').innerText = '10%'; document.getElementById('task_1').classList.add('completed');")
      
      maldi_raw <- readLines(input$file_maldi$datapath)
      maldi_clean <- maldi_raw[!grepl("^#", maldi_raw) & nzchar(trimws(maldi_raw))]
      maldi_text <- paste(maldi_clean, collapse = "\n")
      maldi_data <- tryCatch(
        read.csv(text = maldi_text, sep = ";", stringsAsFactors = FALSE),
        error = function(e) read.csv(text = maldi_text, sep = ",", stringsAsFactors = FALSE)
      )
      
      # Find m/z column with multiple patterns
      mz_col <- grep("m/z|m\\.z|^mz$", names(maldi_data), ignore.case = TRUE, value = TRUE)[1]
      if (is.na(mz_col)) stop("Could not find m/z column. Available: ", paste(names(maldi_data), collapse = ", "))
      
      maldi_mz <- as.numeric(maldi_data[[mz_col]])
      maldi_mz <- maldi_mz[!is.na(maldi_mz) & maldi_mz > 0]
      reactive_data$maldi <- maldi_mz
      
      # Read LC-MS
      reactive_data$progress <- 20
      shinyjs::runjs("document.getElementById('progress_bar_fill').style.width = '20%'; document.getElementById('progress_percent').innerText = '20%'; document.getElementById('task_2').classList.add('completed');")
      
      file_ext <- tolower(tools::file_ext(input$file_lcms$name))
      if (file_ext == "tsv" || file_ext == "txt") {
        lcms_data <- read.delim(input$file_lcms$datapath, stringsAsFactors = FALSE)
      } else {
        lcms_data <- read.csv(input$file_lcms$datapath, stringsAsFactors = FALSE)
      }
      reactive_data$lcms <- lcms_data
      
      # Compute ions
      reactive_data$progress <- 40
      shinyjs::runjs("document.getElementById('progress_bar_fill').style.width = '40%'; document.getElementById('progress_percent').innerText = '40%'; document.getElementById('spinner').style.display = 'inline-block'; document.getElementById('detail_text').innerHTML = '<strong>Computing ion masses and isotope distributions...</strong><br/>This step may take several minutes. Please be patient.';")
      
      mod_col_name <- grep("Modified.Sequence|Modified_Sequence|ModifiedSequence", names(lcms_data), ignore.case = TRUE, value = TRUE)[1]
      if (is.na(mod_col_name)) stop("Could not find Modified Sequence column")
      
      lcms_annotated <- ComputeIon(lcms_data, calc_isotopes = input$calc_isotopes, mod_col = mod_col_name)
      
      # Remove animation and finalize progress
      shinyjs::runjs("")
      reactive_data$lcms_annotated <- lcms_annotated
      
      # Quick final animation to 66%
      for (p in 40:66) {
        Sys.sleep(0.01)
        shinyjs::runjs(paste0("document.getElementById('progress_bar_fill').style.width = '", p, "%'; document.getElementById('progress_percent').innerText = '", p, "%';"))
      }
      shinyjs::runjs("document.getElementById('task_3').classList.add('completed');")
      
      # Matching
      reactive_data$progress <- 70
      shinyjs::runjs("document.getElementById('progress_bar_fill').style.width = '70%'; document.getElementById('progress_percent').innerText = '70%'; document.getElementById('progress_detail').innerText = 'Matching MALDI peaks to LC-MS ions...'; document.getElementById('task_4').classList.add('completed');")
      
      matches <- MALDIMatches(maldi_mz, lcms_annotated, tolerance = input$tolerance_ppm, output_path = NULL)
      reactive_data$matches <- matches
      
      # Venn
      reactive_data$progress <- 90
      shinyjs::runjs("document.getElementById('progress_bar_fill').style.width = '90%'; document.getElementById('progress_percent').innerText = '90%'; document.getElementById('progress_detail').innerText = 'Creating Venn diagram...'; document.getElementById('task_5').classList.add('completed');")
      
      venn_plot_result <- PlotVenn(maldi_mz, lcms_annotated, matches, output_path = NULL, colors = c("#ff8c00", "#1a1a1a"))
      reactive_data$venn_plot <- venn_plot_result
      
      # Complete
      reactive_data$progress <- 100
      shinyjs::runjs("document.getElementById('progress_bar_fill').style.width = '100%'; document.getElementById('progress_percent').innerText = '100%'; document.getElementById('progress_detail').innerText = 'Complete!';")
      
      output$status_message <- renderText("✓ Processing complete!")
      reactive_data$processing <- FALSE
      Sys.sleep(1)
      shinyjs::removeClass("progress_container", "show")
      
    }, error = function(e) {
      output$status_message <- renderText(paste("❌ Error:", e$message))
      reactive_data$processing <- FALSE
      shinyjs::removeClass("progress_container", "show")
    })
  })
  
  output$stat_maldi <- renderText(if (is.null(reactive_data$maldi)) "—" else length(reactive_data$maldi))
  output$stat_lcms <- renderText(if (is.null(reactive_data$lcms_annotated)) "—" else nrow(reactive_data$lcms_annotated))
  output$stat_matches <- renderText(if (is.null(reactive_data$matches)) "—" else nrow(reactive_data$matches))
  output$stat_percentage <- renderText({
    if (is.null(reactive_data$matches) || is.null(reactive_data$maldi)) return("—")
    pct <- (length(unique(reactive_data$matches$mz_MALDI)) / length(reactive_data$maldi)) * 100
    paste0(round(pct, 2), "%")
  })
  
  output$matching_table <- renderDT({
    if (is.null(reactive_data$matches)) return(NULL)
    datatable(as.data.frame(reactive_data$matches), options = list(scrollX = TRUE, pageLength = 25, server = FALSE), rownames = FALSE)
  })
  
  output$venn_plot <- renderPlot({
    if (is.null(reactive_data$venn_plot)) {
      plot.new()
      text(0.5, 0.5, "No data to display.\nPlease process files first.", cex = 1.5, col = "#ff8c00")
      return()
    }
    grid::grid.draw(reactive_data$venn_plot)
  }, bg = "white")
  
  output$btn_download_matches <- downloadHandler(
    filename = function() paste0("MALDImatch_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      if (is.null(reactive_data$matches)) {
        showNotification("No matches to download. Please process files first.", type = "error")
        return()
      }
      write.csv(as.data.frame(reactive_data$matches), file, row.names = FALSE)
    }
  )
  
  output$btn_download_venn <- downloadHandler(
    filename = function() paste0("MALDImatch_venn_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      if (is.null(reactive_data$venn_plot)) {
        showNotification("No Venn diagram to download. Please process files first.", type = "error")
        return()
      }
      png(file, width = 8, height = 8, units = "in", res = 300)
      grid::grid.draw(reactive_data$venn_plot)
      dev.off()
    }
  )
}

shinyApp(ui, server)
