library(shiny)
library(data.table)

# Function to compute Ppk
ppk_fun <- function(x, LSL, USL) {
  m  <- mean(x, na.rm = TRUE)
  s  <- sd(x,  na.rm = TRUE)
  if (is.na(s) || s == 0) return(NA_real_)
  ppk <- min((m - LSL) / (3 * s),
             (USL - m) / (3 * s))
  return(ppk)
}

ui <- fluidPage(
  titlePanel("Bootstrap Ppk Confidence Interval"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload data file (CSV)", accept = c(".csv")),
      uiOutput("col_selector"),
      numericInput("lsl", "Lower Spec Limit (LSL)", value = 0),
      numericInput("usl", "Upper Spec Limit (USL)", value = 1),
      numericInput("n_iter", "Number of bootstrap iterations",
                   value = 1000, min = 10, step = 10),
      
      selectInput("interval_type", "Interval Type",
                  choices = c("Two-sided", "Lower one-sided", "Upper one-sided"),
                  selected = "Two-sided"),
      
      numericInput("alpha", "Total alpha (e.g., 0.10 for 90% CI)",
                   value = 0.10, min = 0.001, max = 0.5, step = 0.01),
      
      actionButton("run", "Run resampling"),
      br(), br(),
      downloadButton("download_csv", "Download Bootstrap Results (CSV)")
    ),
    
    mainPanel(
      plotOutput("ppk_hist")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive: read data
  dat <- reactive({
    req(input$file)
    fread(input$file$datapath)
  })
  
  # Reactive UI: choose numeric column
  output$col_selector <- renderUI({
    req(dat())
    num_cols <- names(dat())[sapply(dat(), is.numeric)]
    selectInput("col", "Select numeric column for Ppk",
                choices = num_cols, selected = num_cols[1])
  })
  
  # Run bootstrap when button clicked
  boot_res <- eventReactive(input$run, {
    req(dat(), input$col, input$lsl, input$usl, input$n_iter)
    
    x <- dat()[[input$col]]
    n <- length(x)
    B <- input$n_iter
    
    ppk_vals <- numeric(B)
    
    for (i in seq_len(B)) {
      set.seed(i)  # seed = loop index for reproducibility
      idx <- sample(seq_len(n), size = n, replace = TRUE)
      x_b <- x[idx]
      ppk_vals[i] <- ppk_fun(x_b, input$lsl, input$usl)
    }
    
    data.table(iteration = seq_len(B), Ppk = ppk_vals)
  })
  
  # Download handler for CSV
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("bootstrap_ppk_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fwrite(boot_res(), file)
    }
  )
  
  # Histogram with interval selection
  output$ppk_hist <- renderPlot({
    dt <- boot_res()
    req(dt)
    
    ppk_vals <- dt$Ppk
    p_med <- median(ppk_vals, na.rm = TRUE)
    
    alpha <- input$alpha
    type  <- input$interval_type
    
    # Determine quantiles based on interval type
    if (type == "Two-sided") {
      q_low  <- quantile(ppk_vals, alpha / 2, na.rm = TRUE)
      q_high <- quantile(ppk_vals, 1 - alpha / 2, na.rm = TRUE)
      label_low  <- paste0(round(alpha/2 * 100, 1), "%")
      label_high <- paste0(round((1 - alpha/2) * 100, 1), "%")
    } else if (type == "Lower one-sided") {
      q_low  <- quantile(ppk_vals, alpha, na.rm = TRUE)
      q_high <- NA
      label_low <- paste0(round(alpha * 100, 1), "%")
      label_high <- NULL
    } else { # Upper one-sided
      q_low  <- NA
      q_high <- quantile(ppk_vals, 1 - alpha, na.rm = TRUE)
      label_low <- NULL
      label_high <- paste0(round((1 - alpha) * 100, 1), "%")
    }
    
    hist(ppk_vals,
         breaks = "FD",
         main = paste("Bootstrap distribution of Ppk\n(", type, " interval)", sep = ""),
         xlab = "Ppk",
         col = "lightgray",
         border = "white")
    
    abline(v = p_med, col = "blue", lwd = 2)
    
    if (!is.na(q_low))  abline(v = q_low,  col = "red", lwd = 2, lty = 2)
    if (!is.na(q_high)) abline(v = q_high, col = "red", lwd = 2, lty = 2)
    
    legend_items <- c(paste0("Median = ", round(p_med, 3)))
    legend_cols  <- c("blue")
    legend_lty   <- c(1)
    legend_lwd   <- c(2)
    
    if (!is.na(q_low)) {
      legend_items <- c(legend_items,
                        paste0(label_low, " = ", round(q_low, 3)))
      legend_cols  <- c(legend_cols, "red")
      legend_lty   <- c(legend_lty, 2)
      legend_lwd   <- c(legend_lwd, 2)
    }
    
    if (!is.na(q_high)) {
      legend_items <- c(legend_items,
                        paste0(label_high, " = ", round(q_high, 3)))
      legend_cols  <- c(legend_cols, "red")
      legend_lty   <- c(legend_lty, 2)
      legend_lwd   <- c(legend_lwd, 2)
    }
    
    legend("topright",
           legend = legend_items,
           col = legend_cols,
           lty = legend_lty,
           lwd = legend_lwd,
           bty = "n")
  })
}

shinyApp(ui, server)