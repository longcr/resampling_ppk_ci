library(shiny)
library(data.table)

# -----------------------------
# Ppk function
# -----------------------------
ppk_fun <- function(x, LSL, USL) {
  m  <- mean(x, na.rm = TRUE)
  s  <- sd(x,  na.rm = TRUE)
  if (is.na(s) || s == 0) return(NA_real_)
  min((m - LSL) / (3 * s),
      (USL - m) / (3 * s))
}

# -----------------------------
# BCa interval helper
# -----------------------------
bca_interval <- function(x, boot_vals, stat_fun, LSL, USL, alpha) {
  
  # Original statistic
  theta_hat <- stat_fun(x, LSL, USL)
  
  # Bias correction
  z0 <- qnorm(mean(boot_vals < theta_hat, na.rm = TRUE))
  
  # Jackknife distribution (leave-one-out on original data)
  n <- length(x)
  jack_vals <- numeric(n)
  
  for (i in seq_len(n)) {
    jack_sample <- x[-i]
    jack_vals[i] <- stat_fun(jack_sample, LSL, USL)
  }
  
  # Remove NA jackknife values
  jack_vals <- jack_vals[!is.na(jack_vals)]
  
  # If jackknife fails, return NA and let caller handle fallback
  if (length(jack_vals) < 5) {
    return(c(low = NA, high = NA))
  }
  
  jack_mean <- mean(jack_vals)
  
  # Acceleration term
  num <- sum((jack_mean - jack_vals)^3)
  den <- 6 * (sum((jack_mean - jack_vals)^2))^(3/2)
  
  if (den == 0) {
    a <- 0
  } else {
    a <- num / den
  }
  
  # Adjusted quantiles
  z_low  <- qnorm(alpha / 2)
  z_high <- qnorm(1 - alpha / 2)
  
  adj_low  <- pnorm(z0 + (z0 + z_low)  / (1 - a * (z0 + z_low)))
  adj_high <- pnorm(z0 + (z0 + z_high) / (1 - a * (z0 + z_high)))
  
  # Final BCa interval
  c(
    low  = quantile(boot_vals, adj_low,  na.rm = TRUE),
    high = quantile(boot_vals, adj_high, na.rm = TRUE)
  )
}

# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  titlePanel("Bootstrap Ppk Confidence Interval"),
  
  tabsetPanel(
    tabPanel("Analysis",
             sidebarLayout(
               sidebarPanel(
                 fileInput("file", "Upload data file (CSV)", accept = c(".csv")),
                 uiOutput("col_selector"),
                 numericInput("lsl", "Lower Spec Limit (LSL)", value = 0),
                 numericInput("usl", "Upper Spec Limit (USL)", value = 1),
                 numericInput("n_iter", "Bootstrap iterations",
                              value = 1000, min = 10, step = 10),
                 
                 selectInput("interval_type", "Interval Type",
                             choices = c("Two-sided", "Lower one-sided", "Upper one-sided"),
                             selected = "Two-sided"),
                 
                 numericInput("alpha", "Total alpha (e.g., 0.10 for 90% CI)",
                              value = 0.10, min = 0.001, max = 0.5, step = 0.01),
                 
                 selectInput("ci_method", "CI Method",
                             choices = c("BCa", "Percentile"),
                             selected = "BCa"),
                 
                 actionButton("run", "Run resampling"),
                 br(), br(),
                 downloadButton("download_csv", "Download Bootstrap Results (CSV)")
               ),
               
               mainPanel(
                 plotOutput("ppk_hist")
               )
             )
    ),
    
    tabPanel("About",
             h3("What This App Does"),
             p("This application estimates confidence intervals for the process capability index Ppk using bootstrap resampling."),
             p("You upload a dataset, select a numeric column, and the app repeatedly resamples the data with replacement. For each resample, it computes Ppk. The distribution of these bootstrap Ppk values is then used to estimate confidence intervals."),
             
             h3("Why This App Was Created"),
             p("Calculations for Ppk assume that the data is normally distributed."),
             p("This app was created for the purpose of providing an alternative to Box-Cox transforms on non-normal data prior to calculating Ppk."),
             
             h3("Bootstrap Approach"),
             p("Bootstrap is a nonparametric resampling method. Instead of assuming a distribution for the data, it repeatedly samples from the observed dataset to approximate the sampling distribution of a statistic."),
             tags$ul(
               tags$li("Resample the data with replacement"),
               tags$li("Compute Ppk for each resample"),
               tags$li("Use the empirical distribution of Ppk values to estimate confidence intervals")
             ),
             p("A good introduction to bootstrap methods:"),
             tags$a(
               href = "https://projecteuclid.org/journals/statistical-science/volume-18/issue-2/An-Introduction-to-the-Bootstrap/10.1214/ss/1063994964.full",
               "Efron & Tibshirani (1993) — An Introduction to the Bootstrap"
             ),
             
             h3("Percentile Interval"),
             p("The percentile interval simply takes the α/2 and 1−α/2 quantiles of the bootstrap distribution. It is easy to compute but can be biased if the statistic is skewed or nonlinear."),
             tags$a(
               href = "https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Percentile_method",
               "Percentile Method (Wikipedia)"
             ),
             
             h3("BCa Interval"),
             p("BCa (Bias-Corrected and Accelerated) intervals adjust for both bias and skewness in the bootstrap distribution. They are generally more accurate than percentile intervals."),
             tags$a(
               href = "https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Bias-corrected_and_accelerated_(BCa)_bootstrap",
               "BCa Method (Wikipedia)"
             )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  dat <- reactive({
    req(input$file)
    fread(input$file$datapath)
  })
  
  output$col_selector <- renderUI({
    req(dat())
    num_cols <- names(dat())[sapply(dat(), is.numeric)]
    selectInput("col", "Select numeric column for Ppk",
                choices = num_cols, selected = num_cols[1])
  })
  
  boot_res <- eventReactive(input$run, {
    req(dat(), input$col)
    
    x <- dat()[[input$col]]
    n <- length(x)
    B <- input$n_iter
    
    ppk_vals <- numeric(B)
    
    for (i in seq_len(B)) {
      set.seed(i)  # seed = loop index for reproducibility
      idx <- sample(seq_len(n), size = n, replace = TRUE)
      ppk_vals[i] <- ppk_fun(x[idx], input$lsl, input$usl)
    }
    
    data.table(iteration = seq_len(B), Ppk = ppk_vals)
  })
  
  output$download_csv <- downloadHandler(
    filename = function() paste0("bootstrap_ppk_", Sys.Date(), ".csv"),
    content = function(file) fwrite(boot_res(), file)
  )
  
  output$ppk_hist <- renderPlot({
    dt <- boot_res()
    req(dt)
    
    x <- dat()[[input$col]]
    ppk_vals <- dt$Ppk
    orig_ppk <- ppk_fun(x, input$lsl, input$usl)
    
    alpha  <- input$alpha
    type   <- input$interval_type
    method <- input$ci_method
    
    q_low  <- NA
    q_high <- NA
    
    # Percentile CI
    if (method == "Percentile") {
      if (type == "Two-sided") {
        q_low  <- quantile(ppk_vals, alpha / 2, na.rm = TRUE)
        q_high <- quantile(ppk_vals, 1 - alpha / 2, na.rm = TRUE)
      } else if (type == "Lower one-sided") {
        q_low  <- quantile(ppk_vals, alpha, na.rm = TRUE)
      } else { # Upper one-sided
        q_high <- quantile(ppk_vals, 1 - alpha, na.rm = TRUE)
      }
    }
    
    # BCa CI
    if (method == "BCa") {
      
      # Two-sided BCa
      if (type == "Two-sided") {
        ci <- bca_interval(
          x = x,
          boot_vals = ppk_vals,
          stat_fun = ppk_fun,
          LSL = input$lsl,
          USL = input$usl,
          alpha = alpha
        )
        q_low  <- ci["low"]
        q_high <- ci["high"]
      }
      
      # Lower one-sided BCa
      else if (type == "Lower one-sided") {
        ci <- bca_interval(
          x = x,
          boot_vals = ppk_vals,
          stat_fun = ppk_fun,
          LSL = input$lsl,
          USL = input$usl,
          alpha = alpha * 2
        )
        q_low  <- ci["low"]
        q_high <- NA
      }
      
      # Upper one-sided BCa
      else {
        ci <- bca_interval(
          x = x,
          boot_vals = ppk_vals,
          stat_fun = ppk_fun,
          LSL = input$lsl,
          USL = input$usl,
          alpha = alpha * 2
        )
        q_low  <- NA
        q_high <- ci["high"]
      }
      
      # Fallback if BCa fails
      if (is.na(q_low) && is.na(q_high)) {
        warning("BCa interval could not be computed; falling back to percentile.")
        if (type == "Two-sided") {
          q_low  <- quantile(ppk_vals, alpha / 2, na.rm = TRUE)
          q_high <- quantile(ppk_vals, 1 - alpha / 2, na.rm = TRUE)
        } else if (type == "Lower one-sided") {
          q_low  <- quantile(ppk_vals, alpha, na.rm = TRUE)
        } else {
          q_high <- quantile(ppk_vals, 1 - alpha, na.rm = TRUE)
        }
      }
    }
    
    # Plot
    hist(ppk_vals,
         breaks = "FD",
         main = paste("Bootstrap Ppk Distribution (", method, " CI)", sep = ""),
         xlab = "Ppk",
         col = "lightgray",
         border = "white")
    
    med <- median(ppk_vals, na.rm = TRUE)
    abline(v = med, col = "blue", lwd = 2)
    
    if (!is.na(q_low))  abline(v = q_low,  col = "red", lwd = 2, lty = 2)
    if (!is.na(q_high)) abline(v = q_high, col = "red", lwd = 2, lty = 2)
    
    legend_items <- c(paste0("Median = ", round(med, 3)))
    legend_cols  <- c("blue")
    legend_lty   <- c(1)
    legend_lwd   <- c(2)
    
    if (!is.na(q_low)) {
      legend_items <- c(legend_items,
                        paste0("Lower CI = ", round(q_low, 3)))
      legend_cols  <- c(legend_cols, "red")
      legend_lty   <- c(legend_lty, 2)
      legend_lwd   <- c(legend_lwd, 2)
    }
    
    if (!is.na(q_high)) {
      legend_items <- c(legend_items,
                        paste0("Upper CI = ", round(q_high, 3)))
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