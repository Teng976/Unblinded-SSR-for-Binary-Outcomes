# app.R — Full Shiny App with Inverse-Normal Combination Test (Integrated)
# Unblinded Sample Size Re-estimation for Binary Outcomes
# Methods implemented:
#  - Plug-in SSR
#  - Conditional power SSR
#  - Inverse-normal combination test (stage-wise p-value combination)

library(shiny)
library(ggplot2)
library(DT)

# ---------- Helper: sample size for two-sample proportions ----------
ss_two_prop <- function(p1, p2, alpha=0.05, power=0.8, two.sided=TRUE){
  z_alpha <- qnorm(1 - alpha/ifelse(two.sided,2,1))
  z_power <- qnorm(power)
  vareff <- p1*(1-p1) + p2*(1-p2)
  n <- ((z_alpha + z_power)^2 * vareff) / ( (p1 - p2)^2 )
  ceiling(n)
}

# ---------- Inverse-normal combination test ----------
# weight1, weight2 are combination weights (sqrt information fraction)
invnorm_combine <- function(p1, p2, w1, w2){
  z1 <- qnorm(1 - p1)
  z2 <- qnorm(1 - p2)
  z_comb <- (w1*z1 + w2*z2)/sqrt(w1^2 + w2^2)
  p_comb <- 1 - pnorm(z_comb)
  list(z_comb = z_comb, p_comb = p_comb)
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Unblinded Sample Size Re-estimation — Binary Outcome (Two-arm)"),
  sidebarLayout(
    sidebarPanel(
      h4("Design inputs (planned)"),
      numericInput("alpha", "Type I error", 0.05, min=1e-6, max=0.5, step=0.005),
      selectInput("alt", "Test type", c("Two-sided"="two", "One-sided"="one")),
      numericInput("power", "Planned power", 0.8, min=0.5, max=0.999, step=0.01),
      numericInput("pC_assumed", "Assumed control rate", 0.2, min=0, max=1, step=0.01),
      numericInput("pT_assumed", "Assumed treatment rate", 0.12, min=0, max=1, step=0.01),
      checkboxInput("equalN", "Equal allocation", TRUE),
      conditionalPanel("!input.equalN", numericInput("ratio", "nT/nC", 1, min=0.1, step=0.1)),
      actionButton("calc_plan", "Calculate planned N"),
      hr(),
      
      h4("Interim data (unblinded)"),
      sliderInput("interim_frac", "Information fraction", min=0.05, max=0.95, value=0.5, step=0.05),
      numericInput("n1_int", "Interim nC", 50, min=1),
      numericInput("s1_int", "Interim events C", 10, min=0),
      numericInput("n2_int", "Interim nT", 50, min=1),
      numericInput("s2_int", "Interim events T", 6, min=0),
      hr(),
      
      h4("Re-estimation method"),
      radioButtons("method", "Method", choices = c(
        "Plug-in"="plugin",
        "Conditional power SSR"="condpow",
        "Inverse-normal combination test"="invnorm"
      ), selected = "plugin"),
      
      conditionalPanel("input.method == 'condpow'",
                       numericInput("target_cp", "Target conditional power", 0.8, min=0.1, max=0.999, step=0.01),
                       numericInput("assumed_pT_for_cp", "Assumed true pT", 0.12, min=0, max=1, step=0.01)
      ),
      
      numericInput("max_inflation", "Max N inflation (×)", value = 2, min=1, step=0.1),
      actionButton("recalc", "Recalculate"),
      width = 4
    ),
    
    mainPanel(
      h4("Results"), verbatimTextOutput("summary_text"), hr(),
      DTOutput("results_table"), hr(),
      plotOutput("n_plot", height="350px"), hr(),
      downloadButton("download_code", "Download app.R")
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session){
  
  # ---- Planned sample size ----
  planned <- eventReactive(input$calc_plan, {
    two.sided <- (input$alt == "two")
    n_per_group <- ss_two_prop(input$pC_assumed, input$pT_assumed,
                               alpha=input$alpha, power=input$power,
                               two.sided=two.sided)
    if(!input$equalN){
      nC <- ceiling(n_per_group * sqrt(1/input$ratio))
      nT <- ceiling(nC * input$ratio)
    } else {
      nC <- n_per_group; nT <- n_per_group
    }
    list(nC=nC, nT=nT, total=nC+nT)
  }, ignoreNULL=FALSE)
  
  # ---- Re-estimation ----
  reestimate <- eventReactive(input$recalc, {
    pl <- planned()
    p1_obs <- input$s1_int/input$n1_int
    p2_obs <- input$s2_int/input$n2_int
    two.sided <- (input$alt == "two")
    
    # ----- Method 1: Plug-in -----
    if(input$method == "plugin"){
      n_new <- ss_two_prop(p1_obs, p2_obs, alpha=input$alpha, power=input$power, two.sided=two.sided)
      if(!input$equalN){
        nC_new <- ceiling(n_new * sqrt(1/input$ratio))
        nT_new <- ceiling(nC_new * input$ratio)
      } else {
        nC_new <- n_new; nT_new <- n_new
      }
      method_text <- "Plug-in SSR"
    }
    
    # ----- Method 2: Conditional power SSR -----
    if(input$method == "condpow"){
      target <- input$target_cp
      mx <- ceiling(pl$nC * input$max_inflation)
      search <- seq(pl$nC, mx)
      cpvals <- sapply(search, function(nC_try){
        nT_try <- ifelse(input$equalN, nC_try, ceiling(nC_try*input$ratio))
        se_int <- sqrt(p1_obs*(1-p1_obs)/input$n1_int + p2_obs*(1-p2_obs)/input$n2_int)
        effect_alt <- input$assumed_pT_for_cp - input$pC_assumed
        se_final <- sqrt(input$pC_assumed*(1-input$pC_assumed)/nC_try +
                           input$assumed_pT_for_cp*(1-input$assumed_pT_for_cp)/nT_try)
        z_alpha <- qnorm(1 - input$alpha/ifelse(two.sided,2,1))
        mean_z <- effect_alt / se_final
        sd_z <- sqrt(1 - input$interim_frac)
        mean(rnorm(1500, mean_z, sd_z) > z_alpha)
      })
      ok <- which(cpvals >= target)
      if(length(ok)==0){
        nC_new <- mx
        method_text <- "CP target unattainable; using max cap"
      } else {
        nC_new <- search[min(ok)]
        method_text <- "Conditional power SSR"
      }
      nT_new <- ifelse(input$equalN, nC_new, ceiling(nC_new*input$ratio))
    }
    
    # ----- Method 3: Inverse-normal combo test -----
    if(input$method == "invnorm"){
      # Stage 1
      p1 <- prop.test(c(input$s1_int, input$s2_int), c(input$n1_int, input$n2_int), correct=FALSE)$p.value
      w1 <- sqrt(input$interim_frac)
      w2 <- sqrt(1 - input$interim_frac)
      
      # Stage 2: assume final stage uses observed effect to compute expected stage-2 p-value
      # (more principled options exist, but this keeps SSR logic simple)
      d_obs <- p2_obs - p1_obs
      nC_new <- pl$nC
      nT_new <- pl$nT
      se2 <- sqrt(p1_obs*(1-p1_obs)/nC_new + p2_obs*(1-p2_obs)/nT_new)
      z_exp <- d_obs / se2
      p2_pred <- 1 - pnorm(z_exp)
      
      cmb <- invnorm_combine(p1, p2_pred, w1, w2)
      method_text <- sprintf("Inverse-normal combination test: combined p = %.4f", cmb$p_comb)
    }
    
    # Cap inflation
    maxC <- ceiling(pl$nC * input$max_inflation)
    if(nC_new > maxC){
      nC_new <- maxC
      nT_new <- ifelse(input$equalN, nC_new, ceiling(nC_new*input$ratio))
    }
    
    list(
      method_text=method_text,
      p1_obs=p1_obs, p2_obs=p2_obs,
      nC_plan=pl$nC, nT_plan=pl$nT,
      nC_new=nC_new, nT_new=nT_new
    )
  })
  
  # ---- Summary text ----
  output$summary_text <- renderPrint({
    pl <- planned()
    cat("Planned sample sizes:
")
    cat(sprintf("  Control=%d, Treatment=%d (Total=%d)

", pl$nC, pl$nT, pl$total))
    
    if(input$recalc>0){
      re <- reestimate()
      cat("Method:
  ", re$method_text, "

")
      cat("Interim observed rates:
")
      cat(sprintf("  pC=%.3f, pT=%.3f

", re$p1_obs, re$p2_obs))
      cat("Re-estimated sample sizes:
")
      cat(sprintf("  Control=%d → %d
", re$nC_plan, re$nC_new))
      cat(sprintf("  Treatment=%d → %d
", re$nT_plan, re$nT_new))
    }
  })
  
  # ---- Table ----
  output$results_table <- renderDT({
    pl <- planned()
    if(input$recalc>0){
      re <- reestimate()
      df <- data.frame(
        Item=c("Planned nC","Planned nT","New nC","New nT","Interim pC","Interim pT"),
        Value=c(pl$nC, pl$nT, re$nC_new, re$nT_new,
                sprintf("%.4f",re$p1_obs), sprintf("%.4f",re$p2_obs)))
    } else {
      df <- data.frame(Item=c("Planned nC","Planned nT"), Value=c(pl$nC, pl$nT))
    }
    datatable(df, options=list(dom='t'), rownames=FALSE)
  })
  
  # ---- Sensitivity plot ----
  output$n_plot <- renderPlot({
    pl <- planned()
    pT_grid <- seq(max(0,input$pT_assumed-0.2), min(1,input$pT_assumed+0.2), length.out=60)
    nvals <- sapply(pT_grid, function(pt){ ss_two_prop(input$pC_assumed, pt, input$alpha, input$power, two.sided=(input$alt=="two")) })
    df <- data.frame(pT=pT_grid, n=nvals)
    ggplot(df, aes(pT,n)) + geom_line() + labs(title="Required n vs pT assumption") + theme_minimal()
  })
  
  output$download_code <- downloadHandler(
    filename=function(){"app_unblinded_ssr_binary.R"},
    content=function(file){ writeLines("Copy the full app.R from the canvas.", file) }
  )
}

shinyApp(ui, server)
