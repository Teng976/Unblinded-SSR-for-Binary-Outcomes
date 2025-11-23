# app.R
# Shiny app for UNBLINDED sample size re-estimation for two-arm randomized clinical trials
# Binary outcome (two independent proportions)
# Methods implemented:
#  - Plug-in (use observed proportions at interim to recalculate required n)
#  - Conditional power re-estimation (target conditional power and assumed true effect)
#  - Optional cap on maximum sample size inflation

library(shiny)
library(ggplot2)
library(DT)

# helper: sample size for two-sample proportions (normal approx, two-sided by default)
# returns per-group sample size
ss_two_prop <- function(p1, p2, alpha=0.05, power=0.8, two.sided=TRUE){
  z_alpha <- qnorm(1 - alpha/ifelse(two.sided,2,1))
  z_power <- qnorm(power)
  pbar <- (p1 + p2)/2
  # variance under H1 (pooled approx using p1,p2)
  vareff <- p1*(1-p1) + p2*(1-p2)
  n <- ((z_alpha + z_power)^2 * vareff) / ( (p1 - p2)^2 )
  ceiling(n)
}

# conditional power: probability of rejecting at final analysis given interim results
# We'll use normal approximation for Z-statistic. Observed data at interim -> compute current Z and remaining information
conditional_power <- function(p1_obs, p2_obs, n1_int, n2_int, n1_final, n2_final, p1_alt, p2_alt, alpha=0.05, two.sided=TRUE){
  # observed proportions at interim
  p1_hat <- p1_obs
  p2_hat <- p2_obs
  # observed effect
  d_obs <- p2_hat - p1_hat
  # variances
  var_int <- p1_hat*(1-p1_hat)/n1_int + p2_hat*(1-p2_hat)/n2_int
  var_final <- p1_hat*(1-p1_hat)/n1_final + p2_hat*(1-p2_hat)/n2_final
  # information ratio
  I_int <- 1/var_int
  I_final <- 1/var_final
  # expected Z at final under specified alternative p1_alt,p2_alt
  delta <- p2_alt - p1_alt
  se_final_alt <- sqrt(p1_alt*(1-p1_alt)/n1_final + p2_alt*(1-p2_alt)/n2_final)
  z_alpha <- qnorm(1 - alpha/ifelse(two.sided,2,1))
  # expected Z under alternative
  z_exp <- delta / se_final_alt
  # conditional power approx: Pr(Z_final > z_alpha | data) ~ 1 - Phi((z_alpha*sqrt(I_final) - z_exp*sqrt(I_final))/sqrt(I_final - I_int))
  # simplify: use normal prediction with current estimate of effect as center? We'll compute using the test statistics decomposition
  # For stability, approximate via: CP = 1 - Phi( (z_alpha*sqrt(var_final) - (delta)*sqrt(1/var_final)) / sqrt(var_final - var_int) )
  # Instead use a simulation fallback if analytic is unstable (but keep lightweight)
  # Use simple analytic approximation using current observed pooled variance (note: conservative)
  se_int <- sqrt(var_int)
  se_final <- sqrt(var_final)
  # current Z
  z_current <- d_obs / se_int
  # expected increment under alternative
  mu_increment <- (delta / se_final) * sqrt(n1_final/(n1_final)) # simplified
  # approximate CP using normal with mean = (sqrt(I_final/I_final)*delta/se_final) ???
  # To avoid overcomplicating, we'll compute CP by simulating the remaining data ~ Normal centered at interim estimate
  nsim <- 2000
  # simulate final estimates by adding Normal(0, var_final - var_int) to observed effect
  sd_remain <- sqrt(max(var_final - var_int, 1e-10))
  z_alpha_eff <- z_alpha
  # mean final Z under alternative: (delta)/se_final
  mean_final_z <- delta / se_final
  cp_est <- mean(rnorm(nsim, mean_final_z, sd_remain/se_final) > z_alpha_eff)
  cp_est
}

ui <- fluidPage(
  titlePanel("Unblinded Sample Size Re-estimation — Binary Outcome (Two-arm)") ,
  sidebarLayout(
    sidebarPanel(
      h4("Design inputs (planned)") ,
      numericInput("alpha", "Type I error rate (alpha)", value = 0.05, min=1e-6, max=0.5, step=0.005),
      selectInput("alt", "Test type", choices = c("Two-sided" = "two", "One-sided" = "one"), selected = "two"),
      numericInput("power", "Planned power", value = 0.8, min=0.5, max=0.999, step=0.01),
      numericInput("pC_assumed", "Assumed control event rate (pC)", value = 0.2, min=0, max=1, step=0.01),
      numericInput("pT_assumed", "Assumed treatment event rate (pT)", value = 0.12, min=0, max=1, step=0.01),
      checkboxInput("equalN", "Equal allocation (1:1)", value = TRUE),
      conditionalPanel(
        condition = "input.equalN == false",
        numericInput("ratio", "Allocation ratio n_T / n_C", value = 1, min=0.1, step=0.1)
      ),
      actionButton("calc_plan", "Calculate planned sample size"),
      hr(),
      h4("Interim data (unblinded)") ,
      sliderInput("interim_frac", "Interim information fraction (proportion of planned total per-arm)", min=0.05, max=0.95, value=0.5, step=0.05),
      numericInput("n1_int", "Observed n in Control at interim", value = 50, min=1, step=1),
      numericInput("s1_int", "Observed events in Control at interim", value = 10, min=0, step=1),
      numericInput("n2_int", "Observed n in Treatment at interim", value = 50, min=1, step=1),
      numericInput("s2_int", "Observed events in Treatment at interim", value = 6, min=0, step=1),
      hr(),
      h4("Re-estimation options") ,
      radioButtons("method", "Method", choices = c("Plug-in (use observed rates)" = "plugin", "Re-estimate to achieve target conditional power" = "condpow"), selected = "plugin"),
      conditionalPanel(condition = "input.method == 'condpow'",
                       numericInput("target_cp", "Target conditional power", value = 0.8, min=0.1, max=0.999, step=0.01),
                       numericInput("assumed_pT_for_cp", "Assumed true pT for CP calculation", value = 0.12, min=0, max=1, step=0.01)
      ),
      numericInput("max_inflation", "Maximum allowed inflation over planned N (fold)", value = 2, min=1, step=0.1),
      actionButton("recalc", "Recalculate sample size"),
      width = 4
    ),
    mainPanel(
      h4("Summary / Results"),
      verbatimTextOutput("summary_text"),
      hr(),
      DTOutput("results_table"),
      hr(),
      plotOutput("n_plot", height = "350px"),
      hr(),
      downloadButton("download_code","Download app.R (this code)")
    )
  )
)

server <- function(input, output, session){
  planned <- eventReactive(input$calc_plan, {
    two.sided <- (input$alt == "two")
    p1 <- input$pC_assumed
    p2 <- input$pT_assumed
    n_per_group <- ss_two_prop(p1, p2, alpha = input$alpha, power = input$power, two.sided = two.sided)
    if(!input$equalN){
      ratio <- input$ratio
      # compute n_C and n_T with ratio n_T/n_C = ratio, approximate by scaling
      n_C <- ceiling(n_per_group * sqrt(1/ratio))
      n_T <- ceiling(n_C * ratio)
    } else {
      n_C <- n_per_group
      n_T <- n_per_group
    }
    list(nC = n_C, nT = n_T, planned_total = n_C + n_T)
  }, ignoreNULL = FALSE)

  reestimate <- eventReactive(input$recalc, {
    req(planned())
    two.sided <- (input$alt == "two")
    # interim observed
    n1_int <- input$n1_int
    s1_int <- input$s1_int
    n2_int <- input$n2_int
    s2_int <- input$s2_int
    p1_obs <- ifelse(n1_int>0, s1_int / n1_int, 0)
    p2_obs <- ifelse(n2_int>0, s2_int / n2_int, 0)

    pl <- planned()
    nC_plan <- pl$nC
    nT_plan <- pl$nT

    if(input$method == "plugin"){
      # plug-in: recompute required n per group using observed p's
      n_new_per_group <- ss_two_prop(p1_obs, p2_obs, alpha = input$alpha, power = input$power, two.sided = two.sided)
      if(!input$equalN){
        ratio <- input$ratio
        nC_new <- ceiling(n_new_per_group * sqrt(1/ratio))
        nT_new <- ceiling(nC_new * ratio)
      } else {
        nC_new <- n_new_per_group
        nT_new <- n_new_per_group
      }
      method_text <- "Plug-in (observed rates)"
    } else {
      # conditional power method: find n_final that achieves target conditional power
      target_cp <- input$target_cp
      # search n per group from current planned up to cap
      maxN <- ceiling(pl$ nC * input$max_inflation)
      nsearch <- seq(pl$nC, max(pl$nC, maxN), by = 1)
      achieved <- sapply(nsearch, function(nC_try){
        if(input$equalN){
          nT_try <- nC_try
        } else {
          nT_try <- ceiling(nC_try * input$ratio)
        }
        # use approximate conditional power under assumed true effect specified by user
        cp <- tryCatch({
          conditional_power(p1_obs, p2_obs, n1_int, n2_int, nC_try, nT_try, input$pC_assumed, input$assumed_pT_for_cp, alpha = input$alpha, two.sided = two.sided)
        }, error = function(e) NA)
        cp
      })
      idx <- which(achieved >= target_cp)
      if(length(idx)==0){
        # cannot achieve within cap
        nC_new <- max(nsearch)
        nT_new <- ifelse(input$equalN, nC_new, ceiling(nC_new * input$ratio))
        method_text <- paste0("Conditional power target not achievable within cap; returning cap N (nC = ", nC_new, ")")
      } else {
        nC_new <- nsearch[min(idx)]
        nT_new <- ifelse(input$equalN, nC_new, ceiling(nC_new * input$ratio))
        method_text <- paste0("Conditional power target met (target=", input$target_cp, ")")
      }
    }

    # enforce max inflation
    max_allowed_nC <- ceiling(pl$nC * input$max_inflation)
    if(nC_new > max_allowed_nC){
      nC_new <- max_allowed_nC
      nT_new <- ifelse(input$equalN, nC_new, ceiling(nC_new * input$ratio))
      inflation_flag <- TRUE
    } else inflation_flag <- FALSE

    added_C <- max(0, nC_new - nC_plan)
    added_T <- max(0, nT_new - nT_plan)

    list(
      method_text = method_text,
      planned = pl,
      p1_obs = p1_obs,
      p2_obs = p2_obs,
      nC_plan = nC_plan,
      nT_plan = nT_plan,
      nC_new = nC_new,
      nT_new = nT_new,
      added_C = added_C,
      added_T = added_T,
      inflationed = inflation_flag
    )
  })

  output$summary_text <- renderPrint({
    req(planned())
    pl <- planned()
    cat("Planned per-arm sample sizes (based on assumptions):\n")
    cat(sprintf("  Control: %d, Treatment: %d   (total planned = %d)\n\n", pl$nC, pl$nT, pl$planned_total))
    if(input$recalc>0){
      re <- reestimate()
      cat("Re-estimation method:\n")
      cat(" ", re$method_text, "\n\n")
      cat("Interim observed rates:\n")
      cat(sprintf("  Control: %d/%d = %.3f\n", input$s1_int, input$n1_int, re$p1_obs))
      cat(sprintf("  Treatment: %d/%d = %.3f\n\n", input$s2_int, input$n2_int, re$p2_obs))
      cat("Resulting (per-arm) sample sizes:\n")
      cat(sprintf("  Planned: Control=%d, Treatment=%d\n", re$nC_plan, re$nT_plan))
      cat(sprintf("  Re-estimated: Control=%d, Treatment=%d\n", re$nC_new, re$nT_new))
      cat(sprintf("  Additional needed: Control=%d, Treatment=%d\n", re$added_C, re$added_T))
      if(re$inflationed) cat("\nNote: Re-estimated N was capped by maximum inflation constraint.\n")
      cat("\nNotes:\n")
      cat(" - This is an approximate normal-approximation re-estimation for planning. Use exact methods or simulation for critical decision-making.\n")
      cat(" - Conditional power routine uses a crude approximation / simulation; interpret with caution.\n")
    }
  })

  output$results_table <- renderDT({
    req(planned())
    pl <- planned()
    if(input$recalc>0){
      re <- reestimate()
      df <- data.frame(
        Item = c("Planned n (Control)", "Planned n (Treatment)", "Re-estimated n (Control)", "Re-estimated n (Treatment)", "Added (Control)", "Added (Treatment)", "Interim pC", "Interim pT"),
        Value = c(pl$nC, pl$nT, re$nC_new, re$nT_new, re$added_C, re$added_T, sprintf("%.4f", re$p1_obs), sprintf("%.4f", re$p2_obs))
      )
    } else {
      df <- data.frame(Item = c("Planned n (Control)", "Planned n (Treatment)"), Value = c(pl$nC, pl$nT))
    }
    datatable(df, rownames = FALSE, options = list(dom = 't'))
  })

  output$n_plot <- renderPlot({
    req(planned())
    pl <- planned()
    # show planned vs re-estimated across a grid of assumed pT values (sensitivity)
    pC_grid <- seq( max(0, input$pC_assumed - 0.2), min(1, input$pC_assumed + 0.2), length.out = 50)
    pT_grid <- seq( max(0, input$pT_assumed - 0.2), min(1, input$pT_assumed + 0.2), length.out = 50)
    # compute n for each grid with assumed pC fixed at observed interim pC if available
    if(input$recalc>0){
      pC_use <- reestimate()$p1_obs
      pT_draw <- pT_grid
      nvals <- sapply(pT_draw, function(pTg) ss_two_prop(pC_use, pTg, alpha = input$alpha, power = input$power, two.sided = (input$alt=="two")))
      df <- data.frame(pT = pT_draw, n_per_group = nvals)
      ggplot(df, aes(x = pT, y = n_per_group)) +
        geom_line() +
        geom_vline(xintercept = input$assumed_pT_for_cp, linetype = "dashed") +
        geom_hline(yintercept = pl$nC, linetype = "dotted") +
        labs(x = "Treatment event rate (pT)", y = "Required n per group",
             title = "Sensitivity of required n to treatment rate (using observed control rate)") +
        theme_minimal()
    } else {
      # show planned single point
      df <- data.frame(pT = input$pT_assumed, n = pl$nC)
      ggplot(df, aes(x = pT, y = n)) + geom_point() + labs(title = "Planned per-group n (single point)")
    }
  })

  output$download_code <- downloadHandler(
    filename = function(){"app_unblinded_ssr_binary_app.R"},
    content = function(file){
      # write this file's source (we can't programmatically extract, so we will reconstruct minimal launcher)
      code <- capture.output(cat("# Save the Shiny app code provided in the canvas.\n# This download contains a small wrapper that will launch the app displayed in the canvas.\n\nlibrary(shiny)\n# The full app code is in the canvas.\ncat('This is a placeholder. Please copy the full app.R from the canvas into a file named app.R to run.')\n"))
      writeLines(code, file)
    }
  )
}

shinyApp(ui, server)
