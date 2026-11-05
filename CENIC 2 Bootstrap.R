
# Calculate empirical power for all four models

# Run bootstrap power calculations for GLM and NB models

power_results <- data.frame()

for (j in sample_sizes) {
  
  message("Running sample size: ", j)
  
  nA <- floor(j / 2)
  nB <- j - nA
  
  # Store p-values for each model
  pvals_glm <- rep(NA_real_, n_boot)
  pvals_any1 <- rep(NA_real_, n_boot)
  pvals_nb <- rep(NA_real_, n_boot)
  
  for (i in seq_len(n_boot)) {
    
    # Bootstrap subjects from each treatment group
    sampleA <- bootstrap_subjects(nb7dayImpA, nA)
    sampleB <- bootstrap_subjects(nb7dayImpB, nB)
    
    boot_data <- bind_rows(sampleA, sampleB)
    
    # Logistic regression: 7-day abstinence
    
    fit_glm <- tryCatch(
      glm(
        abst_7daynoNA ~ Treatment + site,
        data = boot_data,
        family = binomial(link = "logit")
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_glm)) {
      
      sm <- coef(summary(fit_glm))
      
      if ("TreatmentB" %in% rownames(sm)) {
        pvals_glm[i] <- sm["TreatmentB", "Pr(>|z|)"]
      }
    }
    
    
    # Logistic regression: any cigarette-free days
    
    fit_any1 <- tryCatch(
      glm(
        any1day ~ Treatment + site,
        data = boot_data,
        family = binomial(link = "logit")
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_any1)) {
      
      sm <- coef(summary(fit_any1))
      
      if ("TreatmentB" %in% rownames(sm)) {
        pvals_any1[i] <- sm["TreatmentB", "Pr(>|z|)"]
      }
    }
    
    # Negative binomial regression: total cigarette-free days
    
    fit_nb <- tryCatch(
      MASS::glm.nb(
        total_abstinent_days ~ Treatment + site,
        data = boot_data
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_nb)) {
      
      sm <- coef(summary(fit_nb))
      
      if ("TreatmentB" %in% rownames(sm)) {
        pvals_nb[i] <- sm["TreatmentB", "Pr(>|z|)"]
      }
    }
  }
  
  # Calculate power and count successful model fits
  
  result <- data.frame(
    SampleSize = j,
    
    glmpwr = mean(pvals_glm < alpha, na.rm = TRUE),
    
    any1glmpwr = mean(pvals_any1 < alpha, na.rm = TRUE),
    
    nbpwr = mean(pvals_nb < alpha, na.rm = TRUE),
    
    GLM7_Successful = sum(!is.na(pvals_glm)),
    
    Any1_Successful = sum(!is.na(pvals_any1)),
    
    NB_Successful = sum(!is.na(pvals_nb))
  )
  
  power_results <- bind_rows(power_results, result)
  
  # Display results after each sample size
  print(result)
}

# Final GLM and NB power results
power_results

# Event indicator
    abst_7daynoNA = as.integer(!is.na(first(event_day))),
    
    .groups = "drop"
  )

# 7. Cox proportional hazards bootstrap power calculation
# Define sample sizes
sample_sizes <- c(150, 200, 250, 300, 350, 400, 438)

coxpwr <- numeric(length(sample_sizes))

# Separate treatment groups
tbcdataA <- cenic_cox %>%
  filter(Treatment == "A")

tbcdataB <- cenic_cox %>%
  filter(Treatment == "B")


# Run bootstrap power calculations
for (j in seq_along(sample_sizes)) {
  
  n <- sample_sizes[j]
  
  message("Running sample size: ", n)
  
  coxp <- rep(NA_real_, n_boot)
  
  for (i in seq_len(n_boot)) {
    
    # Equal allocation between A and B
    indicesA <- sample(
      seq_len(nrow(tbcdataA)),
      floor(n / 2),
      replace = TRUE
    )
    
    indicesB <- sample(
      seq_len(nrow(tbcdataB)),
      ceiling(n / 2),
      replace = TRUE
    )
    
    data <- bind_rows(
      tbcdataA[indicesA, ],
      tbcdataB[indicesB, ]
    )
    
    
    # Cox proportional hazards model
    fit <- tryCatch(
      coxph(
        Surv(days_from_rand, abst_7daynoNA) ~
          Treatment + strata(site),
        data = data
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit)) {
      
      coxp[i] <- tryCatch(
        summary(fit)$logtest["pvalue"],
        error = function(e) NA_real_
      )
    }
  }
  
  
  # Calculate empirical Cox power
  coxpwr[j] <- mean(coxp < 0.05, na.rm = TRUE)
  
  message("Power = ", round(coxpwr[j], 3))
}


# Final Cox power results
cox_power_results <- data.frame(
  SampleSize = sample_sizes,
  CoxPower = coxpwr
)

cox_power_results

#AG Model Power estimation

sample_sizes <- c(150, 200, 250, 300, 350, 400, 438)
n_boot <- 1000
alpha <- 0.05

AG_power_results <- data.frame(
  SampleSize = sample_sizes,
  Power = NA_real_,
  Successful_Fits = NA_integer_
)

for (j in seq_along(sample_sizes)) {
  
  nA <- sample_sizes[j]
  nB <- sample_sizes[j]
  pvals_AG <- rep(NA_real_, n_boot)
  
  for (i in 1:n_boot) {
    
    # Bootstrap participants within each treatment group
    sampleAGA <- bootstrap_AG(AGdataA, nA, "A", i)
    sampleAGB <- bootstrap_AG(AGdataB, nB, "B", i)
    
    AG_boot_data <- bind_rows(sampleAGA, sampleAGB)
    
    # Fit Andersen-Gill model
    fit_AG <- tryCatch({
      coxph(
        Surv(start, stop, event) ~
          Treatment + strata(site) + cluster(studyid),
        data = AG_boot_data
      )
    }, error = function(e) NULL)
    
    # Extract treatment p-value
    if (!is.null(fit_AG)) {
      sm <- summary(fit_AG)$coefficients
      
      if ("TreatmentB" %in% rownames(sm)) {
        pvals_AG[i] <- sm["TreatmentB", "Pr(>|z|)"]
      }
    }
  }
  
  # Calculate empirical power
  AG_power_results$Power[j] <- mean(pvals_AG < alpha, na.rm = TRUE)
  AG_power_results$Successful_Fits[j] <- sum(!is.na(pvals_AG))
  
  cat("Completed sample size:", sample_sizes[j],
      "| Power:", round(AG_power_results$Power[j], 4),
      "| Successful fits:", AG_power_results$Successful_Fits[j], "\n")
}

print(AG_power_results)

# 8. Reshape GLM and NB results for plotting


longdata <- power_results %>%
  pivot_longer(
    cols = -SampleSize,
    names_to = "MethodCode",
    values_to = "Power"
  )

method_labels <- c(
  glmpwr = "Logistic Regression, 7-Day PPA",
  any1glmpwr = "Logistic Regression, Any Cigarette-Free Days",
  nbpwr = "Negative Binomial, Cigarette-Free Days",
  AGpwr = "Andersen-Gill Model, Time to 7-Day Abstinence"
)

longdata$Method <- method_labels[longdata$MethodCode]
