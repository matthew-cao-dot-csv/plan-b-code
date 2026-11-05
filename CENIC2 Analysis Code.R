#' Quick road map of what I'm Doing Here:
#' 1. Load in the data.
#' 2. Prepare dataset, create survival/ nb/ count variables as needed.
#' 3. Define models for those variables.
#' 4. Use bootstrap-resampling to estimate power.

#1. Load in the data
cenic2_long <- read.csv("~/data and codebook/CENIC2 data/matthew_cenic2_daily_2026_11_02.csv")
cenic2_wide <- read.csv("~/data and codebook/CENIC2 data/matthew_cenic2_wide_2026_11_02.csv")
set.seed(1738)

#2. Prepare Dataset
# Variable names
id_var <- "studyid"
time_var <- "days_from_rand"
cpd_var <- "TotalCPD_ivr"
treatment_var <- "Treatment"

visit_id_var <- names(cenic2_wide)[1]
cutoff_var <- names(cenic2_wide)[4]

# Treatment groups
treatment_levels <- c("B", "A")

#Include site and treatment
cenic2_long <- cenic2_long %>%
  left_join(
    cenic2_wide %>%
      dplyr::select(studyid, site, Treatment),
    by = "studyid"
  )

#Truncate dataset at week 12
cenic_truncated<-data.frame()
cenic_truncated <- cenic2_long %>%
  left_join(
    cenic2_wide %>%
      dplyr::select(studyid, cutoff = 4),
    by = "studyid"
  ) %>%
  filter(days_from_rand < cutoff) %>%
  dplyr::select(-cutoff)


#Create survival variables 

cenic_truncated$abst_7day<-NA

for (i in 7:nrow(cenic_truncated)) {
  ifelse(cenic_truncated$days_from_rand[i]>=7&&0==sum(cenic_truncated[(i-6):i,"TotalCPD_ivr"]),1,0)->cenic_truncated$abst_7day[i]
}

cenic_truncated$abst_3day<-NA

for (i in 3:nrow(cenic_truncated)) {
  ifelse(cenic_truncated$days_from_rand[i]>=3&&0==sum(cenic_truncated[(i-2):i,"TotalCPD_ivr"]),1,0)->cenic_truncated$abst_3day[i]
}

cenic_truncated$abst_7day->cenic_truncated$abst_7daynoNA
cenic_truncated[is.na(cenic_truncated$abst_7day),"abst_7daynoNA"]<-0

cenic_truncated$abst_3day->cenic_truncated$abst_3daynoNA
cenic_truncated[is.na(cenic_truncated$abst_7day),"abst_3daynoNA"]<-0

cenic_truncated$cfd<-ifelse(cenic_truncated$TotalCPD_ivr>0|is.na(cenic_truncated$TotalCPD_ivr),0,1)

#Andersen-Gills variable

cenic_episode_starts <- cenic_truncated %>%
  arrange(studyid, days_from_rand) %>%
  group_by(studyid) %>%
  mutate(
    abstinent = as.integer(
      !is.na(TotalCPD_ivr) & TotalCPD_ivr == 0
    ),
    
    new_run = abstinent == 1 &
      (
        lag(abstinent, default = 0) == 0 |
          days_from_rand != lag(days_from_rand) + 1
      ),
    
    run_id = cumsum(new_run)
  ) %>%
  filter(abstinent == 1) %>%
  group_by(studyid, run_id) %>%
  mutate(
    episode_start = min(days_from_rand),
    episode_length = n()
  ) %>%
  filter(episode_length >= 7) %>%
  summarise(
    event_day = first(episode_start),
    .groups = "drop"
  ) %>%
  distinct(studyid, event_day)

cenic_AG <- cenic_truncated %>%
  arrange(studyid, days_from_rand) %>%
  group_by(studyid) %>%
  mutate(
    start = days_from_rand,
    stop = days_from_rand + 1
  ) %>%
  ungroup() %>%
  left_join(
    cenic_episode_starts %>%
      mutate(event = 1L),
    by = c("studyid" = "studyid",
           "days_from_rand" = "event_day")
  ) %>%
  mutate(
    event = ifelse(is.na(event), 0L, event)
  ) %>%
  arrange(studyid, start)

saveRDS(cenic_truncated, "cenic_truncated.csv")

nb7dayImp_cenic <- cenic_truncated %>%
  group_by(studyid) %>%
  mutate(total_abstinent_days = sum(cfd, na.rm = TRUE)) %>%
  mutate(any1day = ifelse(total_abstinent_days>0,1,0)) %>% 
  mutate(total_PPA = sum(abst_7daynoNA, na.rm = TRUE)) %>%
  mutate(abst_7dayPPA = ifelse(total_PPA>0,1,0)) %>%
  slice_max(days_from_rand,with_ties = FALSE) %>% 
  ungroup()

cenic_cox <- cenic_truncated %>%
  arrange(studyid, days_from_rand) %>%
  group_by(studyid) %>%
  mutate(
    event_day = ifelse(
      any(abst_7daynoNA == 1, na.rm = TRUE),
      min(days_from_rand[abst_7daynoNA == 1], na.rm = TRUE),
      NA_real_
    )
  ) %>%
  summarise(
    Treatment = first(Treatment),
    site = first(site),
    
    # First qualifying abstinence day
    days_from_rand = ifelse(
      !is.na(first(event_day)),
      first(event_day),
      max(days_from_rand)
    ),
    
    # Event indicator
    abst_7daynoNA = as.integer(!is.na(first(event_day))),
    
    .groups = "drop"
  )


# 3.1 Andersen-Gill Model
# Recurrent time-to-event analysis

AG7dayImp <- cenic_AG

AG7_model <- coxph(
  Surv(start, stop, event) ~
    Treatment + strata(site) + cluster(studyid),
  data = AG7dayImp
)

summary(AG7_model)

# 3.2 Logistic Regression: 7 Day Point-prevalent abstinence

glm7_model <- glm(
  abst_7daynoNA ~ Treatment + site,
  data = nb7dayImp_cenic,
  family = binomial(link = "logit")
)

summary(glm7_model)


# 3.3 Logistic Regression: Any Cigarette-Free Days

any1_model <- glm(
  any1day ~ Treatment + site,
  data = nb7dayImp_cenic,
  family = binomial(link = "logit")
)

summary(any1_model)


# 3.4 Negative Binomial Model:
# Number of Cigarette-Free Days

nb_model <- MASS::glm.nb(
  total_abstinent_days ~ Treatment + site,
  data = nb7dayImp_cenic
)

summary(nb_model)

# 5. Prepare treatment groups and reference levels

# Set Treatment A as the reference group
nb7dayImp_cenic <- nb7dayImp_cenic %>%
  mutate(
    Treatment = factor(Treatment, levels = c("A", "B")),
    site = factor(site)
  )

cenic_AG <- cenic_AG %>%
  mutate(
    Treatment = factor(Treatment, levels = c("A", "B")),
    site = factor(site)
  )

# Separate subject-level data by treatment group
nb7dayImpA <- nb7dayImp_cenic %>%
  filter(Treatment == "A")

nb7dayImpB <- nb7dayImp_cenic %>%
  filter(Treatment == "B")

# Separate Andersen-Gill data by treatment group
AGdataA <- cenic_AG %>%
  filter(Treatment == "A")

AGdataB <- cenic_AG %>%
  filter(Treatment == "B")

# 6. Define bootstrap sampling functions

# Subject-level bootstrap for GLM and negative binomial models
bootstrap_subjects <- function(data, n) {
  
  indices <- sample(
    seq_len(nrow(data)),
    size = n,
    replace = TRUE
  )
  
  data[indices, ]
}


# Subject-level bootstrap for Andersen-Gill data
# Preserve all observations belonging to each sampled subject
# and assign a new ID to duplicated subjects.

bootstrap_AG <- function(data, n, arm, iteration) {
  
  subject_ids <- unique(data$studyid)
  
  sampled_ids <- sample(
    subject_ids,
    size = n,
    replace = TRUE
  )
  
  sampled_data <- lapply(seq_along(sampled_ids), function(k) {
    
    temp <- data %>%
      filter(studyid == sampled_ids[k])
    
    temp$studyid <- paste0(
      arm, "_boot_", iteration, "_", k
    )
    
    temp
  })
  
  bind_rows(sampled_data)
}



# 7. Define sample sizes and bootstrap settings


sample_sizes <- c(seq(150, 400, by = 50), 438)

n_boot <- 1000
alpha <- 0.05


# Sample subjects from each treatment arm
sampleA <- bootstrap_subjects(
  nb7dayImpA, nA
)

sampleB <- bootstrap_subjects(
  nb7dayImpB, nB
)

boot_data <- bind_rows(sampleA, sampleB)


#Model 1: Logistic regression for 7-day abstinence

fit_glm <- tryCatch({
  
  glm(
    abst_7daynoNA ~ Treatment + site,
    data = boot_data,
    family = binomial(link = "logit")
  )
  
}, error = function(e) NULL)

if (!is.null(fit_glm)) {
  
  sm <- coef(summary(fit_glm))
  
  if ("TreatmentB" %in% rownames(sm)) {
    pvals_glm[i] <- sm["TreatmentB", "Pr(>|z|)"]
  }
}


# Model 2: Logistic regression for any cigarette-free days

fit_any1 <- tryCatch({
  
  glm(
    any1day ~ Treatment + site,
    data = boot_data,
    family = binomial(link = "logit")
  )
  
}, error = function(e) NULL)

if (!is.null(fit_any1)) {
  
  sm <- coef(summary(fit_any1))
  
  if ("TreatmentB" %in% rownames(sm)) {
    pvals_any1[i] <- sm["TreatmentB", "Pr(>|z|)"]
  }
}


# Model 3: Negative binomial regression
# Number of cigarette-free days
fit_nb <- tryCatch({
  
  MASS::glm.nb(
    total_abstinent_days ~ Treatment + site,
    data = boot_data
  )
  
}, error = function(e) NULL)

if (!is.null(fit_nb)) {
  
  sm <- coef(summary(fit_nb))
  
  if ("TreatmentB" %in% rownames(sm)) {
    pvals_nb[i] <- sm["TreatmentB", "Pr(>|z|)"]
  }
}

# Model 4: Andersen-Gill model

sampleAGA <- bootstrap_AG(
  AGdataA, nA, "A", i
)

sampleAGB <- bootstrap_AG(
  AGdataB, nB, "B", i
)

AG_boot_data <- bind_rows(sampleAGA, sampleAGB)

fit_AG <- tryCatch({
  
  coxph(
    Surv(start, stop, event) ~
      Treatment + strata(site) + cluster(studyid),
    data = AG_boot_data
  )
  
}, error = function(e) NULL)

if (!is.null(fit_AG)) {
  
  sm <- summary(fit_AG)$coefficients
  
  if ("TreatmentB" %in% rownames(sm)) {
    pvals_AG[i] <- sm["TreatmentB", "Pr(>|z|)"]
  }
}


fit_AG <- coxph(
  Surv(start, stop, event) ~ Treatment + strata(site) + cluster(studyid),
  data = cenic_AG
)

# HR, 95% CI, and p-value
summary(fit_AG)
exp(cbind(HR = coef(fit_AG), confint(fit_AG)))

# Model 5: Stratified Log-Rank

library(dplyr)

cenic_cox <- cenic_truncated %>%
  arrange(studyid, days_from_rand) %>%
  group_by(studyid) %>%
  summarise(
    Treatment = first(Treatment),
    site = first(site),
    
    # 7-day abstinence
    time_7day = ifelse(any(abst_7day == 1, na.rm = TRUE),
                       min(days_from_rand[abst_7day == 1], na.rm = TRUE),
                       max(days_from_rand)),
    event_7day = as.integer(any(abst_7day == 1, na.rm = TRUE)),
    
    # 3-day abstinence
    time_3day = ifelse(any(abst_3day == 1, na.rm = TRUE),
                       min(days_from_rand[abst_3day == 1], na.rm = TRUE),
                       max(days_from_rand)),
    event_3day = as.integer(any(abst_3day == 1, na.rm = TRUE)),
    
    # 7-day abstinence, excluding missing values
    time_7daynoNA = ifelse(any(abst_7daynoNA == 1, na.rm = TRUE),
                           min(days_from_rand[abst_7daynoNA == 1], na.rm = TRUE),
                           max(days_from_rand)),
    event_7daynoNA = as.integer(any(abst_7daynoNA == 1, na.rm = TRUE)),
    
    # 3-day abstinence, excluding missing values
    time_3daynoNA = ifelse(any(abst_3daynoNA == 1, na.rm = TRUE),
                           min(days_from_rand[abst_3daynoNA == 1], na.rm = TRUE),
                           max(days_from_rand)),
    event_3daynoNA = as.integer(any(abst_3daynoNA == 1, na.rm = TRUE)),
    
    .groups = "drop"
  )

logrank_results <- list(
  abst_7day = survdiff(Surv(time_7day, event_7day) ~ Treatment + strata(site), data = cenic_cox),
  abst_3day = survdiff(Surv(time_3day, event_3day) ~ Treatment + strata(site), data = cenic_cox),
  abst_7daynoNA = survdiff(Surv(time_7daynoNA, event_7daynoNA) ~ Treatment + strata(site), data = cenic_cox),
  abst_3daynoNA = survdiff(Surv(time_3daynoNA, event_3daynoNA) ~ Treatment + strata(site), data = cenic_cox)
)

do.call(rbind, lapply(logrank_results, function(x) {
  c(Statistic = x$chisq,
    P_value = pchisq(x$chisq, df = length(x$n) - 1, lower.tail = FALSE))
}))

#Stratified Cox proportional hazards
cox_results <- list(
  abst_7day = coxph(Surv(time_7day, event_7day) ~ Treatment + strata(site), data = cenic_cox),
  abst_3day = coxph(Surv(time_3day, event_3day) ~ Treatment + strata(site), data = cenic_cox),
  abst_7daynoNA = coxph(Surv(time_7daynoNA, event_7daynoNA) ~ Treatment + strata(site), data = cenic_cox),
  abst_3daynoNA = coxph(Surv(time_3daynoNA, event_3daynoNA) ~ Treatment + strata(site), data = cenic_cox)
)

# Extract hazard ratios, 95% CIs, and p-values
cox_summary <- do.call(rbind, lapply(names(cox_results), function(model) {
  fit <- cox_results[[model]]
  s <- summary(fit)
  
  data.frame(
    Model = model,
    HR = exp(coef(fit)),
    Lower_95_CI = exp(confint(fit)[1]),
    Upper_95_CI = exp(confint(fit)[2]),
    P_value = s$coefficients[1, "Pr(>|z|)"]
  )
}))

print(cox_summary)
