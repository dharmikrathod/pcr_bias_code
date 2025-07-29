library(tidyverse)
library(driver) # devtools::install_github("jsilve24/driver")
library(fido) # devtools::install_github("jsilve24/stray")
library(phyloseq)
library(NVC) # devtools::install_github("jsilve24/NVC")
library(DECIPHER) # for computing sequence similarity
library(tidybayes)
library(ggforce)
library(ggrepel)
library(Biostrings)



set.seed(18845)
ps <- readRDS("phyloseq.rds")
path.addedmeta <- "2017.08.14MappingFile_addedmetadata.txt"
sample_data(ps) <- phyloseq::import_qiime_sample_data(path.addedmeta)

total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>5000, ps)
sum(sample_sums(ps))/total.reads

ps <- tax_glom(ps, "Genus")
tmp <- taxa_names(ps)
d <- as(otu_table(ps), "matrix")
d <- apply(d, 2, function(x) sum(x > 3) > 0.30*length(x))
tmp <- tmp[!d]
ps <- merge_taxa(ps, tmp, 1)
other <- taxa_names(ps)[taxa_names(ps) %in% tmp]
tmp <- taxa_names(ps)
tmp[tmp==other] <- "other"
taxa_names(ps) <- tmp
rm(tmp, d, other)

otu_table(ps)

taxa_names(ps)
# Prune samples where BiasPCRMachine is not equal to 3
ps <- prune_samples(sample_data(ps)$BiasPCRMachine != 3, ps)

# Convert BiasDonorName to a factor with specified levels
sample_data(ps)$BiasDonorName <- factor(sample_data(ps)$BiasDonorName,
                                        levels = c("16S8921", "908", "Ai96", "J1122526T"))

# fit stray -------------------------------------------------------------

# Prepare Y and X matrices
Y <- as(t(otu_table(ps)), "matrix")
X <- t(as.matrix(model.matrix(~BiasPCRCycleNum + factor(BiasPCRMachine),
                              as(sample_data(ps), "data.frame"))))

# One-hot encode donor names based on the factor levels
X.donor <- onehot(sample_data(ps)$BiasDonorName)
colnames(X.donor) <- c("16S8921", "908", "Ai96", "J1122526T")
X <- rbind(t(X.donor), X[2:nrow(X),])
rm(X.donor)

# Clean Names for PCR Machines
rownames(X)[6:8] <- c("I(Machine2)", "I(Machine4)", "I(Machine5)")

# form priors
A <- create_alr_base(ntaxa(ps), ntaxa(ps))
Xi <- (t(A) %*% diag(ntaxa(ps)) %*% A)/2 # divide by 2 to scale to unit diagonal
Gamma <- diag(c(15, 15, 15, 15, 1, 1, 1, 1)) 
upsilon <- ntaxa(ps)+3
Theta <- matrix(0, ntaxa(ps)-1,ncol(Gamma))





# fit model
fit <- pibble(Y, X, upsilon, Theta, Gamma, Xi)




fit <- to_clr(fit)
fit <- fido:::name.pibblefit(fit)
coord.names <- dimnames(fit$Eta)[[1]]
dimnames(fit$Lambda)[2]
d <- (Y+0.65) %>%
  clr_array(1) %>%
  fido::gather_array(val, coord, sample) %>% # Specify fido::gather_array
  mutate(pcr_cycle = fit$X[5,][sample],
         sample_name = sample_names(ps)[sample],
         PCRMachine = factor(sample_data(ps)$BiasPCRMachine[sample]),
         Donor = sample_data(ps)$BiasDonorName[sample],
         mean=val) %>%
  mutate(coord = coord.names[coord])



fit<- to_proportions(fit)




# Plot Predicted Difference -----------------------------------------------

focus.names <- c("clr_seq_9" = "Parasutterella", 
                 "clr_seq_48" = "Bifidobacterium", 
                 "clr_seq_1" = "Bacteroides", 
                 "clr_seq_192" = "Ruminococcus", 
                 "clr_seq_204" = "Coprococcus_1", 
                 "clr_seq_193" = "Holdemania")


#predicting proportions for the first sample (0th cycle)
LambdaX_zero <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0, 0, 0, 0),
                                           c(0, 1, 0, 0, 0, 0, 0, 0),
                                           c(0, 0, 1, 0, 0, 0, 0, 0),
                                           c(0, 0, 0, 1, 0, 0, 0, 0)
),
pars="LambdaX", use_names=FALSE) 
#predicting proportions for the first sample (35th cycle)
LambdaX_35 <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 35, 0, 0, 0),
                                         c(0, 1, 0, 0, 35, 0, 0, 0),
                                         c(0, 0, 1, 0, 35, 0, 0, 0),
                                         c(0, 0, 0, 1, 35, 0, 0, 0)
),
pars="LambdaX", use_names=FALSE) 



LambdaX_zero <- LambdaX_zero %>% 
  fido::gather_array(val, coord, sample, iter) %>% 
  group_by(coord, sample) %>%
  mutate(
    pcr_cycle = c(0, 0, 0, 0, fit$X[5,])[sample],
    Donor = c("16S8921", "908", "Ai96", "J1122526T", 
              as.character(sample_data(ps)$BiasDonorName))[sample],
    coord = coord.names[coord]
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = iter, values_from = val)

#Shannons Index
LambdaX_normalized <- LambdaX_zero %>%
  pivot_longer(cols = `1`:`2000`, names_to = "iteration", values_to = "proportion") %>%
  group_by(sample, iteration) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  ungroup()

LambdaX_shannon_wide_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = shannon_index
  ) %>%
  arrange(iteration)

# Compute Shannon's index for each sample and posterior iteration

LambdaX_35 <- LambdaX_35 %>% 
  fido::gather_array(val, coord, sample, iter) %>% 
  group_by(coord, sample) %>%
  mutate(
    pcr_cycle = c(0, 0, 0, 0, fit$X[5,])[sample],
    Donor = c("16S8921", "908", "Ai96", "J1122526T", 
              as.character(sample_data(ps)$BiasDonorName))[sample],
    coord = coord.names[coord]
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = iter, values_from = val)

# Compute Shannon's index for each taxa
LambdaX_normalized_35 <- LambdaX_35 %>%
  pivot_longer(cols = `1`:`2000`, names_to = "iteration", values_to = "proportion") %>%
  group_by(sample, iteration) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  
  ungroup()



# Compute Shannon's index for each sample and posterior iteration


LambdaX_shannon_wide_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = shannon_index
  ) %>%
  arrange(iteration)

sample_columns <- intersect(colnames(LambdaX_shannon_wide_35), colnames(LambdaX_shannon_wide_0))
sample_columns <- sample_columns[sample_columns != "iteration"]



# Function to compute ANOVA and R²
compute_anova_r2 <- function(df, group_label) {
  # Convert to long format
  df_long <- df %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "value")
  
  # Fit ANOVA model
  aov_model <- aov(value ~ sample, data = df_long)
  
  # Extract sum of squares
  ss_total   <- sum(summary(aov_model)[[1]][["Sum Sq"]])
  ss_between <- summary(aov_model)[[1]][["Sum Sq"]][1]
  
  # Compute R²
  r_squared <- ss_between / ss_total
  
  # Return summary and R²
  list(
    label = group_label,
    summary = summary(aov_model),
    R2 = r_squared
  )
}

# Apply to both dataframes
result_0  <- compute_anova_r2(LambdaX_shannon_wide_0, "0th cycle")
result_35 <- compute_anova_r2(LambdaX_shannon_wide_35, "35th cycle")

# View results
result_0$R2
result_35$R2


# R² function for ANOVA model
get_r2 <- function(aov_result) {
  ss <- summary(aov_result)[[1]][["Sum Sq"]]
  r_squared <- ss[1] / sum(ss)
  return(r_squared)
}

pairwise_anova_r2 <- function(df, label = "") {
  samples <- setdiff(names(df), "iteration")
  results <- data.frame()
  
  for (i in 1:(length(samples) - 1)) {
    for (j in (i + 1):length(samples)) {
      s1 <- samples[i]
      s2 <- samples[j]
      
      # Reshape to long format for ANOVA
      long_df <- df %>%
        select(iteration, all_of(c(s1, s2))) %>%
        pivot_longer(cols = -iteration, names_to = "sample", values_to = "value") %>%
        mutate(sample = factor(sample))
      
      # Fit ANOVA model
      model <- tryCatch(aov(value ~ sample, data = long_df), error = function(e) NULL)
      
      if (!is.null(model)) {
        anova_table <- summary(model)[[1]]
        ss_sample <- anova_table["sample", "Sum Sq"]
        ss_total <- sum(anova_table[["Sum Sq"]])
        r2 <- ss_sample / ss_total
        
        results <- rbind(results, data.frame(
          Sample1 = s1,
          Sample2 = s2,
          Cycle = label,
          R2 = r2
        ))
      }
    }
  }
  
  return(results)
}

# Function to run pairwise ANOVA and extract R²
pairwise_anova_r2 <- function(df, group_label) {
  sample_cols <- colnames(df)[-1]  # Assumes first column is 'iteration' or grouping
  
  # All pairwise combinations
  results <- combn(sample_cols, 2, function(pair) {
    sample1 <- df[[pair[1]]]
    sample2 <- df[[pair[2]]]
    
    combined_df <- data.frame(
      value = c(sample1, sample2),
      sample = factor(rep(pair, each = length(sample1)))
    )
    
    # Perform ANOVA
    aov_result <- aov(value ~ sample, data = combined_df)
    
    # Extract R²
    r2 <- get_r2(aov_result)
    
    # Return a named data frame
    data.frame(
      pair = paste(pair, collapse = " vs "),
      group = group_label,
      r_squared = r2,
      p_value = summary(aov_result)[[1]][["Pr(>F)"]][1]
    )
  }, simplify = FALSE)
  
  # Combine results into a single data frame
  do.call(rbind, results)
}

str(LambdaX_shannon_wide_0)

# Run for both dataframes
r2_results_0  <- pairwise_anova_r2(LambdaX_shannon_wide_0, "0th cycle")
r2_results_35 <- pairwise_anova_r2(LambdaX_shannon_wide_35, "35th cycle")

# Combine all results
test_df <- LambdaX_shannon_wide_0 %>%
  select(iteration, `1`, `2`)

pairwise_anova_r2(test_df, "0th cycle")

pairwise_anova_results <- rbind(r2_results_0, r2_results_35)

# View
pairwise_anova_results$diff<-r2_results_0$r_squared - r2_results_35$r_squared
pairwise_anova_results




# Simpson Index: 1 - sum(p^2)
LambdaX_simpson_wide_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson_index = sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(names_from = sample, values_from = simpson_index) %>%
  arrange(iteration)

LambdaX_simpson_wide_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson_index = 1 - sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(names_from = sample, values_from = simpson_index) %>%
  arrange(iteration)

# Gini Index: 1 - sum(2 * i * p_i), approx via sorted proportions
LambdaX_gini_index <- function(p) {
  p <- sort(p[!is.na(p)])
  n <- length(p)
  if (n == 0) return(NA)
  gini <- 1 - sum((2 * (1:n) - n - 1) * p) / (n - 1)
  return(gini)
}

LambdaX_gini_wide_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    gini_index = LambdaX_gini_index(proportion),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(names_from = sample, values_from = gini_index) %>%
  arrange(iteration)

LambdaX_gini_wide_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    gini_index = LambdaX_gini_index(proportion),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(names_from = sample, values_from = gini_index) %>%
  arrange(iteration)

# For cycle 0
LambdaX_aitchison_wide_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    aitchison_norm = sqrt(sum((log(proportion) - mean(log(proportion)))^2)),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(names_from = sample, values_from = aitchison_norm) %>%
  arrange(iteration)

# For cycle 35
LambdaX_aitchison_wide_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    aitchison_norm = sqrt(sum((log(proportion) - mean(log(proportion)))^2)),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(names_from = sample, values_from = aitchison_norm) %>%
  arrange(iteration)



# Re-use the same function:

r2_simpson_0  <- pairwise_anova_r2(LambdaX_simpson_wide_0, "0th cycle")
r2_simpson_35 <- pairwise_anova_r2(LambdaX_simpson_wide_35, "35th cycle")

r2_gini_0  <- pairwise_anova_r2(LambdaX_gini_wide_0, "0th cycle")
r2_gini_35 <- pairwise_anova_r2(LambdaX_gini_wide_35, "35th cycle")

r2_aitchison_0  <- pairwise_anova_r2(LambdaX_aitchison_wide_0, "0th cycle")
r2_aitchison_35 <- pairwise_anova_r2(LambdaX_aitchison_wide_35, "35th cycle")


pairwise_aitchison_results <- rbind(r2_aitchison_0, r2_aitchison_35)
pairwise_aitchison_results$diff <- r2_aitchison_0$r_squared - r2_aitchison_35$r_squared

# View






simpson_results <- combine_r2_tables(r2_simpson_0, r2_simpson_35)
gini_results    <- combine_r2_tables(r2_gini_0, r2_gini_35)


pairwise_simpson_results <- rbind(r2_simpson_0, r2_simpson_35)
pairwise_gini_results <- rbind(r2_gini_0, r2_gini_35)

# View
pairwise_gini_results$diff<-r2_gini_0$r_squared - r2_gini_35$r_squared
pairwise_simpson_results$diff<-r2_simpson_0$r_squared - r2_simpson_35$r_squared
pairwise_gini_results
pairwise_simpson_results

