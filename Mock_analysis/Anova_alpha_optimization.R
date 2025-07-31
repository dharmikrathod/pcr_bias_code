library(tidyverse)
library(fido)
library(driver)
library(phyloseq)
library(Biostrings)
library(readxl)
library(ggstance)
library(gt)
library(mvrsquared)

pg <- "2020-03-05_PCRBias_qPCR_Picogreen.xlsx"
pg <- readxl::read_xlsx(pg)

dr <- "DNA dilutions.xlsx"
dr <- readxl::read_xlsx(dr)
dr[,2]/300*pg$PicogreenDNAConc

tmp <- "Mock Communities.xlsx"
mock_ratios <- readxl::read_xlsx(tmp, skip=1)
mock_ratios <- mock_ratios[-11,]
mock_ratios <- mock_ratios[,-1]
mock_ratios <- driver::miniclo_array(as.matrix(mock_ratios), parts=1)
colnames(mock_ratios) <- paste0("Mock_", 1:10)
rownames(mock_ratios) <- c("B.subtilis", "B.longum",
                           "C.hathewayi", "C.innocuum",
                           "C.aerofaciens", "E.faecalis",
                           "L.oris", "L.ruminis",
                           "R.intestinalis", "S.gallolyticus")
rm(tmp)

ps <- readRDS("phyloseq.rds")

# Filter low abundance samples
total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>1000, ps)
sum(sample_sums(ps))/total.reads

# read in seq clusters -- created by sequencing isolates individually
load("seq_clusters.RData")

seq_dictionary <- choices_seq
seq_dictionary[["S.aureus"]] <- NULL
for (n in names(seq_dictionary)) {
  names(seq_dictionary[[n]]) <- rep(n, length(seq_dictionary[[n]]))
}
seq_dict <- seq_dictionary[[1]]
for (n in 2:length(seq_dictionary)){
  seq_dict <- c(seq_dict, seq_dictionary[[n]])
}
rm(seq_dictionary)

# Start Matching
combined_seq <- c(seq_dict, refseq(ps))
seq_dist <-  pwalign::stringDist(combined_seq, method = "levenshtein")
n1 <- length(seq_dict)
n2 <- length(refseq(ps))
seq_dist <- as.matrix(seq_dist)
seq_dist <- seq_dist[1:n1, (n1+1):(n1+n2)]

seq_dist_min <- seq_dist %>% 
  array_branch(2) %>% 
  map(which.min)

seq_mapping <- data.frame(seq=names(seq_dist_min), 
                          name = unlist(map(seq_dist_min, names)), 
                          dist = unlist(seq_dist_min), 
                          tax_sum = taxa_sums(ps)[names(seq_dist_min)])


# Explore mapping
seq_mapping %>% 
  group_by(name) %>% 
  summarize(sum_count = sum(tax_sum), 
            mean_dist = mean(dist))
plot(density(seq_mapping$dist))


# Just double check results... 
hclust(stringDist(combined_seq, method = "levenshtein")) %>% plot()
# seq_17-seq_21 seem off in their own land while everything else falls nicely 
# within a given cluster
# Blasting seq_18 and seq_19 give some weird bacillus species that I have never
# heard of... either way they represent very few counts total I bet. 
ts <- taxa_sums(ps)
100-sum(ts[ts>500])/sum(ts)*100
# Yup only dropping 0.003% of counts
rm(ts)

# drop seq_variants with fewer than 500 counts based on histogram of taxa_counts
# these seem to be observed very infrequently -- also correspond to taxa that are not similar 
# to known isolates -- could be contamination
ps <- filter_taxa(ps, function(x) sum(x) > 500, TRUE)
seq_mapping <- filter(seq_mapping, seq %in% taxa_names(ps))


# Create data structures for modeling -------------------------------------

# Create new count table for modeling
Y <- as(t(otu_table(ps)),"matrix")
Y <- split(Y, seq_mapping$name) %>% 
  map(~matrix(.x, ncol=phyloseq::nsamples(ps)))


Y <- map(Y, colSums) %>% 
  do.call(rbind, .)
colnames(Y) <- sample_names(ps)

# Extract covariates from file names
sn <- colnames(Y)
d <- data.frame(sn = sn, 
                c1 = as.character(str_extract(sn, "^[:alpha:]*")), 
                c2 = as.numeric(str_extract(sn, "[:digit:]*(?=\\.)")), 
                c3 = as.numeric(str_extract(sn, "(?<=\\.)[:digit:]*$"))) %>% 
  mutate(sample_num = if_else(c1 == "cycle", "Calibration", paste0("Mock", c2)), 
         cycle_num = if_else(c1=="cycle", c2, 35), 
         machine = sample_data(ps)[sn,]$machine) %>% 
  mutate(machine = as.factor(machine))

# develop model on just calibration samples -------------------------------

# isolat
# Apply model to unknown mixes --------------------------------------------

X <- t(model.matrix(~sample_num + cycle_num + machine -1, data=d))
fit <- pibble(Y=Y, X=X, Gamma = 10*diag(nrow(X)))
fit <- to_proportions(fit)


newdata <- cbind(diag(11), diag(11))
newdata <- rbind(newdata, rep(c(0,35), each=11), matrix(0, 3, 22))
sn <- str_sub(rownames(X)[1:11], start=11)
rownames(newdata) <- rownames(X)
newdata_meta <- data.frame(sample = rep(sn, times=2),
                           cycle = rep(c(0, 35), each=11))
dimnames(LambdaX)[[1]] <-  names_categories(fit)
cycle0<-newdata[,1:11]
cycle35<-newdata[,12:22]


mockLambdaX_0<- predict(fit, cycle0, pars="LambdaX")
mockLambdaX_35<-predict(fit, cycle35, pars = "LambdaX")

dimnames(mockLambdaX_0)[[1]] <- names_categories(fit)
dimnames(mockLambdaX_35)[[1]] <- names_categories(fit)
dimnames(mockLambdaX_0) <- list(taxon = 1:dim(mockLambdaX_0)[1],
                                sample = 1:dim(mockLambdaX_0)[2],
                                iteration = 1:dim(mockLambdaX_0)[3])

df_long <- as.data.frame.table(mockLambdaX_0, responseName = "proportion") %>%
  mutate(
    sample = as.integer(sample),
    iteration = as.integer(iteration),
    taxon = as.integer(taxon)
  )

# Now compute Shannon index per sample × iteration
LambdaX_shannon_0 <- df_long %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cycle_name = 0) %>%
  arrange(sample, iteration)

dimnames(mockLambdaX_35) <- list(taxon = 1:dim(mockLambdaX_35)[1],
                                sample = 1:dim(mockLambdaX_35)[2],
                                iteration = 1:dim(mockLambdaX_35)[3])

df_long1 <- as.data.frame.table(mockLambdaX_35, responseName = "proportion") %>%
  mutate(
    sample = as.integer(sample),
    iteration = as.integer(iteration),
    taxon = as.integer(taxon)
  )

# Now compute Shannon index per sample × iteration
LambdaX_shannon_35 <- df_long1 %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cycle_name = 0) %>%
  arrange(sample, iteration)

mean_shannon_0 <- df_long %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  summarize(
    mean_shannon = mean(shannon, na.rm = TRUE),
    .groups = "drop"
  )
mean_shannon_35 <- df_long1 %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  summarize(
    mean_shannon = mean(shannon, na.rm = TRUE),
    .groups = "drop"
  )



library(GA)  # Genetic algorithm package
library(vegan)

mean_matrix <- apply(LambdaX, c(1, 2), mean)

mean_shannon_0$cycle <- 0
mean_shannon_35$cycle <- 35

mean_shannon_all <- bind_rows(mean_shannon_0, mean_shannon_35)
# Suppose you have:
# X0: 11x10 matrix of posterior means at cycle 0
# X35: 11x10 matrix of posterior means at cycle 35
library(GA)

# Create lookup tables for cycle 0 and cycle 35
shannon_by_sample <- mean_shannon_all %>%
  pivot_wider(names_from = cycle, values_from = mean_shannon, names_prefix = "cycle_")

# Ensure sorted by sample number
shannon_by_sample <- shannon_by_sample %>% arrange(sample)

# Objective function
objective_shannon <- function(grouping) {
  if (length(unique(grouping)) != 2) return(-Inf)
  
  group <- as.factor(grouping)
  
  r2_cycle <- function(values) {
    summary(lm(values ~ group))$r.squared
  }
  
  r2_0 <- r2_cycle(shannon_by_sample$cycle_0)
  r2_35 <- r2_cycle(shannon_by_sample$cycle_35)
  
  abs(r2_35 - r2_0)
}
ga_result <- ga(
  type = "binary",
  fitness = objective_shannon,
  nBits = 11,
  maxiter = 1000,
  run = 100,
  seed = 123,
  monitor = TRUE
)

# Optimal grouping
best_grouping <- ga_result@fitnessValue



library(dplyr)
library(tidyr)
library(GA)

# === STEP 1: Set dimnames for 3D arrays ===


# === STEP 2: Convert to long format ===
df_long_0 <- as.data.frame.table(mockLambdaX_0, responseName = "proportion") %>%
  mutate(across(c(taxon, sample, iteration), as.integer))

df_long_35 <- as.data.frame.table(mockLambdaX_35, responseName = "proportion") %>%
  mutate(across(c(taxon, sample, iteration), as.integer))

# === STEP 3: Compute Simpson index per sample per iteration ===
LambdaX_simpson_0 <- df_long_0 %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson = 1 - sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  summarize(mean_simpson = mean(simpson, na.rm = TRUE), .groups = "drop")

LambdaX_simpson_35 <- df_long_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson = 1 - sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  summarize(mean_simpson = mean(simpson, na.rm = TRUE), .groups = "drop")

# === STEP 4: Combine and reshape for optimization ===
LambdaX_simpson_0$cycle <- 0
LambdaX_simpson_35$cycle <- 35

mean_simpson_all <- bind_rows(LambdaX_simpson_0, LambdaX_simpson_35)

simpson_by_sample <- mean_simpson_all %>%
  pivot_wider(names_from = cycle, values_from = mean_simpson, names_prefix = "cycle_") %>%
  arrange(sample)

# === STEP 5: Define the optimization objective ===
objective_simpson <- function(grouping) {
  if (length(unique(grouping)) != 2) return(-Inf)
  
  group <- as.factor(grouping)
  
  r2_cycle <- function(values) {
    summary(lm(values ~ group))$r.squared
  }
  
  r2_0 <- r2_cycle(simpson_by_sample$cycle_0)
  r2_35 <- r2_cycle(simpson_by_sample$cycle_35)
  
  abs(r2_35 - r2_0)
}

# === STEP 6: Run the genetic algorithm ===
ga_result_simpson <- ga(
  type = "binary",
  fitness = objective_simpson,
  nBits = 11,  # assuming 11 samples
  maxiter = 1000,
  run = 100,
  seed = 123,
  monitor = TRUE
)

# === STEP 7: Extract best grouping ===
best_grouping_simpson <- ga_result_simpson@solution


library(dplyr)
library(tidyr)
library(GA)
library(ineq)  
install.packages("ineq")# for Gini coefficient


# === STEP 2: Convert to long format ===
df_long_0 <- as.data.frame.table(mockLambdaX_0, responseName = "proportion") %>%
  mutate(across(c(taxon, sample, iteration), as.integer))

df_long_35 <- as.data.frame.table(mockLambdaX_35, responseName = "proportion") %>%
  mutate(across(c(taxon, sample, iteration), as.integer))

# === STEP 3: Compute Gini coefficient per sample × iteration ===
LambdaX_gini_0 <- df_long_0 %>%
  group_by(sample, iteration) %>%
  summarize(
    gini = ineq::Gini(proportion),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  summarize(mean_gini = mean(gini, na.rm = TRUE), .groups = "drop")

LambdaX_gini_35 <- df_long_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    gini = ineq::Gini(proportion),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  summarize(mean_gini = mean(gini, na.rm = TRUE), .groups = "drop")

# === STEP 4: Combine for optimization ===
LambdaX_gini_0$cycle <- 0
LambdaX_gini_35$cycle <- 35

mean_gini_all <- bind_rows(LambdaX_gini_0, LambdaX_gini_35)

gini_by_sample <- mean_gini_all %>%
  pivot_wider(names_from = cycle, values_from = mean_gini, names_prefix = "cycle_") %>%
  arrange(sample)

# === STEP 5: Objective function ===
objective_gini <- function(grouping) {
  if (length(unique(grouping)) != 2) return(-Inf)
  
  group <- as.factor(grouping)
  
  r2_cycle <- function(values) {
    summary(lm(values ~ group))$r.squared
  }
  
  r2_0 <- r2_cycle(gini_by_sample$cycle_0)
  r2_35 <- r2_cycle(gini_by_sample$cycle_35)
  
  abs(r2_35 - r2_0)
}

# === STEP 6: Run GA optimization ===
ga_result_gini <- ga(
  type = "binary",
  fitness = objective_gini,
  nBits = 11,  # adjust if different number of samples
  maxiter = 1000,
  run = 100,
  seed = 123,
  monitor = TRUE
)

# === STEP 7: Best grouping ===
best_grouping_gini <- ga_result_gini@solution



aitchison_norm <- function(p) {
  p <- p + 1e-10  # small pseudocount to avoid log(0)
  clr_p <- compositions::clr(p)
  sqrt(sum(clr_p^2))
}

# Compute norms at cycle 0
norms_0 <- apply(mockLambdaX_0, c(2, 3), function(p) aitchison_norm(p))
mean_norm_0 <- rowMeans(norms_0, na.rm = TRUE)
LambdaX_aitchison_0 <- data.frame(
  sample = 1:length(mean_norm_0),
  mean_aitchison = mean_norm_0,
  cycle = 0
)

# Compute norms at cycle 35
norms_35 <- apply(mockLambdaX_35, c(2, 3), function(p) aitchison_norm(p))
mean_norm_35 <- rowMeans(norms_35, na.rm = TRUE)
LambdaX_aitchison_35 <- data.frame(
  sample = 1:length(mean_norm_35),
  mean_aitchison = mean_norm_35,
  cycle = 35
)

# === STEP 3: Combine for optimization ===
mean_aitchison_all <- bind_rows(LambdaX_aitchison_0, LambdaX_aitchison_35)

aitchison_by_sample <- mean_aitchison_all %>%
  pivot_wider(names_from = cycle, values_from = mean_aitchison, names_prefix = "cycle_") %>%
  arrange(sample)

# === STEP 4: Define objective function ===
objective_aitchison <- function(grouping) {
  if (length(unique(grouping)) != 2) return(-Inf)
  
  group <- as.factor(grouping)
  
  r2_cycle <- function(values) {
    summary(lm(values ~ group))$r.squared
  }
  
  r2_0 <- r2_cycle(aitchison_by_sample$cycle_0)
  r2_35 <- r2_cycle(aitchison_by_sample$cycle_35)
  
  abs(r2_35 - r2_0)
}

# === STEP 5: Run GA ===
ga_result_aitchison <- ga(
  type = "binary",
  fitness = objective_aitchison,
  nBits = 11,  # adjust if sample size changes
  maxiter = 1000,
  run = 100,
  seed = 123,
  monitor = TRUE
)

# === STEP 6: Output grouping ===
best_grouping_aitchison <- ga_result_aitchison@solution
shannon<- round(ga_result@fitnessValue,2)
simpson<- round(ga_result_simpson@fitnessValue,2)
gini<- round(ga_result_gini@fitnessValue,2)
atichison<- round(ga_result_aitchison@fitnessValue,2)
data.frame("Max_R2_difference" = rbind(shannon, simpson, gini, atichison))
