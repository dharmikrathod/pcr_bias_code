
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


# Apply model to unknown mixes --------------------------------------------

X <- t(model.matrix(~sample_num + cycle_num + machine -1, data=d))
fit_mock <- pibble(Y=Y, X=X, Gamma = 10*diag(nrow(X)))
fit_mock <- to_proportions(fit_mock)


newdata <- cbind(diag(11), diag(11))
newdata <- rbind(newdata, rep(c(0,35), each=11), matrix(0, 3, 22))
sn <- str_sub(rownames(X)[1:11], start=11)
rownames(newdata) <- rownames(X)
newdata_meta <- data.frame(sample = rep(sn, times=2),
                           cycle = rep(c(0, 35), each=11))


dimnames(LambdaX)[[1]] <-  names_categories(fit)
cycle0<-newdata[,1:11]
cycle35<-newdata[,12:22]


mockLambdaX_0<- predict(fit_mock, cycle0, pars="LambdaX")
mockLambdaX_35<-predict(fit_mock, cycle35, pars = "LambdaX")

dimnames(mockLambdaX_0)[[1]] <- names_categories(fit)
dimnames(mockLambdaX_35)[[1]] <- names_categories(fit)


library(phyloseq)
library(ape)

# Load the tree
tree <- read.tree("mock_pruned.asv.nwk")

# Confirm the number of taxa
length(tree$tip.label)  # should be 69
# LambdaX_zero: (samples x taxa x draws)
# Transpose to (taxa x samples x draws) if needed
LambdaX_zero <- aperm(mockLambdaX_0, c(1,2, 3))
LambdaX_35 <- aperm(mockLambdaX_35, c(1,2, 3))

# Add taxon names from tree
dimnames(mockLambdaX_0)[[1]] <- tree$tip.label
dimnames(mockLambdaX_35)[[1]] <- tree$tip.label


n_draws <- dim(mockLambdaX_0)[3]
unifrac_zero_list <- vector("list", n_draws)
unifrac_35_list <- vector("list", n_draws)

for (i in 1:n_draws) {
  # Get abundance matrices (taxa x samples)
  abund_zero <- mockLambdaX_0[,,i] + 1e-8  # pseudo-count for numerical safety
  abund_35 <- mockLambdaX_35[,,i] + 1e-8
  
  # Construct phyloseq objects
  ps_zero <- phyloseq(otu_table(abund_zero, taxa_are_rows = TRUE), phy_tree(tree))
  ps_35 <- phyloseq(otu_table(abund_35, taxa_are_rows = TRUE), phy_tree(tree))
  
  # Compute UniFrac distances
  unifrac_zero_list[[i]] <- UniFrac(ps_zero, weighted = TRUE, normalized = TRUE)
  unifrac_35_list[[i]] <- UniFrac(ps_35, weighted = TRUE, normalized = TRUE)
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}

# Sample pairs to extract (as row/column indices)
pair_indices <- combn(1:11, 2, simplify = FALSE)


# Helper function to extract distances from list of matrices
extract_unifrac <- function(unifrac_list, cycle_label) {
  do.call(rbind, lapply(seq_along(unifrac_list), function(i) {
    m <- as.matrix(unifrac_list[[i]])  # convert dist to matrix
    do.call(rbind, lapply(pair_indices, function(pair) {
      data.frame(
        Draw = i,
        Pair = paste0(pair[1], "v", pair[2]),
        Distance = m[pair[1], pair[2]],
        Cycle = cycle_label
      )
    }))
  }))
}

# Combine all into a single data frame
df_unifrac <- rbind(
  extract_unifrac(unifrac_zero_list, "Cycle 0"),
  extract_unifrac(unifrac_35_list, "Cycle 35")
)
bias_list <- lapply(pair_indices, function(pair) {
  pair_label <- paste0(pair[1], "v", pair[2])
  bias_vals <- sapply(1:length(unifrac_35_list), function(i) {
    m_35 <- as.matrix(unifrac_35_list[[i]])
    m_0 <- as.matrix(unifrac_zero_list[[i]])
    m_0[pair[1], pair[2]] - m_35[pair[1], pair[2]]
  })
  data.frame(
    Pair = pair_label,
    Bias = bias_vals
  )
})

# Combine into a long data frame
df_bias <- do.call(rbind, bias_list)


library(dplyr)
library(ggplot2)

# Summarize posterior bias distribution per pair
summary_bias <- df_bias %>%
  group_by(Pair) %>%
  summarize(
    p2.5 = quantile(Bias, 0.025),
    p50 = quantile(Bias, 0.5),
    p97.5 = quantile(Bias, 0.975),
    .groups = "drop"
  )

# Extract Sample1 and Sample2 from Pair name like "1v2"
summary_bias <- summary_bias %>%
  mutate(
    Sample1 = sub("v.*", "", Pair),
    Sample2 = sub(".*v", "", Pair)
  )

# Create dynamic factor levels for all 11 samples
sample_levels <- as.character(1:11)

# Prepare for heatmap (keep only lower triangle: Sample1 < Sample2)
summary_bias_lower <- summary_bias %>%
  mutate(
    Sample1 = factor(Sample1, levels = sample_levels),
    Sample2 = factor(Sample2, levels = sample_levels)
  ) %>%
  filter(as.numeric(as.character(Sample1)) < as.numeric(as.character(Sample2)))

# Plot heatmap with median + credible intervals
ggplot(summary_bias_lower, aes(x = Sample2, y = Sample1, fill = p50)) +
  geom_tile() +
  geom_text(aes(label = signif(p50, 2)),
            color = "white", size = 4, nudge_y = 0.15) +
  geom_text(aes(label = paste0("(", round(signif(p2.5, 2),2), "-", round(signif(p97.5, 2),2), ")")),
            color = "white", size = 3.5, nudge_y = -0.15) +
  scale_fill_gradient(low = "black", high = "#56B4E9", name = "UniFrac bias") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = 'none'
  ) 

ibrary(vegan)

# Use existing: LambdaX_zero and LambdaX_35 (taxa x samples x draws)
n_samples <- dim(LambdaX_zero)[2]
n_draws <- dim(LambdaX_zero)[3]
pair_indices <- combn(1:n_samples, 2, simplify = FALSE)

# Bray-Curtis bias per pair per draw
bias_list <- lapply(pair_indices, function(pair) {
  pair_label <- paste0(pair[1], "v", pair[2])
  
  bias_vals <- sapply(1:n_draws, function(i) {
    samp_0 <- LambdaX_zero[, c(pair[1], pair[2]), i]
    samp_35 <- LambdaX_35[, c(pair[1], pair[2]), i]
    
    # Normalize proportions
    samp_0 <- samp_0 / colSums(samp_0)
    samp_35 <- samp_35 / colSums(samp_35)
    
    # Bray-Curtis
    dist_0 <- vegdist(t(samp_0), method = "bray")[1]
    dist_35 <- vegdist(t(samp_35), method = "bray")[1]
    
    return(dist_35 - dist_0)
  })
  
  data.frame(Pair = pair_label, Bias = bias_vals)
})

# Combine all into a single data frame
df_bias <- do.call(rbind, bias_list)
summary_bias <- df_bias %>%
  group_by(Pair) %>%
  summarize(
    p2.5 = quantile(Bias, 0.025),
    p50  = quantile(Bias, 0.5),
    p97.5 = quantile(Bias, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    Sample1 = sub("v.*", "", Pair),
    Sample2 = sub(".*v", "", Pair)
  )

sample_levels <- as.character(1:n_samples)

summary_bias_lower <- summary_bias %>%
  mutate(
    Sample1 = factor(Sample1, levels = sample_levels),
    Sample2 = factor(Sample2, levels = sample_levels)
  ) %>%
  filter(as.numeric(as.character(Sample1)) < as.numeric(as.character(Sample2)))
library(ggplot2)

ggplot(summary_bias_lower, aes(x = Sample2, y = Sample1, fill = p50)) +
  geom_tile() +
  geom_text(aes(label = signif(p50, 2)),
            color = "white", size = 4, nudge_y = 0.15) +
  geom_text(aes(label = paste0("(", round(signif(p2.5, 2),2), "-", round(signif(p97.5, 2),2), ")")),
            color = "white", size = 3.5, nudge_y = -0.15) +
  scale_fill_gradient(low = "black", high = "#56B4E9", name = "Bray-Curtis bias") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = 'none'
  ) 

