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




print(dim(Y))
print(head(Y[, as.numeric(sample_data(ps)$BiasDonorName) == i]))

# fit model
fit <- pibble(Y, X, upsilon, Theta, Gamma, Xi)


p <- ppc(fit)


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

library(phyloseq)
library(ape)

# Load the tree
tree <- read.tree("pruned_tree.asv.nwk")

# Confirm the number of taxa
length(tree$tip.label)  # should be 69
# LambdaX_zero: (samples x taxa x draws)
# Transpose to (taxa x samples x draws) if needed
LambdaX_zero <- aperm(LambdaX_zero, c(2, 1, 3))
LambdaX_35 <- aperm(LambdaX_35, c(2, 1, 3))

# Add taxon names from tree
dimnames(LambdaX_zero)[[1]] <- tree$tip.label
dimnames(LambdaX_35)[[1]] <- tree$tip.label


n_draws <- dim(LambdaX_zero)[3]
unifrac_zero_list <- vector("list", n_draws)
unifrac_35_list <- vector("list", n_draws)

for (i in 1:n_draws) {
  # Get abundance matrices (taxa x samples)
  abund_zero <- LambdaX_zero[,,i] + 1e-8  # pseudo-count for numerical safety
  abund_35 <- LambdaX_35[,,i] + 1e-8
  
  # Construct phyloseq objects
  ps_zero <- phyloseq(otu_table(abund_zero, taxa_are_rows = TRUE), phy_tree(tree))
  ps_35 <- phyloseq(otu_table(abund_35, taxa_are_rows = TRUE), phy_tree(tree))
  
  # Compute UniFrac distances
  unifrac_zero_list[[i]] <- UniFrac(ps_zero, weighted = TRUE, normalized = TRUE)
  unifrac_35_list[[i]] <- UniFrac(ps_35, weighted = TRUE, normalized = TRUE)
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}

# Sample pairs to extract (as row/column indices)
pair_indices <- list(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(2, 3),
  c(2, 4),
  c(3, 4)
)

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


library(ggplot2)

ggplot(df_unifrac, aes(x = Pair, y = Distance, fill = Cycle)) +
  geom_violin(alpha = 0.6, position = position_dodge(0.8)) +
  labs(
    title = "Posterior Weighted UniFrac Distance by Sample Pair",
    x = "Sample Pair",
    y = "Weighted UniFrac Distance"
  ) +
  scale_fill_manual(values = c("Cycle 0" = "#1b9e77", "Cycle 35" = "#d95f02")) +
  theme_minimal(base_size = 14)


# Define the sample pairs
pair_indices <- list(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(2, 3),
  c(2, 4),
  c(3, 4)
)

# Extract bias values for each pair over all posterior draws
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

# Okabe-Ito colorblind-friendly palette (up to 8 distinct colors)
pair_labels <- unique(df_bias$Pair)

# Choose 6 well-separated, colorblind-friendly colors manually
=
library(ggplot2)

ggplot(df_bias, aes(x = Bias, fill = Pair, color = Pair)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = pair_colors) +
  scale_color_manual(values = rep('black',6)) +
  labs(
    title = "Posterior Density of UniFrac Bias (Cycle 35 − Cycle 0)",
    x = "Change in Weighted UniFrac",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())


library(dplyr)

summary_bias <- df_bias %>%
  group_by(Pair) %>%
  summarize(
    p2.5 = quantile(Bias, 0.025),
    p50 = quantile(Bias, 0.5),
    p97.5 = quantile(Bias, 0.975),
    .groups = "drop"
  )
summary_bias <- summary_bias %>%
  mutate(
    Sample1 = sub("v.*", "", Pair),
    Sample2 = sub(".*v", "", Pair)
  )
# Keep original sample order (not reversed)
sample_levels <- c("1", "2", "3", "4")

summary_bias_lower <- summary_bias %>%
  mutate(
    Sample1 = factor(Sample1, levels = sample_levels),  # x-axis
    Sample2 = factor(Sample2, levels = sample_levels)   # y-axis
  ) %>%
  filter(as.numeric(as.character(Sample1)) < as.numeric(as.character(Sample2)))  # Lower triangle
summary_bias_lower <- summary_bias_lower %>%
  mutate(
    Sample1 = paste("Sample", Sample1),
    Sample2 = paste("Sample", Sample2)
  )


ggplot(summary_bias_lower, aes(x = Sample2, y = Sample1, fill = p50)) +
  geom_tile() +
  geom_text(aes(label = signif(p50, 2)),
            color = "white", size = 5, nudge_y = 0.12) +
  geom_text(aes(label = paste0("(", signif(p2.5, 2), "-", signif(p97.5, 2), ")")),
            color = "white", size = 4, nudge_y = -0.12) +
  scale_fill_gradient(low = "black", high = "#E69F00", name = "UniFrac bias") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    #axis.title = element_blank()
  ) +
  labs(title = "Posterior Change in Weighted UniFrac (Cycle 35 − Cycle 0)")


ggplot(summary_bias_lower, aes(x = Sample2, y = Sample1, fill = p50)) +
  geom_tile() +
  geom_text(aes(label = signif(p50, 2)),
            color = "white", size = 7.5, nudge_y = 0.15) +
  geom_text(aes(label = paste0("(", signif(p2.5, 2), "-", signif(p97.5, 2), ")")),
            color = "white", size = 6.5, nudge_y = -0.15) +
  scale_fill_gradient(low = "black", high = "#56B4E9", name = "UniFrac bias") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Posterior Change in Weighted UniFrac (Cycle 35 − Cycle 0)",
    x = "Sample 2",
    y = "Sample 1"
  )
