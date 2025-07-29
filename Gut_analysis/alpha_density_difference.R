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

Shannnon_diff<- LambdaX_shannon_wide_0 - LambdaX_shannon_wide_35

library(tidyr)
library(dplyr)
library(ggplot2)
colnames(Shannnon_diff)[-1] <- as.character(colnames(Shannnon_diff)[-1])


okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)
# Assume your data frame is named df
df_long <- Shannnon_diff[,c('iteration', '1', '3')] %>%
  pivot_longer(
    cols = -iteration,  # everything except 'iteration'
    names_to = "sample",
    values_to = "shannon_index"
  )
ggplot(df_long, aes(x = shannon_index, color = sample, fill = sample)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Bias in Shannon index",
    x = "Bias in Shannon Index After 35 Cycles of PCR",
    y = 'Density',
    fill = "Sample",
    color = "Sample"
  ) +
  scale_color_manual(values = rep('black',2)) +
  scale_fill_manual(values = rep('grey', 2)) +
  theme_minimal()+
  theme(
    #axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.text.y = element_blank(),
    #axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none'# Remove y-axis tick labels
  )

#for simpsons
# Simpson index: 1 - sum(p^2)
LambdaX_simpson_wide_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson_index = 1 - sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = simpson_index
  ) %>%
  arrange(iteration)

LambdaX_simpson_wide_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson_index = 1 - sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = simpson_index
  ) %>%
  arrange(iteration)

# Difference in Simpson indices
Simpson_diff <- LambdaX_simpson_wide_0 - LambdaX_simpson_wide_35

# Rename columns for clarity
colnames(Simpson_diff)[-1] <- as.character(colnames(Simpson_diff)[-1])

# Reshape for plotting
df_long1 <- Simpson_diff[,c('iteration', '1', '3')] %>%
  pivot_longer(
    cols = -iteration,
    names_to = "sample",
    values_to = "simpson_index"
  )

# Plot
ggplot(df_long1, aes(x = simpson_index, color = sample, fill = sample)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Bias in Simpson index",
    x = "Bias in Simpson Index",
    y = "Density",
    fill = "Sample",
    color = "Sample"
  ) +
  scale_color_manual(values = rep("black", 2)) +
  scale_fill_manual(values = rep('grey', 2)) +
  theme_minimal() +
  theme(
    #axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
   # axis.title.y = element_blank(),
    legend.position = 'none'# Remove y-axis tick labels
  )

LambdaX_gini_wide_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    gini_index = gini_index(proportion),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = gini_index
  ) %>%
  arrange(iteration)

LambdaX_gini_wide_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    gini_index = gini_index(proportion),
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = gini_index
  ) %>%
  arrange(iteration)
# Difference in Gini indices
Gini_diff <- LambdaX_gini_wide_0 - LambdaX_gini_wide_35

# Keep iteration column

# Reshape for selected samples (swap 1 and 4 if needed)
df_long_gini <- Gini_diff[, c('iteration', '1', '3')] %>%
  pivot_longer(
    cols = -iteration,
    names_to = "sample",
    values_to = "gini_index"
  )
ggplot(df_long_gini, aes(x = gini_index, color = sample, fill = sample)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Bias in Gini Index",
    x = "Bias in Gini Index",
    y = "Density",
    fill = "Sample",
    color = "Sample"
  ) +
  scale_color_manual(values = rep("black", 2)) +
  scale_fill_manual(values = rep("grey", 2)) +
  theme_minimal() +
  theme(
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    legend.position = 'none'
  )

LambdaX_aitchison_wide_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    aitchison_norm = {
      p <- proportion
      log_p <- log(p + 1e-10)
      clr_p <- log_p - mean(log_p)
      sqrt(sum(clr_p^2))
    },
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = aitchison_norm
  ) %>%
  arrange(iteration)
LambdaX_aitchison_wide_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    aitchison_norm = {
      p <- proportion
      log_p <- log(p + 1e-10)
      clr_p <- log_p - mean(log_p)
      sqrt(sum(clr_p^2))
    },
    .groups = "drop"
  ) %>%
  mutate(iteration = as.numeric(iteration)) %>%
  pivot_wider(
    names_from = sample,
    values_from = aitchison_norm
  ) %>%
  arrange(iteration)

Aitchison_diff <- LambdaX_aitchison_wide_0 - LambdaX_aitchison_wide_35

# Reshape for selected samples (e.g., "1" and "3" like before)
df_long_aitchison <- Aitchison_diff[, c('iteration', '1', '2')] %>%
  pivot_longer(
    cols = -iteration,
    names_to = "sample",
    values_to = "aitchison_index"
  )
ggplot(df_long_aitchison, aes(x = aitchison_index, color = sample, fill = sample)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Bias in Aitchison Norm",
    x = "Bias in Aitchison Norm",
    y = "Density",
    fill = "Sample",
    color = "Sample"
  ) +
  scale_color_manual(values = rep("black", 2)) +
  scale_fill_manual(values = c('grey', 'grey')) +
  theme_minimal() +
  theme(
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
  
  )



