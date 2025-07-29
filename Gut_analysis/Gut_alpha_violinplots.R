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



focus.names <- c("clr_seq_9" = "Parasutterella", 
                 "clr_seq_48" = "Bifidobacterium", 
                 "clr_seq_1" = "Bacteroides", 
                 "clr_seq_192" = "Ruminococcus", 
                 "clr_seq_204" = "Coprococcus_1", 
                 "clr_seq_193" = "Holdemania")

LambdaX_focused <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0, 0, 0, 0),
                                              c(1, 0, 0, 0, 35, 0, 0, 0)), 
                           pars="LambdaX", use_names=FALSE)
bias <- clrInv_array(LambdaX_focused, 1)
bias <- bias[,2,]/bias[,1,]
bias <- log2(bias)

bias <- bias %>% 
  gather_array(val, coord, iter) %>% 
  group_by(coord) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(coord = coord.names[coord])

p <- bias %>% 
  mutate(sig = ifelse(sign(p2.5)==sign(p97.5), TRUE, FALSE)) %>% 
  arrange(mean) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(rank2 = ifelse((sig & (rank > 66)) | (sig & (rank < 4)), coord, "")) %>% 
  mutate(coord = factor(coord, levels = coord[order(mean)])) %>% 
  mutate(rank2 = focus.names[rank2]) %>% 
  ggplot(aes(x = coord, y = mean)) +
  geom_hline(yintercept = 0) +
  #geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = -2, linetype = "dashed", color = "red") +
  geom_pointrange(aes(ymin = p2.5, ymax = p97.5, color = sig)) +
  geom_label_repel(aes(label = rank2), xlim = c(6, 65), size = 6) +
  coord_cartesian(ylim = c(-4, 2.5)) +
  theme_minimal() +
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 16),         # y-axis tick labels
    axis.title.y = element_text(size = 16),        # y-axis title
    panel.grid.minor = element_blank(), 
    panel.grid.major.x = element_blank(), 
    legend.position = "none"
  )+
  scale_color_manual(values = c("grey", "black")) +
  ylab(expression(log[2]~over(Proportion~~at~~Cycle~~35, Proportion~~at~~Cycle~~0)))


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
LambdaX_shannon_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(iteration = as.numeric(iteration),
         cycle_name = 0  ) %>% # Convert iteration to numeric for proper sorting
  arrange(sample, iteration) # Sort by sample and iteration in ascending order

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
LambdaX_shannon_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(iteration = as.numeric(iteration),
         cycle_name = 35  ) %>% # Convert iteration to numeric for proper sorting
  arrange(sample, iteration) # Sort by sample and iteration in ascending order




# Combine the Shannon indices for 0th and 35th cycles
LambdaX_shannon_combined <- bind_rows(LambdaX_shannon_0, LambdaX_shannon_35)

# Create the violin plot




# Compute Simpson's index for 0th cycle
gini_index <- function(p) {
  n <- length(p)
  if (n == 0) return(NA)
  p <- sort(p)
  G <- sum((2 * seq_along(p) - n - 1) * p) / (n - 1)
  return(G)
}

# Compute Gini index for cycle 0
LambdaX_gini_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    gini_index = gini_index(proportion),
    .groups = "drop"
  ) %>%
  mutate(
    iteration = as.numeric(iteration),
    cycle_name = 0
  ) %>%
  arrange(sample, iteration)

# Compute Gini index for cycle 35
LambdaX_gini_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    gini_index = gini_index(proportion),
    .groups = "drop"
  ) %>%
  mutate(
    iteration = as.numeric(iteration),
    cycle_name = 35
  ) %>%
  arrange(sample, iteration)

# Combine both
LambdaX_gini_combined <- bind_rows(LambdaX_gini_0, LambdaX_gini_35)

# Plot Gini index

# Compute Simpson's index for 0th cycle
LambdaX_simpsons_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson_index = 1 - sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    iteration = as.numeric(iteration),
    cycle_name = 0
  ) %>%
  arrange(sample, iteration)

# Compute Simpson's index for 35th cycle
LambdaX_simpsons_35 <- LambdaX_normalized_35 %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson_index = 1 - sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    iteration = as.numeric(iteration),
    cycle_name = 35
  ) %>%
  arrange(sample, iteration)


LambdaX_simpsons_combined <- bind_rows(LambdaX_simpsons_0, LambdaX_simpsons_35)

LambdaX_aitchison_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    aitchison_norm = {
      p <- proportion
      log_p <- log(p + 1e-10)  # Add small value to avoid log(0)
      clr_p <- log_p - mean(log_p)
      sqrt(sum(clr_p^2))
    },
    .groups = "drop"
  ) %>%
  mutate(
    iteration = as.numeric(iteration),
    cycle_name = 0
  ) %>%
  arrange(sample, iteration)


LambdaX_aitchison_35 <- LambdaX_normalized_35 %>%
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
  mutate(
    iteration = as.numeric(iteration),
    cycle_name = 35
  ) %>%
  arrange(sample, iteration)

LambdaX_aitchison_combined <- bind_rows(LambdaX_aitchison_0, LambdaX_aitchison_35)

# Reorder the sample factor for all gut-related data frames
LambdaX_simpsons_combined$sample <- factor(LambdaX_simpsons_combined$sample,
                                           levels = c("4", "2", "3", "1"))

LambdaX_shannon_combined$sample <- factor(LambdaX_shannon_combined$sample,
                                          levels = c("4", "2", "3", "1"))

LambdaX_gini_combined$sample <- factor(LambdaX_gini_combined$sample,
                                       levels = c("4", "2", "3", "1"))

LambdaX_aitchison_combined$sample <- factor(LambdaX_aitchison_combined$sample,
                                            levels = c("4", "3", "2", "1"))



library(ggplot2)
library(patchwork)

# Simpson index plot
p1 <- ggplot(LambdaX_simpsons_combined, aes(x = factor(sample), y = simpson_index, fill = factor(cycle_name))) +
  geom_violin(trim = FALSE, alpha = 0.5, scale = "width", width = 0.7)+
  labs(
    #title = 'A',
    x = "Samples",
    y = "Simpson's Index",
    fill = "Cycle Name"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "white", "35" = "black")) +
  theme(
    axis.text.y = element_text(size = 20), 
    text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none'

  )

# Shannon index plot
p2 <- ggplot(LambdaX_shannon_combined, aes(x = factor(sample), y = shannon_index, fill = factor(cycle_name))) +
  geom_violin(trim = FALSE, alpha = 0.5, scale = "width", width = 0.7)+
  labs(
    #title = 'A',
    x = "Samples",
    y = "Simpson's Index",
    fill = "Cycle Name"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "white", "35" = "black")) +
  theme(
    axis.text.y = element_text(size = 20), 
    text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none'
    
  )

p3<- ggplot(LambdaX_gini_combined, aes(x = factor(sample), y = gini_index, fill = factor(cycle_name))) +
  geom_violin(
    trim = FALSE,
    alpha = 0.5,
    scale = 'width',
    width = 0.7
  )+
  labs(
    #title = 'A',
    x = "Samples",
    y = "Simpson's Index",
    fill = "Cycle Name"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "white", "35" = "black")) +
  theme(
    axis.text.y = element_text(size = 20), 
    text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none'
    
  )

p4<-ggplot(LambdaX_aitchison_combined, aes(x = factor(sample), y = aitchison_norm, fill = factor(cycle_name))) +
  geom_violin(trim = FALSE, alpha = 0.5, scale = "width", width = 0.7)+
  labs(
    #title = 'A',
    x = "Samples",
    y = "Simpson's Index",
    fill = "Cycle Name"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "white", "35" = "black")) +
  theme(
    axis.text.y = element_text(size = 20), 
    text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none'
    
  )

# Combine side by side
(p1 / p2/p3/p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")











#not for the paper


LambdaXJoint <- predict(fit, newdata=cbind(
                        ## Cycle 0
                        c(1, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 1, 0, 0, 0, 0, 0), 
                        c(0, 0, 0, 1, 0, 0, 0, 0),
                        ## cycle 35
                        c(1, 0, 0, 0, 35, 0, 0, 0),
                        c(0, 1, 0, 0, 35, 0, 0, 0),
                        c(0, 0, 1, 0, 35, 0, 0, 0), 
                        c(0, 0, 0, 1, 35, 0, 0, 0)), 
                        pars="LambdaX", use_names=FALSE)


LambdaXJoint_table <- LambdaXJoint %>%
  fido::gather_array(val, coord, sample, iter) %>%
  group_by(coord, sample)%>%
  mutate(
    pcr_cycle = ifelse(sample <= 4, 0, 35),  # Assign cycle number (0 for first 4 samples, 35 for next 4 samples)
    Donor = c("16S8921", "908", "Ai96", "J1122526T", 
              as.character(sample_data(ps)$BiasDonorName))[sample],
    coord = coord.names[coord]
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = iter, values_from = val) %>%
  mutate(cycle_name = ifelse(pcr_cycle == 0, "0", "35"))


LambdaXJoint_normalized <- LambdaXJoint_table %>%
  # Convert 'coord' and 'sample' to character to avoid type mismatch
  mutate(
    coord = as.character(coord),
    sample = as.character(sample),  # Convert sample to character
    pcr_cycle = as.character(pcr_cycle)  # Convert pcr_cycle to character for consistency
  ) %>%
  # Now pivot the data, specifying the range of columns for iterations
  pivot_longer(cols = `1`:`2000`, names_to = "iteration", values_to = "proportion") %>%
  # Ensure iteration is numeric for proper sorting
  mutate(iteration = as.numeric(iteration)) %>%
  group_by(sample, iteration, cycle_name) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  ungroup()

# Compute Shannon's index for each sample and posterior iteration
LambdaXJoint_shannon <- LambdaXJoint_normalized %>%
  group_by(sample, iteration, cycle_name) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sample, cycle_name, iteration)  # Sort by sample, cycle_name, and iteration



sample_1_data <- LambdaXJoint_shannon %>%
  filter(sample == 1)

sample_5_data <- LambdaXJoint_shannon %>%
  filter(sample == 5)



sample_1_data <- sample_1_data[order(sample_1_data$iteration), ]
sample_5_data <- sample_5_data[order(sample_5_data$iteration), ]

# Merge the two dataframes by `iteration`
merged_data <- merge(sample_1_data, sample_5_data, by = "iteration", suffixes = c("_1", "_5"))

# Calculate the difference in Shannon index
merged_data_1_diff <- merged_data %>%
  mutate(shannon_diff = shannon_index_1 - shannon_index_5)




sample_4_data <- LambdaXJoint_shannon %>%
  filter(sample == 4)

sample_8_data <- LambdaXJoint_shannon %>%
  filter(sample == 8)



sample_4_data <- sample_4_data[order(sample_4_data$iteration), ]
sample_8_data <- sample_8_data[order(sample_8_data$iteration), ]

# Merge the two dataframes by `iteration`
merged_data_2 <- merge(sample_4_data, sample_8_data, by = "iteration", suffixes = c("_4", "_8"))

# Calculate the difference in Shannon index
merged_data_4_diff <- merged_data_2 %>%
  mutate(shannon_diff = shannon_index_4 - shannon_index_8)


merged_data_1_diff <- merged_data_1_diff %>% mutate(group = "Dataset 1")
merged_data_4_diff <- merged_data_4_diff %>% mutate(group = "Dataset 4")

# Combine the datasets
combined_data <- bind_rows(merged_data_1_diff, merged_data_4_diff)

# Plot the combined density
ggplot(combined_data, aes(x = shannon_diff, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Kernel Density Plot of Shannon Index Differences",
    x = "Shannon Index Difference",
    y = "Density"
  ) +
  scale_fill_manual(values = c("blue", "red"), name = "Dataset") +
  theme_minimal()
















