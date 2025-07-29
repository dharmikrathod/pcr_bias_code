library(tidyverse)
library(fido)
library(driver)
library(phyloseq)
library(Biostrings)
library(readxl)
library(ggstance)
library(gt)
library(mvrsquared)
library(pso)  # Load PSO package
library(vegan) 

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
fit <- pibble(Y=Y, X=X, Gamma = 10*diag(nrow(X)))
fit <- to_proportions(fit)


newdata <- cbind(diag(11), diag(11))
newdata <- rbind(newdata, rep(c(0,35), each=11), matrix(0, 3, 22))
sn <- str_sub(rownames(X)[1:11], start=11)
rownames(newdata) <- rownames(X)
newdata_meta <- data.frame(sample = rep(sn, times=2),
                           cycle = rep(c(0, 35), each=11))

LambdaX <- predict(fit, newdata, pars="LambdaX")
dimnames(LambdaX)[[1]] <-  names_categories(fit)

LambdaX_tidy <- LambdaX %>%
  gather_array(val, taxa, sample, iter) %>%
  group_by(taxa, sample) %>%
  summarise_posterior(val) %>%
  ungroup() %>%
  mutate(PCR_Cycle = newdata_meta$cycle[sample],
         sample = newdata_meta$sample[sample],
         taxa = names_categories(fit)[taxa])


LambdaX_long <- LambdaX %>%
  gather_array(val, taxa, sample, iter) %>%
  mutate(PCR_Cycle = newdata_meta$cycle[sample],
         sample = newdata_meta$sample[sample],
         taxa = names_categories(fit)[taxa])

library(GA)  # Genetic algorithm package
library(vegan)

mean_matrix <- apply(LambdaX, c(1, 2), mean)

cycle_0_data<- t(mean_matrix[,1:11])
cycle_35_data<- t(mean_matrix[,12:22])

compute_delta_r2 <- function(group_assignment, cycle_0_data, cycle_35_data) {
  if (length(group_assignment) != 11) {
    stop("Error: group_assignment must have length 11")  
  }
  
  group <- factor(group_assignment)
  
  # Ensure at least two unique group values (0 and 1)
  if (length(unique(group)) < 2) {
    return(-Inf)  # Penalize invalid solutions
  }
  
  # Compute PERMANOVA for cycle 0
  dist_0 <- vegdist(cycle_0_data, method = "bray")
  permanova_0 <- adonis2(dist_0 ~ group)
  r2_0 <- permanova_0$R2[1]  
  
  # Compute PERMANOVA for cycle 35
  dist_35 <- vegdist(cycle_35_data, method = "bray")
  permanova_35 <- adonis2(dist_35 ~ group)
  r2_35 <- permanova_35$R2[1]  
  
  return(abs(r2_0 - r2_35))
}

# PSO Fitness Function
pso_function <- function(params) {
  if (length(params) != 11) {
    stop("Error: params must have length 11")  
  }
  
  group_assignment <- round(params)# Ensure binary assignment
  
  delta_r2 <- compute_delta_r2(group_assignment, cycle_0_data, cycle_35_data)
  
  return(delta_r2)  # Minimize negative R² difference to maximize the original R² difference
}


set.seed(1232)
ga_result <- ga(
  type = "binary",  # Binary representation for group assignments (0 or 1)
  fitness = pso_function,  # Fitness function to maximize R2 difference
  nBits = 11,  # We have 22 samples (11 for cycle 0, 11 for cycle 35)
  maxiter = 500,  # Number of generations to run
  popSize = 200,  # Population size
  pmutation = 0.3,  # Probability of mutation
  pcrossover = 0.7,  # Probability of crossover (between 0 and 1)
  elitism = 2, # Number of best solutions to keep
   run = 10# Number of generations without improvement before stopping
)
array(ga_result@solution)
install.packages("pso")  # Install pso package if not already installed
 # Load vegan for PERMANOVA

# Function to compute R² difference from PERMANOVA


# Run Particle Swarm Optimization
pso_result <- psoptim(
  par = c(rep(0, 5),rep(1,6)),  # Initial solution (11 elements, 0s)
  fn = pso_function,  # Fitness function
  lower = rep(0, 11),  # Lower bound (0)
  upper = rep(1, 11),  # Upper bound (1)
  control = list(
    maxit = 500,  # Maximum iterations
    s = 50,  # Number of particles in the swarm
    p = 0.5,  # Probability of mutation
    trace = 1  # Print optimization progress
  )
)

# Print best solution
print(round(pso_result$par))  # Best group assignment found
print(-pso_result$value)  # Maximum R² difference found




library(ape)

# Load the pruned tree
tree <- read.tree("mock_pruned.asv.nwk")


taxa_in_data <- rownames(mean_matrix)
taxa_in_tree <- tree$tip.label

# Check overlap
setequal(taxa_in_data, taxa_in_tree)  # should return TRUE


rownames(mean_matrix) <- tree$tip.label




library(phyloseq)

# Cycle 0: samples 1 to 11
abund_0 <- mean_matrix[, 1:11]

# Cycle 35: samples 12 to 22
abund_35 <- mean_matrix[, 12:22]

# Create phyloseq objects
ps_0 <- phyloseq(
  otu_table(abund_0, taxa_are_rows = TRUE),
  phy_tree(tree)
)

ps_35 <- phyloseq(
  otu_table(abund_35, taxa_are_rows = TRUE),
  phy_tree(tree)
)

compute_delta_r2_unifrac <- function(group_assignment, ps_0, ps_35) {
  if (length(group_assignment) != 11) stop("group_assignment must have length 11")
  
  group <- factor(round(group_assignment))
  if (length(unique(group)) < 2) return(-Inf)
  
  # Compute UniFrac distances (Weighted)
  dist_0 <- UniFrac(ps_0, weighted = TRUE, normalized = TRUE)
  dist_35 <- UniFrac(ps_35, weighted = TRUE, normalized = TRUE)
  
  # Run PERMANOVA
  adonis_0 <- adonis2(as.dist(dist_0) ~ group)
  adonis_35 <- adonis2(as.dist(dist_35) ~ group)
  
  r2_0 <- adonis_0$R2[1]
  r2_35 <- adonis_35$R2[1]
  
  return(-abs(r2_0 - r2_35))
}

pso_function_unifrac <- function(params) {
  group_assignment <- round(params)
  compute_delta_r2_unifrac(group_assignment, ps_0, ps_35)
}

library(pso)

set.seed(1234)

pso_result_unifrac <- psoptim(
  par = c(rep(0, 5), rep(1, 6)),  # Initial guess
  fn = pso_function_unifrac,
  lower = rep(0, 11),
  upper = rep(1, 11),
  control = list(
    maxit = 500,
    s = 50,
    trace = 1
  )
)

# Output best solution
print(round(pso_result_unifrac$par))     # Optimal group assignment
print(-pso_result_unifrac$value)         # Maximum ΔR² (PERMANOVA) found









































































LambdaX_tidy <- LambdaX %>%
  gather_array(val, taxa, sample, iter) %>%
  group_by(taxa, sample) %>%
  ungroup() %>%
  mutate(PCR_Cycle = newdata_meta$cycle[sample],
         sample = newdata_meta$sample[sample],
         taxa = names_categories(fit)[taxa])%>%
  ungroup() %>%
  pivot_wider(names_from = iter, values_from = val)



Lambda_0<- LambdaX_tidy[1:110,]
Lambda_35<- LambdaX_tidy[111:220,]
LambdaX_normalized <- Lambda_0 %>%
  pivot_longer(cols = `1`:`2000`, names_to = "iteration", values_to = "proportion") %>%
  group_by(sample, iteration) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  ungroup()



LambdaX_shannon_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    shannon_index = -sum(proportion * log(proportion), na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(iteration = as.numeric(iteration),
         cycle_name = 0  ) %>% # Convert iteration to numeric for proper sorting
  arrange(sample, iteration)

LambdaX_normalized_35 <- Lambda_35 %>%
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


LambdaX_shannon_combined <- bind_rows(LambdaX_shannon_0, LambdaX_shannon_35)

# Create the violin plot
ggplot(LambdaX_shannon_combined, aes(x = factor(sample), y = shannon_index, fill = factor(cycle_name))) +
  geom_violin(trim = FALSE) +
  labs(
    title = "Shannon Index Violin Plot",
    x = "Sample",
    y = "Shannon Index",
    fill = "Cycle Name"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "skyblue", "35" = "orange")) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )


LambdaX_simpsons_0 <- LambdaX_normalized %>%
  group_by(sample, iteration) %>%
  summarize(
    simpson_index = sum(proportion^2, na.rm = TRUE),
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
    simpson_index = sum(proportion^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    iteration = as.numeric(iteration),
    cycle_name = 35
  ) %>%
  arrange(sample, iteration)


LambdaX_simpsons_combined <- bind_rows(LambdaX_simpsons_0, LambdaX_simpsons_35)



ggplot(LambdaX_simpsons_combined, aes(x = factor(sample), y = simpson_index, fill = factor(cycle_name))) +
  geom_violin(trim = FALSE) +
  labs(
    title = "Simpson's Index Violin Plot",
    x = "Sample",
    y = "Simpson's Index",
    fill = "Cycle Name"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "skyblue", "35" = "orange")) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )



library(gridExtra)
bray_curtis<- function(indices){
  
  vector10<-c(0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0)
  vector20<-c(0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0)
  vector11<-c(0, 0, 0, 0, 0, 0, 0, 0,0,0,0,35,0,0,0)
  vector21<-c(0, 0, 0, 0, 0, 0, 0, 0,0,0,0,35,0,0,0)
  
  
  vector10[indices[1]]<-1
  vector20[indices[2]]<-1
  vector11[indices[1]]<-1
  vector21[indices[2]]<-1
  
  LambdaX_0 <- predict(fit, newdata=cbind(vector10,vector20),
                       pars="LambdaX", use_names=FALSE)
  LambdaX_35 <- predict(fit, newdata=cbind(vector11,vector21),
                        pars="LambdaX", use_names=FALSE)
  
  sample_10<- LambdaX_0[,1,]
  sample_20<- LambdaX_0[,2,]
  sample_11<- LambdaX_35[,1,]
  sample_21<- LambdaX_35[,2,]
  
  sample_10<- sample_10/colSums(sample_10)
  sample_20<- sample_20/colSums(sample_20)
  sample_11<- sample_11/colSums(sample_11)
  sample_21<- sample_21/colSums(sample_21)
  
  
  
  bray_curtis_values1 <- sapply(1:2000, function(i) {
    vegdist(rbind(sample_10[,i], sample_20[, i]), method = "bray")[1]
  })
  bray_curtis_values2 <- sapply(1:2000, function(i) {
    vegdist(rbind(sample_11[,i], sample_21[, i]), method = "bray")[1]
  })
  
  df<- data.frame(
    bray_curtis = c(bray_curtis_values1, bray_curtis_values2),
    cycles = rep(c("0th cycle","35th cycles"),each = length(bray_curtis_values1))
  )
  
  mean_diff<- mean(bray_curtis_values1 - bray_curtis_values2)
  
  ggplot(data = df, aes(x = bray_curtis, fill = cycles)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c('red', 'blue')) +
    labs(title = paste("Bray-Curtis Dissimilarity btw Sample", indices[1],'and',indices[2], "with mean diff of ",round(mean_diff,2)), 
         x = "Bray-Curtis", 
         y = "Density", 
         fill = "Cycles") +
    theme_minimal() +
    theme(plot.margin = margin(1, 5, 1, 1))  # Extra space on the right margin
  
  return(c(bray_curtis_values1,bray_curtis_values2))
  
}


sample_indices <- combn(11, 2, simplify = FALSE)

# **Run Bray-Curtis function for all pairs**
results_list <- lapply(sample_indices, bray_curtis)

# **Convert results to a data frame**
comparison_names <- sapply(sample_indices, function(x) paste(x[1], "-", x[2]))
bc <- do.call(cbind, results_list)
colnames(bc) <- comparison_names

bc <- data.frame(bc, cycles = rep(c("0th cycle", "35th cycle"), each = nrow(bc) / 2))

# **Transform into long format for plotting**
bc_long <- bc %>%
  pivot_longer(cols = -cycles,  # Select all but the 'cycles' column
               names_to = "comparison", 
               values_to = "value")

# **Plotting the results**
ggplot(bc_long, aes(x = comparison, y = value, fill = cycles)) +
  geom_violin(trim = FALSE) + 
  scale_fill_manual(values = c('red', 'blue')) +  
  labs(title = "Violin Plot of Bray-Curtis Dissimilarity by Cycle", 
       x = "Samples being compared", y = "Bray-Curtis Value", fill = "Cycle") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #




