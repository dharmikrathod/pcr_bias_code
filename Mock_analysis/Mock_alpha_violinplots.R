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

# Same shape/logic as LambdaX_zero
mockLambdaX_0_df <- mockLambdaX_0 %>%
  fido::gather_array(val, coord, sample, iter) %>%
  group_by(coord, sample) %>%
  mutate(coord = names_categories(fit)[coord]) %>%
  ungroup() %>%
  pivot_wider(names_from = iter, values_from = val)

mockLambdaX_35_df <- mockLambdaX_35 %>%
  fido::gather_array(val, coord, sample, iter) %>%
  group_by(coord, sample) %>%
  mutate(coord = names_categories(fit)[coord]) %>%
  ungroup() %>%
  pivot_wider(names_from = iter, values_from = val)
mock_norm_0 <- mockLambdaX_0_df %>%
  pivot_longer(cols = -c(coord, sample), names_to = "iteration", values_to = "proportion") %>%
  group_by(sample, iteration) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  ungroup()

mock_norm_35 <- mockLambdaX_35_df %>%
  pivot_longer(cols = -c(coord, sample), names_to = "iteration", values_to = "proportion") %>%
  group_by(sample, iteration) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  ungroup()


mock_shannon_0 <- mock_norm_0 %>%
  group_by(sample, iteration) %>%
  summarize(shannon_index = -sum(proportion * log(proportion), na.rm = TRUE), .groups = "drop") %>%
  mutate(cycle_name = 0, iteration = as.numeric(iteration))

mock_shannon_35 <- mock_norm_35 %>%
  group_by(sample, iteration) %>%
  summarize(shannon_index = -sum(proportion * log(proportion), na.rm = TRUE), .groups = "drop") %>%
  mutate(cycle_name = 35, iteration = as.numeric(iteration))

mock_shannon_combined <- bind_rows(mock_shannon_0, mock_shannon_35)

gini_index <- function(p) {
  n <- length(p)
  if (n == 0) return(NA)
  p <- sort(p)
  G <- sum((2 * seq_along(p) - n - 1) * p) / (n - 1)
  return(G)
}

mock_gini_0 <- mock_norm_0 %>%
  group_by(sample, iteration) %>%
  summarize(gini_index = gini_index(proportion), .groups = "drop") %>%
  mutate(cycle_name = 0, iteration = as.numeric(iteration))

mock_gini_35 <- mock_norm_35 %>%
  group_by(sample, iteration) %>%
  summarize(gini_index = gini_index(proportion), .groups = "drop") %>%
  mutate(cycle_name = 35, iteration = as.numeric(iteration))

mock_gini_combined <- bind_rows(mock_gini_0, mock_gini_35)
gini_index <- function(p) {
  n <- length(p)
  if (n == 0) return(NA)
  p <- sort(p)
  G <- sum((2 * seq_along(p) - n - 1) * p) / (n - 1)
  return(G)
}

mock_gini_0 <- mock_norm_0 %>%
  group_by(sample, iteration) %>%
  summarize(gini_index = gini_index(proportion), .groups = "drop") %>%
  mutate(cycle_name = 0, iteration = as.numeric(iteration))

mock_gini_35 <- mock_norm_35 %>%
  group_by(sample, iteration) %>%
  summarize(gini_index = gini_index(proportion), .groups = "drop") %>%
  mutate(cycle_name = 35, iteration = as.numeric(iteration))

mock_gini_combined <- bind_rows(mock_gini_0, mock_gini_35)


mock_simpson_0 <- mock_norm_0 %>%
  group_by(sample, iteration) %>%
  summarize(simpson_index = 1 - sum(proportion^2, na.rm = TRUE), .groups = "drop") %>%
  mutate(cycle_name = 0, iteration = as.numeric(iteration))

mock_simpson_35 <- mock_norm_35 %>%
  group_by(sample, iteration) %>%
  summarize(simpson_index = 1 - sum(proportion^2, na.rm = TRUE), .groups = "drop") %>%
  mutate(cycle_name = 35, iteration = as.numeric(iteration))

mock_simpsons_combined <- bind_rows(mock_simpson_0, mock_simpson_35)

mock_aitchison_0 <- mock_norm_0 %>%
  group_by(sample, iteration) %>%
  summarize(
    aitchison_norm = {
      p <- proportion
      log_p <- log(p + 1e-10)  # add pseudo-count to avoid log(0)
      clr_p <- log_p - mean(log_p)
      sqrt(sum(clr_p^2))
    },
    .groups = "drop"
  ) %>%
  mutate(
    cycle_name = 0,
    iteration = as.numeric(iteration)
  )

mock_aitchison_35 <- mock_norm_35 %>%
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
    cycle_name = 35,
    iteration = as.numeric(iteration)
  )

mock_aitchison_combined <- bind_rows(mock_aitchison_0, mock_aitchison_35)


library(patchwork)
  
p5 <- ggplot(mock_simpsons_combined %>% filter(sample != 7),
             aes(x = factor(sample), y = simpson_index, fill = factor(cycle_name))) +
  geom_violin(
    trim = FALSE,
    alpha = 0.5,
    scale = "width",
    width = 0.7     # Wider violins
    #adjust = 2       # Compress tails vertically
  ) +
  labs(x = "Samples", y = "Simpson", fill = NULL) +
  theme_minimal() +
  scale_fill_manual(
    values = c("0" = "white", "35" = "black"),
    labels = c("0" = "Cycle 0", "35" = "Cycle 35")
  ) +
  coord_cartesian(ylim = c(0.6, 1)) +  # Stretch y-axis
  theme(
    axis.text.y = element_text(size = 20), 
    text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none'
  )

  p6<- ggplot(mock_shannon_combined%>% filter(sample!=7), aes(x = factor(sample), y = shannon_index, fill = factor(cycle_name))) +
    geom_violin(
      trim = FALSE,
      alpha = 0.5,
      scale = "width",
      width = 0.7     # Wider violins
      #adjust = 2       # Compress tails vertically
    ) +
    labs(x = "Samples", y = "Simpson", fill = NULL) +
    theme_minimal() +
    scale_fill_manual(
      values = c("0" = "white", "35" = "black"),
      labels = c("0" = "Cycle 0", "35" = "Cycle 35")
    ) +
    coord_cartesian(ylim = c(1.5, 2.5)) +  # Stretch y-axis
    theme(
      axis.text.y = element_text(size = 20), 
      text = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = 'none'
    )
  
  p7 <- ggplot(mock_gini_combined%>% filter(sample!=7), aes(x = factor(sample), y = gini_index, fill = factor(cycle_name))) +
    geom_violin(
      trim = FALSE,
      alpha = 0.5,
      scale = "width",
      width = 0.7     # Wider violins
      #adjust = 2       # Compress tails vertically
    ) +
    labs(x = "Samples", y = "Simpson", fill = NULL) +
    theme_minimal() +
    scale_fill_manual(
      values = c("0" = "white", "35" = "black"),
      labels = c("0" = "Cycle 0", "35" = "Cycle 35")
    ) +
    coord_cartesian(ylim = c(0, 1)) +  # Stretch y-axis
    theme(
      axis.text.y = element_text(size = 20), 
      text = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = 'none'
    )
  
  
  p8 <- ggplot(mock_aitchison_combined%>% filter(sample!=7), aes(x = factor(sample), y = aitchison_norm, fill = factor(cycle_name))) +
    geom_violin(
      trim = FALSE,
      alpha = 0.5,
      scale = "width",
      width = 0.7     # Wider violins
      #adjust = 2       # Compress tails vertically
    ) +
    labs(x = "Samples", y = "Simpson", fill = NULL) +
    theme_minimal() +
    scale_fill_manual(
      values = c("0" = "white", "35" = "black"),
      labels = c("0" = "Cycle 0", "35" = "Cycle 35")
    ) +
    coord_cartesian(ylim = c(1, 7)) +  # Stretch y-axis
    theme(
      axis.text.y = element_text(size = 20), 
      text = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = 'none'
    )
  


#row1<- (p5+p1)
#row2<- (p6+p2)
#row3<- (p7+p3)
#row4<- (p8+p4)
#library(patchwork)
#(p5 +p1) + 
 # plot_layout(widths = c(length(unique(mock_simpsons_combined$sample)),length(unique(LambdaX_simpsons_combined$sample)))) &
 # theme(axis.text.x = element_blank(), axis.title.x = element_blank())
#common_levels <- union(levels(mock_simpsons_combined$sample), levels(LambdaX_simpsons_combined$sample))

