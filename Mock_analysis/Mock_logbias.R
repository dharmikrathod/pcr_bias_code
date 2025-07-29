library(tidyverse)
library(fido)
library(driver)
library(phyloseq)
library(Biostrings)
library(readxl)
library(ggstance)
library(gt)
library(mvrsquared)

# Load Data and Preprocess------------------------------------------------------

# Picogreen 
pg <- "2020-03-05_PCRBias_qPCR_Picogreen.xlsx"
pg <- readxl::read_xlsx(pg)

# Dilution Ratios
dr <- "DNA dilutions.xlsx"
dr <- readxl::read_xlsx(dr)
dr[,2]/300*pg$PicogreenDNAConc
# So it does look like amounts were properly pooled. 


# Read in Mock community ratios
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


# read in phyloseq
ps <- readRDS("phyloseq.rds")

# Filter low abundance samples
total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>1000, ps)
sum(sample_sums(ps))/total.reads

# read in seq clusters -- created by sequencing isolates individually
load("seq_clusters.RData")

# Map seq to isolate clusters
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
seq_dist <- stringDist(combined_seq, method = "levenshtein")
n1 <- length(seq_dict)
n2 <- length(refseq(ps))
seq_dist <- as.matrix(seq_dist)
seq_dist <- seq_dist[1:n1, (n1+1):(n1+n2)]

# Match based on minimum distance
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

# isolate to just calibration samples
calibL <- str_detect(colnames(Y), "^cycle")
Ycalib <- Y[,calibL]
Xcalib <- t(model.matrix(~ cycle_num + machine, data=d[calibL,]))

sigma2 <- seq(.5, 10, by=.5)
lml <- rep(NA, length(sigma2))
for (i in seq_along(sigma2)){
  fit <- pibble(Y=Ycalib, X=Xcalib, Gamma = sigma2[i]*diag(nrow(Xcalib)))
  lml[i] <- fit$logMarginalLikelihood
}
plot(sigma2, lml)

# pick sigma2=10 based on maximizing log-marginal likelihood
fit <- pibble(Y=Ycalib, X=Xcalib, Gamma = 10*diag(nrow(Xcalib)))


# transform to CLR --------------------------------------------------------

fit <- to_clr(fit)

# Magnitude of Bias -------------------------------------------------------

LambdaX_focused <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0),
                                              c(1, 35, 0, 0, 0)), 
                           pars="LambdaX", use_names=FALSE)
bias <- clrInv_array(LambdaX_focused, 1)
bias <- bias[,2,]/bias[,1,]
bias <- log2(bias)

bias %>% 
  gather_array(val, coord, iter) %>% 
  mutate(val = abs(val)) %>% 
  group_by(iter) %>% 
  summarise(val = exp(mean(val))) %>% 
  ungroup() %>% 
  summarise_posterior(val)


bias <- bias %>% 
  gather_array(val, coord, iter) %>% 
  group_by(coord) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(coord = names_categories(fit)[coord])

p <- bias %>% 
  mutate(sig = ifelse(sign(p2.5)==sign(p97.5), TRUE, FALSE)) %>% 
  arrange(mean) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(coord = factor(coord, levels=coord[order(mean)])) %>% 
  ggplot(aes(x = coord, y=mean)) +
  geom_hline(yintercept=0) +
  geom_pointrange(aes(ymin=p2.5, ymax=p97.5, color=sig)) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "black")) +
  ylab(expression(log[2]~over(Proportion~~at~~Cycle~~35,Proportion~~at~~Cycle~~0))) +
  theme(
    axis.title.y = element_text(size = 17), 
    axis.text.y = element_text(size = 17),
    axis.text.x = element_text(size = 17, angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none'
  )
p

