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



fit<- to_proportions(fit)

library(gridExtra)
library(vegan)
bray_curtis<- function(indices){
  
  vector10<-c(0, 0, 0, 0, 0, 0, 0, 0)
  vector20<-c(0, 0, 0, 0, 0, 0, 0, 0)
  vector11<-c(0, 0, 0, 0, 35, 0, 0, 0)
  vector21<-c(0, 0, 0, 0, 35, 0, 0, 0)
  
  
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
  
  
  return(c(bray_curtis_values1,bray_curtis_values2))
 
}

result12<- bray_curtis(c(1,2))
result13<- bray_curtis(c(1,3))
result14<- bray_curtis(c(1,4))
result23<- bray_curtis(c(2,3))
result24<- bray_curtis(c(2,4))
result34<- bray_curtis(c(3,4))


bc<- data.frame("1-2" = result12, "1-3" = result13, "1-4" = result14,
                "2-3" = result23, "2-4" = result24, "3-4" = result34, 
                "cycles" = rep(c("0th cycle","35th cycles"),each = length(result12)/2))

bc_long <- bc %>%
  pivot_longer(cols = 1:6,   # Select columns 1-6 (ignoring the "cycles" column)
               names_to = "comparison",  # Name for the variable
               values_to = "value")     # Name for the values


bc_0<-bc[1:2000,]
bc_35<-bc[2001:4000,]

bc_diff<-bc_0[,-c(7)] - bc_35[,-c(7)]

df_long <- bc_diff %>%
  pivot_longer(cols = everything(), names_to = "Pair", values_to = "Bias")

# Summarize posterior median and 95% CI for each pair
summary_bias <- df_long %>%
  group_by(Pair) %>%
  summarize(
    p2.5 = quantile(Bias, 0.025),
    p50 = quantile(Bias, 0.5),
    p97.5 = quantile(Bias, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    Sample1 = sub("-.*", "", Pair),
    Sample2 = sub(".*-", "", Pair)
  )

# Set sample order
sample_levels <- c("1", "2", "3", "4")

# Filter to lower triangle only
summary_bias_lower <- summary_bias %>%
  mutate(
    Sample1 = factor(Sample1, levels = sample_levels),
    Sample2 = factor(Sample2, levels = sample_levels)
  ) %>%
  filter(as.numeric(as.character(Sample1)) < as.numeric(as.character(Sample2)))
summary_bias <- summary_bias %>%
  mutate(
    Sample1_num = as.numeric(sub("X(\\d)\\..*", "\\1", Pair)),
    Sample2_num = as.numeric(sub("X\\d\\.(\\d)", "\\1", Pair))
  )

summary_bias_lower <- summary_bias %>%
  filter(Sample1_num < Sample2_num)
sample_levels <- c("1", "2", "3", "4")

summary_bias_lower <- summary_bias_lower %>%
  mutate(
    Sample1 = factor(Sample1_num, levels = sample_levels),
    Sample2 = factor(Sample2_num, levels = sample_levels)
  )
# Plot: tile + text with median and 95% CI
ggplot(summary_bias_lower, aes(x = Sample2, y = Sample1, fill = p50)) +
  geom_tile() +
  geom_text(aes(label = signif(p50, 2)),
            color = "white", size = 7.5, nudge_y = 0.15) +  # made text larger
  geom_text(aes(label = paste0("(", signif(p2.5, 2), "-", signif(p97.5, 2), ")")),
            color = "white", size = 6.5, nudge_y = -0.15) +  # made CI text larger
  scale_fill_gradient(low = "black", high = "#56B4E9", name = "Bray-Curtis bias") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  labs(title = "Posterior Change in Bray-Curtis (Cycle 35 − Cycle 0)")



