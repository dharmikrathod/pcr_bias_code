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
