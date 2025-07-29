library(phyloseq)
library(ggtern)
library(driver)
library(compositions)
library(tidyverse)
library(ape)

# ------------ Helper functions ------------

# Weighted UniFrac distance between two compositions
compute_unifrac <- function(x, y, tree) {
  otumat <- rbind(x, y)
  rownames(otumat) <- c("before", "after")
  colnames(otumat) <- tree$tip.label
  ps <- phyloseq(otu_table(otumat, taxa_are_rows = FALSE), phy_tree(tree))
  dist <- UniFrac(ps, weighted = TRUE, normalized = TRUE)
  as.numeric(dist)
}

# ------------ Simulated Tree (3-taxa rooted) ------------

tree_txt <- "((tax1:1,tax2:1):1,tax3:2);"
tree <- read.tree(text = tree_txt)

# ------------ Bias setup ------------

bias <- c(1, 0.85, 1)                     # PCR bias vector
bias <- as.numeric(miniclo(bias))
bias_ilr <- ilr(bias)

x1 <- c(5,3,1)                    # Reference composition
x1 <- as.numeric(miniclo(x1))
x1_ilr <- ilr(x1)

uniform_simplex_grid <- function(n) {
  pts <- list()
  for (i in 0:n) {
    for (j in 0:(n - i)) {
      k <- n - i - j
      pts[[length(pts) + 1]] <- c(i, j, k) / n
    }
  }
  do.call(rbind, pts)
}
# ------------ Grid setup in ILR space ------------\
n_grid <- 70  # larger n = finer grid
x.prop <- uniform_simplex_grid(n_grid)
epsilon <- 1e-8
x2_prop <- x.prop + epsilon
x2_prop <- x2_prop / rowSums(x2_prop)

x2_ilr <- ilr(x2_prop) %>% as.matrix()
y2_ilr <- sweep(x2_ilr, 2, 35 * bias_ilr, FUN = "+")

# ------------ Back-transform to compositions ------------

x1_prop <- ilrInv(x1_ilr)
y1_prop <- ilrInv(x1_ilr + 35 * bias_ilr)


y2_prop <- ilrInv(y2_ilr)

# ------------ Compute UniFrac bias ------------

unifrac1 <- numeric(nrow(x2_ilr))
unifrac2 <- numeric(nrow(x2_ilr))

for (i in 1:nrow(x2_prop)) {
  unifrac1[i] <- compute_unifrac(x1_prop, x2_prop[i,], tree)
  unifrac2[i] <- compute_unifrac(y1_prop, y2_prop[i,], tree)
}

delta_unifrac <- unifrac1 - unifrac2

# ------------ Plotting ------------

df_all <- data.frame(
  v1 = x2_prop[, 1],
  v2 = x2_prop[, 2],
  v3 = x2_prop[, 3],
  val = delta_unifrac
)

x1.plot <- as.data.frame(t(x1_prop))
colnames(x1.plot) <- c("v1", "v2", "v3")

bias.plot <- as.data.frame(t(ilrInv(bias_ilr)))
colnames(bias.plot) <- c("v1", "v2", "v3")

zero.plot <- as.data.frame(t(ilrInv(c(0, 0))))
colnames(zero.plot) <- c("v1", "v2", "v3")

cbPalette <- c(
  "#2166AC", "white","#B2182B"
)
# Use Okabe-Ito palette or other colorblind safe palette
ggtern(df_all, aes(v1, v2, v3, color = val)) +
  geom_point() +
  geom_point(data = bias.plot, aes(v1, v2, v3),
             shape = 21, fill = "#D55E00", color = "black", size = 5.5, stroke = 1.4) +  # orange
  geom_point(data = zero.plot, aes(v1, v2, v3),
             shape = 21, fill = "#009E73", color = "black", size = 4, stroke = 1.4) +    # green
  geom_point(data = x1.plot, shape = 21, fill = "#999", color = "black", size = 4.8, stroke = 1.6) +
  scale_color_gradientn(
    colors = cbPalette,  # should be ordered like: low -> mid -> high (e.g., blue, white, red)
    values = scales::rescale(c(min(df_all$val), 0, max(df_all$val))),  # force zero to midpoint
    name = "Unifrac bias",
    limits = c(min(df_all$val), max(df_all$val))
  )  +
  labs(
    T = "TAXA 1",
    L = "TAXA 2",
    R = "TAXA 3"
  ) +
  scale_size_continuous(range = c(1, 6), name = "UniFrac bias") +
  theme_bw(base_size = 15) +
  theme(
    plot.margin = unit(c(-1, -2, -1, -2), "cm"),
    tern.axis.title.T = element_text(margin = margin(b = 50), size = 18),
    tern.axis.title.L = element_text(margin = margin(t = 50),size = 18),
    tern.axis.title.R = element_text(margin = margin(t = 50),size = 18),
    tern.axis.text = element_text(size = 18),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 20)
    
  )
  
