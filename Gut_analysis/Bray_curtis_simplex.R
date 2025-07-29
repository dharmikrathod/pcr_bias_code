library(tidyverse)

library(ggtern)
library(driver) # jsilve24/driver
bray_curtis <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  stopifnot(length(x) == length(y), is.numeric(x), is.numeric(y))
  num <- sum(abs(x - y))
  denom <- sum(x + y)
  return(num / denom)
}

# ---- PCR bias vector ----
bias <- c(1, 0.85, 1)
bias <- as.numeric(miniclo(bias))
bias_ilr <- ilr(bias)

# ---- Reference composition ----
x1 <- c(1,2,3)
x1 <- as.numeric(miniclo(x1))
x1_ilr <- ilr(x1)

# ---- Uniform grid in simplex ----
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

n_grid <- 70
x.prop <- uniform_simplex_grid(n_grid)
epsilon <- 1e-8
x2_prop <- x.prop + epsilon
x2_prop <- x2_prop / rowSums(x2_prop)

# ---- Transform to ilr and apply bias ----
x2_ilr <- ilr(x2_prop) %>% as.matrix()
y2_ilr <- sweep(x2_ilr, 2, 35 * bias_ilr, FUN = "+")
y2_prop <- ilrInv(y2_ilr)

# ---- Reference before and after PCR ----
x1_prop <- ilrInv(x1_ilr)
y1_prop <- ilrInv(x1_ilr + 35 * bias_ilr)

# ---- Add pseudocount and normalize ----
x2_prop <- x2_prop + 1e-8; x2_prop <- x2_prop / rowSums(x2_prop)
y2_prop <- y2_prop + 1e-8; y2_prop <- y2_prop / rowSums(y2_prop)
x1_prop <- x1_prop + 1e-8; x1_prop <- x1_prop / sum(x1_prop)
y1_prop <- y1_prop + 1e-8; y1_prop <- y1_prop / sum(y1_prop)

# ---- Bray-Curtis distances ----
bray1 <- vapply(1:nrow(x2_prop), function(i) bray_curtis(x1_prop, x2_prop[i, ]), numeric(1))
bray2 <- vapply(1:nrow(y2_prop), function(i) bray_curtis(y1_prop, y2_prop[i, ]), numeric(1))
delta_bray <- bray1 - bray2

# ---- Data for plotting ----
df_all <- data.frame(
  v1 = x2_prop[, 1],
  v2 = x2_prop[, 2],
  v3 = x2_prop[, 3],
  val = delta_bray
)

x1.plot <- as.data.frame(t(x1_prop))
colnames(x1.plot) <- c("v1", "v2", "v3")

bias.plot <- as.data.frame(t(ilrInv(bias_ilr)))
colnames(bias.plot) <- c("v1", "v2", "v3")

zero.plot <- as.data.frame(t(ilrInv(c(0, 0))))
colnames(zero.plot) <- c("v1", "v2", "v3")

# ---- Plot ----
cbPalette <- c(
  "#E69F00", "#56B4E9", "#F0E442", "#009E73",
  "#D55E00", "#CC79A7"
)

cbPalette <- c(
  "#2166AC", "white","#B2182B","#B2182B"
)
# Use Okabe-Ito palette or other colorblind safe palette
ggtern(df_all, aes(v1, v2, v3, color = val)) +
  geom_point() +
  geom_point(data = bias.plot, aes(v1, v2, v3),
             shape = 21, fill = "#D55E00", color = "black", size = 5.5, stroke = 1.4) +  # orange
  geom_point(data = zero.plot, aes(v1, v2, v3),
             shape = 21, fill = "#009E73", color = "black", size = 4, stroke = 1.4) +    # green
  geom_point(data = x1.plot, shape = 21, fill = "#999", color = "black", size = 4.8, stroke = 1.3) +
  scale_color_gradientn(colors = cbPalette, name = "Bray-Curtis bias") +
  labs(T = "TAXA 1", L = "TAXA 2", R = "TAXA 3") +
  scale_size_continuous(range = c(1, 6), name = "Bray-Curtis bias") +
  theme_bw(base_size = 15) +
  theme(
    plot.margin = unit(c(-1, -1, -1, -1), "cm"),
    tern.axis.title.T = element_text(margin = margin(b = 50), size = 18),
    tern.axis.title.L = element_text(margin = margin(t = 50),size = 18),
    tern.axis.title.R = element_text(margin = margin(t = 50),size = 18),
    tern.axis.text = element_text(size = 15),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) 

