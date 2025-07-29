library(tidyverse)

library(ggtern)
library(driver) # jsilve24/driver
aitchison_norm <- function(comp) {
  sqrt(rowSums(comp^2))
}


bias <- c(1,0.85,1)
bias <- as.numeric(miniclo(bias))
bias <- ilr(bias)


n_grid <- 70  # larger n = finer grid
x.prop <- uniform_simplex_grid(n_grid)
epsilon <- 1e-8
x.prop <- x.prop + epsilon
x.prop <- x.prop / rowSums(x.prop)

x<- ilr(x.prop)
x <- as(x, "matrix")
y <- x
for (i in 1:nrow(x)) {
  y[i,] <- x[i,]+35*bias
}
y.prop <- driver::ilrInv(y)
aitchison.x <- aitchison_norm(x)
aitchison.y <- aitchison_norm(y)
delta_aitchison <- aitchison.y - aitchison.x

# Create data frame for plotting
df_all <- data.frame(
  v1 = x.prop[, 1],
  v2 = x.prop[, 2],
  v3 = x.prop[, 3],
  val = delta_aitchison
)

# Identify points near zero change
df_near_zero <- df_all %>% filter(abs(val) < 0.05)

# Reference points
bias.plot <- t(data.frame(as.numeric(ilrInv(bias))))
colnames(bias.plot) <- c("v1", "v2", "v3")
zero.plot <- data.frame(t(as.numeric(ilrInv(c(0, 0)))))
colnames(zero.plot) <- c("v1", "v2", "v3")

cbPalette <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # green
  "#D55E00", 
  "#F0E442",  # yellow
  # blue
  # reddish orange
  "#CC79A7"   # pink
)

cbPalette <- c(
  "#2166AC", "white","white","#B2182B"
)


ggtern(data = df_all, aes(v1, v2, v3, color = val)) +
  geom_point() +
  geom_point(data = bias.plot, aes(v1, v2, v3),
             shape = 21, fill = "#D55E00", color = "black", size = 5.5, stroke = 1.4) +  # orange
  geom_point(data = zero.plot, aes(v1, v2, v3),
             shape = 21, fill = "#009E73", color = "black", size = 4, stroke = 1.4) +    # green
  scale_color_gradientn(
    colors = cbPalette,  # should be ordered like: low -> mid -> high (e.g., blue, white, red)
    values = scales::rescale(c(min(df_all$val), 0, max(df_all$val))),  # force zero to midpoint
    name = "Aitchison bias",
    limits = c(min(df_all$val), max(df_all$val))
  )+
  labs(
    T = "TAXA 1",
    L = "TAXA 2",
    R = "TAXA 3"
  ) +
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

