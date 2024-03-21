pacman::p_load(
  Rvcg, geomorph, MASS, Morpho, colorspace, cluster, DiscriMiner,
  ellipse, ggplot2, mclust, NbClust, shapes, vegan, tidyverse, processx,
  plotly
)

source("src/colLab.r")
source("src/read_multi_nts.R")

label_dc <- list(
  expression("T. occidentalis"), expression("T. triangulum")
)

# Adicionando em uma lista os arquivos de acordo com a nomenclatura

## Todas as populações

list_all_tabanus <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "nts")

## Tabanus occidentalis RS

list_ato <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "ato")

## Tabanus triangulum

list_tt <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "tt")

## Tabanus occidentalis TO

# list_tto <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "tto")

# Lendo as listas

setwd(paste0(getwd(), "/dados_analise/Tabanus"))

## Todas as populações

t_all <- readmulti_nts(list_all_tabanus)

### verificando como estãos agrupados meu array

dim(t_all)

## Tabanus triangulum

tt <- readmulti_nts(list_tt)

## Tabanus occidentalis RS

ato <- readmulti_nts(list_ato)

## Tabanus occidentalis TO

# tto <- readmulti_nts(list_tt)

## voltando para a pasta principal

setwd("../..")

# criando o diretório output

dir.create("output")

# Iniciando as análises Interespecíficas

## GPA

groups_inter <- list(t_all = t_all, tt = tt, ato = ato)

alg_inter <- lapply(groups_inter, gpagen)

alg_inter

## separando os alinhamentos de cada pop

alg_t_all <- alg_inter$t_all$coords

alg_tt <- alg_inter$tt$coords

alg_ato <- alg_inter$ato$coords

# alg_tto <- alg_inter$tto$coords

# Verificando a distribução de todo o dataset

setwd(paste0(getwd(), "/output"))

df_to_plot <- tibble(Csize = alg_inter$t_all$Csize)

## Separando as populações

tab_pop <- as.factor(c(rep("ato", 30), rep("tt", 30)))

## Teste do tamanho

size <- alg_inter$t_all$Csize

sink("output/shapiro_test.txt")
shapiro.test(size)
sink()

## teste t dos tamanhos das pop

teste_t <- t.test(size / 100 ~ tab_pop)
sink("output/t_test.txt")
teste_t
sink()

## Plot
df_to_plot <- tibble(size = size / 100, tab_pop)

ggplot(df_to_plot, aes(y = size, x = tab_pop)) +
  geom_boxplot(
    outlier.colour = NULL,
    aes_string(
      colour = "tab_pop",
      fill = "tab_pop"
    )
  ) +
  stat_summary(
    geom = "crossbar", width = 0.65, fatten = 0,
    color = "white",
    fun.data = function(x) {
      return(c(
        y = median(x),
        ymin = median(x), ymax = median(x)
      ))
    }
  ) +
  scale_color_manual(values = c("#14b8b8", "#BE3455")) +
  scale_fill_manual(values = c("#14b8b8", "#BE3455")) +
  labs(x = "Populations", y = "Centroid Size") +
  theme_bw() +
  theme(legend.position = "none") +
  annotate("text",
    x = 2, y = 10,
    label = paste("P=", round(teste_t$p.value, 4))
  ) +
  scale_x_discrete(labels = c(expression(italic("Tabanus occidentalis")), expression(italic("Tabanus triangulum"))))

ggsave("centroid_size",
  path = "output", width = 8, height = 8,
  device = "tiff", dpi = 700, limitsize = FALSE
)

## Detectando e excluindo outliers

# plotOutliers(alg_t_all)
#
# outliers <- find.outliers(alg_t_all)
#
# outliers$dist.sort

### Correlação do espaço de forma (shape space) com o espaço tangente (tangent space)

regdist(alg_t_all)

# plotRefToTarget(ref, alg_t_all[, , 58], links = links.trach, method = "vector", mag = 5, gridPars = GP1)

# tiff("Shape_species.tiff", units = "in", width = 20, height = 8, res = 300)

# layout(matrix(c(1, 2), 2, 2, byrow = TRUE))

plotRefToTarget(ref, alg_ato[, , 30],
  links = links.trach, method = "vector", mag = 5, gridPars = GP1,
  label = T, axes = T, useRefPts = T
)
title(main = "Tabanus occidentalis")

# plotRefToTarget(ref, alg_tto[, , 28],
#   links = links.trach, method = "vector", mag = 5, gridPars = GP1,
#   label = T, axes = T, useRefPts = T
# )
# title(main = "Tabanus occidentalis of Tocantins")

dev.off()

# PCA

pca_tt <- gm.prcomp(alg_t_all)

pca_tt$x

summary_pca_tt <- summary(pca_tt)

sink("output/PCA.txt")

summary_pca_tt$PC.summary

sink()

a <- cbind(pca_tt$x, tab_pop)

a <- data.frame(a)

a$tab_pop <- tab_pop

df_ell <- data.frame()

for (g in levels(a$tab_pop)) {
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(
    a[a$tab_pop == g, ],
    ellipse::ellipse(cor(Comp1, Comp2),
      scale = c(sd(Comp1), sd(Comp2)),
      centre = c(mean(Comp1), mean(a$Comp2))
    )
  )), tab_pop = g))
}

pca_plot <- ggplot(data = a, aes(x = Comp1, y = Comp2, colour = tab_pop)) +
  geom_point(size = 4, shape = 18) +
  geom_path(data = df_ell, aes(x = x, y = y, colour = tab_pop), linewidth = 0.3, linetype = 3) +
  # scale_color_manual(name = "Populations", labels = label_dc, values = c("#14b8b8", "#BE3455")) +
  # scale_shape_manual(name = "Populations", labels = label_dc, values = c(15, 16)) +
  scale_color_manual(
    name = "Populations",
    # labels = expression(italic(label_dc),
    labels = c(
      expression(italic("T. occidentalis")),
      expression(italic("T. triangulum"))
    ),
    values = c("#14b8b8", "#BE3455")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12), axis.title = element_text(size = 12),
    legend.text = element_text(size = 12), legend.title = element_text(size = 12)
  )

pca_plot

ggsave("pca.tiff",
  path = "output", units = "mm", width = 160, height = 82,
  device = "tiff", dpi = 700, limitsize = FALSE
)

## Plot 2 - PCA

especimes <- rownames(a)
rownames(a) <- NULL
data <- cbind(especimes, a)
data <- data %>%
  select(-tab_pop) %>%
  pivot_longer(!especimes, names_to = "PC", values_to = "value")

lista <- double(26)
pc_names <- list()
for (i in 1:26) {
  lista[[i]] <- summary[["PC.summary"]][[3, i]]
  pc_names[[i]] <- paste0("PC_", i)
}

pc_df <- cbind(pc_names, lista)
pc_df <- as.data.frame(pc_df, row) %>%
  mutate(id = row_number()) %>%
  rename(value = lista)

label_data_2 <- pc_df

# calculate the ANGLE of the labels
number_of_bar <- nrow(pc_df)
angle <- 90 - 360 * (label_data_2$id - 0.5) / number_of_bar # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90

label_data_2$hjust <- ifelse(angle <- 90, 1, 0)

# flip angle BY to make them readable
label_data_2$angle <- ifelse(angle <- 90, angle + 180, angle)
# ----- ------------------------------------------- ---- #

# Start the plot

p <- ggplot(pc_df, aes(x = as.factor(id), y = as.numeric(value))) + # Note that id is a factor. If x is numeric, there is some space between the first bar

  # This add the bars with a blue color
  geom_bar(stat = "identity", fill = alpha("skyblue", 0.7)) +

  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(0, 1) +

  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm"),
    # Adjust the margin to make in sort labels are not truncated!
  ) +

  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +

  # Add the labels, using the label_data dataframe that we have created before
  geom_text(
    data = label_data_2,
    aes(x = id, y = as.numeric(value), label = pc_names, hjust = hjust),
    color = "black", fontface = "bold", alpha = 0.6,
    size = 2.5,
    angle = 90,
    inherit.aes = FALSE
  )

p

# Anova

gdf <- geomorph.data.frame(
  coords = alg_t_all,
  size = df_to_plot$size,
  species = tab_pop
)


fit <- procD.lm(coords ~ size * species,
  data = gdf, iter = 999, turbo = TRUE,
  RRPP = FALSE, print.progress = FALSE
) # randomize raw values
fit2 <- procD.lm(coords ~ size * species,
  data = gdf, iter = 999, turbo = TRUE,
  RRPP = TRUE, print.progress = FALSE
) # randomize residuals



P <- plotAllometry(fit,
  size = gdf$size, logsz = TRUE,
  pch = 19, col = as.numeric(interaction(gdf$species))
)
make_ggplot(P)

sink("output/anova_fit1.txt")
summary(fit)
sink()

sink("output/anova_fit2.txt")
summary(fit2)
sink()
# Comparação par a par

gp <- interaction(tab_pop)

PW <- pairwise(fit2, groups = gp, covariate = NULL)

sink("output/parwise_anova.txt")
summary(PW)
sink()

## MANOVA Wilks's lambda ver

manova_w <- manova(pca_tt$x[, 1:15] ~ tab_pop * df_to_plot$size)

sink("output/manova.txt")
summary(manova_w, test = c("Wilks"))
sink()

## CVA forma matriz de confusão

cva_1 <- CVA(alg_t_all, tab_pop,
  weighting = T,
  tolinv = 1e-10, plot = T,
  rounds = 1000, cv = T,
  p.adjust.method = "bonferroni"
)

sink("output/cva.txt")
cva_1$Dist
sink()

CVscores <- data.frame(cva_1$CVscores)

CVscores$tab_pop <- tab_pop


## Visualização da forma - forma média

# ref <- mshape(alg_t_all)

links.trach <- matrix(c(
  1, 2, 1, 3, 1, 15, 2, 4, 2, 6, 4, 5, 5, 6, 5, 7, 6, 7, 7, 8,
  8, 9, 9, 10, 9, 11, 11, 12, 12, 13, 13, 14, 14, 15
), nrow = 34, ncol = 2, byrow = T)


GP1 <- gridPar(
  pt.bg = "#14b8b8", link.col = "#14b8b8", link.lty = 1, pt.size = 1,
  tar.pt.bg = "#BE3455", tar.pt.size = 1, tar.link.col = "#BE3455",
  tar.link.lty = 1
)

plotRefToTarget(cva_1[["groupmeans"]][, , "ato"], cva_1[["groupmeans"]][, , "tt"],
  links = links.trach, method = "points", mag = 5, gridPars = GP1,
  label = T, axes = T, useRefPts = T
)

tiff("wings_cva.tiff", units = "in", width = 5, height = 5, res = 700)

plotRefToTarget(cva_1[["groupmeans"]][, , "ato"], cva_1[["groupmeans"]][, , "tt"],
  links = links.trach, method = "points", mag = 5, gridPars = GP1,
  label = T, axes = T, useRefPts = T
)

dev.off()

## diferentes plots teste

ggplot(CVscores, aes(CV.1, fill = tab_pop)) +
  geom_histogram() +
  theme_bw() +
  scale_fill_manual(
    name = "Populations",
    # labels = expression(italic(label_dc),
    labels = c(
      expression(italic("T. occidentalis")),
      expression(italic("T. triangulum"))
    ),
    values = c("#14b8b8", "#BE3455")
  ) +
  labs(x = "Canonical Variate 1", y = "Frequency")


ggsave("cva.tiff",
  path = "output", units = "mm", width = 160, height = 82,
  device = "tiff", dpi = 700, limitsize = FALSE
)

# ## plot ggplot
#
df_ell <- data.frame()

for (g in levels(CVscores$tab_pop)) {
  df_ell <- rbind(
    df_ell,
    cbind(as.data.frame(with(
      CVscores[CVscores$tab_pop == g, ],
      ellipse(cor(CV.1, CV.2),
        scale = c(sd(CV.1), sd(CV.2)),
        centre = c(mean(CV.1), mean(CV.2))
      )
    )), tab_pop = g)
  )
}

# cva_plot <- ggplot(data = CVscores, aes(x = CV.1, y = CV.2, colour = tab_pop)) +
#   geom_point(size = 4, shape = 18) +
#   geom_path(data = df_ell, aes(x = x, y = y, colour = tab_pop), size = 0.3, linetype = 3) +
#   scale_color_manual(name = "Populations", labels = label_dc, values = c("#1b9e77", "darkorange4", "blue")) +
#   scale_shape_manual(name = "Populations", labels = label_dc, values = c(15, 16, 17, 18)) +
#   theme_classic() +
#   theme(
#     axis.text = element_text(size = 12), axis.title = element_text(size = 12),
#     legend.text = element_text(size = 12), legend.title = element_text(size = 12)
#   )
#
# cva_plot
