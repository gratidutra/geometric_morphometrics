pacman::p_load(
  Rvcg, geomorph, MASS, Morpho, geomorph, colorspace, cluster, DiscriMiner,
  ellipse, ggplot2, mclust, NbClust, shapes, vegan, tidyverse, processx
)

source("colLab.r")
source("read_multi_nts.R")

# Adicionando em uma lista os arquivos de acordo com a nomenclatura

## Todas as populações

list_all_tabanus <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "nts")

## Tabanus occidentalis RS

list_ato <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "ato")

## Tabanus triangulum

list_tt <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "tt")

## Tabanus occidentalis TO

list_tto <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "tto")

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

tto <- readmulti_nts(list_tt)

## voltando para a pasta principal

setwd("../..")

# criando o diretório output

dir.create("output")

# Iniciando as análises Interespecíficas

## GPA

groups_inter <- list(t_all = t_all, tt = tt, ato = ato, tto = tto)

alg_inter <- lapply(groups_inter, gpagen)

alg_inter

## separando os alinhamentos de cada pop

alg_t_all <- alg_inter$t_all$coords

alg_tt <- alg_inter$tt$coords

alg_ato <- alg_inter$ato$coords

alg_tto <- alg_inter$tto$coords

# Verificando a distribução de todo o dataset

setwd(paste0(getwd(), "/output"))

df_to_plot <- tibble(Csize = alg_inter$t_all$Csize)

tiff("Csize_histogram.tiff", units = "in", width = 5, height = 5, res = 300)

ggplot(df_to_plot, aes(Csize)) +
  geom_histogram() +
  theme_light() +
  ylab("Frequency")

dev.off()

## Separando as populações

tab_pop <- as.factor(c(rep("ato", 30), rep("tt", 30), rep("tto", 28)))

## Teste do tamanho

size <- alg_inter$t_all$Csize

shapiro.test(size)

df_to_plot <- tibble(size, tab_pop)

ggplotly(ggplot(df_to_plot, aes(y = size, x = tab_pop, fill = tab_pop)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic())

## anova dos tamanhos das pop

anova <- aov(size ~ tab_pop)

summary(anova)

TukeyHSD(anova)

## Detectando e excluindo outliers

# plotOutliers(alg_t_all)
#
# outliers <- find.outliers(alg_t_all)
#
# outliers$dist.sort

### Correlação do espaço de forma (shape space) com o espaço tangente (tangent space)

regdist(alg_t_all)

### Visualização da forma - forma média

ref <- mshape(alg_t_all)

links.trach <- matrix(c(
  1, 2, 1, 3, 1, 15, 2, 4, 2, 6, 4, 5, 5, 6, 5, 7, 6, 7, 7, 8,
  8, 9, 9, 10, 9, 11, 11, 12, 12, 13, 13, 14, 14, 15
), nrow = 34, ncol = 2, byrow = T)

GP1 <- gridPar(pt.bg = "#304D63", link.col = "#8FB9AA", link.lty = 10, tar.pt.size = 3)

# plotRefToTarget(ref, alg_t_all[, , 58], links = links.trach, method = "vector", mag = 5, gridPars = GP1)

tiff("Shape_species.tiff", units = "in", width = 20, height = 8, res = 300)

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))

plotRefToTarget(ref, alg_tt[, , 30],
  links = links.trach, method = "vector", mag = 5, gridPars = GP1,
  label = T, axes = T, useRefPts = T
)
title(main = "Tabanus triangulum")

plotRefToTarget(ref, alg_ato[, , 30],
  links = links.trach, method = "vector", mag = 5, gridPars = GP1,
  label = T, axes = T, useRefPts = T
)
title(main = "Tabanus occidentalis of Rio Grande do Sul")

plotRefToTarget(ref, alg_tto[, , 28],
  links = links.trach, method = "vector", mag = 5, gridPars = GP1,
  label = T, axes = T, useRefPts = T
)
title(main = "Tabanus occidentalis of Tocantins")

dev.off()

# PCA

pca_tt <- gm.prcomp(alg_t_all)

pca_tt$x

sink("PCA.txt")

summary(pca_tt)

sink()

label_dc <- list(
  expression("T. occidentalis Rio Grande do Sul"), expression("T. triangulum"),
  expression("T. occidentalis Tocantins")
)

a <- cbind(pca_tt$x, tab_pop)

a <- data.frame(a)

a$tab_pop <- tab_pop

df_ell <- data.frame()

for (g in levels(a$tab_pop)) {
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$tab_pop == g, ], ellipse(cor(Comp1, Comp2),
    scale = c(sd(Comp1), sd(Comp2)),
    centre = c(mean(Comp1), mean(Comp2))
  ))), tab_pop = g))
}

pca_plot <- ggplot(data = a, aes(x = Comp1, y = Comp2, colour = tab_pop)) +
  geom_point(size = 4, shape = 18) +
  geom_path(data = df_ell, aes(x = x, y = y, colour = tab_pop), size = 0.3, linetype = 3) +
  scale_color_manual(name = "Populations", labels = label_dc, values = c("#1b9e77", "darkorange4", "blue")) +
  scale_shape_manual(name = "Populations", labels = label_dc, values = c(15, 16, 17, 18)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12), axis.title = element_text(size = 12),
    legend.text = element_text(size = 12), legend.title = element_text(size = 12)
  )

pca_plot

# Anova

anova <- procD.lm(alg_t_all ~ tab_pop, iter = 9999, RRPP = TRUE)

anova

resultado <- summary(anova)

resultado

# Comparação par a par

gp <- interaction(tab_pop)

PW <- pairwise(anova, groups = gp, covariate = NULL)

summary(PW)

## MANOVA Wilks's lambda ver

manova_w <- manova(pca_tt$x[, 1:26] ~ tab_pop)

summary(manova_w, test = c("Hotelling-Lawley"))


## CVA forma matriz de confusão

cva.1 <- CVA(alg_t_all, tab_pop,
  weighting = T,
  tolinv = 1e-10, plot = T,
  rounds = 1000, cv = T,
  p.adjust.method = "bonferroni"
)

CVscores <- data.frame(cva.1$CVscores)

CVscores$tab_pop <- tab_pop

## diferentes plots teste

ggplot(CVscores, aes(CV.1, fill = tab_pop)) +
  geom_histogram() +
  theme_classic()


## plot ggplot

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

cva_plot <- ggplot(data = CVscores, aes(x = CV.1, y = CV.2, colour = tab_pop)) +
  geom_point(size = 4, shape = 18) +
  geom_path(data = df_ell, aes(x = x, y = y, colour = tab_pop), size = 0.3, linetype = 3) +
  scale_color_manual(name = "Populations", labels = label_dc, values = c("#1b9e77", "darkorange4", "blue")) +
  scale_shape_manual(name = "Populations", labels = label_dc, values = c(15, 16, 17, 18)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12), axis.title = element_text(size = 12),
    legend.text = element_text(size = 12), legend.title = element_text(size = 12)
  )

cva_plot
