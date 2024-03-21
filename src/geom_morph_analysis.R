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

# Creting lists for each landmarks

## All populations

list_all_tabanus <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "nts")

## Tabanus occidentalis RS

list_ato <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "ato")

## Tabanus triangulum

list_tt <- list.files(path = paste0(getwd(), "/dados_analise/Tabanus"), pattern = "tt")

# Reading lists

setwd(paste0(getwd(), "/dados_analise/Tabanus"))

## All populations

t_all <- readmulti_nts(list_all_tabanus)

### verified dimension array 

dim(t_all)

## Tabanus triangulum

tt <- readmulti_nts(list_tt)

## Tabanus occidentalis RS

ato <- readmulti_nts(list_ato)

setwd("../..")

# create path output

dir.create("output")

## GPA

groups_inter <- list(t_all = t_all, tt = tt, ato = ato)

alg_inter <- lapply(groups_inter, gpagen)

alg_inter

## separating the alignment each populations

alg_t_all <- alg_inter$t_all$coords

alg_tt <- alg_inter$tt$coords

alg_ato <- alg_inter$ato$coords

setwd(paste0(getwd(), "/output"))

df_to_plot <- tibble(Csize = alg_inter$t_all$Csize)

## separating species 

tab_pop <- as.factor(c(rep("ato", 30), rep("tt", 30)))

## Shapiro Test

size <- alg_inter$t_all$Csize

sink("output/shapiro_test.txt")
shapiro.test(size)
sink()

## T teste centroid size 
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

## Detecting e excluding outliers

# plotOutliers(alg_t_all)
#
# outliers <- find.outliers(alg_t_all)
#
# outliers$dist.sort

### Correlation shape space with tangent space

regdist(alg_t_all)

plotRefToTarget(ref, alg_ato[, , 30],
                links = links.trach, method = "vector", mag = 5, gridPars = GP1,
                label = T, axes = T, useRefPts = T
)
title(main = "Tabanus occidentalis")

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
  scale_color_manual(
    name = "Populations",
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

# Anova

gdf <- geomorph.data.frame(
  coords = alg_t_all,
  size = df_to_plot$size,
  species = tab_pop
)

fit2 <- procD.lm(coords ~ size * species,
                 data = gdf, iter = 999, turbo = TRUE,
                 RRPP = TRUE, print.progress = FALSE
)

P <- plotAllometry(fit,
                   size = gdf$size, logsz = TRUE,
                   pch = 19, col = as.numeric(interaction(gdf$species))
)

#make_ggplot(P)

sink("output/anova_fit1.txt")
summary(fit)
sink()

sink("output/anova_fit2.txt")
summary(fit2)
sink()

# Comparision peer

gp <- interaction(tab_pop)

PW <- pairwise(fit2, groups = gp, covariate = NULL)

sink("output/parwise_anova.txt")
summary(PW)
sink()

## MANOVA Wilks's lambda

manova_w <- manova(pca_tt$x[, 1:15] ~ tab_pop * df_to_plot$size)

sink("output/manova.txt")
summary(manova_w, test = c("Wilks"))
sink()

## CVA cross-validation

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

## Visualization of shape - mean

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

## Tests plots 

ggplot(CVscores, aes(CV.1, fill = tab_pop)) +
  geom_histogram() +
  theme_bw() +
  scale_fill_manual(
    name = "Populations",
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

#  plot ggplot
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
