library(combinatorixData)
library(ggplot2)
library(reshape2)
library(lme4)
library(lmerTest)
library(combinatorix)

#setwd("~/rhizogenomics/experiments/2017/today/")

################ ANALYZE ################
# Size in pixels
data(binvalidation.sizes)
Sizes <- binvalidation.sizes
binvalidation.sizes <- NULL

res <- analyze_measurement(variable = "npixels",
                           dat = Sizes,
                           trans = "log2plus1",
                           clean = TRUE,
                           contam = FALSE,
                           prefix = "npixels_log_nocontam_nodead_")
AIC(res$m1,res$m2)
BIC(res$m1,res$m2)

# Size in pixels normalized to sticker
res <- analyze_measurement(variable = "Normalized.size",
                           dat = Sizes,
                           trans = "log2plus1",
                           clean = TRUE,
                           contam = FALSE,
                           prefix = "normalized.size_log_nocontam_nodead_")
AIC(res$m1,res$m2)
BIC(res$m1,res$m2)

# Convex hull area in pixels
res <- analyze_measurement(variable = "Hull.area",
                           dat = Sizes,
                           trans = "log2plus1",
                           clean = TRUE,
                           contam = FALSE,
                           prefix = "hull.area_log_nocontam_nodead_")
AIC(res$m1,res$m2)
BIC(res$m1,res$m2)

# Convex hull area in pixels normalized to sticker
res <- analyze_measurement(variable = "Normalized.hull.area",
                           dat = Sizes,
                           trans = "log2plus1",
                           clean = TRUE,
                           contam = FALSE,
                           prefix = "normalized.hull.area_log_nocontam_nodead_")
AIC(res$m1,res$m2)
BIC(res$m1,res$m2)

# Process validation coefs
validation_coefs <- summary(res$m2)$coef
validation_coefs <- data.frame(validation_coefs)
validation_coefs <- validation_coefs[ grep("^Treatment", row.names(validation_coefs)), ]
colnames(validation_coefs) <- c("Validation.Effect", "Validation.SE",
                                "Validation.df","Validation.t","Validation.p")
validation_coefs$Strain <- row.names(validation_coefs)
validation_coefs$Strain <- gsub(pattern = "^Treatment", replacement = "", x = validation_coefs$Strain)
validation_coefs


#######
data(comb.phenotypes)
comb.phenotypes <- subset(comb.phenotypes,Phenotype == "hull.area")
row.names(comb.phenotypes) <- sub(pattern = "^X",replacement = "",x = as.character(comb.phenotypes$Strain))
head(comb.phenotypes)

validation_coefs$Comb.Estimate <- comb.phenotypes[ validation_coefs$Strain, "Coefficient"]
head(validation_coefs)

p1 <- ggplot(validation_coefs, aes(x = Comb.Estimate, y = Validation.Effect)) +
  geom_point(aes(color = Validation.p < 0.01), size = 3) +
  geom_text(aes(label = Strain), vjust = 1.5) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold", color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size= 2))
p1
ggsave("comparison_hull_area.svg", p1, width = 6, height = 4)

subset(comb.phenotypes, p.value < 0.05)

