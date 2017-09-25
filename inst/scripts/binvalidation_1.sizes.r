library(combinatorixData)
library(ggplot2)
library(reshape2)
library(lme4)
library(lmerTest)
############################ FUNCTIONS ##################################
#' Analyze one measurment in binary validation
#' 
#' @param variable String with the name of the variable to study. Must be a column
#' name in dat
#' @param dat data.frame with the measurments and metadata to analyze. Most contain columns:
#' PlateID, Col, Row, PicRep, Bacteria, Strain and Contaminated plus a column corresponding 
#' to variable. Other columns are ignored
#' @param trans String either "log2", "log2plus1" or "none".
#' @param clean logical indicating wheter observations of value zero (dead plants) should be removed
#' @param contam logical indicating whether to keep wells marked as contaminated
#' @param prefix string indicating prefix for created files
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
analyze_measurement <- function(variable, dat, trans = "log2",
                                clean = TRUE, contam = FALSE,
                                prefix = ""){
  # variable <- "Normalized.size"
  # dat <- Sizes
  # log.transform <- TRUE
  # clean <- TRUE
  # contam <- FALSE
  
  # Pre-process
  if(!contam)
    dat <- subset(dat, !Contaminated)
  if(clean)
    dat <- dat[ dat[,variable] > 0, ]
  
  if(trans == "log2plus1"){
    dat[,variable] <- log2(dat[,variable] + 1)
  }else if(trans == "log2" ){
    dat[,variable] <- log2(dat[,variable])
  }else if(trans == "none"){
    # do nothing
  }else{
    stop("ERROR")
  }
    
  ## Check consistency between reps
  dat2 <- subset(dat,PlateID %in% unique(dat$PlateID[ dat$PicRep == "B" ]))
  f1 <- paste("PlateID + Row + Col ~ PicRep",sep = "")
  f1 <- formula(f1)
  f1
  dat2 <- acast(f1,data = dat2,value.var = variable)
  dat2 <- as.data.frame(dat2)
  p1 <- ggplot(dat2,aes(x = A, y = B)) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_point() +
    geom_smooth(method = "lm") +
    ggtitle(label = variable) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA))
  #p1
  filename <- paste(prefix,"picture_consistency.svg", sep = "")
  ggsave(filename,p1,width = 5, height = 5)
  
  ## Aggregate and plot
  f1 <- paste(variable, " ~ PlateID + Col + Row + Bacteria + Strain", sep = "")
  f1 <- formula(f1)
  f1
  dat <- aggregate(f1, data = dat,FUN = mean)
  p1 <- ggplot(dat, aes_string(x = "Strain", y = variable)) +
    geom_boxplot(aes(color = Bacteria)) +
    geom_point(aes(color = Bacteria),position = position_jitterdodge()) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text.x = element_text(color = "black", angle = 90))
  #p1
  filename <- paste(prefix,"per_strain.svg", sep = "")
  ggsave(filename,p1,width = 10, height = 5)
  
  # Models
  f1 <- paste(variable," ~ Treatment + PlateID", sep = "")
  f1 <- formula(f1)
  # f1
  dat$Treatment <- as.character(dat$Strain)
  dat$Treatment[ dat$Bacteria == "No Bacteria" ] <- "No Bacteria"
  dat$Treatment <- relevel(factor(dat$Treatment), ref = "No Bacteria")
  m1 <- lm(f1, data = dat)
  print(summary(m1))
  
  f1 <- paste(variable," ~ Treatment + (1|PlateID)", sep = "")
  f1 <- formula(f1)
  m2 <- lmer(f1, data = dat)
  print(summary(m2))
  # AIC(m1,m2)
  # BIC(m1,m2)
  
  # Plate effects
  f1 <- paste(variable," ~ PlateID + Bacteria", sep = "")
  f1 <- formula(f1)
  #f1
  dat2 <- as.data.frame(acast(PlateID ~ Bacteria,
                              data = aggregate(formula = f1,
                                               data = dat, FUN = mean),
                              fun.aggregate = mean,value.var = variable))
  colnames(dat2) <- c("Bacteria","No.Bacteria")
  p1 <- ggplot(dat2, aes_string(y = "Bacteria", x = 'No.Bacteria')) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_point() +
    geom_smooth(method = "lm") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA))
  #p1
  filename <- paste(prefix,"plate_effects.svg", sep = "")
  ggsave(filename,p1,width = 7, height = 5)
  
  res <- list(m1 = m1, m2 = m2)
  return(res)
}
##################################################################

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

