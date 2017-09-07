library(ggplot2)
library(reshape2)
library(lme4)
library(lmerTest)

Weights <- read.table("FWdata.txt", sep = "\t", header = TRUE)
head(Weights)
Weights$PlateID <- paste(Weights$Condition, Weights$Replicate, sep = ".")
head(Weights)


dat <- Weights
prefix <- "fresh.weight"
variable <- "FreshW"

dat <- subset(dat,FreshW > 0)
dat$FreshW <- log2(dat$FreshW + 1)

p1 <- ggplot(dat, aes(x = Condition, y = FreshW)) +
  geom_boxplot(aes(color = Bacteria)) +
  geom_point(aes(color = Bacteria),position = position_jitterdodge()) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90))
p1
filename <- paste(prefix,"per_strain.svg", sep = "")
ggsave(filename,p1,width = 7, height = 5)

dat$Treatment <- as.character(dat$Condition)
dat$Treatment[ dat$Bacteria == "NB"] <- "NB"
dat$Treatment <- relevel(factor(dat$Treatment),ref = "NB")

f1 <- paste(variable," ~ Treatment + PlateID", sep = "")
f1 <- formula(f1)
f1
m1 <- lm(f1, data = dat)
print(summary(m1))

f1 <- paste(variable," ~ Treatment + (1|PlateID)", sep = "")
f1 <- formula(f1)
m2 <- lmer(f1, data = dat)
print(summary(m2))

AIC(m1,m2)
BIC(m1,m2)


