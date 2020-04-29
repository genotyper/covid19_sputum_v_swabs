#####
#data step (two data.frames made: 1) dat = all data, 2) pos = only those deemed covid positive)
#####
dat = read.csv("data_04062020.csv", stringsAsFactors = FALSE)
dat = dat[dat$Sample.Type == "Sputum" | dat$Sample.Type == "Swab",]

met = read.csv("met.csv", stringsAsFactors = FALSE)
## map samples between metadata and qpcr data using accession numbers, drop rows not present in both.
met = met[met$accession %in% dat$accession,]
dat = dat[dat$accession %in% met$accession,]

## rRT-PCR results of SARS-CoV2 N1-3 regions scored as Ct values
# negative detections coded as zeros, score triple non-zero Ct values as positive
hold = dat$N1.Ct > 0 & dat$N2.Ct > 0 & dat$N3.Ct > 0
dat$Detection = hold
dat$Detection = ifelse(dat$Detection == "TRUE", 1, 0)

## merge dat and met using accession for mapping
dat$accession = as.character(dat$accession)
met$accession = as.character(met$accession)
dat = merge(dat, met, by = "accession")

## cleanup gender
dat$gender = toupper(dat$gender)
dat = dat[dat$gender != "T",]
dat = dat[dat$gender != "",]

## create vectors to store viral copy counts
dat$N1.Abundance = NA
dat$N2.Abundance = NA
dat$N3.Abundance = NA

## calculate copy counts for swabs
swab = dat[dat$Sample.Type == "Swab",]
look = which(swab$Detection == 1)
swab$N1.Abundance[look] = exp((swab$N1.Ct[look]-41.75)/-1.447)
swab$N2.Abundance[look] = exp((swab$N2.Ct[look]-42.646)/-1.482)
swab$N3.Abundance[look] = exp((swab$N3.Ct[look]-43.732)/-1.546)
# repeat for sputum
spit = dat[dat$Sample.Type == "Sputum",]
look = which(spit$Detection == 1)
spit$N1.Abundance[look] = exp((spit$N1.Ct[look]-39.518)/-1.28)
spit$N2.Abundance[look] = exp((spit$N2.Ct[look]-40.427)/-1.238)
spit$N3.Abundance[look] = exp((spit$N3.Ct[look]-41.429)/-1.362)
# merge back together
dat = rbind(spit, swab)
rm(spit, swab)

## average qpcr results across three assays
dat$mean.Ct = (dat$N1.Ct + dat$N2.Ct + dat$N3.Ct)/3
dat$mean.Abundance = (dat$N1.Abundance + dat$N2.Abundance + dat$N3.Abundance)/3

## for only looking at positive detections
pos = dat[dat$Detection == 1, ]

#####
#correlation of Ct of the three assays for those that were positive for covid
#####
cor(pos[,2:4])
cor(pos[,13:15])

## plotting
library(ggplot2)

p = ggplot(pos, aes(N1.Ct, N2.Ct)) + geom_point(alpha = .4) + geom_smooth(method = "lm", color = "black") + 
  theme_bw(base_size = 15)  + scale_y_continuous(breaks=seq(0,40,5)) + scale_x_continuous(breaks=seq(0,40,5))
p1 = p
# png(file = "12Ct.png", width = 4, height = 4, units = "in", res = 600)
# print(p)
# dev.off()

p = ggplot(pos, aes(N1.Ct, N3.Ct)) + geom_point(alpha = .4) + geom_smooth(method = "lm", color = "black") + 
  theme_bw(base_size = 15)  + scale_y_continuous(breaks=seq(0,40,5)) + scale_x_continuous(breaks=seq(0,40,5))
p2 = p
# png(file = "13Ct.png", width = 4, height = 4, units = "in", res = 600)
# print(p)
# dev.off()

p = ggplot(pos, aes(N2.Ct, N3.Ct)) + geom_point(alpha = .4) + geom_smooth(method = "lm", color = "black") + 
  theme_bw(base_size = 15)  + scale_y_continuous(breaks=seq(0,40,5)) + scale_x_continuous(breaks=seq(0,40,5))
p3 = p
# png(file = "23Ct.png", width = 4, height = 4, units = "in", res = 600)
# print(p)
# dev.off()

## merge figure using patchwork
library(patchwork)
# https://patchwork.data-imaginist.com/articles/guides/layout.html
library(grid)
library(ggplotify)

png(file = "assay_correlations.png", width = 9, height = 3, units = "in", res = 600)
print(p1 + p2 + p3)
dev.off()

#####
#viral load as a function of sample type
#####

for(i in 2:4){
  print(names(pos)[i])
  print(shapiro.test(dat[,i]))
}
hist(pos[,2])
hist(pos[,3])
hist(pos[,4])

## calculated copy counts may be better value for analysis than raw Ct - however, different distribution
for(i in 13:15){
  print(names(pos)[i])
  print(shapiro.test(dat[,i]))
}
hist(pos[,13])
hist(pos[,14])
hist(pos[,15])
## to handle skewed distribution we can use two approaches -- log transform before analysis or use poisson (or possibly negative binomial) glm. 
hist(log(pos[,13]))
hist(log(pos[,14]))
hist(log(pos[,15]))
## log transformation helpful, use for plotting

#### formal viral load analysis
library(AER)
library(MASS)

## poisson model evaluated first
lm = glm(mean.Abundance ~ Sample.Type, data = pos, family = "poisson")
dispersiontest(lm) #over dispersed - grounds for using quasipoisson or negative binomial instead
## evaluate each model, and use the one which offers the lowest residual deviance
lm1 = glm(mean.Abundance ~ Sample.Type, data = pos, family = "quasipoisson")
lm2 = glm.nb(mean.Abundance ~ Sample.Type, data = pos)
summary(lm1)
summary(lm2)
anova(lm, lm1, lm2)
## negative binomial appears to be best fit

#take home
lm2 = glm.nb(round(mean.Abundance,0) ~ Sample.Type * gender, data = pos)
summary(lm2) ## no sig interaction; drop interaction term

lm2 = glm.nb(round(mean.Abundance,0) ~ Sample.Type + gender, data = pos)
summary(lm2) 
#                 Estimate    Std.Error z value   Pr(>|z|)    
# (Intercept)      15.2135     0.3639   41.809    < 2e-16 ***
# Sample.TypeSwab  -0.3115     0.3721   -0.837    0.402    
# genderM           1.3330     0.2778   4.798     1.6e-06 ***

est = cbind(Estimate = coef(lm2), confint(lm2))
# Incident ratio can be calculated by exponentiating the model coefficients
exp(est)
# 3.8% viral load increase in males from females

### plotting viral load distributions per sample type and sex
# dropped from paper

pos$sex = pos$gender
pos$sex = gsub("M", "Male", pos$sex)
pos$sex = gsub("F", "Female", pos$sex)
p = ggplot(pos, aes(mean.Abundance)) + geom_histogram(binwidth = 0.5) + facet_grid(sex ~ .) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() + labs(x = "Viral Load (copy count)") + ggtitle(label = "Sex Distribution Comparison") +
  theme(plot.title = element_text(hjust = 0.5, size = 8), axis.text = element_text(size = 7), 
        axis.title = element_text(size = 7, face = "bold"), strip.text = element_text(size = 8))
png(file = "hist_gender_v2.png", width = 3, height = 2, units = "in", res = 600)
print(p)
dev.off()
p1 = p
pos$type = pos$Sample.Type
pos$type = gsub("Sputum", "AST", pos$type)
pos$type = gsub("Swab", "NP Swab", pos$type)

p = ggplot(pos, aes(mean.Abundance)) + geom_histogram(binwidth = 0.5) + facet_grid(type ~ .) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() + labs(x = "Viral Load (copy count)") + ggtitle(label = "Sample Type Distribution Comparison") +
  theme(plot.title = element_text(hjust = 0.5, size = 8), axis.text = element_text(size = 7), 
        axis.title = element_text(size = 7, face = "bold"), strip.text = element_text(size = 8))
png(file = "hist_sampletype_v2.png", width = 3, height = 2, units = "in", res = 600)
print(p)
dev.off()
p2 = p
## merge figure using patchwork
library(patchwork)
# https://patchwork.data-imaginist.com/articles/guides/layout.html
library(grid)
library(ggplotify)

tiff(file = "hist_load_v2.tif", width = 6, height = 2, units = "in", res = 600, compression = 'lzw')
print(p2 + p1)
dev.off()


#####
# chi-sq; compare positive detection rate
#####

chisq.test(dat$Sample.Type, dat$Detection)
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  dat$Sample.Type and dat$Detection
# X-squared = 0.38512, df = 1, p-value = 0.5349

table(dat$Sample.Type, dat$Detection)
#           0    1
# Sputum  719   59
# Swab   3266  297

chisq.test(dat$gender, dat$Detection)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  dat$gender and dat$Detection
# X-squared = 15.668, df = 1, p-value = 7.547e-05

table(dat$gender, dat$Detection)
#     0    1
# F 2253  162
# M 1732  194

spit = dat[dat$Sample.Type == "Sputum",]
swab = dat[dat$Sample.Type == "Swab",]

chisq.test(spit$gender, spit$Detection)
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  spit$gender and spit$Detection
# X-squared = 0.20038, df = 1, p-value = 0.6544

chisq.test(swab$gender, swab$Detection)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  swab$gender and swab$Detection
# X-squared = 16.572, df = 1, p-value = 4.683e-05

mspit = spit[spit$gender == "M",]
sum(mspit$Detection)/nrow(mspit)

fspit = spit[spit$gender == "F",]
sum(fspit$Detection)/nrow(fspit)

male = dat[dat$gender == "M",]
chisq.test(male$Sample.Type, male$Detection)

fmale = dat[dat$gender == "F",]
chisq.test(fmale$Sample.Type, fmale$Detection)

summary(spit$mean.Abundance)
