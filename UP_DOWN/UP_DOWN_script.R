## 8/5/2020

## Get UP/Down regulated genes regarding the mean for termites.
rm(list = ls())
setwd("termites")
library(data.table)

######

## Female_Alates 
## Open the counts
Female_Alates <- read.table("Female_Alates_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Female_Alates$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.00      0.00      0.32      7.64      1.63 104787.61  
'''

boxplot(Female_Alates$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<2)

hist(Female_Alates$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Female_Alates_clean <- Female_Alates[Female_Alates$TPM>0,]
summary(Female_Alates_clean$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.01      0.25      0.74     10.57      2.91 104787.61  
'''
boxplot(Female_Alates_clean$TPM, ylim=c(0,10))
hist(log2(Female_Alates$TPM))
hist(Female_Alates_clean$TPM, breaks =1000000, xlim=c(0,20))

quantile(Female_Alates_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
    1%     5%    10%    25%    30%    90%    95%    99% 
 0.040  0.080  0.120  0.250  0.310  9.470 18.753 91.470 
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Female_Alates_clean$TPM[Female_Alates_clean$TPM > 0.12 & Female_Alates_clean$TPM < 9.47]
m.Female_Alates <- mean(no.ext)
# [1] 1.656052
abline(v=m.Female_Alates, col="red")
md.Female_Alates <- median(no.ext)
# [1] 0.75
abline(v=md.Female_Alates, col="blue")

## FC
fold.Female_Alates <- Female_Alates_clean$TPM/m.Female_Alates
plot(fold.Female_Alates,col = ifelse(fold.Female_Alates < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Female_Alates)
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)

hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (2+(abs(median(FC.log)))))
# 2 = [1] 6017
# 4 = [1] 1454

# Total de down
sum(FC.log < -(2+(abs(median(FC.log)))))
# 2 = [1] 17314
# 4 = [1] 1319

# Total de DE
sum(FC.log > (2+(abs(median(FC.log))))) + sum(FC.log < -(2+(abs(median(FC.log)))))
# 2 = [1] 23331
# 4 = [1] 2773
# 5 = [1] 1034

## Making a table of interesting outputs
Female_Alates.res <- data.table(transcript_id=Female_Alates_clean$transcript_id,TPM=Female_Alates_clean$TPM,FPKM=Female_Alates_clean$FPKM,
                           fold.change=fold.Female_Alates,log2.FC=FC.log)

## Saving outputs

## Expressed transcripts IDs
write.table(Female_Alates_clean$transcript_id, "Female_Alates_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Female_Alates.res, "Female_Alates_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Female_Alates.resUP <- Female_Alates.res[Female_Alates.res$log2.FC > (4+(abs(median(FC.log)))),]
Female_Alates.resDOWN <- Female_Alates.res[Female_Alates.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Female_Alates.resUP, "Female_Alates_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Female_Alates.resDOWN, "Female_Alates_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Female_Alates.resUP$transcript_id, "Female_Alates_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Female_Alates.resDOWN$transcript_id, "Female_Alates_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)


######

## Female_Aspirants
rm(list = ls())
## Open the counts
Female_Aspirants <- read.table("Female_Aspirants_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Female_Aspirants$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.000     0.080     0.490     7.641     2.370 27812.480  
'''

boxplot(Female_Aspirants$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<2)

hist(Female_Aspirants$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Female_Aspirants_clean <- Female_Aspirants[Female_Aspirants$TPM>0,]
summary(Female_Aspirants_clean$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.010     0.320     0.920     9.828     3.680 27812.480 
'''
boxplot(Female_Aspirants_clean$TPM, ylim=c(0,10))
hist(log2(Female_Aspirants$TPM))
hist(Female_Aspirants_clean$TPM, breaks =100000, xlim=c(0,20))

quantile(Female_Aspirants_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
     1%      5%     10%     25%     30%     90%     95%     99% 
 0.0500  0.1000  0.1500  0.3200  0.3900 11.6300 22.1300 97.7304
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Female_Aspirants_clean$TPM[Female_Aspirants_clean$TPM > 0.15 & Female_Aspirants_clean$TPM < 11.63]
m.Female_Aspirants <- mean(no.ext)
# [1] 2.080338
abline(v=m.Female_Aspirants, col="red")
md.Female_Aspirants <- median(no.ext)
# [1] 0.93
abline(v=md.Female_Aspirants, col="blue")

## FC
fold.Female_Aspirants <- Female_Aspirants_clean$TPM/m.Female_Aspirants
plot(fold.Female_Aspirants,col = ifelse(fold.Female_Aspirants < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Female_Aspirants)
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)

hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (4+(abs(median(FC.log)))))
# 2 = [1] 6056
# 4 = [1] 1310

# Total de down
sum(FC.log < -(4+(abs(median(FC.log)))))
# 2 = [1] 17477
# 4 = [1] 1510

# Total de DE
sum(FC.log > (2+(abs(median(FC.log))))) + sum(FC.log < -(2+(abs(median(FC.log)))))
# 2 = [1] 23533
# 4 = [1] 2820
# 5 = [1] 889

## Making a table of interesting outputs
Female_Aspirants.res <- data.table(transcript_id=Female_Aspirants_clean$transcript_id,TPM=Female_Aspirants_clean$TPM,FPKM=Female_Aspirants_clean$FPKM,
                                   fold.change=fold.Female_Aspirants,log2.FC=FC.log)

## Saving outputs

## Expressed transcripts IDs
write.table(Female_Aspirants_clean$transcript_id, "Female_Aspirants_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Female_Aspirants.res, "Female_Aspirants_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Female_Aspirants.resUP <- Female_Aspirants.res[Female_Aspirants.res$log2.FC > (4+(abs(median(FC.log)))),]
Female_Aspirants.resDOWN <- Female_Aspirants.res[Female_Aspirants.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Female_Aspirants.resUP, "Female_Aspirants_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Female_Aspirants.resDOWN, "Female_Aspirants_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Female_Aspirants.resUP$transcript_id, "Female_Aspirants_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Female_Aspirants.resDOWN$transcript_id, "Female_Aspirants_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)



######

## Female_Aspirants_s
rm(list = ls())
## Open the counts
Female_Aspirants_s <- read.table("Female_Aspirants_s_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Female_Aspirants_s$TPM)
'''
    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.000     0.110     0.550      7.641     2.500 30473.060 
'''

boxplot(Female_Aspirants_s$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<3)

hist(Female_Aspirants_s$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Female_Aspirants_s_clean <- Female_Aspirants_s[Female_Aspirants_s$TPM>0,]
summary(Female_Aspirants_s_clean$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.010     0.330     0.960     9.588     3.660 30473.060 
'''
boxplot(Female_Aspirants_s_clean$TPM, ylim=c(0,10))
hist(log2(Female_Aspirants_s$TPM))
hist(Female_Aspirants_s_clean$TPM, breaks =1000000, xlim=c(0,20))

quantile(Female_Aspirants_s_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
     1%      5%     10%     25%     30%     90%     95%     99% 
 0.0500  0.1000  0.1500  0.3300  0.4100 11.6000 22.2140 98.4032 
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Female_Aspirants_s_clean$TPM[Female_Aspirants_s_clean$TPM > 0.15 & Female_Aspirants_s_clean$TPM < 11.60]
m.Female_Aspirants_s <- mean(no.ext)
# [1] 2.080542
abline(v=m.Female_Aspirants_s, col="red")
md.Female_Aspirants_s <- median(no.ext)
# [1] 0.97
abline(v=md.Female_Aspirants_s, col="blue")

## FC
fold.Female_Aspirants_s <- Female_Aspirants_s_clean$TPM/m.Female_Aspirants_s
plot(fold.Female_Aspirants_s,col = ifelse(fold.Female_Aspirants_s < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Female_Aspirants_s)
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)
hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (2+(abs(median(FC.log)))))
# 2 = [1] 6526
# 4 = [1] 1438

# Total de down
sum(FC.log < -(2+(abs(median(FC.log)))))
# 2 = [1] 18364
# 4 = [1] 1630

# Total de DE
sum(FC.log > (2+(abs(median(FC.log))))) + sum(FC.log < -(2+(abs(median(FC.log)))))
# 2 = [1] 24890
# 4 = [1] 3068
# 5 = [1] 947


## Making a table of interesting outputs
Female_Aspirants_s.res <- data.table(transcript_id=Female_Aspirants_s_clean$transcript_id,TPM=Female_Aspirants_s_clean$TPM,FPKM=Female_Aspirants_s_clean$FPKM,
                                   fold.change=fold.Female_Aspirants_s,log2.FC=FC.log)



## Saving outputs

## Expressed transcripts IDs
write.table(Female_Aspirants_s_clean$transcript_id, "Female_Aspirants_s_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Female_Aspirants_s.res, "Female_Aspirants_s_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Female_Aspirants_s.resUP <- Female_Aspirants_s.res[Female_Aspirants_s.res$log2.FC > (4+(abs(median(FC.log)))),]
Female_Aspirants_s.resDOWN <- Female_Aspirants_s.res[Female_Aspirants_s.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Female_Aspirants_s.resUP, "Female_Aspirants_s_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Female_Aspirants_s.resDOWN, "Female_Aspirants_s_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Female_Aspirants_s.resUP$transcript_id, "Female_Aspirants_s_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Female_Aspirants_s.resDOWN$transcript_id, "Female_Aspirants_s_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)



######

## Male_Alates
rm(list = ls())
## Open the counts
Male_Alates <- read.table("Male_Alates_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Male_Alates$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.00e+00 0.00e+00 3.60e-01 7.64e+00 1.72e+00 1.03e+05 
'''

boxplot(Male_Alates$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<2)

hist(Male_Alates$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Male_Alates_clean <- Male_Alates[Male_Alates$TPM>0,]
summary(Male_Alates_clean$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.01      0.31      0.86     10.92      3.13 103041.68 
'''
boxplot(Male_Alates_clean$TPM, ylim=c(0,10))
hist(log2(Male_Alates$TPM))
hist(Male_Alates_clean$TPM, breaks =1000000, xlim=c(0,20))

quantile(Male_Alates_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
     1%      5%     10%     25%     30%     90%     95%     99% 
 0.0500  0.1000  0.1500  0.3100  0.3800  9.7900 19.2860 93.3152
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Male_Alates_clean$TPM[Male_Alates_clean$TPM > 0.15 & Male_Alates_clean$TPM < 9.79]
m.Male_Alates <- mean(no.ext)
# [1] 1.788634
abline(v=m.Male_Alates, col="red")
md.Male_Alates <- median(no.ext)
# [1] 0.87
abline(v=md.Male_Alates, col="blue")

## FC
fold.Male_Alates <- Male_Alates_clean$TPM/m.Male_Alates
plot(fold.Male_Alates,col = ifelse(fold.Male_Alates < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Male_Alates)
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)

hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (4+(abs(median(FC.log)))))
# 2 = [1] 6007
# 4 = [1] 1406

# Total de down
sum(FC.log < -(4+(abs(median(FC.log)))))
# 2 = [1] 14877
# 4 = [1] 1116

# Total de DE
sum(FC.log > (2+(abs(median(FC.log))))) + sum(FC.log < -(2+(abs(median(FC.log)))))
# 2 = [1] 20844
# 4 = [1] 2522
# 5 = [1] 912

## Making a table of interesting outputs
Male_Alates.res <- data.table(transcript_id=Male_Alates_clean$transcript_id,TPM=Male_Alates_clean$TPM,FPKM=Male_Alates_clean$FPKM,
                                     fold.change=fold.Male_Alates,log2.FC=FC.log)


## Saving outputs

## Expressed transcripts IDs
write.table(Male_Alates_clean$transcript_id, "Male_Alates_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Male_Alates.res, "Male_Alates_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Male_Alates.resUP <- Male_Alates.res[Male_Alates.res$log2.FC > (4+(abs(median(FC.log)))),]
Male_Alates.resDOWN <- Male_Alates.res[Male_Alates.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Male_Alates.resUP, "Male_Alates_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Male_Alates.resDOWN, "Male_Alates_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Male_Alates.resUP$transcript_id, "Male_Alates_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Male_Alates.resDOWN$transcript_id, "Male_Alates_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)





######

## Neotenic_Queens_p
rm(list = ls())
## Open the counts
Neotenic_Queens_p <- read.table("Neotenic_Queens_p_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Neotenic_Queens_p$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.000     0.000     0.260     7.641     1.550 28737.140  
'''

boxplot(Neotenic_Queens_p$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<2)

hist(Neotenic_Queens_p$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Neotenic_Queens_p_clean <- Neotenic_Queens_p[Neotenic_Queens_p$TPM>0,]
summary(Neotenic_Queens_p_clean$TPM)
'''
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.01     0.25     0.80    11.28     3.10 28737.14
'''
boxplot(Neotenic_Queens_p_clean$TPM, ylim=c(0,10))
hist(log2(Neotenic_Queens_p$TPM))
hist(Neotenic_Queens_p_clean$TPM, breaks =1000000, xlim=c(0,20))

quantile(Neotenic_Queens_p_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
   1%    5%   10%   25%   30%   90%   95%   99% 
 0.03  0.07  0.11  0.25  0.32 10.15 20.14 96.52 
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Neotenic_Queens_p_clean$TPM[Neotenic_Queens_p_clean$TPM > 0.11 & Neotenic_Queens_p_clean$TPM < 10.15]
m.Neotenic_Queens_p <- mean(no.ext)
# [1] 1.777298
abline(v=m.Neotenic_Queens_p, col="red")
md.Neotenic_Queens_p <- median(no.ext)
# [1] 0.82
abline(v=md.Neotenic_Queens_p, col="blue")

## FC
fold.Neotenic_Queens_p <- Neotenic_Queens_p_clean$TPM/m.Neotenic_Queens_p
plot(fold.Neotenic_Queens_p,col = ifelse(fold.Neotenic_Queens_p < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Neotenic_Queens_p)
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)

hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (4+(abs(mean(FC.log)))))
# 2 = [1] 6676
# 4 = [1] 1569

# Total de down
sum(FC.log < -(4+(abs(mean(FC.log)))))
# 2 = [1] 21142
# 4 = [1] 3024

# Total de DE
sum(FC.log > (5+(abs(mean(FC.log))))) + sum(FC.log < -(5+(abs(mean(FC.log)))))
# 2 = [1] 27818
# 4 = [1] 4593
# 5 = [1] 1284

## Making a table of interesting outputs
Neotenic_Queens_p.res <- data.table(transcript_id=Neotenic_Queens_p_clean$transcript_id,TPM=Neotenic_Queens_p_clean$TPM,FPKM=Neotenic_Queens_p_clean$FPKM,
                              fold.change=fold.Neotenic_Queens_p,log2.FC=FC.log)


## Saving outputs

## Expressed transcripts IDs
write.table(Neotenic_Queens_p_clean$transcript_id, "Neotenic_Queens_p_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Neotenic_Queens_p.res, "Neotenic_Queens_p_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Neotenic_Queens_p.resUP <- Neotenic_Queens_p.res[Neotenic_Queens_p.res$log2.FC > (4+(abs(median(FC.log)))),]
Neotenic_Queens_p.resDOWN <- Neotenic_Queens_p.res[Neotenic_Queens_p.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Neotenic_Queens_p.resUP, "Neotenic_Queens_p_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Neotenic_Queens_p.resDOWN, "Neotenic_Queens_p_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Neotenic_Queens_p.resUP$transcript_id, "Neotenic_Queens_p_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Neotenic_Queens_p.resDOWN$transcript_id, "Neotenic_Queens_p_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)






######

## Primary_Queens
rm(list = ls())
## Open the counts
Primary_Queens <- read.table("Primary_Queens_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Primary_Queens$TPM)
'''
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.00     0.00     0.06     7.64     0.76 73617.35  
'''

boxplot(Primary_Queens$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<1)


hist(Primary_Queens$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Primary_Queens_clean <- Primary_Queens[Primary_Queens$TPM>0,]
summary(Primary_Queens_clean$TPM)
'''
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.01     0.19     0.65    14.36     2.59 73617.35 
'''
boxplot(Primary_Queens_clean$TPM, ylim=c(0,10))
hist(log2(Primary_Queens$TPM))
hist(Primary_Queens_clean$TPM, breaks =1000000, xlim=c(0,20))

quantile(Primary_Queens_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
     1%      5%     10%     25%     30%     90%     95%     99% 
 0.0200  0.0500  0.0800  0.1900  0.2400  8.2400 16.4500 82.6226 
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Primary_Queens_clean$TPM[Primary_Queens_clean$TPM > 0.08 & Primary_Queens_clean$TPM < 8.24]
m.Primary_Queens <- mean(no.ext)
# [1] 1.461006
abline(v=m.Primary_Queens, col="red")
md.Primary_Queens <- median(no.ext)
# [1] 0.67
abline(v=md.Primary_Queens, col="blue")

## FC
fold.Primary_Queens <- Primary_Queens_clean$TPM/m.Primary_Queens
plot(fold.Primary_Queens,col = ifelse(fold.Primary_Queens < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Primary_Queens)
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)

hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (4+(abs(median(FC.log)))))
# 2 = [1] 4393
# 4 = [1] 1069

# Total de down
sum(FC.log < -(4+(abs(median(FC.log)))))
# 2 = [1] 15470
# 4 = [1] 2894

# Total de DE
sum(FC.log > (5+(abs(median(FC.log))))) + sum(FC.log < -(5+(abs(median(FC.log)))))
# 2 = [1] 19863
# 4 = [1] 3963
# 5 = [1] 1336

## Making a table of interesting outputs
Primary_Queens.res <- data.table(transcript_id=Primary_Queens_clean$transcript_id,TPM=Primary_Queens_clean$TPM,FPKM=Primary_Queens_clean$FPKM,
                                    fold.change=fold.Primary_Queens,log2.FC=FC.log)



## Saving outputs

## Expressed transcripts IDs
write.table(Primary_Queens_clean$transcript_id, "Primary_Queens_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Primary_Queens.res, "Primary_Queens_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Primary_Queens.resUP <- Primary_Queens.res[Primary_Queens.res$log2.FC > (4+(abs(median(FC.log)))),]
Primary_Queens.resDOWN <- Primary_Queens.res[Primary_Queens.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Primary_Queens.resUP, "Primary_Queens_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Primary_Queens.resDOWN, "Primary_Queens_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Primary_Queens.resUP$transcript_id, "Primary_Queens_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Primary_Queens.resDOWN$transcript_id, "Primary_Queens_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)






######

## Soldiers 
rm(list = ls())
## Open the counts
Soldiers <- read.table("Soldiers_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Soldiers$TPM)
'''
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.000     0.000     0.220     7.641     1.350 29453.060
'''

boxplot(Soldiers$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<1.5)

hist(Soldiers$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Soldiers_clean <- Soldiers[Soldiers$TPM>0,]
summary(Soldiers_clean$TPM)
'''
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.01     0.19     0.63    10.98     2.65 29453.06 
'''
boxplot(Soldiers_clean$TPM, ylim=c(0,10))
hist(log2(Soldiers$TPM))
hist(Soldiers_clean$TPM, breaks =1000000, xlim=c(0,20))

quantile(Soldiers_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
      1%       5%      10%      25%      30%      90%      95%      99% 
  0.0300   0.0600   0.0900   0.1900   0.2400   8.8400  18.2460 103.0704  
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Soldiers_clean$TPM[Soldiers_clean$TPM > 0.09 & Soldiers_clean$TPM < 8.84]
m.Soldiers <- mean(no.ext)
# [1] 1.514695
abline(v=m.Soldiers, col="red")
md.Soldiers <- median(no.ext)
# [1] 0.65
abline(v=md.Soldiers, col="blue")

## FC
fold.Soldiers <- Soldiers_clean$TPM/m.Soldiers
plot(fold.Soldiers,col = ifelse(fold.Soldiers < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Soldiers)
##FC.log <- log2(Soldiers_clean$TPM)/mean(log2(Soldiers_clean$TPM))
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)

hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (4+(abs(median(FC.log)))))
# 2 = [1] 5666
# 4 = [1] 1529

# Total de down
sum(FC.log < -(4+(abs(median(FC.log)))))
# 2 = [1] 18975
# 4 = [1] 1605

# Total de DE
sum(FC.log > (2+(abs(median(FC.log))))) + sum(FC.log < -(2+(abs(median(FC.log)))))
# 2 = [1] 24641
# 4 = [1] 3134
# 5 = [1] 948

## Making a table of interesting outputs
Soldiers.res <- data.table(transcript_id=Soldiers_clean$transcript_id,TPM=Soldiers_clean$TPM,FPKM=Soldiers_clean$FPKM,
                                 fold.change=fold.Soldiers,log2.FC=FC.log)


## Saving outputs

## Expressed transcripts IDs
write.table(Soldiers_clean$transcript_id, "Soldiers_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Soldiers.res, "Soldiers_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Soldiers.resUP <- Soldiers.res[Soldiers.res$log2.FC > (4+(abs(median(FC.log)))),]
Soldiers.resDOWN <- Soldiers.res[Soldiers.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Soldiers.resUP, "Soldiers_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Soldiers.resDOWN, "Soldiers_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Soldiers.resUP$transcript_id, "Soldiers_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Soldiers.resDOWN$transcript_id, "Soldiers_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)




######

## Workers 
rm(list = ls())
## Open the counts
Workers <- read.table("Workers_bowtie2_RSEM-trinity/RSEM.isoforms.results", stringsAsFactors = FALSE, header = T)

summary(Workers$TPM)
'''
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.00     0.14     0.57     7.64     2.23 53543.06 
'''

boxplot(Workers$TPM, ylim=c(0,10))
## Many transprits have a very low coverage (<2)

hist(Workers$TPM, breaks =1000000, xlim=c(0,20))

## Removing transcripts with 0 TPM
Workers_clean <- Workers[Workers$TPM>0,]
summary(Workers_clean$TPM)
'''
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.01     0.33     0.89     9.31     3.02 53543.06 
'''
boxplot(Workers_clean$TPM, ylim=c(0,10))
hist(log2(Workers$TPM))
hist(Workers_clean$TPM, breaks =1000000, xlim=c(0,20))

quantile(Workers_clean$TPM, probs = c(0.01,0.05,0.1,0.25,0.3, 0.90, 0.95,0.99))
'''
     1%      5%     10%     25%     30%     90%     95%     99% 
 0.0500  0.1000  0.1500  0.3300  0.4000  9.1400 17.4500 83.3464
'''

## Calculating fold change (FC)
## Remove extreme values to get the mean value (between 10 and 90% quantiles)
no.ext <- Workers_clean$TPM[Workers_clean$TPM > 0.15 & Workers_clean$TPM < 9.14]
m.Workers <- mean(no.ext)
# [1] 1.742714
abline(v=m.Workers, col="red")
md.Workers <- median(no.ext)
# [1] 0.9
abline(v=md.Workers, col="blue")

## FC
fold.Workers <- Workers_clean$TPM/m.Workers
plot(fold.Workers,col = ifelse(fold.Workers < 1,"blue","red"), pch = 19)

## Transforming to log2 (more meaninful)
FC.log <- log2(fold.Workers)
plot(FC.log, col = ifelse(FC.log < -(4+(abs(median(FC.log)))),"blue", ifelse(FC.log >(4+(abs(median(FC.log)))), "red", "black")), pch = 19)

hist(FC.log, breaks=100)
abline(v=mean(FC.log), col="red")
abline(v=median(FC.log), col="blue")

# Total de UP correcting for distribution bias
sum(FC.log > (4+(abs(median(FC.log)))))
# 2 = [1] 7082
# 4 = [1] 1598

# Total de down
sum(FC.log < -(4+(abs(median(FC.log)))))
# 2 = [1] 18204
# 4 = [1] 1679

# Total de DE
sum(FC.log > (5+(abs(median(FC.log))))) + sum(FC.log < -(5+(abs(median(FC.log)))))
# 2 = [1] 25286
# 4 = [1] 3277
# 5 = [1] 1078

## Making a table of interesting outputs
Workers.res <- data.table(transcript_id=Workers_clean$transcript_id,TPM=Workers_clean$TPM,FPKM=Workers_clean$FPKM,
                           fold.change=fold.Workers,log2.FC=FC.log)


## Saving outputs

## Expressed transcripts IDs
write.table(Workers_clean$transcript_id, "Workers_full_transcrp_expressed_id.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)
## Full counts
write.table(Workers.res, "Workers_counts.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## DF transcrits
Workers.resUP <- Workers.res[Workers.res$log2.FC > (4+(abs(median(FC.log)))),]
Workers.resDOWN <- Workers.res[Workers.res$log2.FC < -(4+(abs(median(FC.log)))),]
write.table(Workers.resUP, "Workers_UP-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
write.table(Workers.resDOWN, "Workers_DOWN-counts.txt", row.names=FALSE,col.names=TRUE, quote=FALSE)
# IDs
write.table(Workers.resUP$transcript_id, "Workers_UP-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(Workers.resDOWN$transcript_id, "Workers_DOWN-IDs.txt", row.names=FALSE,col.names=FALSE, quote=FALSE)


