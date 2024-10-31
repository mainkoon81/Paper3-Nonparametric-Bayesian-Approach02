################################### Data lab ###################################
#remotes::install_github("dutangc/CASdatasets/pkg")

library(CASdatasets) #how to install?

library(tidyverse)
library(car) # Companion to Applied Regression": to Recodes a numeric, character vector as per specifications.
library(varhandle) # to use "unfactor()"

# check available datasetsssss
data()

#> Data sets in package ‘CASdatasets’:
#PnCdemand::::::::::::NA:::::::::::::Property and casualty insurance demand
#ausautoBI8999                       Automobile bodily injury claim dataset in Australia*****
#ausprivauto0405                     Automobile claim datasets in Australia*****

#brvehins1a :::::::::: NA :::::::::: Two Brazilian datasets for vehicle insurance******
#brvehins1b :::::::::: NA :::::::::: Two Brazilian datasets for vehicle insurance******
#brvehins1c :::::::::: NA :::::::::: Two Brazilian datasets for vehicle insurance******
#brvehins1d :::::::::: NA :::::::::: Two Brazilian datasets for vehicle insurance******
#brvehins1e :::::::::: NA :::::::::: Two Brazilian datasets for vehicle insurance******

#freMTPL2freq                        French Motor Third-Part Liability datasets
#freMTPL2sev                         French Motor Third-Part Liability datasets
#freMTPLfreq                         French Motor Third-Part Liability datasets
#freMTPLsev                          French Motor Third-Part Liability datasets

#usautoBI                            Automobile bodily injury claim dataset


# data(PnCdemand)
# head(PnCdemand)
# summary(PnCdemand) # PecLoss:(0-105), x1, x2

# data(ausautoBI8999 )
# head(ausautoBI8999 )
# summary(ausautoBI8999 ) # (10-4,500k), no NA, x1, x2
# 
# data(ausprivauto0405 )    
# head(ausprivauto0405 )
# summary(ausprivauto0405 ) # (0-50k), no NA, x1, x2

##########################################
# data(brvehins1a)
# head(brvehins1a)
# summary(brvehins1a)    # SumInsAvg:(0-860k), x1, x2
# data(brvehins1b)
# head(brvehins1b)
# summary(brvehins1b)    # (0-970k), x1, x2
# data(brvehins1c)
# head(brvehins1c)
# summary(brvehins1c)    # (0-1100k), x1, x2
# data(brvehins1d)
# head(brvehins1d)
# summary(brvehins1d)    # (0-1307k), x1, x2
# data(brvehins1e)
# head(brvehins1e)
# summary(brvehins1e)    # (0-1250k), x1, x2


# Number of exhibited vehicles (see exhibition concept above)
# PREMIO1 – Sum of the premium values, weighted by the exposure of each
#policy
# EXPOSICAO2 – Field not used
# PREMIO2 – Unused field
# IS_MEDIA - Average Insured Amounts of the policies included in the
#grouping defined by the chosen key, weighted by the exposure of each one
#from them

#(c)
data(brvehins1c)
head(brvehins1c)
sum( brvehins1c$ClaimAmountRob > 0 )
sum( brvehins1c$ClaimAmountPartColl > 0 )
sum( brvehins1c$ClaimNbTotColl > 0 )
sum( brvehins1c$ClaimAmountFire > 0 )
sum( brvehins1c$ClaimAmountOther > 0 )

brvehins1c$AggClaim <- brvehins1c$ClaimAmountRob + brvehins1c$ClaimAmountPartColl + brvehins1c$ClaimNbTotColl + brvehins1c$ClaimAmountFire + brvehins1c$ClaimAmountOther   

head(brvehins1c)

summary(brvehins1c)

sum( brvehins1c$AggClaim > 0 ) # 72831 obv

df.brazil <- brvehins1c[brvehins1c$AggClaim > 0, ]

summary(df.brazil$Gender)
df.brazil$Gender[is.na(df.brazil$Gender)] <- "Corporate"



summary(df.brazil)
cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="pearson" ) # cor: 0.4775715
cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="spearman" ) # cor: 0.2528002
cor( df.brazil$AggClaim, df.brazil$PremTotal, method="pearson" ) # cor: 0.5232885
cor( df.brazil$AggClaim, df.brazil$PremTotal, method="spearman" ) # cor: 0.2615761




# #(a)
# data(brvehins1a)
# 
# brvehins1a$AggClaim <- brvehins1a$ClaimAmountRob + brvehins1a$ClaimAmountPartColl + brvehins1a$ClaimNbTotColl + brvehins1a$ClaimAmountFire + brvehins1a$ClaimAmountOther   
# 
# head(brvehins1a)
# 
# summary(brvehins1a)
# 
# sum( brvehins1a$AggClaim > 0 ) # 73161 obv
# 
# df.brazil <- brvehins1a[brvehins1a$AggClaim > 0, ]
# summary(df.brazil)
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="pearson" ) # cor: 0.4281937
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="spearman" ) # cor: 0.2575091
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="pearson" ) # cor: 0.4591074
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="spearman" ) # cor: 0.4591074
# 
# 
# 
# 
# #(b)
# data(brvehins1b)
# 
# brvehins1b$AggClaim <- brvehins1b$ClaimAmountRob + brvehins1b$ClaimAmountPartColl + brvehins1b$ClaimNbTotColl + brvehins1b$ClaimAmountFire + brvehins1b$ClaimAmountOther   
# 
# head(brvehins1b)
# 
# summary(brvehins1a)
# 
# sum( brvehins1b$AggClaim > 0 ) # 72445 obv
# 
# df.brazil <- brvehins1b[brvehins1b$AggClaim > 0, ]
# summary(df.brazil)
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="pearson" ) # cor: 0.3919972
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="spearman" ) # cor: 0.2596957
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="pearson" ) # cor: 0.4376457
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="spearman" ) # cor: 0.2713166
# 
# 
# 
# #(d)
# data(brvehins1d)
# 
# brvehins1d$AggClaim <- brvehins1d$ClaimAmountRob + brvehins1d$ClaimAmountPartColl + brvehins1d$ClaimNbTotColl + brvehins1d$ClaimAmountFire + brvehins1d$ClaimAmountOther   
# 
# head(brvehins1d)
# 
# summary(brvehins1d)
# 
# sum( brvehins1d$AggClaim > 0 ) # 72479 obv
# 
# df.brazil <- brvehins1d[brvehins1d$AggClaim > 0, ]
# summary(df.brazil)
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="pearson" ) # cor: 0.3809636
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="spearman" ) # cor: 0.2620557
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="pearson" ) # cor: 0.389887
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="spearman" ) # cor: 0.2734049
# 
# 
# 
# #(e)
# data(brvehins1e)
# 
# brvehins1e$AggClaim <- brvehins1e$ClaimAmountRob + brvehins1e$ClaimAmountPartColl + brvehins1e$ClaimNbTotColl + brvehins1e$ClaimAmountFire + brvehins1e$ClaimAmountOther   
# 
# head(brvehins1e)
# 
# summary(brvehins1e)
# 
# sum( brvehins1e$AggClaim > 0 ) # 72160 obv
# 
# df.brazil <- brvehins1e[brvehins1e$AggClaim > 0, ]
# summary(df.brazil)
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="pearson" ) # cor: 0.3521279
# cor( df.brazil$AggClaim, df.brazil$ExposTotal, method="spearman" ) # cor: 0.2558685
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="pearson" ) # cor: 0.3730427
# cor( df.brazil$AggClaim, df.brazil$PremTotal, method="spearman" ) # cor: 0.2666676





##########################################

## [candidate 01.] -------------------------------------------------------------
library(CASdatasets)

## > French Motor Third-Part Liability datasets
data(freMTPL2sev)
head(freMTPL2sev)
summary(freMTPL2sev)

data(freMTPL2freq)
head(freMTPL2freq)
summary(freMTPL2freq)    # (0-970k), no NA, x1, x2

df_MTPL2 <- merge(freMTPL2freq, freMTPL2sev, by = c("IDpol", "IDpol"))

summary(df_MTPL2) # only 26444 obv
head(df_MTPL2)
cor(df_MTPL2$ClaimAmount, df_MTPL2$Exposure, method="pearson") # rho: 0.02
cor(df_MTPL2$ClaimAmount, df_MTPL2$Exposure, method="spearman") # rho: 0.1


data(freMTPLsev)
head(freMTPLsev)
summary(freMTPLsev)
data(freMTPLfreq)
head(freMTPLfreq)
summary(freMTPLfreq)    # (0-970k), no NA, x1, x2

df_MTPL <- merge(freMTPLfreq, freMTPLsev, by = c("PolicyID", "PolicyID"))
summary(df_MTPL) # only 16181 obv
head(df_MTPL)
cor(df_MTPL$ClaimAmount, df_MTPL$Exposure, method="pearson") # rho: 0.03
cor(df_MTPL$ClaimAmount, df_MTPL$Exposure, method="spearman") # rho: 0.1




## [candidate 02. categorical covariate only] 
# library(CASdatasets)
# 
# data(usautoBI)
# head(usautoBI)
# summary(usautoBI)       # (0-1000), x2
# quantile(usautoBI$LOSS, c(0.025, 0.975))
# hist(usautoBI$LOSS) # only 1340 obv
# cor(usautoBI$LOSS, usautoBI$)

# data(ausautoBI8999)
# head(ausautoBI8999)
# summary(ausautoBI8999)



## [candidate 03.] -------------------------------------------------------------
library(insuranceData)
# library(CASdatasets)

## > Automobile claim datasets in Australia
# data(ausprivauto0405) # 2004-2005 one-year vehicle insurance policies  
data(dataCar) # only 67856 obv
head(dataCar)
summary(dataCar)
quantile(dataCar$claimcst0, c(0.025, 0.975))
hist(dataCar$claimcst0)
hist( dataCar[dataCar$claimcst0>0, ]$claimcst0 )
sum( dataCar$claimcst0>0 ) # only 4624 obv

df_car <- dataCar[dataCar$claimcst0>0, ]
cor(df_car$claimcst0, df_car$exposure, method="pearson") # cor:0.13
cor(df_car$claimcst0, df_car$exposure, method="spearman") # cor:0.04
head(df_car)





## [candidate 04.] -------------------------------------------------------------
library(CASdatasets)

## > AON Re Belgian dataset (fire losses)
data(beaonre) # only 1823 obv
head(beaonre)
summary(beaonre)

cor(beaonre$ClaimCost, beaonre$SumInsured, method="pearson") # cor:0.08
cor(beaonre$ClaimCost, beaonre$SumInsured, method="spearman") # cor: 0.33

df_fire <- beaonre





## [candidate 05.] -------------------------------------------------------------
library(CASdatasets)

## > Brazilian premium/claim per region and type of insurance coverage 
data(braggclaim)
head(braggclaim)
data(braggprem)
head(braggprem)

df_brag <- cbind(braggprem, braggclaim$AggClaim) # only 205 obv
head(df_brag)
summary(df_brag)
df_brag$AggClaim <- df_brag$`braggclaim$AggClaim`
df_brag$AggClaim <- as.numeric(df_brag$AggClaim)

cor( df_brag$AggClaim, df_brag$SumInsAvg, method="pearson") # cor:0
cor( df_brag$AggClaim, df_brag$SumInsAvg, method="spearman") # cor:0.1
cor( df_brag$AggClaim, df_brag$PremAvg, method="pearson") # cor:0.46
cor( df_brag$AggClaim, df_brag$PremAvg, method="spearman") # cor:0.43
cor( df_brag$PremAvg, df_brag$SumInsAvg, method="pearson") # cor:0.20
cor( df_brag$PremAvg, df_brag$SumInsAvg, method="spearman") # cor: 0.28




## [candidate 06.] -------------------------------------------------------------
library(CASdatasets)

# Swedish Motor Insurance dataset
# This dataset contains motor insurance data collected in 1977 in Sweden by the Swedish Committee
#on the Analysis of Risk Premium. Records contains individuals characteristics in addition to claim
#counts and severities. 

# > Kilometres: Distance driven by a vehicle, grouped into 5 categories.
# > Zone: Graphic zone of a vehicle, grouped into 7 categories.
# > Bonus: Driver claim experience, grouped into 7 categories.
# > Make: The type of a vehicle, 9 categories
# > Insured: The number of policyholder years. A policyholder year is the fraction of the year that the
#policyholder has a contract with the issuing company.
# > Claims: Number of claims.
# > Payment: Sum of payments.

data(swautoins) # only 2182 obv
head(swautoins)
summary(swautoins)

cor( swautoins$Payment, swautoins$Insured, method="pearson" ) # cor: 0.9
cor( swautoins$Payment, swautoins$Insured, method="spearman" ) # cor: 0.9

df_swed <- swautoins




# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

#<<<< 154 OBV >>>>>

data(PnCdemand) 
head(PnCdemand)
summary(PnCdemand) # PecLoss:(0-105), x1, x2

df11 = select(.data=PnCdemand, Name, LegalSyst, RiskAversion, GenLiab); head(df11)
df11$Name = unfactor(df11$Name)
table(df11$Name)

test_df11 <- df11[df11$Name %in% c("Spain", "Switzerland","United Kingdom", "United States"), ] 
summary(test_df11)
str(test_df11)
head(test_df11)
which(is.na(test_df11$GenLiab))
test_df11[22, ]$GenLiab = 101.313
hist(test_df11$GenLiab)
summary(test_df11)


train_df11 <- df11[df11$Name %in% c("Australia","Austria","Belgium","Canada","Denmark","Finland","France","Greece","Iceland","Ireland","Italy","Japan","Luxembourg","Netherlands","New Zealand","Norway","Portugal","Sweden" ), ] 
summary(train_df11)
str(train_df11)

which(is.na(train_df11$GenLiab))
train_df11[10, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[50, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[53, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[54, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[85, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[92, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[93, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[94, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[95, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[96, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[97, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[98, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[99, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[100, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[101, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[102, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[103, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[104, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[105, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[120, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[121, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[122, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[123, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[124, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[125, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))
train_df11[126, ]$GenLiab = abs(rnorm(n=1, mean=2, sd=2.35))

hist(train_df11$GenLiab)
summary(train_df11)
head(train_df11)

# write.csv(train_df11, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/MK_Demand.train.csv", row.names=FALSE)
# write.csv(test_df11, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/MK_Demand.test.csv", row.names=FALSE)




################################################################################
############################# Medium data for paper3 ###########################
################################################################################
# ---------------------------------------------------------------------------- #
# Read in data
insurance = read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/LGPF_MAR.csv")
head(insurance)


insurance = select(.data=insurance, Year, Total_Losses, Fire5_1.original, Ln.Cov._1, Ln.Cov._2, Error_Indicator, Random_0.7_1.3); head(insurance)
colnames(insurance) = c("year","loss","protect_level","ln_cover","ln_cover_mismeas","err_indicator", "0.7_1.3"); head(insurance)

#table(insurance$err_indicator)
#   0    1 
#3879 1760 



# We use the data in the first 4 years, namely 2006-2009, to develop the model.
# keep the obv in the final year (2010) for validation purposes. OK. now, let's split train and test
train_df <- subset(x=insurance, subset=year<2010); head(train_df)
summary(train_df)

test_df <- subset(x=insurance, subset=year==2010); head(test_df)
summary(test_df)





















################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
############################## BIG data for paper3 #############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
################## Brazilian datasets for vehicle insurance ####################
head(df.brazil) #--------------------------------------------------------------

df1 = select(.data=df.brazil, VehYear, Area, StateAb, PremTotal, Gender, ExposTotal, SumInsAvg, AggClaim); head(df1)
summary(df1)

df1 <- df1[!is.na(df1$VehYear), ]
sum(df1$VehYear<1)
df1 <- df1[df1$VehYear>1, ]

df1$Gender <- recode(df1$Gender, " c('Male','Female')='Individual' ")
names(df1)[names(df1)=='Gender'] <- 'Affiliation'
df1$Affiliation <- ifelse( df1$Affiliation=="Individual", 1, 
                           ifelse(df1$Affiliation=="Corporate", 0, df1$Affiliation) )
#df1$StateAb = unfactor(df1$StateAb)

summary(df1)

# we don't need zero inflation in [ExposTotal], [SumInsAvg]
df1 <- subset(df1, subset=ExposTotal>0) 
summary(df1)

table(df1$Affiliation)
table(df1$Area)
df1$Area <- unfactor(df1$Area)

# df1$Area = factor(labels=c("Acre","Alagoas","Amapa","Amazonas","Bahia",
#                     "Blumenau e demais regioes","Brasilia","Ceara","Demais regioes","Espirito Santo",
#                     "F.Iguatu-Medianeira-Cascavel-Toledo", "Goias","Grande Campinas","Interior","Litoral Norte e Baixada Santista",
#                     "Maranhao","Mato Grosso","Mato Grosso do Sul","Met. Curitiba","Met. de Sao Paulo",
#                     "Met. do Rio de Janeiro","Met. Florianopolis e Sul", "Met. Porto Alegre e Caxias do Sul","Met.BH-Centro Oeste-Zona Mata-C. Vertentes","Oeste",
#                     "Para","Paraiba","Pernambuco","Piaui","Ribeirao Preto e Demais Mun. de Campinas",
#                     "Rio Grande do Norte","Rondonia","Roraima","Sergipe","Sudeste de Goias",
#                     "Sul","Tocantins","Triangulo mineiro","Vale do Aco-Norte-Vale Jequitinhonha","Vale do Paraiba e Ribeira"),
#                   levels=c("grp1","grp2","grp3","grp4","grp5","grp6","grp7","grp8","grp9","grp10",
#                            "grp11","grp12","grp13","grp14","grp15","grp16","grp17","grp18","grp19","grp20",
#                            "grp21","grp22","grp23","grp24","grp25","grp26","grp27","grp28","grp29","grp30",
#                            "grp31","grp32","grp33","grp34","grp35","grp36","grp37","grp38","grp39","grp40"))


df1$Area[df1$Area=="Acre"] <- "grp1"
df1$Area[df1$Area=="Alagoas"] <- "grp2"
df1$Area[df1$Area=="Amapa"] <- "grp3"
df1$Area[df1$Area=="Amazonas"] <- "grp4"
df1$Area[df1$Area=="Bahia"] <- "grp5"

df1$Area[df1$Area=="Blumenau e demais regioes"] <- "grp6"
df1$Area[df1$Area=="Brasilia"] <- "grp7"
df1$Area[df1$Area=="Ceara"] <- "grp8"
df1$Area[df1$Area=="Demais regioes"] <- "grp9"
df1$Area[df1$Area=="Espirito Santo"] <- "grp10"

df1$Area[df1$Area=="F.Iguatu-Medianeira-Cascavel-Toledo"] <- "grp11"
df1$Area[df1$Area=="Goias"] <- "grp12"
df1$Area[df1$Area=="Grande Campinas"] <- "grp13"
df1$Area[df1$Area=="Interior"] <- "grp14"
df1$Area[df1$Area=="Litoral Norte e Baixada Santista"] <- "grp15"

df1$Area[df1$Area=="Maranhao"] <- "grp16"
df1$Area[df1$Area=="Mato Grosso"] <- "grp17"
df1$Area[df1$Area=="Mato Grosso do Sul"] <- "grp18"
df1$Area[df1$Area=="Met. Curitiba"] <- "grp19"
df1$Area[df1$Area=="Met. de Sao Paulo"] <- "grp20"

df1$Area[df1$Area=="Met. do Rio de Janeiro"] <- "grp21"
df1$Area[df1$Area=="Met. Florianopolis e Sul"] <- "grp22"
df1$Area[df1$Area=="Met. Porto Alegre e Caxias do Sul"] <- "grp23"
df1$Area[df1$Area=="Met.BH-Centro Oeste-Zona Mata-C. Vertentes"] <- "grp24"
df1$Area[df1$Area=="Oeste"] <- "grp25"

df1$Area[df1$Area=="Para"] <- "grp26"
df1$Area[df1$Area=="Paraiba"] <- "grp27"
df1$Area[df1$Area=="Pernambuco"] <- "grp28"
df1$Area[df1$Area=="Piaui"] <- "grp29"
df1$Area[df1$Area=="Ribeirao Preto e Demais Mun. de Campinas"] <- "grp30"

df1$Area[df1$Area=="Rio Grande do Norte"] <- "grp31"
df1$Area[df1$Area=="Rondonia"] <- "grp32"
df1$Area[df1$Area=="Roraima"] <- "grp33"
df1$Area[df1$Area=="Sergipe"] <- "grp34"
df1$Area[df1$Area=="Sudeste de Goias"] <- "grp35"

df1$Area[df1$Area=="Sul"] <- "grp36"
df1$Area[df1$Area=="Tocantins"] <- "grp37"
df1$Area[df1$Area=="Triangulo mineiro"] <- "grp38"
df1$Area[df1$Area=="Vale do Aco-Norte-Vale Jequitinhonha"] <- "grp39"
df1$Area[df1$Area=="Vale do Paraiba e Ribeira"] <- "grp40"

table(df1$Area)

table(df1$StateAb)
#AC    AL    AM    AP    BA    CE    DF    ES    GO    MA    MG    MS    MT    PA    PB    PE    PI    PR    RJ 
#153   594   386    84  1908  1162  2170  1526  2208   561  5059   911   933   731   675  1407   487  4428  4818 
#RN    RO    RR    RS    SC    SE    SP    TO 
#883    72    55  4075  3302   606 20268   347 


cor( df1$AggClaim, df1$ExposTotal )
cor( df1$AggClaim, df1$PremTotal )
cor( df1$AggClaim, df1$SumInsAvg )
summary(df1)

set.seed(17)
sample <- sample( c(TRUE, FALSE), nrow(df1), replace=TRUE, prob=c(0.96, 0.04) )
train.df1  <- df1[sample, ]
test.df1   <- df1[!sample, ]


n.breaks = sqrt( nrow(train.df1) ) #******Rule of thumb
hist(train.df1$SumInsAvg, breaks=n.breaks)
hist(train.df1$AggClaim, breaks=n.breaks)


# write.csv(train.df1,
#           "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil_c.train.csv",
#           row.names=FALSE)
# write.csv(test.df1,
#           "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil_c.test.csv",
#           row.names=FALSE)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#################### French Motor Third-Part Liability dataset #################
head(df_MTPL)
summary(df_MTPL)

head(df_MTPL2)
summary(df_MTPL2)

df_fren1 <- select(.data=df_MTPL, PolicyID, Exposure, CarAge, Density, ClaimAmount); head(df_fren1)
df_fren2 <- select(.data=df_MTPL2, IDpol, Exposure, VehAge, Density, ClaimAmount); head(df_fren2)

summary(df_fren1)
summary(df_fren2)


df_fren1$CarAge <- ifelse( df_fren1$CarAge>=10, 1, 0)
table(df_fren1$CarAge)
#     0     1 
# 10342  5839  

df_fren2$VehAge <- ifelse( df_fren2$VehAge>=10, 1, 0)
table(df_fren2$VehAge)
#     0     1 
# 17457  8987  

colnames(df_fren2)[1] <- "PolicyID"
colnames(df_fren2)[3] <- "CarAge"


# let's combine...
df_fren <- rbind(df_fren1, df_fren2); head(df_fren)
summary(df_fren)


set.seed(17)
sample <- sample( c(TRUE, FALSE), nrow(df_fren), replace=TRUE, prob=c(0.95, 0.05) )
train.df1  <- df_fren[sample, ]
test.df1   <- df_fren[!sample, ]


n.breaks = sqrt( nrow(train.df1) ) #******Rule of thumb
hist(train.df1$ClaimAmount, breaks=n.breaks)
cor( df_fren$ClaimAmount, df_fren$Exposure )


# write.csv(train.df1,
#           "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFrench.train.csv",
#           row.names=FALSE)
# write.csv(test.df1,
#           "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFrench.test.csv",
#           row.names=FALSE)











# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#################### Automobile claim datasets in Australia ####################
head(df_car)
summary(df_car)



# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#################### AON Re Belgian dataset (fire losses) ######################
head(df_fire)
summary(df_fire)



# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
###### Brazilian premium/claim per region and type of insurance coverage #######
head(df_brag)
summary(df_brag)

df_brag <- select(.data=df_brag, RegionNb, ExpoAvg, PremAvg, SumInsAvg, StateAb, AggClaim); head(df_brag)



# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
####################### Swedish Motor Insurance dataset ########################

df_swed <- select(.data=swautoins, Kilometres, Zone, Bonus, Insured, Payment)
head(df_swed)
summary(df_swed)

sum( df_swed$Payment==0)

df_swed <- df_swed[df_swed$Payment>0, ]

colnames(df_swed)[3] <- "Experience"
colnames(df_swed)[4] <- "SumInsAvg"
colnames(df_swed)[5] <- "TotalLoss"
colnames(df_swed)[1] <- "Risk"

df_swed$Experience <- ifelse( df_swed$Experience>=4, 1, 0)
table(df_swed$Experience)
#   0    1 
# 744 1053
cor( df_swed$TotalLoss, df_swed$SumInsAvg, method="pearson" )


set.seed(17)
sample <- sample( c(TRUE, FALSE), nrow(df_swed), replace=TRUE, prob=c(0.85, 0.15) )
train.df1  <- df_swed[sample, ]
test.df1   <- df_swed[!sample, ]

n.breaks = sqrt( nrow(train.df1) ) #******Rule of thumb
hist(train.df1$SumInsAvg, breaks=n.breaks)
hist(train.df1$TotalLoss, breaks=n.breaks)


# write.csv(train.df1,
#           "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.train.csv",
#           row.names=FALSE)
# write.csv(test.df1,
#           "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.test.csv",
#           row.names=FALSE)




