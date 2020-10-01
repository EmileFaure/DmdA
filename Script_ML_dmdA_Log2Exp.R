########       Biogeography of dmdA      ############

######## Packages #########

library(ggplot2)
theme_set(theme_minimal()) 
library(data.table)
library(vegan)
library(gridExtra)
library(maps)
library(caret)
library(dplyr)
source("http://peterhaschke.com/Code/multiplot.R")

######## Data #########

Abund_G = read.table('/media/DATA/Emile/DMS/Abund_dmdA_metaG_log2.norm.tsv', sep = '\t', header = T)
Abund_T = read.table('/media/DATA/Emile/DMS/Abund_dmdA_metaT_log2.norm.tsv', sep = '\t', header = T)
Abund_Exp = read.table('/media/DATA/Emile/DMS/Abund_dmdA_exp_log2.norm.tsv', sep = '\t', header = T)
Envi = read.table('/media/DATA/Emile/DMS/Envi_Context_Log2_match.tsv', sep = '\t', header = T, dec = ',')

######## Exploration and preprocessing of envi data #########

summary(Envi)
# lower.size.fraction is constant, we get rid of it, upper.size.fraction should be categorical
Envi = Envi[,-5]
Envi$upper.size.fraction = as.factor(Envi$upper.size.fraction)
# PAR.PC has 92 NAs over 129 values, we get rid of it
Envi = Envi[,-which(names(Envi)=="PAR.PC")]
# Process environmental data for machine learning : center, scale, impute NA values from nearest neighboors, remove highly correlated and near zero variance variables
Envi_mltrans = preProcess(Envi, method = c("center","scale","knnImpute","corr","nzv"))
Envi_ml = predict(Envi_mltrans, Envi)
summary(Envi_ml)
names(Envi)[-which(names(Envi) %in% names(Envi_ml))] # Suppression of Carbon total, NO2NO3 and NO3

####### Exploration of Abundance data ########

# Create an abundance data set :
names(Abund_G)[2] = "K17486_G"
names(Abund_T)[2] = "K17486_T"
names(Abund_Exp)[2] = "K17486_Exp"

dmdA = merge(Abund_G, Abund_T)
dmdA = merge(dmdA, Abund_Exp)

dmdA_Envi = merge(dmdA,Envi)
dmdA_ml = merge(dmdA,Envi_ml)

ggplot() + geom_point(aes(x = 2^K17486_G, y = 2^K17486_T, col = 2^K17486_Exp), data = dmdA, size = 3, alpha = 0.9) +
  labs(x= "Metagenomics sequence abundance", y = "Metatranscriptomics transcript abundance", col = "Gene Expression") +
  scale_color_gradient2(low="indianred", mid="navy", high = "gold")

ggplot() + geom_point(aes(x = 2^K17486_G, y = 2^K17486_T, col = 2^K17486_Exp, shape = Layer), data = dmdA_Envi, size = 3, alpha = 0.9) +
  labs(x= "Metagenomics sequence abundance", y = "Metatranscriptomics transcript abundance", col = "Gene Expression") +
  scale_shape_manual(values=c(15,16,17,18))+
  scale_color_gradient2(low="indianred", mid="navy", high = "gold") +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),  legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))

cor(dmdA_ml$K17486_G,dmdA_ml$K17486_T)

######## Maps of abund data #########

# Plot empty map
world <- map_data('world')
base <- ggplot(world, aes(long,lat)) + geom_map(map=world,aes(map_id=region), fill="grey", color="grey") +
  coord_quickmap() +
  theme_void()

map_G_srf = base + geom_point(data=dmdA_Envi[which(dmdA_Envi$Layer=="SRF"),], aes(x=Longitude, y=Latitude, size=2^K17486_G, col=2^K17486_Exp), alpha = 0.9) +
  labs(x = "Longitude",y = "Latitude", size = "Norm. sequence abund.", col = "Norm. gene expression") +
  scale_color_gradient2(low="indianred", mid="navy", high = "gold")
map_G_srf

map_T_srf = base + geom_point(data=dmdA_Envi[which(dmdA_Envi$Layer=="SRF"),], aes(x=Longitude, y=Latitude, size=2^K17486_T, col=2^K17486_Exp), alpha = 0.9) +
  labs(x = "Longitude",y = "Latitude", size = "Norm. transcript abund.", col = "Norm. gene expression") +
  scale_color_gradient2(low="indianred", mid="navy", high = "gold")
map_T_srf

map_Exp_srf = base + geom_point(data=dmdA_Envi[which(dmdA_Envi$Layer=="SRF"),], aes(x=Longitude, y=Latitude, size=2^K17486_Exp, col=2^K17486_Exp), alpha = 0.9) +
  labs(x = "Longitude",y = "Latitude", size = "Gene expression", col = "Gene expression") +
  scale_color_gradient2(low="indianred", mid="navy", high = "gold")
map_Exp_srf

map_Exp_dcm = base + geom_point(data=dmdA_Envi[which(dmdA_Envi$Layer=="DCM"),], aes(x=Longitude, y=Latitude, size=2^K17486_Exp, col=2^K17486_Exp), alpha = 0.9) +
  labs(x = "Longitude",y = "Latitude", size = "Gene expression", col = "Gene expression") +
  scale_color_gradient2(low="indianred", mid="navy", high = "gold", breaks = c(0.25,0.75,1.25)) +
  scale_size_continuous(breaks=c(0.25,0.75,1.25))
  
map_Exp_dcm

multiplot(map_Exp_srf,map_Exp_dcm,cols=1)

######## Multivariate analysis ############

# Now we need the pre-processed variables and to get rid of variables that are out of our scope of interest:
summary(dmdA_ml)
rownames(dmdA_ml) = dmdA_ml[,1]
dmdA_ml = dmdA_ml[,-c(1,5,9)]

rda = rda(dmdA_ml[,c(1:3)]~., data = dmdA_ml[,-c(1:3)])
summary(rda)
RsquareAdj(rda)
anova.cca(rda) # Significant

# Model selection
set.seed(2020)
rda.step.both = ordistep(rda(dmdA_ml[,c(1:3)]~1, data = dmdA_ml[,-c(1:3)]), scope = formula(rda), direction = "both", pstep = 1000)
formula(rda.step.both)

res.rda = rda(formula = formula(rda.step.both), data = dmdA_ml[,-c(1:3)])
summary(res.rda)
RsquareAdj(res.rda) #72.1%
anova.cca(res.rda, by = "axis") 

# Triplot

# Sites in background
station <- scores(res.rda, scaling=2, choices = c(1:2))$sites
station = as.data.frame(station)
station = merge(station, Envi_ml, by.x = "row.names", by.y = "PANGAEA.sample.id")
rownames(station) = station[,1]
station = station[,-1]

# Abund scores in the background as well
ccscore <- scores(res.rda, scaling=2, choices = c(1:2))$species
ccscore = as.data.frame(ccscore)
rownames(ccscore) = c("MetaG", "MetaT", "Expr.")

g1 <- ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  geom_point(aes(station[,1], station[,2]), size = 3, alpha = 0.3) +
  geom_segment(aes(xend = RDA1, yend = RDA2), x=0, y=0, color = "indianred", data = ccscore, arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(x=RDA1*1.1,y=RDA2*1.1, label=rownames(ccscore)), data = ccscore, color = "indianred", fontface = 2) +
  xlab("RDA1 (66.49%)") +  ylab("RDA2 (7.34%)")

# Adding envi variables 
envi.scores <- scores(res.rda, choices = 1:2, display="bp", scaling=2)
row.names(envi.scores) = c("Depth", "Oxygen", "Temperature", "Max_O2_Depth", "HCO3", "Density", "Chlorophyll_A", "Latitude")

g12 <- g1 + geom_segment(data=as.data.frame(envi.scores),aes(xend = RDA1, yend = RDA2),x=0,y=0,size = 0.5, linetype="F1",color = 'navy',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(envi.scores[,1]*1.1, envi.scores[,2]*1.1, label=rownames(envi.scores)), color="navy", fontface = 2)
g12

######## Elastic net regressions ############

# Now we need the pre-processed variables, we build one dataset per y variable
summary(dmdA_ml)
dmdA_G = dmdA_ml[,-c(2,3)]
dmdA_T = dmdA_ml[,-c(1,3)]
dmdA_Exp = dmdA_ml[,-c(1,2)]

# MetaG

set.seed(2020)
training.samples <- dmdA_G$K17486_G %>%
  createDataPartition(p = 0.8, list = FALSE)

train.data  <- dmdA_G[training.samples, ]
test.data <- dmdA_G[-training.samples, ]

# Predictor variables
x <- model.matrix(K17486_G~., train.data)[,-1]
# Outcome variable
y <- train.data$K17486_G

# Build the model using the training set

control = trainControl(method='repeatedcv', 
                       number = 10, 
                       repeats = 3, 
                       returnResamp = "all", 
                       savePredictions = "all")
set.seed(2020)
model <- train(
  K17486_G ~., 
  data = train.data, 
  method = "glmnet",
  trControl = control,
  tuneLength = 10
)

# When increasing the tunelength, an error message due to no variance in predictions for some lambda/alpha couples appears, eg for tunelength of 20:
model$results[c(250:375),] # Couple 0.8578947 0.5208953925 is problematic
model$pred[which(round(model$pred$lambda, digits=9) == 0.520895392 & round(model$pred$alpha, digits=6) == 0.857895),]
# At each fold, prediction is constant --> Rsquared can't be computed. Prediction is probably constant due to the tree not finding a good split
# using this specific alpha/lambda couple. The rest of models ran fine, and increasing tuning length do not lead to better performing models (RMSE even
# slightly increase in predictions), so we keep 10.

# Best tuning parameter
model$bestTune
coef(model$finalModel, model$bestTune$lambda)
plot(coef(model$finalModel, model$bestTune$lambda))

# Make predictions on the test data
predictions <- model %>% predict(test.data[,-1])

# Model performance metrics
data.frame(
  RMSE = RMSE(predictions, test.data$K17486_G),
  Rsquare = R2(predictions, test.data$K17486_G)
)
# Rsquare of predictions 0.69

ggplot() + geom_point(aes(x=test.data$K17486_G, y = predictions)) +
  geom_abline(intercept=0,slope=1,color="indianred") +
  geom_smooth(aes(x=test.data$K17486_G, y = predictions), method = "lm") +
  labs(x="Metagenomic abundance observations", y="Metagenomic abundance predictions")

# MetaT

set.seed(2020)
training.samples <- dmdA_T$K17486_T %>%
  createDataPartition(p = 0.8, list = FALSE)

train.data  <- dmdA_T[training.samples, ]
test.data <- dmdA_T[-training.samples, ]

# Predictor variables
x <- model.matrix(K17486_T~., train.data)[,-1]
# Outcome variable
y <- train.data$K17486_T

# Build the model using the training set
control = trainControl(method='repeatedcv', 
                       number = 10, 
                       repeats = 3, 
                       returnResamp = "all", 
                       savePredictions = "all")

set.seed(2020)
model <- train(
  K17486_T ~., 
  data = train.data, 
  method = "glmnet",
  trControl = control,
  tuneLength = 10
)
# Best tuning parameter
model$bestTune
coef(model$finalModel, model$bestTune$lambda)
plot(coef(model$finalModel, model$bestTune$lambda))

# Make predictions on the test data
predictions <- model %>% predict(test.data[,-1])

# Model performance metrics
data.frame(
  RMSE = RMSE(predictions, test.data$K17486_T),
  Rsquare = R2(predictions, test.data$K17486_T)
)
# Rsquare of predictions 0.68

ggplot() + geom_point(aes(x=test.data$K17486_T, y = predictions)) +
  geom_abline(intercept=0,slope=1,color="indianred") +
  geom_smooth(aes(x=test.data$K17486_T, y = predictions), method = "lm") +
  labs(x="Metatranscriptomics abundance observations", y="Metatranscriptomics abundance predictions")


# Expression

set.seed(2020)
training.samples <- dmdA_Exp$K17486_Exp %>%
  createDataPartition(p = 0.8, list = FALSE)

train.data  <- dmdA_Exp[training.samples, ]
test.data <- dmdA_Exp[-training.samples, ]

# Predictor variables
x <- model.matrix(K17486_Exp~., train.data)[,-1]
# Outcome variable
y <- train.data$K17486_Exp

# Build the model using the training set
set.seed(2020)
control = trainControl(method='repeatedcv', 
                       number = 10, 
                       repeats = 3, 
                       returnResamp = "all", 
                       savePredictions = "all")

model <- train(
  K17486_Exp ~., 
  data = train.data, 
  method = "glmnet",
  trControl = control,
  tuneLength = 10
)
# Best tuning parameter
model$bestTune
coef(model$finalModel, model$bestTune$lambda)
sort(coef(model$finalModel, model$bestTune$lambda))

# Make predictions on the test data
predictions <- model %>% predict(test.data[,-1])

# Model performance metrics
data.frame(
  RMSE = RMSE(predictions, test.data$K17486_Exp),
  Rsquare = R2(predictions, test.data$K17486_Exp)
)
# Rsquare of predictions 0.74.4%

ggplot() + geom_point(aes(x=test.data$K17486_Exp, y = predictions)) +
  geom_abline(intercept=0,slope=1,color="indianred") +
  geom_smooth(aes(x=test.data$K17486_Exp, y = predictions), method = "lm") +
  labs(x="Gene expression observations", y="Gene expression predictions")

