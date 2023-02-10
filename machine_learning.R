# Load packages
library(tidyverse)
library(caret)
library(party)

# Get data
params <- readr::read_csv("data/params.csv")
dataset <- params



# Clean the data and check for duplicates or missing values
dataset <- params[complete.cases(params),]
#dataset$category.rna.rpf<-as.factor(dataset$category.rna.rpf)
sapply(dataset, function(x) table(is.na(x)))
table(duplicated(dataset))
#dat <- dataset[!duplicated(dataset),]
dat <- dataset

# Visual inspection / descriptive statistics
simp <- dat %>% dplyr::select(GENEID,
                              TXID,
                              TXNAME,
                              TXSTRAND,
                              TXSTART,
                              TXEND,
                              nexon,
                              tx_len,
                              cds_len,
                              utr5_len,
                              utr3_len,
                              orf.length,
                              pqs.length,
                              codon.length,
                              enc,
                              scuo,
                              log2.delta_te,
                              delta_te,
                              level.te) # remove kMers and codon data to simplify analysis

simp <- dat

# Correlation matrix of all predictor variables
cormat <- cor(simp %>% dplyr::select(-c(log2.delta_te, delta_te, level.te)) %>% keep(is.numeric))

cormat %>% as.data.frame %>% mutate(var2=rownames(.)) %>%
  tidyr::pivot_longer(!var2, values_to = "value") %>%
  ggplot(aes(x=name,y=var2,fill=abs(value),label=round(value,2))) +
  geom_tile() + geom_label() + xlab("") + ylab("") +
  ggtitle("Correlation matrix of our numeric predictors") +
  labs(fill="Correlation\n(absolute):")

highcorr <- which(cormat > 0.8, arr.ind = T) # Codon length and 5' UTR length are highly correlated
paste(rownames(cormat)[row(cormat)[highcorr]],
      colnames(cormat)[col(cormat)[highcorr]], sep = " vs. ") %>%
  cbind(cormat[highcorr]) # print a list of all correlations

# Boxplots of the associations between continuous predictors and categorical outcome
simp %>% dplyr::select(-c(GENEID,
                   TXID,
                   TXNAME,
                   TXSTRAND,
                   log2.delta_te,
                   delta_te,
                   codon.length)) %>%
  tidyr::pivot_longer(!level.te, values_to = "value") %>%
  ggplot(aes(x=factor(level.te), y=value, fill=factor(level.te))) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(size = 0.7, width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("steelblue","orangered1")) +
  labs(fill="TE:") +
  theme_minimal() +
  facet_wrap(~name, scales="free")

## CUT OVER FROM CHAID ANALYSIS
# Performing dimensionality-reduction with PCA prior to constructing an LDA model
feats <- simple %>% 
  dplyr::select(-c(GENEID, TXID, TXNAME, TXSTRAND, TXSTART, TXEND,
                   log2.delta_te,
                   delta_te,
                   level.te,
                   log2.te_Normoxic,
                   log2.te_Hypoxic,
                   te_Normoxic,
                   te_Hypoxic)) %>%
  keep(is.numeric)
simp.pr <- prcomp(feats, center = TRUE, scale = TRUE)
summary(simp.pr)

screeplot(simp.pr, type = "l", npcs = 10, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(simp.pr$sdev^2 / sum(simp.pr$sdev^2))
plot(cumpro[0:9], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.97036, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)

plot(simp.pr$x[,1],simp.pr$x[,2], xlab="PC1 (X%)", ylab = "PC2 (Y%)", main = "PC1 / PC2 - plot")

library(factoextra)
fviz_pca_ind(simp.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = simple$level.te, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "TE") +
  ggtitle("2D PCA-plot from feature dataset") +
  theme(plot.title = element_text(hjust = 0.5))

simp.pcst <- simp.pr$x[,1:6]
row.names(simp.pcst) <- simp$GENEID
simp.pcst <- bind_cols(simp.pcst, simp$level.te)
colnames(simp.pcst)[7] <- "level.te"

smp_size <- floor(0.75 * nrow(simp.pcst))
train_ind <- sample(nrow(simp.pcst), size = smp_size)
train.df <- as.data.frame(simp.pcst[train_ind, ])
test.df <- as.data.frame(simp.pcst[-train_ind, ])

simp.lda <- MASS::lda(level.te ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = train.df)
simp.lda.predict <- predict(simp.lda, newdata = test.df)
simp.lda.predict.posteriors <- as.data.frame(simp.lda.predict$posterior) # get the posteriors as a dataframe

library(ROCR)
pred <- ROCR::prediction(simp.lda.predict.posteriors[,2], test.df$level.te) # evaluate the model
roc.perf = ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
auc.train <- ROCR::performance(pred, measure = "auc")
auc.train <- auc.train@y.values

plot(roc.perf)
abline(a=0, b= 1)
text(x = .25, y = .65 ,paste("AUC = ", round(auc.train[[1]],3), sep = ""))

# Partition data into training and test datasets
set.seed(2022)
trainRowNumbers <- createDataPartition(dataset$log2.te_Normoxic, p=0.8, list=FALSE) # get row numbers for the training data
trainData <- simple[trainRowNumbers,] # create the training  dataset
testData <- simple[-trainRowNumbers,] # create the test dataset
x_train <- trainData %>% 
  dplyr::select(-c(GENEID, TXID, TXSTRAND, TXNAME, TXSTART, TXEND,
                   log2.delta_te,
                   delta_te,
                   level.te,
                   log2.te_Normoxic,
                   log2.te_Hypoxic,
                   te_Normoxic,
                   te_Hypoxic)) %>% 
  keep(is.numeric) # store X and Y for later use
y_train <- trainData$log2.te_Normoxic

x_test <- testData %>% 
  dplyr::select(-c(GENEID, TXID, TXSTRAND, TXNAME, TXSTART, TXEND,
                   log2.delta_te,
                   delta_te,
                   level.te,
                   log2.te_Normoxic,
                   log2.te_Hypoxic,
                   te_Normoxic,
                   te_Hypoxic)) %>% 
  keep(is.numeric)
y_test <- testData$log2.te_Normoxic


rmse <- function(y, x){
  sqrt(mean((y - x)^2))
}

rsq <- function(y, x) { 
  1 - sum((y - x) ^ 2) / sum((y - mean(y)) ^ 2) 
}

# Model training

library(caret)
# train the model
model <- caret::train(x_train, y_train, method = "lm",
                      trControl = trainControl(method="cv", number=5, verboseIter = T))


set.seed(2022)
mod <- caret::train(x_train, y_train, method = "rf",
                    trControl = trainControl(method="cv", number=5, verboseIter = T))
plot(varImp(mod), main="Feature importance of random forest model on training data")

set.seed(2022)
mod2 <- caret::train(x_train, y_train, method = "avNNet",
                     preProcess = c("center", "scale", "nzv"),
                     tuneGrid = expand.grid(size = seq(3,21,by=3), decay=c(1e-03, 0.01, 0.1, 0), bag = c(T,F)),
                     trControl= trainControl(method = "cv", number = 5, verboseIter = T),
                     importance = T)
plot(varImp(mod2), main = "Feature importance of neural network classifier on training data")


set.seed(2022)
mod3 <- caret::train(x_train, y_train, method = "xgbTree",
                     tuneGrid = expand.grid(nrounds = c(50,100), max_depth = c(5,7,9),
                                            colsample_bytree = c(0.8,1), subsample = c(0.8,1),
                                            min_child_weight = c(1,5,10), eta = c(0.1, 0.3), gamma = c(0, 0.5)),
                     trControl = trainControl(method="cv", number = 5, verboseIter = T))
plot(varImp(mod3), main = "Feature importance of XGBoost model on training data")

modelSvm <- caret::train(x_train, y_train, method="svmRadial",
                         trControl= trainControl(method="cv", number = 5, verboseIter = T))

mod3 <- caret::train(x_train, y_train, method = "xgbTree",
                     tuneGrid = expand.grid(nrounds = c(50,100), max_depth = c(5,7,9),
                                            colsample_bytree = c(0.8,1), subsample = c(0.8,1),
                                            min_child_weight = c(1,5,10), eta = c(0.1, 0.3), gamma = c(0, 0.5)),
                     trControl = trainControl(method="cv", number = 5, verboseIter = T))
plot(varImp(mod3), main = "Feature importance of XGBoost model on training data")

# Compare the performance of the three algorithms
results <- data.frame(Model = c(mod$method, mod2$method, mod3$method),
                      Accuracy = c(max(mod$results$Accuracy), max(mod2$results$Accuracy), max(mod3$results$Accuracy)))
results %>% ggplot(aes(x=Model, y=Accuracy, label=paste(round(100*Accuracy,1),"%"))) +
  geom_col(fill="steelblue") + theme_minimal() + geom_label() +
  ggtitle("Accuracy in the training data by algorithm")

### Step 7: Model evaluation against the test data
predictions <- predict(mod2, newdata = x_test)
confusionMatrix(predictions, y_test)

precision(predictions, y_test)
recall(predictions, y_test)
F_meas(predictions, y_test)

### Step 8: Model deployment
newpatient <- data.frame(age=62, sex=1, cp=0, trestbps=130, chol=220, fbs=0,
                         restecg=0, thalach=161, exang=0, oldpeak=0, slope=0,
                         ca=0, thal=2)

preprocess_new_data <- function(df){
  
  # Convert features to int like the original dataset
  df[,names(df) != "oldpeak"] <- purrr::map_df(df[,names(df) != "oldpeak"], as.integer)
  
  df <- df %>% mutate(restecg = recode(restecg, `1` = 0L),
                      thal = recode(thal, `1` = 0L),
                      ca = recode(ca, `4` = 0L))
  
  # Nominal variables
  # NOTE: We do not have all the values for the dummies in the new dataset
  existing_cols <- names(x_train)[names(x_train) %in% names(df)]
  new_cols <- names(x_train)[!names(x_train) %in% names(df)]
  df[new_cols] <- 0
  nomvars <- c("cp","ca","thal","restecg","slope")
  
  for (i in 1:nrow(df)){
    for (j in length(nomvars)){
      df[i,paste0(nomvars[j],df[nomvars[j]][i])] <- 1
    }
  }
  
  df <- df[,names(df) %in% c(existing_cols, new_cols)]
  
  df$hr_age <- df$thalach / df$age
  df$chol_age <- df$chol / df$age
  df$st <- ifelse(df$oldpeak>0,1,0)
  
  return(df)
}

save(mod2, x_train, preprocess_new_data, file="Heart_disease_prediction.RData")

predict(mod2, newdata = preprocess_new_data(newpatient))
predict(mod2, newdata = preprocess_new_data(newpatient), type = "prob")
