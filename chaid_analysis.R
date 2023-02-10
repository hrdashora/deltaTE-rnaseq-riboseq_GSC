# Load libraries
require(rsample) # for dataset and splitting also loads broom and tidyr
require(dplyr)
require(ggplot2)
library(partykit)
require(CHAID)
require(purrr)
require(caret)

# Load dataset
dataset <- readr::read_csv("data/params.csv")
# Remove rows with missing values
dataset <- dataset[complete.cases(dataset),]
# Remove kMers and codon data to simplify the matrix of features
simple <- dataset %>% dplyr::select(GENEID,
                                    TXID,
                                    #TXNAME,
                                    #TXSTRAND,
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
                                    level.te,
                                    log2.te_Normoxic,
                                    log2.te_Hypoxic,
                                    te_Normoxic,
                                    te_Hypoxic) %>%
  mutate_at(vars(level.te), as.ordered)

# Correlation matrix of all predictor variables
cormat <- cor(simple %>% dplyr::select(-c(log2.delta_te,
                                          delta_te,
                                          level.te,
                                          log2.te_Normoxic,
                                          log2.te_Hypoxic,
                                          te_Normoxic,
                                          te_Hypoxic)) %>% keep(is.numeric))

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
simple %>% dplyr::select(-c(GENEID,
                            TXID,
                            TXNAME,
                            TXSTRAND,
                            log2.delta_te,
                            delta_te,
                            log2.te_Normoxic,
                            log2.te_Hypoxic,
                            te_Normoxic,
                            te_Hypoxic,
                            codon.length)) %>%
  tidyr::pivot_longer(!level.te, values_to = "value") %>%
  ggplot(aes(x=factor(level.te), y=value, fill=factor(level.te))) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(size = 0.7, width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("steelblue","orangered1")) +
  labs(fill="TE:") +
  theme_minimal() +
  facet_wrap(~name, scales="free")

library(mRMRe)
dd <- mRMR.data(data = data.frame(dataset %>% dplyr::select(-c(GENEID, TXID, TXNAME))))
spearman_mim <- mim(dd, continuous_estimator = "spearman")
pearson_mim <- mim(dd, continuous_estimator = "pearson")
mRMR.classic(data = dd, target_indices = c(1),
             feature_count = 30)


set.seed(2022)
library(mlbench)
library(caret)
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(x_train, y_train, rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))

# Partition data into training and test sets
set.seed(2022)
trainRowNumbers <- createDataPartition(dataset$log2.te_Normoxic, p=0.8, list=FALSE) # get row numbers for the training data
trainData <- dataset[trainRowNumbers,] # create the training  dataset
testData <- dataset[-trainRowNumbers,] # create the test dataset
x_train <- trainData %>% 
  dplyr::select(-c(GENEID, TXID, TXSTRAND, TXNAME, TXSTART, TXEND,
                   log2.delta_te,
                   delta_te,
                   level.te,
                   log2.te_Normoxic,
                   log2.te_Hypoxic,
                   te_Normoxic,
                   te_Hypoxic,
                   symbol)) %>% 
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
                   te_Hypoxic,
                   symbol)) %>% 
  keep(is.numeric)
y_test <- testData$log2.te_Normoxic

# train the model
model <- caret::train(x_train, y_train, method = "lm",
                      trControl = trainControl(method="cv", number=5, verboseIter = T))
