## R code for publication: "The use of text-mining and machine learning algorithms in systematic reviews: 
## reducing workload in preclinical biomedical sciences and reducing human screening error" - Alexandra Bannach-Brown 

## load necessary packages
install.packages("caret")
install.packages("e1071")
library(lattice)
library(ggplot2)
library(e1071)
library(caret)
library(tidyverse)

## read in data files
library(readr)

classifier1 <- read_csv("~/classifier1.csv")
classifier2 <- read_csv("~/classifier2.csv")
classifier3 <- read_csv("~/classifier3.csv")
classifier4 <- read_csv("~/classifier4.csv")
classifier5 <- read_csv("~/classifier5.csv")

cutoffs <- read_csv("~/cutoffs.csv", col_names = FALSE)
cutoffs2 <- read_csv("~/cutoffs2.csv", col_names = FALSE)

## datasets for the error analysis
classifier6updated <- read_csv("~/classifier6-newvalidation.csv")

classifier5updated <- classifier5

# applying the corrected validation set results to classifier 5 datasheet
classifier5updated$`Inclusion Status`  <- classifier6updated$`Inclusion Status`


## remove empty colums
classifier1 <- classifier1 %>% na.omit()
classifier2 <- classifier2 %>% na.omit()
classifier3 <- classifier3 %>% na.omit()
classifier4 <- classifier4 %>% na.omit()
classifier5 <- classifier5 %>% na.omit()

classifier5updated <- classifier5updated %>% na.omit()
classifier6updated <- classifier6updated %>% na.omit()

## make inclusion status to factors
classifier1$`Incl(1)/Excl(0)`<- as.factor(classifier1$`Incl(1)/Excl(0)`)
classifier2$`Incl(1)/Excl(0)`<- as.factor(classifier2$`Incl(1)/Excl(0)`)
classifier3$`Incl(1)/Excl(0)`<- as.factor(classifier3$`Incl(1)/Excl(0)`)
classifier4$`Inclusion Status`<- as.factor(classifier4$`Inclusion Status`)
classifier5$`Inclusion Status`<- as.factor(classifier5$`Inclusion Status`)

classifier5updated$`Inclusion Status`<- as.factor(classifier5updated$`Inclusion Status`)
classifier6updated$`Inclusion Status`<- as.factor(classifier6updated$`Inclusion Status`)

## make all other columns levels
#levels(classifier1) <- c(0,1)
#levels(classifier2) <- c(0,1)
#levels(classifier3) <- c(0,1)
#levels(classifier4) <- c(0,1)
#levels(classifier5) <- c(0,1)


##create output matrices for for loops
outputClassifier1 <- data.frame(matrix(0, nrow = nrow(cutoffs), ncol=5))
  colnames(outputClassifier1) <- c("cutoff", "sens", "spec", "precision", "accuracy")
outputClassifier2 <- data.frame(matrix(0, nrow = nrow(cutoffs2), ncol=5))
  colnames(outputClassifier2) <- c("cutoff", "sens", "spec", "precision", "accuracy")
outputClassifier3 <- data.frame(matrix(0, nrow = nrow(cutoffs), ncol=5))
  colnames(outputClassifier3) <- c("i", "sens", "spec", "precision", "accuracy")
outputClassifier4 <- data.frame(matrix(0, nrow = nrow(cutoffs), ncol=5))
  colnames(outputClassifier4) <- c("i", "sens", "spec", "precision", "accuracy")
outputClassifier5 <- data.frame(matrix(0, nrow = nrow(cutoffs), ncol=5))
  colnames(outputClassifier5) <- c("i", "sens", "spec", "precision", "accuracy")
  
  
outputClassifier5updated <- data.frame(matrix(0, nrow = nrow(cutoffs), ncol=5))
  colnames(outputClassifier5updated) <- c("i", "sens", "spec", "precision", "accuracy")
outputClassifier6updated <- data.frame(matrix(0, nrow = nrow(cutoffs), ncol=5))
  colnames(outputClassifier6updated) <- c("i", "sens", "spec", "precision", "accuracy")
  
## for each column in dataset compare that column to the reference data, the correct included/excluded
  
df1 <- function(classifiers, outputClassifier, cutoffs) {
    
    for (i in 4:ncol(classifiers)){ 
      classifier <- classifiers[,i]
      classifier <- factor(unlist(classifier))
      levels(classifier) <- c(0,1)
      ## calculate sensitivity and specificity for each cut off value
      
      class1Conf <- confusionMatrix(data = classifier, reference = classifier1$`Incl(1)/Excl(0)`, positive = "1", dnn = c("Classifier 1", "Gold Standard"))
      
      outputClassifier[i-3,2] <- class1Conf$byClass["Sensitivity"]
      outputClassifier[i-3,3] <- class1Conf$byClass["Specificity"]
      outputClassifier[i-3,4] <- class1Conf$byClass["Precision"]
      outputClassifier[i-3,5] <- class1Conf$overall["Accuracy"]
    }
    
    ## append sens and spec values along with cut off value to a new table
    outputClassifier[,1] <- cutoffs
    
    return(outputClassifier)
  }
  
  class1Conf <- df1(classifier1, outputClassifier1, cutoffs)

  
## classifiers 2 
df2 <- function(classifiers, outputClassifier, cutoffs2) {
    
    for (i in 4:ncol(classifiers)){ 
      classifier <- classifiers[,i]
      classifier <- factor(unlist(classifier))
      levels(classifier) <- c(0,1)
      
      ## calculate sensitivity and specificity for each cut off value
      class2Conf <- confusionMatrix(data = classifier, reference = classifier2$`Incl(1)/Excl(0)`, positive = "1", dnn = c("Classifier 2", "Gold Standard"))
      
      outputClassifier[i-3,2] <- class2Conf$byClass["Sensitivity"]
      outputClassifier[i-3,3] <- class2Conf$byClass["Specificity"]
      outputClassifier[i-3,4] <- class2Conf$byClass["Precision"]
      outputClassifier[i-3,5] <- class2Conf$overall["Accuracy"]
    }
    
    ## append sens and spec values along with cut off value to a new table
    outputClassifier[,1] <- cutoffs2
    
    return(outputClassifier)
  }
  
  class2Conf <- df2(classifier2, outputClassifier2, cutoffs2)  

  
## classifier 3
  df3 <- function(classifiers, outputClassifier, cutoffs) {
    
    for (i in 4:ncol(classifiers)){ 
      classifier <- classifiers[,i]
      classifier <- factor(unlist(classifier))
      levels(classifier) <- c(0,1)
      
      ## calculate sensitivity and specificity for each cut off value
      class3Conf <- confusionMatrix(data = classifier, reference = classifier3$`Incl(1)/Excl(0)`, positive = "1", dnn = c("Classifier 3", "Gold Standard"))
      
      outputClassifier[i-3,2] <- class3Conf$byClass["Sensitivity"]
      outputClassifier[i-3,3] <- class3Conf$byClass["Specificity"]
      outputClassifier[i-3,4] <- class3Conf$byClass["Precision"]
      outputClassifier[i-3,5] <- class3Conf$overall["Accuracy"]
    }
    
    ## append sens and spec values along with cut off value to a new table
    outputClassifier[,1] <- cutoffs
    
    return(outputClassifier)
  }
  
  class3Conf <- df3(classifier3, outputClassifier3, cutoffs)  
  
  
## classifier 4
df4 <- function(classifiers, outputClassifier, cutoffs) {
    
    for (i in 4:ncol(classifiers)){ 
      classifier <- classifiers[,i]
      classifier <- factor(unlist(classifier))
      levels(classifier) <- c(0,1)
      
      ## calculate sensitivity and specificity for each cut off value
      class4Conf <- confusionMatrix(data = classifier, reference = classifier4$`Inclusion Status`, positive = "1", dnn = c("Classifier 4", "Gold Standard"))
      
      outputClassifier[i-3,2] <- class4Conf$byClass["Sensitivity"]
      outputClassifier[i-3,3] <- class4Conf$byClass["Specificity"]
      outputClassifier[i-3,4] <- class4Conf$byClass["Precision"]
      outputClassifier[i-3,5] <- class4Conf$overall["Accuracy"]
    }
    
    ## append sens and spec values along with cut off value to a new table
    outputClassifier[,1] <- cutoffs
    
    return(outputClassifier)
  }
  
  class4Conf <- df4(classifier4, outputClassifier4, cutoffs)

  
##classifier 5
df5 <- function(classifiers, outputClassifier, cutoffs) {
    
    for (i in 4:ncol(classifiers)){ 
      classifier <- classifiers[,i]
      classifier <- factor(unlist(classifier))
      levels(classifier) <- c(0,1)
      
      ## calculate sensitivity and specificity for each cut off value
      class5Conf <- confusionMatrix(data = classifier, reference = classifier5$`Inclusion Status`, positive = "1", dnn = c("Classifier 5", "Gold Standard"))
      
      outputClassifier[i-3,2] <- class5Conf$byClass["Sensitivity"]
      outputClassifier[i-3,3] <- class5Conf$byClass["Specificity"]
      outputClassifier[i-3,4] <- class5Conf$byClass["Precision"]
      outputClassifier[i-3,5] <- class5Conf$overall["Accuracy"]
    }
    
    ## append sens and spec values along with cut off value to a new table
    outputClassifier[,1] <- cutoffs
    
    return(outputClassifier)
  }
  
  class5Conf <- df5(classifier5, outputClassifier5, cutoffs)

  
  
  
##classifier 5 updated
df5up <- function(classifiers, outputClassifier, cutoffs) {
    
    for (i in 4:ncol(classifiers)){ 
      classifier <- classifiers[,i]
      classifier <- factor(unlist(classifier))
      levels(classifier) <- c(0,1)
      
      ## calculate sensitivity and specificity for each cut off value
      class5updatedConf <- confusionMatrix(data = classifier, reference = classifier5updated$`Inclusion Status`, positive = "1", dnn = c("Classifier 5", "Gold Standard"))
      
      outputClassifier[i-3,2] <- class5updatedConf$byClass["Sensitivity"]
      outputClassifier[i-3,3] <- class5updatedConf$byClass["Specificity"]
      outputClassifier[i-3,4] <- class5updatedConf$byClass["Precision"]
      outputClassifier[i-3,5] <- class5updatedConf$overall["Accuracy"]
    }
    
    ## append sens and spec values along with cut off value to a new table
    outputClassifier[,1] <- cutoffs
    
    return(outputClassifier)
  }
  
class5updatedConf <- df5up(classifier5updated, outputClassifier5updated, cutoffs)  
  
  
  
##classifier 6 
df6up <- function(classifiers, outputClassifier, cutoffs) {
    
    for (i in 4:ncol(classifiers)){ 
      classifier <- classifiers[,i]
      classifier <- factor(unlist(classifier))
      levels(classifier) <- c(0,1)
      
      ## calculate sensitivity and specificity for each cut off value
      class6Confupdated <- confusionMatrix(data = classifier, reference = classifier6updated$`Inclusion Status`, positive = "1", dnn = c("Classifier 6", "Gold Standard"))
      
      outputClassifier[i-3,2] <- class6Confupdated$byClass["Sensitivity"]
      outputClassifier[i-3,3] <- class6Confupdated$byClass["Specificity"]
      outputClassifier[i-3,4] <- class6Confupdated$byClass["Precision"]
      outputClassifier[i-3,5] <- class6Confupdated$overall["Accuracy"]
    }
    
    ## append sens and spec values along with cut off value to a new table
    outputClassifier[,1] <- cutoffs
    
    return(outputClassifier)
  }

class6Confupdated <- df6up(classifier6updated, outputClassifier6updated, cutoffs)
  


## write classifier sensitivity and specificity to file
write.csv(class1Conf, file = "class1Conf.csv")
write.csv(class2Conf, file = "class2Conf.csv")
write.csv(class3Conf, file = "class3Conf.csv")
write.csv(class4Conf, file = "class4Conf.csv")
write.csv(class5Conf, file = "class5Conf.csv")
write.csv(class6Conf, file = "class6Conf.csv")

write.csv(class5updatedConf, file = "class5ConfUpdated.csv")
write.csv(class6Confupdated, file = "class6ConfUpdated.csv")


## graph each classifier output on sens/spec by cutoff
install.packages("plotly")
library(plotly)

## prepare data to be graphed
graphData <- list(cutoff1 = class1Conf$cutoff, class1sens = class1Conf$sens, class1spec = class1Conf$spec, 
#classifier 2 ommitted due to different row length # cutoff2 = class2Conf$X1, class2sens = class2Conf$X2, class2spec = class2Conf$X3, 
                  class3sens = class3Conf$sens, class3spec = class3Conf$spec, 
                  class4sens = class4Conf$sens, class4spec = class4Conf$spec, 
                  class5sens =  class5Conf$sens, class5spec = class5Conf$spec, 
          class5updatedsens =  class5updatedConf$sens, class5updatedspec = class5updatedConf$spec, 
          class6updatedsens = class6Confupdated$sens, class6updatedspec = class6Confupdated$spec)

dataGraph<- as.data.frame(graphData)
row.names(dataGraph) <- graphData$cutoff1
row.names(class2Conf) <- class2Conf$cutoff


## build plot 
p <- plot_ly(dataGraph, x = ~class1spec[5:95], y = ~class1sens[5:95], 
name = 'Classifier 1', type = 'scatter', mode = 'lines', 
text = paste("Cut-Off Value: ", row.names(dataGraph[5:95,]), 
             "<br>Sensitivity: ", dataGraph$class1sens[5:95],
             "<br>Specificity: ", dataGraph$class1spec[5:95]), 
hoverinfo='text'
) %>%
  #Classifier 2 data added manuualy here due to different row length
  add_trace(x=~class2Conf$sens[5:285], y = ~class2Conf$spec[5:285], 
            name = 'Classifier 2', mode='lines',
            text = paste("Cut-Off Value: ", row.names(class2Conf[5:285,]), 
                         "<br>Sensitivity: ", class2Conf$sens[5:285],
                         "<br>Specificity: ", class2Conf$spec[5:285]), 
            hoverinfo='text')%>%
add_trace(x=~class3spec[5:95], y = ~class3sens[5:95], 
          name = 'Classifier 3', mode='lines', 
          text = paste("Cut-Off Value: ", row.names(dataGraph[5:95,]), 
 "<br>Sensitivity: ", dataGraph$class3sens[5:95],
 "<br>Specificity: ", dataGraph$class3spec[5:95]), 
          hoverinfo='text') %>%
add_trace(x=~class4spec[5:95], y = ~class4sens[5:95], 
          name = 'Classifier 4', mode='lines',
          text = paste("Cut-Off Value: ", row.names(dataGraph[5:95,]), 
                       "<br>Sensitivity: ", dataGraph$class4sens[5:95],
                       "<br>Specificity: ", dataGraph$class4spec[5:95]), 
          hoverinfo='text')%>%
add_trace(x=~class5spec[2:95], y = ~class5sens[2:95], 
          name = 'Classifier 5', mode='lines',
          text = paste("Cut-Off Value: ", row.names(dataGraph[2:95,]), 
                       "<br>Sensitivity: ", dataGraph$class5sens[2:95],
                       "<br>Specificity: ", dataGraph$class5spec[2:95]), 
          hoverinfo='text')%>%
  layout(title = 'Performance of Machine Learning Classifiers',
         xaxis = list(title = 'Specificity'),
         yaxis = list (title = 'Sensitivity'))

##print plot
p



##plot comparing class. 5 & 6
plot56 <- plot_ly(dataGraph, x=~class5updatedspec[2:100], y = ~class5updatedsens[2:100], 
                  name = 'Classifier 5', mode='lines', type = 'scatter',
                  text = paste("Cut-Off Value: ", row.names(dataGraph[2:100,]), 
                               "<br>Sensitivity: ", dataGraph$class5sens[2:100],
                               "<br>Specificity: ", dataGraph$class5spec[2:100]), 
                  hoverinfo='text'
)%>%
  add_trace(x=~class6updatedspec[2:100], y = ~class6updatedsens[2:100], 
            name = 'Classifier 6', mode='lines',
            text = paste("Cut-Off Value: ", row.names(dataGraph[2:100,]), 
                         "<br>Sensitivity: ", dataGraph$class6updatedsens[2:100],
                         "<br>Specificity: ", dataGraph$class6updatedspec[2:100]), 
            hoverinfo='text')%>%
  layout(title = 'Performance of Machine Learning Classifiers on Corrected Validation Set',
         xaxis = list(title = 'Specificity'),
         yaxis = list (title = 'Sensitivity'))

plot56


## ROC curve & AUC comparison

#install necessary packages

install.packages("pROC")
library(pROC)
library(ggplot2)

rocClass5 <- roc(classifier5updated$`Inclusion Status`, classifier5updated$`0.1`, ci=TRUE, of="auc")
rocClass5
ci.auc(rocClass5)

rocClass6 <- roc(classifier6updated$`Inclusion Status`, classifier6updated$`0.1`, ci=TRUE, of="auc")
ci.auc(rocClass6)

roc.test(rocClass5, rocClass6, reuse.auc = FALSE)

#plotting Roc curves
ROCplot <- ggroc(list(rocClass5, rocClass6))
ROCplot

#sample size and power calcs
powerROCclass5 <- power.roc.test(rocClass5, sig.level = 0.001, power= NULL, method = delong)
powerROCclass5

powerROCclass6 <- power.roc.test(rocClass6, sig.level = 0.001, power= NULL, method = delong)
powerROCclass6
