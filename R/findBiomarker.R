#' Identify biomarkers
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param condition Which condition is the target condition
#' @param useAssay useAssay Indicate which assay to use
#' @param gene_subset gene names that will be used as biomarker pool, default is first 100 genes
#' @param nfolds number of splits in CV
#' @param nrepeats number of CVs with different random splits
#' @param seed for repeatable research
#' @param percent_top_biomarker Top importance percentage to pick biomarker
#' @param model_name one of "logistic regression", "random forest"
#'
#' @return A R list containing biomarker info, importantce score plot, and ROC plot
#'
#' @import caret
#' @import plotROC
#' @import forcats
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import glmnet
#' @import DMwR
#' @importFrom ggplot2 geom_col aes coord_flip theme_bw coord_equal annotate
#'
#' @examples
#' data("mouseBrainSubsetSCE")
#' p <- findBiomarker(mouseBrainSubsetSCE,
#'                     condition = "level1class",
#'                     useAssay = "logcounts",
#'                     gene_subset = seq(100),
#'                     nfolds = 3,
#'                     nrepeats = 3,
#'                     seed = 99,
#'                     percent_top_biomarker = 0.2,
#'                     model_name = "logistic regression")
#' 
#'

#'
#' @export

findBiomarker <- function(inSCE,
                          condition,
                          useAssay,
                          gene_subset = seq(100),
                          nfolds = 3,
                          nrepeats = 3,
                          seed = 99,
                          percent_top_biomarker = 0.2,
                          model_name = c("logistic regression", "random forest"),
                          level_selected_1va = NULL,
                          level_selected_1v1 = NULL) {

    ## SEED
    # bioC not suggesst add set seed function in R code
    # set.seed(seed)

  
    
    # retrieve data from SCE
    cnts <- SummarizedExperiment::assay(inSCE, useAssay)
    annotData <-
    SingleCellExperiment::colData(inSCE)[, c(condition),
                                            drop = FALSE]
    # subset of genes
    if (sum(suppressWarnings(gene_subset == seq(100)))){
      cnts <- cnts[gene_subset,]
    } else{
      cnts <- cnts[which(gene_subset %in% rownames(cnts)),]
    }
    # print(dim(cnts))
    
    # filter out zero variance genes
    gene_index_keep <- c()
    for (i in seq(nrow(cnts))){
      if (sd(cnts[i,]) > 0.1){
        gene_index_keep <- c(gene_index_keep, i)
      }
    }
    cnts <- cnts[gene_index_keep,]
    # print(dim(cnts))
    
    # transpose
    cnts %<>%
      base::t() %>%
      base::as.data.frame()


    # add target variable
    cnts[,'y'] <- annotData[[1]]
  
    # check levels
    y <- as.character(cnts[,'y'])
    if (!is.null(level_selected_1va)){
      y[which(y == level_selected_1va)] <- "positive_1"
      y[which(y != "positive_1")] <- "negative_0"
      cnts[,'y'] <- as.factor(y)
    } else if (!is.null(level_selected_1v1)){
      y[which(y == level_selected_1va[1])] <- "positive_1"
      y[which(y == level_selected_1va[2])] <- "negative_0"
      cnts <- cnts[which(y %in% c("positive_1", "negative_0")),]
      y <- y[which(y %in% c("positive_1", "negative_0"))]
      cnts[,'y'] <- as.factor(y)
    } else {
      y[which(y == y[1])] <- "positive_1"
      y[which(y != "positive_1")] <- "negative_0"
      cnts[,'y'] <- as.factor(y)
    }
    
    # check if level number is enough for modeling
    num_positive <- sum(y == "positive_1")
    num_negative <- sum(y == "negative_0")
    if (num_positive < 6 | num_negative < 6){
      stop("One level of dataset has less than 6 samples, need more data for cross-validation!")
    }
    
    
    
    # remove "-" from gene names
    colnames(cnts) <- gsub("-", "", colnames(cnts))

    # set up classification model prameters
    fitControl <- caret::trainControl(## n1-fold CV
                               method = "repeatedcv",
                               number = nfolds,
                               ## repeated n2 times
                               repeats = nrepeats,
                               classProbs = TRUE,
                               summaryFunction = twoClassSummary,
                               sampling = "smote",
                               savePredictions = TRUE)

    # choose different model
    if (model_name == "logistic regression"){
        model_fit <- caret::train(y ~ .,
                    data = cnts,
                    method = "glmnet",
                    tuneLength = 5,
                    trControl = fitControl,
                    metric = "ROC")
    } else if (model_name == "gbm"){
        model_fit <- caret::train(y ~ .,
                     data = cnts,
                     method = "gbm",
                     trControl = fitControl,
                     tuneLength = 5,
                     metric = "ROC",
                     ## This last option is actually one
                     ## for gbm() that passes through
                     verbose = FALSE)
    } else if (model_name == "random forest"){
        model_fit <- caret::train(y ~ .,
                    data = cnts,
                    method = "ranger",
                    trControl = fitControl,
                    tuneLength = 5,
                    metric = "ROC",
                    # ranger specific parameter
                    importance = "impurity")
    }


    biomarker <- caret::varImp(model_fit)$importance %>%
                      base::as.data.frame() %>%
                      tibble::rownames_to_column() %>%
                      dplyr::rename(importance = Overall) %>%
                      dplyr::rename(biomarker = rowname) %>%
                      dplyr::arrange(importance) %>%
                      dplyr::filter(importance > quantile(importance, 1-percent_top_biomarker)) %>%
                      dplyr::select(biomarker) %>%
                      .$biomarker

    importance_plot <- caret::varImp(model_fit)$importance %>%
                      base::as.data.frame() %>%
                      tibble::rownames_to_column() %>%
                      dplyr::rename(importance = Overall) %>%
                      dplyr::rename(biomarker = rowname) %>%
                      dplyr::arrange(importance) %>%
                      dplyr::filter(importance > quantile(importance, 1-percent_top_biomarker)) %>%
                      dplyr::mutate(biomarker = forcats::fct_inorder(biomarker)) %>%
                      ggplot2::ggplot()+
                        geom_col(aes(x = biomarker, y = importance))+
                        coord_flip()+
                        theme_bw()

    # retrain the model using the biomarker
    cnts <- cnts %>%
                        dplyr::select(biomarker,y)

    # choose different model
    if (model_name == "logistic regression"){
        model_fit <- caret::train(y ~ .,
                    data = cnts,
                    method = "glmnet",
                    tuneLength = 5,
                    trControl = fitControl,
                    metric = "ROC")
    } else if (model_name == "gbm"){
        model_fit <- caret::train(y ~ .,
                     data = cnts,
                     method = "gbm",
                     trControl = fitControl,
                     tuneLength = 5,
                     metric = "ROC",
                     ## This last option is actually one
                     ## for gbm() that passes through
                     verbose = FALSE)
    } else if (model_name == "random forest"){
        model_fit <- caret::train(y ~ .,
                    data = cnts,
                    method = "ranger",
                    trControl = fitControl,
                    tuneLength = 5,
                    metric = "ROC",
                    # ranger specific parameter
                    importance = "impurity")
    }



    # print the biomarker CV performance
    # biomarker_cv_performance <- model_fit$results %>%
    #     dplyr::select(ROC, Sens, Spec) %>%
    #     dplyr::filter(ROC == max(ROC))

    prob_pred <- as.numeric(model_fit$pred$obs)
    prob_pred[prob_pred == 1] <- 0
    prob_pred[prob_pred == 2] <- 1
    df_roc <- data.frame(m = model_fit$pred[,which(colnames(model_fit$pred)
                                                   == levels(model_fit$pred$obs)[2])],
                         d = prob_pred,
                         stringsAsFactors = FALSE)

    g <- ggplot(df_roc, aes(m=m, d=d)) +
    geom_roc(n.cuts=0) +
    coord_equal() +
    style_roc()

    roc_plot <- g + annotate("text", x=0.75, y=0.25,
                             label=paste("AUC =", round((calc_auc(g))$AUC, 4)))

    # output a list
    list_output <- list(biomarker = biomarker,
                        importance_plot = importance_plot,
                        roc_plot = roc_plot)

    return(list_output)

}
