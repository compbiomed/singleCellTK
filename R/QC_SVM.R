# to implement the SVM for the quality control 
# (Ref: Classification of low quality cells from single-cell RNA-seq data
#       https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0888-1#Sec12  )

# ToDo: need to (1)select parameters for SVM, need to think about whether let users choose the parameters
#       ,or(2) use cross-validation to choose the parameters for the users, ----ing 
#        
#         (3) show the test error rate if test_label is included ?     


#' To implement SVM for quality control 
#' 
#' @param train_set The training dataset (cells x features) 
#' @param train_label   The training label annotation specifying the high/low quality of the training dataset.
#' @param test_set The test dataset (cells x features)
#' @param var feature names used for predicting cell quality ## 
#' @param tune_para boolean values determining if doing the parameter tuning 
#' @param gamma_range vector of gamma paramters provided if tune_para is TRUE
#' @param cost_range vecor of cost patameters provided if tune_para is TRUE
#' @export
singlecell_SVM <- function(train_set, train_label, test_set, var,
                           tune_para=FALSE, gamma_range, cost_range ){
  # train data frame preparation
  ## should examine train_set is data.frame
  train_df <- train_set[, colnames(train_set) %in% var]
  train_df <- data.frame(train_df, l=as.factor(train_label))

  # test data frame preparation
  test_df <- test_set[, colnames(test_set) %in% var]

 # fit SVM model using the train dataset
  sfit <- e1071::svm(l~., data=train_df,
                     scale= FALSE,
                     kernel="radial",   # shoudl also make this changable but set the "radial" as default
                     gamma = 1,  # gamma should not be fixed later
                     cost = 0.1 # cost should not be fixed later
                     # also consider including weight later
                     # also consider c-v here
                     )
  test.pred = predict(sfit, test_set)

  # parameter tuning for the SVM model
    # if the set_tune = T, then also provide a data.frame of the parameters(para =),
    # and the result would be an object, which gives the parameters for the best model
    # as well as the prediction for the test dataset
  sfit.tune <- NULL
  if(tune_para ==TRUE){
    if(!is.null(gamma_range) & !is.null(cost_range)){  # both gamma & cost are given
      sfit.tune <- e1071::tune.svm(l~., data = train_df,
                                   sampling="fix",
                                   gamma = gamma_range,   # gamma should be given by users
                                   cost = cost_range      # cost should be given by users
                                   )
    }else if(!is.null(gamma_range) & is.null(cost_range)){ # only gamma is given
      sfit.tune <- e1071::tune.svm(l~., data = train_df,
                                   sampling="fix",
                                   gamma = gamma_range   # gamma should be given by users
                                   )
    }else if(is.null(gamma_range) & !is.null(cost_range)){
      sfit.tune <- e1071::tune.svm(l~., data = train_df,
                                  sampling="fix",
                                  cost = cost_range   # cost should be given by users
                                  )
    }else{ # neither gamma or cost is given
      sfit.tune <- e1071::tune.svm(l~., data = train_df,
                                   sampling="fix")
    }
  }  # should has output for this later

 #if(!is.null(sfit.tune)){
   #cat("which part of the sfit.tune should be returned to the user")
   #}
  return(test.pred)
}
