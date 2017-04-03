# to implement the SVM for the quality control 
# (Ref: Classification of low quality cells from single-cell RNA-seq data
#       https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0888-1#Sec12  )

# ToDo: need to (1)select parameters for SVM, need to think about whether let users choose the parameters
#       ,or(2) use cross-validation to choose the parameters for the users,
#        
#         (3) show the test error rate if test_label is included ?


 #' To implement SVM for quality control 
 #' 
 #' @param train_set The training dataset (cells x features) 
 #' @param train_label   The training label annotation specifying the high/low quality of the training dataset.
 #' @param test_set The test dataset (cells x features)
 #' @param var feature names used for predicting cell quality ## 
 #' @example  ToDo
 #' @export
 singlecell_SVM <- function(train_set, train_label, test_set, var){
   # train data frame preparation 
        ## should examine train_set is data.frame 
   train_df <- train_set[, colnames(train_set) %in% var]
   train_df <- data.frame(train_df, l=as.factor(train_label))
   
   # test data frame preparation 
   test_df <- test_set[, colnames(test_set) %in% var]
   
   # fit SVM model using the train dataset 
   sfit <- e1071::svm(l~., data=train_df, 
                      scale= FALSE,
                      kernel="radial", 
                      gamma = 1,  # gamma should not be fixed later 
                      cost = 0.1 # cost should not be fixed later 
                      # also consider including weight later 
                      # also consider c-v here
                      )
   
   
   test.pred = predict(sfit, test_set)
   return(test.pred)
   
 } 