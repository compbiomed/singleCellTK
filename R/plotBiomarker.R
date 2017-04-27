
#' get Biomarker
#'   
#' Use this function to get expression or binary data of gene list
#'
#' @param count_data A SCE object
#' @param gene gene list
#' @param binary binary or gradient
#' @param visual visualization type of plot (PCA or tSNE)
#' @param shape visualization shape
#' @param axis_df df of PC or tSNE components
#' @param x x coordinate for PCA
#' @param y y coordinate for PCA
#' 
#' @return A Biomarker plot
#' @export plotBiomarker
#'

plotBiomarker <- function(count_data, gene, binary="Binary", visual="PCA",shape=NULL,axis_df=NULL,x="PC1", y="PC2"){
  if(shape == "No Shape"){
    shape <- NULL
  }
  if(visual=="PCA"){
    if(is.null(axis_df)){
      axis_df <- getPCA(count_data)
    }
    variances <- attr(axis_df,"percentVar")
  }
  if(visual=="tSNE"){
    if(is.null(axis_df)){
      axis_df <- getTSNE(count_data)
    }
  }
  
  bio_df <- getBiomarker(count_data,gene,binary)
  l <- axis_df
  if(!is.null(shape)){
    l$shape <- factor(eval(parse(text = paste("pData(count_data)$",shape,sep=""))))
  }
  gene_name <- colnames(bio_df)[2]
  colnames(bio_df)[2] <- "expression"
  l$Sample <- bio_df$sample
  l$expression <- bio_df$expression
  c <- counts(count_data)[c(gene_name),]
  percent <- round(100*sum(c>0)/length(c),2)
  if(visual == "PCA"){
    if(binary == "Binary"){
      g <- ggplot(l, aes_string(x, y, label="Sample")) + aes(color=ifelse(expression=="TRUE","blue","grey"))+
        geom_point() +
        scale_color_manual(labels = c("Yes", "No"), values=c("Blue", "Grey")) + 
        labs(color = "Expression")
    }
    else if (binary == "Continuous"){
      if(min(round(l$expression, 6)) == max(round(l$expression, 6))){
        g <- ggplot(l, aes_string(x, y, label="Sample")) +
          geom_point(color="grey") 
      } else{
        g <- ggplot(l, aes_string(x, y,label="Sample", color="expression"))+
          scale_colour_gradient(limits=c(min(l$expression), max(l$expression)), low="grey", high="blue")+
          geom_point()
      }
      g <- g + labs(color = "Expression")
    }
    g <- g+
      ggtitle(paste(gene_name," - ",percent,"%"," cells",sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))+
      labs(x = paste(x,toString(round(variances[strtoi(strsplit(x,"PC")[[1]][-1])]*100,2)),"%"),y = paste(y,toString(round(variances[strtoi(strsplit(y,"PC")[[1]][-1])]*100,2)),"%"))
  } else if(visual == "tSNE"){
    if(binary == "Binary"){
      g <- ggplot(l, aes(X1, X2, label=Sample, color=ifelse(expression=="TRUE","blue","grey")))+
        geom_point() +
        scale_color_manual(labels = c("Yes", "No"), values=c("Blue", "Grey")) + 
        labs(color = "Expression")
    }
    else if (binary == "Continuous"){
      if(min(round(l$expression, 6)) == max(round(l$expression, 6))) {
        g <- ggplot(l, aes(X1, X2, label=Sample)) +
          geom_point(color="grey")
      } else{
      g <- ggplot(l, aes(X1, X2, label=Sample, color=expression))+
        scale_colour_gradient(limits=c(min(l$expression), max(l$expression)), low="grey", high="blue") +
        geom_point()
      }
      g <- g + labs(color = "Expression")
    }
    g <- g +
      ggtitle(paste(gene_name," - ",percent,"%"," cells",sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  if(!is.null(shape)){
    g <- g + aes_string(shape="shape") +
         labs(shape = shape)
  }
  return(g)
}