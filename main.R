#Rachel Payne 
#BF591 Assignment 3 

library(tidyverse)
library(RColorBrewer)
library(gplots)
#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, sep) {
  data <- read.csv(intensity_data, sep = " ") %>%
  return()
}


#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  #find variance of each pca -- variance is equal to the square of standard deviation
  pca_var <- pca_results$sdev**2
  # variance explained -divide the variance of each PC by the sum total of the variance of all PCs
  explained_var <- pca_var/sum(pca_var) %>%
  return()
}


#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PCA names, variance explained by each PC, and the
#' cumulative sum of variance explained
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#'
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained and the cumulative variance explained
#' @export
#'
#' @examples
make_variance_tibble <- function(pca_ve, pca_results) {
  pca_var_tib <- tibble(
    PC = factor(str_c("PC", 1:35), str_c("PC", 1:35)),
    variance_explained = pca_ve,
    cumulative = cumsum(pca_ve)) %>%
  return()
}


#' Define a function that creates a bar plot of the variance explained by each
#' PC along with a scatter plot showing the cumulative sum of variance explained
#' using ggplot2
#' @param variance_tibble (tibble): the tibble gnerated in the previous function
#' that contains each PC label, the variance explained by each PC, and the
#' cumulative sum of variance explained
#'
#' @return A ggplot with a barchart representing individual variance
#'   explained and a scatterplot (connected with a line) that represents the
#'   cumulative sum of PCs
#' @export
#'
#' @examples

plot_pca_variance <- function(variance_tibble) {
  plot <- ggplot(variance_tibble) +
    geom_bar(aes(x=PC, y=variance_explained, fill = "Variance Explained"),
             size = 0.3, 
             stat = "identity", color = "black") +
    geom_point(aes(x=PC, y=cumulative, 
                   group=1, color = "Cumulative")) +
    geom_line(aes(x=PC, y=cumulative,
                  group=1, color = "Cumulative")) +
    scale_color_manual(name = "Cumulative",
                       values = "black") +
    scale_fill_manual(name = "Variance Explained",
                      values = "skyblue") +
      theme_classic() +
      theme(axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5), plot.title = element_text(hjust = 0.5)) +
      labs( x = "PC", y = "% variance", title = "Proportion of Variance Explained", title.position = 0.5)

    ## need to add caption to each graph 
    
    return(plot)
}



#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples

make_biplot <- function(metadata, pca_results) {
  #load metadata file
  read_metadata <- read.csv(metadata, sep = ",")
  #make pca_results$x a df, add geo_accession column
  pca_result_df <- as.data.frame(pca_results$x) %>%
    rownames_to_column("geo_accession")
  #join df together -- pca, geo_accession and sixsubtypeclass
  df_joined <- left_join(pca_result_df, select(read_metadata, SixSubtypesClassification, geo_accession), by = "geo_accession")

  ggplot(data=df_joined, mapping = aes(x = PC1, y=PC2, color=SixSubtypesClassification)) +
    theme_classic()+
    geom_point() %>%
    return()
}



#' #' Define a function to return a list of probeids filtered by signifiance
#' #'
#' #' @param diff_exp_csv (str): The path to the differential expression results
#' #'   file we have provided
#' #' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#' #'   value of .01. This is the column "padj" in the CSV.
#' #'
#' #' @return A list with the names of the probeids passing the fdr_threshold
#' #' @export
#' #'
#' #' @examples
#' 

list_significant_probes <- function(diff_exp_csv, fdr_threshold) {
  filt_DE <- read.delim(diff_exp_csv, header = T, sep = ",") %>%
    #filter for probeids <0.01 padj 
    dplyr::filter(padj < fdr_threshold)
  return(rownames(filt_DE))
  }
  


#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  intensity %>%
    #filter probeids present in significant ids list 
    dplyr::filter(rownames(intensity) %in% sig_ids_list) %>%
    as.matrix() %>%
    return()
}

#' #' Define a function that takes the intensity values for significant probes and
#' #' creates a color-blind friendly heatmap
#' #'
#' #' @param de_intensity (matrix): The matrix of intensity values for significant
#' #'   differentially expressed probes returned in part 7
#' #' @param num_colors (int): The number of colors in a specificed RColorBrewer
#' #'   palette
#' #' @param palette (str): The name of the chosen RColorBrewer palette
#' #'
#' #' @return A heatmap displaying the intensity values for the differentially
#' #'   expressed probes
#' #' @export
#' #'
#' #' @examples
plot_heatmap <- function(de_intensity, num_colors, palette){
  ht <- heatmap.2(de_intensity, scale="none",
                       col = brewer.pal(num_colors, palette),
                       ylab="probe ids", trace = "none", density.info = "none",
                       main=" Normalized Intensity of DE Probes")
  return(ht)
}


