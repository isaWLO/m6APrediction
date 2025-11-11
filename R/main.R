#' @title Encode DNA Strings into a Factor Data Frame
#' @description Converts a character vector of DNA strings into a data frame where each nucleotide becomes a factor column. This is useful for preparing sequence data for machine learning models.
#' @param dna_strings A character vector where each element is a DNA sequence (e.g., "ATCGG").
#' @return A data frame. Each column corresponds to a nucleotide position (named "nt_pos1", "nt_pos2", ...) and contains factors with levels "A", "T", "C", "G".
#' @examples
#' \dontrun{
#' dna_sequences <- c("ATCG", "GCTA")
#' encoded_df <- dna_encoding(dna_sequences)
#' print(encoded_df)
#' }
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}



#' @title Predict m6A Sites for Multiple Sequences
#' @description Predicts m6A modification status for multiple sequences provided in a data frame. It internally encodes DNA sequences and sets correct factor levels before prediction.
#' @param ml_fit A trained machine learning model (e.g., from randomForest) that can be used with the predict() function.
#' @param feature_df A data frame containing features for each sequence. Must include columns: "gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", and "DNA_5mer".
#' @param positive_threshold A numeric value between 0 and 1. Probabilities above this threshold are classified as "Positive". Default is 0.5.
#' @return The original feature_df augmented with two new columns: "predicted_m6A_prob" (the probability of being positive) and "predicted_m6A_status" ("Positive" or "Negative").
#' @export
#' @import randomForest
#' @examples
#' \dontrun{
#' # Assuming 'rf_fit.rds' and 'example.csv' are in inst/extdata
#' model_file <- system.file("extdata", "rf_fit.rds", package="m6APrediction")
#' data_file <- system.file("extdata", "example.csv", package="m6APrediction")
#' ml_model <- readRDS(model_file)
#' my_data <- read.csv(data_file)
#'
#' predictions_df <- prediction_multiple(ml_model, my_data)
#' head(predictions_df)
#' }
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame

  # Encode DNA sequence and set factor levels for categorical variables

  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  prob_test <- predict(ml_fit,newdata=feature_df,type="prob")[,"Positive"]
  prob_status <- ifelse(prob_test>positive_threshold,"Positive","Negative")


  feature_df <- cbind(feature_df,predicted_m6A_prob= prob_test,predicted_m6A_status=prob_status)


  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}




#' @title Predict m6A Site for a Single Sequence
#' @description Predicts m6A modification status for a single sequence based on its features. It creates a one-row data frame, processes it, and returns a named vector with the prediction.
#' @param ml_fit A trained machine learning model.
#' @param gc_content Numeric. The GC content of the sequence.
#' @param RNA_type Character. The type of RNA (e.g., "mRNA").
#' @param RNA_region Character. The region of the RNA (e.g., "CDS").
#' @param exon_length Numeric. Length of the exon.
#' @param distance_to_junction Numeric. Distance to the nearest splice junction.
#' @param evolutionary_conservation Numeric. A score for evolutionary conservation.
#' @param DNA_5mer Character. The 5-mer DNA sequence.
#' @param positive_threshold A numeric value between 0 and 1. Probabilities above this threshold are classified as "Positive". Default is 0.5.
#' @return A named vector with two elements: "predicted_m6A_prob" (the probability) and "predicted_m6A_status" ("Positive" or "Negative").
#' @export
#' @import randomForest
#' @examples
#' \dontrun{
#' # Assuming 'rf_fit.rds' is in inst/extdata
#' model_file <- system.file("extdata", "rf_fit.rds", package="m6APrediction")
#' ml_model <- readRDS(model_file)
#'
#' single_prediction <- prediction_single(ml_model,
#'                                       gc_content = 0.5,
#'                                       RNA_type = "mRNA",
#'                                       RNA_region = "CDS",
#'                                       exon_length = 100,
#'                                       distance_to_junction = 50,
#'                                       evolutionary_conservation = 0.8,
#'                                       DNA_5mer = "ATCGG")
#' print(single_prediction)
#' }
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  #Complete this function by writing code at the `___`

  feature_df <- data.frame(gc_content,RNA_type,RNA_region,exon_length,distance_to_junction, evolutionary_conservation,DNA_5mer)
  feature_df <- cbind(feature_df,dna_encoding(feature_df$DNA_5mer))

  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  prob_test <- predict(ml_fit,newdata=feature_df,type="prob")[,"Positive"]
  prob_status <- ifelse(prob_test>positive_threshold,"Positive","Negative")


  returned_vector <- c(predicted_m6A_prob=prob_test,predicted_m6A_status=prob_status)

  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}
