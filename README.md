README
================

*M6APrediction*

**Package Overview**

The **m6APrediction** package provides tools for predicting **m6A RNA
modification sites**  
based on sequence and structural features using a **Random Forest
model**.

It includes functions for: - Model training and validation  
- Batch prediction on multiple samples (`prediction_multiple()`)  
- Single-sequence prediction (`prediction_single()`)  
- Performance visualization (ROC and PRC curves)

------------------------------------------------------------------------

*Installation*

You can install the development version of **m6APrediction** from
GitHub:

``` r
install.packages("devtools")
devtools::install_github("isaWLO/m6APrediction")
```

*Example: Using m6APrediction package*

library(m6APrediction)

``` r
# Load pre-trained model
model_file <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
ml_model <- readRDS(model_file)

# ---- Single sequence prediction ----
single_result <- prediction_single(
  ml_fit = ml_model,
  gc_content = 0.5,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 100,
  distance_to_junction = 50,
  evolutionary_conservation = 0.8,
  DNA_5mer = "ATCGG"
)

print(single_result)

# ---- Multiple sequences prediction ----
data_file <- system.file("extdata", "example.csv", package = "m6APrediction")
example_df <- read.csv(data_file)

multi_result <- prediction_multiple(ml_fit = ml_model, feature_df = example_df)

head(multi_result)
```
