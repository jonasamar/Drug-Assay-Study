This file is here to explain the order in which one should run the multiple Rmd files (when values are equals it means it doesn't matter which one is first):

1 - Onset of Labor preprocessing
2 - Onset of Labor curves classification
3 - Onset of Labor feature index and patient clustering (except section : Simulated drug effect visualization)

4 - Drug assay preprocessing
5 - Drug assay data visualization
6 - Drug assay sigmoids fitting
7 - Drug assay drug tensors

8 - Drug assay scores and visualization
9 - Onset of Labor feature index and patient clustering (all sections) : switch predicting.simulated.data to TRUE
