## Requirements

The R scripts in this folder have been written and should be executed in the order described below and with the following R version :

platform       x86_64-apple-darwin20       
arch           x86_64                      
os             darwin20                    
system         x86_64, darwin20            
status                                     
major          4                           
minor          3.0                         
year           2023                        
month          04                          
day            21                          
svn rev        84292                       
language       R                           
version.string R version 4.3.0 (2023-04-21)
nickname       Already Tomorrow

You can check your own R version by running the command "version" in you R console.
You can doawnload this version of R for windows here : https://cran.r-project.org/bin/windows/base/old/4.3.0/
You can doawnload this version of R for macOS here : https://cran.r-project.org/bin/macosx/

We strongly advice to also use R Studio to go through all the scripts and their comments which is easy to doawnload from this website : https://posit.co/downloads/

## R scripts order

This file is here to explain the order in which one should run the multiple R/Rmd files (when values are equals it means it doesn't matter which one is first):

0 - install.R

1 - Onset of Labor preprocessing
2 - Onset of Labor curves classification
3 - Onset of Labor feature index and patient clustering (except section Simulated drug effect visualization)

4 - Drug assay preprocessing
5 - Drug assay data visualization
6 - Drug assay sigmoids fitting
7 - Drug assay drug tensors

8 - Drug assay scores and visualization
9 - Onset of Labor feature index and patient clustering (all sections)
