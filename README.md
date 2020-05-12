# E-Py-Genetics

![logo](gif/sequencing_plot.gif)

Detection of epigenetic DNA modification is important for several research applications, including cancer research and analysis of evolutionary mechanisms. However, currently established methods are not powerful enough. To tackle this problem, an AI-based analysis of nanopore sequencing data was implemented in this project.

Data was received from https://github.com/tleonardi/nanocompore/

[Google Slides](https://docs.google.com/presentation/d/1usy1A-_OmuS3tVYgu_4owRPxRHVGXeNvdjLj_fFP6Yc/edit?usp=sharing) version of Presentation

# Structure of this repository:

## Jupyter Notebooks:

### 1_Data_Mining_EDA

* Data import
* Feature engineering
* Data exploration
* Yeast data import
* Data export

### 2_Predictive_Modeling_Uncombined_DataFrame

* Data import
* Train test split
* Modeling

### 3_Predictive_Modeling_Combined_DataFrame

* Data import
* Train test split
* Modeling

### 4_RNN (work in progress)

* Data import
* Packages and Train Test Split
* Separate Data
* Separate Data with reduced Feature set

### 5_Visualizations

* Data Import
* UMAP and tSNE
* Dwell Time to Mod Status
* Effect of Base Identity
* Yeast Data Prediction
* DNA Feature GIF

## data Folder

* RAW Data received from Nanocompore repository

## figures Folder

* figures generated in notebook 5

## GIF

* gif generated in notebook 5
* explanatory gif for nanopore sequencing ([source](https://www.youtube.com/watch?v=E9-Rm5AoZGw&t=90s))
