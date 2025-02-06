# circadian-rhythm-
Codes for JTK analysis and harmonic regression to check and fit rhythm in the data

# Harmonic Regression Analysis  

## 📌 Overview  
This repository contains the **R script and datasets** used to analyze **Any gene expression over time** using **harmonic regression**. This method models periodic fluctuations in **GLUT4 expression** to study its circadian rhythm.  

## 📁 Repository Structure  

📂 Harmonic_Regression/ │── 📜 README.md # Project documentation
│── 📜 harmonic_regression.R # Main R script for harmonic regression analysis
│── 📁 data folder # Input dataset containing GLUT4 expression over time
│── 📁 output # Generated plot of expression
  
## 🛠️ Requirements & Installation  
To run this analysis, you need **R (v4.2.1)** and the following R packages:  

```r
install.packages(c("ggplot2", "minpack.lm", "dplyr"))

git clone (https://github.com/RashmiSivasengh/circadian-rhythm)

source("harmonic_regression.R")


Harmonic Regression Model
The model follows the equation:

𝑦
=
𝐴
⋅
sin
⁡
(
𝜔
𝑡
+
𝜙
)
+
𝐶
y=A⋅sin(ωt+ϕ)+C
where:

A = Amplitude
ω = Frequency component
φ = Phase shift
C = Baseline expression
The script uses nlsLM (Nonlinear Least Squares) for curve fitting.

# JTK_CYCLE Analysis of PER3 Gene

## 📌 Overview
This repository contains **JTK_CYCLE analysis** for detecting **circadian rhythmicity** in ** gene expression** under **Insulin treatment conditions**. The analysis identifies **periodic patterns**, estimating **adjusted p-values (BH.Q), period length (PER), lag phase (LAG), and amplitude (AMP).**

## 📁 Repository Structure

## 🛠️ Dependencies & Installation  
To run this analysis, you need **R (v4.2.1 or later)** with the following dependencies:  

```r
options(stringsAsFactors=FALSE)
source("JTK_CYCLEv3.1.R")

JTK_Analysis/ │── 📜 README.md # Project documentation │── 📜 JTK_CYCLEv3.1.R # JTK_CYCLE core script │── 📜 JTK_PER3_Analysis.R # Main R script for analysis │──
📄 X(Anydate).txt # Input dataset (gene expression over time) │── 📄 results.txt # Processed JTK results (ordered by adj.P & AMP) │── 📂 results/ # Directory for output files


## 🛠️ Dependencies & Installation  
To run this analysis, you need **R (v4.2.1 or later)** with the following dependencies:  

```r
options(stringsAsFactors=FALSE)
source("JTK_CYCLEv3.1.R")

git clone 
cd JTK_PER3_Analysis https://github.com/RashmiSivasengh/circadian-rhythm

source("JTK_PER3_Analysis.R")
A JTK results.txt file containing ranked circadian genes by BH-adjusted p-values (BH.Q) and amplitude (AMP).
The script estimates optimal period (PER) and phase shift (LAG).


# Code for third analysis, re-analysis of public data
library(GEOquery)
library(EnsDb.Hsapiens.v79)
library(dplyr)
library(ggplot2)
library(gridExtra)




# Read Batch-Corrected logCPM tsv file from 
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE182121 
# GEO accession code GSE182121 

ccount <- read.table("/Users/andrewscott/Downloads/GSE182117_logCPM_batchCorrected.tsv")

# restructure dataframe for ease of use
colnames(ccount) <- ccount[1,]
ccount <- ccount[-1,]
colnames(ccount)[1] <- "EnsembleIDs"

# extract Ensemble IDs from df
ensembleIDs <- ccount[,1]

# obtain corresponding HGNC Gene Symbols
geneSymbols <- AnnotationDbi::select(EnsDb.Hsapiens.v79, 
                      key= ensembleIDs, 
                      columns=c("SYMBOL"), 
                      keytype="GENEID"
                      )

# add gene symbols to the dataframe
ccount <- left_join(ccount, 
                    geneSymbols,
                    join_by(EnsembleIDs == GENEID)
                    )

# reorder the dataframe for ease of use
ccount <- dplyr::select(ccount,
                        "SYMBOL",
                        everything()
                        )

# extract the salient dataframe subset
genes_of_interest <- c("PER3", "ARNTL", "HOXB5", "TSSK6")
ccount_sub <- dplyr::filter(ccount,
                            ccount[,1] %in% genes_of_interest
                            )

# convert the logCPM values the salient subset to Z-scores
# and splice back together

# remove the "non-values columns"
df <- ccount_sub[,-c(1:2)]

# make these values doubles 
df <- mutate_all(df, function(x) as.numeric(as.character(x)))

# transpose for z-score calculations
df = t(df)

# convert df to z scores using the scale() function
df_z <- scale(df, scale = TRUE, center = TRUE)

# re-transpose to original form
df_z <- t(df_z)

# add back in the gene names
ccount_sub_z <- cbind(ccount_sub[,1:2],
                      df_z
                      )

# transpose table to make it easier to add disease condition and time information
ccount_sub_z <- t(ccount_sub_z)


# create a vector of values by searching the rownames for T2D status indicator
T2D_status <- rep("NGT", length(ccount_sub_z))
T2D_status[which(grepl("T2D",rownames(ccount_sub_z)))] <- "T2D"

# attach to ccount_sub_z df
ccount_sub_z <- cbind(ccount_sub_z, T2D_status)

# repeat this for the participant number
participant_number <- c(rep(NA,   2), 
                        rep("5",  8),
                        rep("5",  8),
                        rep("6",  8),
                        rep("6",  8),
                        rep("1",  7),
                        rep("1",  7),
                        rep("2",  8),
                        rep("2",  8),
                        rep("3",  8),
                        rep("3",  8),
                        rep("4",  8),
                        rep("4",  8),
                        rep("12", 7),
                        rep("12", 7),
                        rep("7",  8),
                        rep("7",  8),
                        rep("8",  8),
                        rep("8",  8),
                        rep("9",  7),
                        rep("9",  7),
                        rep("10", 8),
                        rep("10", 8),
                        rep("11", 8),
                        rep("11", 8)
                        )

ccount_sub_z <- cbind(ccount_sub_z, participant_number)

## add the glucose high or control condition
glucose_condition <-  c(rep(NA,         2), 
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  7),
                        rep("High",     7),
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  7),
                        rep("High",     7),
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  7),
                        rep("High",     7),
                        rep("Control",  8),
                        rep("High",     8),
                        rep("Control",  8),
                        rep("High",     8)
                        )

ccount_sub_z <- cbind(ccount_sub_z, glucose_condition)

# add biopsy time information
# NB spaces show where gaps occur in the data. 
# the number after the hash is the number of times there 
# are for that row/participant

biopsy_times      <-    c(rep(NA, 2), 
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,   30,36,42,48,54), #7
                        c(12,18,   30,36,42,48,54), #7
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,   24,30,36,42,48,54), #7
                        c(12,18,24,30,36,42,48   ), #7
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,   48,54), #7
                        c(12,   24,30,36,42,48,54), #7
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54), #8
                        c(12,18,24,30,36,42,48,54)  #8
                        )

ccount_sub_z <- cbind(ccount_sub_z, biopsy_times)

# tidy up the matrix and convert to df
colnames(ccount_sub_z)[1:4] <- c("PER3", "ARNTL", "HOXB5", "TSSK6")
ccount_sub_z <- ccount_sub_z[-(1:2),]
ccount_sub_z <- as.data.frame(ccount_sub_z)

# ensure that logCPM values and times are coded as doubles
ccount_sub_z$PER3         <- as.numeric(ccount_sub_z$PER3)
ccount_sub_z$ARNTL        <- as.numeric(ccount_sub_z$ARNTL)
ccount_sub_z$HOXB5        <- as.numeric(ccount_sub_z$HOXB5)
ccount_sub_z$TSSK6        <- as.numeric(ccount_sub_z$TSSK6)
ccount_sub_z$biopsy_times <- as.numeric(ccount_sub_z$biopsy_times)

# make a data subset with only the "Control" glucose condition
ccount_sub_z_c<- dplyr::filter(ccount_sub_z, glucose_condition == "Control")

# custom function to perform harmonic regression and plot
plot_harmonic_regression <- function(data, gene, Tau = 24) {
  # filter data for the specific gene
  gene_data <- data %>% filter(dataset == gene)
  
  # perform harmonic regression for NGT
  model_NGT <- lm(y ~ sin(2*pi/Tau * x) + cos(2*pi/Tau * x), data = gene_data[gene_data$z == "NGT", ])
  
  # perform harmonic regression for T2D
  model_T2D <- lm(y ~ sin(2*pi/Tau * x) + cos(2*pi/Tau * x), data = gene_data[gene_data$z == "T2D", ])
  
  # create a data frame for smooth predictions
  pred_data <- data.frame(x = seq(min(gene_data$x, na.rm = TRUE), max(gene_data$x, na.rm = TRUE), length.out = 100))
  
  # predict for NGT
  pred_data$y_NGT <- predict(model_NGT, newdata = pred_data)
  
  # predict for T2D
  pred_data$y_T2D <- predict(model_T2D, newdata = pred_data)
  
  # custom plotting
  p <- ggplot(gene_data, aes(x = x, y = y, color = z)) +
    geom_point(alpha = 0.5) +
    geom_line(data = pred_data, aes(x = x, y = y_NGT), color = "black") +
    geom_line(data = pred_data, aes(x = x, y = y_T2D), color = "red") +
    scale_color_manual(values = c("T2D" = "red", "NGT" = "black")) +
    scale_x_continuous(breaks=seq(12,54,6)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey80"),
          legend.position = "right") +
    labs(x="Time (Hours)",
         y = "Expression (logCPM Z-Scores)",
         title = "") +
    ylim(-3,3)
  
  return(p)
}

# prepare data for plotting
plot_data <- list(
  transmute(ccount_sub_z_c, x=biopsy_times, y=ARNTL, z = T2D_status, dataset="ARNTL"),
  transmute(ccount_sub_z_c, x=biopsy_times, y=HOXB5, z = T2D_status, dataset="HOXB5"),
  transmute(ccount_sub_z_c, x=biopsy_times, y=PER3,  z = T2D_status, dataset="PER3"),
  transmute(ccount_sub_z_c, x=biopsy_times, y=TSSK6, z = T2D_status, dataset="TSSK6")
) %>%
bind_rows()

# plot each gene separately
plots <- lapply(c("ARNTL", "HOXB5", "PER3", "TSSK6"), function(gene) {
  plot_harmonic_regression(plot_data, gene)
})

# arrange plots in a grid
grid.arrange(grobs = plots, nrow = 2, top = "Circadian Gene Expression")

📌 Citation
If you use this code in your research, please cite:

[Live cell GLUT4 translocation assay reveals Per3 as a novel regulator of circadian insulin sensitivity  in skeletal muscle cells], [Rashmi Sivasengh.,Andrew scott., Brendan Gabriel*], [2025]

📬 Contact
For questions, please contact [iamrashmi96@gmail.com] or open an issue in this repository.
