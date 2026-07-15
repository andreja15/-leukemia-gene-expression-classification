# leukemia-gene-expression-classification

**Predicting Leukemia Type from Gene Expression Profiles Using Machine Learning**

**Overview**

This project compares five supervised learning methods to classify leukemia type, ALL vs. AML, using the Golub et al. (1999) gene expression dataset. The dataset includes 7,129 gene expression measurements across 72 patients. The goal is to evaluate predictive performance (AUC, sensitivity, specificity, accuracy) and identify genes that are consistently important predictors across methods.

**Methods**

Ridge Regression, an L2-penalized logistic regression that retains all predictors
Lasso Regression, an L1-penalized logistic regression that performs variable selection
Elastic Net Regression, combining L1 and L2 penalties (alpha = 0.9)
Sure Independence Screening + Lasso (SIS + Lasso), a two-stage screening and selection approach
Random Forest, a tree-based ensemble method that captures non-linear relationships

**Data**

Source: Golub et al. (1999) leukemia dataset, publicly available on Kaggle.
Training set: 38 patients (27 ALL, 11 AML)
Independent test set: 34 patients (20 ALL, 14 AML)
Features: 7,129 gene expression measurements per patient, reduced to 5,064 after filtering
Outcome: leukemia type (ALL = 0, AML = 1)

**Preprocessing**

Transpose expression matrices so patients are rows and genes are columns
Remove 58 Affymetrix control probes (AFFX prefix)
Filter out low-variability genes using training set SD (SD of 100 or less removed), reducing genes from 7,071 to 5,064
Standardize all features to mean 0 and SD 1 using training set parameters, applied to both train and test sets

**Results**

Ridge Regression achieved the highest AUC (0.9964) but the lowest sensitivity (0.5714), indicating a conservative model, as the L2 penalty pulls predictions toward the center.

Lasso and Elastic Net both achieved an AUC of 0.9750, a sensitivity of 0.7143, and an accuracy of 0.8824. Lasso selected 13 genes, and Elastic Net selected 15.

Random Forest achieved an AUC of 0.9929, a sensitivity of 0.5714, and an accuracy of 0.8235.

SIS + Lasso achieved an AUC of 0.9786, the highest sensitivity (0.7857), and the highest accuracy (0.9118), selecting 14 genes.

All five methods achieved a specificity of 1.0, meaning every ALL patient in the test set was correctly identified. All methods reached an AUC of at least 0.975, showing that gene expression profiles are highly informative for distinguishing ALL from AML. SIS + Lasso is arguably the most clinically relevant method given its strength in correctly identifying AML cases.

**Consensus Genes**

Y12670_at and X95735_at were selected as important predictors by all five methods, making them the most robust markers for distinguishing ALL from AML. Several other genes, including U50136_rna1_at, D49950_at, U82759_at, M19507_at, and M23197_at, were selected by four of the five methods.

**Limitations**

The dataset includes only 72 patients, which limits statistical power. The data were collected in 1999 using Affymetrix microarray technology, which may not generalize well to modern sequencing platforms. All methods were evaluated on a single dataset from one source, which limits generalizability.
