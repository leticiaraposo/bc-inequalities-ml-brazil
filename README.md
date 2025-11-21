
# Geographical and Socioeconomic Inequalities in Breast Cancer Mortality: A Machine Learning Analysis in Brazil

This repository contains the data and reproducible R code supporting the article:

**Geographical and Socioeconomic Inequalities in Breast Cancer Mortality: A Machine Learning Analysis in Brazil**

---

## Overview

This study investigates how geographical, socioeconomic, demographic, and health service–related factors shape the spatial distribution of breast cancer mortality in Brazil. Using Immediate Geographic Regions (IGRs) as the unit of analysis, we apply machine learning methods — including Random Forest (RF), Geographical Random Forest (GRF), and Shapley-value interpretability — to capture both global and spatially varying predictor importance.

The repository provides the complete workflow for data preparation, processing, and modeling.

---

## Repository Structure

### **data/raw/**
Contains the original dataset:
- `dados_CM_Ferrete_Raposo.xlsx`

### **data/processed/**
Contains the cleaned and processed dataset used in the analyses:
- `base_cancer_mama_imputados_sem_correlacao_alta_ri.rds`

### **R/**
Contains R scripts used in the study:
- `Dataset_assembly_BC_Ferrete_Raposo.R`
- `Data_aggregation_ASMR_BC_Ferrete_Raposo.R`
- `Multicollinearity_Imputation_BC_Ferrete_Raposo.R`
- `Modeling_BC_Ferrete_Raposo.R`
---

## Data Description

The dataset compiles indicators from multiple Brazilian sources:

- **Ipeadata (IPEA)** – socioeconomic indicators  
- **Proadess (Fiocruz)** – health services and infrastructure  
- **SIDRA (IBGE)** – demographic and socioeconomic census indicators  
- **DATASUS / Tabnet (SUS)** – mortality and health system data  

Variables include sociodemographic composition, socioeconomic structure, urbanization metrics, and availability of health professionals and services.

---

## Methodological Summary

### **Machine Learning Models**
- **Random Forest (RF)** for global variable importance  
- **Geographical Random Forest (GRF)** for spatially varying importance  
- **Shapley values** for model-agnostic interpretability  

### **Performance Evaluation**
- Out-of-bag (OOB) estimates  
- RMSE and R² for global models  
- Pseudo-R² for GRF local submodels  

### **Spatial Scale**
- Immediate Geographic Regions (IGRs) as recommended by IBGE for regional analysis

