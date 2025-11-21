
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

---


```r
dados <- readRDS("data/processed/base_cancer_mama_imputados_sem_correlacao_alta_ri.rds")
