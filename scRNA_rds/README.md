# Download scRNA-seq Seurat RDS Files

This document describes how to download single-cell RNA-seq datasets in **Seurat RDS format** used in this project.

## 1. Datasets from scPlantDB

Most scRNA-seq Seurat objects can be downloaded directly from **scPlantDB**:

🔗 https://biobigdata.nju.edu.cn/scplantdb/dataset

Please note that **all datasets are available from scPlantDB except**:

- `GSE152766`
- `GSE297576`

These two datasets must be downloaded from **NCBI GEO**, as described below.

---

## 2. Datasets from NCBI GEO

### 2.1 GSE152766 (Root Atlas)

- **Accession**: GSE152766  
- **Format**: Seurat `.rds.gz`

**Download link**:  
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152766/suppl/GSE152766_Root_Atlas.rds.gz

**Command-line download example**:
```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152766/suppl/GSE152766_Root_Atlas.rds.gz
gunzip GSE152766_Root_Atlas.rds.gz
