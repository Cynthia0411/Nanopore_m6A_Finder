# Nanopore_m6A_Finder

Nanopore_m6A_Finder is a machine‑learning based pipeline to identify **m6A sites** and estimate **m6A levels** from Nanopore Direct RNA Sequencing (DRS) data.

It supports detecting m6A in both **coding regions (CDS)** and **poly(A) tails** by combining basecalling features with trained **XGBoost** and **Random Forest** models.

---

## Key Features

- **End‑to‑end m6A calling pipeline** for Nanopore DRS datasets
- **Two regions supported**:
  - Coding regions (CDS)
  - Poly(A) regions
- **Supervised learning models**:
  - XGBoost model
  - Random Forest model
- **Site‑level m6A level estimation**:
  - Per‑site read counts, coverage and m6A levels
  - WT vs writer‑KO comparison to find **reliable m6A sites**
- **Scripted workflow**:
  - Shell scripts + Python scripts for each step
  - Easy to integrate into existing analysis pipelines

---

## Repository Structure (Main Scripts)

Below is a brief overview of important scripts in this repository:

### Core m6A calling pipeline

- `get_pos_CDs.py`  
  Extract coding region (CDS) positions from aligned SAM files.

- `get_pos_polyA.py`  
  Extract poly(A) positions from aligned SAM files.

- `Extract_feature_per_5mer.py`  
  Extract signal / basecalling features per 5‑mer from DRS reads (using fast5 + positional info).

- `Extract_feature_per_read.py`  
  Extract features at the read level (alternative / complementary feature representation).

- `basecall_feature.py`  
  Encode the extracted per‑5‑mer features into numeric feature matrices for ML models.

- `xgb_model_pred.py`  
  Apply the trained **XGBoost** model to 5‑mer features to obtain methylation scores.

- `RF_combine_xgb_predict.py`  
  Combine **Random Forest** and **XGBoost** predictions to obtain final m6A calls (mainly for poly(A) pipeline).

- `get_chr_pos.py`  
  Convert read‑level predictions into genomic coordinates (chr, position, strand).

- `get_interval_and_read_count.py`  
  Aggregate predictions at the site level:
  - Identify m6A candidate sites
  - Count supporting reads / intervals

- `calculate_coverage.sh`  
  Shell script to compute coverage per candidate m6A site.

- `merge_to_get_coverage.py`  
  Merge site‑level predictions with coverage to compute **m6A levels**.

- `merge_WT_KO.py`  
  Combine WT and writer‑KO datasets to compare m6A levels per site.

- `bootstrap_new.py` + `config_bootstrap.py`  
  Perform bootstrap‑based statistical testing to identify **reliable m6A sites** whose WT m6A levels are significantly higher than in KO.

- `calculate_overlap.py`, `calculate_overlap_bed_tsv.py`, `calculate_overlap_form_bed.py`, `calculate_overlap_from_bed_to_csv.py`  
  Utilities to compute overlaps between different site sets (BED / TSV / CSV formats).

### Data processing and utilities

- `basecall_feature_predict.py`  
  Variant of prediction step on encoded basecall features.

- `get_AAAAA_pos.py`, `get_level_70_AAAAA_pos.py`  
  Utilities related to poly(A) related motifs / thresholds.

- `calculate_coverage.sh`, `tsv_to_bed.sh`, `bed_to_gene.sh`  
  Shell utilities for coverage calculation, format conversion and gene annotation.

### Model training & evaluation

- `train_imbalanced_xgb.py`  
  Train an imbalanced **XGBoost** model on 5‑mer features.

- `train_random_forest_model.py`  
  Train the **Random Forest** model.

- `train_SVM_model.py`, `train_bayes_model.py`, `train_bayes_model_cn.py`, `train_hnsw_model.py`  
  Alternative model training scripts (SVM, Bayesian, HNSW‑based methods).

- `draw_roc_curve_for_the_two_models_xgb_rf.py`  
  Plot ROC curves to compare XGBoost and Random Forest models.

- `hardvoting_xgb_rf_eval.py`  
  Hard‑voting ensemble evaluation between XGBoost and Random Forest.

### File format & helper scripts

- `csv_to_bed.py`  
  Convert CSV site tables to BED format.

- `merge_WT_KO_site_coverage.py`, `merge_3_tsv_files.py`, `merge_to_get_coverage.py`  
  Various merge utilities for site coverage and level calculation.

- `extract_read_kmers_from_bed.py`, `select_reads.py`, `select_1500_fastq.py`  
  Read / k‑mer selection helpers.

- `change_file_name.py`, `whether_pos_in_1308.py`  
  Small helper scripts for file management and site filtering.

---

## Pretrained Models & Training Data

The pretrained XGBoost and Random Forest models, as well as the training data, are provided via Google Drive:

- Models and training data:  
  `https://drive.google.com/drive/folders/1w7Ff1Ay76DVXK-Mk_Mi9LF1csrQqcGgq?usp=drive_link`

Training data were derived from in vitro transcribed 5‑mers with As or m6As (e.g. `81_1_1500_top100.sam`, `82_1_1500.sam`).

---

## Installation

### Requirements

- Python 3.x
- Typical scientific Python stack:
  - `numpy`, `scipy`, `pandas`
  - `scikit-learn`
  - `xgboost`
  - `h5py` (for fast5‑related operations, depending on your environment)
- Nanopore tools and dependencies:
  - Aligned **SAM/BAM** files from Nanopore DRS
  - Access to raw **fast5** files
- Shell environment with common tools:
  - `bash`, `samtools`, etc. (depending on your existing pipeline)

> **Note:** Exact version requirements may depend on your environment.  
> We recommend using a dedicated **conda** or **virtualenv** environment.

### Example: create a conda environment

```bash
conda create -n nanopore_m6a python=3.10
conda activate nanopore_m6a

pip install numpy scipy pandas scikit-learn xgboost h5py
```

---

## Workflow Overview

### 1. m6A calling in coding regions (CDS)

1. **Extract CDS positions from aligned reads**

   ```bash
   python get_pos_CDs.py \
       --input your_aligned.sam \
       --output cds_positions.tsv
   ```

2. **Extract features per 5‑mer using fast5 and positional info**

   ```bash
   python Extract_feature_per_5mer.py \
       --pos cds_positions.tsv \
       --fast5_dir path/to/fast5 \
       --out_dir features/
   ```

3. **Encode basecalling features**

   ```bash
   python basecall_feature.py \
       --feature_dir features/ \
       --output encoded_features.npy
   ```

4. **Predict m6A scores with pretrained XGBoost model**

   ```bash
   python xgb_model_pred.py \
       --model xgb_model.pkl \
       --features encoded_features.npy \
       --output xgb_pred.tsv
   ```

5. **Convert predictions to genomic coordinates**

   ```bash
   python get_chr_pos.py \
       --pred xgb_pred.tsv \
       --output pred_with_chr_pos.tsv
   ```

6. **Get per‑site interval and read counts**

   ```bash
   python get_interval_and_read_count.py \
       --input pred_with_chr_pos.tsv \
       --output m6a_sites_with_counts.tsv
   ```

7. **Compute coverage for each candidate site**

   ```bash
   bash calculate_coverage.sh
   ```

8. **Merge predictions with coverage to estimate m6A level**

   ```bash
   python merge_to_get_coverage.py \
       --sites m6a_sites_with_counts.tsv \
       --coverage coverage.tsv \
       --output m6a_sites_with_level.tsv
   ```

9. **Merge WT and KO datasets**

   ```bash
   python merge_WT_KO.py \
       --wt wt_m6a_sites.tsv \
       --ko ko_m6a_sites.tsv \
       --output merged_WT_KO.tsv
   ```

10. **Bootstrap‑based filtering of reliable m6A sites**

    ```bash
    python bootstrap_new.py \
        --config config_bootstrap.py \
        --input merged_WT_KO.tsv \
        --output reliable_m6a_sites.tsv
    ```

### 2. m6A calling in poly(A) regions

1. Extract poly(A) positions:

   ```bash
   python get_pos_polyA.py \
       --input your_aligned.sam \
       --output polyA_positions.tsv
   ```

2. Run the same **feature extraction + encoding** as for CDS:

   ```bash
   python Extract_feature_per_5mer.py ...
   python basecall_feature.py ...
   ```

3. Combine Random Forest and XGBoost predictions:

   ```bash
   python RF_combine_xgb_predict.py \
       --xgb_model xgb_model.pkl \
       --rf_model rf_model.pkl \
       --features encoded_features.npy \
       --output polyA_m6a_pred.tsv
   ```

4. Downstream steps (coordinate conversion, coverage and level calculation) are analogous to the CDS pipeline.

---

## Notes & Recommendations

- This pipeline assumes **Nanopore DRS** data with corresponding **aligned SAM/BAM** and **fast5**.  
- For large datasets, we strongly recommend running the pipeline on a compute node or server with sufficient memory and storage.
- Training scripts (`train_imbalanced_xgb.py`, `train_random_forest_model.py`, etc.) can be used to:
  - Retrain models with your own in vitro datasets
  - Adapt the pipeline to other organisms or experimental conditions

---

## Citation

If you use Nanopore_m6A_Finder in your work, please cite the corresponding paper (add bibliographic information here when available).

---

## License

This project is licensed under the **MIT License**.

See the [`LICENSE`](./LICENSE) file for the full license text.

---

## Contact

For questions, issues, or suggestions, please open an issue on this GitHub repository:  
`https://github.com/Cynthia0411/Nanopore_m6A_Finder`

