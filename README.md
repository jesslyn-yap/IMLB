# Prediction of transcription factor binding sites

# Project Overview

This project predicts CTCF transcription factor binding on 200 bp DNA sequences in the GM12878 cell line, comparing CNN model with LightGBM.

> **Status:** Work in Progress  
> **Note:** This project is part of a course assignment. The code, analysis, and results are still evolving.

---

## 1. Machine Learning Task

Binary classification task:

> Predict whether CTCF binds to a given 200 bp genomic bin.

Two model types are implemented:

* **CNN** — learns directly from raw DNA sequence
* **LightGBM** — uses engineered genomic features for interpretability

---

## 2. Dataset

This dataset construction is inspired by:
- [ENCODE-DREAM Challenge](https://www.synapse.org/Synapse:syn6131484/wiki/402026)
- [DeepGRN paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03952-1)


Data sources:

* CTCF optimal (IDR) peaks: [ENCFF710VEH](https://www.encodeproject.org/files/ENCFF710VEH/)
* CTCF replicate 1: [ENCFF256QBB](https://www.encodeproject.org/files/ENCFF256QBB/)
* CTCF replicate 2: [ENCFF307KRD](https://www.encodeproject.org/files/ENCFF307KRD/)
* Reference genome (hg19): UCSC download

Labeling strategy:

* Positive = bins intersecting optimal (IDR) peaks
* Ambiguous = intersect replicate peaks but not IDR → discarded
* Negative = bins with no overlap with any replicate

Genomic binning:

* Fixed bin size: 200 bp
* Only chromosomes chr1–chr22

Train/validation/test split (chromosome-based):

* Train: 2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,20,22
* Validation: 14, 19
* Test: 1, 21

---

## 3. Models

### Convolutional Neural Network (CNN)

* Input: one-hot encoded 200 bp sequence
* Learns DNA motifs & local patterns

### LightGBM

* Input: engineered features

  * One-hot encoding (summaries, not flat 800-dim one-per-base)
  * GC content
  * CpG density
  * (optional) PWM motif scores
  * (optional) DNase-seq signal (averaged in the bin)

LightGBM is used for interpretability and speed on tabular genomic features.

---

## 4. Feature Engineering

### CNN

* Raw DNA sequence → one-hot encoded
* No handcrafted features

### LightGBM Features

* Base composition features
  * Fraction A/C/G/T
  * GC content
  * CpG density
* PWM scores (optional)

---

## 5. Regularization

### CNN

* Dropout
* Early stopping

### LightGBM

* max_depth
* num_leaves
* min_data_in_leaf
* lambda_l1, lambda_l2
* feature_fraction
* bagging_fraction, bagging_freq

Hyperparameter tuning improves generalization, especially due to heavy class imbalance.

---

## 6. Evaluation Metrics

* ROC-AUC
* PR-AUC (more informative for imbalanced data)

Evaluation is done separately on validation and held-out test chromosomes.

---

## 7. Possible Extensions

If time permits:

* Single-task learning approach to predict other TF in GM12878
* Incorporate k-mer embeddings into LightGBM
* Add DNase-seq features
