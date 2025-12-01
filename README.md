# GammaAnalyzer

A medical physics tool for visualizing CT + RTDOSE datasets and computing 3D Gamma Index
between two dose distributions. Designed for clinical research and QA comparison between
OriginalCT-Plan vs. SynCT-Plan (or other dose sets).

---

## 📂 Required Folder Structure

The root folder you select in the GUI must follow this format:


✔ Each scan subfolder must include:
- CT Image Storage DICOM files
- One (or more) RT Dose DICOM file (largest used automatically)

---

## 🔍 Features

- Automatic detection of patient folders + scan types
- CT + Dose loading into memory
- Interactive slice viewer (**CT + dose overlay**)
- 3D Gamma Index calculation
- Export of results:
  - `.npy`
  - `.csv`
  - `.png`
  - Summary `.txt`

---

## 🖥 Running the App

```bash
python GammaAnalyzer.py
