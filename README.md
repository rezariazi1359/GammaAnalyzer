GammaAnalyzer – Usage Instructions
---------------------------------------------

1) Root Folder Structure
 • Choose a ROOT folder that contains one folder per PATIENT.
 • Inside each patient folder, create numbered SCAN subfolders (Scan1, Scan2, … ScanN).
 • Each scan folder should contain all DICOM files directly in its root:
   - CT series (DICOM with CT Image Storage UID)
   - RT Dose file (DICOM with RT Dose Storage UID)

 Example:
 MyData/
 ├── P001/
 │   ├── Scan1/
 │   │   CT*.dcm
 │   │   RTDOSE*.dcm
 │   ├── Scan2/
 │   │   CT*.dcm
 │   │   RTDOSE*.dcm
 │   └── Scan3/
 │       CT*.dcm
 │       RTDOSE*.dcm
 ├── P002/
 │   ├── Scan1/
 │   │   CT*.dcm
 │   │   RTDOSE*.dcm
 │   └── Scan2/
 │       CT*.dcm
 │       RTDOSE*.dcm
 └── ...

2) Workflow
 • Launch GammaAnalyzer.exe.
 • Click 'Load Patient Data' and select the ROOT folder.
 • Choose a patient ID from the dropdown.
 • Select two scan types (A = reference, B = evaluation).
 • Use 'View Slices' to visually verify CT + dose overlay.
 • Adjust gamma criteria if needed (e.g. 2% / 2 mm / 10% threshold).
 • Click 'Calculate Gamma Index' to run the 3D gamma analysis.

3) Outputs
 • Flattened gamma values (.csv)
 • Text summary report (.txt) with pass rate and criteria
 • Debug log file (.log) showing processing steps
