#!/usr/bin/env python
# coding: utf-8

import os, sys, pydicom, pickle, time, logging, pymedphys, threading

import tkinter as tk
import tkinter.ttk as ttk

from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import simpledialog
# Global GUI variables
selected_index_var = None
selected_index_menu = None


# Global data holders
data_structure = {}
dose_cache = {}  # currently unused


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)




def calculate_gamma_index(pid, scanA, scanB, gamma_inputs):
    logger.info("Gamma calculation started")
    entryA = data_structure[pid][scanA]
    entryB = data_structure[pid][scanB]

    doseA, xA, yA, zA = entryA["dose"]
    doseB, xB, yB, zB = entryB["dose"]

    # --- Sanitize arrays before gamma ---
    def ensure_numpy(arr, name="array"):
        arr_np = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
        logger.debug(
            f"{name}: shape={arr_np.shape}, dtype={arr_np.dtype}, "
            f"contiguous={arr_np.flags['C_CONTIGUOUS']}"
        )
        return arr_np

    doseA = ensure_numpy(doseA, "doseA")
    doseB = ensure_numpy(doseB, "doseB")
    xA, yA, zA = (ensure_numpy(xA, "xA"),
                  ensure_numpy(yA, "yA"),
                  ensure_numpy(zA, "zA"))
    xB, yB, zB = (ensure_numpy(xB, "xB"),
                  ensure_numpy(yB, "yB"),
                  ensure_numpy(zB, "zB"))

    # Geometry validation
    dx = np.max(np.abs(xA - xB))
    dy = np.max(np.abs(yA - yB))
    dz = np.max(np.abs(zA - zB))

    logger.debug(f"Geometry deltas: dX={dx:.6f}, dY={dy:.6f}, dZ={dz:.6f}")

    if dx > 1e-6 or dy > 1e-6 or dz > 1e-6:
        return {
            "error": f"Geometry Mismatch:\n"
                     f"dX: {dx:.3f} mm | dY: {dy:.3f} mm | dZ: {dz:.3f} mm"
        }

    # Gamma Criteria from UI
    dose_percent  = gamma_inputs["dose_percent"]
    dta_mm        = gamma_inputs["dta_mm"]
    dose_threshold= gamma_inputs["dose_threshold"]
    local_flag    = gamma_inputs["local_flag"]

    gamma_options = {
        "dose_percent_threshold": dose_percent,
        "distance_mm_threshold":  dta_mm,
        "lower_percent_dose_cutoff": dose_threshold,
        "local_gamma": local_flag,
        "max_gamma": 2.0,
        "interp_fraction": 10,
    }

    for name, arr in [("xA", xA), ("yA", yA), ("zA", zA), ("doseA", doseA),
                  ("xB", xB), ("yB", yB), ("zB", zB), ("doseB", doseB)]:
        logger.debug(f"{name}: type={type(arr)}, dtype={arr.dtype}, contiguous={arr.flags['C_CONTIGUOUS']}")

    
    logger.debug(f"xA range: {xA[0]} → {xA[-1]}, spacing={np.mean(np.diff(xA))}")
    logger.debug(f"xB range: {xB[0]} → {xB[-1]}, spacing={np.mean(np.diff(xB))}")
    logger.debug(f"yA range: {yA[0]} → {yA[-1]}, spacing={np.mean(np.diff(yA))}")
    logger.debug(f"yB range: {yB[0]} → {yB[-1]}, spacing={np.mean(np.diff(yB))}")
    logger.debug(f"zA range: {zA[0]} → {zA[-1]}, spacing={np.mean(np.diff(zA))}")
    logger.debug(f"zB range: {zB[0]} → {zB[-1]}, spacing={np.mean(np.diff(zB))}")
    
    logger.debug(f"doseA max={np.max(doseA)}, min={np.min(doseA)}")
    logger.debug(f"doseB max={np.max(doseB)}, min={np.min(doseB)}")
    logger.debug(f"Expected dose array shape: ({len(xA)}, {len(yA)}, {len(zA)})")
    logger.debug(f"Actual doseA shape: {doseA.shape}")
    logger.debug(f"Actual doseB shape: {doseB.shape}")
    logger.debug(f"Types: xA={type(xA)}, yA={type(yA)}, zA={type(zA)}, doseA={type(doseA)}")
    
    logger.debug(f"Gamma options: {gamma_options}")
    
    logger.info("Geometry check passed")
    try:
        gamma = pymedphys.gamma(
            (xA.astype(np.float64), yA.astype(np.float64), zA.astype(np.float64)), doseA.astype(np.float64),
            (xB.astype(np.float64), yB.astype(np.float64), zB.astype(np.float64)), doseB.astype(np.float64),
            **gamma_options
        )
        logger.debug(f"Gamma result stats: min={np.nanmin(gamma)}, max={np.nanmax(gamma)}")
    except Exception as e:
        logger.error(f"Gamma calculation failed: {e}")
        return {"error": str(e)}
        
    logger.info("Gamma calculation completed successfully")

    valid_voxels = ~np.isnan(gamma)
    pass_mask = (gamma[valid_voxels] <= 1)
    pass_rate = 100 * np.sum(pass_mask) / np.sum(valid_voxels)

    logger.info(f"Gamma pass rate: {pass_rate:.2f}% "
                f"(valid voxels={np.sum(valid_voxels)})")

    # File save base name
    saveBase = f"gamma_{pid}_{scanA}_vs_{scanB}_{int(dose_percent*100)}pc_{int(dta_mm)}mm"
    paths = {
        "npy": os.path.join(root_path, saveBase + ".npy"),
        "csv": os.path.join(root_path, saveBase + ".csv"),
        "txt": os.path.join(root_path, saveBase + "_summary.txt"),
        "png": os.path.join(root_path, saveBase + ".png"),
    }

    # Save gamma map
    np.save(paths["npy"], gamma)
    np.savetxt(paths["csv"],
               gamma.reshape(-1, gamma.shape[2]), delimiter=",")

    # Summary text file
    with open(paths["txt"], "w") as f:
        f.write(f"Gamma Report\nPID: {pid}\n{scanA} vs {scanB}\n\n")
        f.write(f"Criteria: {dose_percent*100:.1f}% / {dta_mm:.1f} mm\n")
        f.write(f"Threshold: {dose_threshold*100:.1f}%\n")
        f.write(f"Local Gamma: {local_flag}\n\n")
        f.write(f"Pass Rate: {pass_rate:.2f}%\n")
        f.write(f"Valid Voxels: {np.sum(valid_voxels)}\n")

    # Save slice image
    mid_idx = gamma.shape[2] // 2
    fig = plt.figure(figsize=(7,7))
    plt.imshow(gamma[:,:,mid_idx], cmap="coolwarm", vmin=0, vmax=2)
    plt.colorbar(label="Gamma")
    plt.title(f"Mid Slice Gamma — {pid}\n{scanA} vs {scanB}")
    plt.savefig(paths["png"])
    plt.close(fig)

    
    # Debug BEFORE shutdown and return
    logger.debug(f"Gamma calculation finished, result shape={gamma.shape}")

    # --- Cleanup logging handlers to release file lock ---
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)
    logging.shutdown()   # <-- ensures Windows releases the file
    
    return {
        "pass_rate": pass_rate,
        "paths": paths
    }



def load_patient_data():
    global root_path, data_structure, selected_index_var
    
    root_path = filedialog.askdirectory(title="Select Patient Root Folder")
    if not root_path:
        messagebox.showwarning("No Folder", "Please select a valid folder.")
        return
    
    # 👉 Initialize logger here, pointing to root_path
    log_file = os.path.join(root_path, "gamma_debug.log")
    if logger.hasHandlers():
        logger.handlers.clear()
    fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info(f"Logger initialized for patient root: {root_path}")    
    logger.info("Starting patient data load...")

    # existing code...
    if not root_path:
        logger.warning("No folder selected for patient data")
        return

    # after successful load
    logger.info(f"Patient data loaded successfully from {root_path}")

    data_structure = {}

    

    patient_ids = sorted(
        d for d in os.listdir(root_path)
        if os.path.isdir(os.path.join(root_path, d))
    )

    if not patient_ids:
        messagebox.showerror("Error", "No patient folders found.")
        return

    missing = []

    for pid in patient_ids:
        patient_path = os.path.join(root_path, pid)
        scan_types = sorted(
            d for d in os.listdir(patient_path)
            if os.path.isdir(os.path.join(patient_path, d))
        )

        data_structure[pid] = {}

        for scan in scan_types:
            scan_path = os.path.join(patient_path, scan)

            ct_slices, dose_files = find_dicom_files(scan_path)

            if not ct_slices:
                missing.append((pid, scan, "NO CT FOUND"))
                continue

            if not dose_files:
                missing.append((pid, scan, "NO RTDOSE FOUND"))
                continue

            dose_file = dose_files[-1]  # simple + safe

            ct_3d, x_ct, y_ct, z_ct = load_ct_data(ct_slices)
            dose_3d, x_ds, y_ds, z_ds = load_dose_data(dose_file)

            data_structure[pid][scan] = {
                "ct":   (ct_3d, x_ct, y_ct, z_ct),
                "dose": (dose_3d, x_ds, y_ds, z_ds)
            }
    
    # Select the first patient
    selected_index_menu.configure(values=patient_ids)
    selected_index_var.set(patient_ids[0])
    _update_scan_type_dropdowns()



    messagebox.showinfo(
        "Import Complete",
        f"{len(patient_ids)} Patients Loaded\n"
        f"Missing entries: {len(missing)}\n"
        f"Check console for detail"
    )

    if missing:
        print("\n--- Missing Data Report ---")
        for pid, scan, issue in missing:
            print(f" PID: {pid:<6} | Scan: {scan:<10} | {issue}")
    #Only touch root if it exists (running inside full GUI)
    if "root" in globals():
        root.focus_force()
        root.after(10, lambda: root.attributes('-topmost', False))
    

def load_dose_data(dose_file):
    logger.info(f"Loading RTDOSE file: {dose_file}")

    """Load RTDOSE DICOM file into 3D dose + coordinate arrays"""
        
    ds = pydicom.dcmread(dose_file)

    # Extract dose matrix and scaling
    dose_3d = ds.pixel_array * ds.DoseGridScaling

    # Orientation fix: RTDose stored as (z,y,x), convert to (x,y,z)
    dose_3d = np.transpose(dose_3d, (2, 1, 0))

    # Coordinates
    pixel_spacing = np.array(ds.PixelSpacing)  # [row, col] = [y,x]
    image_pos = np.array(ds.ImagePositionPatient)
    z_positions = np.array(ds.GridFrameOffsetVector) + image_pos[2]

    x_coords = np.linspace(
        image_pos[0],
        image_pos[0] + pixel_spacing[1]*(dose_3d.shape[0]-1),
        dose_3d.shape[0]
    )
    
    y_coords = np.linspace(
        image_pos[1],
        image_pos[1] + pixel_spacing[0]*(dose_3d.shape[1]-1),
        dose_3d.shape[1]
    )
    
    logger.info("RTDOSE file loaded and dose matrix extracted")

    return dose_3d, x_coords, y_coords, z_positions

def find_dicom_files(scan_path):
    """Walk all subfolders to find CT slices and RTDOSE file properly."""
    ct_slices = []
    dose_files = []

    for root, _, files in os.walk(scan_path):
        for fname in files:
            if not fname.lower().endswith('.dcm'):
                continue

            f = os.path.join(root, fname)

            try:
                ds = pydicom.dcmread(f, stop_before_pixels=True)
                sop = ds.SOPClassUID

                if sop == pydicom.uid.CTImageStorage:
                    ct_slices.append(f)

                elif sop == pydicom.uid.RTDoseStorage:
                    dose_files.append(f)

            except Exception:
                continue

    return ct_slices, dose_files


def _update_scan_type_dropdowns(*args):
    """Refresh scan selectors when patient selection changes"""
    pid = selected_index_var.get()
    if pid not in data_structure:
        return

    available_scans = list(data_structure[pid].keys())
    if not available_scans:
        return

    scan_type_A_menu['values'] = available_scans
    scan_type_B_menu['values'] = available_scans

    scan_type_A_var.set(available_scans[0])
    scan_type_B_var.set(available_scans[-1])



def load_ct_data(ct_files):
    slices = []
    positions = []

    for f in ct_files:
        ds = pydicom.dcmread(f)
        slices.append(ds.pixel_array)
        positions.append(ds.ImagePositionPatient[2])

    order = np.argsort(positions)
    ct_3d = np.stack([slices[i] for i in order])

    # Extract spacing + coords
    ds0 = pydicom.dcmread(ct_files[0])
    px, py = ds0.PixelSpacing
    x0, y0, z0 = ds0.ImagePositionPatient

    z = np.array(sorted(positions))
    x = np.linspace(x0, x0 + px*(ct_3d.shape[2]-1), ct_3d.shape[2])
    y = np.linspace(y0, y0 + py*(ct_3d.shape[1]-1), ct_3d.shape[1])

    return ct_3d, x, y, z

def _get_scan_data(pid, scan_type):
    entry = data_structure[pid][scan_type]
    return entry["ct"], entry["dose"]
    

def view_slices():
    """Split view: Scan A and Scan B dose overlay on CT for sanity check."""
    pid = selected_index_var.get()
    scanA = scan_type_A_var.get()
    scanB = scan_type_B_var.get()

    if pid not in data_structure:
        messagebox.showerror("Selection Error", "Please select a valid patient.")
        return

    if scanA not in data_structure[pid] or scanB not in data_structure[pid]:
        messagebox.showerror("Selection Error",
                             "Selected scan types not available for this patient.")
        return

    # Load CT + dose from memory (NO MORE DISK LOAD)
    (ctA, xA, yA, zA) = data_structure[pid][scanA]["ct"]
    (doseA, dxA, dyA, dzA) = data_structure[pid][scanA]["dose"]

    (ctB, xB, yB, zB) = data_structure[pid][scanB]["ct"]
    (doseB, dxB, dyB, dzB) = data_structure[pid][scanB]["dose"]

    # Check grid match for dose overlay
    if doseA.shape != doseB.shape:
        messagebox.showerror("Geometry Error",
                             f"Dose grid mismatch:\n{scanA}: {doseA.shape}\n{scanB}: {doseB.shape}")
        return

    nz = doseA.shape[2]
    maxdose = max(np.max(doseA), np.max(doseB))
    doseA_norm = doseA / maxdose
    doseB_norm = doseB / maxdose

    viewer = tk.Toplevel()
    viewer.title(f"Viewer — {pid}: {scanA} vs {scanB}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    fig.tight_layout()
    canvas = FigureCanvasTkAgg(fig, viewer)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    plt.close(fig)  
    plt.ioff()  # prevent Jupyter from auto-rendering


    slice_idx = tk.IntVar(value=nz // 2)

    def update_slice(*_):
        idx = slice_idx.get()
        idx = max(0, min(idx, nz-1))

        # Correct orientation — CT[z,y,x], Dose[x,y,z]
        ctA_img = np.rot90(ctA[idx, :, :], 1)
        ctB_img = np.rot90(ctB[idx, :, :], 1)

        doseA_img = np.rot90(doseA_norm[:, :, idx], 1)
        doseB_img = np.rot90(doseB_norm[:, :, idx], 1)

        # Normalize CT for display
        def normalize(img):
            return (img - img.min()) / (img.max() - img.min() + 1e-6)

        ax[0].clear()
        ax[0].imshow(normalize(ctA_img), cmap='gray')
        ax[0].imshow(doseA_img, cmap='jet', alpha=(doseA_img>0.05)*0.5)
        ax[0].set_title(f"{scanA} — Slice {idx+1}/{nz}")
        ax[0].axis('off')

        ax[1].clear()
        ax[1].imshow(normalize(ctB_img), cmap='gray')
        ax[1].imshow(doseB_img, cmap='jet', alpha=(doseB_img>0.05)*0.5)
        ax[1].set_title(f"{scanB} — Slice {idx+1}/{nz}")
        ax[1].axis('off')

        fig.suptitle(f"Patient: {pid}", fontsize=14)
        canvas.draw()

    slider = tk.Scale(viewer, from_=0, to=nz-1,
                      orient=tk.HORIZONTAL, variable=slice_idx,
                      command=update_slice, length=500)
    slider.pack(pady=10)

    update_slice()
    if "root" in globals():
        root.focus_force()
        root.after(10, lambda: root.attributes('-topmost', False))

def start_gamma_thread():
    pid = selected_index_var.get()
    scanA = scan_type_A_var.get()
    scanB = scan_type_B_var.get()

    if pid == "" or scanA == "" or scanB == "":
        messagebox.showwarning("Missing Data", "Load data and select patient/scans first")
        return
    
    if scanA == scanB:
        messagebox.showerror("Error", "Choose two different scan types")
        return

    # Get gamma settings safely in main thread
    try:
        gamma_inputs = {
            "dose_percent": float(dose_diff_entry.get()) / 100.0,
            "dta_mm": float(dta_entry.get()),
            "dose_threshold": float(threshold_entry.get()) / 100.0,
            "local_flag": (gamma_type_var.get() == "Local")
        }
    except ValueError:
        messagebox.showerror("Error", "Invalid gamma settings input")
        return

    # Disable button during compute
    calc_btn.config(state="disabled")

    # Start threaded computation
    threading.Thread(
        target=calculate_gamma_index_threaded,
        args=(pid, scanA, scanB, gamma_inputs)
    ).start()


def calculate_gamma_index_threaded(pid, scanA, scanB, gamma_inputs):
    try:
        result = calculate_gamma_index(pid, scanA, scanB, gamma_inputs)

        def handle_result():
            if result is None:
                messagebox.showerror("Gamma Error", "Unknown error occurred.")
                return
            
            if "error" in result:
                messagebox.showerror("Geometry Error", result["error"])
                return
            
            # Success — update UI
            pass_rate = result["pass_rate"]
            paths = result["paths"]

            messagebox.showinfo(
                "Gamma Complete",
                f"Pass Rate: {pass_rate:.2f}%\n\n"
                "Saved:\n" +
                "\n".join(paths.values())
            )

        root.after(0, handle_result)

    except Exception as e:
        root.after(0, lambda err=e: messagebox.showerror("Gamma Error", str(err)))
    finally:
        root.after(0, finish_gamma_ui_updates)



def finish_gamma_ui_updates():
    if "calc_btn" in globals():
        calc_btn.config(state="normal")

    if "root" in globals():
        try:
            root.focus_force()
            root.lift()
            root.focus_force()
            root.after(10, lambda: root.attributes('-topmost', False))
        except tk.TclError:
            # In case Tk is already destroyed or not in a GUI context
            pass

root = tk.Tk()
root.withdraw()  # Hide until UI builds completely
root.title("GammaAnalyzer – Gamma Index GUI")
root.geometry("600x600")
root.attributes('-topmost', True)

# --- Title ---
ttk.Label(root, text="GammaAnalyzer – Gamma Index GUI",
          font=("Arial", 16, "bold")).grid(
    row=0, column=0, columnspan=2, pady=(10, 20)
)

# --- Section: Patient Data ---
ttk.Label(root, text="Patient Data", font=("Arial", 14, "bold")).grid(
    row=1, column=0, sticky="w", padx=10
)
ttk.Button(root, text="Load Patient Data",
           command=load_patient_data).grid(
    row=2, column=0, padx=10, pady=5, sticky="ew"
)
ttk.Button(root, text="Pre-Processed (disabled)",
           state="disabled").grid(
    row=2, column=1, padx=10, pady=5, sticky="ew"
)

ttk.Separator(root, orient="horizontal").grid(
    row=3, column=0, columnspan=2, sticky="ew", pady=10
)

# --- Section: Gamma Settings ---
ttk.Label(root, text="Gamma Settings", font=("Arial", 14, "bold")).grid(
    row=7, column=0, sticky="w", padx=10
)

ttk.Label(root, text="Dose Difference (%):").grid(
    row=8, column=0, sticky="w", padx=10
)
dose_diff_entry = ttk.Entry(root)
dose_diff_entry.insert(0, "2")  # default: 2%
dose_diff_entry.grid(row=8, column=1, padx=10, pady=2, sticky="ew")

ttk.Label(root, text="DTA (mm):").grid(
    row=9, column=0, sticky="w", padx=10
)
dta_entry = ttk.Entry(root)
dta_entry.insert(0, "2")  # default: 2 mm
dta_entry.grid(row=9, column=1, padx=10, pady=2, sticky="ew")

ttk.Label(root, text="Dose Threshold (%):").grid(
    row=10, column=0, sticky="w", padx=10
)
threshold_entry = ttk.Entry(root)
threshold_entry.insert(0, "10")  # default: 10%
threshold_entry.grid(row=10, column=1, padx=10, pady=2, sticky="ew")

ttk.Label(root, text="Gamma Type:").grid(
    row=11, column=0, sticky="w", padx=10
)
gamma_type_var = tk.StringVar(value="Global")
gamma_type_menu = ttk.Combobox(
    root, textvariable=gamma_type_var,
    values=["Global", "Local"], state="readonly"
)
gamma_type_menu.grid(row=11, column=1, padx=10, pady=2, sticky="ew")

# --- Section: Gamma Comparison Setup ---
ttk.Label(root, text="Gamma Comparison Setup",
          font=("Arial", 14, "bold")).grid(
    row=12, column=0, sticky="w", padx=10, pady=(10, 5)
)

ttk.Label(root, text="Select Patient:").grid(
    row=13, column=0, sticky="w", padx=10
)
selected_index_var = tk.StringVar(value='')
selected_index_menu = ttk.Combobox(
    root, textvariable=selected_index_var, state="readonly"
)
selected_index_menu.grid(row=13, column=1, padx=10, pady=2, sticky="ew")
selected_index_menu.bind("<<ComboboxSelected>>", lambda *_: _update_scan_type_dropdowns())

# Scan Type A Selector
ttk.Label(root, text="Scan Type A (Reference):").grid(
    row=14, column=0, sticky="w", padx=10
)
scan_type_A_var = tk.StringVar(value='')
scan_type_A_menu = ttk.Combobox(
    root, textvariable=scan_type_A_var, state="readonly"
)
scan_type_A_menu.grid(row=14, column=1, padx=10, pady=2, sticky="ew")

# Scan Type B Selector
ttk.Label(root, text="Scan Type B (Evaluation):").grid(
    row=15, column=0, sticky="w", padx=10
)
scan_type_B_var = tk.StringVar(value='')
scan_type_B_menu = ttk.Combobox(
    root, textvariable=scan_type_B_var, state="readonly"
)
scan_type_B_menu.grid(row=15, column=1, padx=10, pady=2, sticky="ew")

# --- Visualization & QA Check ---
vis_row = 16
ttk.Separator(root, orient="horizontal").grid(
    row=vis_row, column=0, columnspan=2, sticky="ew", pady=10
)
ttk.Label(root, text="Visualization & QA Check",
          font=("Arial", 14, "bold")).grid(
    row=vis_row + 1, column=0, sticky="w", padx=10
)

ttk.Button(root, text="View Slices",
           command=view_slices).grid(
    row=vis_row + 2, column=0, columnspan=2,
    padx=10, pady=5, sticky="ew"
)

# --- Gamma Calculation ---
gamma_row = vis_row + 3
ttk.Separator(root, orient="horizontal").grid(
    row=gamma_row, column=0, columnspan=2, sticky="ew", pady=10
)
ttk.Label(root, text="Gamma Calculation",
          font=("Arial", 14, "bold")).grid(
    row=gamma_row + 1, column=0, sticky="w", padx=10
)

calc_btn = ttk.Button(
    root, text="Calculate Gamma Index",
    command=start_gamma_thread
)
calc_btn.grid(
    row=gamma_row + 2, column=0, columnspan=2,
    padx=10, pady=10, sticky="ew"
)

# Make columns expand evenly
root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=1)

def on_close():
    root.quit()
    root.destroy()

root.protocol("WM_DELETE_WINDOW", on_close)


if __name__ == "__main__":
    root.deiconify()  # show the window now that it's built
    root.mainloop()