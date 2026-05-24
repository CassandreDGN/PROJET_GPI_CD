# Projet GPI - Cassandre DEGEN M1BBS

**This is my M1's GPI project.**
---
```   _______ .______    __     .______   .______        ______          __   _______   ______ .___________.
 /  _____||   _  \  |  |    |   _  \  |   _  \      /  __  \        |  | |   ____| /      ||           |
|  |  __  |  |_)  | |  |    |  |_)  | |  |_)  |    |  |  |  |       |  | |  |__   |  ,----'`---|  |----`
|  | |_ | |   ___/  |  |    |   ___/  |      /     |  |  |  | .--.  |  | |   __|  |  |         |  |     
|  |__| | |  |      |  |    |  |      |  |\  \----.|  `--'  | |  `--'  | |  |____ |  `----.    |  |     
 \______| | _|      |__|    | _|      | _| `._____| \______/   \______/  |_______| \______|    |__|   
  ``` 

 # Cryo-EM Particle Picker — Cassandre DEGEN M1BBS

A Python script that automatically detects particles in a cryo-EM micrograph (`.mrc` file) using a reference protein structure (`.pdb` file) as a template. It detects two orientations — **top views** (the protein seen from above, circular shape) and **side views** (the protein seen from the side, rectangular shape) — and saves the results as images.

---

## What you need before running

- Python 3.9 or newer
- The three files below all in the **same folder** on your computer

---

## Setup

**1. Download the script and requirements**

Download these two files and put them in a folder on your computer:

- [`particle_picker.py`](https://raw.githubusercontent.com/CassandreDGN/PROJET_GPI_CD/master/particle_picker.py) — right click → Save link as
- [`requirements.txt`](https://raw.githubusercontent.com/CassandreDGN/PROJET_GPI_CD/master/requirements.txt) — right click → Save link as

**2. Download the PDB structure**

- [`6BDF.pdb`](https://raw.githubusercontent.com/CassandreDGN/PROJET_GPI_CD/master/6BDF.pdb) — right click → Save link as

Put it in the **same folder** as the script.

**3. Download the MRC micrograph**

The `.mrc` file is too large to host on GitHub. Download it here:

👉 [Download map.mrc from Google Drive](https://drive.google.com/uc?export=download&id=1Qj30jSXcHEpkzE04cisbP6ljtnQ2Ausr)

Save it as `map.mrc` in the **same folder** as the script.

Your folder should look like this before running:
```
your-folder/
├── particle_picker.py
├── requirements.txt
├── map.mrc
└── 6BDF.pdb
```

**4. Install the dependencies**

Open a terminal in your folder and run:
```bash
pip install -r requirements.txt
```
---

## How to run

```bash
python particle_picker.py --mrc map.mrc --pdb 6BDF.pdb
```

---

## Outputs

All output files are saved in the same folder as your `.mrc` file by default.

| Output | Description |
|---|---|
| `particle_map.png` | The full micrograph with detected particles boxed in cyan (top views) and lime green (side views) |
| `gallery_top/` | Folder containing one cropped image per top-view particle |
| `gallery_top_overview.png` | A single image showing all top-view particles in a grid |
| `gallery_side/` | Folder containing one cropped image per side-view particle |
| `gallery_side_overview.png` | A single image showing all side-view particles in a grid |

To save outputs to a different folder, use `--outdir`:
```bash
python particle_picker.py --mrc map.mrc --pdb 6BDF.pdb --outdir results/
```

---

## Optional parameters

You don't need to touch these to get started, but they let you fine-tune the detection if results aren't great.

| Parameter | Default | Description |
|---|---|---|
| `--bin` | `5` | Downscaling factor applied to the image before processing (higher = faster but less precise) |
| `--thresh-top` | `0.23` | Detection sensitivity for top views — lower catches more particles but also more false positives |
| `--thresh-side` | `0.12` | Same for side views |
| `--min-dist` | `20` | Minimum distance in pixels between two detected particles |
| `--border` | `70` | How many pixels to ignore around the image edges |
| `--max-gallery` | `64` | Maximum number of particles shown in the gallery overview images |

**Example with custom parameters:**
```bash
python particle_picker.py --mrc map.mrc --pdb 6BDF.pdb --thresh-top 0.30 --bin 4
```

---

## Troubleshooting

**`ModuleNotFoundError`** — a dependency is missing. Run `pip install -r requirements.txt` again.

**`file not found`** — check that your `.mrc` and `.pdb` files are in the same folder as the script and that the names match exactly.

**Too many / too few particles detected** — adjust `--thresh-top` and `--thresh-side`. Increase the threshold to detect fewer (more strict), decrease it to detect more (more permissive).
