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

  --- 

#Particle Picker

A Python script that automatically detects particles in a cryo-EM micrograph (`.mrc` file) using a reference protein structure (`.pdb` file) as a template. It detects two orientations : **top views** (the protein seen from above, circular shape) and **side views** (the protein seen from the side, rectangular shape) — and saves the results as images.
---
## What you need before running

- Python 3.9 or newer
- The two input files in the same folder:
  - your micrograph file (`.mrc`)
  - your protein structure file (`.pdb`)

---

## Installation

**1. Clone or download this repository**

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

**2. Install the dependencies**

```bash
pip install -r requirements.txt
```
---

## How to run

```bash
python particle_picker.py --mrc your_map.mrc --pdb your_structure.pdb
```

Replace `your_map.mrc` and `your_structure.pdb` with the actual filenames.

**Example for the files provided you should use :**
```bash
python particle_picker.py --mrc map.mrc --pdb 6BDF.pdb
```

---

## Outputs

All output files are saved in the **same folder as your `.mrc` file** by default.

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

You don't need to change any of these to get started, especially since its just a M1 project but nonetheless  they let you fine-tune the detection if results aren't great and change settings as you please (i did not try to change them that way though).

| Parameter | Default | Description |
|---|---|---|
| `--bin` | `5` | Downscaling factor applied to the image before processing (higher = faster but less precise) |
| `--thresh-top` | `0.23` | Detection sensitivity for top views — lower catches more particles but also more false positives |
| `--thresh-side` | `0.12` | Same for side views |
| `--min-dist` | `20` | Minimum distance in pixels between two detected particles |
| `--border` | `70` | How many pixels to ignore around the image edges |
| `--max-gallery` | `64` | Maximum number of particles shown in the gallery overview images |

---

## Troubleshooting

**`ModuleNotFoundError`** — a dependency is missing. Run `pip install -r requirements.txt` again.

**`file not found`** — check that the path to your `.mrc` or `.pdb` file is correct, and that both files are in the same folder.

**Too many / too few particles detected** — adjust `--thresh-top` and `--thresh-side`. Increase the threshold to detect fewer (more strict), decrease it to detect more (more permissive).

