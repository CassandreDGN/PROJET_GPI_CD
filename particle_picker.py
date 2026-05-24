#!/usr/bin/env python3
"""
particle_picker.py — Cryo-EM particle picker (GPI project)

Usage:
    python particle_picker.py --mrc map.mrc --pdb 6BDF.pdb [options]

Inputs:
    --mrc   Path to the .mrc micrograph file
    --pdb   Path to the .pdb structure file

Outputs (written next to the input files by default, override with --outdir):
    particle_map.png          — micrograph with top/side bounding boxes
    gallery_top/              — individual top-view particle crops (.png)
    gallery_side/             — individual side-view particle crops (.png)
    gallery_top_overview.png  — full top-view gallery montage
    gallery_side_overview.png — full side-view gallery montage
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")          # telling plt to not try and open the outputs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import cv2
import mrcfile
from scipy import ndimage
from scipy.ndimage import rotate, maximum_filter
from Bio.PDB import PDBParser

np.random.seed(42)

# Arguments that can be modified

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cryo-EM particle picker: takes an MRC micrograph and a "
                    "PDB reference structure, detects top-view (circle) and "
                    "side-view (rectangle) particles, and saves the results."
    )
    parser.add_argument("--mrc", required=True,
                        help="Path to the input .mrc micrograph file")
    parser.add_argument("--pdb", required=True,
                        help="Path to the input .pdb structure file")
    parser.add_argument("--outdir", default=None,
                        help="Output directory (default: same folder as --mrc)")
    parser.add_argument("--bin", type=int, default=5,
                        help="Binning factor applied to the micrograph (default: 5)")
    parser.add_argument("--thresh-top", type=float, default=0.23,
                        help="NCC threshold for top-view detection (default: 0.23)")
    parser.add_argument("--thresh-side", type=float, default=0.12,
                        help="NCC threshold for side-view detection (default: 0.12)")
    parser.add_argument("--min-dist", type=int, default=20,
                        help="Minimum pixel distance between picked particles (default: 20)")
    parser.add_argument("--border", type=int, default=70,
                        help="Border width (px) to ignore around the image edge (default: 70)")
    parser.add_argument("--max-gallery", type=int, default=64,
                        help="Maximum number of particles shown in the gallery montage (default: 64)")
    return parser.parse_args()


#ETAPE 1 : LOADING THE MRC

def load_mrc(path):
    print(f"\n[1/6] Loading MRC file: {path}")
    if not os.path.isfile(path):
        sys.exit(f"  ERROR: file not found — {path}")

    with mrcfile.open(path, permissive=True) as mrc:
        if mrc.data is None:
            sys.exit("  ERROR: MRC file contains no readable data.")
        data = mrc.data.copy()

    print(f"  Shape : {data.shape}")
    print(f"  Dtype : {data.dtype}")
    print(f"  Range : [{data.min():.4f}, {data.max():.4f}]")

    if data.ndim == 3:
        print("  3-D volume detected — using XY (sum-Z) projection as micrograph.")
        data = np.sum(data, axis=0).astype(np.float32)
    elif data.ndim == 2:
        print("  2-D micrograph detected.")
    else:
        sys.exit(f"  ERROR: unexpected MRC dimensionality ({data.ndim}D).")

    return data.astype(np.float32)


# STEP 2 = Pre-process image (micrograph)


def preprocess(data, bin_factor=5, pad_width=70):
    print(f"\n[2/6] Pre-processing (median filter → bin ×{bin_factor} → normalise → invert → pad)")

    # Median filter
    filtered = ndimage.median_filter(data, size=3)
    filtered = np.clip(filtered, 0, 255).astype(np.uint8)

    # Binning
    h, w = filtered.shape
    new_h, new_w = h // bin_factor, w // bin_factor
    binned = (filtered[:new_h * bin_factor, :new_w * bin_factor]
              .reshape(new_h, bin_factor, new_w, bin_factor)
              .mean(axis=(1, 3))
              .astype(np.uint8))

    # Normalise
    normalised = (binned - binned.mean()) / (binned.std() + 1e-8)

    # Invert (particles are bright in raw EM → dark after inversion for NCC)
    inverted = normalised * -1.0

    # Edge-pad
    padded = np.pad(inverted,
                    ((pad_width, pad_width), (pad_width, pad_width)),
                    mode="edge")

    print(f"  Final padded shape: {padded.shape}")
    return padded.astype(np.float32)


# STEP 3 = Load PDB and compute 2D projections

def load_pdb_projections(path):
    print(f"\n[3/6] Loading PDB file: {path}")
    if not os.path.isfile(path):
        sys.exit(f"  ERROR: file not found — {path}")

    parser = PDBParser(QUIET=True)
    pdb_id = os.path.splitext(os.path.basename(path))[0]
    structure = parser.get_structure(pdb_id, path)

    model = structure[0]
    all_coords = []
    for chain in model.get_chains():
        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                all_coords.append(atom.get_coord())

    coords = np.array(all_coords)
    print(f"  Total atoms: {len(coords)}")
    print(f"  X range: [{coords[:,0].min():.1f}, {coords[:,0].max():.1f}] Å")
    print(f"  Y range: [{coords[:,1].min():.1f}, {coords[:,1].max():.1f}] Å")
    print(f"  Z range: [{coords[:,2].min():.1f}, {coords[:,2].max():.1f}] Å")

    coords_c = coords - coords.mean(axis=0)
    padding_ang = 5.0  

    total_span = coords_c.max() - coords_c.min()
    ang_per_px = total_span / 48.0

    def make_proj(c1, c2):
        mn1, mx1 = c1.min() - padding_ang, c1.max() + padding_ang
        mn2, mx2 = c2.min() - padding_ang, c2.max() + padding_ang
        span = max(mx1 - mn1, mx2 - mn2)
        bins = max(8, int(np.ceil(span / ang_per_px)))
        H, _, _ = np.histogram2d(c1, c2, bins=bins,
                                 range=[[mn1, mx1], [mn2, mx2]])
        return (H / H.max()).astype(np.float32)

    proj_xy = make_proj(coords_c[:, 0], coords_c[:, 1])  # top view
    proj_xz = make_proj(coords_c[:, 0], coords_c[:, 2])  # side view

    print(f"  Top-view template (XY): {proj_xy.shape}")
    print(f"  Side-view template (XZ): {proj_xz.shape}")
    return proj_xy, proj_xz


# STEP 4 — Template matching helpers

def _ncc(image, template):
    """Normalised cross-correlation"""
    img_h, img_w = image.shape
    tmp_h, tmp_w = template.shape
    out_h = img_h - tmp_h + 1
    out_w = img_w - tmp_w + 1
    output = np.zeros((out_h, out_w), dtype=np.float32)

    tmpl_norm = template - template.mean()
    tmpl_std = np.sqrt((tmpl_norm ** 2).sum())
    if tmpl_std == 0:
        return output

    for i in range(out_h):
        for j in range(out_w):
            region = image[i:i + tmp_h, j:j + tmp_w]
            reg_norm = region - region.mean()
            reg_std = np.sqrt((reg_norm ** 2).sum())
            if reg_std == 0:
                output[i, j] = 0.0
            else:
                output[i, j] = (reg_norm * tmpl_norm).sum() / (reg_std * tmpl_std)
    return output


def template_matching(image, template, method="ncc"):
    """Run template matching with half-template padding so the result is
    the same spatial size as the input image."""
    tmp_h, tmp_w = template.shape
    pad_h, pad_w = tmp_h // 2, tmp_w // 2
    padded = np.pad(image,
                    ((pad_h, pad_h), (pad_w, pad_w)),
                    mode="constant", constant_values=0)
    if method == "ncc":
        result = _ncc(padded, template)
    else:
        raise ValueError(f"Unknown method: {method}")
    return result[:image.shape[0], :image.shape[1]]


def extract_objects(response_map, threshold=0.7, min_distance=15):
    """Peak detection: threshold + non-maximum suppression."""
    above = response_map >= threshold
    local_max = maximum_filter(response_map, size=min_distance) == response_map
    peaks = above & local_max
    y_coords, x_coords = np.where(peaks)
    return list(zip(x_coords.tolist(), y_coords.tolist()))


# STEP 5 — Particle picking

def pick_particles(padded_image, proj_xy, proj_xz,
                   thresh_top=0.23, thresh_side=0.12,
                   min_dist=20, border=70):
    print(f"\n[4/6] Picking particles")

    box_size = max(max(proj_xy.shape), max(proj_xz.shape)) + 10
    half_box = box_size // 2
    img_h, img_w = padded_image.shape

    def _valid(pts):
        out = []
        for x, y in pts:
            if (y - half_box >= border and y + half_box < img_h - border and
                    x - half_box >= border and x + half_box < img_w - border):
                out.append((x, y))
        return out

    # --- Top views (circles) via NCC ---
    print("  Scanning for Top Views (Circles) via NCC…")
    hm_top = template_matching(padded_image, proj_xy, method="ncc")
    hm_top[:border, :] = 0;  hm_top[-border:, :] = 0
    hm_top[:, :border] = 0;  hm_top[:, -border:] = 0
    detected_tops = _valid(extract_objects(hm_top, threshold=thresh_top,
                                           min_distance=min_dist))
    print(f"  → {len(detected_tops)} valid top-view particles")

    # --- Side views (rectangles) via OpenCV multi-angle ---
    print("  Scanning for Side Views (Rectangles) via multi-angle OpenCV…")
    img_cv = padded_image.astype(np.float32)
    hm_side = np.zeros_like(padded_image, dtype=np.float32)

    for angle in np.arange(0, 360, 15):
        rot_tmpl = rotate(proj_xz, angle, reshape=True,
                          mode="constant", cval=0.0).astype(np.float32)
        mr = cv2.matchTemplate(img_cv, rot_tmpl, cv2.TM_CCOEFF_NORMED)
        dh = img_h - mr.shape[0]
        dw = img_w - mr.shape[1]
        pm = np.pad(mr,
                    ((dh // 2, dh - dh // 2), (dw // 2, dw - dw // 2)),
                    mode="constant")
        hm_side = np.maximum(hm_side, pm[:img_h, :img_w])

    hm_side[:border, :] = 0;  hm_side[-border:, :] = 0
    hm_side[:, :border] = 0;  hm_side[:, -border:] = 0

    # Mutual exclusion: zero out already-picked top-view locations
    for (xc, yc) in detected_tops:
        cv2.circle(hm_side, (int(xc), int(yc)), int(min_dist), 0, -1)

    detected_sides = _valid(extract_objects(hm_side, threshold=thresh_side,
                                            min_distance=min_dist))
    print(f"  → {len(detected_sides)} valid side-view particles")

    # --- Crop galleries ---
    top_gallery  = [padded_image[y - half_box:y + half_box,
                                 x - half_box:x + half_box]
                    for x, y in detected_tops]
    side_gallery = [padded_image[y - half_box:y + half_box,
                                 x - half_box:x + half_box]
                    for x, y in detected_sides]

    return detected_tops, detected_sides, top_gallery, side_gallery, proj_xy, proj_xz


# STEP 6 — Save outputs

def save_particle_map(padded_image, detected_tops, detected_sides,
                      proj_xy, proj_xz, outdir):
    path = os.path.join(outdir, "particle_map.png")
    print(f"\n[5/6] Saving particle map → {path}")

    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax.imshow(padded_image, cmap="gray",
              vmin=np.percentile(padded_image, 5),
              vmax=np.percentile(padded_image, 95))
    ax.set_title(
        f"Micrograph Particle Map\n"
        f"Cyan: Top view ({len(detected_tops)}) | "
        f"Lime: Side view ({len(detected_sides)})",
        fontsize=14, fontweight="bold"
    )

    h_t, w_t = proj_xy.shape
    for x, y in detected_tops:
        ax.add_patch(patches.Rectangle(
            (x - w_t // 2, y - h_t // 2), w_t, h_t,
            linewidth=1.5, edgecolor="cyan", facecolor="none"
        ))

    h_s, w_s = proj_xz.shape
    for x, y in detected_sides:
        ax.add_patch(patches.Rectangle(
            (x - w_s // 2, y - h_s // 2), w_s, h_s,
            linewidth=1.5, edgecolor="lime", facecolor="none"
        ))

    ax.axis("off")
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved.")


def save_gallery(particle_list, label, outdir, num_cols=8, max_images=64):
    """Save each particle as an individual PNG, plus a montage overview."""
    if not particle_list:
        print(f"  No particles for '{label}' — skipping.")
        return

    folder = os.path.join(outdir, label)
    os.makedirs(folder, exist_ok=True)

    for i, crop in enumerate(particle_list):
        fname = os.path.join(folder, f"{label}_{i + 1:04d}.png")
        fig, ax = plt.subplots(1, 1, figsize=(2, 2))
        ax.imshow(crop, cmap="gray")
        ax.axis("off")
        fig.savefig(fname, dpi=72, bbox_inches="tight", pad_inches=0)
        plt.close(fig)

    # Overview montage
    num_show  = min(len(particle_list), max_images)
    num_rows  = (num_show + num_cols - 1) // num_cols
    fig, axes = plt.subplots(num_rows, num_cols,
                             figsize=(num_cols * 1.6, num_rows * 1.6))
    fig.suptitle(f"{label.replace('_', ' ').title()} ({num_show} shown)",
                 fontsize=14, fontweight="bold", y=1.01)

    for i, ax in enumerate(np.array(axes).flat):
        if i < num_show:
            ax.imshow(particle_list[i], cmap="gray")
        ax.axis("off")

    plt.tight_layout()
    overview_path = os.path.join(outdir, f"{label}_overview.png")
    fig.savefig(overview_path, dpi=120, bbox_inches="tight")
    plt.close(fig)

    print(f"  {len(particle_list)} crops  → {folder}/")
    print(f"  Gallery overview         → {overview_path}")



def main():
    args = parse_args()
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.mrc))
    os.makedirs(outdir, exist_ok=True)
    print(f"\nOutput directory: {outdir}")

    
    data          = load_mrc(args.mrc)
    padded_image  = preprocess(data, bin_factor=args.bin, pad_width=args.border)
    proj_xy, proj_xz = load_pdb_projections(args.pdb)

    (detected_tops, detected_sides,
     top_gallery, side_gallery,
     proj_xy, proj_xz) = pick_particles(
        padded_image, proj_xy, proj_xz,
        thresh_top=args.thresh_top,
        thresh_side=args.thresh_side,
        min_dist=args.min_dist,
        border=args.border,
    )

    print(f"\n[5/6] Saving outputs…")
    save_particle_map(padded_image, detected_tops, detected_sides,
                      proj_xy, proj_xz, outdir)

    print(f"\n[6/6] Saving galleries…")
    save_gallery(top_gallery,  "gallery_top",  outdir, max_images=args.max_gallery)
    save_gallery(side_gallery, "gallery_side", outdir, max_images=args.max_gallery)

    print("\n✓ Done.")
    print(f"  particle_map.png")
    print(f"  gallery_top/          ({len(top_gallery)} files)")
    print(f"  gallery_top_overview.png")
    print(f"  gallery_side/         ({len(side_gallery)} files)")
    print(f"  gallery_side_overview.png")


if __name__ == "__main__":
    main()
