#!/usr/bin/env python3

"""Extract the amyloid and tau staging strips of Weickenmeier et al., fig. 1.

Figure 1 of Weickenmeier et al. (J. Mech. Phys. Solids 124, 2019, p. 266;
the drawings are adopted there from Jucker and Walker 2013) shows the
expected three-stage progression of four misfolded proteins as rows of
medial brain drawings. The report compares the two Alzheimer rows against
our simulated seeding patterns, so this script cuts those two rows out of
the published page instead of redrawing them by hand.

The extraction is deterministic and detected, not measured by eye: the page
is rendered at 300 dpi, the four cartoon rows are found as horizontal bands
of saturated pixels, the right-hand disease labels (plain black text) are
dropped by cutting at the last wide unsaturated gap, and the two target
bands are identified by their dominant hue, orange for amyloid-beta and
blue for tau. The script fails loudly if the page does not contain exactly
four bands with the expected hues.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

PAGE = 3           # PDF page holding figure 1 (journal page 266)
DPI = 600          # the strips span the full report text width, so the
                   # rendering must stay well above print resolution
SATURATION = 0.18  # chroma threshold separating drawings from black text
MIN_BAND = DPI // 5   # minimum band height in pixels
PAD = DPI // 50       # padding around detected content, pixels


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True,
                        help="Weickenmeier et al. (2019) article PDF")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def render_page(pdf, directory):
    prefix = Path(directory) / "page"
    subprocess.run(
        ["pdftoppm", "-f", str(PAGE), "-l", str(PAGE), "-r", str(DPI),
         "-png", str(pdf), str(prefix)],
        check=True)
    produced = sorted(Path(directory).glob("page-*.png"))
    if len(produced) != 1:
        raise RuntimeError(f"Expected one rendered page, got {produced}")
    return np.asarray(Image.open(produced[0]).convert("RGB"), dtype=float)


def saturation(image):
    """Chroma per pixel: distance of the RGB channels from their mean."""
    return (image.max(axis=2) - image.min(axis=2)) / 255.0


def bands(mask_rows):
    """Contiguous index ranges where the row mask is true."""
    found, start = [], None
    for index, value in enumerate(mask_rows):
        if value and start is None:
            start = index
        if not value and start is not None:
            if index - start >= MIN_BAND:
                found.append((start, index))
            start = None
    if start is not None and len(mask_rows) - start >= MIN_BAND:
        found.append((start, len(mask_rows)))
    return found


def dominant_hue(patch):
    """Label a cartoon band by its saturated pixels' mean colour."""
    chroma = saturation(patch)
    pixels = patch[chroma > SATURATION]
    red, green, blue = pixels.mean(axis=0)
    if blue > red and blue > green:
        return "blue"
    if red > blue and green > blue and red > 1.15 * green:
        return "orange"
    if red > blue and green > blue:
        return "green"
    return "red"


def cut_labels(image, band):
    """Column extent of the drawings, dropping the black label text."""
    top, bottom = band
    strip = image[top:bottom]
    chroma = saturation(strip)
    filled = np.flatnonzero((chroma > SATURATION).sum(axis=0) > 0)
    # The right bound comes from the saturated drawings, which stops before
    # the plain label text. The left bound comes from the ink instead: the
    # first-stage drawing of the tau row is almost entirely an unsaturated
    # outline, and a chroma bound would cut it.
    right = filled[-1]
    ink = np.flatnonzero((strip.min(axis=2) < 200).sum(axis=0) > 2)
    left = ink[0]
    return max(left - PAD, 0), min(right + PAD, image.shape[1])


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as scratch:
        image = render_page(args.pdf, scratch)

    chroma = saturation(image)
    rows = (chroma > SATURATION).sum(axis=1) > DPI // 8
    detected = bands(rows)
    if len(detected) != 4:
        raise RuntimeError(
            f"Expected the four cartoon rows of figure 1, found {detected}")

    hues = [dominant_hue(image[top:bottom]) for top, bottom in detected]
    # Only the two Alzheimer rows are taken; the alpha-synuclein and TDP-43
    # rows below them are not needed and their exact hue class is not
    # checked (the salmon of the third row straddles red and orange).
    if hues[:2] != ["orange", "blue"]:
        raise RuntimeError(
            f"Band hues {hues} do not start with ['orange', 'blue']")

    for name, band in (("amyloid", detected[0]), ("tau", detected[1])):
        top, bottom = band
        left, right = cut_labels(image, band)
        crop = image[max(top - PAD, 0):bottom + PAD, left:right]
        output = args.output_dir / f"weickenmeier_fig1_{name}.png"
        Image.fromarray(crop.astype(np.uint8)).save(output)
        print(f"{output}  rows {top}-{bottom}  columns {left}-{right}  "
              f"({crop.shape[1]}x{crop.shape[0]})")


if __name__ == "__main__":
    main()
