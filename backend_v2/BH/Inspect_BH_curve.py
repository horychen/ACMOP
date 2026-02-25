"""
Inspect all BH curve data in the directory of this script.
- Files not ending with .BH are converted to .BH format (content normalized, saved as *.BH).
- Plots all BH curves for comparison.
"""

from pathlib import Path
import numpy as np

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# Extensions treated as BH curve data (besides .BH)
DATA_EXTENSIONS = {".tab", ".csv", ".txt", ".dat"}


def script_dir():
    return Path(__file__).resolve().parent


def is_bh_curve_file(path: Path) -> bool:
    """Return True if path looks like a BH curve data file (by extension)."""
    ext = path.suffix.lower()
    return ext == ".bh" or ext in DATA_EXTENSIONS


def read_bh_data(path: Path) -> np.ndarray:
    """
    Read H (A/m) and B (T) from a file. Returns (N, 2) array [H, B].
    Handles .BH (no header, space/tab separated) and .tab (optional header).
    """
    raw = path.read_text(encoding="utf-8", errors="replace").strip()
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]

    if not lines:
        return np.array([]).reshape(0, 2)

    # Detect header: first line might be column names (e.g. "H (A_per_meter)" \t "B (tesla)")
    first = lines[0]
    if first.lstrip().startswith('"') or first.lower().startswith("h ") or "a_per_meter" in first.lower() or "tesla" in first.lower():
        lines = lines[1:]

    rows = []
    for ln in lines:
        parts = ln.replace("\t", " ").split()
        if len(parts) >= 2:
            try:
                h = float(parts[0])
                b = float(parts[1])
                rows.append([h, b])
            except ValueError:
                continue
    return np.array(rows) if rows else np.array([]).reshape(0, 2)


def write_bh_file(path: Path, data: np.ndarray) -> None:
    """Write H, B data to a .BH file (no header, space-separated)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for row in data:
            f.write(f"{row[0]}\t{row[1]}\n")


def ensure_bh_file(source: Path) -> Path:
    """
    If source does not end with .BH, read it, normalize to H,B and save as same stem with .BH.
    Returns the path to the .BH file (same as source if already .BH, else new path).
    """
    if source.suffix.lower() == ".bh":
        return source
    data = read_bh_data(source)
    if data.size == 0:
        return source
    bh_path = source.with_suffix(".BH")
    write_bh_file(bh_path, data)
    return bh_path


def collect_and_convert_bh_curves(dir_path: Path):
    """
    Find all BH curve data files in dir_path. Convert non-.BH to .BH.
    Return list of (label, data_array) and set of .BH paths.
    """
    results = []  # (label, data)
    bh_paths = set()

    for path in sorted(dir_path.iterdir()):
        if path.is_dir() or path.name.startswith("."):
            continue
        if path.suffix.lower() == ".py":
            continue
        if not is_bh_curve_file(path):
            continue

        data = read_bh_data(path)
        if data.size == 0:
            continue

        bh_path = ensure_bh_file(path)
        bh_paths.add(bh_path)
        label = bh_path.stem
        results.append((label, data))

    return results, bh_paths


def plot_bh_curves(labels_and_data, out_path: Path | None = None):
    """Plot all BH curves on one figure for comparison."""
    if not labels_and_data:
        print("No BH curve data to plot.")
        return
    if not HAS_MATPLOTLIB:
        print("Optional: install matplotlib to show BH curve comparison plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    for label, data in labels_and_data:
        H = data[:, 0]
        B = data[:, 1]
        ax.plot(H, B, "-", label=label, linewidth=1.5)

    ax.set_xlabel("H (A/m)")
    ax.set_ylabel("B (T)")
    ax.set_title("BH curves comparison")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.tight_layout()

    if out_path:
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved figure: {out_path}")
    plt.show()


def main():
    base = script_dir()
    print(f"Inspecting BH curve data in: {base}")

    labels_and_data, bh_paths = collect_and_convert_bh_curves(base)

    if not labels_and_data:
        print("No BH curve data files found.")
        return

    print(f"Found {len(labels_and_data)} BH curve(s):")
    for label, data in labels_and_data:
        print(f"  - {label}: {len(data)} points, H in [{data[:, 0].min():.2f}, {data[:, 0].max():.2f}] A/m, B in [{data[:, 1].min():.4f}, {data[:, 1].max():.4f}] T")

    plot_bh_curves(labels_and_data, out_path=base / "BH_curves_comparison.png")


if __name__ == "__main__":
    main()
