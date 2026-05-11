#!/usr/bin/env python3
"""
CellDynamicST: Cell Type Distribution Heatmap
==============================================
Refactored from: cell_type_distribution_heatmap_enhanced_v6.py

Creates a publication-quality multi-dimensional heatmap showing cell type
distribution across brain regions with:
  - Main heatmap: circle size = cell count, color = % change (treatment/control)
  - Circle border = significance (chi-squared test)
  - Outer ring = variance (CV across samples)
  - Inner glyph = dominant receptor family
  - Side annotations: stacked area charts showing composition changes
  - Mesostructure color bar, Gini coefficient, Shannon diversity

Usage from R:
  cdst_plot_distribution_heatmap(seurat_obj, genotype = "AA",
                                  control_suffix = "SAL", treatment_suffix = "MOR")

Usage from command line:
  python celltype_distribution_heatmap.py --input metadata.csv --output ./output/
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Circle, Rectangle
import warnings
import argparse
import os
import sys
from dataclasses import dataclass, field

warnings.filterwarnings("ignore")


def _scalar(val):
    """Extract scalar from potential Series (duplicate index)."""
    if hasattr(val, 'iloc'):
        return str(val.iloc[0])
    return str(val)


# =============================================================================
# CONFIGURABLE COLUMN MAPPING
# =============================================================================
@dataclass
class ColMap:
    """Column name mapping for generic use. Users override these to match
    their metadata column names."""
    cell_type: str = "cell_type_annotation"
    region: str = "brain_region"
    meso_region: str = "meso_structure"
    cell_group: str = "cell_type"
    nt_type: str = "nt_class"
    experiment_group: str = "experiment_group"
    ap_location: str = "AP_location"
    receptor_criteria: str = "receptor_criteria"

    # Legacy Allen Brain Atlas column names (from original paper)
    @classmethod
    def legacy(cls):
        """Column mapping matching the original Oprm1A118G paper data."""
        return cls(
            cell_type="cell_type_annotation",
            region="level_7_parent_roi_acronym",
            meso_region="level_3_parent_roi_acronym",
            cell_group="high_level_cell_type",
            nt_type="neurotransmitter_type",
            experiment_group="experiment_group",
            ap_location="AP_location",
            receptor_criteria="receptor_criteria",
        )


DEFAULT_COLMAP = ColMap()


# =============================================================================
# CONFIGURABLE COLOR PALETTES
# =============================================================================
DEFAULT_NT_COLORS = {
    "Glut": "#E78AC3", "GABA": "#8DA0CB", "Dopa": "#FC8D62",
    "Sero": "#66C2A5", "Hist": "#FFD92F", "Ach": "#A6D854",
    "Glutamatergic": "#E78AC3", "GABAergic": "#8DA0CB",
    "Dopaminergic": "#FC8D62", "Serotonergic": "#66C2A5",
    "Cholinergic": "#A6D854",
    "Other": "#B3B3B3",
}

DEFAULT_CELL_GROUP_COLORS = {
    "Glut": "#E78AC3", "GABA": "#8DA0CB", "Astrocyte": "#B3DE69",
    "Endothelial": "#BC80BD", "Ependymal": "#FFFFB3", "Microglia": "#BEBADA",
    "Oligodendrocyte": "#FDB462", "Neuron": "#E78AC3", "OPC": "#D9D9D9",
    "Pericyte": "#CCEBC5", "Other": "#D9D9D9",
}

DEFAULT_MESO_COLORS = {
    "Cortex": "#66C2A5", "Striatum": "#FC8D62", "Pallidum": "#8DA0CB",
    "Thalamus": "#E78AC3", "Hypothalamus": "#A6D854", "Midbrain": "#FFD92F",
    "Pons": "#E5C494", "Medulla": "#B3B3B3", "Cerebellum": "#FDBF6F",
    "Hippocampus": "#FB9A99", "Olfactory": "#CAB2D6",
    "Isocortex": "#66C2A5", "OLF": "#CAB2D6", "HPF": "#FB9A99",
    "CTXsp": "#8DD3C7", "STR": "#FC8D62", "PAL": "#8DA0CB",
    "TH": "#E78AC3", "HY": "#A6D854", "MB": "#FFD92F",
    "P": "#E5C494", "MY": "#B3B3B3", "CB": "#FDBF6F",
    "root": "#DDDDDD",
}

RECEPTOR_GLYPHS = {1: "\u25CF", 2: "\u25C6", 3: "\u25B2", 4: "\u2726"}
SIG_COLORS = {0: "#CCCCCC", 1: "#888888", 2: "#333333", 3: "#000000"}

HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "cdst_heatmap",
    ["#3B4CC0", "#6B8DF0", "#AAC7FD", "#F7F7F7", "#F7A889", "#E26952", "#B40426"],
)


# =============================================================================
# DATA PROCESSING FUNCTIONS
# =============================================================================
def filter_metadata(meta, cm=None, min_cells_per_type=20, min_cells_per_region=20):
    """Filter metadata to retain well-represented cell types and regions."""
    if cm is None:
        cm = DEFAULT_COLMAP
    ct_counts = meta[cm.cell_type].value_counts()
    valid_cts = ct_counts[ct_counts >= min_cells_per_type].index
    meta = meta[meta[cm.cell_type].isin(valid_cts)].copy()

    roi_counts = meta[cm.region].value_counts()
    valid_rois = roi_counts[roi_counts >= min_cells_per_region].index
    meta = meta[meta[cm.region].isin(valid_rois)].copy()

    return meta


def calculate_counts(meta, experiment_group, cm=None):
    """Calculate cell counts per cell_type x region for a given group."""
    if cm is None:
        cm = DEFAULT_COLMAP
    subset = meta[meta[cm.experiment_group] == experiment_group]
    counts = (
        subset.groupby([cm.cell_type, cm.region])
        .size()
        .reset_index(name="count")
    )
    return counts


def calculate_percentage_change(treatment_counts, control_counts, cm=None):
    """Calculate percentage change: (treatment - control) / control * 100."""
    if cm is None:
        cm = DEFAULT_COLMAP
    merged = pd.merge(
        treatment_counts,
        control_counts,
        on=[cm.cell_type, cm.region],
        suffixes=("_treat", "_ctrl"),
        how="outer",
    ).fillna(0)
    merged["pct_change"] = np.where(
        merged["count_ctrl"] > 0,
        (merged["count_treat"] - merged["count_ctrl"]) / merged["count_ctrl"] * 100,
        np.where(merged["count_treat"] > 0, 100, 0),
    )
    return merged


def get_ordering(meta, cell_types, regions, cm=None):
    """Order rows/columns by mesostructure and AP location."""
    if cm is None:
        cm = DEFAULT_COLMAP
    clean = meta[
        (meta[cm.cell_type].isin(cell_types))
        & (meta[cm.region].isin(regions))
    ].copy()

    # Column ordering: by mesostructure, then AP
    if cm.ap_location in clean.columns:
        col_means = clean.groupby(cm.region)[cm.ap_location].mean()
        if cm.meso_region in clean.columns:
            col_data = (
                clean[[cm.region, cm.meso_region]]
                .drop_duplicates()
            )
            col_data["ap_mean"] = col_data[cm.region].map(col_means)
            meso_means = clean.groupby(cm.meso_region)[cm.ap_location].mean()
            col_data["meso_mean"] = col_data[cm.meso_region].map(meso_means)
            col_data = col_data.sort_values(["meso_mean", "ap_mean"])
        else:
            col_data = pd.DataFrame({cm.region: list(regions)})
            col_data["ap_mean"] = col_data[cm.region].map(col_means)
            col_data = col_data.sort_values("ap_mean")
        col_order = col_data[cm.region].tolist()
    else:
        col_order = sorted(regions)

    # Row ordering: by cell group, then AP
    if cm.cell_group in clean.columns and cm.ap_location in clean.columns:
        row_means = clean.groupby(cm.cell_type)[cm.ap_location].mean()
        group_means = clean.groupby(cm.cell_group)[cm.ap_location].mean()
        row_data = (
            clean[[cm.cell_type, cm.cell_group]]
            .drop_duplicates()
        )
        row_data["ap_mean"] = row_data[cm.cell_type].map(row_means)
        row_data["group_mean"] = row_data[cm.cell_group].map(group_means)
        row_data = row_data.sort_values(["group_mean", "ap_mean"])
        row_order = row_data[cm.cell_type].tolist()
    else:
        row_order = sorted(cell_types)

    return row_order, col_order


def prepare_matrices(pct_change, ctrl_counts, row_order, col_order, cm=None):
    """Build percentage-change and size matrices aligned to row/col order."""
    if cm is None:
        cm = DEFAULT_COLMAP
    n_rows, n_cols = len(row_order), len(col_order)
    pct_matrix = np.zeros((n_rows, n_cols))
    size_matrix = np.zeros((n_rows, n_cols))

    pct_pivot = pct_change.pivot_table(
        index=cm.cell_type,
        columns=cm.region,
        values="pct_change",
        fill_value=0,
    )
    size_pivot = ctrl_counts.pivot_table(
        index=cm.cell_type,
        columns=cm.region,
        values="count",
        fill_value=0,
    )

    for i, ct in enumerate(row_order):
        for j, roi in enumerate(col_order):
            if ct in pct_pivot.index and roi in pct_pivot.columns:
                pct_matrix[i, j] = pct_pivot.loc[ct, roi]
            if ct in size_pivot.index and roi in size_pivot.columns:
                size_matrix[i, j] = size_pivot.loc[ct, roi]

    return pct_matrix, size_matrix


def calculate_significance(meta, row_labels, col_labels, ctrl_group, treat_group, cm=None):
    """Chi-squared test for each cell type x region combination."""
    if cm is None:
        cm = DEFAULT_COLMAP
    from scipy import stats as sp_stats

    n_rows, n_cols = len(row_labels), len(col_labels)
    sig_matrix = np.zeros((n_rows, n_cols))

    for i, ct in enumerate(row_labels):
        for j, roi in enumerate(col_labels):
            ctrl_n = len(
                meta[
                    (meta[cm.experiment_group] == ctrl_group)
                    & (meta[cm.cell_type] == ct)
                    & (meta[cm.region] == roi)
                ]
            )
            treat_n = len(
                meta[
                    (meta[cm.experiment_group] == treat_group)
                    & (meta[cm.cell_type] == ct)
                    & (meta[cm.region] == roi)
                ]
            )
            total = ctrl_n + treat_n
            if total >= 10:
                expected = total / 2
                chi2 = ((ctrl_n - expected) ** 2 + (treat_n - expected) ** 2) / expected
                p_val = 1 - sp_stats.chi2.cdf(chi2, df=1)
                if p_val < 0.001:
                    sig_matrix[i, j] = 3
                elif p_val < 0.01:
                    sig_matrix[i, j] = 2
                elif p_val < 0.05:
                    sig_matrix[i, j] = 1

    return sig_matrix


def calculate_proportions(meta, group_col, value_col, order, sep="-"):
    """Calculate proportions of value_col categories within each group."""
    if value_col not in meta.columns:
        # Return empty proportions if column missing
        return pd.DataFrame(0, index=order, columns=["Other"])

    if sep:
        expanded = meta.assign(
            **{value_col: meta[value_col].astype(str).str.split(sep)}
        ).explode(value_col)
    else:
        expanded = meta.copy()

    counts = expanded.groupby([group_col, value_col]).size().reset_index(name="n")
    totals = counts.groupby(group_col)["n"].transform("sum")
    counts["prop"] = counts["n"] / totals

    pivot = counts.pivot_table(
        index=group_col, columns=value_col, values="prop", fill_value=0
    )

    valid_order = [o for o in order if o in pivot.index]
    missing = [o for o in order if o not in pivot.index]
    if missing:
        missing_df = pd.DataFrame(0, index=missing, columns=pivot.columns)
        pivot = pd.concat([pivot, missing_df])
    pivot = pivot.loc[valid_order] if valid_order else pivot

    return pivot


def normalize_sizes(size_matrix, min_r=0.15, max_r=0.85):
    """Normalize size matrix to [min_r, max_r] with log transform."""
    if size_matrix.max() == 0:
        return np.full_like(size_matrix, min_r)

    lower = np.percentile(size_matrix[size_matrix > 0], 5) if np.any(size_matrix > 0) else 0
    upper = np.percentile(size_matrix, 95)
    clipped = np.clip(size_matrix, lower, upper)
    normed = (clipped - lower) / (upper - lower + 1e-10)
    log_t = np.log1p(normed)
    scaled = log_t / (np.max(log_t) + 1e-10) * (max_r - min_r) + min_r
    return scaled


# =============================================================================
# STACKED AREA CHART HELPERS
# =============================================================================
def _get_set2_colors(n):
    """Generate Set2-like colors."""
    base = [
        "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
        "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
    ]
    if n <= len(base):
        return base[:n]
    import colorsys
    colors = base.copy()
    while len(colors) < n:
        for bc in base:
            if len(colors) >= n:
                break
            r, g, b = int(bc[1:3], 16) / 255, int(bc[3:5], 16) / 255, int(bc[5:7], 16) / 255
            h, s, v = colorsys.rgb_to_hsv(r, g, b)
            h = (h + 0.1) % 1.0
            r2, g2, b2 = colorsys.hsv_to_rgb(h, s, v)
            colors.append(f"#{int(r2*255):02x}{int(g2*255):02x}{int(b2*255):02x}")
    return colors[:n]


def _get_colors(categories, color_dict):
    """Map categories to colors, using Set2 for unknowns."""
    colors = []
    unknown = []
    for cat in categories:
        if cat in color_dict:
            colors.append(color_dict[cat])
        else:
            unknown.append(cat)
            colors.append(None)
    if unknown:
        unk_colors = _get_set2_colors(len(unknown))
        unk_map = dict(zip(unknown, unk_colors))
        colors = [unk_map.get(cat, c) if c is None else c for cat, c in zip(categories, colors)]
    return colors


def _draw_stacked_area_h(ax, sal_props, mor_props, colors, n_items, cell_w=1.0):
    """Draw horizontal stacked area charts (SAL left, MOR right) as trapezoid polygons."""
    from matplotlib.patches import Polygon
    for i in range(n_items):
        x_left = i * cell_w + 0.05 * cell_w
        x_right = (i + 1) * cell_w - 0.05 * cell_w

        sal_p = sal_props[i] if i < sal_props.shape[0] else np.zeros(sal_props.shape[1])
        mor_p = mor_props[i] if i < mor_props.shape[0] else np.zeros(mor_props.shape[1])

        sal_sum = sal_p.sum()
        mor_sum = mor_p.sum()
        if sal_sum > 0:
            sal_p = sal_p / sal_sum
        if mor_sum > 0:
            mor_p = mor_p / mor_sum

        sal_cumsum = np.cumsum(sal_p)
        mor_cumsum = np.cumsum(mor_p)

        prev_sal = 0
        prev_mor = 0
        for k in range(len(colors)):
            curr_sal = sal_cumsum[k] if k < len(sal_cumsum) else 1
            curr_mor = mor_cumsum[k] if k < len(mor_cumsum) else 1
            color = colors[k] if k < len(colors) else "#AAAAAA"

            polygon = Polygon([
                (x_left, prev_sal),
                (x_left, curr_sal),
                (x_right, curr_mor),
                (x_right, prev_mor),
            ], facecolor=color, edgecolor='white', linewidth=0.3, alpha=0.9)
            ax.add_patch(polygon)
            prev_sal = curr_sal
            prev_mor = curr_mor


def _draw_stacked_area_v(ax, sal_props, mor_props, colors, n_items, cell_h=1.0):
    """Draw vertical stacked area charts (SAL top, MOR bottom) as trapezoid polygons."""
    from matplotlib.patches import Polygon
    for i in range(n_items):
        y_bottom = (n_items - i - 1) * cell_h + 0.05 * cell_h
        y_top = (n_items - i) * cell_h - 0.05 * cell_h

        sal_p = sal_props[i] if i < sal_props.shape[0] else np.zeros(sal_props.shape[1])
        mor_p = mor_props[i] if i < mor_props.shape[0] else np.zeros(mor_props.shape[1])

        sal_sum = sal_p.sum()
        mor_sum = mor_p.sum()
        if sal_sum > 0:
            sal_p = sal_p / sal_sum
        if mor_sum > 0:
            mor_p = mor_p / mor_sum

        sal_cumsum = np.cumsum(sal_p)
        mor_cumsum = np.cumsum(mor_p)

        prev_sal = 0
        prev_mor = 0
        for k in range(len(colors)):
            curr_sal = sal_cumsum[k] if k < len(sal_cumsum) else 1
            curr_mor = mor_cumsum[k] if k < len(mor_cumsum) else 1
            color = colors[k] if k < len(colors) else "#AAAAAA"

            polygon = Polygon([
                (prev_sal, y_bottom),
                (curr_sal, y_bottom),
                (curr_mor, y_top),
                (prev_mor, y_top),
            ], facecolor=color, edgecolor='white', linewidth=0.3, alpha=0.9)
            ax.add_patch(polygon)
            prev_sal = curr_sal
            prev_mor = curr_mor


# =============================================================================
# MAIN HEATMAP PLOTTING
# =============================================================================
def create_heatmap(
    pct_matrix, size_matrix, row_labels, col_labels,
    row_groups, col_mesostructures,
    nt_by_region_sal, nt_by_region_mor,
    rc_by_cell_sal, rc_by_cell_mor,
    nt_by_cell_sal, nt_by_cell_mor,
    ct_by_region_sal, ct_by_region_mor,
    title, output_path,
    significance_matrix=None,
    variance_matrix=None,
    receptor_matrix=None,
    heatmap_range=None,
    nt_colors=None, cell_group_colors=None, meso_colors=None,
):
    """Create the full multi-dimensional heatmap with all annotations."""
    if nt_colors is None:
        nt_colors = DEFAULT_NT_COLORS
    if cell_group_colors is None:
        cell_group_colors = DEFAULT_CELL_GROUP_COLORS
    if meso_colors is None:
        meso_colors = DEFAULT_MESO_COLORS

    n_rows, n_cols = len(row_labels), len(col_labels)
    if n_rows == 0 or n_cols == 0:
        print("Warning: Empty row or column labels. Skipping heatmap.")
        return

    # Layout proportions
    fig = plt.figure(figsize=(max(24, n_cols * 0.4), max(16, n_rows * 0.35)))

    # Axes positions [left, bottom, width, height]
    w_cg = 0.015
    w_area = 0.04
    w_main = 0.55
    h_main = 0.55
    h_area = 0.04
    h_meso = 0.015

    x_cg = 0.06
    x_rc = x_cg + w_cg + 0.005
    x_main = x_rc + w_area + 0.005
    y_main = 0.18
    y_meso = y_main + h_main + 0.005
    y_nt_top = y_meso + h_meso + 0.005
    y_ct_top = y_nt_top + h_area + 0.005

    # Heatmap range
    if heatmap_range is None:
        vmin = np.percentile(pct_matrix, 10)
        vmax = np.percentile(pct_matrix, 90)
    else:
        vmin, vmax = heatmap_range
    if vmin >= 0:
        vmin = -1  # Ensure diverging colormap works
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=max(vmax, 1))
    sizes = normalize_sizes(size_matrix)

    # ---- MAIN HEATMAP ----
    ax_main = fig.add_axes([x_main, y_main, w_main, h_main])
    for i in range(n_rows):
        for j in range(n_cols):
            x, y = j + 0.5, n_rows - i - 0.5
            radius = sizes[i, j] * 0.45
            if radius < 0.06:
                continue

            pct = pct_matrix[i, j]
            fill = HEATMAP_CMAP(norm(pct))
            sig = int(significance_matrix[i, j]) if significance_matrix is not None else 0
            sig = min(sig, 3)
            var_val = variance_matrix[i, j] if variance_matrix is not None else 0
            rec = int(receptor_matrix[i, j]) if receptor_matrix is not None else 0

            # Variance ring
            if var_val > 0.15:
                ring_w = min(var_val * 0.15, 0.08)
                ax_main.add_patch(
                    Circle((x, y), radius + ring_w, facecolor="none",
                           edgecolor="#555555", linewidth=var_val * 3, alpha=0.6)
                )

            # Main circle with significance border
            bw = 0.8 + sig * 0.6
            ax_main.add_patch(
                Circle((x, y), radius, facecolor=fill,
                       edgecolor=SIG_COLORS[sig], linewidth=bw)
            )

            # Receptor glyph
            if rec > 0 and radius > 0.12:
                glyph = RECEPTOR_GLYPHS.get(rec, "")
                if glyph:
                    ax_main.text(x, y, glyph, fontsize=max(radius * 32, 7),
                                ha="center", va="center", color="white",
                                fontweight="bold", alpha=0.95)

    ax_main.set_xlim(0, n_cols)
    ax_main.set_ylim(0, n_rows)
    ax_main.set_aspect("equal")
    ax_main.set_title(title, fontsize=15, fontweight="bold", pad=10)
    for spine in ax_main.spines.values():
        spine.set_visible(False)

    # Row labels (right)
    ax_main.set_yticks([n_rows - i - 0.5 for i in range(n_rows)])
    ax_main.set_yticklabels(row_labels, fontsize=9)
    # Column labels (bottom)
    ax_main.set_xticks([j + 0.5 for j in range(n_cols)])
    ax_main.set_xticklabels(col_labels, fontsize=9, rotation=90, ha="center")

    # ---- LEFT: Cell Group bar ----
    ax_cg = fig.add_axes([x_cg, y_main, w_cg, h_main])
    boundaries = [0]
    current = row_groups[0] if row_groups else ""
    for i, g in enumerate(row_groups):
        if g != current:
            boundaries.append(i)
            current = g
    boundaries.append(n_rows)

    for i in range(len(boundaries) - 1):
        group = row_groups[boundaries[i]]
        color = cell_group_colors.get(group, "#AAAAAA")
        ax_cg.add_patch(
            Rectangle((0, n_rows - boundaries[i + 1]), 1,
                      boundaries[i + 1] - boundaries[i],
                      facecolor=color, edgecolor="white", linewidth=0.5)
        )
    ax_cg.set_xlim(0, 1)
    ax_cg.set_ylim(0, n_rows)
    ax_cg.axis("off")

    # ---- LEFT: RC by Cell Type stacked area ----
    ax_rc = fig.add_axes([x_rc, y_main, w_area, h_main])
    rc_colors = _get_set2_colors(rc_by_cell_sal.shape[1])
    _draw_stacked_area_v(ax_rc, rc_by_cell_sal.values, rc_by_cell_mor.values,
                         rc_colors, n_rows)
    ax_rc.set_xlim(0, 1)
    ax_rc.set_ylim(0, n_rows)
    ax_rc.axis("off")

    # ---- TOP: Mesostructure bar ----
    ax_meso = fig.add_axes([x_main, y_meso, w_main, h_meso])
    meso_boundaries = [0]
    current_meso = col_mesostructures[0] if col_mesostructures else ""
    for j, m in enumerate(col_mesostructures):
        if m != current_meso:
            meso_boundaries.append(j)
            current_meso = m
    meso_boundaries.append(n_cols)

    for j in range(len(meso_boundaries) - 1):
        meso = col_mesostructures[meso_boundaries[j]]
        color = meso_colors.get(meso, "grey")
        ax_meso.add_patch(
            Rectangle((meso_boundaries[j], 0),
                      meso_boundaries[j + 1] - meso_boundaries[j], 1,
                      facecolor=color, edgecolor="white", linewidth=0.5)
        )
    ax_meso.set_xlim(0, n_cols)
    ax_meso.set_ylim(0, 1)
    ax_meso.axis("off")

    # ---- TOP: NT by Region stacked area ----
    ax_nt_top = fig.add_axes([x_main, y_nt_top, w_main, h_area])
    nt_cols_list = list(nt_by_region_sal.columns)
    nt_col_colors = _get_colors(nt_cols_list, nt_colors)
    _draw_stacked_area_h(ax_nt_top, nt_by_region_sal.values,
                         nt_by_region_mor.values, nt_col_colors, n_cols)
    ax_nt_top.set_xlim(0, n_cols)
    ax_nt_top.set_ylim(0, 1)
    ax_nt_top.axis("off")

    # ---- TOP: CT by Region stacked area ----
    ax_ct_top = fig.add_axes([x_main, y_ct_top, w_main, h_area])
    ct_cols_list = list(ct_by_region_sal.columns)
    ct_col_colors = _get_colors(ct_cols_list, cell_group_colors)
    _draw_stacked_area_h(ax_ct_top, ct_by_region_sal.values,
                         ct_by_region_mor.values, ct_col_colors, n_cols)
    ax_ct_top.set_xlim(0, n_cols)
    ax_ct_top.set_ylim(0, 1)
    ax_ct_top.axis("off")

    # ---- RIGHT: NT by Cell Type stacked area ----
    x_nt_right = x_main + w_main + 0.005
    ax_nt_right = fig.add_axes([x_nt_right, y_main, w_area, h_main])
    _draw_stacked_area_v(ax_nt_right, nt_by_cell_sal.values,
                         nt_by_cell_mor.values, nt_col_colors, n_rows)
    ax_nt_right.set_xlim(0, 1)
    ax_nt_right.set_ylim(0, n_rows)
    ax_nt_right.axis("off")

    # ---- BOTTOM: CT stacked area ----
    y_bottom = y_main - h_area - 0.01
    ax_bottom = fig.add_axes([x_main, y_bottom, w_main, h_area])
    _draw_stacked_area_h(ax_bottom, ct_by_region_sal.values,
                         ct_by_region_mor.values, ct_col_colors, n_cols)
    ax_bottom.set_xlim(0, n_cols)
    ax_bottom.set_ylim(0, 1)
    ax_bottom.axis("off")

    # ---- ANNOTATION LABELS ----
    label_x = x_cg - 0.04
    for label_text, label_y in [
        ("Cell\nGroup", y_main + h_main / 2),
        ("RC", y_main + h_main / 2),
        ("NT by Region\n(SAL\u2192MOR)", y_nt_top + h_area / 2),
        ("CT by Region\n(SAL\u2192MOR)", y_ct_top + h_area / 2),
        ("Meso", y_meso + h_meso / 2),
    ]:
        fig.text(label_x, label_y, label_text, fontsize=10, fontweight="bold",
                 ha="right", va="center", rotation=0)

    fig.text(x_nt_right + w_area + 0.01, y_main + h_main / 2,
             "NT by Cell\n(SAL\u2192MOR)", fontsize=10, fontweight="bold",
             ha="left", va="center")

    # ---- COLORBAR ----
    cbar_ax = fig.add_axes([x_main + w_main + w_area + 0.04,
                            y_main, 0.015, h_main * 0.4])
    sm = plt.cm.ScalarMappable(cmap=HEATMAP_CMAP, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, cax=cbar_ax, label="% Change (Treatment / Control)")

    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved heatmap: {output_path}")


# =============================================================================
# MAIN PROCESSING PIPELINE
# =============================================================================
def process_and_visualize(
    meta_data, output_dir, genotype="AA",
    control_suffix="SAL", treatment_suffix="MOR",
    nt_colors=None, cell_group_colors=None, meso_colors=None,
    colmap=None,
):
    """Full pipeline: filter -> compute -> plot.

    Parameters
    ----------
    meta_data : pd.DataFrame
        Cell-level metadata.
    output_dir : str
        Output directory.
    genotype : str
        Genotype to analyze (e.g., "AA").
    control_suffix, treatment_suffix : str
        Group suffixes.
    colmap : ColMap or None
        Column name mapping. If None, auto-detected.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Auto-detect column mapping
    if colmap is None:
        if "level_7_parent_roi_acronym" in meta_data.columns:
            cm = ColMap.legacy()
        else:
            cm = DEFAULT_COLMAP
    else:
        cm = colmap

    meta_data = filter_metadata(meta_data, cm)

    # Detect group name format: "AA SAL" vs "AA_SAL"
    groups = meta_data[cm.experiment_group].unique()
    ctrl_group = f"{genotype}_{control_suffix}"
    treat_group = f"{genotype}_{treatment_suffix}"
    if ctrl_group not in groups:
        ctrl_group = f"{genotype} {control_suffix}"
        treat_group = f"{genotype} {treatment_suffix}"
    if ctrl_group not in groups:
        print(f"Warning: Control group '{ctrl_group}' not found. Available: {list(groups)}")
        return None

    ctrl_counts = calculate_counts(meta_data, ctrl_group, cm)
    treat_counts = calculate_counts(meta_data, treat_group, cm)

    common_cts = set(ctrl_counts[cm.cell_type]) & set(treat_counts[cm.cell_type])
    common_rois = set(ctrl_counts[cm.region]) & set(treat_counts[cm.region])

    if not common_cts or not common_rois:
        print("Warning: No common cell types or regions between groups.")
        return None

    ctrl_counts = ctrl_counts[
        ctrl_counts[cm.cell_type].isin(common_cts)
        & ctrl_counts[cm.region].isin(common_rois)
    ]
    treat_counts = treat_counts[
        treat_counts[cm.cell_type].isin(common_cts)
        & treat_counts[cm.region].isin(common_rois)
    ]

    pct_change = calculate_percentage_change(treat_counts, ctrl_counts, cm)
    row_order, col_order = get_ordering(meta_data, common_cts, common_rois, cm)

    pct_matrix, size_matrix = prepare_matrices(pct_change, ctrl_counts, row_order, col_order, cm)

    print("  Computing significance...")
    sig_matrix = calculate_significance(meta_data, row_order, col_order, ctrl_group, treat_group, cm)

    # Proportion annotations
    ctrl_meta = meta_data[meta_data[cm.experiment_group] == ctrl_group]
    treat_meta = meta_data[meta_data[cm.experiment_group] == treat_group]

    def _aligned_props(sal_m, mor_m, group_col, val_col, order, sep="-"):
        sal_p = calculate_proportions(sal_m, group_col, val_col, order, sep)
        mor_p = calculate_proportions(mor_m, group_col, val_col, order, sep)
        all_cols = sorted(set(sal_p.columns) | set(mor_p.columns))
        for c in all_cols:
            if c not in sal_p.columns:
                sal_p[c] = 0
            if c not in mor_p.columns:
                mor_p[c] = 0
        return sal_p[all_cols], mor_p[all_cols]

    nt_reg_s, nt_reg_m = _aligned_props(ctrl_meta, treat_meta,
                                          cm.region,
                                          cm.nt_type, col_order)
    rc_cell_s, rc_cell_m = _aligned_props(ctrl_meta, treat_meta,
                                            cm.cell_type,
                                            cm.receptor_criteria, row_order)
    nt_cell_s, nt_cell_m = _aligned_props(ctrl_meta, treat_meta,
                                            cm.cell_type,
                                            cm.nt_type, row_order)

    # CT by region (use primary cell type)
    def _primary(x):
        s = str(x).split("-")[0].strip()
        return s if s else "Other"

    ctrl_ct = ctrl_meta.copy()
    treat_ct = treat_meta.copy()
    if cm.cell_group in ctrl_ct.columns:
        ctrl_ct["primary_ct"] = ctrl_ct[cm.cell_group].apply(_primary)
        treat_ct["primary_ct"] = treat_ct[cm.cell_group].apply(_primary)
    else:
        ctrl_ct["primary_ct"] = "Other"
        treat_ct["primary_ct"] = "Other"
    ct_reg_s = calculate_proportions(ctrl_ct, cm.region, "primary_ct", col_order, sep=None)
    ct_reg_m = calculate_proportions(treat_ct, cm.region, "primary_ct", col_order, sep=None)
    all_ct_cols = sorted(set(ct_reg_s.columns) | set(ct_reg_m.columns))
    for c in all_ct_cols:
        if c not in ct_reg_s.columns:
            ct_reg_s[c] = 0
        if c not in ct_reg_m.columns:
            ct_reg_m[c] = 0
    ct_reg_s, ct_reg_m = ct_reg_s[all_ct_cols], ct_reg_m[all_ct_cols]

    # Row groups and mesostructures
    if cm.cell_group in meta_data.columns:
        row_data = (
            meta_data[meta_data[cm.cell_type].isin(row_order)]
            [[cm.cell_type, cm.cell_group]]
            .drop_duplicates()
            .set_index(cm.cell_type)
        )
        row_groups = [
            _scalar(row_data.loc[r, cm.cell_group]) if r in row_data.index else "Other"
            for r in row_order
        ]
    else:
        row_groups = ["Other"] * len(row_order)

    if cm.meso_region in meta_data.columns:
        col_data = (
            meta_data[meta_data[cm.region].isin(col_order)]
            [[cm.region, cm.meso_region]]
            .drop_duplicates()
            .set_index(cm.region)
        )
        col_meso = [
            _scalar(col_data.loc[c, cm.meso_region]) if c in col_data.index else "root"
            for c in col_order
        ]
    else:
        col_meso = ["root"] * len(col_order)

    heatmap_range = (np.percentile(pct_matrix, 10), np.percentile(pct_matrix, 90))

    out_path = os.path.join(output_dir, f"{genotype}_{treatment_suffix}_vs_{control_suffix}_heatmap.png")
    create_heatmap(
        pct_matrix=pct_matrix, size_matrix=size_matrix,
        row_labels=row_order, col_labels=col_order,
        row_groups=row_groups, col_mesostructures=col_meso,
        nt_by_region_sal=nt_reg_s, nt_by_region_mor=nt_reg_m,
        rc_by_cell_sal=rc_cell_s, rc_by_cell_mor=rc_cell_m,
        nt_by_cell_sal=nt_cell_s, nt_by_cell_mor=nt_cell_m,
        ct_by_region_sal=ct_reg_s, ct_by_region_mor=ct_reg_m,
        title=f"Cell Type Distribution: {genotype} {treatment_suffix} vs {control_suffix}",
        output_path=out_path,
        significance_matrix=sig_matrix,
        heatmap_range=heatmap_range,
        nt_colors=nt_colors, cell_group_colors=cell_group_colors,
        meso_colors=meso_colors,
    )

    return {
        "pct_matrix": pct_matrix, "size_matrix": size_matrix,
        "significance_matrix": sig_matrix,
        "row_labels": row_order, "col_labels": col_order,
    }


# =============================================================================
# CLI ENTRY POINT
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="CellDynamicST: Cell Type Distribution Heatmap"
    )
    parser.add_argument("--input", "-i", required=True, help="Metadata CSV")
    parser.add_argument("--output", "-o", default="./output/", help="Output dir")
    parser.add_argument("--genotypes", "-g", nargs="+", default=["AA", "GG"])
    parser.add_argument("--control", default="SAL", help="Control suffix")
    parser.add_argument("--treatment", default="MOR", help="Treatment suffix")
    # Column mapping overrides
    parser.add_argument("--col-celltype", default=None, help="Cell type column name")
    parser.add_argument("--col-region", default=None, help="Region column name")
    parser.add_argument("--col-meso", default=None, help="Mesostructure column name")
    parser.add_argument("--col-cellgroup", default=None, help="Cell group column name")
    parser.add_argument("--col-nt", default=None, help="NT type column name")
    parser.add_argument("--col-group", default=None, help="Experiment group column name")
    parser.add_argument("--legacy-columns", action="store_true",
                        help="Use legacy Allen Brain Atlas column names")
    args = parser.parse_args()

    # Build column mapping
    if args.legacy_columns:
        cm = ColMap.legacy()
    else:
        cm = ColMap()
        if args.col_celltype:
            cm.cell_type = args.col_celltype
        if args.col_region:
            cm.region = args.col_region
        if args.col_meso:
            cm.meso_region = args.col_meso
        if args.col_cellgroup:
            cm.cell_group = args.col_cellgroup
        if args.col_nt:
            cm.nt_type = args.col_nt
        if args.col_group:
            cm.experiment_group = args.col_group

    meta = pd.read_csv(args.input, index_col=0)

    # Auto-detect if legacy columns present
    if "level_7_parent_roi_acronym" in meta.columns and not args.legacy_columns:
        print("Detected legacy column names. Using legacy column mapping.")
        cm = ColMap.legacy()

    for genotype in args.genotypes:
        print(f"\nProcessing {genotype}...")
        try:
            process_and_visualize(meta, args.output, genotype,
                                  args.control, args.treatment, colmap=cm)
        except Exception as e:
            print(f"Error processing {genotype}: {e}")
            import traceback
            traceback.print_exc()
    print("\nDone!")


if __name__ == "__main__":
    main()
