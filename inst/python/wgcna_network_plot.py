#!/usr/bin/env python3
"""
CellDynamicST: WGCNA Differential Gene Coexpression Network Visualization
==========================================================================
Refactored from: WGCNA_module_differential_gene_coexpression_network.py
                 WGCNA_circos.py
                 Interregional_WCGNA_Module_hub_gene.R
                 gene_coexpression_network_atlas_plot.py

Creates publication-quality network visualizations:
  1. Brain atlas overlay: nodes = brain regions, edges = module correlations,
     overlaid on Allen Brain Atlas annotation volume slice (sagittal/coronal/horizontal)
  2. Circos plot: 5 enrichment rings (OUD, 4 DEG effects) + module color ring
  3. Hub gene subnetwork: gold hub nodes + sky blue neighbor nodes

Usage from command line:
  python wgcna_network_plot.py atlas  --nodes n.csv --edges e.csv --output out.png --plane sagittal --slice 510
  python wgcna_network_plot.py circos --modules m.csv --output out.png
  python wgcna_network_plot.py hub    --hubs h.csv --edges e.csv --output out.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, Normalize
from matplotlib.patches import FancyArrowPatch, Circle, Wedge, Ellipse, Polygon, Patch
import matplotlib.patheffects as pe
import warnings
import argparse
import os

warnings.filterwarnings("ignore")


MODULE_COLORS = {
    "blue": "#1F77B4", "brown": "#8C564B", "green": "#2CA02C",
    "red": "#D62728", "turquoise": "#17BECF", "yellow": "#BCBD22",
    "black": "#333333", "pink": "#E377C2", "magenta": "#9467BD",
    "purple": "#7B4173", "greenyellow": "#B5CF6B", "tan": "#C49C94",
    "salmon": "#FF9896", "cyan": "#9EDAE5", "midnightblue": "#003366",
    "lightcyan": "#E0FFFF", "lightgreen": "#90EE90",
    "lightyellow": "#FFFFE0", "royalblue": "#4169E1", "darkred": "#8B0000",
    "darkgreen": "#006400", "darkturquoise": "#00CED1", "darkgrey": "#A9A9A9",
    "orange": "#FF7F0E", "darkorange": "#FF8C00", "white": "#F0F0F0",
    "skyblue": "#87CEEB", "saddlebrown": "#8B4513", "steelblue": "#4682B4",
    "paleturquoise": "#AFEEEE", "violet": "#EE82EE", "darkolivegreen": "#556B2F",
    "grey": "#CCCCCC", "grey60": "#999999",
}


def _norm(series):
    """Min-max normalize a series to [0, 1]."""
    s = series.fillna(series.min())
    rng = s.max() - s.min()
    if rng < 1e-9:
        return pd.Series(0.5, index=s.index)
    return (s - s.min()) / rng


# =============================================================================
# ALLEN BRAIN ATLAS HELPERS
# =============================================================================
def _load_allen_annotation(cache_dir=None, resolution=10):
    """
    Download (or load from cache) the Allen Mouse Brain annotation volume
    and structure tree via AllenSDK.

    Returns
    -------
    annotation : np.ndarray  shape (AP, DV, ML) at given resolution
    structure_tree : StructureTree
    """
    try:
        from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
    except ImportError:
        raise ImportError(
            "AllenSDK is required for atlas overlay. "
            "Install with: pip install allensdk"
        )
    if cache_dir is None:
        cache_dir = os.path.join(os.path.expanduser("~"), ".allen_cache")
    os.makedirs(cache_dir, exist_ok=True)
    mcc = MouseConnectivityCache(
        resolution=resolution,
        manifest_file=os.path.join(cache_dir, "manifest.json")
    )
    annotation, _ = mcc.get_annotation_volume()
    structure_tree = mcc.get_structure_tree()
    return annotation, structure_tree


def _get_region_centers_3d(annotation, structure_tree):
    """
    Compute 3D centroid of each brain structure from the FULL annotation volume.
    This matches the original script's approach: centers are computed globally,
    then projected onto the chosen slice plane.

    Returns
    -------
    region_centers : dict  {structure_id: np.array([AP, DV, ML])}
    """
    from tqdm import tqdm

    region_centers = {}
    unique_structures = np.unique(annotation)

    for struct_id in tqdm(unique_structures, desc="Computing 3D centers"):
        if struct_id != 0:
            coords = np.argwhere(annotation == struct_id)
            if coords.size > 0:
                center = coords.mean(axis=0)  # [AP, DV, ML]
                # Store as [DV, AP] (flipped first two) — matches original np.flip(center[:2])
                region_centers[int(struct_id)] = np.flip(center[:2])
            else:
                print(f"Warning: No coordinates found for structure ID '{struct_id}'.")

    return region_centers


def _map_acronyms_to_ids(structure_tree, acronyms):
    """Map a list of region acronyms to Allen structure IDs."""
    acronym_to_id = {}
    for acr in acronyms:
        try:
            structs = structure_tree.get_structures_by_acronym([acr])
            if structs:
                acronym_to_id[acr] = structs[0]["id"]
        except (KeyError, IndexError):
            pass
    return acronym_to_id


def _get_parent_center(structure_tree, parent_id, region_centers):
    """
    For a parent region (e.g., Isocortex), compute the centroid by averaging
    the centroids of all descendant leaf structures that have centers.
    """
    try:
        desc_ids = set(structure_tree.descendant_ids([parent_id])[0])
    except Exception:
        desc_ids = {parent_id}

    coords = []
    for did in desc_ids:
        if did in region_centers:
            coords.append(region_centers[did])

    if coords:
        return np.mean(coords, axis=0)
    return None


# =============================================================================
# 1. BRAIN ATLAS OVERLAY NETWORK (AllenSDK-based)
# =============================================================================
def plot_atlas_network(
    nodes_df, edges_df, output_path,
    plane="sagittal",
    slice_idx=None,
    resolution=10,
    cache_dir=None,
    title=None,
    node_size_col="n_genes",
    node_color_col="module_color",
    edge_weight_col="correlation",
    edge_color_col=None,
    min_edge_weight=0.0,
    figsize=(20, 10),
):
    """
    Plot WGCNA network overlaid on Allen Brain Atlas annotation volume.
    Coordinate system matches the original gene_coexpression_network_atlas_plot.py:
    - Region centers computed from FULL 3D volume (not just the slice)
    - Slice drawn as np.flipud(np.transpose(annotation[:, :, z]))
    - Node coordinates: flipud(transpose(center)), then y = shape[1] - y

    Parameters
    ----------
    nodes_df : pd.DataFrame
        Must have 'region' column (Allen acronyms), plus size/color columns.
    edges_df : pd.DataFrame
        Must have 'source', 'target' (region acronyms), and weight column.
    plane : str
        'sagittal', 'coronal', or 'horizontal'.
    slice_idx : int or None
        Slice coordinate along the chosen axis. None = midpoint.
        For sagittal: ML coordinate (z). For coronal: AP coordinate (x).
        For horizontal: DV coordinate (y).
    resolution : int
        Allen annotation volume resolution in microns (10 or 25).
    """
    import networkx as nx
    from sklearn.preprocessing import MinMaxScaler
    try:
        from adjustText import adjust_text
        has_adjust = True
    except ImportError:
        has_adjust = False

    # Load Allen annotation
    print(f"Loading Allen Brain Atlas (resolution={resolution}um)...")
    annotation, structure_tree = _load_allen_annotation(cache_dir, resolution)
    print(f"  Annotation volume shape: {annotation.shape} (AP, DV, ML)")

    # Compute region centers from FULL 3D volume (matching original script)
    print("Computing region centers from full 3D volume...")
    region_centers = _get_region_centers_3d(annotation, structure_tree)
    print(f"  Computed centers for {len(region_centers)} structures")

    # Determine slice index
    shape = annotation.shape  # (AP, DV, ML)
    if plane == "sagittal":
        if slice_idx is None:
            slice_idx = shape[2] // 2
        slice_idx = min(max(0, slice_idx), shape[2] - 1)
        # Original: annotation_slice = np.flipud(np.transpose(annotation[:, :, z]))
        annotation_slice = np.flipud(np.transpose(annotation[:, :, slice_idx]))
    elif plane == "coronal":
        if slice_idx is None:
            slice_idx = shape[0] // 2
        slice_idx = min(max(0, slice_idx), shape[0] - 1)
        annotation_slice = np.flipud(np.transpose(annotation[slice_idx, :, :]))
    elif plane == "horizontal":
        if slice_idx is None:
            slice_idx = shape[1] // 2
        slice_idx = min(max(0, slice_idx), shape[1] - 1)
        annotation_slice = np.flipud(np.transpose(annotation[:, slice_idx, :]))
    else:
        raise ValueError(f"Unknown plane '{plane}'. Use sagittal/coronal/horizontal.")

    print(f"  Plane: {plane}, slice index: {slice_idx}")
    print(f"  Annotation slice shape: {annotation_slice.shape}")

    # Map acronyms to IDs
    all_acronyms = set(nodes_df["region"].unique())
    if "source" in edges_df.columns:
        all_acronyms |= set(edges_df["source"].unique())
    if "target" in edges_df.columns:
        all_acronyms |= set(edges_df["target"].unique())
    acronym_to_id = _map_acronyms_to_ids(structure_tree, all_acronyms)
    print(f"  Mapped {len(acronym_to_id)}/{len(all_acronyms)} acronyms to IDs")

    # For parent regions not directly in region_centers, aggregate descendants
    for acr, sid in acronym_to_id.items():
        if sid not in region_centers:
            parent_center = _get_parent_center(structure_tree, sid, region_centers)
            if parent_center is not None:
                region_centers[sid] = parent_center
                print(f"  Aggregated descendants for '{acr}' (id={sid})")
            else:
                print(f"  No descendant centers for '{acr}'")

    # Build figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Draw brain structure contours from annotation slice
    unique_structs = np.unique(annotation_slice)
    for sid in unique_structs:
        if sid == 0:
            continue
        mask = (annotation_slice == sid).astype(float)
        ax.contour(mask, colors="grey", levels=[0.5], linewidths=0.8, alpha=0.2)

    # Build networkx graph
    G = nx.Graph()
    pos = {}
    texts = []

    for _, row in nodes_df.iterrows():
        acr = row["region"]
        rid = acronym_to_id.get(acr)
        if rid and rid in region_centers:
            center = np.array(region_centers[rid])
            # Original transformation:
            #   center = np.flipud(np.transpose(center))
            #   center[1] = annotation.shape[1] - center[1]
            # center is [DV, AP] from _get_region_centers_3d
            # flipud(transpose) on a 1D array: effectively swaps the two values -> [AP, DV]
            center = np.flipud(np.transpose(center))  # [AP, DV] -> swap -> [DV, AP] ... 
            # Actually for a 1D array, transpose is no-op, flipud reverses: [DV, AP] -> [AP, DV]
            # So center is now [AP, DV]
            center[1] = annotation.shape[1] - center[1]  # flip DV axis
            
            cx, cy = center[0], center[1]
            print(f"  center for: {acr} ({cx:.1f}, {cy:.1f})")
            
            G.add_node(acr, size=row.get(node_size_col, 50),
                       color=row.get(node_color_col, "grey"))
            pos[acr] = (cx, cy)
        else:
            print(f"  Warning: '{acr}' not found in atlas")

    # Normalize edge weights for line width
    if edge_weight_col in edges_df.columns:
        edge_weights = edges_df[edge_weight_col].values
        scaler = MinMaxScaler(feature_range=(2, 12))
        norm_weights = scaler.fit_transform(
            np.abs(edge_weights).reshape(-1, 1)
        ).flatten()
    else:
        norm_weights = np.full(len(edges_df), 4.0)

    # Edge colormap
    import matplotlib.cm as cm_mod
    import matplotlib.colors as mcolors
    if edge_weight_col in edges_df.columns:
        w_vals = edges_df[edge_weight_col].values
        enorm = mcolors.Normalize(vmin=w_vals.min(), vmax=w_vals.max())
        ecmap = cm_mod.Greys
    else:
        enorm = mcolors.Normalize(0, 1)
        ecmap = cm_mod.Greys

    for idx, (_, row) in enumerate(edges_df.iterrows()):
        src = str(row.get("source", row.get("region", "")))
        tgt = str(row.get("target", row.get("region2", "")))
        if src in pos and tgt in pos:
            w = row.get(edge_weight_col, 0.5)
            if abs(w) < min_edge_weight:
                continue
            color = ecmap(enorm(abs(w)))
            G.add_edge(src, tgt, weight=norm_weights[idx], color=color)

    # Draw edges
    for u, v, d in G.edges(data=True):
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        ax.plot([x1, x2], [y1, y2], color=d.get("color", "grey"),
                linewidth=d.get("weight", 3) * 0.8, alpha=0.5, zorder=2)

    # Draw nodes (using networkx for consistency with original)
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_size=[G.nodes[n].get("size", 50) * 10 for n in G.nodes],
        node_color=[MODULE_COLORS.get(str(G.nodes[n].get("color", "grey")),
                                       str(G.nodes[n].get("color", "grey")))
                    for n in G.nodes],
        alpha=0.8
    )
    nx.draw_networkx_edges(
        G, pos, ax=ax,
        width=[G.edges[e].get("weight", 3) * 0.8 for e in G.edges],
        edge_color=[G.edges[e].get("color", "grey") for e in G.edges],
        alpha=0.5
    )

    # Add text labels
    for node in G.nodes:
        if node not in pos:
            continue
        x, y = pos[node]
        txt = ax.text(x, y, node, fontsize=14, ha="center", va="center",
                      fontweight="bold", zorder=5,
                      path_effects=[pe.withStroke(linewidth=3, foreground="white")])
        texts.append(txt)

    # Adjust text to avoid overlaps
    if has_adjust and texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))

    ax.axis("off")

    if title is None:
        title = f"Network on Brain Atlas at {plane.title()} Plane z={slice_idx}"
    ax.set_title(title, fontsize=20, fontweight="bold", pad=20)

    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved atlas network: {output_path}")


# =============================================================================
# 1b. FALLBACK: Simple brain silhouette (no AllenSDK required)
# =============================================================================
DEFAULT_REGION_COORDS = {
    "CTXpl": (0.25, 0.75), "CTXsp": (0.30, 0.60), "STR": (0.35, 0.55),
    "PAL": (0.38, 0.50), "TH": (0.45, 0.48), "HY": (0.42, 0.35),
    "MB": (0.55, 0.55), "P": (0.62, 0.50), "MY": (0.70, 0.45),
    "CB": (0.72, 0.65), "HPF": (0.40, 0.65), "OLF": (0.18, 0.45),
    "Isocortex": (0.28, 0.80), "CNU": (0.36, 0.52), "IB": (0.45, 0.45),
    "HB": (0.65, 0.48), "CBX": (0.73, 0.68),
}


def plot_atlas_network_simple(
    nodes_df, edges_df, output_path,
    region_coords=None,
    title="WGCNA Module Network",
    node_size_col="n_genes",
    node_color_col="module_color",
    edge_weight_col="correlation",
    min_edge_weight=0.3,
    figsize=(14, 10),
):
    """
    Fallback atlas network using hand-drawn brain silhouette.
    Used when AllenSDK is not available.
    """
    if region_coords is None:
        region_coords = DEFAULT_REGION_COORDS

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")

    # Draw brain sagittal silhouette
    cerebrum = Ellipse((0.38, 0.62), 0.55, 0.45, angle=-8,
                       facecolor="#F5F5F0", edgecolor="#B0B0B0",
                       linewidth=2.5, alpha=0.4, zorder=0)
    ax.add_patch(cerebrum)
    cerebellum = Ellipse((0.73, 0.62), 0.18, 0.22, angle=10,
                         facecolor="#F0F0E8", edgecolor="#B0B0B0",
                         linewidth=2.0, alpha=0.4, zorder=0)
    ax.add_patch(cerebellum)
    brainstem = Ellipse((0.62, 0.38), 0.30, 0.15, angle=-30,
                        facecolor="#EEEEEA", edgecolor="#B0B0B0",
                        linewidth=2.0, alpha=0.4, zorder=0)
    ax.add_patch(brainstem)

    import networkx as nx
    G = nx.Graph()
    pos = {}

    for _, row in nodes_df.iterrows():
        acr = row["region"]
        if acr in region_coords:
            G.add_node(acr, size=row.get(node_size_col, 50),
                       color=row.get(node_color_col, "grey"))
            pos[acr] = region_coords[acr]

    for _, row in edges_df.iterrows():
        src = str(row.get("source", row.get("region", "")))
        tgt = str(row.get("target", row.get("region2", "")))
        if src in pos and tgt in pos:
            w = abs(row.get(edge_weight_col, 0.5))
            if w >= min_edge_weight:
                G.add_edge(src, tgt, weight=w)

    nx.draw_networkx_edges(G, pos, ax=ax, width=2, alpha=0.4, edge_color="grey")
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_size=[G.nodes[n].get("size", 50) * 8 for n in G.nodes],
        node_color=[MODULE_COLORS.get(str(G.nodes[n].get("color", "grey")),
                                       str(G.nodes[n].get("color", "grey")))
                    for n in G.nodes],
        alpha=0.85, edgecolors="white", linewidths=2
    )
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=13, font_weight="bold")

    ax.set_title(title, fontsize=18, fontweight="bold", pad=15)
    ax.axis("off")
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved simple atlas network: {output_path}")


# =============================================================================
# 2. CIRCOS PLOT
# =============================================================================
def plot_circos(modules_df, output_path, figsize=(14, 14)):
    """
    Multi-ring circos plot for WGCNA module enrichment.
    Outer ring = module color, inner rings = enrichment values.
    """
    ring_cols = [c for c in modules_df.columns
                 if c not in ("module", "module_color", "n_genes")]
    n_rings = len(ring_cols)
    n_mod = len(modules_df)
    if n_mod == 0:
        print("No modules to plot"); return

    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw={"polar": True})
    angles = np.linspace(0, 2 * np.pi, n_mod, endpoint=False)
    width = 2 * np.pi / n_mod * 0.85

    r_outer = 1.0
    ring_width = 0.10
    gap = 0.02

    # Module color ring (outermost)
    for i, (_, row) in enumerate(modules_df.iterrows()):
        color = MODULE_COLORS.get(str(row.get("module_color", "grey")), "#CCCCCC")
        ax.bar(angles[i], ring_width, width=width, bottom=r_outer,
               color=color, edgecolor="white", linewidth=0.5, zorder=3)

    # Module name labels
    for i, (_, row) in enumerate(modules_df.iterrows()):
        label_r = r_outer + ring_width + 0.06
        angle_deg = np.degrees(angles[i])
        rotation = angle_deg - 90 if angle_deg < 180 else angle_deg + 90
        ha = "left" if angle_deg < 180 else "right"
        ax.text(angles[i], label_r, str(row["module"]),
                ha=ha, va="center", fontsize=11, fontweight="bold",
                rotation=rotation, rotation_mode="anchor")

    # Enrichment rings (inside)
    cmap_pos = LinearSegmentedColormap.from_list("pos", ["#FFFFFF", "#D62728"])
    cmap_neg = LinearSegmentedColormap.from_list("neg", ["#1F77B4", "#FFFFFF"])

    for ri, col in enumerate(ring_cols):
        r_base = r_outer - (ri + 1) * (ring_width + gap)
        vals = modules_df[col].fillna(0).values
        vmax = max(abs(vals.max()), abs(vals.min()), 1e-9)
        for i, v in enumerate(vals):
            if v >= 0:
                color = cmap_pos(min(v / vmax, 1.0))
            else:
                color = cmap_neg(min(abs(v) / vmax, 1.0))
            ax.bar(angles[i], ring_width, width=width, bottom=r_base,
                   color=color, edgecolor="white", linewidth=0.3, zorder=2)

    # Ring labels
    for ri, col in enumerate(ring_cols):
        r_label = r_outer - (ri + 1) * (ring_width + gap) + ring_width / 2
        label = col.replace("_", " ").title()
        ax.text(np.pi, r_label, label, ha="center", va="center",
                fontsize=9, fontstyle="italic", color="#555555")

    ax.set_ylim(0, r_outer + ring_width + 0.25)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["polar"].set_visible(False)
    ax.set_title("WGCNA Module Enrichment", fontsize=18, fontweight="bold",
                 pad=25, y=1.05)

    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved circos: {output_path}")


# =============================================================================
# 3. HUB GENE SUBNETWORK
# =============================================================================
def plot_hub_network(hubs_df, edges_df, output_path, figsize=(14, 14)):
    """
    Hub gene subnetwork: gold hub nodes + sky blue neighbor nodes.
    """
    import networkx as nx

    G = nx.Graph()
    hub_genes = set(hubs_df["gene"].unique()) if "gene" in hubs_df.columns else set()

    for _, row in edges_df.iterrows():
        src = str(row.get("source", row.get("gene1", "")))
        tgt = str(row.get("target", row.get("gene2", "")))
        w = float(row.get("weight", row.get("correlation", 0.5)))
        if src and tgt:
            G.add_edge(src, tgt, weight=abs(w))

    if len(G) == 0:
        print("Empty graph"); return

    pos = nx.spring_layout(G, k=2.5 / np.sqrt(max(len(G), 1)),
                           iterations=80, seed=42)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    node_colors = ["#FFD700" if n in hub_genes else "#87CEEB" for n in G.nodes]
    node_sizes = [800 if n in hub_genes else 300 for n in G.nodes]

    edge_weights = [G.edges[e].get("weight", 0.5) for e in G.edges]
    max_w = max(edge_weights) if edge_weights else 1
    edge_widths = [1 + 4 * (w / max_w) for w in edge_weights]

    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                           alpha=0.3, edge_color="#999999")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color=node_colors, edgecolors="white",
                           linewidths=2, alpha=0.9)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=11, font_weight="bold")

    hub_patch = Patch(facecolor="#FFD700", edgecolor="white", label=f"Hub genes ({len(hub_genes)})")
    nbr_patch = Patch(facecolor="#87CEEB", edgecolor="white",
                      label=f"Neighbors ({len(G) - len(hub_genes)})")
    edge_patch = Patch(facecolor="none", edgecolor="#999999",
                       label=f"Edges ({G.number_of_edges()})")
    ax.legend(handles=[hub_patch, nbr_patch, edge_patch],
              loc="upper left", fontsize=13, framealpha=0.9)

    ax.set_title("Hub Gene Network", fontsize=18, fontweight="bold", pad=15)
    ax.axis("off")
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved hub network: {output_path}")


# =============================================================================
# CLI
# =============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CellDynamicST WGCNA Network Visualization"
    )
    sub = parser.add_subparsers(dest="cmd")

    # atlas
    p_atlas = sub.add_parser("atlas", help="Brain atlas overlay network")
    p_atlas.add_argument("--nodes", required=True, help="Nodes CSV (region, n_genes, module_color)")
    p_atlas.add_argument("--edges", required=True, help="Edges CSV (source, target, correlation)")
    p_atlas.add_argument("--output", required=True)
    p_atlas.add_argument("--plane", default="sagittal", choices=["sagittal", "coronal", "horizontal"])
    p_atlas.add_argument("--slice", type=int, default=None, help="Slice index along chosen axis")
    p_atlas.add_argument("--resolution", type=int, default=10)
    p_atlas.add_argument("--title", default=None)
    p_atlas.add_argument("--min-edge-weight", type=float, default=0.0)

    # circos
    p_circ = sub.add_parser("circos", help="Circos enrichment plot")
    p_circ.add_argument("--modules", required=True, help="Modules CSV")
    p_circ.add_argument("--output", required=True)

    # hub
    p_hub = sub.add_parser("hub", help="Hub gene subnetwork")
    p_hub.add_argument("--hubs", required=True, help="Hub genes CSV (gene column)")
    p_hub.add_argument("--edges", required=True, help="Edges CSV (source, target, weight)")
    p_hub.add_argument("--output", required=True)

    args = parser.parse_args()

    if args.cmd == "atlas":
        nodes = pd.read_csv(args.nodes)
        edges = pd.read_csv(args.edges)
        try:
            plot_atlas_network(
                nodes, edges, args.output,
                plane=args.plane,
                slice_idx=args.slice,
                resolution=args.resolution,
                title=args.title,
                min_edge_weight=args.min_edge_weight,
            )
        except ImportError as e:
            print(f"AllenSDK not available: {e}")
            print("Falling back to simple brain silhouette...")
            plot_atlas_network_simple(nodes, edges, args.output)
    elif args.cmd == "circos":
        modules = pd.read_csv(args.modules)
        plot_circos(modules, args.output)
    elif args.cmd == "hub":
        hubs = pd.read_csv(args.hubs)
        edges = pd.read_csv(args.edges)
        plot_hub_network(hubs, edges, args.output)
    else:
        parser.print_help()
