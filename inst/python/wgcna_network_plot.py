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


def _get_region_centers(annotation, structure_tree, plane="sagittal", slice_idx=None):
    """
    Compute 2D centroid of each brain structure in the chosen slice plane.
    For parent regions, aggregates all descendant structure voxels to compute
    the centroid.

    Parameters
    ----------
    annotation : 3D ndarray (AP, DV, ML)
    structure_tree : AllenSDK StructureTree
    plane : str  'sagittal' | 'coronal' | 'horizontal'
    slice_idx : int  index along the slicing axis (None = midpoint)

    Returns
    -------
    region_centers : dict  {structure_id: np.array([x, y])}
    slice_2d : 2D ndarray  the annotation slice for contour drawing
    """
    shape = annotation.shape  # (AP, DV, ML)
    if plane == "sagittal":
        if slice_idx is None:
            slice_idx = shape[2] // 2
        slice_idx = min(max(0, slice_idx), shape[2] - 1)
        slice_2d = annotation[:, :, slice_idx]  # (AP, DV)
    elif plane == "coronal":
        if slice_idx is None:
            slice_idx = shape[0] // 2
        slice_idx = min(max(0, slice_idx), shape[0] - 1)
        slice_2d = annotation[slice_idx, :, :]  # (DV, ML)
    elif plane == "horizontal":
        if slice_idx is None:
            slice_idx = shape[1] // 2
        slice_idx = min(max(0, slice_idx), shape[1] - 1)
        slice_2d = annotation[:, slice_idx, :]  # (AP, ML)
    else:
        raise ValueError(f"Unknown plane '{plane}'. Use sagittal/coronal/horizontal.")

    # Build a lookup: leaf_id -> set of all ancestor IDs (including self)
    # This lets us compute centroids for parent regions by aggregating descendants
    unique_ids_in_slice = set(int(x) for x in np.unique(slice_2d) if x != 0)

    # Compute centroids for each unique leaf-level structure
    leaf_centers_raw = {}  # {leaf_id: coords_array}
    for sid in unique_ids_in_slice:
        coords = np.argwhere(slice_2d == sid)
        if coords.size > 0:
            leaf_centers_raw[sid] = coords

    # Also build parent -> descendant mapping for requested regions
    region_centers = {}
    # Store leaf-level centroids too
    for sid, coords in leaf_centers_raw.items():
        center = coords.mean(axis=0)
        region_centers[sid] = np.array([center[1], slice_2d.shape[0] - center[0]])

    return region_centers, slice_2d, slice_idx, leaf_centers_raw


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

    # Get region centers and 2D slice
    region_centers, slice_2d, actual_slice, leaf_raw = _get_region_centers(
        annotation, structure_tree, plane=plane, slice_idx=slice_idx
    )
    print(f"  Plane: {plane}, slice index: {actual_slice}")
    print(f"  Found {len(region_centers)} leaf structures in slice")

    # Map acronyms to IDs
    all_acronyms = set(nodes_df["region"].unique())
    if "source" in edges_df.columns:
        all_acronyms |= set(edges_df["source"].unique())
    if "target" in edges_df.columns:
        all_acronyms |= set(edges_df["target"].unique())
    acronym_to_id = _map_acronyms_to_ids(structure_tree, all_acronyms)
    print(f"  Mapped {len(acronym_to_id)}/{len(all_acronyms)} acronyms to IDs")

    # For parent regions, compute centroid from all descendant voxels in slice
    unique_ids_in_slice = set(leaf_raw.keys())
    for acr, sid in acronym_to_id.items():
        if sid in region_centers:
            continue  # already a leaf-level match
        # Get all descendant IDs
        try:
            desc_ids = set(structure_tree.descendant_ids([sid])[0])
        except Exception:
            desc_ids = {sid}
        # Collect all voxel coordinates from descendants present in slice
        all_coords = []
        for did in desc_ids:
            if did in leaf_raw:
                all_coords.append(leaf_raw[did])
        if all_coords:
            merged = np.vstack(all_coords)
            center = merged.mean(axis=0)
            region_centers[sid] = np.array([center[1], slice_2d.shape[0] - center[0]])
            print(f"  Aggregated {len(all_coords)} sub-regions for '{acr}' (id={sid})")
        else:
            print(f"  No descendant voxels for '{acr}' in this slice")

    # Build figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Draw brain structure contours from annotation slice
    flipped_slice = np.flipud(slice_2d.T) if plane == "sagittal" else slice_2d
    unique_structs = np.unique(flipped_slice)
    for sid in unique_structs:
        if sid == 0:
            continue
        mask = (flipped_slice == sid).astype(float)
        ax.contour(mask, colors="grey", levels=[0.5], linewidths=0.8, alpha=0.2)

    # Build networkx graph
    G = nx.Graph()
    pos = {}
    texts = []

    for _, row in nodes_df.iterrows():
        acr = row["region"]
        rid = acronym_to_id.get(acr)
        if rid and rid in region_centers:
            center = region_centers[rid]
            # For sagittal, transform coordinates to match flipped slice
            if plane == "sagittal":
                cx = center[0]
                cy = flipped_slice.shape[0] - center[1]
            else:
                cx, cy = center[0], center[1]
            G.add_node(acr, size=row.get(node_size_col, 50),
                       color=row.get(node_color_col, "grey"))
            pos[acr] = (cx, cy)
        else:
            print(f"  Warning: '{acr}' not found in this slice")

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
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    if edge_weight_col in edges_df.columns:
        w_vals = edges_df[edge_weight_col].values
        enorm = mcolors.Normalize(vmin=w_vals.min(), vmax=w_vals.max())
        ecmap = cm.Greys
    else:
        enorm = mcolors.Normalize(0, 1)
        ecmap = cm.Greys

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

    # Draw nodes
    for node in G.nodes:
        if node not in pos:
            continue
        x, y = pos[node]
        size = G.nodes[node].get("size", 50) * 10
        color_name = str(G.nodes[node].get("color", "grey"))
        color = MODULE_COLORS.get(color_name, color_name)
        ax.scatter(x, y, s=size, c=color, edgecolors="white",
                   linewidths=2.5, zorder=4, alpha=0.85)
        txt = ax.text(x, y, node, fontsize=16, ha="center", va="center",
                      fontweight="bold", zorder=5,
                      path_effects=[pe.withStroke(linewidth=3, foreground="white")])
        texts.append(txt)

    # Adjust text to avoid overlaps
    if has_adjust and texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))

    if title is None:
        title = f"Network on Brain Atlas — {plane.title()} Plane (slice={actual_slice})"
    ax.set_title(title, fontsize=20, fontweight="bold", pad=15)
    ax.axis("off")

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
    olf = Ellipse((0.12, 0.50), 0.10, 0.08, angle=0,
                  facecolor="#F0F0E8", edgecolor="#B0B0B0",
                  linewidth=1.5, alpha=0.4, zorder=0)
    ax.add_patch(olf)

    # Edge colormap
    cmap_edge = LinearSegmentedColormap.from_list(
        "edge_cmap", ["#3B4CC0", "#F7F7F7", "#B40426"]
    )
    if edge_weight_col in edges_df.columns:
        edges_plot = edges_df[edges_df[edge_weight_col].abs() >= min_edge_weight].copy()
        vmin = edges_plot[edge_weight_col].min()
        vmax = edges_plot[edge_weight_col].max()
        if vmin >= 0: vmin = -0.01
        if vmax <= 0: vmax = 0.01
        edge_norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    else:
        edges_plot = edges_df.copy()
        edge_norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

    for _, row in edges_plot.iterrows():
        src = str(row["source"])
        tgt = str(row["target"])
        if src not in region_coords or tgt not in region_coords:
            continue
        x1, y1 = region_coords[src]
        x2, y2 = region_coords[tgt]
        weight = row.get(edge_weight_col, 0.5)
        color = cmap_edge(edge_norm(weight))
        lw = abs(weight) * 5 + 0.5
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw,
                alpha=0.5, zorder=1, solid_capstyle="round")

    for _, row in nodes_df.iterrows():
        region = str(row.get("region", ""))
        if region not in region_coords:
            continue
        x, y = region_coords[region]
        size = row.get(node_size_col, 50)
        color_name = str(row.get(node_color_col, "grey"))
        color = MODULE_COLORS.get(color_name, color_name)
        radius = np.clip(np.sqrt(size) / 40, 0.02, 0.07)
        circle = Circle((x, y), radius, facecolor=color,
                        edgecolor="white", linewidth=2.5, zorder=3)
        ax.add_patch(circle)
        ax.text(x, y - radius - 0.025, region, fontsize=14,
                ha="center", va="top", fontweight="bold",
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])

    sm = plt.cm.ScalarMappable(cmap=cmap_edge, norm=edge_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.02, aspect=20)
    cbar.set_label("Module Correlation", fontsize=15, fontweight="bold")

    ax.set_title(title, fontsize=20, fontweight="bold", pad=15)
    ax.axis("off")
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved atlas network (simple): {output_path}")


# =============================================================================
# 2. CIRCOS PLOT: Multi-ring enrichment (matching manuscript Fig 4B)
# =============================================================================
def plot_circos(
    module_df, output_path,
    comparison_cols=None,
    module_col="module_color",
    enrichment_col_prefix="NES_",
    title="Consensus Module Enrichment Circos Plot",
    figsize=(12, 12),
    ring_colormaps=None,
):
    """
    Create a multi-ring circos plot where each ring represents a different
    enrichment type (OUD, DEG effects), with the outermost ring showing
    module colors.
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw={"projection": "polar"})

    if comparison_cols is None:
        comparison_cols = [c for c in module_df.columns
                          if c.startswith(enrichment_col_prefix)]

    n_modules = len(module_df)
    n_rings = len(comparison_cols)
    if n_modules == 0 or n_rings == 0:
        print("No data for circos plot.")
        return

    if ring_colormaps is None:
        default_cmaps = [plt.cm.Oranges, plt.cm.Blues, plt.cm.Purples,
                         plt.cm.Greens, plt.cm.Reds]
        ring_colormaps = [default_cmaps[i % len(default_cmaps)]
                          for i in range(n_rings)]

    for col in comparison_cols:
        module_df[col + "_norm"] = _norm(module_df[col])

    angles = np.linspace(0, 2 * np.pi, n_modules + 1)
    ring_base = 0.85
    ring_width = 0.06
    radii = [ring_base + ring_width * i for i in range(n_rings + 1)]

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    ax.set_ylim(0, radii[-1] + 0.20)

    for j, (comp_col, cmap) in enumerate(zip(comparison_cols, ring_colormaps)):
        for idx in range(n_modules):
            row = module_df.iloc[idx]
            theta0, theta1 = angles[idx], angles[idx + 1]
            norm_val = row[comp_col + "_norm"]
            color = cmap(norm_val)
            ax.bar(
                x=(theta0 + theta1) / 2,
                height=radii[j + 1] - radii[j],
                width=theta1 - theta0,
                bottom=radii[j],
                color=color, linewidth=0, align="center"
            )

    for idx in range(n_modules):
        row = module_df.iloc[idx]
        theta0, theta1 = angles[idx], angles[idx + 1]
        color_name = str(row.get(module_col, "grey"))
        color = MODULE_COLORS.get(color_name, color_name)
        ax.bar(
            x=(theta0 + theta1) / 2,
            height=0.055,
            width=theta1 - theta0,
            bottom=radii[-1],
            color=color, edgecolor="black", linewidth=1.2, align="center"
        )

    for idx in range(n_modules):
        row = module_df.iloc[idx]
        theta0, theta1 = angles[idx], angles[idx + 1]
        mid_angle = (theta0 + theta1) / 2
        label = str(row.get(module_col, f"M{idx}"))
        rotation = (np.degrees(mid_angle) + 270) % 360
        if 90 < rotation < 270:
            rotation += 180
        ax.text(mid_angle, radii[-1] + 0.09, label,
                fontsize=12, fontweight="bold",
                ha="center", va="center",
                rotation=rotation, rotation_mode="anchor")

    ring_labels = [c.replace(enrichment_col_prefix, "").replace("_", " ")
                   for c in comparison_cols]
    legend_patches = []
    for lbl, cmap in zip(ring_labels, ring_colormaps):
        legend_patches.append(Patch(facecolor=cmap(0.85), edgecolor="k", label=lbl))
    legend_patches.append(Patch(facecolor="#888888", edgecolor="k", label="Module Color"))
    ax.legend(handles=legend_patches, bbox_to_anchor=(0.5, 1.08),
              loc="lower center", ncol=min(3, len(legend_patches)),
              fontsize=13, frameon=False)

    subtitle = "(inner\u2192outer: " + ", ".join(ring_labels) + ", Module Color)"
    ax.set_title(title + "\n" + subtitle, fontsize=18, fontweight="bold",
                 pad=20, y=1.15)

    plt.savefig(output_path, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved circos plot: {output_path}")


# =============================================================================
# 3. HUB GENE SUBNETWORK (matching manuscript Fig 4F-I)
# =============================================================================
def plot_hub_network(
    hub_df, edge_df, output_path,
    gene_col="gene",
    module_col="module_color",
    kme_col="kME",
    is_hub_col="is_hub",
    source_col="gene1",
    target_col="gene2",
    weight_col="weight",
    top_n=50,
    title="Hub Gene Subnetwork",
    figsize=(14, 14),
):
    """
    Create a force-directed hub gene subnetwork with gold hub nodes and
    sky blue neighbor nodes.
    """
    try:
        import networkx as nx
    except ImportError:
        print("networkx required. Install with: pip install networkx")
        return

    hub_df = hub_df.nlargest(top_n, kme_col).copy()

    if is_hub_col in hub_df.columns:
        hub_genes = set(hub_df[hub_df[is_hub_col].astype(str).isin(
            ["TRUE", "True", "true", "1", "yes"])][gene_col])
    else:
        n_hubs = max(3, int(len(hub_df) * 0.2))
        hub_genes = set(hub_df.nlargest(n_hubs, kme_col)[gene_col])

    all_genes = set(hub_df[gene_col])
    edge_df = edge_df[
        edge_df[source_col].isin(all_genes) & edge_df[target_col].isin(all_genes)
    ].copy()

    G = nx.Graph()
    for _, row in hub_df.iterrows():
        gene = row[gene_col]
        is_hub = gene in hub_genes
        G.add_node(gene, is_hub=is_hub, kme=row[kme_col])

    for _, row in edge_df.iterrows():
        if row[source_col] in G and row[target_col] in G:
            G.add_edge(row[source_col], row[target_col],
                       weight=abs(row.get(weight_col, 0.5)))

    if len(G.nodes) == 0:
        print("No nodes in hub network.")
        return

    isolated = [n for n in G.nodes if G.degree(n) == 0]
    G.remove_nodes_from(isolated)
    if len(G.nodes) == 0:
        print("No connected nodes in hub network.")
        return

    pos = nx.spring_layout(G, k=2.5 / np.sqrt(max(len(G.nodes), 1)),
                           iterations=150, seed=42)

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    fig.patch.set_facecolor("white")

    max_w = max((d.get("weight", 0.5) for _, _, d in G.edges(data=True)), default=1)
    for u, v, d in G.edges(data=True):
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        w = d.get("weight", 0.5)
        lw = 1.0 + 4.0 * abs(w / max(max_w, 1e-9))
        ax.plot([x1, x2], [y1, y2], color="#808080",
                linewidth=lw, alpha=0.35, zorder=1)

    for node in G.nodes:
        x, y = pos[node]
        is_hub = G.nodes[node].get("is_hub", False)
        kme = G.nodes[node].get("kme", 0.5)

        if is_hub:
            color = "#DAA520"
            size = 350 + kme * 500
            edge_color = "#B8860B"
            fontsize = 13
            fontweight = "bold"
        else:
            color = "#87CEEB"
            size = 120 + kme * 200
            edge_color = "#4682B4"
            fontsize = 11
            fontweight = "normal"

        ax.scatter(x, y, s=size, c=color, edgecolors=edge_color,
                   linewidths=2.0, zorder=3)
        ax.text(x, y + 0.035, node, fontsize=fontsize, ha="center",
                va="bottom", fontweight=fontweight, color="black",
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#DAA520",
               markeredgecolor="#B8860B", markersize=14, label="Hub gene"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#87CEEB",
               markeredgecolor="#4682B4", markersize=10, label="Neighbor gene"),
    ]
    ax.legend(handles=legend_elements, loc="lower left", fontsize=15,
              frameon=True, fancybox=True, shadow=True)

    n_nodes = len(G.nodes)
    n_edges = len(G.edges)
    n_hubs_shown = sum(1 for n in G.nodes if G.nodes[n].get("is_hub", False))
    ax.set_title(f"{title}\n({n_hubs_shown} hubs, {n_nodes - n_hubs_shown} neighbors, "
                 f"{n_edges} edges)", fontsize=18, fontweight="bold")
    ax.axis("off")

    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved hub network: {output_path}")


# =============================================================================
# CLI ENTRY POINT
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="CellDynamicST: WGCNA Network Visualization"
    )
    sub = parser.add_subparsers(dest="command")

    # Atlas network
    atlas = sub.add_parser("atlas", help="Brain atlas overlay network")
    atlas.add_argument("--nodes", required=True, help="Nodes CSV")
    atlas.add_argument("--edges", required=True, help="Edges CSV")
    atlas.add_argument("--output", default="wgcna_atlas_network.png")
    atlas.add_argument("--title", default=None)
    atlas.add_argument("--plane", default="sagittal",
                       choices=["sagittal", "coronal", "horizontal"],
                       help="Atlas slice plane (default: sagittal)")
    atlas.add_argument("--slice", type=int, default=None,
                       help="Slice index along the chosen axis (default: midpoint)")
    atlas.add_argument("--resolution", type=int, default=10,
                       help="Allen annotation resolution in microns (10 or 25)")
    atlas.add_argument("--simple", action="store_true",
                       help="Use simple silhouette fallback (no AllenSDK)")

    # Circos
    circos = sub.add_parser("circos", help="Module enrichment circos plot")
    circos.add_argument("--modules", required=True, help="Module enrichment CSV")
    circos.add_argument("--output", default="wgcna_circos.png")
    circos.add_argument("--title", default="Consensus Module Enrichment Circos")

    # Hub network
    hub = sub.add_parser("hub", help="Hub gene force-directed network")
    hub.add_argument("--hubs", required=True, help="Hub genes CSV")
    hub.add_argument("--edges", required=True, help="Gene-gene edges CSV")
    hub.add_argument("--output", default="wgcna_hub_network.png")
    hub.add_argument("--top-n", type=int, default=50)
    hub.add_argument("--title", default="Hub Gene Subnetwork")

    args = parser.parse_args()

    if args.command == "atlas":
        nodes = pd.read_csv(args.nodes)
        edges = pd.read_csv(args.edges)
        if args.simple:
            plot_atlas_network_simple(nodes, edges, args.output, title=args.title or "WGCNA Module Network")
        else:
            plot_atlas_network(
                nodes, edges, args.output,
                plane=args.plane,
                slice_idx=args.slice,
                resolution=args.resolution,
                title=args.title,
            )

    elif args.command == "circos":
        modules = pd.read_csv(args.modules)
        plot_circos(modules, args.output, title=args.title)

    elif args.command == "hub":
        hubs = pd.read_csv(args.hubs)
        edges = pd.read_csv(args.edges)
        plot_hub_network(hubs, edges, args.output,
                        top_n=args.top_n, title=args.title)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
