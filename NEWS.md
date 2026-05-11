# CellDynamicST 1.0.0

## Initial Release

### New Features

* **Generalized Clustering** (`cdst_cluster`): Reference-based clustering with Random Forest label propagation, including automatic technical bias PC detection and exclusion.
* **Neurotransmitter Classification** (`cdst_classify_nt`): GMM-based neurotransmitter type classification supporting 8 NT types with species-configurable marker gene lists.
* **Glial Classification** (`cdst_classify_glia`): Dual-module GMM classification for astrocytes, microglia, oligodendrocytes, OPCs, ependymal, and endothelial cells.
* **Hierarchical Annotation** (`cdst_annotate_hierarchy`): Multi-level cell type hierarchy construction from cluster to superclass.
* **Marker Validation** (`cdst_validate_markers`): Cross-validation against CellMarker2 database.
* **Atlas Visualization** (`cdst_plot_atlas`): Brain atlas projection via Allen CCFv3 with Python backend.
* **Spatial Cell Dynamics** (`cdst_plot_heatmap`, `cdst_plot_sankey`, `cdst_plot_dimred`): Distribution heatmaps, Sankey diagrams, and dimensionality reduction plots.
* **DEG Analysis** (`cdst_compute_deg`): Differential expression analysis with configurable pairwise comparisons.
* **GO/GSEA Enrichment** (`cdst_compute_enrichment`): Gene Ontology and Gene Set Enrichment Analysis.
* **Spatial Enrichment** (`cdst_compute_spatial_enrichment`): Single-cell level spatial GO enrichment via weighted GSEA.
* **WGCNA** (`cdst_detect_modules`): Consensus and per-group weighted gene co-expression network analysis with automatic soft power selection.
* **Hub Gene Detection** (`cdst_detect_hubs`): ARACNe/minet-based hub gene identification.
* **Differential Networks** (`cdst_compute_diff_network`): Differential co-expression network analysis between experimental conditions.
* **Full Pipeline** (`cdst_run_pipeline`): One-command execution of the complete analysis pipeline.

### Infrastructure

* YAML-based configuration system with sensible defaults.
* S4 class design (`CdstConfig`, `CdstResult`) with type-safe validation.
* Comprehensive input validation with informative error messages via `checkmate` and `cli`.
* Provenance tracking for FAIR4RS compliance.
* GitHub Actions CI/CD for multi-OS R CMD check and code coverage.
* Docker images for reproducible environments.
* pkgdown documentation website.
