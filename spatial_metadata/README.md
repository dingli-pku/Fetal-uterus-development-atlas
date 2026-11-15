# Uterus spatial clustering metadata

These tables contain cell-level metadata for spatial clustering of uterus samples. The BMK files cover stromal spatial annotations for the INN and smooth muscle cell (SMC) analysis, with epithelial subsets and combined epithelial–stromal versions where available. The STAGATE files contain broad cell-group annotations from the merged spatial datasets. In file names, `19w` means gestational week 19, and so on.

## BMK spatial clustering

Source directory on the analysis server: `/Data6/chenzl/uterus/inn_new_dl/BMK_stromal/`.

- [`BMK_stromal/`](BMK_stromal/) contains seven existing tab-delimited `.txt` exports copied from the server. These tables include the source `id` column.
- [`BMK_stromal/from_RData/`](BMK_stromal/from_RData/) contains 11 TSV files extracted from the Seurat objects in the source directory. Each TSV preserves all columns in the object's `meta.data` and adds `cell_id` from its row names.

Each item below is a source file stem. Append `__metadata.tsv` to get the exported file name.

| Week | File stems and cell counts |
| --- | --- |
| 19 | `location_19w` (28,101), `location_19w_epi` (2,550), `location_19w_stromal_epi` (30,651) |
| 20 | `location_20w` (27,653) |
| 22 | `location_22w` (65,982), `location_22w_epi` (4,220), `location_22w_epi_stromal` (70,202) |
| 26 | `location_26w` (151,820), `location_26w_epi` (3,066), `location_26w_epi_stromal` (151,820), `location_26w_stromal` (148,754) |

The suffixes `_epi` and `_stromal` identify epithelial and stromal subsets. Both `_epi_stromal` and `_stromal_epi` denote combined versions. Unsuffixed files retain the original object naming. Files for the same week can contain overlapping cells, so use the version that matches the intended subset before counting cells.

## STAGATE broad cell-group clustering

[`STAGATE/`](STAGATE/) contains metadata exported from same-named RData copies under `/Data6/chenzl/uterus/new_cellbin/`. The original project path was unavailable on the analysis server, so equivalence to the original RData files was unverified.

| Week | Metadata file | Cells |
| --- | --- | ---: |
| 19 | [`cellNewBin_19w_merge_clean__metadata.tsv`](STAGATE/cellNewBin_19w_merge_clean__metadata.tsv) | 44,081 |
| 20 | [`cellNewBin_20w_merge_clean__metadata.tsv`](STAGATE/cellNewBin_20w_merge_clean__metadata.tsv) | 51,996 |
| 22 | [`cellNewBin_22w_merge_clean__metadata.tsv`](STAGATE/cellNewBin_22w_merge_clean__metadata.tsv) | 98,853 |
| 26 | [`cellNewBin_26w_merge_clean__metadata.tsv`](STAGATE/cellNewBin_26w_merge_clean__metadata.tsv) | 209,246 |

## Reading the tables

Each row in the new TSVs represents one cell in the source Seurat object. Use `cell_id` to match rows by Seurat cell name within the corresponding sample. Where present, `id` is a separate annotation field carried over from the source. The available columns vary across weeks and objects.

| Column | Meaning |
| --- | --- |
| `orig.ident` | Original Seurat sample or project identifier |
| `x`, `y` | Spatial coordinates from the source object |
| `nCount_RNA`, `nFeature_RNA` | Total RNA counts and number of detected RNA features per cell |
| `nCount_SCT`, `nFeature_SCT` | Count and feature summaries for the SCT assay, where available |
| `cellType5` | Broad cell-group annotation used in the project analysis |
| `cellType`, `cellType_subType`, `inn`, `area` | Additional project annotations, where available |
| `louvain`, `seurat_clusters`, `RNA_snn_res.*`, `SCT_snn_res.*` | Clustering assignments from the source analysis |

In R, run the following from this folder to read one export.

```r
metadata <- read.delim(
  "STAGATE/cellNewBin_19w_merge_clean__metadata.tsv",
  check.names = FALSE,
  na.strings = "NA"
)
```

All 15 newly exported TSVs were checked against their server copies using SHA-256.
