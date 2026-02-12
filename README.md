# MashClust

An automated workflow for large-scale genomic proximity analysis.

## Workflow Overview

### 1. Genome Downloader
Able to handle massive datasets (>120.000 genomes) without connection failures:
* **Batch Processing:** Downloads genomes in groups of 500 to prevent NCBI `Gateway Timeout` errors.
* **Retry Logic:** Built-in automatic retries with strategic pauses to bypass API rate limits and connection resets.
* **Auto-Organization:** Extracts and renames files using official metadata: `genus_species_taxid_accession`.

### 2. Mash Distance Analysis
* **Fast Sketching:** Efficient MinHash-based genomic fingerprinting.
* **Distance Calculation:** Computes Jaccard distances and mutation estimates across all pairs.
* **Representative Selection:** Identifies centroid genomes to reduce redundancy in large-scale datasets.

### 3. Advanced Visualization
Interactive HTML outputs powered by **Plotly**:
* **Dual-Clustered Heatmap:** Features synchronized top and left dendrograms via UPGMA clustering.
* **Species Boundary Focus:** Custom color scale highlighting the **0.06 (Blue) to 0.07 (Yellow)** range to visually distinguish the 94-95% ANI species boundary.
* **PCoA / MDS:** Interactive 2D projection of the genomic space for population structure analysis.


## 🛠️ Requirements

* **Tools:** `datasets` (NCBI CLI), `mash`, `quicktree`.
* **Python Libraries:** `numpy`, `scipy`, `plotly`, `scikit-learn`, `click`, `jsonlines`.

## 📊 Key Outputs

| File | Description |
| :--- | :--- |
| `representatives_clustered_heatmap.html` | Interactive dual-dendrogram heatmap with custom species-boundary scale. |
| `representatives_pcoa.html` | Multi-dimensional scaling plot for genomic space visualization. |
| `genome_list.txt` | Formatted list of paths for all successfully processed genomes. |
| `representatives_tree.nwk` | Rapid phylogenetic tree in Newick format for downstream analysis. |
