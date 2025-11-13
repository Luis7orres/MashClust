#!/usr/bin/env python3
"""
Script 3: Create interactive visualizations from clustering data
Generates HTML (Plotly) or SVG (matplotlib) visualizations
"""


import json
import numpy as np
from pathlib import Path
import argparse
import sys
import subprocess
import re
import shutil
import os



def load_clustering_data(data_file):
    """Load clustering data from JSON file"""
    print(f"[INFO] Loading clustering data from: {data_file}")
    
    try:
        with open(data_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"[ERROR] File not found: {data_file}")
        print(f"[ERROR] Run mash_distances.py first!")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"[ERROR] Invalid JSON file: {e}")
        sys.exit(1)
    
    # Convert back to appropriate formats
    clusters = [set(c) for c in data['clusters']]
    neighbors = {k: set(v) for k, v in data['neighbors'].items()}
    genome_names = {k: Path(v) for k, v in data['genome_names'].items()}
    identity_threshold = data['identity_threshold']
    
    print(f"[INFO] Loaded {len(clusters)} clusters")
    print(f"[INFO] Total genomes: {len(genome_names)}")
    print(f"[INFO] Identity threshold: {identity_threshold*100:.2f}%")
    print(f"[INFO] Data from: {data.get('timestamp', 'unknown')}")
    
    return clusters, neighbors, genome_names, identity_threshold



def find_binary(name, provided_path=None):
    """Find binary in PATH or conda environment"""
    if provided_path:
        bin_path = Path(provided_path)
        if bin_path.exists():
            return str(bin_path)
    
    # Try shutil.which
    found = shutil.which(name)
    if found:
        return found
    
    # Try conda environment
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        conda_bin = Path(conda_prefix) / 'bin' / name
        if conda_bin.exists():
            return str(conda_bin)
    
    # Try common conda locations
    home = Path.home()
    common_paths = [
        home / 'anaconda3' / 'bin' / name,
        home / 'miniconda3' / 'bin' / name,
        home / '.conda' / 'envs' / 'base' / 'bin' / name,
        Path('/opt/anaconda3/bin') / name,
        Path('/opt/miniconda3/bin') / name,
    ]
    
    for path in common_paths:
        if path.exists():
            return str(path)
    
    return None



def extract_genome_name(path_str):
    """
    Extract full genome name from path
    Example: /path/to/mycobacterium_tuberculosis_1773_GCF_001323605_1/GENOME.fna
    Returns: mycobacterium_tuberculosis_1773_GCF_001323605_1
    """
    path = Path(path_str)
    
    # Get parent directory name (where the genome file is located)
    parent_dir = path.parent.name
    
    # If it looks like a genome directory, return it
    if 'GCF_' in parent_dir or 'GCA_' in parent_dir:
        return parent_dir
    
    # Otherwise, try to extract from filename
    filename = path.name
    
    # Try to match pattern like: mycobacterium_tuberculosis_1773_GCF_001323605_1
    match = re.search(r'([a-z_]+_GC[FA]_\d{9}_\d+)', filename, re.IGNORECASE)
    if match:
        return match.group(1)
    
    # Fallback: return filename without extension
    return filename.split('.')[0] if '.' in filename else filename



def generate_quicktree_from_distances(output_dir, quicktree_path, dist_file):
    """
    Generate phylogenetic tree using QuickTree from Mash distance matrix
    Uses RELAXED PHYLIP format to preserve full genome names
    """
    
    print(f"[3/5] Generating phylogenetic tree with QuickTree...")
    
    output_dir = Path(output_dir)
    representatives_file = output_dir / "representatives.txt"
    dist_file = Path(dist_file)
    
    if not representatives_file.exists():
        print(f"      [ERROR] Representatives file not found: {representatives_file}")
        return None
    
    if not dist_file.exists():
        print(f"      [ERROR] Distance file not found: {dist_file}")
        return None
    
    # Check binary
    quicktree_bin = find_binary('quicktree', quicktree_path)
    
    if not quicktree_bin:
        print(f"      [ERROR] QuickTree not found!")
        print(f"      [INFO] Install: conda install -c bioconda quicktree")
        print(f"      [INFO] Or specify path: --quicktree /path/to/quicktree")
        return None
    
    print(f"      Found QuickTree: {quicktree_bin}")
    
    # Read representatives
    with open(representatives_file, 'r') as f:
        representatives = [line.strip() for line in f if line.strip()]
    
    print(f"      Found {len(representatives)} representative genomes")
    
    # Convert representatives paths to full genome names
    rep_names = set(extract_genome_name(rep) for rep in representatives)
    
    # Parse distance matrix and extract representatives subset
    print(f"      Extracting representative distances from matrix...")
    
    genome_order = []
    distance_matrix = []
    
    with open(dist_file, 'r') as f:
        # Read header
        header_line = f.readline().strip()
        header_parts = header_line.split('\t')
        all_genomes = header_parts[1:]
        
        # Find representative indices
        rep_indices = []
        for i, genome_path in enumerate(all_genomes):
            genome_name = extract_genome_name(genome_path)
            if genome_name in rep_names:
                rep_indices.append(i)
                genome_order.append(genome_name)
        
        # Read distance rows
        for row_idx, line in enumerate(f):
            if not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            query_path = parts[0]
            query_name = extract_genome_name(query_path)
            
            if query_name in rep_names:
                distances = parts[1:]
                # Extract only representative columns
                rep_distances = [float(distances[i]) for i in rep_indices if i < len(distances)]
                distance_matrix.append(rep_distances)
    
    if len(genome_order) != len(distance_matrix):
        print(f"      [ERROR] Matrix dimension mismatch")
        return None
    
    print(f"      Matrix size: {len(genome_order)} x {len(genome_order)}")
    
    # Write RELAXED PHYLIP format (no name length restriction)
    phylip_file = output_dir / "representatives_phylip.dist"
    n = len(genome_order)
    
    print(f"      Writing RELAXED PHYLIP format with FULL genome names...")
    with open(phylip_file, 'w') as f:
        f.write(f"{n}\n")
        for i, genome_name in enumerate(genome_order):
            # RELAXED PHYLIP: no truncation, names separated by spaces/tabs
            dist_line = '\t'.join(f"{d:.6f}" for d in distance_matrix[i])
            f.write(f"{genome_name}\t{dist_line}\n")
    
    print(f"      Phylip format: {phylip_file.name}")
    print(f"      Names preserved: {genome_order[0]} (example)")
    
    # Run QuickTree
    tree_file = output_dir / "representatives_tree.nwk"
    
    try:
        # QuickTree command:
        # -in m : input is a distance matrix
        # -out t : output a tree
        
        cmd = [str(quicktree_bin), "-in", "m", "-out", "t", str(phylip_file)]
        
        print(f"      Command: {' '.join(cmd)}")
        
        with open(tree_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, check=True, 
                                  stderr=subprocess.PIPE, text=True)
        
        print(f"      Saved: {tree_file.name}")
        print(f"      ✓ Tree contains FULL genome names!")
        
        # Verify tree has full names
        with open(tree_file, 'r') as f:
            tree_content = f.read()
            if genome_order[0] in tree_content:
                print(f"      ✓ Verified: Full names present in tree")
            else:
                print(f"      ⚠ Warning: Names may have been truncated")
        
        return tree_file
        
    except subprocess.CalledProcessError as e:
        print(f"      [ERROR] QuickTree failed: {e.stderr}")
        return None



def create_plotly_visualizations(clusters, neighbors, genome_names, output_dir, 
                                 identity_threshold, quicktree_path, dist_file):
    """Create interactive HTML visualizations using Plotly"""
    
    try:
        import plotly.graph_objects as go
        import plotly.express as px
        from plotly.subplots import make_subplots
    except ImportError:
        print(f"[ERROR] Plotly not installed!")
        print(f"[INFO] Install with: pip install plotly kaleido")
        sys.exit(1)
    
    print(f"\n{'='*70}")
    print(f"CREATING INTERACTIVE VISUALIZATIONS")
    print(f"{'='*70}\n")
    
    output_dir = Path(output_dir)
    cluster_sizes = [len(c) for c in clusters]
    
    # ========== 1. CLUSTER SIZE DISTRIBUTION ==========
    print(f"[1/5] Creating cluster size distribution...")
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Histogram (Linear Scale)', 
            'Histogram (Log Scale)', 
            'Box Plot Statistics', 
            'Cumulative Distribution'
        ),
        specs=[[{"type": "histogram"}, {"type": "histogram"}],
               [{"type": "box"}, {"type": "scatter"}]],
        vertical_spacing=0.12,
        horizontal_spacing=0.1
    )
    
    # Linear histogram
    fig.add_trace(
        go.Histogram(
            x=cluster_sizes, 
            nbinsx=min(50, len(set(cluster_sizes))),
            name='Linear',
            marker_color='steelblue',
            hovertemplate='Size: %{x}<br>Count: %{y}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # Log histogram
    fig.add_trace(
        go.Histogram(
            x=cluster_sizes, 
            nbinsx=min(50, len(set(cluster_sizes))),
            name='Log',
            marker_color='darkgreen',
            hovertemplate='Size: %{x}<br>Count: %{y}<extra></extra>'
        ),
        row=1, col=2
    )
    
    # Box plot
    fig.add_trace(
        go.Box(
            y=cluster_sizes, 
            name='Statistics',
            marker_color='lightcoral',
            boxmean='sd',
            hovertemplate='Size: %{y}<extra></extra>'
        ),
        row=2, col=1
    )
    
    # Cumulative
    sorted_sizes = sorted(cluster_sizes, reverse=True)
    cumsum = np.cumsum(sorted_sizes)
    fig.add_trace(
        go.Scatter(
            x=list(range(1, len(cumsum)+1)), 
            y=cumsum,
            mode='lines',
            name='Cumulative',
            line=dict(color='purple', width=2),
            fill='tozeroy',
            hovertemplate='Clusters: %{x}<br>Genomes: %{y}<extra></extra>'
        ),
        row=2, col=2
    )
    
    # Update axes
    fig.update_xaxes(title_text="Cluster Size", row=1, col=1)
    fig.update_yaxes(title_text="Frequency", row=1, col=1)
    fig.update_xaxes(title_text="Cluster Size", row=1, col=2)
    fig.update_yaxes(title_text="Frequency (log)", type="log", row=1, col=2)
    fig.update_yaxes(title_text="Cluster Size", row=2, col=1)
    fig.update_xaxes(title_text="Number of Clusters", row=2, col=2)
    fig.update_yaxes(title_text="Cumulative Genomes", row=2, col=2)
    
    fig.update_layout(
        title_text=f"Cluster Size Distribution<br><sub>{len(clusters)} clusters, {len(genome_names)} genomes, threshold={identity_threshold*100:.2f}%</sub>",
        title_x=0.5,
        showlegend=False,
        height=800,
        template='plotly_white'
    )
    
    output_file = output_dir / 'cluster_distribution.html'
    fig.write_html(output_file)
    print(f"      Saved: {output_file.name}")
    
    # ========== 2. CLUSTER HIERARCHY ==========
    print(f"[2/5] Creating cluster hierarchy...")
    
    clusters_sorted = sorted(clusters, key=len, reverse=True)
    cluster_sizes_sorted = [len(c) for c in clusters_sorted]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        y=[f'C{i+1}' for i in range(len(clusters_sorted))],
        x=cluster_sizes_sorted,
        orientation='h',
        marker=dict(
            color=cluster_sizes_sorted,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Size", thickness=15)
        ),
        text=cluster_sizes_sorted,
        textposition='outside',
        hovertemplate='<b>Cluster %{y}</b><br>Size: %{x} genomes<extra></extra>'
    ))
    
    if len(clusters_sorted) > 100:
        fig.update_yaxes(showticklabels=False, title_text="Clusters (sorted by size)")
    else:
        fig.update_yaxes(title_text="Cluster ID")
    
    fig.update_xaxes(title_text="Number of Genomes")
    
    fig.update_layout(
        title=f"Hierarchical View of {len(clusters)} Clusters<br><sub>Sorted by size</sub>",
        title_x=0.5,
        height=max(600, len(clusters_sorted) * 6),
        template='plotly_white',
        showlegend=False
    )
    
    output_file = output_dir / 'cluster_hierarchy.html'
    fig.write_html(output_file)
    print(f"      Saved: {output_file.name}")
    
    # ========== 3. PHYLOGENETIC TREE (QUICKTREE) ==========
    generate_quicktree_from_distances(output_dir, quicktree_path, dist_file)
    
    # ========== 4. NETWORK GRAPH ==========
    print(f"[4/5] Creating network graph...")
    
    try:
        import networkx as nx
        
        # Decide on subset
        if len(genome_names) > 1000:
            print(f"      Large dataset ({len(genome_names)} genomes) - plotting top 10 clusters")
            clusters_to_plot = clusters_sorted[:10]
            title_suffix = "Top 10 Clusters"
        elif len(genome_names) > 500:
            print(f"      Medium dataset - plotting top 20 clusters")
            clusters_to_plot = clusters_sorted[:20]
            title_suffix = "Top 20 Clusters"
        else:
            print(f"      Small dataset - plotting all clusters")
            clusters_to_plot = clusters
            title_suffix = f"All {len(clusters)} Clusters"
        
        # Build graph
        G = nx.Graph()
        node_cluster = {}
        
        for cluster_idx, cluster in enumerate(clusters_to_plot):
            for node in cluster:
                G.add_node(node)
                node_cluster[node] = cluster_idx
                
                # Add edges
                for neighbor in neighbors.get(node, set()):
                    if neighbor in cluster:
                        G.add_edge(node, neighbor)
        
        print(f"      Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        # Layout
        if G.number_of_nodes() > 500:
            pos = nx.kamada_kawai_layout(G)
        else:
            pos = nx.spring_layout(G, k=2/np.sqrt(G.number_of_nodes()), iterations=50, seed=42)
        
        # Edges
        edge_x = []
        edge_y = []
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines',
            showlegend=False
        )
        
        # Nodes
        node_x = []
        node_y = []
        node_color = []
        node_text = []
        node_size = []
        
        color_map = px.colors.qualitative.Set3
        
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            
            cluster_id = node_cluster[node]
            node_color.append(color_map[cluster_id % len(color_map)])
            
            degree = G.degree(node)
            node_size.append(5 + degree * 2)
            
            node_text.append(
                f"<b>{str(node)[:50]}</b><br>"
                f"Cluster: {cluster_id+1}<br>"
                f"Connections: {degree}"
            )
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            text=node_text,
            marker=dict(
                color=node_color,
                size=node_size,
                line=dict(width=0.5, color='white')
            ),
            showlegend=False
        )
        
        fig = go.Figure(
            data=[edge_trace, node_trace],
            layout=go.Layout(
                title=f"Network Graph - {title_suffix}<br><sub>{G.number_of_nodes()} genomes, {G.number_of_edges()} connections</sub>",
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=80),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                height=800,
                template='plotly_white'
            )
        )
        
        output_file = output_dir / 'network_graph.html'
        fig.write_html(output_file)
        print(f"      Saved: {output_file.name}")
        
    except ImportError:
        print(f"      Skipped (NetworkX not installed)")
        print(f"      Install with: pip install networkx")
    
    # ========== 5. SUMMARY TABLE ==========
    print(f"[5/5] Creating summary statistics...")
    
    stats = {
        'Total Genomes': [len(genome_names)],
        'Total Clusters': [len(clusters)],
        'Singleton Clusters': [sum(1 for c in clusters if len(c) == 1)],
        'Mean Cluster Size': [f"{np.mean(cluster_sizes):.2f}"],
        'Median Cluster Size': [f"{np.median(cluster_sizes):.0f}"],
        'Largest Cluster': [max(cluster_sizes)],
        'Smallest Cluster': [min(cluster_sizes)],
        'Identity Threshold': [f"{identity_threshold*100:.2f}%"],
        'Genomes with Neighbors': [len([g for g in genome_names if g in neighbors])],
        'Isolated Genomes': [len([g for g in genome_names if g not in neighbors])]
    }
    
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=['<b>Metric</b>', '<b>Value</b>'],
            fill_color='steelblue',
            font=dict(color='white', size=14),
            align='left',
            height=35
        ),
        cells=dict(
            values=[list(stats.keys()), [str(v[0]) for v in stats.values()]],
            fill_color='lavender',
            font=dict(size=12),
            align='left',
            height=30
        )
    )])
    
    fig.update_layout(
        title="Clustering Summary Statistics",
        title_x=0.5,
        height=450,
        margin=dict(l=20, r=20, t=60, b=20)
    )
    
    output_file = output_dir / 'summary_stats.html'
    fig.write_html(output_file)
    print(f"      Saved: {output_file.name}")
    
    print(f"\n[SUCCESS] All visualizations created!")



def create_matplotlib_svg(clusters, neighbors, genome_names, output_dir, 
                         identity_threshold, quicktree_path, dist_file):
    """Create SVG visualizations using matplotlib"""
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print(f"[ERROR] Matplotlib not installed!")
        print(f"[INFO] Install with: pip install matplotlib seaborn")
        sys.exit(1)
    
    print(f"\n{'='*70}")
    print(f"CREATING SVG VISUALIZATIONS")
    print(f"{'='*70}\n")
    
    output_dir = Path(output_dir)
    cluster_sizes = [len(c) for c in clusters]
    
    # 1. Distribution
    print(f"[1/3] Creating distribution plot...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    axes[0, 0].hist(cluster_sizes, bins=min(50, len(set(cluster_sizes))), 
                    edgecolor='black', alpha=0.7, color='steelblue')
    axes[0, 0].set_xlabel('Cluster Size', fontsize=12)
    axes[0, 0].set_ylabel('Frequency', fontsize=12)
    axes[0, 0].set_title('Histogram (Linear)', fontsize=14, fontweight='bold')
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    axes[0, 1].hist(cluster_sizes, bins=min(50, len(set(cluster_sizes))), 
                    edgecolor='black', alpha=0.7, color='darkgreen')
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_xlabel('Cluster Size', fontsize=12)
    axes[0, 1].set_ylabel('Frequency (log)', fontsize=12)
    axes[0, 1].set_title('Histogram (Log)', fontsize=14, fontweight='bold')
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    axes[1, 0].boxplot(cluster_sizes, vert=True, patch_artist=True,
                       boxprops=dict(facecolor='lightcoral', alpha=0.7))
    axes[1, 0].set_ylabel('Cluster Size', fontsize=12)
    axes[1, 0].set_title('Statistics', fontsize=14, fontweight='bold')
    axes[1, 0].grid(axis='y', alpha=0.3)
    
    sorted_sizes = sorted(cluster_sizes, reverse=True)
    cumsum = np.cumsum(sorted_sizes)
    axes[1, 1].plot(range(1, len(cumsum)+1), cumsum, linewidth=2, color='purple')
    axes[1, 1].fill_between(range(1, len(cumsum)+1), cumsum, alpha=0.3)
    axes[1, 1].set_xlabel('Number of Clusters', fontsize=12)
    axes[1, 1].set_ylabel('Cumulative Genomes', fontsize=12)
    axes[1, 1].set_title('Cumulative Distribution', fontsize=14, fontweight='bold')
    axes[1, 1].grid(alpha=0.3)
    
    plt.suptitle(f'Cluster Distribution - {len(clusters)} clusters, {len(genome_names)} genomes', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    output_file = output_dir / 'cluster_distribution.svg'
    plt.savefig(output_file, format='svg', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"      Saved: {output_file.name}")
    
    # 2. Hierarchy
    print(f"[2/3] Creating hierarchy plot...")
    
    clusters_sorted = sorted(clusters, key=len, reverse=True)
    widths = [len(c) for c in clusters_sorted]
    
    fig, ax = plt.subplots(figsize=(16, max(10, len(clusters_sorted) * 0.1)))
    
    colors = [plt.cm.viridis(w / max(widths)) for w in widths]
    bars = ax.barh(range(len(widths)), widths, color=colors, edgecolor='black', linewidth=0.3)
    
    # Add size labels
    for i, (bar, width) in enumerate(zip(bars, widths)):
        if width > np.median(widths) or i < 20:
            ax.text(width, i, f'  {width}', va='center', fontsize=8)
    
    ax.set_ylabel('Cluster ID', fontsize=12)
    ax.set_xlabel('Number of Genomes', fontsize=12)
    ax.set_title(f'Cluster Hierarchy - {len(clusters)} Clusters', fontsize=14, fontweight='bold')
    
    if len(clusters_sorted) <= 50:
        ax.set_yticks(range(len(widths)))
        ax.set_yticklabels([f'C{i+1}' for i in range(len(widths))], fontsize=8)
    
    plt.tight_layout()
    output_file = output_dir / 'cluster_hierarchy.svg'
    plt.savefig(output_file, format='svg', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"      Saved: {output_file.name}")
    
    # 3. Phylogenetic Tree
    print(f"[3/3] Generating phylogenetic tree...")
    generate_quicktree_from_distances(output_dir, quicktree_path, dist_file)
    
    print(f"\n[SUCCESS] All SVG visualizations created!")



def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Create visualizations from clustering results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # HTML interactive (auto-detect QuickTree)
  python mash_visualize.py mash_output
  
  # Specify QuickTree location explicitly
  python mash_visualize.py mash_output --quicktree /path/to/quicktree
  
  # SVG format
  python mash_visualize.py mash_output --format svg

Installation:
  # Install QuickTree (conda)
  conda install -c bioconda quicktree
  
  # Or compile from source
  git clone https://github.com/khowe/quicktree
  cd quicktree
  make
  # Binary will be at: ./quicktree

Pipeline:
  1. Uses Mash distances already calculated
  2. QuickTree builds tree from distance matrix
  3. Tree preserves FULL genome names (e.g., mycobacterium_tuberculosis_1773_GCF_001323605_1)
  4. FAST: No need for alignment step!

Output files (HTML format):
  - cluster_distribution.html     : Size distribution plots
  - cluster_hierarchy.html        : Hierarchical bar chart
  - network_graph.html            : Interactive network
  - summary_stats.html            : Statistics table
  - representatives_tree.nwk      : Phylogenetic tree with FULL genome names
  - representatives_phylip.dist   : Distance matrix in PHYLIP format
        """
    )
    
    parser.add_argument("output_dir", help="Directory with clustering_data.json")
    parser.add_argument("--format", choices=['html', 'svg'], default='html',
                       help="Output format (default: html)")
    parser.add_argument("--quicktree", default=None,
                       help="Path to QuickTree binary (default: search in PATH/conda)")
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    
    if not output_dir.exists():
        print(f"[ERROR] Directory not found: {output_dir}")
        sys.exit(1)
    
    data_file = output_dir / 'clustering_data.json'
    
    if not data_file.exists():
        print(f"[ERROR] Clustering data not found: {data_file}")
        print(f"[ERROR] Run mash_distances.py first!")
        sys.exit(1)
    
    dist_file = output_dir / 'distances.txt'
    
    if not dist_file.exists():
        print(f"[ERROR] Distance file not found: {dist_file}")
        print(f"[ERROR] Run mash_distances.py first!")
        sys.exit(1)
    
    # Load data
    clusters, neighbors, genome_names, identity_threshold = load_clustering_data(data_file)
    
    # Create visualizations
    if args.format == 'html':
        create_plotly_visualizations(clusters, neighbors, genome_names, output_dir, 
                                     identity_threshold, args.quicktree, dist_file)
    else:
        create_matplotlib_svg(clusters, neighbors, genome_names, output_dir, 
                             identity_threshold, args.quicktree, dist_file)
    
    # Summary
    print(f"\n{'='*70}")
    print("VISUALIZATION COMPLETE")
    print(f"{'='*70}")
    print(f"Output directory: {output_dir}")
    print(f"Format: {args.format.upper()}")
    print(f"\nGenerated files:")
    for f in sorted(output_dir.glob(f'*.{args.format}')):
        print(f"  - {f.name}")
    
    tree_file = output_dir / 'representatives_tree.nwk'
    if tree_file.exists():
        print(f"  - {tree_file.name} (distance-based, FULL names)")
    
    phylip_file = output_dir / 'representatives_phylip.dist'
    if phylip_file.exists():
        print(f"  - {phylip_file.name} (for inspection)")
    
    print(f"{'='*70}\n")



if __name__ == "__main__":
    main()
