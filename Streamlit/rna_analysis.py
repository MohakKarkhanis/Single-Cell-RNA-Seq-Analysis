# -*- coding: utf-8 -*-
"""RNA_Analysis.ipynb"""

import streamlit as st
import scanpy as sc
import os
import matplotlib.pyplot as plt
import pandas as pd

model_dir = "Streamlit"
DATA_PATH = "C:/Users/mohak/Desktop/adata.h5ad"

st.set_page_config(layout="wide", page_title="scRNA-Seq Explorer")
st.title("Single-Cell RNA-Seq Data Explorer")

@st.cache_resource
def load_data(path):
    try:
        return sc.read(path)
    except FileNotFoundError:
        st.error(f"Error: AnnData file not found at {path}.")
        return None

adata = load_data(DATA_PATH)

if adata is not None:
    st.sidebar.header("Explore Data")
    color_options = [col for col in adata.obs.columns if adata.obs[col].dtype == 'category' or adata.obs[col].dtype == 'object']
    if 'leiden' in color_options:
        color_options.remove('leiden')
        color_options = ['leiden'] + color_options

    selected_color_by = st.sidebar.selectbox("Color UMAP by:", options=color_options, index=0)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=selected_color_by, legend_loc='on data',
               title=f"UMAP colored by {selected_color_by}", show=False, ax=ax)
    st.pyplot(fig)
    plt.close(fig) # Close specific figure to free memory

    if 'rank_genes_groups' in adata.uns:
        st.header("Marker Gene Exploration")
        available_clusters = adata.obs['leiden'].cat.categories.tolist()
        selected_cluster_for_markers = st.sidebar.selectbox("Select Cluster to view Markers:", options=available_clusters)
        n_top_genes_display = st.sidebar.slider("Number of top marker genes to display:", min_value=5, max_value=50, value=10)

        st.subheader(f"Top {n_top_genes_display} Marker Genes for {selected_cluster_for_markers}")
        try:
            marker_names = pd.DataFrame(adata.uns['rank_genes_groups']['names'])[selected_cluster_for_markers].head(n_top_genes_display).tolist()
            marker_logfc = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])[selected_cluster_for_markers].head(n_top_genes_display).tolist()
            marker_pvals = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])[selected_cluster_for_markers].head(n_top_genes_display).tolist()

            marker_df = pd.DataFrame({
                'Gene': marker_names,
                'Log2FC': marker_logfc,
                'Adjusted P-value': marker_pvals
            })
            st.dataframe(marker_df)
        except KeyError:
            st.warning(f"Marker genes for '{selected_cluster_for_markers}' not found. Check if rank_genes_groups was run correctly or if cluster names were changed post-ranking.")
        except Exception as e:
            st.error(f"An error occurred: {e}")

else:
    st.warning("Please ensure the AnnData file is saved correctly at the specified DATA_PATH.")

plt.close('all') # Ensure all matplotlib figures are closed

