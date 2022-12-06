#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc


# In[84]:


os.getcwd()
os.chdir('/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/')
data_path = '/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/'

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=300, figsize=(10,8), format='png')
sc.settings.figdir = 'Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/'



# In[9]:





# In[40]:


adata = sc.read_h5ad(data_path + 'iec_controls.h5ad')

#adata = sc.read(data_path + 'iec_controls.h5ad')


# In[85]:


adata.X = adata.X.astype('float64')  # this is not required and results will be comparable without it
adata


# In[86]:


sc.pp.neighbors(adata, n_pcs=18)
sc.tl.draw_graph(adata)
sc.tl.draw_graph(adata)
sc.pl.umap(adata, color='cells', legend_loc='on data')


# In[83]:


sc.pl.draw_graph(adata, color='cells', legend_loc='on data')


# In[69]:


## Change datatype to category
#adata.obs['cells'] = cells[0].astype('category').values
adata.obs['cells'] = adata.obs['cells'].astype('category').values
#sc.tl.paga(adata, groups='clusters')
## Rename Clusters
adata.obs['cells'].cat.categories = ["EPL1", "MEP1", "MEP2","TA","MED","IEP", "MEP3","EPE1", "ISC",
       "EPE2","EPL2","Goblet","Tuft","EEC","MEP4","Paneth"]


# In[70]:


adata.obs['cells']


# In[71]:


sc.tl.paga(adata, groups='cells')


# In[76]:


sc.pl.paga(adata, threshold=0.03, show="TRUE",save='_Control_1')


# In[77]:


sc.pl.paga(adata, threshold=0.03, show="TRUE",layout='eq_tree',save='_Control_2')
sc.pl.paga(adata, threshold=0.03, show="TRUE",layout='eq_tree',root=8,save='_Control_3' )
#sc.tl.umap(adata, init_pos='paga')


# In[87]:


sc.pl.paga(adata, threshold=0.03, show="TRUE",layout='fr')


# In[80]:


sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
adata.uns['iroot'] = np.flatnonzero(adata.obs['cells']  == 'ISC')[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color=['cells', 'dpt_pseudotime'], legend_loc='on data')


# In[82]:


sc.tl.diffmap(adata)
#sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')


# In[ ]:




