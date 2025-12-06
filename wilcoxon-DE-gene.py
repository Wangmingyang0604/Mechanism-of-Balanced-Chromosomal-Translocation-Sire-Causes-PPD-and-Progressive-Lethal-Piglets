#singularity exec -B /home/tmpdir/zhiyan/Polytoe/Duanjunhui/:/home/tmpdir/zhiyan/Polytoe/Duanjunhui/ /home/tmpdir/zhiyan/Polytoe/Duanjunhui/software/singcell_python.sif /opt/conda/envs/scRNA/bin/python
import scanpy as sc
import scanpy.external as sce
import doubletdetection
#basic modules
import os
import random
from tqdm import tqdm
import glob
#stats/data modules
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation as mad
from scipy.stats import iqr
import matplotlib.pyplot as plt
import seaborn as sns

###画umap图
##需要使用的值是normalized的数值，而且是经历过log1p的。
##Expects logarithmized data.
##https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html?utm_source


os.chdir('/home/tmpdir/zhiyan/Polytoe/Duanjunhui/single-cell/prepare/cellbender_output/20250718_all_restart/wilcoxon-DE-gene')
os.chdir('/home/tmpdir/zhiyan/Polytoe/Duanjunhui/single-cell/prepare/cellbender_output/20250718_all_restart/Enterocytes/wilxocon-DE')
os.chdir('/home/tmpdir/zhiyan/Polytoe/Duanjunhui/single-cell/prepare/cellbender_output/20250718_all_restart/immune/immune_wilcoxon_DE')
os.getcwd()
adata = sc.read("merged_annotated_finally.h5ad")
adata=sc.read("Enterocytes_annotated_finally.h5ad")


import numpy as np
import scipy

X_norm = adata.layers["normlized"]
vals = X_norm.data[:10000] if scipy.sparse.issparse(X_norm) else X_norm.flatten()

print("normlized layer 最大值：", vals.max())
print("normlized layer 最小非零值：", vals[vals > 0].min())
###检查一下数据集是不是normalized的数值，而且是经历过log1p的。


adata
adata.obs.columns
sc.tl.rank_genes_groups(adata, groupby='group', method='wilcoxon', layer='normlized')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('all_celltype_markerGeneInfo_all.csv')


###########这个是赖晨靓师姐的#################################################################################################################
sc.tl.rank_genes_groups(dat,'leiden0.3',method='wilcoxon',key_added='rank_genes_groups_leiden0.3',pts=True)
result = dat.uns['rank_genes_groups_leiden0.3']
groups = result['names'].dtype.names
df1 = pd.DataFrame({group+'_'+ key:result[key][group] for group in groups for key in ['names','logfoldchanges','pvals',]})
df2 = pd.DataFrame({group+'_'+ key:result[key][group] for group in groups for key in ['pts','pts_rest']})
df=pd.concat(
    [
        df1[[f"{group}_names", f"{group}_logfoldchanges", f"{group}_pvals"]]
        .merge(df2[[f"{group}_pts",f"{group}_pts_rest"]], how="left", left_on=f"{group}_names", right_index=True)
        for group in groups
    ],
    axis=1
)

result = []
for i in range(0, df.shape[1], 5):
    n_col = df.columns[i]
    l_col = df.columns[i + 1]
    p_col = df.columns[i + 2]
    pts_col = df.columns[i + 3]
    pts_rest_col = df.columns[i + 4]
    filtered = df[(df[l_col].abs() > 0.01) & (df[p_col] < 0.05) & (df[pts_col] > 0.01)]
    filtered = filtered.dropna(subset=[n_col, l_col, p_col,pts_col,pts_rest_col])
    filtered = filtered[[n_col, l_col, p_col,pts_col,pts_rest_col]]
    result.append(filtered)

DEG = pd.concat(result, axis=1)
DEG.to_csv('chorion_ct_marker(logfc0.01pval0.05pts0.01).csv',index=False)
##########################################################################################################################################







print(adata.obs['cell_type'].unique().tolist())

adata_Enterocytes = adata[adata.obs['cell_type'] == 'ENterocytes'].copy()
print(adata_Enterocytes.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Enterocytes,'group',method='wilcoxon',layer='normlized')
result = adata_Enterocytes.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Enterocytes_celltype_markerGeneInfo_all.csv')




adata_Immune = adata[adata.obs['cell_type'] == 'Immune'].copy()
print(adata_Immune.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Immune,'group',method='wilcoxon',layer='normlized')
result = adata_Immune.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Immune_celltype_markerGeneInfo_all.csv')


adata_Mesothelial = adata[adata.obs['cell_type'] == 'Mesothelial'].copy()
print(adata_Mesothelial.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Mesothelial,'group',method='wilcoxon',layer='normlized')
result = adata_Mesothelial.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Mesothelial_celltype_markerGeneInfo_all.csv')




adata_Neurons = adata[adata.obs['cell_type'] == 'Neurons'].copy()
print(adata_Neurons.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Neurons,'group',method='wilcoxon',layer='normlized')
result = adata_Neurons.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Neurons_celltype_markerGeneInfo_all.csv')



adata_Tuft = adata[adata.obs['cell_type'] == 'Tuft'].copy()
print(adata_Tuft.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Tuft,'group',method='wilcoxon',layer='normlized')
result = adata_Tuft.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Tuft_celltype_markerGeneInfo_all.csv')




adata_Venous_EC = adata[adata.obs['cell_type'] == 'Venous EC'].copy()
print(adata_Venous_EC.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Venous_EC,'group',method='wilcoxon',layer='normlized')
result = adata_Venous_EC.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Venous_EC_celltype_markerGeneInfo_all.csv')


adata_SMC = adata[adata.obs['cell_type'] == 'SMC'].copy()
print(adata_SMC.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_SMC,'group',method='wilcoxon',layer='normlized')
result = adata_SMC.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('SMC_celltype_markerGeneInfo_all.csv')



print(adata.obs['cell_type'].unique().tolist())

adata_ENterocytes_BEST4 = adata[adata.obs['cell_type'] == 'ENterocytes BEST4+ '].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_ENterocytes_BEST4.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_ENterocytes_BEST4,'group',method='wilcoxon',layer='normlized')
result = adata_ENterocytes_BEST4.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('ENterocytes_BEST4_celltype_markerGeneInfo_all.csv')





print(adata.obs['cell_type'].unique().tolist())

adata_EEC = adata[adata.obs['cell_type'] == 'EEC'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_EEC.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_EEC,'group',method='wilcoxon',layer='normlized')
result = adata_EEC.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('EEC_celltype_markerGeneInfo_all.csv')


print(adata.obs['cell_type'].unique().tolist())

adata_Tuft = adata[adata.obs['cell_type'] == 'Tuft'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_Tuft.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Tuft,'group',method='wilcoxon',layer='normlized')
result = adata_Tuft.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Tuft_celltype_markerGeneInfo_all.csv')




print(adata.obs['cell_type'].unique().tolist())

adata_TA = adata[adata.obs['cell_type'] == 'TA'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_TA.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_TA,'group',method='wilcoxon',layer='normlized')
result = adata_TA.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('TA_celltype_markerGeneInfo_all.csv')




print(adata.obs['cell_type'].unique().tolist())

adata_LEC = adata[adata.obs['cell_type'] == 'LEC'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_LEC.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_LEC,'group',method='wilcoxon',layer='normlized')
result = adata_LEC.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('LEC_celltype_markerGeneInfo_all.csv')


print(adata.obs['cell_type'].unique().tolist())

adata_Fibroblasts = adata[adata.obs['cell_type'] == 'Fibroblasts'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_Fibroblasts.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Fibroblasts,'group',method='wilcoxon',layer='normlized')
result = adata_Fibroblasts.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Fibroblasts_celltype_markerGeneInfo_all.csv')



print(adata.obs['cell_type'].unique().tolist())

adata_Goblet = adata[adata.obs['cell_type'] == 'Goblet'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_Goblet.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Goblet,'group',method='wilcoxon',layer='normlized')
result = adata_Goblet.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Goblet_celltype_markerGeneInfo_all.csv')



print(adata.obs['cell_type'].unique().tolist())

adata_Stem = adata[adata.obs['cell_type'] == 'Stem'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_Stem.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Stem,'group',method='wilcoxon',layer='normlized')
result = adata_Stem.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Stem_celltype_markerGeneInfo_all.csv')




print(adata.obs['cell_type'].unique().tolist())

adata_Glial = adata[adata.obs['cell_type'] == 'Glial'].copy() ###这个ENterocytes BEST4+ 需要保留一个空格。
print(adata_Glial.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Glial,'group',method='wilcoxon',layer='normlized')
result = adata_Glial.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Glial_celltype_markerGeneInfo_all.csv')





os.chdir('/home/tmpdir/zhiyan/Polytoe/Duanjunhui/single-cell/prepare/cellbender_output/20250718_all_restart/Enterocytes/wilxocon-DE')
os.getcwd()
adata=sc.read("Enterocytes_annotated_finally.h5ad")
adata
adata.obs.columns

print(adata.obs['cell_type'].unique().tolist())

ILCs_Immune = adata[adata.obs['cell_type'] == 'ILCs'].copy()
print(ILCs_Immune.obs['group'].value_counts())
sc.tl.rank_genes_groups(ILCs_Immune,'group',method='wilcoxon',layer='normlized')
result = ILCs_Immune.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('ILCs_Immune_celltype_markerGeneInfo_all.csv')


print(adata.obs['cell_type'].unique().tolist())
IFIT1_ENterocytes = adata[adata.obs['cell_type'] == 'IFIT1+ ENterocytes'].copy()
print(IFIT1_ENterocytes.obs['group'].value_counts())
sc.tl.rank_genes_groups(IFIT1_ENterocytes,'group',method='wilcoxon',layer='normlized')
result = IFIT1_ENterocytes.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('IFIT1_ENterocytes_celltype_markerGeneInfo_all.csv')




print(adata.obs['cell_type'].unique().tolist())
ACE_ENterocytes = adata[adata.obs['cell_type'] == 'ACE+ ENterocytes'].copy()
print(ACE_ENterocytes.obs['group'].value_counts())
sc.tl.rank_genes_groups(ACE_ENterocytes,'group',method='wilcoxon',layer='normlized')
result = ACE_ENterocytes.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('ACE_ENterocytes_celltype_markerGeneInfo_all.csv')


print(adata.obs['cell_type'].unique().tolist())
APOC3_ENterocytes = adata[adata.obs['cell_type'] == 'APOC3+ ENterocytes'].copy()
print(APOC3_ENterocytes.obs['group'].value_counts())
sc.tl.rank_genes_groups(APOC3_ENterocytes,'group',method='wilcoxon',layer='normlized')
result = APOC3_ENterocytes.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('APOC3_ENterocytes_celltype_markerGeneInfo_all.csv')





print(adata.obs['cell_type'].unique().tolist())
NEIL1_ENterocytes = adata[adata.obs['cell_type'] == 'NEIL1+ ENterocytes'].copy()
print(NEIL1_ENterocytes.obs['group'].value_counts())
sc.tl.rank_genes_groups(NEIL1_ENterocytes,'group',method='wilcoxon',layer='normlized')
result = NEIL1_ENterocytes.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('NEIL1_ENterocytes_celltype_markerGeneInfo_all.csv')


print(adata.obs['cell_type'].unique().tolist())

adata_Enterocytes = adata[adata.obs['cell_type'] == 'ENterocytes'].copy()
print(adata_Enterocytes.obs['group'].value_counts())
sc.tl.rank_genes_groups(adata_Enterocytes,'group',method='wilcoxon',layer='normlized')
result = adata_Enterocytes.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores','logfoldchanges','pvals']})
markerGeneRank_all.to_csv('Enterocytes_celltype_markerGeneInfo_all.csv')




import scanpy as sc

print(adata.obs['cell_type'].unique().tolist())

AQP8_ENterocytes = adata[adata.obs['cell_type'] == 'AQP8+ ENterocytes'].copy()
# 设置保存路径为当前目录
sc.settings.figdir = "."

# 目标基因列表
genes = [
    "TBC1D4", "G6PC", "STAT3", "IKBKB", "SLC2A2", "PRKAA2",
    "CREB3L3", "MGEA5", "INSR", "PTPN1", "PPARGC1A", "TRIB3"
]

# 绘图并保存为 linshi.pdf
sc.pl.dotplot(
    AQP8_ENterocytes,
    var_names=genes,
    groupby="group",
    layer="normlized",
    standard_scale="var",
    cmap="Spectral_r",
    save="linshi.pdf"
)




