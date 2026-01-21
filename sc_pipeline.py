import scanpy as sc
import pandas as pd
import numpy as np
import h5py
import anndata
from scipy import sparse
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
# Fix imports
import matplotlib.patches as mpatches 
from matplotlib.patches import Rectangle, Polygon
import matplotlib.path as mpath
import gseapy as gp

class SCPipeline:
    def __init__(self, config, inspect_data = False):
        self.config = config
        self.adata = None
        self.marker_df = None
        self.data_dir, self.fig_dir = self._setup_dirs()
        
        sc.settings.figdir = self.fig_dir
        sc.settings.set_figure_params(dpi=150, frameon=False)
        
        if inspect_data:
            self._check_h5(self.config['input_path'])
        
        # --- Mode Switching Logic Check ---
        # This line confirms the mode based on the config
        mode = "🧪 Experiment Mode" if self.config.get('use_experiment') else "📦 Standard Mode"
        print(f"\n🚀 [Pipeline] Initialized: {self.config['dataset_name']} | {mode}")

    def _setup_dirs(self):
        base = self.config['output_base']
        name = self.config['dataset_name']
        # Separate output directory for Experiment Mode
        if self.config.get('use_experiment'):
            name += "_Experiment"
        data_dir = os.path.join(base, name, 'Data')
        fig_dir = os.path.join(base, name, 'Figures')
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(fig_dir, exist_ok=True)
        return data_dir, fig_dir

    def _inspect_node(self, name, obj):
        # 缩进显示层级
        indent = "  " * name.count('/')
        
        if isinstance(obj, h5py.Group):
            print(f"{indent}📂 Group: {name}")
            # 检查 Group 的属性 (attrs)
            for key, val in obj.attrs.items():
                print(f"{indent}  🔹 Attr: {key} = {val}")
                
        elif isinstance(obj, h5py.Dataset):
            print(f"{indent}📄 Dataset: {name}")
            print(f"{indent}   Shape: {obj.shape}, Dtype: {obj.dtype}")
            
            # 1. 如果是一维数组，且比较短（可能是类别名），打印出来看看
            if obj.ndim == 1 and obj.shape[0] < 100:
                try:
                    data = np.array(obj)
                    # 处理字节字符串
                    if data.dtype.kind == 'S':
                        data = [x.decode('utf-8') for x in data]
                    print(f"{indent}   👀 内容预览: {data}")
                except:
                    pass
                    
            # 2. 检查 Dataset 的属性 (attrs) - 有时候标签映射藏在这里！
            for key, val in obj.attrs.items():
                print(f"{indent}  🔹 Attr: {key} = {val}")

    def _check_h5(self, file_path):
        print(f"正在深度侦察文件: {file_path}")
        print("-" * 60)
        try:
            with h5py.File(file_path, 'r') as f:
                # 递归遍历整个文件
                f.visititems(self._inspect_node)
        except FileNotFoundError:
            print("❌ 找不到文件，请检查路径。")
        except Exception as e:
            print(f"❌ 发生错误: {e}")

        print("-" * 60)
        print("侦察结束。")
        
    # ==============================================================================
    # 1. Data Loading
    # ==============================================================================
    def load_data(self):
        file_path = self.config['input_path']
        if not os.path.exists(file_path): raise FileNotFoundError(f"❌ File not found: {file_path}")
        print(f"📖 [Loader] Reading data...")
        with h5py.File(file_path, 'r') as f:
            X = self._read_matrix_smart(f)
            n_cells = X.shape[0]
            obs_df = self._read_obs_smart(f, n_cells)
            var_names = self._read_var_names(f, X.shape[1])
            if len(obs_df) != n_cells and not obs_df.empty: self.adata = anndata.AnnData(X=X)
            else: self.adata = anndata.AnnData(X=X, obs=obs_df)
            if var_names is not None: self.adata.var_names = var_names
            self.adata.var_names_make_unique()
            self._decode_labels(f)
        
        # --- Experiment Mode Logic ---
        if self.config.get('use_experiment'):
            print("   -> Loading experiment embeddings & predictions...")
            try:
                self.custom_embedding = np.loadtxt(self.config['experiment_files']['embeddings'])
                pred_path = self.config['experiment_files']['pred']
                try: self.custom_prediction = np.loadtxt(pred_path, dtype=str)
                except: self.custom_prediction = pd.read_csv(pred_path, header=None).iloc[:,0].astype(str).values
            except Exception as e:
                raise ValueError(f"Failed to load experiment files: {e}")
                
        print(f"✅ [Loader] Done: {self.adata.shape}")

    def _read_matrix_smart(self, f):
        candidates = ['exprs', 'X', 'data', 'counts', 'matrix']
        target_key = next((c for c in candidates if c in f.keys()), None)
        obj = f[target_key]
        if isinstance(obj, h5py.Group):
            mapper = {}
            for k in obj.keys():
                k_low = k.lower()
                if k_low in ['data', 'x']: mapper['data'] = k
                if k_low in ['indices', 'i']: mapper['indices'] = k
                if k_low in ['indptr', 'p']: mapper['indptr'] = k
                if k_low in ['shape', 'dims']: mapper['shape'] = k
            data = np.array(obj[mapper['data']])
            indices = np.array(obj[mapper['indices']])
            indptr = np.array(obj[mapper['indptr']])
            shape = tuple(np.array(obj[mapper['shape']])) if 'shape' in mapper else (indices.max() + 1, len(indptr) - 1)
            if len(indptr) == shape[1] + 1: matrix = sparse.csc_matrix((data, indices, indptr), shape=shape)
            else: matrix = sparse.csc_matrix((data, indices, indptr), shape=(shape[1], shape[0]))
        else: matrix = sparse.csr_matrix(np.array(obj))
        if matrix.shape[0] > matrix.shape[1] and matrix.shape[0] > 5000: return matrix.T
        return matrix

    def _read_obs_smart(self, f, n_cells):
        obs_df = pd.DataFrame()
        if 'obs_names' in f:
            names = np.array(f['obs_names']).astype(str)
            if len(names) == n_cells: obs_df.index = names
        if 'obs' in f and isinstance(f['obs'], h5py.Group):
            for key in f['obs'].keys():
                try:
                    col_data = np.array(f['obs'][key])
                    if len(col_data) == n_cells:
                        if col_data.dtype.kind == 'S': col_data = col_data.astype(str)
                        obs_df[key] = col_data
                except: pass
        label_keys = ['Y', 'label', 'cell_type', 'cluster']
        for k in label_keys:
            if k in f.keys():
                labels = np.array(f[k]).flatten()
                if len(labels) == n_cells:
                    if labels.dtype.kind == 'S': labels = labels.astype(str)
                    obs_df['cell_type'] = labels
                    break
        return obs_df

    def _read_var_names(self, f, n_genes):
        candidates = ['var_names', 'gene_names', 'genes']
        if 'var' in f and 'var_names' in f['var']:
            names = np.array(f['var']['var_names']).astype(str)
            if len(names) == n_genes: return names
        for k in candidates:
            if k in f.keys():
                names = np.array(f[k]).astype(str)
                if len(names) == n_genes: return names
        return None

    def _decode_labels(self, f):
        target_col = 'cell_type'
        if target_col not in self.adata.obs: return
        labels = self.adata.obs[target_col]
        if not (np.issubdtype(labels.dtype, np.number) or labels.str.isnumeric().all()): return
        n_classes = len(np.unique(labels))
        candidates = ['names', 'label_names']
        for k in candidates:
            if k in f.keys():
                names = np.array(f[k]).astype(str)
                if len(names) == n_classes:
                    print(f"   ℹ️ [Loader] Found numeric label mapping: {k}")
                    mapping = {i+int(labels.min()): n for i, n in enumerate(names)}
                    mapping.update({str(i+int(labels.min())): n for i, n in enumerate(names)})
                    self.adata.obs[target_col] = self.adata.obs[target_col].map(mapping).astype(str)
                    break

    # ==============================================================================
    # 2. Preprocessing
    # ==============================================================================
    def run_preprocessing(self, n_top_genes=2000):
        print("⚙️ [Preprocess] Starting preprocessing...")
        sc.pp.filter_genes(self.adata, min_cells=3)
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        self.adata.raw = self.adata 
        
        # --- Mode Switching Logic ---
        if self.config.get('use_experiment'):
            # Experiment Mode: Inject external embeddings and predictions
            print("   -> 🧪 Experiment Mode: Injecting Embeddings & Predictions")
            self.adata.obsm['X_deep_emb'] = self.custom_embedding
            self.adata.obs['grouper'] = self.custom_prediction
            self.adata.obs['grouper'] = self.adata.obs['grouper'].astype('category')
            # Compute neighbors on custom embedding
            sc.pp.neighbors(self.adata, use_rep='X_deep_emb', n_neighbors=15)
            sc.tl.umap(self.adata)
        else:
            # Standard Mode: Standard Scanpy pipeline
            print("   -> 📦 Standard Mode: HVG -> PCA -> Neighbors -> Leiden")
            target = self.config.get('target_col', 'cell_type')
            if target not in self.adata.obs: target = 'leiden'
            
            sc.pp.highly_variable_genes(self.adata, n_top_genes=n_top_genes)
            self.adata = self.adata[:, self.adata.var.highly_variable]
            sc.pp.scale(self.adata, max_value=10)
            sc.tl.pca(self.adata)
            sc.pp.neighbors(self.adata)
            sc.tl.umap(self.adata)
            
            if target == 'leiden': sc.tl.leiden(self.adata)
            self.adata.obs['grouper'] = self.adata.obs[target]
        print("✅ [Preprocess] Done")

    # ==============================================================================
    # 3. Differential Expression Analysis
    # ==============================================================================
    def find_markers(self):
        print(f"🧮 [Analysis] Finding markers...")
        if self.adata.obs['grouper'].dtype.name != 'category':
            self.adata.obs['grouper'] = self.adata.obs['grouper'].astype('category')
        try:
            sc.tl.rank_genes_groups(self.adata, groupby='grouper', method='wilcoxon', pts=True)
        except:
            sc.tl.rank_genes_groups(self.adata, groupby='grouper', method='t-test', pts=True)

        raw_df = sc.get.rank_genes_groups_df(self.adata, group=None)
        df = raw_df[(raw_df['pvals_adj'] < 0.05) & (raw_df['logfoldchanges'] > 1)].copy()
        if len(df) == 0: df = raw_df.nlargest(200, 'logfoldchanges').copy()

        df = df.rename(columns={'names': 'Gene', 'logfoldchanges': 'Log2FC', 'pvals_adj': 'Adj_P', 'group': 'Cell_Type'})
        self.marker_df = df
        self.marker_df.to_csv(os.path.join(self.data_dir, f"{self.config['dataset_name']}_markers.csv"), index=False)
        print(f"💾 [Analysis] Markers saved")

    # ==============================================================================
    # 4. Plotting
    # ==============================================================================
    def plot_markers(self):
        print(f"🎨 [Plotting] Generating plots...")
        prefix = self.config['dataset_name']
        specific_genes = self.config.get('specific_genes', None)
        
        try: sc.pl.rank_genes_groups_dotplot(self.adata, n_genes=5, groupby='grouper', standard_scale='var', save=f'_{prefix}_dotplot.png')
        except: pass
        try: sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=3, groupby='grouper', save=f'_{prefix}_violin.png')
        except: pass
        sc.pl.umap(self.adata, color='grouper', title=f'{prefix} Clusters', save=f'_{prefix}_umap.png')

        if specific_genes:
            self._plot_manual_heatmap(specific_genes, prefix)
            self._plot_manual_features(specific_genes, prefix)
        else:
            try:
                n_groups = len(self.adata.obs['grouper'].unique())
                dynamic_height = max(8, n_groups * 5 * 0.25)
                sc.pl.rank_genes_groups_heatmap(self.adata, n_genes=5, groupby='grouper', cmap='plasma', dendrogram=True, standard_scale='var', figsize=(12, dynamic_height), save=f'_{prefix}_heatmap.png')
            except: pass
            self._plot_auto_features(prefix)
        print(f"✅ [Plotting] Done")

    def _plot_manual_heatmap(self, specific_genes, prefix):
        final_var_names = specific_genes
        use_raw = False
        flat_genes = []
        if isinstance(specific_genes, dict):
            for genes in specific_genes.values(): flat_genes.extend(genes)
        else: flat_genes = specific_genes
        if any(g not in self.adata.var_names for g in flat_genes) and self.adata.raw: use_raw = True
        try:
            dynamic_height = max(8, len(flat_genes) * 0.3)
            sc.pl.heatmap(self.adata, var_names=final_var_names, groupby='grouper', use_raw=use_raw, cmap='magma', standard_scale='var', figsize=(10, dynamic_height), save=f'_{prefix}_paper_heatmap.png')
        except: pass

    def _plot_manual_features(self, specific_genes_dict, prefix):
        genes, titles = [], []
        if isinstance(specific_genes_dict, dict):
            for grp, glist in specific_genes_dict.items():
                if glist: genes.append(glist[0]); titles.append(f"{glist[0]}\n({grp})")
        if genes:
            sc.pl.umap(self.adata, color=genes, title=titles, ncols=4, use_raw=(self.adata.raw is not None), color_map='viridis', save=f'_{prefix}_paper_features.png')
            self._save_single_features(genes, titles, prefix, (self.adata.raw is not None))

    def _plot_auto_features(self, prefix):
        genes, titles = [], []
        for g in self.adata.obs['grouper'].unique():
            gl = self.marker_df[self.marker_df['Cell_Type'] == g]['Gene'].head(1).tolist()
            if gl: genes.append(gl[0]); titles.append(f"{gl[0]}\n({g})")
        if genes:
            sc.pl.umap(self.adata, color=genes, title=titles, ncols=4, use_raw=(self.adata.raw is not None), color_map='viridis', save=f'_{prefix}_auto_features.png')
            self._save_single_features(genes, titles, prefix, (self.adata.raw is not None))

    def _save_single_features(self, genes, titles, prefix, use_raw):
        sdir = os.path.join(self.fig_dir, 'Single_Features')
        os.makedirs(sdir, exist_ok=True)
        old = sc.settings.figdir
        sc.settings.figdir = sdir
        for i, g in enumerate(genes):
            try: sc.pl.umap(self.adata, color=g, title=titles[i], use_raw=use_raw, color_map='viridis', show=False, save=f'_{prefix}_{g.replace("/","_")}.png')
            except: pass
        sc.settings.figdir = old

    # ==============================================================================
    # 5. Enrichment Analysis (GO + KEGG + Sankey)
    # ==============================================================================
    def run_enrichment(self, top_n=100):
        print(f"🧬 [Enrichment] Starting enrichment analysis...")
        go_dir = os.path.join(self.data_dir, 'Enrichment_GO')
        kegg_dir = os.path.join(self.data_dir, 'Enrichment_KEGG')
        net_dir = os.path.join(self.fig_dir, 'Enrichment_Network') 
        for d in [go_dir, kegg_dir, net_dir]: os.makedirs(d, exist_ok=True)
        
        organism = self.config.get('organism', 'Mouse')
        kegg_lib = 'KEGG_2019_Mouse' if organism.lower() == 'mouse' else 'KEGG_2021_Human'
        go_libs = ['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021']
        go_labels = ['BP', 'CC', 'MF']

        for ct in self.marker_df['Cell_Type'].unique():
            sub_marker_df = self.marker_df[self.marker_df['Cell_Type'] == ct]
            genes = sub_marker_df['Gene'].head(top_n).tolist()
            genes_upper = [g.upper() for g in genes]
            if len(genes) < 5: continue
            
            safe_name = str(ct).replace('/', '_').replace(' ', '_')

            # 1. GO
            go_results = []
            for lib, label in zip(go_libs, go_labels):
                try:
                    enr = gp.enrichr(gene_list=genes_upper, gene_sets=lib, organism=organism, outdir=None)
                    res = enr.results
                    res = res[res['Adjusted P-value'] < 0.05].copy()
                    if not res.empty:
                        if 'Overlap' in res.columns: res['Count'] = res['Overlap'].apply(lambda x: int(str(x).split('/')[0]))
                        else: res['Count'] = res['Genes'].apply(lambda x: len(str(x).split(';')))
                        res = res.sort_values('Adjusted P-value', ascending=True).head(10)
                        res['Group'] = label
                        go_results.append(res)
                except: pass

            if go_results:
                combined_go = pd.concat(go_results)
                combined_go.to_csv(os.path.join(go_dir, f"{safe_name}_GO_combined.csv"), index=False)
                self._plot_go_faceted(combined_go, ct, os.path.join(go_dir, f"{safe_name}_GO_barplot.png"))

            # 2. KEGG & 3. Sankey
            try:
                enr_kegg = gp.enrichr(gene_list=genes_upper, gene_sets=kegg_lib, organism=organism, outdir=None)
                res_kegg = enr_kegg.results
                res_kegg = res_kegg[res_kegg['Adjusted P-value'] < 0.05].copy()
                res_kegg = res_kegg[~res_kegg['Term'].str.contains('disease|infection|cancer|cardiomyopathy', case=False, regex=True)]
                
                if not res_kegg.empty:
                    res_kegg['Hit'] = res_kegg['Overlap'].apply(lambda x: int(str(x).split('/')[0]))
                    res_kegg['Total'] = res_kegg['Overlap'].apply(lambda x: int(str(x).split('/')[1]))
                    res_kegg['GeneRatio'] = res_kegg['Hit'] / res_kegg['Total']
                    res_kegg_top = res_kegg.sort_values('Adjusted P-value', ascending=True).head(15)
                    
                    res_kegg_top.to_csv(os.path.join(kegg_dir, f"{safe_name}_KEGG.csv"), index=False)
                    self._plot_kegg_bubble(res_kegg_top, ct, os.path.join(kegg_dir, f"{safe_name}_KEGG_bubble.png"))
                    
                    self._plot_kegg_sankey(res_kegg_top.head(8), sub_marker_df, ct, os.path.join(net_dir, f"{safe_name}_Sankey.png"))
                    
            except: pass
        print(f"✅ [Enrichment] Done")

    def _plot_go_faceted(self, df, title, save_path):
        try:
            df['Term'] = df['Term'].apply(lambda x: x.split(' (GO:')[0])
            df['Group'] = pd.Categorical(df['Group'], categories=['BP', 'CC', 'MF'], ordered=True)
            df = df.sort_values(['Group', 'Adjusted P-value'], ascending=[True, False])
            groups = df['Group'].unique()
            fig, axes = plt.subplots(nrows=len(groups), ncols=1, figsize=(8, max(4, 3 * len(groups))), sharex=False, constrained_layout=True)
            if len(groups) == 1: axes = [axes]
            norm = mcolors.LogNorm(vmin=df['Adjusted P-value'].min(), vmax=df['Adjusted P-value'].max())
            cmap = plt.cm.Spectral_r
            for ax, grp in zip(axes, groups):
                sub = df[df['Group'] == grp]
                ax.barh(sub['Term'], sub['Count'], color=cmap(norm(sub['Adjusted P-value'].values)))
                ax.set_title(grp, loc='right', fontsize=12, fontweight='bold', color='gray')
                ax.grid(axis='x', linestyle='--', alpha=0.3)
            fig.supxlabel('Gene Count', fontsize=12)
            fig.suptitle(f'{title}\nGO Enrichment', fontsize=14)
            cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=axes, orientation='vertical', fraction=0.05, pad=0.02)
            cbar.set_label('Adjusted P-value')
            fig.savefig(save_path, bbox_inches='tight', dpi=300)
            plt.close(fig)
        except: plt.close()

    def _plot_kegg_bubble(self, df, title, save_path):
        try:
            df['Term'] = df['Term'].apply(lambda x: x.split(' (')[0] if '(' in x else x)
            df = df.sort_values('GeneRatio', ascending=True)
            fig, ax = plt.subplots(figsize=(8, max(5, len(df) * 0.4)))
            norm = mcolors.LogNorm(vmin=df['Adjusted P-value'].min(), vmax=df['Adjusted P-value'].max())
            scatter = ax.scatter(x=df['GeneRatio'], y=df['Term'], s=df['Hit'] * 70, c=df['Adjusted P-value'], cmap=plt.cm.RdYlBu_r, norm=norm, alpha=0.9, edgecolors='black', linewidth=0.5)
            for i in range(len(df)):
                row = df.iloc[i]
                ax.text(row['GeneRatio'], row['Term'], str(int(row['Hit'])), ha='center', va='center', fontsize=8, fontweight='bold', color='black')
            ax.set_xlabel("GeneRatio")
            ax.set_title(f"{title}\nKEGG")
            plt.colorbar(scatter, ax=ax).set_label('Adj P')
            fig.savefig(save_path, bbox_inches='tight', dpi=300)
            plt.close(fig)
        except: plt.close()

    def _plot_kegg_sankey(self, kegg_df, marker_df, title, save_path):
        try:
            pathways = []
            all_genes_set = set()
            gene_fc_dict = dict(zip(marker_df['Gene'].str.upper(), marker_df['Log2FC']))
            
            for _, row in kegg_df.iterrows():
                p_name = row['Term'].split(' (')[0]
                pval = row['Adjusted P-value']
                genes = [g.upper() for g in str(row['Genes']).split(';') if g.upper() in gene_fc_dict]
                if not genes: continue
                pathways.append({'name': p_name, 'pval': pval, 'genes': genes, 'count': len(genes)})
                for g in genes: all_genes_set.add(g)
            
            pathways.sort(key=lambda x: x['pval'])
            gene_list = sorted(list(all_genes_set), key=lambda x: gene_fc_dict.get(x, 0), reverse=True)
            if not pathways or not gene_list: return

            fig_height = max(8, len(gene_list) * 0.25)
            fig, ax = plt.subplots(figsize=(12, fig_height))
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 100)
            
            p_col_x, p_col_w = 5, 20
            g_col_x, g_col_w = 75, 20
            ribbon_start_x = p_col_x + p_col_w
            ribbon_end_x = g_col_x
            
            p_vals = [p['pval'] for p in pathways]
            norm_p = mcolors.LogNorm(vmin=min(p_vals), vmax=max(p_vals))
            cmap_p = plt.cm.RdYlBu_r

            fc_vals = [gene_fc_dict[g] for g in gene_list]
            max_fc = max(abs(min(fc_vals)), abs(max(fc_vals))) if fc_vals else 1
            norm_fc = mcolors.Normalize(vmin=-max_fc, vmax=max_fc)
            cmap_fc = plt.cm.bwr

            total_gene_slots = sum(p['count'] for p in pathways)
            total_genes_unique = len(gene_list)
            y_avail = 90
            y_start = 5
            spacing = 1
            slot_height_unit = (y_avail - (len(pathways)-1)*spacing) / total_gene_slots
            g_box_h = (y_avail - (total_genes_unique-1)*spacing*0.5) / total_genes_unique
            
            g_y_coords = {}
            current_g_y = y_start
            for g in gene_list:
                g_y_coords[g] = {'bottom': current_g_y, 'top': current_g_y + g_box_h, 'center': current_g_y + g_box_h/2}
                current_g_y += g_box_h + spacing*0.5

            def sigmoid(t): return 1 / (1 + np.exp(-t))
            x_ribbon = np.linspace(ribbon_start_x, ribbon_end_x, 100)
            t = np.linspace(-6, 6, 100)
            sig = sigmoid(t)
            current_y = y_start
            
            for p in pathways:
                p_h = p['count'] * slot_height_unit
                p_internal_y = current_y
                ribbon_color = cmap_p(norm_p(p['pval']))
                
                for g_name in p['genes']:
                    y_start_top = p_internal_y + slot_height_unit
                    y_end_top = g_y_coords[g_name]['top']
                    y_end_bottom = g_y_coords[g_name]['bottom']
                    
                    y_top_curve = y_start_top + (y_end_top - y_start_top) * sig
                    y_bottom_curve = p_internal_y + (y_end_bottom - p_internal_y) * sig
                    
                    ax.fill_between(x_ribbon, y_bottom_curve, y_top_curve, color=ribbon_color, alpha=0.3, lw=0)
                    p_internal_y += slot_height_unit
                
                rect = Rectangle((p_col_x, current_y), p_col_w, p_h, facecolor=ribbon_color, edgecolor='black', lw=1, zorder=5)
                ax.add_patch(rect)
                ax.text(p_col_x - 1, current_y + p_h/2, p['name'], ha='right', va='center', fontsize=10, fontweight='bold')
                current_y += p_h + spacing

            for g_name in gene_list:
                coords = g_y_coords[g_name]
                color = cmap_fc(norm_fc(gene_fc_dict[g_name]))
                rect = Rectangle((g_col_x, coords['bottom']), g_col_w, g_box_h, facecolor=color, edgecolor='black', lw=0.5, zorder=5)
                ax.add_patch(rect)
                ax.text(g_col_x + g_col_w + 1, coords['center'], g_name, ha='left', va='center', fontsize=9)

            ax.axis('off')
            plt.title(f"{title} - Gene-Pathway Sankey", fontsize=16, pad=20)
            
            sm_p = plt.cm.ScalarMappable(cmap=cmap_p, norm=norm_p)
            cbar_p_ax = fig.add_axes([0.15, 0.05, 0.2, 0.02])
            plt.colorbar(sm_p, cax=cbar_p_ax, orientation='horizontal').set_label('Pathway Adj. P-value')
            
            sm_fc = plt.cm.ScalarMappable(cmap=cmap_fc, norm=norm_fc)
            cbar_fc_ax = fig.add_axes([0.65, 0.05, 0.2, 0.02])
            plt.colorbar(sm_fc, cax=cbar_fc_ax, orientation='horizontal').set_label('Gene Log2FC')

            plt.savefig(save_path, bbox_inches='tight', dpi=300)
            plt.close()
        except Exception as e:
            print(f"      ⚠️ Sankey plot failed: {e}")
            plt.close()

    # ==============================================================================
    # 6. Trajectory Analysis (Updated for Larger Fonts)
    # ==============================================================================
    def run_trajectory(self):
        root = self.config.get('root_cell_type')
        if not root: return
        print(f"⏳ [Trajectory] Starting trajectory analysis (Root: {root})...")
        try:
            matches = np.flatnonzero(self.adata.obs['grouper'].astype(str) == str(root))
            if len(matches) > 0: self.adata.uns['iroot'] = matches[0]
            else: return
        except: return
        
        # Ensure diffmap is run before dpt
        sc.tl.diffmap(self.adata)
        sc.tl.dpt(self.adata)
        
        cat_map = {name: str(i) for i, name in enumerate(self.adata.obs['grouper'].cat.categories)}
        legend_map = {i: name for i, name in enumerate(self.adata.obs['grouper'].cat.categories)}
        self.adata.obs['grouper_code'] = self.adata.obs['grouper'].map(cat_map).astype('category')
        
        sc.tl.paga(self.adata, groups='grouper_code')
        sc.settings.figdir = self.fig_dir
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # 【关键修改】使用 rc_context 临时增大 PAGA 图节点上的数字字体
        with plt.rc_context({'font.size': 14, 'font.weight': 'bold'}):
            try:
                sc.pl.paga(self.adata, threshold=0.01, layout='fr', random_state=42, 
                           edge_width_scale=2.0, show=False, ax=ax, fontoutline=2, node_size_scale=1.5)
            except:
                sc.pl.paga(self.adata, threshold=0.01, show=False, ax=ax)
            
        colors = self.adata.uns.get('grouper_code_colors', sc.pl.palettes.default_20)
        patches = [mpatches.Patch(color=colors[i] if i < len(colors) else 'gray', label=f"{i}: {name}") for i, name in legend_map.items()]
        
        # 【关键修改】增大图例字体和标题字体
        ax.legend(handles=patches, title="Cell Types", 
                  bbox_to_anchor=(1.02, 0.5), loc='center left', 
                  frameon=False, 
                  fontsize=12,       # 图例项字体大小
                  title_fontsize=14) # 图例标题字体大小
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.fig_dir, f'_{self.config["dataset_name"]}_paga_clean.png'), bbox_inches='tight', dpi=300)
        plt.close(fig)
        
        sc.pl.umap(self.adata, color=['dpt_pseudotime', 'grouper'], save=f'_{self.config["dataset_name"]}_traj.png')
        del self.adata.obs['grouper_code']
        print(f"✅ [Trajectory] Done")

    # ==============================================================================
    # 7. Saving Data
    # ==============================================================================
    def save_data(self):
        path = os.path.join(self.data_dir, f"{self.config['dataset_name']}_processed.h5ad")
        self.adata.write(path)
        print(f"💾 [Save] {path}")