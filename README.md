

# SCPipeline: 全流程单细胞转录组分析与可视化框架

SCPipeline：仅需8行代码，跑通单细胞转录组数据分析基础实验！

**SCPipeline** 是一个基于 [Scanpy](https://scanpy.readthedocs.io/) 封装的高级单细胞分析框架。它专为生成 **出版级Publication-Ready** 图表而设计，支持从原始数据加载到拟时序分析的全流程处理。

该框架具有独特的**双模式架构**，既支持标准的 Scanpy 分析流程，也能无缝集成深度学习聚类模型（Deep Clustering）的输出结果。

<div align="center"> <img src="./abs.png" width="800"/></div>
---

## 🌟 核心特性 (Key Features)

### 1. 双模式架构 (Dual-Mode Architecture)

* **📦 标准模式 (Standard Mode)**:
* 执行标准的单细胞分析流程：HVG 选择 -> PCA 降维 -> Neighbors 计算 -> Leiden 聚类。
* 适用于常规数据集的探索性分析。


* **🧪 实验模式 (Experiment Mode)**:
* 专为深度学习模型设计。直接加载外部计算的 **Embeddings** 和 **预测标签 (Predictions)**。
* 跳过 PCA 和 Leiden，直接基于深度特征构建邻接图（Neighborhood Graph）并生成 UMAP/PAGA。
* 用于验证自研算法在下游生物学分析中的表现。



### 2. 高级可视化 (Advanced Visualization)

* **桑基图 (Sankey Diagram)**: 全新设计的基因-通路映射图。摒弃传统的点状图，使用贝塞尔曲线（Bezier Curves）连接左侧通路与右侧基因，直观展示 Crosstalk。
* **论文级热力图 (Paper-Style Heatmap)**: 自动根据基因数量动态调整画布高度，支持从 Raw 数据回溯丢失的 Marker 基因。
* **PAGA 轨迹图**: 优化布局算法，消除多余留白，自动放大节点标签和图例字体，确保视觉协调。
* **Feature Plots**: 支持批量自动导出每个 Marker 基因的独立 UMAP 图（保存至 `Single_Features` 子目录），同时也支持生成组合大图。

### 3. 自动化与模块化

* **配置驱动**: 所有参数通过 `main.py` 中的配置字典控制，无需修改底层代码。
* **自动 Marker 筛选**: 内置 Wilcoxon 检验，自动计算并导出差异基因 CSV。
* **功能富集**: 集成 `GSEApy`，自动进行 GO 和 KEGG 分析，并生成气泡图和柱状图。

---

## 🛠️ 安装指南 (Installation)

### 1. 环境要求

建议使用 Anaconda 创建独立的虚拟环境：

```bash
conda create -n sc_pipeline python=3.9
conda activate sc_pipeline

```

### 2. 依赖包安装

请确保安装以下 Python 库：

```bash
pip install scanpy pandas numpy matplotlib seaborn h5py scipy gseapy leidenalg anndata

```

*注意：`leidenalg` 是 Scanpy 进行 Leiden 聚类所必需的。*

---

## 📂 目录结构 (Directory Structure)

在使用本 Pipeline 前，建议的数据组织方式如下：

```text
Project_Root/
├── main.py                 # 入口文件 (用户配置)
├── sc_pipeline.py          # 核心代码库 (无需修改)
├── 数据集/
│   └── Adam/
│       ├── data.h5         # 原始表达矩阵 (H5格式)
│       ├── embeddings.txt  # (可选) 实验模式所需的嵌入向量
│       └── pred.txt        # (可选) 实验模式所需的聚类标签
└── Results/                # 输出目录 (自动生成)
    └── Adam/
        ├── Data/           # 存放处理后的 h5ad 和 CSV 表格
        └── Figures/        # 存放所有图片
            └── Single_Features/  # 存放单基因 UMAP 图

```

---

## 🚀 使用指南 (Usage)

### 1. 配置 `main.py`

在 `main.py` 中，通过修改 `datasets` 列表来定义分析任务。

#### 场景 A：运行标准流程

```python
datasets = [
    {
        'dataset_name': 'Adam',
        'input_path': '数据集/Adam/data.h5', 
        'output_base': 'Results',
        'target_col': 'cell_ontology_class', # 指定原本的细胞类型列（如果有）
        'root_cell_type': 'ureteral cell',   # 指定拟时序分析的起点细胞
        'organism': 'Mouse',                 # 用于富集分析 ('Mouse' or 'Human')
        'use_experiment': False              # 关闭实验模式
    }
]

```

#### 场景 B：运行实验模式 (深度聚类)

```python
datasets = [
    {
        'dataset_name': 'Adam_DeepClust',
        'input_path': '数据集/Adam/data.h5', 
        'output_base': 'Results',
        'use_experiment': True,              # 开启实验模式
        'experiment_files': {
            'embeddings': '数据集/Adam/embeddings.txt', # 纯数字矩阵
            'pred': '数据集/Adam/pred.txt'              # 预测标签 (单列)
        },
        'root_cell_type': 'ureteral cell', 
        'organism': 'Mouse'
    }
]

```

#### 场景 C：手动指定 Marker 基因 (复刻论文图)

如果你有感兴趣的特定基因列表，可以添加到配置中。**如果不添加此项，程序会自动计算 Top 基因并绘图。**

```python
        # ... 在 config 字典中添加 ...
        'specific_genes': {
            'Group_A': ['Gene1', 'Gene2'],
            'Group_B': ['Gene3', 'Gene4']
        }

```

### 2. 运行

```bash
python main.py

```

---

## 📊 输出结果说明 (Output)

运行结束后，`Results/dataset_name/` 目录下将生成：

### 📁 Data/

* `*_processed.h5ad`: 保存了所有分析结果（降维、聚类、拟时序）的 Anndata 对象。
* `*_markers.csv`: 差异基因分析结果表。
* `Enrichment_GO/*.csv` & `Enrichment_KEGG/*.csv`: 富集分析原始数据。

### 📁 Figures/

* **基础图**:
* `_dotplot.png`: 气泡图。
* `_violin.png`: 小提琴图。
* `_umap.png`: 聚类 UMAP 图。


* **高级图**:
* `_heatmap.png`: 动态高度的热力图。
* `_auto_features.png` / `_paper_features.png`: Marker 基因在 UMAP 上的分布大图。


* **富集分析**:
* `_KEGG_bubble.png`: KEGG 气泡图（数字标注在气泡内）。
* `_GO_barplot.png`: GO 分面柱状图。
* `Enrichment_Network/*_Sankey.png`: **基因-通路桑基图**。


* **拟时序**:
* `_paga_clean.png`: 优化布局的 PAGA 拓扑图。
* `_trajectory.png`: 伪时间（Pseudotime）UMAP 映射图。


* **Single_Features/**: 包含每个 Marker 基因单独的 UMAP 分布图。

---

## ⚙️ 核心类设计 (Design Details)

代码逻辑封装在 `SCPipeline` 类中 (`sc_pipeline.py`)：

1. **Loader**: 智能识别 H5 内部结构（`data`, `indices`, `indptr`），兼容多种存储格式。
2. **Preprocessing**:
* `use_experiment=True`: 将外部 Embedding 注入 `adata.obsm['X_deep_emb']`，并以此计算 Neighbors。
* `use_experiment=False`: 执行标准 PCA + Neighbors。


3. **Plotting Engine**:
* 使用了 `matplotlib.patches` (Rectangle, Polygon) 手绘桑基图，不再依赖复杂的第三方网络库，确保了极高的定制化自由度（如颜色映射、布局控制）。
* 使用 `plt.rc_context` 动态控制 PAGA 图的字体大小。



---

✉️ 作者信息 (Author)
开发者: Zhongyang Zhou (周中阳)

单位: Chongqing Normal University (重庆师范大学)

研究方向: AI4Bios、数据挖掘、深度学习

Email: zhouzhongyang@163.com