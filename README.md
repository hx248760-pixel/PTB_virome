# PTB Virome Project / 肺结核病毒组项目

**English**  
Reproducibility package for the Pulmonary Tuberculosis (PTB) virome project.
This pipeline covers bacterial differential abundance, virome diversity and
composition, viral differential abundance, viral–host linkage, functional
enrichment, co-occurrence network analysis, random forest classification,
and external cohort validation.

**中文**  
肺结核（PTB）病毒组项目的可复现分析流程。内容包括细菌差异丰度、病毒组
多样性与组成、病毒差异丰度、病毒–宿主关联、功能富集、共现网络分析、
随机森林分类以及外部队列验证。

---

## 1. Repository structure / 仓库结构

**English**
.
├── 01_bac_analysis.R # Bacterial differential abundance (MaAsLin2)
├── 02_virome_analysis.R # Alpha / Beta diversity, family barplot
├── 03_differential_analysis.R # Viral differential abundance, volcano, pie, host
├── 04_functional_analysis.R # KO functional enrichment (Fisher)
├── 05_prepare_inputf.R # Prepare correlation inputs + run calc_corr.py
├── 05_network_analysis.R # Co-occurrence network analysis & visualization
├── 06_random_forest.R # Random forest (importance, K-fold, ROC, Top20)
├── 06_random_forest_external.R # External validation of RF models
├── environment.yml # Conda environment for R + Python
├── data/
│ ├── metadata.tsv
│ ├── votu_abundance.tsv
│ ├── virus_family_abundance.tsv
│ ├── bacterial_species_abundance.tsv
│ ├── votu_ko_abundance.tsv
│ ├── annotations/
│ │ ├── votu_ko_mapping.tsv
│ │ ├── ko_description.tsv
│ │ ├── vOTU.family.tsv
│ │ ├── vOTU.host.tax.genus.tsv
│ │ └── network/
│ │ ├── hc_node.csv
│ │ └── tb_node.csv
│ └── external/ # (Optional) External validation inputs
│ ├── votu.profile
│ ├── species_relative_abundance.tsv
│ └── sample_info.txt
├── script/
│ ├── config.R # Centralized paths & parameters
│ ├── calculate_alpha.R
│ ├── plot_aplha.R
│ ├── plot_beta.R
│ ├── plot_family_barplot.R
│ ├── plot_virf_host.R
│ ├── plot_roc.R
│ ├── calc_corr.py # Spearman correlation (multi-thread)
│ ├── rf.py # Random forest engine
│ ├── run_rf.sh # RF pipeline shell script
│ └── run_rf_external.sh # RF external validation shell script
└── results/ # Auto-generated
├── figures/
├── tables/
├── maaslin/
└── random_forest/

---

## 2. Environment setup / 环境安装

**English**

We recommend using conda / mamba:

```bash
conda env create -f environment.yml
conda activate ptb-virome
The environment includes:

R 4.3 with tidyverse, vegan, pROC, ggpubr, patchwork,
igraph, ggraph, ggtext, ggrepel, and Maaslin2.

Python 3.10 with pandas, numpy, scipy, scikit-learn, dask.

GNU parallel for the random forest shell pipeline.
3. Configuration / 配置
English

All file paths and parameters are centralized in script/config.R.
Key parameters:

Variable	Description
METADATA_FILE	Sample metadata (must contain Sample, Group, Sex, Age)
VOTU_FILE	vOTU abundance matrix (features × samples)
BACTERIAL_FILE	Bacterial species abundance matrix
VIRUS_FAMILY_FILE	Viral family abundance matrix
MAASLIN_*	MaAsLin2 parameters (min abundance, prevalence, covariates, reference)
OCCURRENCE_THRESHOLD	Prevalence filter for correlation analysis
CORRELATION_THRESHOLD	Spearman |r| cutoff
FDR_THRESHOLD	FDR cutoff
RF_*	Random forest configuration (threads, K, top-N grid, seeds)
RF_EXTERNAL_*	External validation file paths
GROUP_LEVELS should be c("Health", "Patient") and match the values
used in metadata.tsv.

中文

所有文件路径与参数均集中在 script/config.R。关键参数：

变量	说明
METADATA_FILE	样本元数据（需包含 Sample、Group、Sex、Age）
VOTU_FILE	vOTU 丰度矩阵（特征 × 样本）
BACTERIAL_FILE	细菌物种丰度矩阵
VIRUS_FAMILY_FILE	病毒科丰度矩阵
MAASLIN_*	MaAsLin2 参数（最小丰度、流行率、协变量、参考水平）
OCCURRENCE_THRESHOLD	相关性分析的出现率阈值
CORRELATION_THRESHOLD	Spearman |r| 阈值
FDR_THRESHOLD	FDR 阈值
RF_*	随机森林配置（线程数、K 折、top-N 网格、随机种子）
RF_EXTERNAL_*	外部验证文件路径
GROUP_LEVELS 应为 c("Health", "Patient")，并与 metadata.tsv 中的取值一致。

4. Usage / 使用流程
Run the scripts in numerical order. All results are written under
results/ (figures, tables, MaAsLin2 outputs, RF outputs).

按编号顺序运行脚本，所有结果输出到 results/ 下（图片、表格、MaAsLin2 输出、随机森林输出）。
