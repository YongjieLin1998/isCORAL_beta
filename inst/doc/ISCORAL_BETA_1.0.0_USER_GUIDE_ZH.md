# isCORAL beta 1.0.0 中文使用说明

## 1. 软件定位

isCORAL beta 1.0.0 用于从单细胞转录组中识别稳定的转录记忆信号，并推断细胞间的 pseudoclone（伪克隆）关联。该框架可用于探索仅依赖转录组的 **in-silico lineage tracing**。本测试版以 isCORAL 6.1.5 的功能实现为基线，并将正式主扫描入口统一为 `isCORAL()`。

需要区分以下概念：

- **pseudoclone**：由转录组推断得到的连续、概率性细胞关联，不等同于实验测得的真实 clone；
- **`CORAL_cell_embed`**：由 memory genes 投影得到的细胞表达空间；
- **`CORAL_clone_embed`**：由 pseudoclone 概率加权得到的克隆关联空间；
- **Metaclone**：在 `CORAL_clone_embed` 中构建 SNN 图并进行社区发现后得到的离散群落，不等同于 pseudoclone 或实验 lineage label。

论文在其明确的数据集、指标和比较方法范围内报告了领先的 clone-reconstruction benchmark。本文档用于说明软件接口和标准工作流，不构成对该 benchmark 或 SOTA 结论的新一轮独立验证。

## 2. 系统要求与安装

### 2.1 平台要求

核心扫描通过 `callr` 启动独立 R 进程，并在该进程中使用 fork-based `doMC` 后端。因此：

- 主扫描面向 Linux、macOS、HPC 或其他类 Unix 环境；
- beta 1.0.0 没有为原生 Windows 提供等价的 PSOCK 主扫描后端；
- 已经包含扫描结果的 Seurat 对象，可在依赖满足的环境中继续运行 embedding、Metaclone 和绘图函数。

主扫描和 embedding 都会构建细胞×细胞矩阵。数据集较大时，应提前评估内存需求；`max_cells_lasso` 只减少 Elastic Net 的训练量，不减少前面的全细胞能量距离矩阵计算。

### 2.2 加载或安装

```r
BASE_DIR <- "/path/to/project"
PKG_DIR <- file.path(BASE_DIR, "packages", "isCORAL_beta_1.0.0")

devtools::load_all(PKG_DIR)
utils::packageVersion("isCORAL")
```

正式安装：

```r
devtools::install(
  PKG_DIR,
  dependencies = TRUE,
  upgrade = "never"
)

library(isCORAL)
utils::packageVersion("isCORAL")
```

`leidenbase` 是位于 `Suggests` 中的可选依赖。如需默认的 Leiden Metaclone 聚类，建议预先安装：

```r
install.packages("leidenbase")
```

`FindMetaclone()` 不会在运行过程中联网或安装软件。若以 `algorithm = 4` 请求 Leiden，但无法加载 `leidenbase`，或当前 Seurat 早于 5.2.0、尚未采用 `leidenbase` 后端，函数会明确警告并以 `algorithm = 1` 回退到 Louvain。

## 3. 输入数据

输入必须是 Seurat 对象，并具有唯一的细胞条形码。推荐使用 `RNA` assay。

### 3.1 原始 counts

`input_mode = "raw_counts"` 是默认路径。函数读取所选 assay 的 `counts`，可选执行 UMI downsampling，调用 `Seurat::NormalizeData()`，并使用生成的 `data` 进行主扫描。

### 3.2 预归一化数据

仅有 log-normalized `data`、缺少 counts 的对象应显式使用：

```r
sample <- isCORAL(
  seurat_obj = sample,
  assay = "RNA",
  input_mode = "pre_normalized_data",
  pre_normalized_layer = "data",
  species = "human"
)
```

该路径不再次调用 `NormalizeData()`，也不会验证输入是否真的经过恰当归一化。使用者应自行确认矩阵来源和数值含义。此模式下如果把 `downsampleUMI` 设为非 `NULL`，该参数会被忽略并产生警告。

### 3.3 物种

请显式设置 `species = "human"` 或 `species = "mouse"`。该参数只控制线粒体基因前缀和 cell-cycle gene symbol 的处理，不会转换表达矩阵的物种。

## 4. 标准工作流

```text
Seurat object
    ↓
isCORAL()
    ├── PlotMemoryGenes()                   # memory-gene 排名和候选基因选择
    ├── Assign_deterministic_pseudoclone()  # 共识矩阵确定性层次聚类
    │       └── simulate_clonal_expression()# 指定基因的 isCORAL 克隆表达图
    ├── PlotPseudoclone()                   # 原始 UMAP 上的 seed-level 可视化
    └── CORALembedding()
            ├── PlotPseudoclone()          # 可选的 seed-level 可视化
            └── FindMetaclone()
                    └── Metaclone 可视化与验证
```

所有返回 Seurat 对象的函数都必须将结果重新赋给对象。下文各节给出了包含参数的完整调用示例；请使用 `sample <- function_name(seurat_obj = sample, ...)` 的赋值形式，不要只调用函数而丢弃返回值。

## 5. 主扫描：`isCORAL()`

### 5.1 推荐的 raw-counts 调用

```r
set.seed(20251009)

sample <- isCORAL(
  seurat_obj = sample,
  assay = "RNA",
  input_mode = "raw_counts",
  check_raw_counts = TRUE,
  results_dir = "iscoral_results_raw",
  species = "human",
  random_seed = 20251009,
  downsampleUMI = NULL,
  bin_size = 200,
  n_repeats = 100,
  parallel_cores = 25,
  max_cells = NULL,
  genes_to_remove = NULL,
  gene_min_mean_for_ranking = 0.05,
  nmf_feature_selection_ratio = 0.20,
  nmf_nrun = 1,
  nmf_feature_number = 5,
  nmf_k_value = 50,
  downsample_gemli_ratio = 1,
  specc_kernel = "laplacedot",
  nmf_eval_downsample_cells = NULL
)
```

对于小鼠数据，将 `species` 改为 `"mouse"`。

示例中的 `parallel_cores = 25` 必须根据作业调度器实际分配的 CPU 和可用内存调整。主扫描不会自动把显式请求值限制到可用核心数。

### 5.2 完整参数

| 参数 | 默认值 | 说明 |
| --- | ---: | --- |
| `assay` | `"RNA"` | 主扫描读取的 assay |
| `input_mode` | `"raw_counts"` | 原始 counts 或预归一化输入 |
| `pre_normalized_layer` | `"data"` | 预归一化模式读取的 layer/slot |
| `check_raw_counts` | `TRUE` | 对 raw-count 输入执行诊断性检查 |
| `results_dir` | `"iscoral_results_raw"` | 运行日志目录 |
| `species` | `"human"` | 决定线粒体和 cell-cycle genes 的处理 |
| `random_seed` | `20251009` | 控制后台重复扫描的 RNG |
| `bin_size` | `200` | 目标每个临时 bin 的细胞数，不是 clone 数 |
| `n_repeats` | `100` | 共识扫描重复次数 |
| `parallel_cores` | `25` | 后台请求的 CPU workers |
| `downsampleUMI` | `NULL` | raw-count 模式中的可选 UMI 深度下采样 |
| `max_cells` | `NULL` | 主扫描前的可选全局细胞下采样 |
| `gemli_mean_q1` | `0.98` | 第一条候选基因规则的 mean quantile |
| `gemli_mean_q2` | `0.90` | 第二条候选基因规则的 mean quantile |
| `gemli_fano_q2` | `0.40` | 第二条候选基因规则的 Fano quantile |
| `gemli_mean_q3` | `0.60` | 第三条候选基因规则的 mean quantile |
| `gemli_fano_q3` | `0.90` | 第三条候选基因规则的 Fano quantile |
| `genes_to_remove` | `NULL` | 除内置过滤外额外排除的基因 |
| `gene_min_mean_for_ranking` | `0.05` | component validation 排名使用的最低平均表达 |
| `nmf_feature_selection_ratio` | `0.2` | 名义保留的 NMF components 比例；实际至少保留两个 |
| `nmf_nrun` | `1` | 每次特征评估中的 NMF 运行次数 |
| `nmf_feature_number` | `5` | 每个 NMF component 选择的 top genes |
| `nmf_k_value` | `50` | NMF rank |
| `downsample_gemli_ratio` | `1` | 每次重复使用的候选基因比例 |
| `nmf_eval_downsample_cells` | `NULL` | NMF component validation 使用的可选细胞数 |
| `specc_kernel` | `"laplacedot"` | spectral clustering 的 kernlab kernel 名称 |
| `custom_log_file` | `NULL` | 可选的自定义进度日志路径 |

`bin_size` 可以是向量，例如 `c(50, 100, 200)`，相应结果保存为 `binsize_50`、`binsize_100` 和 `binsize_200`。内部目标 bin 数为 `max(2, floor(number_of_cells / bin_size))`。

### 5.3 四类下采样

- `downsampleUMI`：降低 raw counts 的 UMI 深度；
- `max_cells`：在整个分析前减少细胞数；
- `downsample_gemli_ratio`：每次重复只使用部分候选基因；
- `nmf_eval_downsample_cells`：仅在 NMF component validation 中使用部分细胞。

其中，`max_cells` 会对子集化返回的 Seurat 对象；raw-count 模式下的 `downsampleUMI` 会把下采样 counts 写回所选 assay，再运行 `NormalizeData()`。因此这两项不只是内部临时计算设置。

使用全部信息时，应在第 5.1 节的主函数调用中显式保持 `downsampleUMI = NULL`、`max_cells = NULL`、`downsample_gemli_ratio = 1` 和 `nmf_eval_downsample_cells = NULL`。

`check_raw_counts = TRUE` 只会发出诊断警告，不会自动修复输入。

## 6. 主扫描输出

在某个 bin size 至少得到两个有效重复的标准情况下，核心结果位于 `sample@misc$CORAL_scan_results`：

```text
CORAL_scan_results
├── run_parameters
└── binsize_<value>
    ├── ss_filter_average
    ├── S_stab
    ├── n_repeats_requested
    ├── n_repeats_valid
    └── n_repeats_failed
```

检查示例：

```r
names(sample@misc$CORAL_scan_results)

run_info <- sample@misc$CORAL_scan_results$run_parameters
run_info[c(
  "iscoral_version",
  "input_mode",
  "normalization_applied",
  "expression_layer_used"
)]

scan <- sample@misc$CORAL_scan_results[["binsize_200"]]

c(
  requested = scan$n_repeats_requested,
  valid = scan$n_repeats_valid,
  failed = scan$n_repeats_failed
)

memory_table <- scan$ss_filter_average
memory_table <- memory_table[
  order(memory_table$omega_square_infer, decreasing = TRUE),
]
head(memory_table)
```

`S_stab` 是尚未行归一化的对称共关联矩阵；`omega_square_infer` 是跨有效重复汇总并将负值截断为零的 memory score；`CORALembedding()` 会将 `S_stab` 行归一化为 `P_norm`。若某个 bin size 少于两个有效重复，函数会警告，且不会生成该 `binsize_<value>` key。

再次运行主函数会覆盖 `@misc$CORAL_scan_results`。`results_dir` 主要保存日志，不会自动保存返回的 Seurat 对象。

## 7. Memory genes 可视化：`PlotMemoryGenes()`

`PlotMemoryGenes()` 从指定扫描结果读取 `ss_filter_average`，按平均表达量过滤基因，并依据 `omega_square_infer` 降序绘图。默认使用严格条件 `mean > 0.1`，显示排名最高的 30 个基因。

```r
memory_gene_result <- PlotMemoryGenes(
  seurat_obj = sample,
  binsize_key = "binsize_200",
  min_mean = 0.1,
  top_n = 30,
  return_data = TRUE
)

memory_gene_table <- memory_gene_result$memory_genes
p_memory_genes <- memory_gene_result$plot

print(p_memory_genes)

memory_gene_plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(
  filename = "memory_genes_mean_gt_0.1.pdf",
  plot = p_memory_genes,
  width = memory_gene_plot_size[1],
  height = memory_gene_plot_size[2],
  units = "in",
  device = grDevices::cairo_pdf
)
```

`memory_gene_table` 保留全部通过过滤的基因及原有统计列；`top_n` 只限制图中显示的数量。设置 `top_n = NULL` 可绘制全部通过过滤的基因。该函数不会修改 `@misc$CORAL_scan_results`。

## 8. 确定性伪克隆与克隆表达可视化

```r
sample <- Assign_deterministic_pseudoclone(
  seurat_obj = sample,
  bin_size = 200,
  k_values = c(10, 20, 30),
  hclust_method = "ward.D2"
)

deterministic_group <- "pseudoclone_k10_bs200"
features_to_plot <- utils::head(memory_gene_table$name, 2)

clonal_expression_result <- simulate_clonal_expression(
  seurat_obj = sample,
  features = features_to_plot,
  group_by = deterministic_group,
  assay = "RNA",
  layer = "data",
  reduction = "umap",
  min_cells = 2,
  n_pseudoclones = 12,
  return_data = TRUE
)

p_clonal_expression <- clonal_expression_result$plot
pseudoclone_expression_table <- clonal_expression_result$pseudoclone_summary

print(p_clonal_expression)

clonal_expression_plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(
  filename = "iscoral_clonal_expression.pdf",
  plot = p_clonal_expression,
  width = clonal_expression_plot_size[1],
  height = clonal_expression_plot_size[2],
  units = "in",
  device = grDevices::cairo_pdf
)
```

`Assign_deterministic_pseudoclone()` 创建 `pseudoclone_k10_bs200` 等 metadata 列。`pseudoclone_` 前缀专门表示由 `S_stab` 经过层次聚类得到的确定性伪克隆分组，以避免与 embedding-based `FindMetaclone()` 产生的 Metaclone 概念混淆。该步骤没有随机抽样；在 `S_stab`、细胞顺序、`k_values` 和 `hclust_method` 相同的条件下，分组是确定性的。

`simulate_clonal_expression()` 是从 clone-reconstruction benchmark 中提取的 isCORAL-only 可视化 API。它在左侧绘制指定基因的原始 UMAP expression，在右侧绘制同一基因在确定性 pseudoclones 中的单细胞 violin distributions；不读取实验 clone ground truth、GEMLI 或 SNN 结果，也不生成合成表达矩阵。默认保留至少有两个细胞的最大 12 个 pseudoclones，并针对每个基因按平均表达从高到低排列。

`features` 必须显式提供，可以来自 `PlotMemoryGenes(..., return_data = TRUE)` 返回的 `memory_gene_table$name`，也可以是生物学上预先指定的基因。`assay` 和 `layer` 明确决定表达来源；默认 `layer = "data"`，可兼容只有预归一化 `data`、没有 counts 的对象。`return_data = TRUE` 返回作图对象、cell-level expression data、实际用于 violin 的数据和 pseudoclone summary。

`bin_size` 必须对应现有的 `binsize_<value>`。确定性分组函数要求规范缓存路径 `@misc$CORAL_scan_results`，且其中的 `S_stab` 必须与当前 Seurat 对象覆盖完全相同的细胞。

## 9. pCORAL embedding：`CORALembedding()`

```r
sample <- CORALembedding(
  seurat_obj = sample,
  binsize_key = "binsize_200",
  n_dims = 30,
  top_n_genes = 3000,
  max_cells_lasso = NULL,
  alpha_elasticnet = 0.05,
  workers_glmnet = 5,
  random_seed = 20251009
)
```

`binsize_key` 必须与主扫描输出对应。主扫描使用 `bin_size = 100` 时，这里必须使用 `"binsize_100"`。

重要行为：

- 固定读取已经归一化的 `RNA/data`；
- 不读取 counts，也不调用 `NormalizeData()`；
- `max_cells_lasso = NULL` 表示 Elastic Net 使用全部有效细胞；
- `max_cells_lasso` 不会减少前面的全细胞能量距离和 MDS 计算；
- 至少需要两个匹配的 memory genes；
- 实际 embedding 维数为 `min(n_dims, n_valid_cells - 5)`，因此至少需要 7 个有效细胞；
- 固定五折 Elastic Net 至少需要五个训练细胞；
- `workers_glmnet` 是请求值；实际 worker 数至少为 1、至多为 5，并依据可用核心数限制；核心检测失败或少于 2 核时使用 1 个 worker；
- 建议新分析显式设置 `random_seed`。

如果主扫描使用了非 `RNA` assay，仍需确保 embedding 所需的适当归一化表达已经存在于 `RNA/data`。这里不存在 assay 自动继承：使用者必须确认 `RNA/data` 与主扫描输入具有预期的一致生物学含义和基因命名，否则两个阶段会基于不同的表达空间。

输出包括：

```r
Seurat::Embeddings(sample, "CORAL_cell_embed")
Seurat::Embeddings(sample, "CORAL_clone_embed")
cache <- sample@misc$isCORAL_embedding
```

检查请求和实际参数：

```r
cache$run_parameters
cache$effective_parameters
cache$reproducibility
cache$preprocessing
cache$provenance
cache$outputs

cache$effective_parameters$effective_dimensions
cache$effective_parameters$actual_workers

c(
  valid_cells = cache$effective_parameters$n_valid_cells,
  training_cells = cache$effective_parameters$n_training_cells,
  subsampling_applied = cache$effective_parameters$subsampling_applied
)
```

再次运行 `CORALembedding()` 会覆盖 `CORAL_cell_embed`、`CORAL_clone_embed` 和 `@misc$isCORAL_embedding`；该 cache 不自动保留运行历史。

可选平滑：

```r
sample <- SmoothEmbedding(sample, lambda = 0.5)
```

`lambda = 0` 对应纯 cell embedding，`lambda = 1` 对应完整 pseudoclone diffusion。

## 10. Pseudoclone 可视化

### 10.1 单个 seed cell

```r
seed_cells_to_plot <- colnames(sample)[1:2]

p_pseudoclone <- PlotPseudoclone(
  seurat_obj = sample,
  seed_cell = seed_cells_to_plot[1],
  binsize_key = "binsize_200",
  reduction = "CORAL_clone_embed",
  show_ground_truth = FALSE
)

print(p_pseudoclone)

pseudoclone_size <- grDevices::dev.size("in")
ggplot2::ggsave(
  filename = "pseudoclone_seed_cell.pdf",
  plot = p_pseudoclone,
  width = pseudoclone_size[1],
  height = pseudoclone_size[2],
  units = "in",
  device = grDevices::cairo_pdf
)
```

### 10.2 两个 seed cells

```r
p_compare <- PlotPseudoclone(
  seurat_obj = sample,
  seed_cell = seed_cells_to_plot,
  binsize_key = "binsize_200",
  reduction = "CORAL_clone_embed",
  show_ground_truth = FALSE
)

print(p_compare)

compare_size <- grDevices::dev.size("in")
ggplot2::ggsave(
  filename = "pseudoclone_two_seeds.pdf",
  plot = p_compare,
  width = compare_size[1],
  height = compare_size[2],
  units = "in",
  device = grDevices::cairo_pdf
)
```

如有实验 clone label，可提供 `clone_col` 和一个或两个 `clone_ids`。若扫描结果包含多个 bin-size key，必须显式设置 `binsize_key`。在没有显式提供 seed cell 时，函数会随机选择一个细胞；此时建议设置 `random_seed`。设置 `return_data = TRUE` 可同时返回绘图数据、实际 seed cells、reduction 和解析后的 bin-size key。

## 11. Metaclone：`FindMetaclone()`

建议根据实际 reduction 维数构建 `dims`：

```r
metaclone_dims <- seq_len(
  ncol(Seurat::Embeddings(sample, "CORAL_clone_embed"))
)

sample <- FindMetaclone(
  seurat_obj = sample,
  resolution = 0.7,
  reduction = "CORAL_clone_embed",
  dims = metaclone_dims,
  graph_prefix = "pcoral",
  cluster_col_name = "metaclone",
  algorithm = 4,
  random_seed = 42,
  verbose = TRUE
)
```

默认输出包括 `sample$metaclone`、`sample@graphs$pcoral_nn`、`sample@graphs$pcoral_snn` 和 `sample@misc$isCORAL_metaclone`，并将 `metaclone` 设为 active identity。

检查实际聚类后端：

```r
sample@misc$isCORAL_metaclone$latest$effective_parameters$clustering_engine
```

`algorithm = 4` 请求 Leiden，所支持的可选后端为 `leidenbase`，该后端由 Seurat 5.2.0 起采用。函数只检查该 namespace 与 Seurat 版本是否满足要求，不会自动安装依赖；若无法使用该后端，函数会明确警告并将有效算法改为 Louvain（`algorithm = 1`）。不同后端的结果不能仅因为 seed 相同就视为数值等价。

`resolution` 控制社区划分粒度，不存在适用于所有数据集的唯一正确值。使用相同的 `graph_prefix` 或 `cluster_col_name` 重新运行会覆盖同名 graphs 或 metadata，而 `@misc$isCORAL_metaclone$history` 会追加运行记录。

### 11.1 Metaclone landscape

```r
p_metaclone <- PlotMetacloneLandscape(
  seurat_obj = sample,
  reduction = "CORAL_clone_embed",
  cluster_col = "metaclone"
)

print(p_metaclone)

metaclone_size <- grDevices::dev.size("in")
ggplot2::ggsave(
  filename = "metaclone_landscape.pdf",
  plot = p_metaclone,
  width = metaclone_size[1],
  height = metaclone_size[2],
  units = "in",
  device = grDevices::cairo_pdf
)
```

## 12. 有实验 lineage labels 时的验证

如果需要另外计算实验 lineage labels 对应的 gene-wise omega-squared 及 permutation background，可选运行 `calculate_heritable_genes()`：

```r
Seurat::DefaultAssay(sample) <- "RNA"

sample <- calculate_heritable_genes(
  seurat_obj = sample,
  group_by = "barcode",
  invalid_barcodes = c("0", "NA", "Unknown", "Unassigned", ""),
  min_cells_per_group = 2,
  permutation_repeats = 100,
  n_cores = 20
)
```

示例中的 `n_cores = 20` 必须按实际分配的 CPU 和内存调整。

该函数固定读取当前 `DefaultAssay` 的 `data`，并将结果写入 `@misc$groundtruth`：

```r
sample@misc$groundtruth$heritable_genes_df
sample@misc$groundtruth$heritable_genes_threshold
sample@misc$groundtruth$omega_squared_null_distribution
```

完成 isCORAL 核心扫描后，可直接调用 `CompareMemoryScores()`，将实验 lineage barcode 对应的 ground-truth omega-squared 与 isCORAL memory score 进行比较：

```r
memory_comparison <- CompareMemoryScores(
  seurat_obj = sample,
  binsize_key = "binsize_200",
  assay = "RNA",
  min_mean = 0.1,
  top_n = 2000
)

memory_comparison$summary
head(memory_comparison$comparison)
```

该 API 使用 log-normalized `RNA/data` 自动完成基因对齐和平均表达过滤，将实验 omega-squared 排名前 2,000 的基因定义为阳性集，并计算 AUROC、AUPRC、AUPRC enrichment、Max F1、Spearman/Pearson 相关性及 RMSE。它同时返回逐基因比较表和 ROC/PR 曲线，且不会修改输入的 Seurat 对象。

函数不会缓存过滤后的表达矩阵或标签。每次调用都会重新读取当前 `DefaultAssay` 的 log-normalized `data` layer/slot，并按细胞名重新解析本次 `group_by` 与过滤设置；旧对象中遗留的 `@misc$groundtruth$internal_data` 会被自动移除。轻量输入记录保存在 `@misc$groundtruth$input_provenance`，因此可以核对本次使用的 assay、标签列、过滤参数、细胞数和并行参数。当前函数没有 `random_seed` 参数，也未使用 `doRNG`，因此 permutation null distribution 的严格复现还取决于外部 RNG 状态和并行后端。该函数还可能改变当前注册的 foreach 后端，且当前实现不会恢复调用前的后端；后续并行步骤不应假定原后端仍然有效。

这一步不是 `PlotValidationSuite()` 的前置条件。验证图直接读取实验标签 metadata、Metaclone 标签、指定 reduction 和 `@misc$isCORAL_embedding$P_norm`：

```r
p_validation <- PlotValidationSuite(
  seurat_obj = sample,
  reduction = "CORAL_clone_embed",
  space_title = "pCORAL clone space",
  barcode_col = "experimental_clone",
  cluster_col = "metaclone",
  na_barcodes = c("Unassigned")
)

print(p_validation)

validation_size <- grDevices::dev.size("in")
ggplot2::ggsave(
  filename = "pseudoclone_validation.pdf",
  plot = p_validation,
  width = validation_size[1],
  height = validation_size[2],
  units = "in",
  device = grDevices::cairo_pdf
)
```

`PlotValidationSuite()` 需要 `ggrastr`，并假设 Metaclone 标签可转换为数字；标准 `FindMetaclone()` 产生的标签符合这一要求。字符标签 `"0"` 固定按无效实验 barcode 处理，`na_barcodes` 只能在内置无效标签集合上追加，不能移除其中的值。`PlotValidationSuite()` 和 `PlotInferenceSuite()` 建议仅用于未 subset、未重排的全细胞标准流程。使用前应确认 Seurat 对象、reduction 和 `P_norm` 覆盖同一批且顺序一致的细胞。一般 seed-level 展示优先使用会按条形码对齐的 `PlotPseudoclone()`。

`PlotMetacloneSeedHeatmap()`、`PlotMetacloneInferredHeatmap()` 和 `PlotMetacloneGroundTruthHeatmap()` 返回列表，其中包含 `$table` 和 `$plot`。这些函数要求至少两个可聚类的 Metaclones，并要求相应的 cell-type 或 barcode metadata 已存在。

## 13. 随机性与可复现性

主扫描、embedding 和 Metaclone 三个核心函数分别暴露独立的 `random_seed` 参数。建议分别在相应的完整函数调用中显式设置：

| 阶段 | 参数 | 本文示例值 |
| --- | --- | ---: |
| 主扫描 | `isCORAL(random_seed = ...)` | `20251009` |
| pCORAL embedding | `CORALembedding(random_seed = ...)` | `20251009` |
| Metaclone | `FindMetaclone(random_seed = ...)` | `42` |

注意：

1. 主扫描的 `random_seed` 在 callr/doRNG 后台阶段控制重复扫描；运行 UID、`max_cells` 抽样和 raw-count 模式的 `SampleUMI()` 在此之前发生。启用这些下采样时，应在调用前另行 `set.seed()`。
2. `CORALembedding(random_seed = NULL)` 不会设置种子；IRLBA、可选训练细胞抽样和五折划分仍然具有随机性。
3. `CORALembedding()` 会推进全局 RNG；即使显式设置 `random_seed`，退出时也不会恢复调用前的 RNG 状态。
4. `max_cells_lasso = NULL` 只表示全部有效细胞参与 Elastic Net 训练，不表示整个 embedding 没有随机性。
5. `PlotPseudoclone()` 在未提供 seed cell 时会随机选取细胞，应设置它自己的 `random_seed`；显式提供 seed 后不存在这一步随机选择。
6. `calculate_heritable_genes()` 的标签置换也涉及随机性，但当前没有 seed 参数或 doRNG 记录；其严格复现需额外固定并记录外部 RNG 和并行后端。
7. 完整复现还依赖输入矩阵、细胞和基因顺序、R/包/BLAS版本及聚类后端。

建议保存对象、session 信息和 `results_dir/logs`：

```r
saveRDS(sample, "sample_isCORAL_beta_1.0.0.rds")

writeLines(
  capture.output(sessionInfo()),
  "sample_isCORAL_beta_1.0.0_sessionInfo.txt"
)
```

## 14. 兼容入口

新脚本应使用 `isCORAL()`、`PlotMemoryGenes()`、`Assign_deterministic_pseudoclone()`、`simulate_clonal_expression()`、`CORALembedding()` 和 `FindMetaclone()`。beta 1.0.0 不再提供任何含内部版本号的主扫描入口。

以下非主扫描名称仅用于旧脚本调用兼容：

- `FindMetaclones()`：弃用包装器，转发到 `FindMetaclone()`；
- `SimulateClone()`、`Simulate_clonal_expression_simple()` 和 `SimulateMetaClone()`：仅为旧脚本保留的弃用别名，调用时发出警告并转发到 `Assign_deterministic_pseudoclone()`。
- `PlotPseudoclone()` 的 `seed_cells` 和 `seed` 参数：仅为旧调用保留的弃用参数；新代码分别使用 `seed_cell` 和 `random_seed`。

这些入口不会模拟归档版本的历史行为。严格复现历史分析时，应使用对应归档包，并匹配输入、参数、RNG、依赖版本和聚类后端。

当前 beta 1.0.0 的导出接口中没有 `isCORAL_scale()`；不要依据旧文档将其作为本版本的可用函数。

## 15. 常见问题

| 问题 | 原因与处理 |
| --- | --- |
| 找不到 `binsize_200` | 检查 `names(sample@misc$CORAL_scan_results)` 并设置正确的 `binsize_key` |
| 找不到 `RNA/data` | `CORALembedding()` 不执行归一化；先确保 `RNA/data` 存在适当的 log-normalized expression |
| counts 缺失 | 主扫描显式使用 `input_mode = "pre_normalized_data"` |
| 至少需要两个 memory genes | 检查基因名与 `ss_filter_average` 是否一致 |
| 五折 CV 训练细胞不足 | 将 `max_cells_lasso` 设为 `NULL`、`Inf` 或至少 5 |
| probability matrix barcode 报错 | `S_stab` 两轴细胞集合不一致；不要只删除或重排一侧 |
| 找不到 `CORAL_clone_embed` | 先运行 `CORALembedding()` |
| `dims` 超出范围 | 使用实际 reduction 维数构建 `dims` |
| Leiden 自动变成 Louvain | 运行 `install.packages("leidenbase")` 后重新请求 Leiden，或显式使用 `algorithm = 1` 并记录后端 |
| `metaclone` 没有保存 | 必须执行 `sample <- FindMetaclone(sample, ...)` |
| `PlotPseudoclone()` 报多个结果 key | 显式指定 `binsize_key` |
| validation suite 长度或下标错误 | 核对对象、reduction 与 `P_norm` 的细胞集合和顺序 |

## 16. 使用边界

- pseudoclone 是转录组推断结果，不应在没有实验验证时直接称为真实 clone；
- Metaclone 是 clone embedding 中的图社区，不是实验 lineage identity；
- beta 1.0.0 保留 6.1.5 基线的核心分析定义、输入校验和运行参数记录，但移除了含内部版本号的主扫描入口；
- 不应仅凭静态源码审查宣称不同版本逐位一致、所有 Seurat v4/v5 组合均已全面测试，或所有数据规模均可无条件处理；
- `@misc` 中的 provenance 是扩展运行记录，但不包含原始输入副本或输入 checksum。

更多细节见同目录中的：

- `BETA_1.0.0_COMPATIBILITY.md`
