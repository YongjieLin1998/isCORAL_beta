# isCORAL 6.1.0–6.1.5：每个版本的差异和更新
# isCORAL 6.1.0–6.1.5: Version-by-Version Differences and Updates

> **归档说明 / Archive notice:** 本文记录 beta 1.0.0 之前的内部开发版本，
> 其中的版本化函数名仅用于解释历史源码，不代表 beta 1.0.0 的公开 API。
> beta 1.0.0 仅公开未带版本号的主扫描入口 `isCORAL()`。
> This document records internal development versions that predate beta 1.0.0.
> Version-numbered function names are retained here only to explain archived
> source code and are not part of the beta 1.0.0 public API. Beta 1.0.0 exposes
> only the unversioned core entry point `isCORAL()`.
>
> **当前 beta 下游接口 / Current beta downstream API:** beta 1.0.0 使用
> `Assign_deterministic_pseudoclone()` 作为确定性层次分组的正式入口，并使用
> `simulate_clonal_expression()` 绘制指定基因的 isCORAL-only clonal-expression
> 图。本文随后出现的 `Simulate_clonal_expression_simple()` 仍用于准确记录归档
> 6.1.5 的历史接口；在当前 beta 中，`SimulateClone()`、
> `Simulate_clonal_expression_simple()` 和 `SimulateMetaClone()` 均为弃用别名。
> Beta 1.0.0 uses `Assign_deterministic_pseudoclone()` as the canonical
> deterministic hierarchical-assignment entry point and
> `simulate_clonal_expression()` for selected-gene, isCORAL-only
> clonal-expression visualization. Later references to
> `Simulate_clonal_expression_simple()` accurately describe the archived 6.1.5
> API. In the current beta, `SimulateClone()`,
> `Simulate_clonal_expression_simple()`, and `SimulateMetaClone()` are all
> deprecated aliases.

本文面向使用 isCORAL 开展单细胞转录组分析的使用者、生物学家和方法开发者，
完整说明当前归档中 6.1.0、6.1.1、6.1.2、6.1.3、6.1.4 和 6.1.5 六个版本
的功能、默认行为、兼容性边界和迁移方式。

This document is written for isCORAL users, biologists, and method developers.
It provides a self-contained account of the functions, default behavior,
compatibility boundaries, and migration paths of the six archived releases:
6.1.0, 6.1.1, 6.1.2, 6.1.3, 6.1.4, and 6.1.5.

本文件以各版本现存包源码、`DESCRIPTION`、`NAMESPACE`、帮助页和 Seurat
缓存结构的直接比较为依据，不依赖 Git 提交记录。除非明确写为“已运行验证”，
文中的一致性结论均属于源码级确认，而不是对所有 R、Seurat、glmnet、IRLBA、
BLAS 或 Leiden 组合的运行保证。

The account is based on direct comparison of the package sources,
`DESCRIPTION`, `NAMESPACE`, help pages, and Seurat cache schemas present in
each archived release; it does not depend on Git commit history. Unless a
result is explicitly described as runtime-validated, compatibility statements
refer to source-level findings rather than a guarantee across every
combination of R, Seurat, glmnet, IRLBA, BLAS, or Leiden.

---

## 1. 科学定位与术语 / Scientific scope and terminology

### 中文

isCORAL 的核心目标，是仅从单细胞转录组中提取可遗传或具有谱系记忆意义的
表达结构，并推断细胞之间的 pseudoclone（伪克隆）关联。该概率性推断可作为
**in-silico lineage tracing**：它不要求预先测得 DNA barcode、CellTag 或
体细胞突变，而是从转录组的结构相关性中重建潜在谱系关系。

论文及其 Figure 3 已经在明确的数据集、指标和比较方法范围内完成 clone
reconstruction benchmark，并报告了领先结果。因此，在论文所覆盖的
benchmark 范围内，
pseudoclone inference 可以被描述为仅基于转录组进行克隆重建的
state-of-the-art（SOTA）模型。该表述应始终限定于论文报告的实验范围；软件
版本比较本身不是一次新的独立 benchmark，也不能据此宣称每个后续版本都会
进一步提高生物学准确率。

本文区分以下概念：

- **isCORAL core scan**：识别 memory/GEMLI genes，并生成细胞间共识结构
  `S_stab`。
- **pseudoclone inference**：由 `S_stab` 得到每个细胞的概率性克隆关联
  分布，而不是实验测得的真实 clone label。
- **pCORAL embedding**：`CORALembedding()` 使用平方能量距离、经典 MDS
  和多响应 Elastic Net，建立 `CORAL_cell_embed`；再通过概率期望
  `P_norm %*% Y_ref_raw` 建立 `CORAL_clone_embed`。
- **Metaclone**：在 `CORAL_clone_embed` 上建立 SNN 图并执行 Leiden 或
  Louvain 聚类后得到的离散社区。
- **Simulate_clonal_expression_simple**：对 `1 - S_stab` 进行层次聚类并按指定 `k` 切树。
  它不是 embedding-based Metaclone，也不是实验 clone 的生成模型。

### English

The central purpose of isCORAL is to extract heritable or lineage-memory-like
expression structure from single-cell transcriptomes and infer probabilistic
pseudoclone relationships among cells. This inference can serve as
**in-silico lineage tracing**: it does not require an experimentally measured
DNA barcode, CellTag, or somatic mutation, but instead reconstructs latent
lineage relationships from transcriptomic structural correlation.

The manuscript and its Figure 3 include an explicit clone-reconstruction
benchmark and report leading results within the stated datasets, metrics, and
comparator set.
Within that benchmarked scope, pseudoclone inference may therefore be
described as a state-of-the-art model for transcriptome-only clone
reconstruction. The claim should remain bounded by the experiments reported
in the manuscript. A software-version comparison is not a new independent
benchmark and does not imply that every later version further improves
biological accuracy.

This document distinguishes the following concepts:

- **isCORAL core scan**: identifies memory/GEMLI genes and produces the
  cell-cell consensus structure `S_stab`.
- **Pseudoclone inference**: derives a probabilistic clone-association
  distribution for each cell from `S_stab`; this is not an experimentally
  observed ground-truth clone label.
- **pCORAL embedding**: `CORALembedding()` combines squared energy distance,
  classical MDS, and multi-response Elastic Net to construct
  `CORAL_cell_embed`, then obtains `CORAL_clone_embed` through the probability
  expectation `P_norm %*% Y_ref_raw`.
- **Metaclone**: a discrete community detected by Leiden or Louvain after
  building an SNN graph in `CORAL_clone_embed`.
- **Simulate_clonal_expression_simple**: hierarchical clustering of `1 - S_stab`, followed by
  tree cutting at specified values of `k`. It is neither an embedding-based
  Metaclone nor a generative model of experimentally observed clones.

---

## 2. 一页式版本总览 / One-page version overview

| 版本 / Version | 主要定位 / Main role | 主扫描默认输入 / Default scan input | 关键更新或边界 / Principal update or boundary |
| --- | --- | --- | --- |
| 6.1.0 | 该系列的早期完整工作流 / Early complete workflow in this series | 无独立输入模式；总是调用 `NormalizeData()` / No separate input mode; always calls `NormalizeData()` | 默认 mouse、`bin_size = 50`；导出复数 `FindMetaclones()`；从GEMLI/NMF候选集合硬编码排除大写 `MALAT1` / Defaults to mouse and `bin_size = 50`; exports plural `FindMetaclones()`; hard-coded exclusion of uppercase `MALAT1` from the GEMLI/NMF candidate set |
| 6.1.1 | 小型预处理修订 / Focused preprocessing revision | 与6.1.0相同 / Same as 6.1.0 | 仅取消GEMLI/NMF候选集合中的内置 `MALAT1` 排除；其他算法、默认值和公开函数与6.1.0相同 / Only removes the built-in `MALAT1` exclusion from the GEMLI/NMF candidate set; other algorithms, defaults, and exports are unchanged |
| 6.1.2 | 更新科研默认值的过渡版本 / Transitional release with updated scientific defaults | 仍总是归一化 / Still always normalizes | 默认改为 human、scan-worker `random_seed = 20251009`、`bin_size = 200`；存在 `downsampleUMI` 源码错误和公开导出不一致 / Defaults switch to human, scan-worker `random_seed = 20251009`, and `bin_size = 200`; contains a `downsampleUMI` source defect and export inconsistencies |
| 6.1.3 | 预归一化数据优先的版本 / Pre-normalized-data-first release | `pre_normalized_data` | 增加输入模式、Seurat v4/v5版本引导读取和 `PlotPseudoclone()`；`downsampleUMI` 被保留但忽略；旧 `isCORAL_metaclone.R`、embedding-based Find函数及validation suite缺失，`SimulateClone()`仍为未导出源码 / Adds input modes, version-guided Seurat v4/v5 reads, and `PlotPseudoclone()`; retains but ignores `downsampleUMI`; the old `isCORAL_metaclone.R`, embedding-based Find functions, and validation suite are absent, while `SimulateClone()` remains source-only |
| 6.1.4 | 功能重新整合与大数据版本 / Reintegrated full-feature and large-data release | `raw_counts` | 完成Seurat v4/v5 layer/slot读写；恢复可视化；正式使用单数 `FindMetaclone()`；增加 `isCORAL_scale()` / Completes Seurat v4/v5 layer/slot read/write handling; restores visualization; canonical singular `FindMetaclone()`; adds `isCORAL_scale()` |
| 6.1.5 | 兼容性、审计和可复现性版本 / Compatibility, auditability, and reproducibility release | `raw_counts` | 核心扫描定义继承6.1.4；增强 embedding 与 Metaclone 记录；将简单层次聚类正式命名为 `Simulate_clonal_expression_simple()` 并保留旧入口；恢复弃用的复数兼容入口；清理 namespace 冲突 / Retains the 6.1.4 core scan; expands embedding and Metaclone records; names the simple hierarchical method `Simulate_clonal_expression_simple()` while retaining old entry points; restores a deprecated plural wrapper; removes namespace masking |

最重要的用户结论是：6.1.0–6.1.5 不是一条“版本号越高、算法精度必然越高”
的序列。版本间最大的实际差异来自默认输入预处理、默认参数、公开 API 的完整
程度、异常输入处理和可复现性记录。pCORAL 的主要 embedding 数学路径在
6.1.0–6.1.5 中没有被有意改写。

The most important practical conclusion is that versions 6.1.0–6.1.5 do not
form a sequence in which a higher version number necessarily means greater
analytical accuracy. The largest differences concern default preprocessing,
default parameters, completeness of the public API, handling of invalid
inputs, and reproducibility records. The principal pCORAL embedding
mathematics was not intentionally rewritten across versions 6.1.0–6.1.5.

---

## 3. 主扫描默认值的版本变化 / Changes in core-scan defaults

| 版本 / Version | 可公开调用的正式入口 / Public canonical entry | `input_mode` 的实际默认值 / Effective default | `species` | `random_seed` | `bin_size` | `parallel_cores` | `downsample_gemli_ratio` | `nmf_eval_downsample_cells` |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| 6.1.0 | 无；源码定义未导出的 `isCORAL_v6_1_0()` / None; source defines unexported `isCORAL_v6_1_0()` | 无；总是 `NormalizeData()` / None; always normalizes | mouse | `NULL` | 50 | 30 | 0.75 | 2000 |
| 6.1.1 | 无；仍为未导出的 `isCORAL_v6_1_0()` / None; still unexported `isCORAL_v6_1_0()` | 与6.1.0相同 / Same as 6.1.0 | mouse | `NULL` | 50 | 30 | 0.75 | 2000 |
| 6.1.2 | 无；仍为未导出的 `isCORAL_v6_1_0()` / None; still unexported `isCORAL_v6_1_0()` | 无；总是 `NormalizeData()` / None; always normalizes | human | 20251009 | 200 | 25 | 1 | `NULL` |
| 6.1.3 | `isCORAL_v6_1_3()` | `pre_normalized_data` | human | 20251009 | 200 | 25 | 1 | `NULL` |
| 6.1.4 | `isCORAL_v6_1_4()` | `raw_counts` | human | 20251009 | 200 | 25 | 1 | `NULL` |
| 6.1.5 | `isCORAL_v6_1_5()` | `raw_counts` | human | 20251009 | 200 | 25 | 1 | `NULL` |

表中未列出的 GEMLI、NMF 和谱聚类默认参数在这六个归档版本中保持一致。
尤其是 `n_repeats = 100`、`gemli_mean_q1 = 0.98`、
`gemli_mean_q2 = 0.90`、`gemli_fano_q2 = 0.40`、
`gemli_mean_q3 = 0.60`、`gemli_fano_q3 = 0.90`、
`gene_min_mean_for_ranking = 0.05`、`nmf_feature_selection_ratio = 0.2`、
`nmf_nrun = 1`、`nmf_feature_number = 5`、`nmf_k_value = 50` 和
`specc_kernel = "laplacedot"` 没有改变。

The GEMLI, NMF, and spectral-clustering defaults not shown in the table remain
the same in the six archived releases. In particular, `n_repeats = 100`,
`gemli_mean_q1 = 0.98`, `gemli_mean_q2 = 0.90`,
`gemli_fano_q2 = 0.40`, `gemli_mean_q3 = 0.60`,
`gemli_fano_q3 = 0.90`, `gene_min_mean_for_ranking = 0.05`,
`nmf_feature_selection_ratio = 0.2`, `nmf_nrun = 1`,
`nmf_feature_number = 5`, `nmf_k_value = 50`, and
`specc_kernel = "laplacedot"` do not change.

---

## 4. isCORAL 6.1.0

### 中文

6.1.0提供了该系列较早的完整分析组合：核心结构相关性扫描、层次聚类版
`SimulateClone()`、pCORAL embedding、embedding-based
`FindMetaclones()`、ground-truth 辅助计算和多种 validation 图形。

需要注意以下行为：

- 主扫描源码函数名为 `isCORAL_v6_1_0()`，但没有被 `NAMESPACE` 导出，
  也没有对应帮助页；它不是一个完整包装的公开入口。
- 主扫描默认 `species = "mouse"`、`random_seed = NULL`、
  `bin_size = 50`。
- 没有 `input_mode`。函数总是调用 `NormalizeData()`，随后读取归一化
  `data`；因此它没有受支持的“仅使用已有 data、跳过归一化”路径。
- `downsampleUMI` 非空且为 numeric 时，函数从形参 `assay` 指定的 counts
  layer读取计数；counts缺失或为空时停止。结果写回该 assay 后，函数调用未
  显式传入 `assay` 的 `NormalizeData()`。因此，只有
  `DefaultAssay(seurat_obj) == assay` 时，才能确认刚下采样的 counts 被
  同一次调用重新归一化。该旧实现也没有版本引导的 v4 `slot` 回退或 setter。
- 大写 `"MALAT1"` 被硬编码排除出 GEMLI/NMF候选基因集合；它并没有从表达
  矩阵、全部表达统计或后续 embedding特征中被全局删除。
- `CORALembedding()` 的默认 `binsize_key = "binsize_200"`，但核心扫描
  默认产生 `binsize_50`。按默认核心参数运行后，embedding 时应显式指定
  `binsize_key = "binsize_50"`。
- 公开的 embedding 社区入口是复数 `FindMetaclones()`。它会把标签写入
  metadata 并设为 active identity，但调用者必须接收返回的 Seurat 对象。
- `CORALembedding()` 默认
  `n_dims = 30`、`top_n_genes = 3000`、
  `max_cells_lasso = NULL`、`alpha_elasticnet = 0.05`、
  `workers_glmnet = 5`。`max_cells_lasso = NULL` 表示 Elastic Net
  使用全部有效细胞，而不是自动下采样。

`species` 只控制线粒体及细胞周期基因的命名匹配，不执行物种转换。标准小鼠
符号通常是 `Malat1`，因此6.1.0硬编码的大写 `"MALAT1"` 对标准小鼠数据
通常本来就不生效；它主要影响使用标准大写人类符号的数据。即使在6.1.0中，
`ss_filter_average` 仍可包含 `MALAT1`，而 `CORALembedding()` 也不会重新
应用该排除集。

### English

Version 6.1.0 provides an early complete combination of the structural-
correlation scan, hierarchical `SimulateClone()`, pCORAL embedding,
embedding-based `FindMetaclones()`, ground-truth helpers, and multiple
validation plots.

Important behavior includes:

- The main scan is defined in source as `isCORAL_v6_1_0()`, but it is not
  exported in `NAMESPACE` and has no help page. It is therefore not a fully
  packaged public entry point.
- The scan defaults to `species = "mouse"`, `random_seed = NULL`, and
  `bin_size = 50`.
- There is no `input_mode`. The function always calls `NormalizeData()` and
  then reads normalized `data`; it has no supported path that uses existing
  data while explicitly skipping normalization.
- If `downsampleUMI` is non-`NULL` and numeric, counts are read from the
  layer selected by the formal `assay` argument, and the function stops if
  that layer is missing or empty. The result is written back to that assay,
  after which `NormalizeData()` is called without an explicit `assay`.
  The newly downsampled counts are therefore guaranteed to be renormalized by
  that call only when `DefaultAssay(seurat_obj) == assay`. This old
  implementation also has no version-guided v4 `slot` fallback or setter.
- Uppercase `"MALAT1"` is hard-coded out of the GEMLI/NMF candidate-gene
  set; it is not globally removed from the expression matrix, all expression
  statistics, or later embedding features.
- `CORALembedding()` defaults to `binsize_key = "binsize_200"`, whereas the
  default scan produces `binsize_50`. After a default scan, embedding should
  therefore be called with `binsize_key = "binsize_50"`.
- The public embedding-community entry point is plural
  `FindMetaclones()`. It writes labels to metadata and sets the active
  identity, but the caller must assign the returned Seurat object.
- `CORALembedding()` defaults to `n_dims = 30`,
  `top_n_genes = 3000`, `max_cells_lasso = NULL`,
  `alpha_elasticnet = 0.05`, and `workers_glmnet = 5`.
  `max_cells_lasso = NULL` means that all valid cells, rather than an
  automatically subsampled set, are used for Elastic Net training.

`species` only controls gene-symbol matching for mitochondrial and cell-cycle
genes; it does not convert species. The canonical mouse symbol is generally
`Malat1`, so the hard-coded uppercase `"MALAT1"` in version 6.1.0 would
normally not affect standard mouse data. Its main effect is on data using
standard uppercase human symbols. Even in version 6.1.0,
`ss_filter_average` can still contain `MALAT1`, and `CORALembedding()` does
not reapply this exclusion set.

---

## 5. isCORAL 6.1.1

### 中文

6.1.1是一个范围很小但可能具有生物学影响的预处理更新。除版本号和一处基因
排除规则外，6.1.0与6.1.1的包文件基本相同：核心算法、公开函数、全部默认
参数、pCORAL embedding、Metaclone 聚类和输出结构均没有改变。

唯一的实质变化是取消6.1.0中从 GEMLI/NMF候选集合排除大写 `"MALAT1"` 的
规则：

- 对标准人类符号，`MALAT1` 在6.1.1中可以重新进入 GEMLI/NMF候选集合；
  只有当它通过表达量/Fano等筛选时，才可能直接参与候选特征和NMF，并由此
  间接影响 `S_stab`、embedding 和 Metaclone。
- 对标准小鼠符号 `Malat1`，6.1.0的原规则通常本来就匹配不到，因此仅就
  这一规则而言，6.1.0和6.1.1通常等效。
- 如果希望在6.1.1或后续版本中明确排除它，应使用与表达矩阵完全同名且大小写
  一致的 `genes_to_remove`。

该排除规则并非全局删除：6.1.0的 `ss_filter_average` 仍然可以保存
`MALAT1` 的统计量，`CORALembedding()` 也不重新应用
`genes_to_remove`。因此，两个版本的直接差异是它能否进入 GEMLI/NMF候选
集合，而不是它能否在任何下游对象中出现。

6.1.1的主函数在源码中仍名为 `isCORAL_v6_1_0()`，仍未导出，控制台信息
仍显示6.1.0。因此，仅根据函数名、控制台输出或早期 `@misc` 记录，无法可靠
区分一次扫描究竟来自6.1.0还是6.1.1。

### English

Version 6.1.1 is a narrowly scoped preprocessing update with potential
biological consequences. Apart from the version field and one gene-exclusion
rule, versions 6.1.0 and 6.1.1 are effectively the same package: the core
algorithms, exported functions, all defaults, pCORAL embedding, Metaclone
clustering, and output schemas are unchanged.

The only substantive change is removal of the version 6.1.0 rule that
excluded uppercase `"MALAT1"` from the GEMLI/NMF candidate set:

- With standard human symbols, `MALAT1` can re-enter the candidate set in
  version 6.1.1. Only if it passes the expression/Fano filters can it
  directly participate in candidate-feature selection and NMF, with possible
  indirect effects on `S_stab`, embedding, and Metaclone.
- With the standard mouse symbol `Malat1`, the original 6.1.0 rule would
  normally fail to match, so the two versions are generally equivalent with
  respect to this rule alone.
- To exclude the gene explicitly in version 6.1.1 or later, use
  `genes_to_remove` with exactly the same spelling and case as the expression
  matrix.

This is not a global gene removal. Version 6.1.0 can still store `MALAT1`
statistics in `ss_filter_average`, and `CORALembedding()` does not reapply
`genes_to_remove`. The direct difference between the two releases is
eligibility for the GEMLI/NMF candidate set, not whether the gene can appear
anywhere in downstream objects.

The source-level main function in version 6.1.1 is still named
`isCORAL_v6_1_0()`, remains unexported, and still prints a 6.1.0 console
message. A scan from version 6.1.0 therefore cannot be reliably distinguished
from one from version 6.1.1 using only that function name, console text, or
the early `@misc` record.

---

## 6. isCORAL 6.1.2

### 中文

6.1.2把主扫描的科研默认值调整为后来版本沿用的方向：

- `species` 从 mouse 改为 human；
- `random_seed` 从 `NULL` 改为 `20251009`；
- `bin_size` 从50改为200；
- `parallel_cores` 从30改为25；
- `downsample_gemli_ratio` 从0.75改为1；
- `nmf_eval_downsample_cells` 从2000改为 `NULL`。

它仍没有独立的 `input_mode`，仍会调用 `NormalizeData()`。主函数仍名为
`isCORAL_v6_1_0()` 且未导出。因此，6.1.2不是只通过更新包版本号就能直接
替代6.1.0/6.1.1的版本：即使显式参数相同，公开入口和默认预处理仍需核对。

6.1.2还存在两个重要的实现/包装边界：

1. `downsampleUMI` 分支中读取 counts 的源码表达式被路径文本污染。默认
   `downsampleUMI = NULL` 时该分支不执行；一旦传入 numeric 值，该版本
   预计会在此处失败。因此6.1.2不能作为可靠的 UMI 下采样版本。
2. `NAMESPACE` 只导出四个名称，其中导出了源码中不存在的
   `SimulateMetaClone`，却没有导出实际存在的 `SimulateClone()`、
   `FindMetaclones()`、`SmoothEmbedding()` 和 validation suite。
   Metaclone及绘图源码仍在包内，但普通使用者无法把它们视为受支持的公开 API。

该版本的 validation 源码增加了栅格化，并把字符串 `"0"` 加入默认无效
barcode；但是函数未导出、缺少完整帮助页，而且缺少 `ggrastr` 时会尝试在
运行过程中调用 `install.packages()`。这属于包封装问题，不应解释为新的
clone reconstruction 算法。

`CORALembedding()` 的计算源码与6.1.1保持一致。其真实源码默认仍为
`max_cells_lasso = NULL`；该版本帮助页把它写成3000属于陈旧文档，不代表
运行时会自动抽取3000个细胞。

### English

Version 6.1.2 changes the scientific defaults of the main scan toward the
values retained by later releases:

- `species` changes from mouse to human;
- `random_seed` changes from `NULL` to `20251009`;
- `bin_size` changes from 50 to 200;
- `parallel_cores` changes from 30 to 25;
- `downsample_gemli_ratio` changes from 0.75 to 1;
- `nmf_eval_downsample_cells` changes from 2000 to `NULL`.

There is still no separate `input_mode`, and the function still calls
`NormalizeData()`. The main function remains named `isCORAL_v6_1_0()` and
is not exported. Version 6.1.2 is therefore not a drop-in replacement for
6.1.0/6.1.1 based on its package version alone: both the public entry point
and default preprocessing need to be checked.

Version 6.1.2 also has two important implementation and packaging boundaries:

1. The source expression that reads counts in the `downsampleUMI` branch is
   corrupted by embedded path text. The branch remains dormant under the
   default `downsampleUMI = NULL`, but a numeric value is expected to fail at
   that line. Version 6.1.2 is therefore not a reliable UMI-downsampling
   release.
2. `NAMESPACE` exports only four names. It exports a nonexistent
   `SimulateMetaClone` symbol while failing to export the implemented
   `SimulateClone()`, `FindMetaclones()`, `SmoothEmbedding()`, and validation
   suite. Metaclone and plotting source remains in the package tree, but
   ordinary users should not regard it as a supported public API.

The validation source adds rasterization and includes the string `"0"` among
default invalid barcodes. However, the function is not exported, lacks
complete help-page packaging, and attempts to call `install.packages()` at
runtime if `ggrastr` is absent. These are packaging concerns, not a new
clone-reconstruction algorithm.

The computational source of `CORALembedding()` remains the same as in version
6.1.1. Its actual source default is still `max_cells_lasso = NULL`; the value
3000 shown in the help page of this release is stale documentation and does
not cause runtime subsampling to 3000 cells.

---

## 7. isCORAL 6.1.3

### 中文

6.1.3首次把预归一化输入明确纳入主扫描接口，并导出正式入口
`isCORAL_v6_1_3()`：

```r
input_mode = c("pre_normalized_data", "raw_counts")
pre_normalized_layer = "data"
check_raw_counts = TRUE
```

由于 `match.arg()` 选择第一个值，实际默认是
`input_mode = "pre_normalized_data"`。在该模式中，函数跳过
`NormalizeData()`，并把指定 layer视为预归一化表达；源码不验证该层是否
确实经过归一化、是否非负或是否适合本分析。在
`input_mode = "raw_counts"` 中，它尝试对能够取得的 counts 进行诊断性
检查，然后调用 `NormalizeData()`；该检查本身不保证 counts 存在。这使
6.1.3特别适合已经有适当规范化 `RNA/data`、但缺少 counts 的 CellTag对象。
assay读取已经根据 SeuratObject 主版本在 v5 `layer` 和 v4 `slot` 之间选择；
但是该版本不需要实现 counts 写回，因为 `downsampleUMI` 被禁用。

`downsampleUMI` 虽然仍在形参中，但6.1.3在两个输入模式下都明确忽略它并
给出警告。该决定避免默认 rescue workflow 对 counts 的依赖，但意味着
6.1.3无法执行 UMI 深度下采样。若研究设计必须下采样 counts，应使用有 counts
的对象和能够执行该分支的版本，或在包外完成并重新归一化。

6.1.3还带来以下变化：

- 增加并导出 `PlotPseudoclone()`，用于在现有 reduction 上显示一个或两个
  seed cell 对应的 pseudoclone probability field。
- 主扫描的 `run_parameters` 开始明确记录 isCORAL 版本、输入模式、是否执行
  归一化和所读取的表达层。
- `isCORAL_v6_1_0(...)` 只是转发到6.1.3；它保留旧名称的调用入口，却不会
  自动恢复6.1.0的 mouse、bin 50或 raw-count-style 默认行为。
- `CORALembedding()` 与6.1.2的计算源码相同。
- 与6.1.2相同，6.1.3帮助页把 `max_cells_lasso` 写成3000，但函数源码的
  实际默认值是 `NULL`；运行行为应以源码为准。

同时，6.1.3包中没有 `isCORAL_metaclone.R`，因此没有可用的
`FindMetaclone(s)`、`PlotValidationSuite()` 或整套 Metaclone 可视化。
`NAMESPACE` 仍导出了不存在的 `SimulateMetaClone`，而实际源码中的
`SimulateClone()` 和 `SmoothEmbedding()` 没有正确导出。这是公开 API
不完整的过渡版本。

### English

Version 6.1.3 is the first release to make pre-normalized input an explicit
part of the main scan interface and to export a canonical
`isCORAL_v6_1_3()` entry point:

```r
input_mode = c("pre_normalized_data", "raw_counts")
pre_normalized_layer = "data"
check_raw_counts = TRUE
```

Because `match.arg()` selects the first value, the effective default is
`input_mode = "pre_normalized_data"`. In this mode the function skips
`NormalizeData()` and treats the selected layer as pre-normalized expression;
the source does not verify that the layer is truly normalized, non-negative,
or otherwise suitable for the analysis. In `input_mode = "raw_counts"`, it
attempts a diagnostic check of any available counts and then calls
`NormalizeData()`; that check does not itself guarantee that counts exist.
This makes version 6.1.3 particularly suitable for CellTag objects that have
appropriate normalized `RNA/data` but no counts.
Assay reads are already selected between the v5 `layer` and v4 `slot` APIs
according to the SeuratObject major version. A version-aware counts setter is
not needed in this release because `downsampleUMI` is disabled.

Although `downsampleUMI` remains in the formal interface, version 6.1.3
explicitly ignores it with a warning in both input modes. This avoids making
the default rescue workflow depend on counts, but it also means that version
6.1.3 cannot perform UMI-depth downsampling. Studies requiring count
downsampling should use an object with counts and a version that implements
that path, or perform downsampling and renormalization outside the package.

Version 6.1.3 also introduces the following changes:

- It adds and exports `PlotPseudoclone()` for displaying pseudoclone
  probability fields associated with one or two seed cells in an existing
  reduction.
- The main scan's `run_parameters` begins to record the isCORAL version,
  input mode, whether normalization was performed, and the expression layer
  used.
- `isCORAL_v6_1_0(...)` only forwards to version 6.1.3. It preserves an old
  call name but does not restore the 6.1.0 defaults of mouse, bin 50, or
  raw-count-style preprocessing.
- The computational source of `CORALembedding()` is unchanged from version
  6.1.2.
- As in version 6.1.2, the 6.1.3 help page shows
  `max_cells_lasso = 3000`, whereas the actual function-source default is
  `NULL`; runtime behavior follows the source.

At the same time, `isCORAL_metaclone.R` is absent from the version 6.1.3
package. Consequently, there is no packaged `FindMetaclone(s)`,
`PlotValidationSuite()`, or complete Metaclone visualization suite.
`NAMESPACE` still exports the nonexistent `SimulateMetaClone`, while the
implemented `SimulateClone()` and `SmoothEmbedding()` are not correctly
exported. Version 6.1.3 should therefore be understood as a transitional
release with an incomplete public API.

---

## 8. isCORAL 6.1.4

### 中文

6.1.4重新整合了扫描、pseudoclone 可视化、Metaclone 和大数据分析，并把默认
输入改回：

```r
input_mode = c("raw_counts", "pre_normalized_data")
```

因此实际默认是 `raw_counts`。有 counts 的标准对象会读取 counts，可选执行
`downsampleUMI`，再调用 `NormalizeData()`；仅有预归一化 data 的 CellTag
对象必须显式设置：

```r
input_mode = "pre_normalized_data"
pre_normalized_layer = "data"
```

在预归一化模式下，`downsampleUMI` 会被忽略并警告，因为没有原始 counts
就无法定义 UMI 下采样。6.1.4的主扫描增加了 SeuratObject v5 `layer` 与
Seurat v4 `slot` 的版本兼容写回路径，并沿用6.1.3已经引入的版本引导读取。

该版本的重要功能包括：

- 正式入口 `isCORAL_v6_1_4()`，并保留转发到当前实现的旧名称
  `isCORAL_v6_1_3()` 和 `isCORAL_v6_1_0()`。
- 恢复并修正 `SimulateClone()`；`SimulateMetaClone()` 仅作为发出弃用
  提示后转发的可选别名。
- 将 embedding-based 社区函数正式命名为单数 `FindMetaclone()`。
  它在 `CORAL_clone_embed` 上建立 SNN 图，默认用
  `resolution = 0.7`、`dims = 1:30`、`algorithm = 4` 和
  `random_seed = 42` 聚类。
- `FindMetaclone()` 已经把结果写入默认 metadata 列 `metaclone`，并将其
  设为 active identity。正确用法必须接收返回值：

```r
sample <- FindMetaclone(sample)
```

  若只调用 `FindMetaclone(sample)` 而不重新赋值，外部的原对象不会保留
  新列。6.1.4缺少的是详细 `@misc` 运行记录，而不是 metadata 写入本身。
- 恢复并导出 `PlotValidationSuite()`、`PlotInferenceSuite()` 和其他
  Metaclone 图形；保留 `PlotPseudoclone()`。
- 增加 `isCORAL_scale()`，通过局部分块推断和全局 embedding 支持更大数据。
- 补齐源码显式使用的包依赖，并声明核心 doMC 扫描需要 Unix-like 系统。
- 对 `S_stab` 的细胞名和维度进行更严格检查，并在 `SimulateClone()` 中按
  Seurat barcode 对齐。

6.1.4的 `CORALembedding()` 与6.1.3在主要计算上相同，仍没有专属
`random_seed` 参数，也只保存七个历史缓存字段。加载时出现的
`RcppML::predict/head` 被 `stats::predict` 和 `utils::head` 替换警告来自
宽泛 namespace import。静态namespace解析表明这些名称最终指向预期函数，
且NMF继续显式调用 `RcppML::nmf()`；没有发现该警告引起的有意算法变化。
当前环境尚未完成成对运行回归。该冗余在6.1.5中被清理。

### English

Version 6.1.4 reintegrates the scan, pseudoclone visualization, Metaclone
workflow, and large-data support, while changing the default input back to:

```r
input_mode = c("raw_counts", "pre_normalized_data")
```

The effective default is therefore `raw_counts`. A standard object with
counts is read from counts, optionally processed by `downsampleUMI`, and then
passed through `NormalizeData()`. A CellTag object containing only
pre-normalized data must explicitly use:

```r
input_mode = "pre_normalized_data"
pre_normalized_layer = "data"
```

In pre-normalized mode, `downsampleUMI` is ignored with a warning because UMI
downsampling is not defined without original counts. The 6.1.4 core scan adds
version-compatible write-back through SeuratObject v5 `layer` or Seurat v4
`slot`, while retaining the version-guided read path introduced in 6.1.3.

Major functionality in this release includes:

- The canonical `isCORAL_v6_1_4()` entry point, with the old names
  `isCORAL_v6_1_3()` and `isCORAL_v6_1_0()` forwarding to the current
  implementation.
- Restoration and correction of `SimulateClone()`; `SimulateMetaClone()` is
  retained only as an optional deprecated alias that forwards to it.
- The canonical singular name `FindMetaclone()` for embedding-based
  community detection. It constructs an SNN graph in `CORAL_clone_embed` and
  defaults to `resolution = 0.7`, `dims = 1:30`, `algorithm = 4`, and
  `random_seed = 42`.
- `FindMetaclone()` already writes the result to the default `metaclone`
  metadata column and sets it as the active identity. The returned value must
  be assigned:

```r
sample <- FindMetaclone(sample)
```

  Calling `FindMetaclone(sample)` without reassignment does not preserve the
  new column in the caller's original object. What version 6.1.4 lacks is a
  detailed `@misc` run record, not metadata writing itself.
- Restoration and export of `PlotValidationSuite()`, `PlotInferenceSuite()`,
  and the other Metaclone plots, while retaining `PlotPseudoclone()`.
- Addition of `isCORAL_scale()` for larger datasets through local,
  partitioned inference and a global embedding.
- Completion of package metadata for explicitly referenced dependencies and
  declaration that the core doMC scan requires a Unix-like system.
- Stricter cell-name and dimension checks for `S_stab`, with Seurat-barcode
  alignment in `SimulateClone()`.

The main computation in `CORALembedding()` remains the same as in version
6.1.3. It still has no embedding-specific `random_seed` and stores only seven
historical cache fields. The load-time warnings in which
`RcppML::predict/head` are replaced by `stats::predict` and `utils::head`
arise from a broad namespace import. The calls resolve to the intended
functions in static namespace inspection, and NMF continues to call
`RcppML::nmf()` explicitly; no intentional analytical change caused by the
warning was found. Paired runtime regression has not been completed in the
current environment. This namespace redundancy is removed in version 6.1.5.

---

## 9. isCORAL 6.1.5

### 中文

6.1.5以6.1.4的核心分析定义为基础，重点增强旧脚本兼容、输入安全性、参数审计
和可复现性。主扫描 `isCORAL_v6_1_5()` 的参数表、默认输入模式和核心
`S_stab` 计算与6.1.4相同；原 `SimulateClone()` 的层次聚类计算路径在
`Simulate_clonal_expression_simple()` 下保留，`PlotPseudoclone()` 和主要
Metaclone 图形也没有被有意改写。

主要更新如下：

1. **历史入口兼容**

   - `FindMetaclone()` 仍是唯一正式的 embedding-based 聚类实现。
   - 恢复 `FindMetaclones()`，但它只发出弃用提示并转发到单数函数，不包含
     第二套算法。
   - `isCORAL_v6_1_4()`、`isCORAL_v6_1_3()` 和
     `isCORAL_v6_1_0()` 均转发到6.1.5核心。它们保证普通调用入口兼容，
     不是历史行为模拟器；其形参也变为 `...`，依赖 `formals()` 的程序需
     注意。
   - `Simulate_clonal_expression_simple()` 成为简单层次聚类的正式名称；
     `SimulateClone()` 以相同完整形参静默转发，`SimulateMetaClone()` 保留为
     弃用别名。算法、默认参数和 metadata 列名不变。

2. **`CORALembedding()` 可复现性**

   - 在末尾新增可选参数 `random_seed = NULL`，不破坏旧位置参数。
   - `NULL` 表示不调用 `set.seed()`，但 IRLBA、可选训练细胞抽样和五折
     分配仍会消耗随机数并推进全局 `.Random.seed`。指定整数时，在
     embedding入口设置一次种子。函数退出时不恢复调用者原先的RNG状态。
   - 五折 `foldid` 被显式生成、按训练细胞命名、传入 `cv.glmnet()` 并保存。
   - `max_cells_lasso = NULL` 仍表示全部有效细胞参与训练；缓存会同时记录
     `subsampling_applied = FALSE` 和实际训练细胞数。
   - 增加少于两个 memory genes、少于五个训练细胞、非法参数和不一致概率
     矩阵 dimnames 的提前错误。
   - worker 数安全限制为至少1，并在临时 doParallel cluster 关闭后尝试恢复
     调用者原先的 foreach 后端状态。

3. **Seurat v4/v5访问记录**

   - SeuratObject 5.x优先读取 `layer = "data"`，失败后回退到 `slot`。
   - 4.x或无法识别版本时优先读取 `slot = "data"`，失败后回退到 `layer`。
   - 实际使用的 API 和回退策略写入 `@misc$isCORAL_embedding`。
   - 与之前一样，embedding 固定读取 `RNA/data`，不读取 counts，也不调用
     `NormalizeData()`。

4. **完整但向后兼容的缓存**

   以下七个历史字段保留原路径：

```text
P_norm
original_MDS_coords
Y_ref_raw
valid_cells
top_genes
source_binsize
lasso_model
```

   同时增加 `run_parameters`、`effective_parameters`、
   `reproducibility`、`preprocessing`、`provenance` 和 `outputs`，
   保存调用参数、有效维度、训练细胞、CV folds、缩放统计、RNG状态、assay
   读取方式、lambda path、软件版本和 reduction 信息。该记录不包含原始输入
   矩阵副本或输入 checksum，因此不是自包含的完整复现实验包。

5. **Metaclone安全性和记录**

   - 聚类标签按 cell barcode 校验、重排后写入 metadata。
   - 精确历史默认值 `dims = 1:30` 遇到低于30维的 reduction 时，会警告并
     使用现有维度；其他显式越界值仍报错。
   - `@misc$isCORAL_metaclone$latest/history` 保存 requested/effective
     参数、图名、聚类数、实际 Leiden/Louvain 引擎和软件版本。
   - Leiden不可用时回退Louvain的行为在6.1.4已经存在。6.1.5记录实际引擎，
     但不同后端即使种子相同也不能假定结果相同。

6. **namespace清理**

   删除宽泛的 `import(RcppML)`，但保留依赖并继续显式调用
   `RcppML::nmf()`。这消除了 `predict` 和 `head` 的加载替换警告，而NMF
   调用路径仍然保留。

### English

Version 6.1.5 builds on the analytical definitions of 6.1.4 and focuses on
old-script compatibility, input safety, parameter auditing, and
reproducibility. The formal parameters, default input mode, and core `S_stab`
calculation of `isCORAL_v6_1_5()` are the same as in version 6.1.4. The
hierarchical computation formerly exposed through `SimulateClone()` is
retained under `Simulate_clonal_expression_simple()`. `PlotPseudoclone()` and
the principal Metaclone plots are likewise not intentionally rewritten.

The principal updates are:

1. **Historical entry-point compatibility**

   - `FindMetaclone()` remains the sole canonical embedding-based clustering
     implementation.
   - `FindMetaclones()` is restored, but only as a deprecated forwarding
     wrapper; it does not contain a second algorithm.
   - `isCORAL_v6_1_4()`, `isCORAL_v6_1_3()`, and `isCORAL_v6_1_0()` all
     forward to the 6.1.5 core. They preserve ordinary call-entry
     compatibility but do not simulate historical behavior. Their formal
     interface is also reduced to `...`, which matters to programs that
     inspect `formals()`.
   - `Simulate_clonal_expression_simple()` becomes the canonical name for the
     simple hierarchical method. `SimulateClone()` silently forwards with its
     full formal interface, while `SimulateMetaClone()` remains deprecated.
     The algorithm, defaults, and metadata-column names are unchanged.

2. **`CORALembedding()` reproducibility**

   - An optional final `random_seed = NULL` argument is added without
     disrupting old positional arguments.
   - `NULL` means that `set.seed()` is not called, but IRLBA, optional
     training-cell sampling, and five-fold assignment still consume random
     numbers and advance the global `.Random.seed`. A supplied integer sets
     the RNG once at embedding entry. The caller's prior RNG state is not
     restored on exit.
   - Five-fold `foldid` is generated explicitly, named by training cell,
     supplied to `cv.glmnet()`, and saved.
   - `max_cells_lasso = NULL` still means that all valid cells are used for
     training. The cache also records `subsampling_applied = FALSE` and the
     actual training-cell count.
   - Early errors are added for fewer than two memory genes, fewer than five
     training cells, invalid parameters, and inconsistent probability-matrix
     dimnames.
   - The worker count is safely limited to at least one, and the function
     attempts to restore the caller's prior foreach backend after closing its
     temporary doParallel cluster.

3. **Recorded Seurat v4/v5 access**

   - SeuratObject 5.x first reads `layer = "data"` and falls back to `slot`.
   - Version 4.x, or an unavailable version, first reads `slot = "data"` and
     falls back to `layer`.
   - The actual API and fallback strategy are written to
     `@misc$isCORAL_embedding`.
   - As before, embedding always reads `RNA/data`; it does not read counts or
     call `NormalizeData()`.

4. **Expanded but backward-compatible cache**

   The following seven historical fields remain at their original paths:

```text
P_norm
original_MDS_coords
Y_ref_raw
valid_cells
top_genes
source_binsize
lasso_model
```

   Additional `run_parameters`, `effective_parameters`, `reproducibility`,
   `preprocessing`, `provenance`, and `outputs` records store call arguments,
   effective dimensions, training cells, CV folds, scaling statistics, RNG
   states, assay access, the lambda path, software versions, and reduction
   information. The record does not include a copy of the original input
   matrix or an input checksum, so it is not a self-contained experiment
   archive.

5. **Metaclone safety and records**

   - Cluster labels are validated and reordered by cell barcode before being
     written to metadata.
   - If the exact historical default `dims = 1:30` is used with a reduction
     containing fewer than 30 dimensions, the function warns and uses the
     available dimensions. Other explicitly out-of-range requests remain
     errors.
   - `@misc$isCORAL_metaclone$latest/history` stores requested and effective
     parameters, graph names, cluster count, the actual Leiden/Louvain engine,
     and software versions.
   - The Leiden-to-Louvain fallback already existed in version 6.1.4.
     Version 6.1.5 records the effective engine, but results from different
     backends should not be assumed identical even with the same seed.

6. **Namespace cleanup**

   The broad `import(RcppML)` is removed, while the dependency and explicit
   `RcppML::nmf()` call remain. This eliminates the load-time `predict` and
   `head` replacement warnings while retaining the NMF call path.

---

## 10. 公开函数兼容性矩阵 / Public-function compatibility matrix

| 功能 / Function | 6.1.0 | 6.1.1 | 6.1.2 | 6.1.3 | 6.1.4 | 6.1.5 |
| --- | --- | --- | --- | --- | --- | --- |
| 正式主扫描入口 / Canonical scan entry | 源码有、未导出 / Source only | 源码有、未导出 / Source only | 源码有、未导出 / Source only | `isCORAL_v6_1_3` | `isCORAL_v6_1_4` | `isCORAL_v6_1_5` |
| `CORALembedding` | 导出 / Exported | 导出 / Exported | 导出 / Exported | 导出 / Exported | 导出 / Exported | 导出并扩展记录 / Exported, expanded records |
| `Simulate_clonal_expression_simple` | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 正式入口 / Canonical entry |
| `SimulateClone` | 导出 / Exported | 导出 / Exported | 源码有、未导出 / Source only | 源码有、未导出 / Source only | 导出 / Exported | 静默兼容别名 / Silent compatibility alias |
| `SimulateMetaClone` | 不存在 / Absent | 不存在 / Absent | 错误导出、无实现 / Exported but missing | 错误导出、无实现 / Exported but missing | 弃用别名 / Deprecated alias | 弃用别名 / Deprecated alias |
| `FindMetaclones` | 正式复数入口 / Canonical plural entry | 正式复数入口 / Canonical plural entry | 源码有、未导出 / Source only | 不存在 / Absent | 不存在 / Absent | 弃用包装器 / Deprecated wrapper |
| `FindMetaclone` | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 正式单数入口 / Canonical singular entry | 正式单数入口 / Canonical singular entry |
| `PlotPseudoclone` | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 导出 / Exported | 导出 / Exported | 导出 / Exported |
| `PlotValidationSuite` | 导出 / Exported | 导出 / Exported | 源码有、未导出 / Source only | 不存在 / Absent | 导出 / Exported | 导出 / Exported |
| `PlotInferenceSuite` | 不存在 / Absent | 不存在 / Absent | 源码有、未导出 / Source only | 不存在 / Absent | 导出 / Exported | 导出 / Exported |
| `isCORAL_scale` | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 不存在 / Absent | 导出 / Exported | 导出 / Exported |

“源码有、未导出”不等于稳定的公开函数；“错误导出、无实现”则表示该名称可能
出现在 `NAMESPACE` 中，但正常调用不会得到相应函数。6.1.2和6.1.3的这些问题
属于包接口缺陷，不表示其内部 scan 或 embedding 公式发生了同名变化。

“Source only” is not equivalent to a stable public function. “Exported but
missing” means that a name may appear in `NAMESPACE` even though no callable
implementation is present. These issues in versions 6.1.2 and 6.1.3 are
package-interface defects; they do not imply a correspondingly named change
to the internal scan or embedding formulas.

---

## 11. `CORALembedding()` 的跨版本稳定性 / Cross-version stability of `CORALembedding()`

### 中文

6.1.0至6.1.4的 `CORALembedding()` 具有相同的实际默认值：

```r
CORALembedding(
  seurat_obj,
  binsize_key = "binsize_200",
  n_dims = 30,
  top_n_genes = 3000,
  max_cells_lasso = NULL,
  alpha_elasticnet = 0.05,
  workers_glmnet = 5
)
```

其主要数学步骤均为：

1. 按 `omega_square_infer` 选择 top memory genes；
2. 将 `S_stab` 行归一化为 `P_norm`；
3. 在归一化表达的 memory-gene 空间计算细胞欧氏距离；
4. 由 pseudoclone probability distributions 构造平方能量距离；
5. 对距离矩阵执行经典 MDS；
6. 对表达矩阵执行细胞级 L2-normalization 和基因级 Z-score；
7. 使用 `family = "mgaussian"`、默认 `alpha = 0.05` 的多响应
   Elastic Net 将 MDS 坐标投影回纯表达空间；
8. 对正则化收缩后的各维进行方差回标；
9. 写入 `CORAL_cell_embed`，并以 `P_norm %*% Y_ref_raw` 写入
   `CORAL_clone_embed`。

6.1.5保留上述核心数学路径，同时增加末尾 `random_seed`、显式五折、版本
引导的assay访问、严格输入校验、worker下限、foreach后端恢复和详细记录。
6.1.0–6.1.4由 `cv.glmnet()` 内部生成 folds，且没有 embedding专属
种子参数。即使 `max_cells_lasso = NULL`、全部细胞参与训练，五折划分仍可能
涉及随机性；“没有训练细胞下采样”和“没有随机交叉验证”不是同一件事。

### English

Versions 6.1.0 through 6.1.4 have the same effective
`CORALembedding()` defaults:

```r
CORALembedding(
  seurat_obj,
  binsize_key = "binsize_200",
  n_dims = 30,
  top_n_genes = 3000,
  max_cells_lasso = NULL,
  alpha_elasticnet = 0.05,
  workers_glmnet = 5
)
```

Their principal mathematical steps are:

1. select top memory genes by `omega_square_infer`;
2. row-normalize `S_stab` to form `P_norm`;
3. calculate cell distances in normalized memory-gene expression space;
4. construct squared energy distances between pseudoclone probability
   distributions;
5. apply classical MDS to the distance matrix;
6. perform cell-wise L2 normalization and gene-wise Z-scoring;
7. project MDS coordinates back into pure expression space with
   multi-response Elastic Net using `family = "mgaussian"` and default
   `alpha = 0.05`;
8. recalibrate each regularized dimension to the MDS variance scale;
9. write `CORAL_cell_embed`, then write `CORAL_clone_embed` as
   `P_norm %*% Y_ref_raw`.

Version 6.1.5 retains this principal mathematical path while adding the final
`random_seed`, explicit five-fold assignment, version-guided assay access,
strict input validation, a worker lower bound, foreach-backend restoration,
and expanded records. Versions 6.1.0–6.1.4 allow `cv.glmnet()` to generate
folds internally and have no
embedding-specific seed argument. Even when `max_cells_lasso = NULL` and all
cells are used for training, five-fold assignment can remain stochastic.
“No training-cell subsampling” is not the same as “no random
cross-validation.”

---

## 12. 输入数据兼容性 / Input-data compatibility

### 中文

**有原始 counts 的标准数据**

- 6.1.0–6.1.2总是调用 `NormalizeData()`。
- 6.1.3需显式使用 `input_mode = "raw_counts"` 才进入该路径。
- 6.1.4–6.1.5默认就是 `input_mode = "raw_counts"`。
- 6.1.2的 numeric `downsampleUMI` 分支损坏，不应使用。
- 6.1.3完全忽略 `downsampleUMI`。
- 6.1.0、6.1.1、6.1.4和6.1.5可以在相应 raw-count 路径执行 UMI
  下采样，但前提是 counts存在且相应版本的assay读写接口可用。对
  6.1.0–6.1.1，还应保证形参 `assay` 与 DefaultAssay一致；同时需要固定
  数据、参数、细胞顺序和前置 RNG 条件。

**缺少 counts、只有预归一化 data 的 CellTag 数据**

- 6.1.3的默认输入模式可直接表达这一设计。
- 6.1.4–6.1.5必须显式设置
  `input_mode = "pre_normalized_data"`。
- 6.1.0–6.1.2的核心扫描没有正式的跳过归一化路径，不应作为此类输入的首选。
- 所有版本的 `CORALembedding()` 都读取归一化的 `RNA/data`，不需要 counts
  并且不会再次调用 `NormalizeData()`；但它仍需要先前生成的 scan cache。

**species 参数**

`species` 只影响线粒体和细胞周期基因的排除列表，不切换数学模型，也不执行
人鼠同源基因转换。把小鼠对象误写成 `"human"` 可能导致小鼠格式的
`mt-` 和细胞周期基因未被正确排除，从而改变候选特征，但不是把整个分析变成
“人类模型”。对 benchmark 的严重程度取决于这些基因是否进入候选集以及结果
对其敏感程度，应通过正确参数的成对重跑量化，而不应仅凭标签错误作结论。

**非 RNA assay**

六个版本的主扫描都接收 `assay` 参数。6.1.0–6.1.2调用
`NormalizeData()` 时没有传入该参数，因此归一化的是当时的 DefaultAssay，
随后读取的却是形参 `assay` 指定的 `data`；二者不同时，本次调用可能读取
一个没有被同步归一化的assay。6.1.3–6.1.5会把 `assay` 传给
`NormalizeData()`。但是，所有版本的 `CORALembedding()` 都固定读取
`RNA/data`。如果 scan 使用其他assay，embedding前必须确认相应表达空间已
合理地存在于RNA assay，否则两步分析可能不一致。

### English

**Standard data with raw counts**

- Versions 6.1.0–6.1.2 always call `NormalizeData()`.
- Version 6.1.3 requires explicit `input_mode = "raw_counts"` to enter this
  path.
- Versions 6.1.4–6.1.5 default to `input_mode = "raw_counts"`.
- The numeric `downsampleUMI` path in version 6.1.2 is defective and should
  not be used.
- Version 6.1.3 ignores `downsampleUMI` entirely.
- Versions 6.1.0, 6.1.1, 6.1.4, and 6.1.5 can perform UMI downsampling in
  their corresponding raw-count paths, provided that counts exist and the
  relevant version's assay access and write path are available. For versions
  6.1.0–6.1.1, the formal `assay` should also equal the DefaultAssay. Data,
  parameters, cell order, and the preceding RNG conditions must all be
  controlled for reproducibility.

**CellTag data with pre-normalized data but no counts**

- The default input mode in version 6.1.3 directly represents this design.
- Versions 6.1.4–6.1.5 require explicit
  `input_mode = "pre_normalized_data"`.
- The core scans in versions 6.1.0–6.1.2 have no formal skip-normalization
  path and should not be preferred for this input.
- `CORALembedding()` in every version reads normalized `RNA/data`, does not
  require counts, and does not call `NormalizeData()` again, although it
  still requires a previously generated scan cache.

**The `species` argument**

`species` only changes the exclusion lists for mitochondrial and cell-cycle
genes. It does not switch the mathematical model or perform human-mouse
ortholog conversion. Mistakenly labeling a mouse object as `"human"` can
leave mouse-formatted `mt-` and cell-cycle genes unfiltered and thereby alter
candidate features, but it does not turn the complete analysis into a “human
model.” Its seriousness for a benchmark depends on whether those genes enter
the candidate set and whether the results are sensitive to them. It should be
quantified by paired reruns with the correct setting rather than inferred from
the label error alone.

**Non-RNA assays**

The core scan in all six releases accepts an `assay` argument. Versions
6.1.0–6.1.2 call `NormalizeData()` without passing that argument, so
normalization is applied to the current DefaultAssay while the subsequent
`data` read uses the assay named by the formal argument. If they differ, the
scan may read an assay not normalized by that call. Versions 6.1.3–6.1.5 pass
`assay` to `NormalizeData()`. However, `CORALembedding()` in every version is
fixed to `RNA/data`. If the scan uses a different assay, users must confirm
that the corresponding expression space is appropriately represented in the
RNA assay before embedding; otherwise the two stages may be inconsistent.

---

## 13. Metaclone与模拟标签的迁移 / Migration of Metaclone and simulated labels

### 中文

推荐在6.1.5中使用：

```r
# 正式的 embedding-based Metaclone
sample <- FindMetaclone(
  sample,
  reduction = "CORAL_clone_embed",
  dims = 1:30,
  resolution = 0.7,
  random_seed = 42
)

# 正式的简单层次聚类版本
sample <- Simulate_clonal_expression_simple(
  sample,
  bin_size = 200,
  k_values = c(5, 10),
  hclust_method = "ward.D2"
)
```

历史名称的解释如下：

- 6.1.0–6.1.1脚本中的 `FindMetaclones()` 可在6.1.5继续运行，但会收到弃用
  提示并转发到 `FindMetaclone()`。
- 6.1.4脚本应使用单数 `FindMetaclone()`。
- 6.1.5中的 `SimulateClone()` 是静默兼容别名，`SimulateMetaClone()` 是弃用
  别名；新代码应使用 `Simulate_clonal_expression_simple()`。
- `FindMetaclone()` 和 `Simulate_clonal_expression_simple()` 不可互换：前者是 embedding
  上的图社区发现，后者是 `S_stab` 的层次切树。

### English

The recommended version 6.1.5 calls are:

```r
# Canonical embedding-based Metaclone
sample <- FindMetaclone(
  sample,
  reduction = "CORAL_clone_embed",
  dims = 1:30,
  resolution = 0.7,
  random_seed = 42
)

# Canonical simple hierarchical-clustering method
sample <- Simulate_clonal_expression_simple(
  sample,
  bin_size = 200,
  k_values = c(5, 10),
  hclust_method = "ward.D2"
)
```

Historical names should be interpreted as follows:

- `FindMetaclones()` calls from version 6.1.0–6.1.1 scripts remain accepted
  by version 6.1.5, but emit a deprecation warning and forward to
  `FindMetaclone()`.
- Version 6.1.4 scripts should use the singular `FindMetaclone()`.
- In version 6.1.5, `SimulateClone()` is a silent compatibility alias and
  `SimulateMetaClone()` is deprecated; new code should use
  `Simulate_clonal_expression_simple()`.
- `FindMetaclone()` and `Simulate_clonal_expression_simple()` are not interchangeable: the first
  detects graph communities in an embedding, whereas the second cuts a
  hierarchical tree derived from `S_stab`.

---

## 14. 历史分析复现建议 / Recommendations for reproducing historical analyses

### 中文

参数映射可以恢复历史版本的**预期分析配置**，但不能自动恢复每个历史版本的
所有实现细节。若目标是严格复现已发表结果，优先使用当时归档的包版本，并保留
原始输入对象、细胞和基因顺序、RNG入口状态、R/依赖版本、BLAS和聚类后端。

在6.1.5中近似映射历史主扫描配置时，可采用：

| 目标配置 | 6.1.5中需要显式设置的关键参数 |
| --- | --- |
| 6.1.0 | `input_mode="raw_counts"`, `species="mouse"`, `random_seed=NULL`, `bin_size=50`, `parallel_cores=30`, `downsample_gemli_ratio=0.75`, `nmf_eval_downsample_cells=2000`；在现有 `genes_to_remove` 基础上显式加入大写 `"MALAT1"`，以恢复GEMLI候选排除规则。标准小鼠符号 `Malat1` 通常不会被该规则匹配 |
| 6.1.1 | 与6.1.0的数值默认值相同，但不要自动加入 `"MALAT1"`；只有原历史调用本身显式传入时才保留 |
| 6.1.2 | `input_mode="raw_counts"`, `species="human"`, `random_seed=20251009`, `bin_size=200`, `parallel_cores=25`, `downsample_gemli_ratio=1`, `nmf_eval_downsample_cells=NULL`; 不要把6.1.2的损坏 `downsampleUMI` 行为当作分析定义 |
| 6.1.3 | `input_mode="pre_normalized_data"`, `pre_normalized_layer="data"`；如需模拟该版行为，不执行 UMI 下采样 |
| 6.1.4 | 6.1.5主扫描默认值已经匹配；仍需显式记录实际输入模式和所有非默认参数 |

这些映射不能保证逐位一致，原因包括：

- 6.1.5中的旧版本函数名都转发到当前核心；
- 6.1.5 embedding显式保存并传入CV folds；
- 早期版本未完整保存 RNG、训练细胞和依赖信息；
- glmnet、IRLBA、Seurat、BLAS及Leiden/Louvain后端版本可能改变数值；
- 并行执行和细胞/基因顺序也可能影响结果。

使用6.1.5近似6.1.0–6.1.2时，应默认令
`DefaultAssay(seurat_obj) == assay`。否则，早期版本“归一化DefaultAssay、
随后读取指定assay”的错配行为无法由当前参数表准确复现；表中的
`input_mode = "raw_counts"` 只代表最接近的正常分析配置。

还需区分核心 scan 和 embedding 的种子。主扫描的 `random_seed` 在 callr
子进程内注册给 doRNG；运行UID生成、`max_cells` 抽样以及 raw-count模式下的
`SampleUMI()` 都发生在此之前。因此，仅设置主扫描的 `random_seed` 不能
单独保证这些前置随机步骤完全复现。6.1.5的
`CORALembedding(random_seed = ...)` 是另一个独立控制项，不会继承
`FindMetaclone(random_seed = ...)` 的含义。

建议至少比较：

```text
CORAL_scan_results$run_parameters
ss_filter_average$omega_square_infer
S_stab
P_norm
original_MDS_coords
Y_ref_raw
CORAL_cell_embed
CORAL_clone_embed
lambda.min and glmnet coefficients
Metaclone labels and the effective clustering engine
```

### English

Parameter mapping can recover the **intended analytical configuration** of a
historical version, but it cannot automatically restore every implementation
detail of that release. Strict reproduction of a published result should
prefer the archived package actually used, together with the original input
object, cell and gene order, entry RNG state, R and dependency versions, BLAS,
and clustering backend.

The following settings can be used to approximate historical core-scan
profiles in version 6.1.5:

| Target profile | Key explicit settings in version 6.1.5 |
| --- | --- |
| 6.1.0 | `input_mode="raw_counts"`, `species="mouse"`, `random_seed=NULL`, `bin_size=50`, `parallel_cores=30`, `downsample_gemli_ratio=0.75`, `nmf_eval_downsample_cells=2000`; explicitly append uppercase `"MALAT1"` to the existing `genes_to_remove` vector to restore the GEMLI-candidate exclusion rule. The canonical mouse symbol `Malat1` normally does not match this rule |
| 6.1.1 | The same numerical defaults as 6.1.0, but do not add `"MALAT1"` automatically unless the historical call itself supplied it |
| 6.1.2 | `input_mode="raw_counts"`, `species="human"`, `random_seed=20251009`, `bin_size=200`, `parallel_cores=25`, `downsample_gemli_ratio=1`, `nmf_eval_downsample_cells=NULL`; do not treat the defective 6.1.2 `downsampleUMI` path as an analytical definition |
| 6.1.3 | `input_mode="pre_normalized_data"`, `pre_normalized_layer="data"`; omit UMI downsampling when emulating this release |
| 6.1.4 | Version 6.1.5 core-scan defaults already match; still record the actual input mode and every non-default argument |

These mappings do not guarantee bitwise identity because:

- old versioned function names in version 6.1.5 all forward to the current
  core;
- version 6.1.5 embedding explicitly stores and supplies CV folds;
- early versions did not fully record RNG, training cells, or dependencies;
- glmnet, IRLBA, Seurat, BLAS, and Leiden/Louvain versions can alter numerical
  results;
- parallel execution and cell/gene ordering can also matter.

When version 6.1.5 is used to approximate versions 6.1.0–6.1.2, the normal
assumption should be `DefaultAssay(seurat_obj) == assay`. Otherwise, the early
behavior of normalizing the DefaultAssay and then reading a separately
specified assay cannot be reproduced accurately by a current parameter
profile. The `input_mode = "raw_counts"` entries in the table represent the
closest normal analytical configuration.

The core-scan seed must also be distinguished from the embedding seed. The
scan's `random_seed` is registered with doRNG inside the callr subprocess,
whereas run-UID generation, `max_cells` sampling, and raw-count
`SampleUMI()` occur earlier. Setting only the scan's `random_seed` therefore
does not by itself make those preceding stochastic steps fully reproducible.
The `CORALembedding(random_seed = ...)` control added in version 6.1.5 is
independent and does not inherit the meaning of
`FindMetaclone(random_seed = ...)`.

At minimum, compare:

```text
CORAL_scan_results$run_parameters
ss_filter_average$omega_square_infer
S_stab
P_norm
original_MDS_coords
Y_ref_raw
CORAL_cell_embed
CORAL_clone_embed
lambda.min and glmnet coefficients
Metaclone labels and the effective clustering engine
```

---

## 15. 可视化兼容性与已知边界 / Visualization compatibility and known boundaries

### 中文

`PlotPseudoclone()` 从6.1.3开始提供，6.1.4与6.1.5的实现相同。它可以接收
一个或两个 seed、clone ID选种或显式 seed cell。在概率矩阵具有可用
dimnames时，函数按barcode对齐矩阵和reduction；若行名被判定为纯数字且矩阵
与Seurat细胞数完全相同，函数会按Seurat顺序重设行列名，原顺序不同时存在
条件性位置错配风险。该函数适合展示 transcriptome-only pseudoclone
inference产生的连续概率场，不应把色彩强度解释为实验clone的直接测量值。

`PlotValidationSuite()` 在6.1.0–6.1.1公开可用，6.1.2仅有未导出源码，
6.1.3缺失，6.1.4–6.1.5恢复公开。6.1.2以后默认把字符串 `"0"` 当作无效
barcode。6.1.4–6.1.5的 `na_barcodes` 只能向固定无效集合增加标签，不能
把内置 `"0"` 恢复为有效值。如果某个数据集把 `"0"` 用作真实clone label，
必须在绘图前重编码该标签，或修改函数实现。不要与 `PlotPseudoclone()` 中可
整体覆盖的 `invalid_clones` 参数混淆。

6.1.4–6.1.5的 `PlotValidationSuite()` 和 `PlotInferenceSuite()` 仍假定
Seurat对象、reduction和 `P_norm` 包含同一组、同一顺序的细胞。在标准的
全细胞、未重排流程中该假定成立；如果对象被 subset、reduction 被重排或缓存
只覆盖部分细胞，可能报长度/下标错误，或在重排情况下发生静默错配。正式比较
图生成前，应先核对：

```r
identical(
  rownames(Seurat::Embeddings(sample, "CORAL_clone_embed")),
  sample@misc$isCORAL_embedding$valid_cells
)
```

并检查 `P_norm` 两轴的 barcode。`PlotPseudoclone()` 自身具有更明确的
barcode交集对齐，不应与两个旧suite的实现边界混为一谈。

### English

`PlotPseudoclone()` is available from version 6.1.3 onward, and its
implementation is the same in versions 6.1.4 and 6.1.5. It accepts one or two
seeds, clone-ID-based selection, or explicit seed cells. When usable
probability-matrix dimnames are available, the matrix and reduction are
aligned by barcode. If row names are classified as purely numeric and the
matrix exactly matches the Seurat cell count, the function renames both axes
using Seurat order; a different original order creates a conditional
positional-misalignment risk. The function is designed to display the
continuous probability field produced by transcriptome-only pseudoclone
inference; color intensity should not be interpreted as a direct measurement
of an experimental clone.

`PlotValidationSuite()` is public in versions 6.1.0–6.1.1, source-only in
6.1.2, absent in 6.1.3, and restored as a public function in
6.1.4–6.1.5. Starting with 6.1.2, the string `"0"` is treated as an invalid
barcode. In versions 6.1.4–6.1.5, `na_barcodes` can only add values to the
fixed invalid set; it cannot restore the built-in `"0"` as valid. If a
dataset uses `"0"` as a genuine clone label, that label must be recoded before
plotting or the function implementation must be changed. This is distinct
from the fully replaceable `invalid_clones` argument of
`PlotPseudoclone()`.

In versions 6.1.4–6.1.5, `PlotValidationSuite()` and
`PlotInferenceSuite()` still assume that the Seurat object, reduction, and
`P_norm` contain the same cells in the same order. This assumption holds in a
standard full-cell, unreordered workflow. If the object is subsetted, the
reduction reordered, or the cache covers only a subset of cells, length or
index errors may occur; a reordered case may also be silently misaligned.
Before producing a formal comparison figure, check:

```r
identical(
  rownames(Seurat::Embeddings(sample, "CORAL_clone_embed")),
  sample@misc$isCORAL_embedding$valid_cells
)
```

and verify the barcodes on both axes of `P_norm`. `PlotPseudoclone()` itself
uses a more explicit barcode-intersection alignment and should not be
conflated with the implementation boundary of the two older suites.

---

## 16. 依赖、namespace与平台 / Dependencies, namespace, and platform

### 中文

包元数据的完整性也随版本发生变化：

| 版本 / Version | 静态包装状态 / Static packaging status |
| --- | --- |
| 6.1.0–6.1.1 | 源码显式使用的多个包没有完整列入 `DESCRIPTION`，包括核心 scan、embedding 和部分绘图所需依赖；主 scan 固定使用 doMC / Several packages explicitly referenced by the scan, embedding, and plots are not fully declared; the main scan uses doMC |
| 6.1.2 | 延续上述问题，并增加未声明的 `ggrastr` 使用；validation代码甚至会尝试运行时安装包，顶层 `library(ggdendro)` 也可能使干净环境加载失败 / Retains the earlier issues, adds undeclared `ggrastr`, attempts runtime installation, and may fail to load in a clean environment because of top-level `library(ggdendro)` |
| 6.1.3 | 删除旧Metaclone文件后绘图依赖减少，但核心 scan/embedding 的若干显式依赖仍未完整声明 / Removing the old Metaclone file reduces plotting dependencies, but several explicit scan/embedding dependencies remain undeclared |
| 6.1.4 | 静态 `pkg::function` 扫描未发现未声明的显式依赖；正式声明 Unix-like doMC平台要求，但宽泛导入 `RcppML` 会产生名称替换警告 / Static `pkg::function` inspection finds no undeclared explicit dependency; Unix-like doMC requirement is declared, although broad `RcppML` import produces masking warnings |
| 6.1.5 | 延续6.1.4的依赖清单，移除宽泛 `RcppML` import并消除 `predict/head` 警告 / Retains the 6.1.4 dependency list, removes the broad `RcppML` import, and eliminates the `predict/head` warnings |

6.1.0–6.1.3的依赖元数据不完整，意味着在开发者已经预装大量包的环境中
`load_all()`能够运行，并不能证明干净安装环境可复现。6.1.4–6.1.5改善的是
包安装和加载的可预期性，而不是改变 pseudoclone 数学模型。

所有版本的主 scan 实质上都采用 Unix/Linux fork-oriented doMC路径，没有
Windows原生PSOCK替代。`CORALembedding()` 的 glmnet阶段使用
doParallel/PSOCK风格，但这不使整个主工作流自动成为Windows兼容实现。
6.1.4–6.1.5已在 `SystemRequirements` 中明确这一边界。

依赖只被声明而没有精确锁定版本。即使 isCORAL源码相同，不同Seurat、
SeuratObject、glmnet、IRLBA、BLAS、RcppML或Leiden版本仍可能造成数值、
速度或标签差异。

### English

Package-metadata completeness also changes across releases:

| Version | Static packaging status |
| --- | --- |
| 6.1.0–6.1.1 | Several packages explicitly referenced by the core scan, embedding, and plots are not fully declared in `DESCRIPTION`; the main scan is fixed to doMC |
| 6.1.2 | The earlier issues remain, with additional undeclared use of `ggrastr`; validation code attempts runtime package installation, and top-level `library(ggdendro)` can prevent loading in a clean environment |
| 6.1.3 | Removal of the old Metaclone file reduces plotting dependencies, but several explicitly referenced scan and embedding dependencies remain undeclared |
| 6.1.4 | Static `pkg::function` inspection finds no undeclared explicit dependency; the Unix-like doMC requirement is declared, although broad `RcppML` import produces masking warnings |
| 6.1.5 | The 6.1.4 dependency list is retained, while the broad `RcppML` import is removed and the `predict/head` warnings are eliminated |

Incomplete dependency metadata in versions 6.1.0–6.1.3 means that successful
`load_all()` execution in a developer environment containing many preinstalled
packages does not establish reproducibility in a clean installation.
Versions 6.1.4–6.1.5 improve installation and loading predictability; they do
not thereby change the pseudoclone mathematical model.

The main scan in every release is effectively a Unix/Linux, fork-oriented
doMC workflow without a native Windows PSOCK alternative.
`CORALembedding()` uses a doParallel/PSOCK-style glmnet stage, but this does
not make the complete scan workflow automatically Windows-compatible.
Versions 6.1.4–6.1.5 state this boundary in `SystemRequirements`.

Dependencies are declared but not pinned to exact versions. Even with
identical isCORAL source, different versions of Seurat, SeuratObject, glmnet,
IRLBA, BLAS, RcppML, or Leiden can still alter numerical values, runtime, or
cluster labels.

---

## 17. 分析性能与兼容性可以确认到什么程度 / What can be claimed about analytical performance

### 中文

从源码可以确认：

- 6.1.1相对6.1.0没有改变除 GEMLI/NMF候选集合中的 `MALAT1` 排除规则以外
  的算法。
- 6.1.3、6.1.4之间的 pCORAL embedding主要计算没有改变。
- 6.1.4、6.1.5之间的核心 scan 和主要可视化计算没有被有意改变；原
  `SimulateClone()` 的层次聚类路径在6.1.5中以
  `Simulate_clonal_expression_simple()` 作为正式入口保留。
- 6.1.5保留6.1.4的 pCORAL核心公式和七个历史 embedding缓存字段。

不能仅从源码确认：

- 所有版本在任意数据集上产生完全相同的输出；
- 6.1.4和6.1.5逐位一致；
- benchmark数值、运行时间或内存使用完全不变；
- 所有Seurat v4/v5、glmnet、IRLBA、BLAS和聚类后端组合均等价。

输入模式和默认参数的变化本身就可能造成分析差异。6.1.5新增的严格校验会把
某些早期版本可能静默接受的异常输入改为明确错误；显式 fold记录也改变了可审计
方式。扩展 `@misc` 会增加 Seurat对象大小和RDS序列化时间，但不代表核心
pseudoclone公式发生改变。

因此，最严谨的表述是：

> 6.1.5以6.1.4的核心分析定义为基础，目标是在合法、同序且参数匹配的输入上
> 保持数值兼容，同时加强接口兼容、输入校验和可复现性记录。逐位一致、运行
> 性能完全相同和 benchmark完全不变仍需在目标R环境中进行成对回归验证。

### English

Source inspection confirms that:

- Version 6.1.1 does not alter the 6.1.0 algorithms apart from the `MALAT1`
  exclusion rule in the GEMLI/NMF candidate set.
- The principal pCORAL embedding computation is unchanged between versions
  6.1.3 and 6.1.4.
- The core scan and principal visualization calculations are not intentionally
  changed between versions 6.1.4 and 6.1.5. The hierarchical path formerly
  exposed through `SimulateClone()` is retained under the canonical
  `Simulate_clonal_expression_simple()` name in version 6.1.5.
- Version 6.1.5 retains the 6.1.4 pCORAL formulas and all seven historical
  embedding-cache fields.

Source inspection alone cannot establish that:

- every version produces identical output on every dataset;
- versions 6.1.4 and 6.1.5 are bitwise identical;
- benchmark values, runtime, or memory use are completely unchanged;
- every Seurat v4/v5, glmnet, IRLBA, BLAS, and clustering-backend combination
  is equivalent.

Changes in default input mode and parameters can themselves change an
analysis. The stricter checks in version 6.1.5 convert some invalid inputs that
an earlier release might have accepted silently into explicit errors, and
stored folds change auditability. Expanded `@misc` records increase Seurat
object size and RDS serialization time, but do not imply a change to the core
pseudoclone formulas.

The most rigorous statement is:

> Version 6.1.5 is based on the core analytical definitions of version 6.1.4
> and is designed to remain numerically compatible for valid, identically
> ordered, parameter-matched inputs, while improving interface compatibility,
> input validation, and reproducibility records. Bitwise identity, identical
> runtime performance, and completely unchanged benchmark results still
> require paired regression testing in the target R environment.

---

## 18. 版本选择建议 / Which version should be used

### 中文

- **新的常规分析**：优先使用6.1.5，并明确记录 `input_mode`、`species`、
  `bin_size`、所有下采样参数、embedding seed和实际聚类后端。
- **缺少 counts 的 CellTag 数据**：使用6.1.5并显式设置
  `input_mode = "pre_normalized_data"`；确认规范化表达确实位于
  `RNA/data`。
- **大规模数据**：使用6.1.5中的 `isCORAL_scale()`。若需保留完整局部扫描
  缓存，应显式设置 `save_local_results = TRUE` 并保存 `results_dir`；
  返回对象只保留compact subset summaries和全局结果，不包含完整局部
  `CORAL_scan_results`。
- **旧 `FindMetaclones()` 脚本迁移**：6.1.5可以直接接收旧调用，但新代码
  应改为 `FindMetaclone()`。
- **严格复现已发表的历史运行**：使用当时的归档版本，而不是依赖6.1.5中的
  旧函数名；同时恢复原输入、参数、顺序、随机状态和软件环境。
- **解释结果**：把 pseudoclone看作转录组推断的概率性谱系关联，把
  Metaclone看作该空间中的离散社区。二者均不应在没有 ground-truth验证时
  直接等同于实验测得的 clone。

### English

- **New routine analyses**: prefer version 6.1.5 and explicitly record
  `input_mode`, `species`, `bin_size`, every subsampling parameter, the
  embedding seed, and the effective clustering backend.
- **CellTag data without counts**: use version 6.1.5 with explicit
  `input_mode = "pre_normalized_data"` and confirm that normalized expression
  is actually present in `RNA/data`.
- **Large datasets**: use `isCORAL_scale()` in version 6.1.5. To retain full
  local scan caches, explicitly set `save_local_results = TRUE` and preserve
  `results_dir`. The returned object contains compact subset summaries and
  global results, not complete local `CORAL_scan_results`.
- **Migration of old `FindMetaclones()` scripts**: version 6.1.5 accepts the
  old call directly, but new code should use `FindMetaclone()`.
- **Strict reproduction of a published historical run**: use the archived
  release actually used rather than relying on an old function name inside
  version 6.1.5; restore the original input, parameters, ordering, RNG state,
  and software environment.
- **Interpretation**: treat a pseudoclone as a probabilistic lineage
  association inferred from transcriptomes and a Metaclone as a discrete
  community in that space. Neither should be equated directly with an
  experimentally measured clone without ground-truth validation.

---

## 19. 当前审查状态 / Current audit status

### 中文

本版本文档所依据的检查包括逐版本源码和导出比较、UTF-8文本检查、接口默认值
核对、核心计算路径比较、缓存字段比较、依赖声明静态扫描以及 Metaclone
metadata写入审查。

当前构建环境没有可用的 R/Rscript，因此尚未在此环境执行 `testthat`、
`R CMD check`、Seurat v4/v5运行矩阵、6.1.4/6.1.5成对数值回归或
benchmark重跑。测试文件的存在不等于测试已经执行。完成这些运行检查之前，
6.1.5应称为“通过源码静态审计的候选版本”，不应称为“已证明与6.1.4逐位
一致且性能完全相同”。

另有两项不影响分析公式、但会阻止正式公开发布的元数据问题：

- `DESCRIPTION` 中的 Author 和 Maintainer 仍是
  `Your Name <your.email@example.com>` 模板值；
- `License: MIT` 尚未配套填写MIT模板所需的实际版权持有人和年份。

项目负责人应在公开分发或提交软件仓库前提供真实信息；本次审查不猜测作者身份
或版权归属。

核心 scan及调用它的 `isCORAL_scale()` 仍依赖 Unix-like系统上的 fork-based
doMC后端。依赖版本没有被精确锁定，因此正式复现还应保存
`sessionInfo()`、包版本和外部聚类后端信息。

### English

The evidence supporting this document includes direct version-by-version
comparison of source and exports, UTF-8 text checks, verification of interface
defaults, comparison of core computational paths and cache fields, static
dependency inspection, and review of Metaclone metadata writing.

No usable R/Rscript is available in the current build environment.
`testthat`, `R CMD check`, a Seurat v4/v5 runtime matrix, paired 6.1.4/6.1.5
numerical regression, and benchmark reruns have therefore not been executed
here. The presence of test files is not evidence that the tests have run.
Until those runtime checks are completed, version 6.1.5 should be described as
a release candidate that has passed static source auditing, not as a release
proven bitwise identical to version 6.1.4 with identical performance.

Two metadata issues do not affect the analytical formulas but prevent a
formal public release:

- Author and Maintainer in `DESCRIPTION` remain the template value
  `Your Name <your.email@example.com>`;
- `License: MIT` is not yet accompanied by the actual copyright holder and
  year required to complete the MIT template.

The project owner should supply the real information before public
distribution or repository submission. This audit does not guess authorship
or copyright ownership.

The core scan, and consequently `isCORAL_scale()` when it invokes that scan,
still depends on a fork-based doMC backend on a Unix-like system. Dependency
versions are not pinned exactly, so formal reproduction should also preserve
`sessionInfo()`, package versions, and external clustering-backend
information.
