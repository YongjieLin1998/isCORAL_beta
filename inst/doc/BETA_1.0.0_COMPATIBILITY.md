# isCORAL beta 1.0.0 兼容性与可复现性说明
# isCORAL beta 1.0.0 Compatibility and Reproducibility Notes

## 中文

### 1. 版本定位

isCORAL beta 1.0.0 以归档的 isCORAL 6.1.5 源代码为功能基线。重打包没有有意修改主扫描、pseudoclone 共识计算、memory score、pCORAL embedding、成功执行的 Leiden/Louvain 聚类或绘图函数的数学路径。对外发布层面的变化为：

- 包版本设为 `1.0.0`，发布状态标记为 beta；
- 正式主扫描入口统一为 `isCORAL()`；
- 删除所有含内部开发版本号的主扫描入口，不保留同名别名；
- 将运行记录中的包版本写为 `1.0.0`，并增加 `release_status = "beta"`；
- 新增 `PlotMemoryGenes()` 作为 memory-gene 排名与可视化 API，不改变主扫描保存的统计量；
- 将 `FindMetaclone()` 的可选 Leiden 后端统一为 `leidenbase`；缺少时明确回退到 Louvain，且运行过程中不安装软件；
- README、帮助页、用户指南和测试使用新的公开入口与包名。

因此，beta 1.0.0 应被理解为 **6.1.5 功能实现的公开测试版封装**，而不是一套新的分析算法。

### 2. 公开接口兼容性

新的主扫描代码应使用：

```r
sample <- isCORAL(
  seurat_obj = sample,
  assay = "RNA",
  input_mode = "raw_counts"
)
```

旧脚本若调用含版本号的主扫描函数，需要将函数名改为 `isCORAL()`。除函数名外，beta 1.0.0 保留 6.1.5 主扫描函数的参数表、参数顺序和默认值。

以下下游与兼容入口可用：

- `PlotMemoryGenes()`：读取并可视化主扫描保存的 memory-gene 统计；
- `FindMetaclone()`：正式的 embedding-based Metaclone 聚类函数；
- `FindMetaclones()`：弃用包装器，转发到 `FindMetaclone()`；
- `Assign_deterministic_pseudoclone()`：正式的确定性层次聚类伪克隆分组函数；
- `simulate_clonal_expression()`：只绘制 isCORAL 结果的指定基因克隆表达 API；
- `SimulateClone()`、`Simulate_clonal_expression_simple()` 和 `SimulateMetaClone()`：仅为旧脚本保留的弃用兼容别名，均转发到 `Assign_deterministic_pseudoclone()`。
- `PlotPseudoclone()` 的 `seed_cells` 和 `seed` 参数同样为弃用别名；新代码分别使用 `seed_cell` 和 `random_seed`。

这些别名只保证调用迁移，不模拟某个历史包版本的全部行为。严格复现已完成的历史分析时，应保存并使用当时的包、输入对象、参数、随机数状态和依赖环境。

### 3. 输入数据兼容性

主扫描支持两种明确的输入模式：

| 模式 | 数据来源 | 主要行为 |
| --- | --- | --- |
| `input_mode = "raw_counts"` | 原始 counts | 检查 counts，并在当前 assay 上执行 Seurat `NormalizeData()`；UMI 下采样仅在该路径生效 |
| `input_mode = "pre_normalized_data"` | 已归一化的表达层 | 直接读取 `pre_normalized_layer`，不重复归一化；`downsampleUMI` 会被忽略并给出警告 |

缺少 counts 的数据必须显式使用 `pre_normalized_data`。这不是物种参数可以替代的选择。`species` 仍控制物种相关的基因过滤规则，必须与实际数据一致。

### 4. Seurat v4/v5 兼容性

6.1.5 基线中的表达矩阵访问逻辑被保留：

- 优先使用 SeuratObject 5.x 的 layer 接口；
- 在适用环境中回退到 Seurat v4 风格的 slot 接口；
- `CORALembedding()` 在 `@misc$isCORAL_embedding` 中记录请求的访问方式、实际使用的访问方式和回退信息；
- 各函数仍要求对象中存在所需的 assay、层、reduction、metadata 或 `@misc` 内容。

“兼容 Seurat v4/v5”表示源码实现了两类访问路径，不等同于所有 Seurat、SeuratObject、Matrix 和 R 版本组合都已经完成运行矩阵验证。

### 5. CORALembedding 兼容性

`CORALembedding()` 的主要默认值保持为：

```r
n_dims = 30
max_cells_mds = 10000
max_cells_lasso = NULL
workers_glmnet = 2
alpha_elasticnet = 0.05
random_seed = NULL
```

其中 `max_cells_lasso = NULL` 表示 Elastic Net 阶段不进行细胞下采样。`random_seed = NULL` 不重置调用者的随机数生成器；若需要可审计的重复运行，应显式传入整数种子。

embedding 结果继续写入 Seurat reductions，并在 `@misc$isCORAL_embedding` 中保存原有缓存字段和扩展 provenance，包括请求参数、有效参数、细胞和基因选择、训练细胞、交叉验证折分、lambda、缩放统计、并行后端、软件版本和随机数信息。beta 1.0.0 仅更新其中的发布版本标识。

### 6. FindMetaclone 兼容性

`FindMetaclone()` 保持 6.1.5 的正式接口与默认值：

```r
FindMetaclone(
  seurat_obj,
  resolution = 0.7,
  reduction = "CORAL_clone_embed",
  dims = 1:30,
  graph_prefix = "pcoral",
  cluster_col_name = "metaclone",
  algorithm = 4,
  random_seed = 42,
  verbose = TRUE
)
```

函数将聚类标签写入指定 metadata 列，将该列设为 active identity，并把请求参数、有效参数、聚类引擎、图名称、cluster 数量和软件版本写入 `@misc$isCORAL_metaclone`。当默认 `dims = 1:30` 超过实际 reduction 维数时，函数沿用 6.1.5 行为：警告后使用现有维数；其他显式越界请求仍报错。

`leidenbase` 保留在 `Suggests`，不会随普通依赖被强制安装。若要使用默认的 `algorithm = 4`，应在分析前运行 `install.packages("leidenbase")`，并使用已采用该后端的 Seurat >= 5.2.0。`FindMetaclone()` 本身不会联网或安装依赖；若无法加载 `leidenbase`，或 Seurat 版本过旧，它会明确警告，将有效算法设为 `algorithm = 1`，并把 Louvain 回退记录到 `@misc$isCORAL_metaclone`。这一依赖解析策略取代了旧的 Python `leiden` 尝试分支；Seurat v4 对象仍可使用 Louvain，此处不影响前述 layer/slot 输入兼容性。

### 7. 数值复现边界

在输入对象、细胞与基因顺序、参数、随机种子、依赖版本、并行环境和聚类后端全部一致时，beta 1.0.0 预期保持 6.1.5 的分析行为。以下因素仍可能导致结果不同：

- `input_mode`、assay 或 layer 不同；
- UMI、整细胞、候选基因或模块验证下采样设置不同；
- 随机种子或进入函数时的 RNG 状态不同；
- `cv.glmnet()` 折分、训练细胞或并行调度不同；
- Leiden/Louvain 可用性或 Seurat 默认行为不同；
- R、Seurat、SeuratObject、Matrix、glmnet、kernlab 或底层 BLAS 版本不同；
- 输入矩阵的细胞/基因顺序、稀疏表示或数值精度不同。

因此，“算法路径没有有意改变”不等同于“在任意环境中逐位一致”。发布前应在目标环境中完成 `R CMD check`、testthat、Seurat v4/v5 代表性对象测试和与归档 6.1.5 的成对数值回归。

### 8. 当前验证状态

本重打包过程执行源码级检查，包括文件级差异、公开导出、帮助页别名、UTF-8、文档链接和压缩包目录结构。若目标机器没有可用的 R 运行时，这些检查不能替代实际安装、加载和数值回归；在运行验证完成前，应将该包描述为测试版。

---

## English

### 1. Release identity

isCORAL beta 1.0.0 uses the archived isCORAL 6.1.5 source tree as its functional baseline. The repackaging does not intentionally change the mathematics of the core scan, pseudoclone consensus calculation, memory score, pCORAL embedding, successful Leiden/Louvain clustering, or plotting functions. The release-level changes are:

- the package version is `1.0.0` and the release status is beta;
- the canonical core entry point is `isCORAL()`;
- all core entry points containing internal development version numbers are removed without aliases;
- provenance records report version `1.0.0` and `release_status = "beta"`;
- `PlotMemoryGenes()` is added as a memory-gene ranking and visualization API without changing the statistics stored by the core scan;
- `FindMetaclone()` uses `leidenbase` as its supported optional Leiden backend, explicitly falls back to Louvain when it is unavailable, and never installs software at runtime;
- README files, help pages, user documentation, and tests use the public package name and entry point.

Beta 1.0.0 should therefore be understood as a **public beta packaging of the 6.1.5 functional implementation**, not as a new analytical method.

### 2. Public API compatibility

New core-scan code should use:

```r
sample <- isCORAL(
  seurat_obj = sample,
  assay = "RNA",
  input_mode = "raw_counts"
)
```

Scripts that call a version-numbered core function must rename that call to `isCORAL()`. Apart from the function name, beta 1.0.0 retains the 6.1.5 core-scan formal arguments, argument order, and defaults.

The following downstream and compatibility entries are available:

- `PlotMemoryGenes()` reads and visualizes memory-gene statistics stored by the core scan;
- `FindMetaclone()` is the canonical embedding-based Metaclone function;
- `FindMetaclones()` is a deprecated forwarding wrapper;
- `Assign_deterministic_pseudoclone()` is the canonical deterministic hierarchical pseudoclone-assignment function;
- `simulate_clonal_expression()` is the isCORAL-only selected-gene clonal-expression visualization API;
- `SimulateClone()`, `Simulate_clonal_expression_simple()`, and `SimulateMetaClone()` are deprecated compatibility aliases retained only for old scripts; all forward to `Assign_deterministic_pseudoclone()`.
- The `seed_cells` and `seed` arguments of `PlotPseudoclone()` are likewise deprecated aliases; new code should use `seed_cell` and `random_seed`, respectively.

These aliases support call migration only; they do not reproduce every behavior of an archived package. Strict reproduction of a completed historical analysis requires the original package, input object, parameters, RNG state, and dependency environment.

### 3. Input compatibility

The core scan supports two explicit input modes:

| Mode | Data source | Main behavior |
| --- | --- | --- |
| `input_mode = "raw_counts"` | Raw counts | Validates counts and runs Seurat `NormalizeData()` on the selected assay; UMI downsampling is available only on this path |
| `input_mode = "pre_normalized_data"` | A normalized expression layer | Reads `pre_normalized_layer` without repeating normalization; `downsampleUMI` is ignored with a warning |

Datasets without counts must explicitly use `pre_normalized_data`. The `species` argument controls species-dependent gene filtering and must match the biological source.

### 4. Seurat v4/v5 compatibility

The 6.1.5 matrix-access strategy is retained: SeuratObject 5.x layer access is preferred, with an applicable Seurat v4-style slot fallback. `CORALembedding()` records the requested and effective access paths in `@misc$isCORAL_embedding`.

Source support for both access styles does not establish that every combination of R, Seurat, SeuratObject, and Matrix versions has been runtime-tested.

### 5. CORALembedding compatibility

The principal defaults remain:

```r
n_dims = 30
max_cells_mds = 10000
max_cells_lasso = NULL
workers_glmnet = 2
alpha_elasticnet = 0.05
random_seed = NULL
```

`max_cells_lasso = NULL` means that the Elastic Net stage does not subsample cells. `random_seed = NULL` does not reset the caller's RNG; pass an integer seed when a self-contained reproducibility record is required.

The embedding still writes Seurat reductions and stores both legacy cache fields and expanded provenance in `@misc$isCORAL_embedding`. Beta 1.0.0 changes only the release identity recorded there.

### 6. FindMetaclone compatibility

`FindMetaclone()` retains the 6.1.5 interface and defaults. It writes labels to the requested metadata column, makes that column the active identity, and stores requested/effective parameters, clustering engine, graph names, cluster count, and software versions in `@misc$isCORAL_metaclone`.

If the historical default `dims = 1:30` exceeds the available reduction dimensions, the function warns and uses the available dimensions. Other explicit out-of-range requests remain errors.

`leidenbase` remains in `Suggests` and is not forced by ordinary dependency installation. To use the default `algorithm = 4`, install it before analysis with `install.packages("leidenbase")` and use Seurat >= 5.2.0, which adopts that backend. `FindMetaclone()` never downloads or installs dependencies. If `leidenbase` cannot be loaded or Seurat is too old, it warns, changes the effective algorithm to Louvain (`algorithm = 1`), and records the fallback in `@misc$isCORAL_metaclone`. This deterministic policy replaces the former attempt to use the Python-backed `leiden` package. Seurat v4 objects remain usable with Louvain, and this restriction does not alter the layer/slot input compatibility described above.

### 7. Numerical reproducibility boundary

Beta 1.0.0 is expected to preserve 6.1.5 analytical behavior when the input object, cell and gene order, parameters, seed, dependency versions, parallel environment, and clustering backend are identical. Changes in any of these factors can alter the result. In particular, input mode, downsampling, `cv.glmnet()` folds, Leiden/Louvain availability, package versions, BLAS, and matrix ordering must be controlled.

An unchanged intended algorithmic path is not a guarantee of bitwise equality in arbitrary environments. Release validation should include `R CMD check`, testthat, representative Seurat v4/v5 objects, and paired numerical regression against the archived 6.1.5 package.

### 8. Current verification status

The repackaging process performs source-level checks of file differences, public exports, help aliases, UTF-8 encoding, documentation links, and archive layout. If no usable R runtime is available on the packaging machine, these checks do not substitute for installation, loading, or numerical regression. The package should remain labeled beta until those runtime checks are completed.
