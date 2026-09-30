# LipidGate

LipidGate 用于同一 LC 方法、单一离子模式下的脂质鉴定。

## 启动与流程

```powershell
cd "D:\Vscode Projects\LipidGate"
& "C:\Python313\python.exe" -m lipidgate gui
```

GUI 顺序：选择/创建项目 → 导入文件 → MS1 特征参数 → MS2 匹配、加合物及脂质类型 → 分数/ECN 过滤 → 结果查看（含导出）。

Windows 用户可直接运行 [exe 发布说明](docs/windows_exe.md)中的 `LipidGate.exe`，无需自行安装 Python。分析页会显示当前阶段、已完成文件数和已用时间，并按已完成文件的实际速度估算剩余时间。打包版包含预解析谱库快照，首次进入 MS2 匹配时会显示谱库准备状态，后续启动复用稳定缓存。

输入文件登记在项目中，不复制原始数据。当前分析引擎读取 mzML；运行前检查类型及离子模式。MS-DIAL 特征表在 MS1 页指定。每次运行使用开始时的参数快照，并写入新的 `项目/runs/时间戳/` 目录。项目参数保存在 `lipidgate.project.json`，重新打开可恢复。

## 固定核心（2026-09-20）

- 谱库：`libraries/ms2/current_positive.msp.gz`、`current_negative.msp.gz`；PS `[M+NH4]+` 已移除。
- 默认 MS1：pyOpenMS，5 ppm，noise 1000，S/N 5，半高全宽（FWHM）5–60 秒；最低特征峰高可独立设置。XCMS 的 5–60 秒为色谱峰宽，并可设置跨样本的最小比例和最小样本数。多文件检测保留各文件原始特征及 OpenMS 对齐成员关系；MS2 先关联同文件原始峰，再把注释回填到跨样本对齐特征。结果页按对齐特征显示一行，未对齐峰和 MS2-only 记录另列。
- 默认 MS2：母离子 5 ppm、碎片 15 ppm、最低相对强度 0.2%、Top 3。参数可改为 10/10 ppm；本轮不改变已确认的评分和门控规则。
- 分数过滤默认开启，阈值 50；ECN 默认关闭，可单独开启。关闭 MS1 特征检测时，输出为 MS2-only、低置信度。
- ECN 使用高置信度候选和已确认的有序系列规则；低置信度只回捞，不反向拟合。通过阈值默认 ±0.5 分钟。拟合点不足的亚类可选择保留，保留时 ECN 状态仍为“无法判断”。
- 图像：固定 DB 颜色、白底、每个 DB/C 只展示最贴近曲线的实际 RT；可选阈值带，默认 300 DPI。不把预测值伪装成观测点。
- 结果页的“导出结果与图像”集中选择当前筛选或全部特征、CSV / Excel 以及 ECN 图像和 DPI；特征表一行一个对齐特征或 MS2-only 谱簇，逐谱证据表保留全部扫描。分析时的原始审计 CSV 自动存入 `audit/`。
- 真假峰识别已移出产品代码，GUI/CLI 不再提供入口，也不再依赖 torch/torchvision。

详见 [核心审查与固定基线](docs/core_review_20260920.md)。哈希与库结构检查结果保存在 `config/core_baseline_20260920.json`。后续改动须运行回归检查，并明确更新基线，不能无记录地更换库、规则或图像参数。

## 命令行

```powershell
python -m lipidgate project-run --project "D:\LipidProjects\example"
python -m lipidgate filter-results --input "D:\LipidProjects\example\runs\时间戳\ms2\audit\ms2_candidates.csv" --output results\filtered
python -m lipidgate ms2-search --mode negative --mzml "D:\data\sample.mzML" --output results\ms2
python -m lipidgate detect --algo pyopenms --input "D:\data\mzml" --output results\ms1
```

`project-run` 与 GUI 使用同一流水线；`filter-results` 是当前普通项目筛选入口。旧 `ecn-filter` 保留为兼容命令，使用旧算法，不作为新项目入口。`filter-results` 需要原始候选审计表，不能用删减后的 17 列展示表代替。

多文件项目可在「MS2 参数」选择 1–4 个处理进程，默认 1。每个进程各自加载一次谱库并处理独立的 mzML 文件；结果按文件顺序合并，之后统一执行 MS1 特征关联、分数/ECN 筛选和导出。单文件始终串行。运行前检查可用内存，不足时提示减少进程数；实测及设计依据见 [MS2 并行处理说明](docs/parallel_ms2_20260924.md)。

LargeIntestine 正模式样本的 XCMS 双进程更慢，pyOpenMS 峰提取已很快且保持串行；对照数据见 [MS1 并行评估](docs/ms1_parallel_20260924.md)。跨样本对齐成员关系、9.06 分钟附近峰的修正及结果页速度见 [真实样本核对](docs/results_alignment_20260927.md)。Windows 中文目录的 mzML 读取也已兼容，不需移动原始文件。

正模式库对象已收紧内存布局；单个真实样本新流程的工作集峰值为约 5.49 GiB，进程容量按 6.5 GiB 估算。5.475 分钟附近三个 MS2-only 扫描现归为一行，仍可切换查看各原始谱图；正模式 LPS 的 164 条 `[M+NH4]+` 库记录已移除。实测与谱图复核见 [内存、ECN 与结果导出核对](docs/memory_ecn_export_20260928.md)。

结果查看默认仅显示已鉴定特征；正模式 LPE-P、LPE-O 门控及 PI 铵加合物碎片库已核对。指定 PI(18:0_20:3) 在样本 04 的重搜结果和与本机 MS-DIAL/LipidIN 参考库的比较见 [正模式 LPE 与 PI 库复核](docs/positive_lyso_pi_review_20260928.md)。

运行结果浏览默认只显示通过分析时分数/ECN 筛选的逐谱候选；明确打开审计表才显示被淘汰候选。实验谱图、轴缩放和脂质类别双击筛选的核对见 [结果页分数与谱图复核](docs/result_browser_score_visual_review_20260928.md)。

结果图的轴上右键连续缩放、导航图左键平移和按需显示的 MS1 EIC 见 [图谱交互与 EIC 核对](docs/result_axes_eic_20260929.md)。

导航图的新配色、无网格坐标轴、图标工具栏及完整代码位置见 [结果图视觉整理](docs/navigation_plot_polish_20260929.md)。

下拉框和脂质类别树的细线箭头、导航工具栏与两列详情布局见 [结果页控件样式](docs/chevron_details_ui_20260929.md)。

MS1 峰宽、峰高与 XCMS 样本门槛的含义，以及结果页最终交互见 [MS1 参数和结果页收尾](docs/ms1_params_results_axes_20260930.md)。无 MS1 扫描的 mzML 会跳过 MS1 提取，继续 MS2 匹配；对应的 MS1 图显示无可用扫描。

## 安装与发布范围

```powershell
python -m pip install -e .[dev]
git lfs pull
python -m pytest -q
```

正式运行代码在 `src/lipidgate/`，MS1 复用 `src/lipidbench/`。Asari、R/XCMS 需另外安装对应运行环境。谱库由 Git LFS 管理，直接读取 gzip MSP；解析缓存位于 `%LOCALAPPDATA%\LipidGate\Cache`（可用 `LIPIDGATE_CACHE_DIR` 覆盖）。

二维实验的分段、归一化和正负模式合并属于本地专项脚本，不进入普通项目流程。历史二维结果不因本轮 GUI 更新而自动重算。旧真假峰源码在 `archive/peak_truth/`，历史模型在 `models/peak_truth/`，均不属于产品运行依赖。发布导出边界见 [release_scope.md](docs/release_scope.md)。
