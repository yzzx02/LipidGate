# LipidGate

LipidGate 是用于 LC–MS 脂质数据的桌面软件，提供 MS1 特征提取与跨样本对齐、MS2 规则检索、可选 ECN 保留时间筛选，以及 MS1 EIC / MS2 谱图复核。

LipidGate is a desktop workbench for LC–MS lipid identification, aligned MS1 features, MS2 evidence review and optional ECN filtering. The Windows release includes Python, pyOpenMS and both reference libraries.

## 下载与运行

从 [v1.0.0 Release](https://github.com/yzzx02/LipidGate/releases/tag/v1.0.0) 下载 **LipidGate-v1.0.0-Windows-x64.zip**，解压后双击 `LipidGate.exe`。默认流程无需安装 Python、pyOpenMS 或另外导入谱库。正、负离子谱库均已内置；需要自己的 MSP 时才选择自定义库。首次启动会解压依赖，需要等待片刻。程序面向 Windows 64 位。

1. 创建项目，导入同一 LC 方法、同一离子模式的 mzML 文件。
2. 设置 MS1 参数；默认使用 pyOpenMS，也可导入 MS-DIAL 特征表。
3. 选择 MS2 模式、质量容差和分数门槛。需要时限制前体 m/z 范围，并在搜索窗口勾选加合物和脂质类型。
4. 选择是否使用 ECN 筛选，保存参数并开始分析。
5. 在结果页复核 MS1 EIC、逐次 MS2 谱图、匹配碎片和置信度，导出 CSV、Excel 或图像。

仪器 RAW 数据需先转换为 mzML。可选 XCMS 需要另外安装 R 及相应 R 包；Asari 需要独立环境。详见 [使用说明](docs/user_guide.md)。

## 特征与鉴定证据

跨样本注释按色谱峰归组。本样本未检出 MS1 特征时，同一峰的 MS2 仍可归入已有对齐特征；有谷底分隔的邻峰保留为不同特征。程序不会给未检出的样本借用其他样本的积分面积或峰边界。

`MS1 feature` 描述检测峰关联，`MS1 evidence` 描述原始扫描中的前体证据。`Confirmed precursor` 表示原始 MS1 扫描在质量和时间容差内确认了前体信号，不等于检测器已经输出可积分的峰。`High / Low` 同时考虑 MS1 证据、MS2 门控与结构分辨要求；得分高也可能是 Low。详见 [关联、积分与置信度](docs/identification_display_and_confidence.md)。

MS1 EIC 按需读取，默认显示 RT ±1 min；未关联的注释也能查看原始 MS1 信号。MS2 以 m/z 为横轴显示质谱。右拖 / 上拖放大，左拖 / 下拖缩小，强度轴从 0 开始缩放。

## 资源使用

内置谱库采用只读磁盘索引和有限缓存，只还原当前候选，不将百万条记录全部展开到内存。m/z 范围、加合物和类型限制进一步减少候选。默认 MS2 进程数为 1，多文件最多可设置 4 个；实际内存还取决于原始文件和候选数量。

Windows 包仅携带运行索引。源码保留无损压缩 MSP，供检查和重建。自定义大型 MSP 的首次解析可能需要更多内存。详见 [架构与资源策略](docs/architecture.md)。

## 从源码运行

开发环境使用 Python 3.13 构建并验证；项目声明支持 Python 3.10 及以上，其他版本需自行验证。Windows 普通用户优先下载程序包。

可下载 Release 中的 **LipidGate-v1.0.0-Source.zip**，它包含实际谱库文件；也可通过 Git LFS 克隆：

```powershell
git lfs install
git clone https://github.com/yzzx02/LipidGate.git
cd LipidGate
git lfs pull
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev,build]"
python scripts/build_prebuilt_libraries.py
python -m lipidgate gui
```

首次从 MSP 生成索引是开发构建步骤，需要额外时间和内存；后续运行复用索引。GitHub 自动生成的 `Source code (zip/tar.gz)` 可能只有 Git LFS 指针，请使用提供的 Source 包或执行 `git lfs pull`。

```powershell
python -m lipidgate --help
python -m lipidgate project-run --project "path/to/project"
python -m pytest -q
```

发布步骤见 [Windows 构建说明](docs/windows_exe.md)。谱库规则见 [MS2 evidence policy](docs/ms2_final_policy.md)，ECN 方法见 [ECN filtering](docs/ecn_filter.md)。

## 项目目录

| 路径 | 用途 |
| --- | --- |
| `src/lipidgate/` | 项目、GUI、MS1/MS2 工作流、ECN 和结果导出 |
| `src/lipidbench/` | 检测算法桥接及特征表工具 |
| `libraries/ms2/` | 正负模式压缩 MSP，使用 Git LFS |
| `scripts/` | 桌面打包、结构审计和谱库维护工具 |
| `tests/` | 检索、关联、积分、界面与导出的回归测试 |
| `docs/` | 使用和开发说明 |
| `config/` | 发布版本与验证指纹 |
| `assets/` | 字体及图标 |

## 许可证

LipidGate 项目代码使用 [MIT License](LICENSE)。第三方运行库和 Inter 字体保留各自许可证，见 [第三方说明](THIRD_PARTY_NOTICES.md) 和 `assets/fonts/LICENSE.txt`。
