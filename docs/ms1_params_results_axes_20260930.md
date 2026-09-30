# MS1 参数和结果页收尾（2026-09-30）

## 参数含义

pyOpenMS 的 `ElutionPeakDetection.min_fwhm/max_fwhm` 表示峰的半高全宽（FWHM），单位为秒；XCMS `CentWaveParam.peakwidth` 表示近似色谱峰宽范围，单位也是秒。界面按算法分别标明含义并原值传入，没有使用固定系数换算。参考 [OpenMS 参数文档](https://openms.de/documentation/classOpenMS_1_1ElutionPeakDetection.html)和 [XCMS 手册](https://www.bioconductor.org/packages/release/bioc/manuals/xcms/man/xcms.pdf)。

噪声阈值是 MS1 信号检测的强度门槛；最低特征峰高是检测完成后按峰顶强度筛选的独立门槛，默认 0 表示不额外筛选。pyOpenMS 使用特征的 `max_height` 元数据，不把积分面积当峰高。XCMS 在 centWave 检测后按峰表的 `maxo` 筛选峰顶强度，独立于 `snthresh`。XCMS `PeakDensityParam` 的 `minFraction` 与 `minSamples` 分别为跨样本比例和绝对样本数，仅选择 XCMS 时启用。参考 [XCMS 分组参数](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_density.html)。

如果输入 mzML 没有 MS1 扫描，流水线跳过该文件的 MS1 提取，仍对其 MS2 谱图执行匹配。若整个项目都没有 MS1 扫描，结果以 MS2-only 形式显示；MS1 图在对应 RT 范围提示没有 MS1 扫描。已有 MS1 的其他文件仍正常检测并对齐。

## 结果页

详情卡标题为“注释详情”，两列字段采用一致字号；主概览省去“支持样本”和“Source / scan”，完整审计信息仍可从“完整详情…”查看。置信度显示 High/Low；ECN 未判定显示 `-`。右下角直接使用 `MS1`、`MS2` 两个标签，默认 MS2；碎片图例固定单行，谱图纵轴标题保持在图内。

导航图和谱图横轴上按住右键从左向右拖动为缩小，从右向左为放大；在坐标轴上按住左键可平移。导航图画布内左键拖动仍可平移，点选仍可选择特征。

## 核对

使用 `D:\脂质匹配算法软件\dataset\LargeIntestine\Pos` 的真实正模式 mzML 确认扫描级 MS1 检测；用已有结果在 Windows 界面核对两列详情、标签、图例和坐标轴。截图与日志存于 `.test_outputs/parallel_tests_20260924/axis_eic_20260929/`。完整测试为 631 项及 3,118 个子测试通过；重新打包的 `dist/LipidGate.exe --self-test` 退出码为 0。没有重跑整个真实数据集。
