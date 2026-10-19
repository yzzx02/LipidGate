# MS1 gap filling

源码开发版在跨样本对齐后，从原始 MS1 对缺失样本的特征窗口积分。GUI 的「MS1 特征」页提供
「启用 Gap filling」开关，默认开启；`lipidgate project-run` 使用同一组保存参数。
已完成的旧运行不会自动重新补峰，发布的 v1.0.3 EXE 尚未包含此功能。

## 方法与参考实现

核对了以下固定提交的源码：

- [XCMS `ChromPeakAreaParam` 补峰](https://github.com/sneumann/xcms/blob/b781f8b0c39a3310d74aebe5ae59d24f7e818431/R/methods-XCMSnExp.R#L2637)：
  分组后，仅对缺失的样本/特征组合提取原始信号，并保存 `is_filled`。其默认区域由各
  检出峰边界的分位数构成；采用校正后的扫描时间进行区域提取。
- [MS-DIAL LC-MS gap filler](https://github.com/systemsomicslab/MsdialWorkbench/blob/349d7ef6edbd09e34c0eb38b849bd9fd8af79e7e/src/MSDIAL5/MsdialLcMsApi/Algorithm/Alignment/LcmsGapFiller.cs#L42)：
  根据已有峰确定目标质量、时间和宽度；逆变换到原始时间后寻找附近的峰顶和边界，
  对边界内信号积分，保留补峰状态；另有强制插入分支。该提交的
  [`IsForceInsertForGapFilling` 默认值为 true](https://github.com/systemsomicslab/MsdialWorkbench/blob/349d7ef6edbd09e34c0eb38b849bd9fd8af79e7e/src/MSDIAL5/MsdialCore/Parameter/ParameterBase.cs#L1383)。

LipidGate 的 `feature_window_auc_v2` 是参考上述思路的独立实现，没有复制第三方代码，
也不声称与 XCMS 或 MS-DIAL 的数值完全相同：

1. 保留原始检测表和对齐表；仅检查面积为空或为 0 的样本单元格。该样本已有原生
   检测峰时不补，已有正面积不覆盖。无效负面积、无穷值保留原值并记录原因。
2. 在对齐时间轴上，取各原生检出峰左边界的 25% 分位数、右边界的 75% 分位数作为
   积分区域，参考 XCMS `ChromPeakAreaParam`。没有可用原生边界时依次使用特征表
   边界、中心 ± 备用 RT 半窗，不使用各样本峰顶的离散范围替代峰宽。
3. 有 **raw→aligned RT 映射** 时，将该窗口逆变换到各样本原始时间；只在映射覆盖
   范围内插值。不提供映射时，使用参考窗口并在两侧增加备用 RT 余量，审计记录
   `reference_window_with_rt_allowance`。已提供但无效的映射不参与补峰。
4. 提取质量窗口内最近的 centroid m/z 信号。**不要求再次找出合格的峰，不设 S/N、
   FWHM、最低峰高门槛**；弱信号、平坦背景、尖峰均可得到窗口积分。S/N 仅作诊断。
   平滑仅辅助确定邻峰谷底，不参与积分。
5. 明确的相邻峰以同文件实测谷底分界；两个局部峰同等接近目标时间时保留缺失。
   不占用已属于其他原生峰的信号，同一实测峰不能填到两个特征行。邻峰判断优先
   使用已记录的实测 EIC 峰顶，旧表才回退到检测器 RT。
6. 用未平滑的原始强度和实际扫描时间做梯形积分，单位为 `intensity*seconds`。
   同时记录直线边界基线扣除后的面积。主矩阵使用零基线面积，以保持与原生
   pyOpenMS 的原始 EIC 积分口径一致。至少需要两张实测扫描；全零、读取失败、
   明显采集断档等仍保留缺失。采集范围只覆盖部分窗口时只积分实际覆盖的部分，
   不补造扫描点，也不按名义窗口宽度放大面积。无峰形的非零背景也可能产生面积，
   这是取消再次找峰门槛后的预期行为。

这一步只恢复 MS1 定量。它不生成原生 Feature ID，不改变检测峰边界、MS2 关联、
分数、gate、候选顺序、逐谱置信度或合格 Top1 的样本投票。只有某张 MS2 自身具有
对应证据时，它才可能支持鉴定；补出的面积不能独立证明某样本已经鉴定该脂质。

旧版 `targeted_eic_v1` 的 S/N、FWHM、峰形和必须具有 RT 映射等硬门槛已撤销。
旧运行仍按原来保存的面积打开，不会自动重新积分。

## Fill %

主表、正式 CSV / Excel 和 GUI 导出增加 `Fill %`，候选表不增加此列：

`Fill % = 原始检出该特征的样本数 / 总 MS1 样本数`

与核对的 MS-DIAL 源码一致，范围为 0–1，浏览器保留两位小数：`0.50` 表示 50%。
MS-DIAL 在 [DataObjConverter.cs](https://github.com/systemsomicslab/MsdialWorkbench/blob/349d7ef6edbd09e34c0eb38b849bd9fd8af79e7e/src/MSDIAL5/MsdialCore/Utility/DataObjConverter.cs#L113)
先筛选 `PeakID >= 0` 的检测峰，再除以全部样本数；其
[导出器](https://github.com/systemsomicslab/MsdialWorkbench/blob/349d7ef6edbd09e34c0eb38b849bd9fd8af79e7e/src/MSDIAL5/MsdialCore/Export/IMetadataAccessor.cs#L119)
按 `F2` 输出。补峰不增加此值，补出面积也不会增加 MS2 鉴定支持票数。

原生检测成员按样本去重，检测到但面积为零的峰也算原始检出。总样本数包含参与
MS1 分析但没有检出峰的文件，排除仅含 MS2 的文件。外部 MS-DIAL 表保留已有的
`Fill %`；缺少原生成员及该字段的旧表，回退到补峰前面积矩阵的有效正值覆盖率。
没有完整 MS1 样本信息的纯 MS2 记录留空。

## 参数与后端

| 参数 | 默认 | 含义 |
| --- | ---: | --- |
| `enabled` | `true` | 对新运行启用补峰；关闭 MS1 时跳过 |
| `mz_ppm` | 5 | 原始 MS1 质量容差 |
| `rt_tolerance_sec` | 10 | 缺少边界时的备用 RT 半窗；缺少映射时两侧的 RT 余量；邻峰搜索上下文 |

GUI 暴露开关、质量误差和备用 RT 半窗。项目 JSON 的
`settings.ms1.gap_filling` 例如：

```json
"gap_filling": {
  "enabled": true,
  "mz_ppm": 5.0,
  "rt_tolerance_sec": 10.0
}
```

pyOpenMS 保存实际线性校正映射；XCMS 保存扫描级的原始/校正 RT 及原生峰成员。
XCMS 桥接脚本原来隐式调用 `fillChromPeaks()`，现由上述共享补峰阶段替代，以免丢失
补峰状态或重复填充。独立 `lipidgate detect` 只做检测/对齐，补峰由项目流水线执行。

旧项目中的 `min_snr`、`min_fwhm_sec`、`max_fwhm_sec`、`min_height` 可正常读取，
新补峰阶段忽略这些退休参数，GUI 保存时不再写回。MS1 原始检测的门槛不受影响。

Asari 或外部 MS-DIAL 表没有 RT 映射时也可使用上述参考窗口回退。该回退依赖样本
的原始 RT 足够接近参考；大幅漂移需要先进行可靠的 RT 校正。单样本原生峰表没有
跨样本缺口，本步骤不适用。Profile MS1 需要先质心化。

## 保存、展示与内存

每次运行在 `ms1/gap_filling/` 保存：

- `gap_filled_features.csv`：一行一个特征的完整面积矩阵及补峰样本。
- `gap_filling.csv`：每个尝试的样本/特征组合、成功/缺失/跳过/读取失败原因、
  目标坐标、参考窗口来源、RT 映射方式、实测最大值位置、积分边界、诊断 S/N、
  扫描数及积分。FWHM 不再作为门槛或计算结果，兼容旧审计时保留空字段。
- `filled_peaks.csv`：仅成功补峰的实测记录，供 EIC 按需读取。
- `gap_filling.json`：方法、参数、参考源码、文件身份和表/映射/审计 SHA-256。

正式 CSV / Excel 和 GUI 导出只在原有样本列中使用补回的面积，主表不增加补峰样本
列表或计数；GUI 详情也不额外列出这些字段。补峰来源、状态和坐标保留在上述项目
内部文件中。仍为一行一个特征，候选表不复制面积矩阵。选择相应样本查看 EIC 时
使用其实际最大值位置和积分边界。MS1 检出样本数仍只按原始检测矩阵统计，加合物关系与
MS2 关联继续读取原始检测表。

每个原始文件流式读取一次，只保留仍处于提取窗口内的紧凑 EIC 数组，不将整个
项目的扫描数组载入内存。提取窗口结束便释放 EIC；审计表、目标列表和输出矩阵
本身仍需要内存。浏览器仅缓存最多 128 个补峰坐标查询，不常驻全部补峰详情。

## 验证范围

回归覆盖多个数量级强度、窄/宽峰、RT 偏移、非线性逆映射、原生边界分位数、
不均匀扫描、部分采集窗口、邻峰谷底、重复占峰、原生峰保护、无信号、背景噪声、
尖峰、采集缺口、缺失/非单调 RT 映射、profile 数据、读取失败及 CSV/Excel
两条导出路径。实际 pyOpenMS 项目流水线对照开启/关闭补峰，检查面积恢复，同时
确认原始文件、原生检测成员、MS2 分数、gate、置信度和候选次序不受补峰影响。
还验证了空检测样本计入分母、重复检测去重、MS-DIAL 旧比例保留、旧参数迁移和
旧审计重新打开。

2026-10-19 的本地验证另有两组：

- 实际 R/XCMS 4.8.0 检测、Obiwarp 和 `fillChromPeaks(ChromPeakAreaParam())` 对照：
  三个合成弱峰的积分区域均一致，LipidGate 面积与区域内高斯解析面积相差约 0.10%。
  与 XCMS 面积相差约 +7.94%；XCMS 按强度和、区域宽度和扫描间隔计算，LipidGate
  按实际原始秒数做梯形积分，校正时间轴和端点处理也不同，不能声称数值相同。
- 11 个真实 centroid mzML，先实际检测/对齐，再以固定随机种子遮掉 500 个已检出
  样本/特征的面积及其检测成员，保留其余成员作为参考。补回 428/500 个，补回部分
  相对原面积的绝对相对偏差中位数为 28.7%；其中原本已有明确 EIC 谷底边界的
  243 个补回 219 个，偏差中位数为 14.7%。未补回的 72 个为原生邻峰占用 41 个、
  局部双峰歧义 30 个、采集断档 1 个。与其他原生峰冲突时仍不强行填充。
  这次同时尝试了所选特征中原本缺失的样本单元格，累计得到 2,784 个非零面积；
  这些原缺失单元格没有已知真值，不能据此计算正确恢复率。

实测中，边界未解析的原生峰仍使用检测器面积，和当前 EIC 窗口面积的口径差异较大；
该组遮峰偏差中位数为 46.5%。特征窗口也可能因 RT 校正残差截去部分信号。验证未
覆盖所有仪器、色谱条件、空白或标准品，不能把上述遮峰统计视为实际漏检的恢复率
或误填率。原始 mzML、原生峰表和对齐表均通过前后 SHA-256 检查确认未修改。

原始数据、派生 CSV、日志、图和验证脚本保留在忽略的 `.test_outputs/`，不上传
GitHub；源码检查点记录验证报告的指纹及汇总。
