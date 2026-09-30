# 氧化脂质碎片分池与正模式谱库修订

2026-09-07。本次为 `oxpe_chain_localization_20260907.md` 的后续修订。

## 全库检查范围

检查正式正、负谱库全部 Ox 前缀记录：OxPC、OxPE、OxPG、OxPI、OxPS、OxPGCN、OxTG、OxFA，共 114,046 条。

| 范围 | 发现与处理 |
| --- | --- |
| 负模式 OxPE | 酸／ketene 丢失仍为 other，保留上一轮修正。 |
| 负模式 OxPGCN | 已正确归 other，无需改库。 |
| 负模式 OxPC | 链丢失以 Neutral_Loss 存储，已在 other；补充带链标签的酸／ketene 丢失参与氧化定位与分链门控，计分仍仅在 other。 |
| 负模式 OxPG、OxPI、OxPS | 当前没有相应链丢失条目；采用母类分池继承，今后导入的 Diagnostic_FA_Loss 进入 other。 |
| 正模式 OxTG | 49,487 条记录的无链标签 M+H-H2O 原为 Diagnostic_FA_Loss，改为 Common。它不再充当必须命中的链丢失，也不计链池。真实链损失与三链专用门控保留。 |
| 正模式 OxPC [M+H]+ | 两条链的酸／ketene 丢失共四个，重建并保留链标签；移除原多余的氧化链 RCO 支持条目，加入两种母离子脱水支持峰。 |
| 正模式 OxPE [M+H]+ | 两条链各一个脱 HG 后 ketene 丢失和一个 RCO+，共四个链碎片；消除氧化链别名重复条目，RCO+ 按母类 PE 进入链池。 |
| OxFA、OxPC/OxPE 钠加合物 | 未套用质子化 PC/PE 的碎片重建规则。 |

## 正模式 PC/PE 衍生类

对 OxPC/OxPE [M+H]+，新增 `[M+H-H2O]+`、`[M+H-2H2O]+` 为 Common/other。
OxPE 另加 `[M+H-HG-H2O]+`，其中 HG = C2H8NO4P = 141.019094 Da。
这些是可匹配的支持离子，不是必须出现的峰，也不定位氧化链。氧数本身不能说明具体官能团或保证脱水发生。

酸 FA 按链的元素组成计算，FAk = FA - H2O；RCO+ = FA - OH - electron。
保留项目现有等价名称 `[M-(ROOH)+H]+(chain)`、`[M-(R=O)+H]+(chain)`、`[M-R=O-C2H8O4NP+H]+(chain)` 和 `(R=O)+(chain)`，避免破坏已有链标签识别。

氧化衍生类正模式继承相应母类的 HG 主导及支持池策略；OxPC/OxPE 的 HG/链/other 为 60/20/20，包括其现有钠加合物。正模式链丢失不按负模式规则放到 other。
双链氧化脂质的链池仍取两个最强独立实测峰的饱和强度质量分平均，不要求这两个评分峰属于不同链。

高置信度仍受独立门控约束：1O 两条链各至少一个独立峰；2O 以上氧化链至少两个，另一链至少一个；氧化链必须有非脱水的定位证据，同一实测峰不能同时补足两条链。因此两个同链峰本身不能满足完整氧化结构的确认门控。MS1 与 HG 要求保留。

## 验证与文件状态

- 549 tests、3118 subtests 通过；17 项为既有警告。
- 覆盖 1O/2O/3O，OxPC/OxPE 分池、准确质量、脱水不能补链门控、负模式母类继承、OxTG 脱水改类、负 OxPC 链丢失定位，以及归一化幂等性。
- 正模式 1,253,691 条总记录数不变，更新 64,487 条：OxTG 49,487、质子化 OxPC 7,500、质子化 OxPE 7,500。
- `scripts/curate_oxidized_library.py` 对输出逐条验证原记录顺序、指定规则转换结果、无额外记录及幂等性；非目标记录保留。
- 正式 `current_positive.msp.gz` 与 `current_positive.msp` 已同步。负模式磁盘谱库不需要重写；运行时分池与门控已更新。
- 库缓存版本升至 35，旧缓存不会继续加载。整批二维工作簿尚未重跑。
- gzip SHA256: `143826bcaa73974ab7287f31b8b3b7cc0bf0eea6843a7f364d6b28d552700bd8`。
- MSP SHA256: `c63be8c5acfd94df2c281857c363638f910e6a382e473b2eccaea5f643bc68e5`。

## 碎裂依据

用户指定的规则结合原有母类模板实施；可匹配离子不等于每一种氧化异构体必然产生该离子。

- [氧化 PC 的酸、ketene 及脱水证据原始研究](https://www.mdpi.com/1422-0067/16/4/8351)
- [氧化 PE 正负模式碎裂原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC11778248/)
- [氧化磷脂的水与过氧化氢中性丢失研究](https://sfrbm.org/site/assets/files/1240/lc-ms_oxidized-phospholipids_pitt_frbm2011.pdf)
