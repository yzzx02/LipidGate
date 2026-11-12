# 其他溶血磷脂的正模式碎片池（2026-11-12，Unreleased）

正模式常见 LPA/LPC/LPE/LPG/LPI/LPS 及现有 O/P 形式统一为 **HG60 / Other40**。
原 FAH 的评分碎片并入 Other，但保留源类型、名称、质量和脂链证据身份。
HG默认至少命中一半，所需数为 `ceil(HG通道数/2)`：2中1、3中2、4中2。
没有额外指定某个HG必须命中；只有一个HG时须命中该一个。Other不能替代HG，
不增加 Other 准入条件。本表权重是 LipidGate 的方法选择，
不能称为文献证明的最佳权重。LPC/LPE 细节见 `lpc_lpe_positive_pools.md`。

## 当前通道

| 类型 / 加合物 | 源记录数 | 保留的 HG 通道 | 有效权重 |
| --- | ---: | --- | --- |
| LPA `[M+H]+` | 164 | `P−H3PO4`、`P−H3PO4−H2O` | HG60 / Other40 |
| LPA-O `[M+H]+` | 54 | `P−H3PO4` | HG60 / Other40 |
| LPG `[M+H]+` | 164 | `P−C3H9O6P`、再失H2O | HG60 / Other40 |
| LPG `[M+NH4]+` | 164 | 源库 `P−C3H9O6P`，保留铵形式的产物 | HG60 / Other40 |
| LPG-O `[M+H]+` | 54 | `P−C3H9O6P`、再失H2O | HG60 / Other40 |
| LPI `[M+NH4]+` | 164 | `[M+H−C6H13O9P]+`，相对铵前体净损失277.0563 | HG60 / Other40 |
| LPI-O `[M+NH4]+` | 54 | 同上 | HG60 / Other40 |
| LPS `[M+H]+` | 164 | `P−C3H8NO6P`185.0089；Other仅保留失水、NL105.0426和前体 | HG60 / Other40 |
| LPS-O `[M+H]+` | 54 | `P−C3H8NO6P`185.0089、`P−C3H5NO2`87.0320 | HG60 / Other40 |

`P` 指该行的单电荷前体。LPI 铵前体还释放 NH3，不能将相对于 `[M+H]+`
的260.0297损失直接当作相对于 `[M+NH4]+` 的损失。C3H9O6P172.0137可写成
“脱水甘油 C3H6O2 + H3PO4”，不是游离甘油加 H3PO4；该拆写不规定反应先后。

## Other 和给分

- 原 LPA-H、LPG-H、LPI-NH4 的 FAH 通道全部移入 Other。现有的酰鎓离子、
  脂链损失、醚链支持、失水及前体仍按各自模板保留，没有借本次改权重增删质量。
- HG 沿用最强匹配峰的强度质量函数；Other 沿用实际物理匹配数和既有饱和上限。
  两池均充分命中为100分，满强度 HG-only 为60分。总分过滤独立于 HG 准入。
- 撤销这些正模式lyso记录的FAH-only低置信度准入例外；只有Other而没有HG
  不准入。最终输出FAH分数为0，不恢复第三个评分池。
- 脂链确认与准入分开。原来通过完整HG集合才能升到脂链层级的模板仍需完整
  集合；LPE-P命中一个HG可准入，但不会因此获得原需两个HG的脂链确认。
  最弱HG的额外强度限制不构成新的准入条件，正常实验峰阈值及总分过滤保留。

## 保留的 LPG 铵模板疑点

LPG(18:1) `[M+NH4]+`528.3296 的原 HG 为356.3159，即净损失172.0137，
名称为 `[M−C3H9O6P+NH4]+`。如果先失 NH3 再脱甘油磷酸组成，产物应为
339.2894，净损失189.0402。本次只改HG准入比例，保留原HG成员和质量，
不把质子化 LPG 的文献通道直接移作铵加合物的新强制条件。

## 可定位的原始研究图表

| 研究 | 可以核对的内容 | 适用限制 |
| --- | --- | --- |
| [Chen等，2017，Talanta162:530–539](https://doi.org/10.1016/j.talanta.2016.10.045) | **Fig.3上三行是LPC16:0、LPE16:0、LPS18:1的正模式标准品谱**。LPS为524.3→339.3（NL185.0089）、419.3（NL105.0426）、506.3（失水），另有265.3/247.2支持；3.2节正文明确LPG正模式NL172.0137 | Fig.3下三行LPA/LPI/LPG是负模式。LPS本次按用户指定只保留HG185损失及失水、105损失、前体；移除旧HG166.9984和两个酰链产物通道 |
| [Huynh等，2019，CellChem.Biol.26:71–84.e4](https://researchmgt.monash.edu/ws/portalfiles/portal/323259211/287614400_oa.pdf)，doi:10.1016/j.chembiol.2018.10.008 | **Table1的LPI为[M+NH4]+、NL277.0**；STAR Methods“Expanding the Analysis of Lipid Species in Human Plasma”用LPI13:0标准品确认NL277.0，并列LPE-P18:0标准品NL171.9 | 正模式标准品验证的特征转移，不是LPI完整产物离子谱图。277.0563为NH3+肌醇磷酸组成净损失的元素计算质量 |
| [Chen与Li，1994，J.Biochem.116:811–817](https://www.jstage.jst.go.jp/article/biochemistry1922/116/4/116_4_811/_pdf)，doi:10.1093/oxfordjournals.jbchem.a124600 | Fig.5a普通LPE16:0、5b烯醚LPE18:0；Fig.6 LPS18:0，526→341及失水508 | 正离子LSIMS-CID，非现代ESI-HCD；不照搬相对强度。其PA图为二钠加合物，不能称为质子化LPA标准谱 |
| [Morel等，2021，Anal.Chim.Acta1186:339088](https://pmc.ncbi.nlm.nih.gov/articles/PMC8581249/)，doi:10.1016/j.aca.2021.339088 | Fig.1C LPE-P、1D普通LPE；3.1节和补充SchemeS1讨论甘油磷酸组成损失与重排 | 正模式参考；通道角色与60/40仍由本项目定义 |
| [Lísa等，2011，J.Chromatogr.A1218:5146–5156](https://holcapek.upce.cz/reprints/RE_JCA_1218_2011_5146_lipidomics_2D_LCMS.pdf)，doi:10.1016/j.chroma.2011.05.081 | Table1列质子化LPG损失172、LPI损失260；Fig.6B为LPE馏分，Fig.7B为LPC馏分 | Table1支持特征损失；不是LPG/LPI各自完整标准品碎裂谱图 |
| [Wang等，2015，Anal.Chem.87:4879–4887](https://pmc.ncbi.nlm.nih.gov/articles/PMC4657943/)，doi:10.1021/acs.analchem.5b00410 | Fig.1标准LPI16:0/17:1、LPS16:0/17:1碎裂谱；Fig.4多类扫描 | Fig.1为负模式。LPC/LPE的部分正模式实验为钠加合物，不混入H/NH4模板 |

LPA正模式完整标准品碎裂谱本轮仍缺少可直接核查的原始图页。Murphy的
《Tandem Mass Spectrometry of Lipids: Molecular Analysis of Complex Lipids》，
[第5章](https://doi.org/10.1039/9781782626350-00130)，5.8.3节/Fig.5.12是进一步
检索线索；未取得原图核对前，不据此给LPA新增或强制任何通道。

LPC-O/P的正负模式参考另见 `lpc_ether_negative_fragments.md`。本次没有将其他
负模式溶血磷脂改为60/40，也没有用未定位到完整谱图的资料新增强制通道。
