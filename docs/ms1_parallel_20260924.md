# MS1 峰提取并行评估（2026-09-24）

用户数据：`D:\脂质匹配算法软件\dataset\LargeIntestine\Pos`，共 46 份正模式 mzML，每份约 123–142 MB。本机为 Windows、16 逻辑处理器、32 GB 内存。基准文件通过同卷硬链接放在临时英文路径，源文件未移动或修改。

## pyOpenMS

对四份真实文件，默认参数下串行检测、对齐、分组共 6.02 秒。两线程试验约 5.64–6.77 秒，未表现出稳定提速；两进程试验为 4.31 秒，但进程之间必须用 featureXML 中转原生 `FeatureMap`。每文件原始特征表完全一致，对齐表在中转后出现强度、质量等字段的舍入。因此没有把两进程 pyOpenMS 峰提取接入正式工作流，保留原精度和默认串行流程。[OpenMS 并发说明](https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/THREAD_SAFETY.md)也提醒，释放 Python GIL 并不等于所有原生对象可安全共享。

## XCMS

本项目在 R 脚本中固定注册 `SerialParam`。试验中用 BiocParallel 的 `SnowParam` 配置 Windows 双进程，并显式传给 `findChromPeaks`。XCMS [官方手册](https://bioconductor.org/packages/release/bioc/manuals/xcms/man/xcms.pdf)支持 `BPPARAM`，BiocParallel [官方手册](https://bioconductor.org/packages/release/bioc/manuals/BiocParallel/man/BiocParallel.pdf)说明了 Windows 上的 `SnowParam`。

| 真实文件数 | XCMS 串行 | XCMS 双进程 | 导出核对 |
| ---: | ---: | ---: | --- |
| 2 | 87.84 秒 | 138.80 秒 | 1939 行完全一致 |
| 4 | 120.23 秒 | 204.71 秒 | 2350 行完全一致 |

对这批样本，双进程明显更慢，因此没有把 MS1 进程设置加入产品。

## 中文路径

本机 pyOpenMS 无法直接打开中文目录下的 mzML。现在原生读取前会在同卷临时英文目录创建硬链接，结束后清理；原始文件不复制。已从用户目录中的原始路径完成一份真实文件的项目全流程：pyOpenMS 得到 737 条特征，MS2 候选审计 2166 行，其中 800 行有 MS1 支持，最终导出 1535 行。MS2 读取也使用同一办法；输出仍保留原始 mzML 文件名。
