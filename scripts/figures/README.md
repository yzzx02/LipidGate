# 脂质鉴定结果绘图

`generate_lipid_figures.py` 直接读取最终 Excel 工作簿，不依赖中间 CSV、JSON 或预览缓存。

默认只生成两张最常用的图：

```powershell
C:\Python313\python.exe scripts\figures\generate_lipid_figures.py
```

只生成指定图片：

```powershell
C:\Python313\python.exe scripts\figures\generate_lipid_figures.py `
  --figures ecn mass-defect iteration
```

生成全部五类图片：

```powershell
C:\Python313\python.exe scripts\figures\generate_lipid_figures.py --all
```

可选图片名称：

- `ecn`：保留时间过滤后的 ECN 分布，以及 RT 模型删除点。
- `mass-defect`：PC、PE、TG、FA、Cer、SM 六类脂质的原始/过滤后质量亏损对比。
- `iteration`：正负模式迭代 DDA 累计鉴定。
- `class-ring`：Top1 脂质类别层级环形图。
- `collision-energy`：正负模式不同碰撞能量的鉴定重叠图。

默认只输出 PNG，不生成中间数据文件。需要矢量图时可明确指定：

```powershell
C:\Python313\python.exe scripts\figures\generate_lipid_figures.py `
  --figures ecn `
  --formats png svg
```

可通过 `--raw`、`--filtered` 和 `--output` 指向下一批鉴定结果及输出目录。
