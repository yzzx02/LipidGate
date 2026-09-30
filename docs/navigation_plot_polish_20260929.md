# 结果图视觉整理（2026-09-29）

当前结果页在 `result_browser.py` 中实例化 `qt_navigation.NavigationPlot` 和 `qt_spectrum.SpectrumPlot`；它们用 PySide6 的 `QPainter` 按屏幕像素比例绘制。`result_plots.py` 中还保留旧的 Matplotlib 类，但结果页不使用那两类。因此本次直接改当前画布，没有引入 `NavigationToolbar2QT`。

## 改动

- 导航散点按 GP、SP、GL、FA、ST、Other 使用固定配色：`#1f77b4`、`#ff7f0e`、`#2ca02c`、`#d62728`、`#9467bd`、`#7f7f7f`。普通点的透明度为 0.70、大小统一；先画数量多的类别，让少数类别能出现在上层。选中点有白色留白、深色边线和彩色中心。
- 图例置于绘图区外的右上方，窄窗口会自动折行。导航控件使用图标和悬停提示；框选后显示框内特征数。历史箭头在没有对应历史视图时禁用。随后工具栏移到左上，控件放大并移除导航图标题，见 [结果页控件样式](chevron_details_ui_20260929.md)。
- 导航图、MS/MS 图及 EIC 图保留白色背景，绘制明确的左轴、下轴及向外的刻度短线，去掉图内的横纵网格线；轴标题加粗。原有的十进制刻度、边界限制、左键平移和轴上右键拖动缩放继续使用。顺便修复了导航图轴上右键拖动没有记录按下状态的问题。

## 完整实现位置

- 导航图画布、工具栏和交互：`src/lipidgate/gui/qt_navigation.py`
- 家族色板：`src/lipidgate/gui/result_data.py`
- 图标：`src/lipidgate/gui/ui_icons.py`
- MS/MS 轴：`src/lipidgate/gui/qt_spectrum.py`
- EIC 轴：`src/lipidgate/gui/qt_eic.py`

实际 Windows 界面和导航画布截图、运行脚本集中在 `.test_outputs/parallel_tests_20260924/axis_eic_20260929/plot_polish_*`；推荐查看 `plot_polish_final_page_windows.png` 和 `plot_polish_final_navigation_windows.png`。截图使用 LargeIntestine/Pos 的已有最终运行结果，没有重新处理 mzML。完整回归为 627 项测试和 3,118 个子测试通过；后续补充的框选计数交互在 21 项定向测试中通过。
