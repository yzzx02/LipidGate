# 下拉箭头与结果页布局（2026-09-29）

LipidGate 使用 PySide6；全局控件样式位于 `MainWindow._apply_style`，结果页另有局部样式。参数页、导出弹窗和结果页继续使用原生 `QComboBox` / `QTreeWidget`，没有改动选项及展开逻辑。

Qt 样式表可以定位 `QComboBox::down-arrow` 和 `QTreeWidget::branch`，但不支持 CSS 的 `transform: rotate(...)` 来直接画细线 V 形。Qt 官方示例使用 `image` 替换箭头。本项目原有 `combo_down.svg`，本次将其笔画收细到 1.45，并补齐悬停和树节点向右状态的三个极小 SVG；没有添加位图。

## 完整 QSS 示例

下面路径为示意。运行时代码会使用 `resource_path` / `icon_path` 拼出绝对路径，兼容源码运行和打包的 exe。

```css
QComboBox {
    padding: 0 34px 0 10px;
    border: 1px solid #cbd5e1;
    border-radius: 7px;
    background: #ffffff;
}
QComboBox:hover { border-color: #60a5fa; }
QComboBox::drop-down {
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 30px;
    border: none;
    background: transparent;
    border-top-right-radius: 7px;
    border-bottom-right-radius: 7px;
}
QComboBox::drop-down:hover { background: #f3f6fa; }
QComboBox::down-arrow {
    image: url("assets/icons/combo_down.svg");
    width: 18px;
    height: 18px;
}
QComboBox::down-arrow:hover {
    image: url("assets/icons/combo_down_hover.svg");
}

QTreeWidget#lipidClassTree::branch:closed:has-children {
    image: url("assets/icons/tree_right.svg");
    width: 18px;
    height: 18px;
}
QTreeWidget#lipidClassTree::branch:open:has-children {
    image: url("assets/icons/combo_down.svg");
    width: 18px;
    height: 18px;
}
QTreeWidget#lipidClassTree::branch:closed:has-children:hover {
    image: url("assets/icons/tree_right_hover.svg");
}
QTreeWidget#lipidClassTree::branch:open:has-children:hover {
    image: url("assets/icons/combo_down_hover.svg");
}
```

`QComboBox` 规则放在 `src/lipidgate/gui/app.py` 的全局主窗口 QSS 中；结果页的局部 QSS 在 `src/lipidgate/gui/result_browser.py`，也包含下拉规则与只针对 `lipidClassTree` 的树节点规则。两个作用域使用同一组资源，避免局部样式重新露出默认箭头。

## 同轮结果页改动

导航图去掉“RT–m/z 特征导航”标题，工具栏移至左上、按钮增至 36 像素，并收窄画布边距。MS/MS 图例中的“头基 / 类别”改为“头基”。注释详情由原来单列、可滚动的 `QTextBrowser` 改为两列字段概览；长值可悬停查看，“完整详情…”弹窗保留完整审计字段及碎片证据。2026-09-30 收尾时，两列字段字体统一，主概览移除了“支持样本”和“Source / scan”，置信度显示为 High/Low，ECN 无法判断显示为 `-`。

真实 Windows 界面、树节点和参数页截图集中在 `.test_outputs/parallel_tests_20260924/axis_eic_20260929/chevron_details_*`。截图使用 LargeIntestine/Pos 的已有最终结果，没有重跑 mzML。

回归验证：627 项测试、3,118 个子测试通过；在实际 Windows 窗口中核对了结果页、参数页和树节点展开/折叠两种状态。Qt 离屏测试不稳定地显示系统弹出层，因此用选项切换和模型项测试下拉框功能，用点击树节点分支测试展开/收起。
