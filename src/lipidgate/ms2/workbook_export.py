"""Shared result number formats and workbook writer."""
from __future__ import annotations
from datetime import datetime
from pathlib import Path
import pandas as pd
from lipidgate.result_file_io import group_header_values

def format_workbook_sheet(worksheet, *, header_row=1) -> None:
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter
    worksheet.freeze_panes = f'E{header_row + 1}'
    worksheet.auto_filter.ref = f'A{header_row}:{get_column_letter(worksheet.max_column)}{worksheet.max_row}'
    worksheet.sheet_view.showGridLines = False
    for row in worksheet:
        for cell in row:
            cell.fill = PatternFill('solid', fgColor='FFFFFF')
            cell.font = Font(name='Arial', size=11, color='000000', bold=cell.row <= header_row)
            cell.alignment = Alignment(vertical='top', wrap_text=True)
            if isinstance(cell.value,str) and cell.value.startswith('='):
                cell.data_type = 's'
    for i,cell in enumerate(worksheet[header_row],start=1):
        width = (70 if cell.value in {'matched_fragments', 'MS/MS 碎片注释', '同位素'} else
                 40 if cell.value in {'matched_name', 'source_file', 'Annotation', 'MS/MS 来源文件', '鉴定提示'} else 20)
        worksheet.column_dimensions[get_column_letter(i)].width = width
    header_to_index = {cell.value: index for index, cell in enumerate(worksheet[header_row], start=1)}
    for column_name, number_format in [
        ("rt_minutes", "0.000"),
        ("selected_ms2_rt", "0.000"),
        ("feature_rt", "0.000"),
        ("feature_rtmin", "0.000"),
        ("feature_rtmax", "0.000"),
        ("precursor_mz", "0.0000"),
        ("feature_mz", "0.0000"),
        ("mz", "0.0000"),
        ("ppm_error", "0.00"),
        ("mz_error_to_feature_ppm", "0.00"),
        ("rt_delta_sec", "0.0"),
        ("final_score", "0.00"),
        ("RT", "0.000"),
        ("m/z", "0.0000"),
        ("Score", "0.00"),
        ("同位素 MS1 RT (min)", "0.000"),
    ]:
        column_index = header_to_index.get(column_name)
        if column_index is None:
            continue
        for row in worksheet.iter_rows(
            min_row=header_row + 1,
            max_row=worksheet.max_row,
            min_col=column_index,
            max_col=column_index,
        ):
            row[0].number_format = number_format



def write_workbook(path: Path, sheets: dict[str, pd.DataFrame], *, permission_fallback: bool = False) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    def write(target):
        with pd.ExcelWriter(target, engine="openpyxl") as writer:
            for name, frame in sheets.items():
                groups = group_header_values(frame)
                frame.to_excel(writer, sheet_name=name, index=False, startrow=1 if groups else 0)
                worksheet = writer.sheets[name]
                if groups:
                    for column, value in enumerate(groups, 1):
                        worksheet.cell(1, column, value)
                format_workbook_sheet(worksheet, header_row=2 if groups else 1)
        return target
    try:
        return write(path)
    except PermissionError:
        if not permission_fallback:
            raise
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        return write(path.with_name(f"{path.stem}_{stamp}{path.suffix}"))
