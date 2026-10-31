#!/usr/bin/env python
# 47b_supp_tables_docx.py —— 把 Supplementary_Tables.md 转成投稿用 DOCX
#   版式：A4 纵向 / 左右 2.0 cm / Times New Roman / 三线表 / 表格 8 pt
#   复用 46_md2docx.py 的字体与 Markdown 解析函数，保证与主稿风格一致
import os, re, importlib.util
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(BASE, "10_manuscript", "Supplementary_Tables.md")
OUT  = os.path.join(BASE, "10_manuscript", "Supplementary_Tables.docx")

# 复用主稿脚本的排版原语
_spec = importlib.util.spec_from_file_location("m2d", os.path.join(BASE, "00_code", "46_md2docx.py"))
m2d = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m2d)
add_runs = m2d.add_runs
set_font  = m2d.set_font
BODY = 10
TBL  = 8
TBL_WIDE = 7      # 列数 >= 10 的宽表


def three_line(t):
    """把 Table Grid 改成三线表：上下粗线，无竖线无内横线。"""
    tblPr = t._tbl.tblPr
    for old in tblPr.findall(qn("w:tblBorders")):
        tblPr.remove(old)
    b = OxmlElement("w:tblBorders")
    for edge, sz in (("top", "12"), ("bottom", "12"), ("left", "0"),
                     ("right", "0"), ("insideH", "0"), ("insideV", "0")):
        el = OxmlElement("w:" + edge)
        el.set(qn("w:val"), "none" if sz == "0" else "single")
        el.set(qn("w:sz"), sz)
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), "000000")
        b.append(el)
    tblPr.append(b)
    # 表头下细线
    for cell in t.rows[0].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tb = OxmlElement("w:tcBorders")
        bt = OxmlElement("w:bottom")
        bt.set(qn("w:val"), "single"); bt.set(qn("w:sz"), "6")
        bt.set(qn("w:space"), "0");    bt.set(qn("w:color"), "000000")
        tb.append(bt)
        tcPr.append(tb)


def repeat_header(t):
    """跨页时重复表头行。"""
    trPr = t.rows[0]._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader"); el.set(qn("w:val"), "true")
    trPr.append(el)


def set_col_widths(t, rows, avail_cm=17.0, min_cm=1.0):
    """按内容长度比例分配列宽，宽表不至于把长字符串挤成多行。"""
    ncol = len(t.columns)
    lens = []
    for j in range(ncol):
        mx = max((len(re.sub(r"\*\*|<sup>|</sup>", "", r[j])) if j < len(r) else 0)
                 for r in rows)
        lens.append(max(mx, 3))
    raw = [max(avail_cm * l / sum(lens), min_cm) for l in lens]
    k = avail_cm / sum(raw)
    widths = [w * k for w in raw]
    t.autofit = False
    for j, w in enumerate(widths):
        t.columns[j].width = Cm(w)
        for cell in t.columns[j].cells:
            cell.width = Cm(w)


def add_table(doc, rows):
    ncol = max(len(r) for r in rows)
    size = TBL_WIDE if ncol >= 10 else TBL
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(rows):
        for j in range(ncol):
            cell = t.cell(i, j)
            cell.text = ""
            p = cell.paragraphs[0]
            pf = p.paragraph_format
            pf.space_after = Pt(1); pf.space_before = Pt(1)
            pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
            txt = row[j] if j < len(row) else ""
            add_runs(p, txt, size=size, bold=(i == 0))
    three_line(t)
    repeat_header(t)
    set_col_widths(t, rows, avail_cm=17.0, min_cm=1.4 if ncol >= 10 else 1.0)
    return t


def para(doc, text="", size=BODY, bold=False, italic=False, align=None,
         space_after=4, space_before=0, keep=False, indent=None):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
    pf.space_after = Pt(space_after); pf.space_before = Pt(space_before)
    if indent is not None:
        pf.left_indent = Cm(indent)
    if align is not None:
        p.alignment = align
    if keep:
        pf.keep_with_next = True
    if text:
        add_runs(p, text, size=size, bold=bold, italic=italic)
    return p


def main():
    blocks = m2d.split_md(SRC)
    doc = Document()

    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.0)
    sec.top_margin = sec.bottom_margin = Cm(2.2)

    cp = doc.core_properties
    cp.author = "Mingjing Wang; Yi Chen; Richeng Quan; Lv Yan; Xiaomei Hu"
    cp.last_modified_by = "Xiaomei Hu"
    cp.title = "Supplementary Tables - Germline CHIP liability and myeloproliferative neoplasms"

    ntab = nrow = 0
    for kind, payload in blocks:
        if kind == "title":
            para(doc, payload, size=14, bold=True, space_after=8)
        elif kind == "h2":
            # 表题：粗体，与表同页
            para(doc, payload, size=10, bold=True, space_before=10,
                 space_after=3, keep=True)
            ntab += 1
        elif kind == "h3":
            para(doc, payload, size=10, bold=True, space_before=6, space_after=2, keep=True)
        elif kind == "table":
            add_table(doc, payload)
            nrow += len(payload) - 1
            para(doc, "", size=4, space_after=0)
        elif kind == "hr":
            pass
        else:
            para(doc, payload, size=BODY, space_after=4)

    doc.save(OUT)
    print(f"[OK] {OUT}")
    print(f"     表 {ntab} 张，数据行 {nrow} 行，文件 {os.path.getsize(OUT)} 字节")


if __name__ == "__main__":
    main()
