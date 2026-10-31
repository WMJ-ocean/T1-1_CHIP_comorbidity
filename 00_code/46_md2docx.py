#!/usr/bin/env python
# 46_md2docx.py —— 把 Leukemia 格式手稿 Markdown 转成投稿用 DOCX
#   版式：A4 / Times New Roman 12pt / 双倍行距 / 连续行号 / 页码 / 表格文末
import os, re, sys
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING, WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(BASE, "10_manuscript", "Manuscript_Leukemia.md")
OUT  = os.path.join(BASE, "10_manuscript", "Manuscript_Leukemia_Submission.docx")

BODY_FONT = "Times New Roman"
BODY_SIZE = 12


# ---------- 低层排版工具 ----------
def set_font(run, size=BODY_SIZE, bold=False, italic=False, sup=False):
    run.font.name = BODY_FONT
    run._element.rPr.rFonts.set(qn("w:eastAsia"), BODY_FONT)
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if sup:
        run.font.superscript = True


TOKEN = re.compile(r"(\*\*.+?\*\*|<sup>.+?</sup>|\*[^*\n]+?\*)")


def add_runs(par, text, size=BODY_SIZE, bold=False, italic=False):
    text = text.replace("\\*", "*")
    pos = 0
    for m in TOKEN.finditer(text):
        if m.start() > pos:
            r = par.add_run(text[pos:m.start()])
            set_font(r, size, bold, italic)
        tok = m.group(0)
        if tok.startswith("**"):
            r = par.add_run(tok[2:-2]); set_font(r, size, True, italic)
        elif tok.startswith("<sup>"):
            r = par.add_run(tok[5:-6]); set_font(r, size, bold, italic, sup=True)
        else:
            r = par.add_run(tok[1:-1]); set_font(r, size, bold, True)
        pos = m.end()
    if pos < len(text):
        r = par.add_run(text[pos:])
        set_font(r, size, bold, italic)


def para(doc, text="", size=BODY_SIZE, bold=False, italic=False, align=None,
         space_after=0, space_before=0, first_indent=None, double=True, keep=False):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE if double else WD_LINE_SPACING.SINGLE
    pf.space_after = Pt(space_after)
    pf.space_before = Pt(space_before)
    if first_indent is not None:
        pf.first_line_indent = Cm(first_indent)
    if align is not None:
        p.alignment = align
    if keep:
        pf.keep_with_next = True
    if text:
        add_runs(p, text, size=size, bold=bold, italic=italic)
    return p


def add_line_numbers(section):
    sectPr = section._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "360")
    sectPr.append(ln)


def add_page_number_footer(section):
    footer = section.footer
    p = footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    fldChar = OxmlElement("w:fldChar"); fldChar.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText"); instr.set(qn("xml:space"), "preserve"); instr.text = "PAGE"
    fldEnd = OxmlElement("w:fldChar"); fldEnd.set(qn("w:fldCharType"), "end")
    run._r.append(fldChar); run._r.append(instr); run._r.append(fldEnd)
    set_font(run, 10)


def add_table(doc, rows, header=True, font_size=8):
    ncol = len(rows[0])
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    for i, row in enumerate(rows):
        for j in range(ncol):
            cell = t.cell(i, j)
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            txt = row[j] if j < len(row) else ""
            add_runs(p, txt, size=font_size, bold=(i == 0 and header))
    return t


# ---------- 解析 Markdown ----------
def split_md(path):
    raw = open(path, encoding="utf-8").read()
    lines = raw.split("\n")
    blocks = []          # (kind, payload)
    i = 0
    while i < len(lines):
        ln = lines[i]
        s = ln.rstrip()
        if not s.strip():
            i += 1; continue
        if s.startswith("|"):
            tbl = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in row if c):
                    tbl.append(row)
                i += 1
            blocks.append(("table", tbl)); continue
        if s.startswith("# "):
            blocks.append(("title", s[2:].strip())); i += 1; continue
        if s.startswith("## "):
            blocks.append(("h2", s[3:].strip())); i += 1; continue
        if s.startswith("### "):
            blocks.append(("h3", s[4:].strip())); i += 1; continue
        if s.startswith("---"):
            blocks.append(("hr", "")); i += 1; continue
        if re.match(r"^\d+\.\s", s):
            blocks.append(("ref", s)); i += 1; continue
        if s.startswith("**Figure") or s.startswith("**Table S"):
            blocks.append(("boldblock", s)); i += 1; continue
        # 普通段落：合并到下一个空行
        buf = [s]
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "|", "---")):
            if re.match(r"^\d+\.\s", lines[i].strip()):
                break
            buf.append(lines[i].strip())
            i += 1
        blocks.append(("p", " ".join(buf)))
    return blocks


def main():
    blocks = split_md(SRC)
    doc = Document()

    # 文档属性（避免投稿时显示 "python-docx"）
    cp = doc.core_properties
    cp.author = "Mingjing Wang; Yi Chen; Richeng Quan; Lv Yan; Xiaomei Hu"
    cp.last_modified_by = "Xiaomei Hu"
    cp.title = ("Germline CHIP liability causes myeloproliferative neoplasms "
                "independently of JAK2")
    cp.subject = "Mendelian randomisation and colocalisation"
    cp.comments = "Corresponding author: Xiaomei Hu (huxiaomei_2@163.com)"

    # 页面设置
    sec = doc.sections[0]
    sec.page_width  = Cm(21.0)
    sec.page_height = Cm(29.7)
    for attr in ("left_margin", "right_margin"):
        setattr(sec, attr, Cm(3.0))
    sec.top_margin = Cm(2.5); sec.bottom_margin = Cm(2.5)
    add_line_numbers(sec)
    add_page_number_footer(sec)

    st = doc.styles["Normal"]
    st.font.name = BODY_FONT
    st.font.size = Pt(BODY_SIZE)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), BODY_FONT)

    # ---------------- 标题页 ----------------
    idx_title = next(k for k, b in enumerate(blocks) if b[0] == "title")
    title = blocks[idx_title][1]

    para(doc, "Leukemia — Article (original research)", size=10, italic=True, space_after=12, double=False)
    para(doc, title, size=14, bold=True, space_after=12)

    # 作者与单位块：从 title 之后到第一个 h2 之前
    j = idx_title + 1
    meta = []
    while j < len(blocks) and blocks[j][0] != "h2":
        if blocks[j][0] == "p":
            meta.append(blocks[j][1])
        j += 1
    for m in meta:
        para(doc, m, size=11, space_after=6)

    # 分页 → 正文
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

    # ---------------- 正文 ----------------
    body_sections = {"Abstract", "Introduction", "Methods", "Results", "Discussion",
                     "Declarations", "References", "Tables", "Figure legends",
                     "Supplementary material"}
    in_refs = False
    for k, (kind, payload) in enumerate(blocks):
        if k <= j - 1:
            continue
        if kind == "h2":
            in_refs = (payload == "References")
            if payload in ("Tables", "Figure legends", "Supplementary material"):
                doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
            para(doc, payload, size=13, bold=True, space_before=12, space_after=6, keep=True)
            continue
        if kind == "h3":
            para(doc, payload, size=12, bold=True, italic=True, space_before=8, space_after=4, keep=True)
            continue
        if kind == "p":
            if payload.startswith("**Keywords:**"):
                para(doc, payload, space_after=10)
            elif payload.startswith("**Word counts:**"):
                continue
            else:
                para(doc, payload, space_after=6, first_indent=0.75)
            continue
        if kind == "table":
            add_table(doc, payload)
            para(doc, "", space_after=6, double=False)
            continue
        if kind == "ref":
            num, _, rest = payload.partition(".")
            p = para(doc, "", space_after=2)
            p.paragraph_format.left_indent = Cm(0.8)
            p.paragraph_format.first_line_indent = Cm(-0.8)
            add_runs(p, payload, size=10.5)
            continue
        if kind == "boldblock":
            para(doc, payload, space_before=8, space_after=4)
            continue
        if kind == "hr":
            continue

    doc.save(OUT)
    print("[OK]", OUT)


if __name__ == "__main__":
    main()
