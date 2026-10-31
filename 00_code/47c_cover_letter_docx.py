#!/usr/bin/env python
# 47c_cover_letter_docx.py —— Cover Letter 转 docx / rtf / txt 三格式
#   与 P2 投稿包保持同一套交付习惯（docx 正式版、rtf 富文本、txt 可直接粘贴）
import os, re, importlib.util
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(BASE, "10_manuscript", "Cover_Letter_Leukemia.md")
OUTD = os.path.join(BASE, "10_manuscript")
OUT_DOCX = os.path.join(OUTD, "Cover_Letter_Leukemia_Hu.docx")
OUT_RTF  = os.path.join(OUTD, "Cover_Letter_Leukemia_Hu.rtf")
OUT_TXT  = os.path.join(OUTD, "Cover_Letter_Leukemia_Hu.txt")

_spec = importlib.util.spec_from_file_location("m2d", os.path.join(BASE, "00_code", "46_md2docx.py"))
m2d = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m2d)
add_runs = m2d.add_runs

SIZE = 11


# ---------- Markdown → 段落序列 ----------
def parse(path):
    """返回 [(kind, text)]，kind ∈ {h1, p, li, sign}"""
    out = []
    for raw in open(path, encoding="utf-8").read().split("\n"):
        s = raw.rstrip()
        if not s.strip():
            continue
        if s.startswith("# "):
            out.append(("h1", s[2:].strip()))
        elif re.match(r"^\d+\.\s", s):
            out.append(("li", s))
        elif s.startswith("**Xiaomei Hu") or s.startswith("Department") \
                or s.startswith("No. 1 Xiyuan") or s.startswith("E-mail"):
            out.append(("sign", s))
        else:
            out.append(("p", s))
    return out


def plain(text):
    t = text.replace("\\*", "*")
    t = re.sub(r"\*\*(.+?)\*\*", r"\1", t)
    t = re.sub(r"\*(.+?)\*", r"\1", t)
    t = re.sub(r"<sup>(.+?)</sup>", r"^\1", t)
    return t


# ---------- DOCX ----------
def build_docx(blocks):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.5)
    sec.top_margin = sec.bottom_margin = Cm(2.5)

    cp = doc.core_properties
    cp.author = "Xiaomei Hu"
    cp.last_modified_by = "Xiaomei Hu"
    cp.title = "Cover letter - Germline CHIP liability and myeloproliferative neoplasms"

    for kind, text in blocks:
        p = doc.add_paragraph()
        pf = p.paragraph_format
        pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
        pf.space_after = Pt(10 if kind != "sign" else 0)
        pf.space_before = Pt(0)
        if kind == "li":
            pf.left_indent = Cm(0.75)
            pf.first_line_indent = Cm(-0.75)
        add_runs(p, text, size=SIZE if kind != "h1" else 13,
                 bold=(kind == "h1"))
    doc.save(OUT_DOCX)
    return OUT_DOCX


# ---------- RTF ----------
def rtf_escape(s):
    out = []
    for ch in s:
        o = ord(ch)
        if ch in "\\{}":
            out.append("\\" + ch)
        elif o < 128:
            out.append(ch)
        elif o <= 32767:
            out.append("\\u%d?" % o)
        else:
            out.append("\\u%d?" % (o - 65536))
    return "".join(out)


def rtf_inline(text):
    """**bold** / *italic* → RTF 控制字"""
    parts, pos = [], 0
    for m in re.finditer(r"(\*\*.+?\*\*|\*[^*\n]+?\*)", text):
        if m.start() > pos:
            parts.append(rtf_escape(text[pos:m.start()]))
        tok = m.group(0)
        if tok.startswith("**"):
            parts.append("{\\b " + rtf_escape(tok[2:-2]) + "}")
        else:
            parts.append("{\\i " + rtf_escape(tok[1:-1]) + "}")
        pos = m.end()
    if pos < len(text):
        parts.append(rtf_escape(text[pos:]))
    return "".join(parts)


def build_rtf(blocks):
    head = (r"{\rtf1\ansi\ansicpg1252\deff0"
            r"{\fonttbl{\f0\fnil\fcharset0 Times New Roman;}}"
            r"\paperw11906\paperh16838\margl1418\margr1418\margt1418\margb1418"
            r"\f0\fs22")
    body = []
    for kind, text in blocks:
        size = r"\fs26" if kind == "h1" else r"\fs22"
        if kind == "li":
            body.append(
                r"\pard\fi-425\li425\sb0\sa200" + size + " " + rtf_inline(text) + r"\par")
        elif kind == "sign":
            body.append(r"\pard\sb0\sa0" + size + " " + rtf_inline(text) + r"\par")
        else:
            body.append(r"\pard\sb0\sa200" + size + " " + rtf_inline(text) + r"\par")
    open(OUT_RTF, "w", encoding="ascii", errors="ignore").write(head + "\n" + "\n".join(body) + "\n}")
    return OUT_RTF


# ---------- TXT ----------
def build_txt(blocks):
    lines = []
    for kind, text in blocks:
        lines.append(plain(text))
        if kind in ("h1", "p", "li"):
            lines.append("")
    with open(OUT_TXT, "w", encoding="utf-8-sig", newline="\r\n") as f:
        f.write("\n".join(lines).rstrip() + "\n")
    return OUT_TXT


def main():
    blocks = parse(SRC)
    for fn in (build_docx, build_rtf, build_txt):
        p = fn(blocks)
        print(f"[OK] {p}  ({os.path.getsize(p)} 字节)")
    words = len(plain(open(SRC, encoding="utf-8").read()).split())
    print(f"     段 {len(blocks)} 个，约 {words} 词")


if __name__ == "__main__":
    main()
