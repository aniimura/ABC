# -*- coding: utf-8 -*-
"""`ResearchPaper/0_Source/*.pdf` から `.txt` を作る **唯一の道具**。

★★★なぜ要るか(2026-09-06、メタ第 6 回 M11/M16 の実測)

`0_Source/*.txt` は **7 本の道具が読んでいる**
(hedge-index / bibmap / cycle-probe / full-graph / frdi-progress /
 genell-progress / paper-items)。うち `hedge-index` は CLAUDE.md が
「着手前に必ず数える」と定めている道具である。

ところが **`.txt` を書く道具はリポジトリに 1 本も無かった**。
`===== [page` を grep すると **読む側 7 本・書く側 0 本**。手順が外部にしか無く、
その結果「`.txt` は `pdftotext -layout` で作る」という**誤った契約**が
`paper-items.mjs` に書かれていた(2026-09-06 に訂正済み)。

★★実体は **PyMuPDF(`fitz`)の `page.get_text()` + 合字正規化 + `===== [page N] =====` の包み**
であることが第 6 回の逐語照合で判明した:

* FrdI p.25 を Xpdf 4.00 で抽出すると `φ` 0 / `∈` 0 / `→` 0 個。`.txt` には 17 / 6 / 10 個ある
* 6 通り(Xpdf/poppler × 既定/-layout/-raw)どれも不一致
* PyMuPDF が一致(標本 458 頁中 244 頁が完全一致、31 本は標本 5 頁すべて一致)

★**この道具は `.txt` を再現可能にするためのものである。**
`0_Source` は gitignore なので成果物は git に入らないが、**手順は入る**。

使い方:
    python tools/source-text.py "<pdf の絶対パスか 0_Source からの相対パス>"
    python tools/source-text.py --all          # .txt の無い PDF をすべて
    python tools/source-text.py --check <pdf>  # 既存 .txt と突き合わせるだけ(書かない)

★出力の様式(読む側 7 本がこれを前提にしている):
    ===== [page 1] =====
    (1 頁目の本文)
    ===== [page 2] =====
    ...

★合字の正規化: PyMuPDF は `ﬁ` `ﬂ` 等の合字をそのまま出すので、
既存の `.txt` に合わせて ASCII へ戻す。第 6 回が「差は合字だけ」と測った箇所である。
"""

import io
import os
import sys
import tempfile

# ★合字。第 6 回が「PyMuPDF と既存 .txt の差は合字だけ」と測った分。
LIGATURES = {
    "ﬀ": "ff",
    "ﬁ": "fi",
    "ﬂ": "fl",
    "ﬃ": "ffi",
    "ﬄ": "ffl",
    "ﬅ": "st",
    "ﬆ": "st",
}

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "ResearchPaper", "0_Source")


def normalize(s):
    for a, b in LIGATURES.items():
        s = s.replace(a, b)
    return s


def extract(pdf_path):
    """PDF から `===== [page N] =====` 包みのテキストを作る。"""
    import fitz  # PyMuPDF

    doc = fitz.open(pdf_path)
    out = []
    for i, page in enumerate(doc, start=1):
        out.append("===== [page %d] =====" % i)
        # ★`get_text()` は末尾に改行を付けるので、`join` と合わせると
        #   区切りの前に空行が余分に入る。既存の .txt に合わせて落とす。
        out.append(normalize(page.get_text()).rstrip("\n"))
    doc.close()
    return "\n".join(out) + "\n"


def write_atomic(path, text):
    """★`io.open(path, "w")` を直に使わない —— 開いた時点で切り詰めるので、
    符号化に失敗するとファイルが 0 行になる(2026-09-06 に 1 件やらかした)。
    先に encode し、一時ファイル経由で置き換える。"""
    data = text.encode("utf-8")
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    with os.fdopen(fd, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


def resolve(arg):
    if os.path.isabs(arg):
        return arg
    p = os.path.join(SRC, arg)
    return p if os.path.exists(p) else arg


def report(pdf, txt, text):
    pages = text.count("===== [page ")
    print("  %s" % os.path.basename(pdf))
    print("    頁 %d / %d 文字 / %d 行" % (pages, len(text), text.count("\n") + 1))
    if os.path.exists(txt):
        old = io.open(txt, encoding="utf-8", errors="replace").read()
        same = old == text
        print("    既存 .txt: %s (%d 文字)" % ("★一致" if same else "☆不一致", len(old)))
        return same
    return None


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0

    check_only = False
    if argv[0] == "--check":
        check_only = True
        argv = argv[1:]

    if argv and argv[0] == "--all":
        targets = []
        for e in sorted(os.listdir(SRC)):
            if not e.lower().endswith(".pdf"):
                continue
            pdf = os.path.join(SRC, e)
            txt = pdf[:-4] + ".txt"
            if not os.path.exists(txt):
                targets.append(pdf)
        print("★.txt の無い PDF: %d 本" % len(targets))
    else:
        targets = [resolve(a) for a in argv]

    for pdf in targets:
        if not os.path.exists(pdf):
            print("  見つからない: %s" % pdf)
            continue
        txt = pdf[:-4] + ".txt"
        text = extract(pdf)
        report(pdf, txt, text)
        if not check_only:
            write_atomic(txt, text)
            print("    → 書き出した: %s" % os.path.basename(txt))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
