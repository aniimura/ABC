---
name: txt-cache-vs-live-pdftotext-differs
description: 0_Source/*.txt(静的キャッシュ)とcheck.mjsのライブpdftotext呼び出しが、同じシェル・同じpoppler版でも同じページで食い違うことがある(プライム記号・マイナス記号)。
metadata:
  type: feedback
---

CorrHyp の `Theorem 4.2`(物理 p.10)で実測(2026-09-04、ABC3b)。
`ResearchPaper/0_Source/Correspondences on Hyperbolic Curves.txt`(既存の静的キャッシュ)を
読んで `(g′, r′)` や `2g′−2+r′`(U+2032 プライム、U+2212 マイナス)をそのまま Lean の
逐語引用に写したが、`tools/check.mjs` が呼ぶ**ライブの** `pdftotext -f 10 -l 10 …` は
**同じページで** `(g, r)`(プライム脱落)・`2g-2+r`(素の ASCII ハイフン)を返した。
Bash から実行・poppler 版も同じなので [[gate-shell-pdftotext-differs]](シェル間差)とは別の罠。

**Why:** `.txt` キャッシュがいつ・どのオプションで作られたか分からない
(本プロジェクトでも複数回、複数セッションが再生成している)。フォントの
ligature/kerning 由来の抽出ゆらぎは同一ファイル内の**同じ単語の別の出現**でも
文字が変わりうる——1箇所だけ様子が違っても驚かない。

**How to apply:**
- 逐語引用を書いたら `.txt` を読んだだけで満足せず、実際に
  `node tools/check.mjs`(または対象ページだけ `pdftotext -f N -l N` を直接呼ぶ)で
  ライブ照合すること。[[frdi-verbatim-ascii-only]] の「ASCII 断片に切る」は
  この手の揺れも一緒に潰す実務上の理由がある。
- 揺れが出た箇所は諦めてプライム/マイナス記号を素の文字(`g,r` / `-`)に落とし、
  「この箇所だけ pdftotext で脱落する」とコメントに残す(CorrHyp §4 Section4.lean の実例)。
- 大量の項目を一気に書いたときほど、この手の1〜2文字のズレを見落としやすい——
  [[report-progress-not-shortfall]] の逆で、**検査を後回しにしない**([[lean-build-check-discipline]])。
