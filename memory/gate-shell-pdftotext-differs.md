---
name: gate-shell-pdftotext-differs
description: tools/check.mjs の引用照合はシェルによって結果が変わる。PowerShell は poppler 25、Git Bash は xpdf 4.00 の pdftotext を拾う。ゲートは Bash で回すこと。
metadata:
  type: feedback
---

`node tools/check.mjs` の S4 / 引用照合は `pdftotext` を子プロセスで呼ぶ。
★★**この `pdftotext` は起動したシェルによって別物になる**(2026-08-18 実測):

| シェル | 解決先 | 版 |
|---|---|---|
| PowerShell | `…\WinGet\Packages\oschwartz10612.Poppler_…\poppler-25.07.0\Library\bin\pdftotext.exe` | poppler 25.07.0 |
| Git Bash | `/mingw64/bin/pdftotext` | xpdf/Glyph & Cog 4.00 |

★poppler 25 は `→` `≅` などの記号の出し方が違うため、**正しい引用が
軒並み「逐語が見つからない」で落ちる**。PowerShell で回すと
`ElementaryFrobenioid` `MorphismTypes` `Prop16` `Prop19` `GenEll/Lemma31`
`NCBelyi/*` `pGC section-1.html` が一斉に NG になり、
自己診断まで「D6 正しい入力は通る → ★落ちた(偽陽性)」と出る。
**Bash で同じコミットを回すと NG は 0 になる。**

**Why:** 一度これで「自分が壊した」と誤読しかけた。NG の中身が
**自分が触っていないファイルばかり**なら、まずシェルを疑う。

**How to apply:**

- ★ゲートは **Bash tool から** `node tools/check.mjs` で回す。PowerShell では回さない。
- NG の一覧に自分の変更ファイルが 1 つも無いときは、(a) シェル、(b) 平行セッション
  ABC3b の未追跡 WIP、の順に疑う。`git status --porcelain` で切り分く。
- 引用そのものの書き方は [[frdi-verbatim-ascii-only]]。判定の規律は
  [[lean-build-check-discipline]]。平行セッションの件は
  [[parallel-session-sweeps-my-files]]。
