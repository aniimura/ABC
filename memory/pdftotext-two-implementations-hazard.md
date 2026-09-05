---
name: pdftotext-two-implementations-hazard
description: check.mjs は PDF を pdftotext で読む(2実装の罠・キャッシュが本命、第1024で対策済)。★0_Source の .txt は pdftotext 製ではなく PyMuPDF 製——別物である
metadata:
  type: reference
---

## ★★★まず区別すること(2026-09-06 の訂正。ここを混ぜると測定が壊れる)

**この木には原文を読む経路が 2 本あり、読んでいる文字列が違う。**

| 経路 | 誰が使うか | 実体 |
|---|---|---|
| **PDF 直読** | `check.mjs` **1 本だけ** | `pdftotext`(較正済み **Xpdf 4.00**) |
| **`.txt` 経由** | **他の 7 本** —— `hedge-index` / `bibmap` / `cycle-probe` / `full-graph` / `frdi-progress` / `genell-progress` / `paper-items` | **PyMuPDF 1.27.2(`fitz`)の `page.get_text()`** + 合字正規化(`ﬁ`→`fi`)+ `===== [page N] =====` の包み |

★**`.txt` は `pdftotext` 製ではない。** `paper-items.mjs` が書いていた契約
「各自 `pdftotext -layout` で作る」は**誤り**だった(2026-09-06、メタ第 6 回の実測)。
FrdI p.25 を Xpdf 4.00 で抽出すると `φ` 0 / `∈` 0 / `→` 0 個だが、`.txt` の同じ頁には
`φ` 17 / `∈` 6 / `→` 10 ある。6 通り(Xpdf/poppler × 既定/`-layout`/`-raw`)どれも不一致で、
PyMuPDF が一致した(標本 458 頁中 244 頁が完全一致、31 本は標本 5 頁すべて一致)。

★**したがって `CLAUDE.md` の「着手前に `hedge-index` で数える」の分母は PyMuPDF 産である。**
`check.mjs` の逐語照合が見ている文字列とは別物なので、
「`check.mjs` が通ったから `.txt` も正しい」とは**言えない**。

★★★**撤回**: 2026-09-05 のメタ第 3 回が書いた「`.txt` は Xpdf 風 111 本 / poppler 風 3 本」は
**原典 PDF 側の符号位置を見ていただけ**で、抽出実装の指紋ではなかった。数字も再現しない。
実際の混ざり方は**世代**である —— `page-marker` 式 114 本(PyMuPDF 世代)/ `formfeed` 式 22 本。
★**UTF-8 として復号できない `.txt` が 21 本**あり、formfeed 式とほぼ一致する
(EGA1 に U+FFFD 15,758 個、EGA2 14,206、BC 1,520)。

## `check.mjs` 側の話(ここは正しいまま。第 1024 で対策済み)

この機械には `pdftotext` が 2 つある(2026-09-05 実測):

| | 実装 | 場所 |
|---|---|---|
| Git Bash の PATH 先頭 | **Xpdf 4.00**(較正済み) | `C:\Program Files\Git\mingw64\bin\pdftotext.exe` |
| PowerShell の PATH 先頭 | poppler 25.07.0 | WinGet の Poppler |

2,157 頁のうち 1,718 頁(79.6%)で違うテキストが出る。最多は `´etale`(Xpdf・分解)対
`étale`(poppler・合成)。`Ŝ` は poppler では `U+0002`+`S`。

**★本当の罠はシェルではなくキャッシュだった。** `.cache/pdf-pages.json` の鍵が
`パス#頁#mtime#size` で実装を持っていなかったため、**正しいシェルで正しい実装を使っていても
poppler 産のキャッシュが残っていれば NG 175 件**になり、しかも 2 秒で返るので
「速い＝正常」に見えた。第 1024 で (1) `-v` から実装を同定(★poppler は `-v` を **stderr** に
終了コード 0 で出すので `spawnSync` が要る)、(2) PATH に無い既知の設置場所も含めて**版で選ぶ**、
(3) **キャッシュの鍵に実装を入れる**、の 3 点で塞いだ。

**How to apply:**

- 「原文にこう書いてある」と言う前に、**どちらの経路で読んだか**を意識する。
- `.txt` を作り直すときは **PyMuPDF** を使う。`pdftotext` で作り直すと
  7 本の道具の分母が黙って変わる。★再現する道具は**まだ無い**
  (`===== [page` を grep すると読む側 7 本・書く側 0 本。手順が外部にしか無い)。
- ★`0_Source/*.txt` の壊れ(孤立括弧・UTF-8 不能)を測る道具は
  メタ第 6 回が `tools/source-health.mjs` として書いたが、**取り込みは保留中**
  (`decisions-pending.md` の D12)。[[mathlib-index-nonascii-truncation]] と同じ
  「測定器そのものが壊れていた」系である。
