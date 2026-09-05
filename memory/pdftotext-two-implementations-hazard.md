---
name: pdftotext-two-implementations-hazard
description: この機械には pdftotext が 2 つある。Git Bash=Xpdf、PowerShell=poppler。check.mjs は必ず Git Bash から
metadata:
  type: feedback
---

**`node tools/check.mjs` は必ず Git Bash(Bash ツール)から走らせること。PowerShell から走らせてはならない。**

**Why:** この機械には `pdftotext` が 2 つ入っている(2026-09-05 実測):

| シェル | 実装 |
|---|---|
| Git Bash | **Xpdf 4.00** |
| PowerShell | **poppler 25.07.0** |

出力が違うため、**同じ木で NG 13 件(Bash)対 176 件(PowerShell)**になる。
差 163 件はすべて S4 逐語照合と引用照合。

いま鳴っていないのは `.cache/pdf-pages.json` に Xpdf の出力が入っているからで、
★**キャッシュの鍵に `check.mjs` 自身のハッシュが含まれる**ので、
**`check.mjs` を編集した次の実行でキャッシュは必ず捨てられ、作り直される**。
そのとき PowerShell から走らせると NG が 176 に跳ねる。

**How to apply:** `check.mjs` を触った直後は特に注意。作り直しは Bash から。
メタ agent がこれを踏んで 1 時間使った(`ResearchPaper/meta-backlog.md` M8)。
隔離 worktree はさらに `core.autocrlf=true` で CRLF 検出され逐語照合がバイトで
食い違うので、worktree の NG 件数を本体と直接比べないこと(同 M10)。
