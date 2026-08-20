---
name: node-e-eats-backticks
description: Bash ツールの node -e "..." はバックティックをコマンド置換として実行する。識別子が黙って消える。
metadata:
  type: feedback
---

`node -e "..."`(二重引用符)を Bash ツールから走らせると、**バックティックが
コマンド置換として実行される**。Bash ツールは `bash -c '...'` で包むが、
その内側の二重引用符では `` ` `` が生きているためである。

★実害(2026-08-20): 依存グラフの note に `` `mprec_effSub_iff` `` のような
識別子を書いたら、**識別子がまるごと消えた JSON がコミットされた**。
bash は `mprec_effSub_iff: command not found` を吐いていたが、
`node` 自体は成功して終了するので**ゲートも通ってしまう**。

**Why:** JSON の note は「なぜそう決めたか」の唯一の記録で、識別子が消えると
後から辿れない。しかも失敗が exit code に出ない(静かに壊れる)。

**How to apply:** 識別子をバックティックで囲む文章を JSON/Markdown に書くときは、

* クォート付き heredoc の Python(`<<'PYEOF'` … `PYEOF`)で `json.load` / `json.dumps` する、または
* Write ツールで `.mjs` を書いてから `node <file>` で走らせる。

★同じ理由で `$(...)` と `$VAR` も展開される。
★[[heredoc-eats-backslash]] は別の罠(バックスラッシュ)で、こちらは
クォート付き heredoc でも起きる。Python 側の文字列は `\uXXXX` で書くと安全。
★書いたあとは必ず読み返して確かめること([[verify-insertion-not-just-ok]])。
