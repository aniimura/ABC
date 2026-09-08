---
name: lean-check-via-mcp-repl
description: Lean の検査は node tools/leanfile.mjs（8〜13 秒）。★MCP lean_check は 2026-09-08 に使用禁止になった。確定は必ず build.mjs。
metadata:
  type: feedback
---

★**2026-09-08 訂正: MCP（`abc3-lean` の `lean_check`）は使わない。**
理由は速度ではなく事故 —— ★**基準環境がエラー文を 1 つも出さずにすり替わる**。
実測で 11 件 / agent 9 体、うち 9 件で「同じ imports を頼んだ別の agent が同時に走っていた」と
相手を名指しできた。同じ木に「REPL は処理中」が 19 件。
⇒ ★**`node tools/leanfile.mjs <file>`（8〜13 秒/往復）を使う。**olean を書かないので並行して安全。

**形式化の律速は検査である**（この観察自体は今も有効、2026-08-17 実測）:
`lake env lean` は**ファイル全体を再検査する**ので、ファイルが育つほど重くなる
（3,804 行で数分 / 322 行で 9 秒）。★だから**ファイルを小さく保ち、抽象核を別宣言に切り出す**
（[[lean-build-check-discipline]]）。

**How to apply:**
1. 書きかけは `leanfile.mjs` で通す。既定は診断ブロックだけ（上限 60 行）、全文は `.cache/leanfile-*.log`。
   ★後から `--errors` / `--grep` / `--full` で切り出す（`lean` を呼ばない）。
2. ★**節目では必ず `node tools/build.mjs <target>`。** 書いた順序・`variable` の効き方・
   リンタは build でしか出ない。「1 ファイル通った」を「木が通った」と読まない。
3. コミット前にゲート `node tools/check.mjs --brief`。

関連: [[lean-build-check-discipline]] [[agent-budget-2026-09-08]]
