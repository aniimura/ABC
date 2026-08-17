---
name: lean-check-via-mcp-repl
description: Lean の検査は abc3-lean MCP の lean_check（0.01 秒）で 1 宣言ずつ。確定は必ず lake build。
metadata:
  type: feedback
---

**形式化の律速は検査**である。`lake env lean` は**ファイル全体を再検査する**ので、
ファイルが育つほど 1 往復が重くなる（2026-08-17 実測）:

| 検査のしかた | 1 往復 |
|---|---|
| `lake env lean` on 3,804 行 | **数分** |
| `lake env lean` on 322 行 | **9 秒** |
| `abc3-lean` MCP の `lean_check` | **0.01〜0.02 秒** |

**How to apply:**

1. 書きかけの定理は `lean_check` で 1 個ずつ通す。`lake env lean` でファイル全体を回さない。
2. `lean_start(imports)` はセッション中 1 回（冷 90 秒／温 8 秒）。
   **`olean` を作り直したら `lean_start` を呼び直す**。
3. ★★**節目では必ず `lake build`。** REPL は `olean` しか読まないので、
   **書いた順序・`variable` の効き方・リンタは `lake build` でしか出ない**。
   `lean_check` が通ったことを「ファイルが通った」と読まないこと。
4. コミット前にゲート `node tools/check.mjs`。

**Why:** 今日 1 日で `lake env lean` を 30 回以上回し、その大半が 3,800 行の再検査だった。
同じ編集の検査に 20〜40 倍の差が出ていた。

**罠（実装時に踏んだ）:** Windows では `repl.exe` を直に spawn してはならない。
`libleanshared.dll` の解決にツールチェインの `PATH` が要る。直起動は無応答のまま
600 秒でタイムアウトした。`lake env <repl>` 経由で起動すること。

道具: `tools/mcp-lean/server.mjs`（依存なし）＋ `tools/lean-repl/`（`repl` の `v4.31.0`、gitignore）。

関連: [[lean-build-check-discipline]] [[leaf-first-with-graph-feedback]]
