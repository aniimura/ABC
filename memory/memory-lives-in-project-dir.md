---
name: memory-lives-in-project-dir
description: このプロジェクトのメモリは D:\Math_ABC3\memory\ に置く(Claude の既定の隠しフォルダではない)。
metadata:
  type: feedback
---

メモリファイルは `D:\Math_ABC3\memory\` に保存する。`C:\Users\Aruta\.claude\projects\D--Math-ABC3\memory\` には、ここを指すポインタ1行だけを置く(ユーザー指示、2026-08-14)。

**Why:** プロジェクトの前提・判断の記録が、隠しフォルダではなくプロジェクト本体と同じ場所にあれば、ユーザーが直接読め、git でバージョン管理でき、フォルダごと移動・複製しても失われない。本プロジェクトはメモリの内容(北極星・争点の所在・過去世代の扱い)が計画そのものと不可分なので、`PLAN.md` と同じ場所にあるのが自然。

**How to apply:** 新しいメモリは `D:\Math_ABC3\memory\<slug>.md` に書き、索引は `D:\Math_ABC3\memory\MEMORY.md` に1行追加する。既定の隠しフォルダ側には書かない(ポインタ1行のみ維持)。セッション開始時に自動で読み込まれるのはポインタだけなので、**プロジェクトの作業を始める前に `D:\Math_ABC3\memory\MEMORY.md` を読むこと**。
