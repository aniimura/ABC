---
name: orchestration-graph-first
description: 依存グラフを外部の真実とし、agent には持ち場だけを渡す。★2026-09-08 に役割は lean-prover のみに絞られた。
metadata:
  type: project
---

2026-09-05 に「1 セッションが順番に作業する」形をやめ、**依存グラフを外部の真実**とし、
agent には**持ち場(brief)だけ**を渡す形にした。

- 「次にどこを叩くか」: `node tools/frontier.mjs`（`startable` かつ `downstream` の大きいノードから）
- 「持ち場を切り出す」: `node tools/brief.mjs --node <rel>`
- 役割: `.claude/agents/` に **lean-prover / lean-search / lean-verifier** の 3 本

★**2026-09-08 の変更**: `math-planner` と `meta-optimizer` は**廃止**した（[[agent-budget-2026-09-08]]）。
方針立て・検証・測定は**本体が inline でやる**。同時起動は **1 波 1 体**
（旧 [[orchestration-parallel-cap-5]] の上限 5 も、その後の上限 2 も無効）。

**Why:** 形式化は一本道にならず、作業中に依存 DAG が成長する。全体を 1 つのコンテキストに
入れることはできず、単独ノードだけを渡すと「なぜこの補題が要るのか分からない」になる。

**How to apply:** 持ち場は**短く**（実パスと禁止事項だけ。数学の見込みは書かない ——
本体の名指しは 6 波連続で外れ、実装者は毎回「その 1 つ手前の部品」で通した）。
docstring は渡す（過去の監査の結論は道具の在庫に依存するので、結論だけでは同じ誤りを繰り返す）。
