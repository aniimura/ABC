---
name: orchestration-graph-first
description: 進め方は「外部依存グラフが真実、agent には持ち場だけ」。設計は ResearchPaper/orchestration.md
metadata:
  type: project
---

2026-09-05 に進め方を切り替えた。「1 セッションが順番に作業する」形をやめ、
**依存グラフを外部の真実**とし、agent には**持ち場(brief)だけ**を渡す。

- 設計: `ResearchPaper/orchestration.md`
- 「次にどこを叩くか」: `node tools/frontier.mjs`
  (`sorry` ノードごとに downstream / dsItems / blockers / startable を計算)
- 「持ち場を切り出す」: `node tools/brief.mjs --node <rel>`
  (sorry の宣言・docstring・`.src`・`.needs`・依存/被依存・許された手順)
- 役割: `.claude/agents/`(lean-search / math-planner / lean-prover /
  lean-verifier / meta-optimizer)

**Why:** 形式化は一本道にならず、作業中に依存 DAG が成長する。全体を 1 つの
コンテキストに入れることはできず、単独ノードだけを渡すと「なぜこの補題が
要るのか分からない」になる。全体像を Claude の外に置くのが答え。

**How to apply:** 着手前に必ず `frontier.mjs` を引き、`startable` かつ
`downstream` の大きいノードから配る。`startable` でないノードに agent を
割かない(上流の `sorry` で止まる)。同時起動は [[orchestration-parallel-cap-5]]。
docstring を必ず渡すこと——過去の監査の結論は**道具の在庫に依存する**ので、
結論だけでなく経緯ごと渡さないと同じ誤りを繰り返す。
