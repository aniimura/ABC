---
name: orchestration-parallel-cap-5
description: 並列 agent 数の上限は 5 程度。増やしすぎない(2026-09-05 ユーザー指示)
metadata:
  type: feedback
---

サブエージェントの同時起動は **最大 5 程度**に抑える。増やしすぎない。

**Why:** 2026-09-05、外部依存グラフ中心のオーケストレーションを準備した直後に
ユーザーが明示。この木の前線は `sorry` ノード 24 / うち着手可能 13 しかなく、
着手可能でないノードに agent を割いても上流の `sorry` に当たって止まる。
並列度を上げてもスループットは上がらず、共有 worktree(`.git` と作業ツリーを
ABC3b/c と共有)への同時書き込みでビルドが不安定になる副作用だけが増える
——同日に `lake build` の過渡的失敗を 3 回観測している。

**How to apply:** `tools/frontier.mjs` の既定 `--limit` は 5。
1 波で配る持ち場は 5 件まで。足りなければ**次の波**にする。
Workflow を使う場合も同時 agent 数を 5 以下に保つ。
並列度は結果であって目標ではない([[orchestration-graph-first]])。
