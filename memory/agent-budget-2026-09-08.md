---
name: agent-budget-2026-09-08
description: 2026-09-08 のユーザー判断 —— サブエージェントは lean-prover のみ。meta-optimizer / math-planner / general-purpose は廃止。
metadata:
  type: feedback
---

2026-09-08、ユーザーの判断で**サブエージェントは `lean-prover` だけ**に絞った。
**`meta-optimizer` / `math-planner` / `general-purpose` は廃止。配ってはならない。**
`lean-search` / `Explore` は安いので必要なら短く使ってよい。

**Why:** トークン消費が Max プランの週間上限に対して速すぎる（ユーザー計測で
「約 15 時間で週間上限の 30%、1.5 日で使い切る速度」）。
実測（通知 238 件、`node tools/agent-timing.mjs --limit 0`）ではエージェント総トークン
46,658,603 のうち lean-prover 61.2% / meta-optimizer 15.0% / general-purpose 10.1% /
math-planner 7.7% で、**廃止した 3 種で 32.8%**。
さらにユーザーの指摘どおり、meta-optimizer 36 回の実態は「改善」ではなく
**「測定器の修理と本体の数字の訂正」= 検証**だった。検証は外注する仕事ではない。
速度改善の証拠も無い（1 ノードあたり中央 27.3 分 → 31.2 分、172k → 243k と**上がっている**）。
ユーザーは「進捗速度も 8 月から変わっている実感がない」とも述べており、私の測定と矛盾しない。

**How to apply:**
- 検証・測定・自分の数字の検算は**本体が inline で**やる。数字を報告に出す前に必ず期間を切って数える
  （[[measurement-numbers-must-be-measured]]）。
- 方針立ても本体が原典を直接読む。
- **持ち場（brief）は従来の 1/3 の長さに削る。** 貼るのは道具の実パスと禁止事項だけで、
  数学の見込みは要らない（本体の見込みは 6 波連続で外れ、実装者は毎回「その 1 つ手前の部品」で通した）。
- `ResearchPaper/decisions-pending.md` の記録も短くする。残すのは
  ①真偽の判定 ②在庫の測定 ③次の 1 点 ④VERDICT/COST の 4 つだけ。
- 規約本体は `ResearchPaper/autonomy-policy.md` の「3.5 エージェント構成」に書いた。
