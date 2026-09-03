---
name: corrhyp-track-goal
description: CorrHyp(Correspondences on Hyperbolic Curves)のLean形式化ゴールを2026-09-04に節ごとに設定した。
metadata:
  type: project
---

`ResearchPaper/corrhyp-goal.md` に CorrHyp(Mochizuki, *Correspondences on Hyperbolic Curves*, 18頁)の
節ごとの形式化ゴールを設定した: `§1 0/5, §2 0/6, §3 0/3, §4 0/2, §5 0/7, §6 0/1`(合計 0/24、Skeleton 未着手)。
§0(Theorem A/B/C)は本文中で後続定理の再掲と明記されているため対象外——[[genell-track-b]] の §1 起点と同じ扱い。

**Why:** 論文は `ResearchPaper/papers.json` に登記済みだが notationRisk が unmeasured(誰も目視していない)。
着手(Skeleton化)の前に p.3-4 の 260dpi 目視が要る。

**How to apply:** 次の一手は §1(Def 1.1-1.5、依存0の葉)から Skeleton を立てること。
まだどのセッションもこのトラックの担当を宣言していない([[genell-track-b]]=ABC3b が GenEll、
otricomm/Prop44 チェーンを別セッション ABC3c が担当中、と 2026-09-04 時点で確認済み)。
着手する前に他セッションと重複がないか `ListAgents`/直接メッセージで確認すること。
