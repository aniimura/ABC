---
name: abc3-independent-of-prior-generations
description: Math_ABC3 は D:\Math_ABC / D:\Math_ABC2 とは無関係の独立プロジェクト。過去世代の記録を根拠に使わない。
metadata: 
  node_type: memory
  type: project
  originSessionId: fbb27ed5-3b09-4626-a97a-bdf85c79f517
  modified: 2026-08-14T01:53:09.000Z
---

`D:\Math_ABC3` は `D:\Math_ABC`(第1〜2期)・`D:\Math_ABC2`(第3期)とは**無関係**の独立プロジェクトとして扱う(ユーザー明言、2026-08-14)。過去世代の設計・失敗記録(`KNOWLEDGE.md`・`PIPELINE_REVIEW.md`・9工程パイプライン・JSON台帳等)を、ABC3の計画の根拠として引用しない。

**Why:** ABC3 は `ResearchPaper/`(0_Source・1_Structured)と `Donburi_v2.0.exe` を ABC2 から複製して作られているため、**残存ファイルが過去世代の作業体系を参照している**。実際、`1_Structured` の全21ファイルが存在しないフォルダ(`2_LocatorMap` 38箇所・`MAP-` 103箇所・`4_CalibrationPlan` 13箇所・`IUT_DependencyMap` 12箇所)を参照しており、これを手がかりに `D:\Math_ABC2` を読んで計画を組み立て、ユーザーから訂正を受けた(`1_Structured/README.md` はその後ユーザーが削除)。

**How to apply:** 残存する過去世代への参照を見つけても、`D:\Math_ABC`・`D:\Math_ABC2` を読みに行かない。ABC3 の計画の根拠は、(a)`idea.md`・`idea2.md`、(b)`ResearchPaper/0_Source` の原文(PDFが権威)、(c)mathlib の実測、(d)Lean で実際に実行して確かめた結果、に限る。現在の計画は [[abc3-plan-two-track]] を参照。
