---
name: frdi-s1-two-blockers
description: [FrdI] §1 は 13/15 で止まっている。残る 2 件は実装不足でなく原文の仮定から出ない。捻れ積が反例の道。
metadata:
  type: project
---

`[FrdI] §1` は **13/15**（2026-08-16 時点）。残る 2 件は**実装が足りないのではなく、
原文の仮定から結論が出ない**ことが式で特定されている。`lean/ABC3/Gap/FrdI/Section1.lean`
に `GapRecord`（分類 ② `missingMath`）として記録済み。

**1. `Proposition 1.14, (iii)` の `⟸`** — 不足は「prime-Frobenius 射が FSM 射」。
`Definition 1.3` で mono を与える条項は `preStepMono` **1 つだけ**で pre-step にしか効かない。
穴の正体は `Found/FrdI/Prop114.lean` の `frobType_cancel_invariants`（Frobenius 型で消去すると
不変量は全部一致する）と `mono_of_frobType_of_faithful`（`𝒞 → 𝔽_Φ` が忠実なら mono）で式にした。
つまり**穴はちょうど単元の分だけ**開いている（`Definition 1.3` は `faithfulUpToUnits` しか与えない）。
迂回不可：(ii) の特徴づけが「pre-step でない」を含んで循環しており、(iii) がそれを切るので両向き要る。
使用箇所（FrdI p.63 `Theorem 3.4`）でも追加仮定は無い。

**2. `Proposition 1.6, (v)` の `base-trivial` / `metrically trivial` の `⟸`** — 不足は `Aut-ample`。
`𝒞'` の同型を作ると四角形が `f` の `Base` を 1 つに指定するが、base-trivial が与えるのは
**ある**同型であって底を指定した同型ではない。原文は (vi) を「if」（片向き）としか書いておらず、
著者は (v) と (vi) で向きを書き分けている。`⟹` 向きは両方とも実装済み。

**反例の道**：`lean/ABC3/Check/FrdI/TwistedFrobenioid.lean` に**捻れ積 `𝔽_Φ ⋉ G`** を構成した。
合成を Frobenius 次数で捻る（`(f,g) ≫ (f′,g′) = (f ≫ f′, g^{degFr f′}·g′)`）のが要点で、
単なる直積だと `Frobenius-normalized` が `G` を自明にしてしまう。
**圏であること・`PreFrobenioid` であること・次数 `d` の射が `G` の `d`-捻れで mono を失うこと**
（`not_mono_twist`）は証明済み。残るのは `Definition 1.3` の全条件（`FrobenioidCore` の約 21 条）。
検算では `preStepMono`・`faithfulUpToUnits`・圏同値（三角形が `G` 成分を一意に決めるので前順序のまま）
・`frobDegSurj` すべて通る見通し。

**判断**：`.src`（条なし = その原典項目を完全に実装したという主張）を追加仮定つきで付けることは
しない。[[challenger-audit-without-context]] の Prop 1.10 取り下げで定めた規律。
15/15 に到達する道は「捻れ積の検証を終えて原文の主張が偽だと確定させる」か
「規律を変える」かのどちらかで、後者は勝手にはやらない。
