---
name: frdi-s1-two-blockers
description: [FrdI] §1 は 13/15 で止まる。Prop 1.14 は主張が偽と反例で確定したので 14/15 が上限。Prop 1.6 (v) は Aut-ample 不足。
metadata:
  type: project
---

`[FrdI] §1` は **13/15**(2026-08-16 時点)。★**残る 2 件は実装不足ではない。**
`lean/ABC3/Gap/FrdI/Section1.lean` に `GapRecord` として記録済み。

## 1. `Proposition 1.14, (iii)` の `⟸` —— ★★**主張が偽であることを反例で確定した**

分類 ③ `sourceGap`。★**したがってこの項目に `.src`(完全実装の主張)を付ける道は無い。
§1 の上限は 14/15 である。**

反例は `lean/ABC3/Check/FrdI/TwistedFrobenioid.lean` の **捻れ積** `𝔽_Φ ⋉ G`:
合成を Frobenius 次数で捻る `(f,g) ≫ (f′,g′) = (f ≫ f′, g^{degFr f′}·g′)`。
単なる直積だと `Frobenius-normalized` が `G` を自明にしてしまうのが要点。
`twFrobenioidCore`(21 条)＋ `twIsFrobenioid` で `Definition 1.3` を**全条件**証明済み。

★具体化 `cx2P`: `𝒟 = Discrete PUnit`、`Φ = ℕ`(定数)、`G = ∏_n ℤ/n`。
- `G` に**すべての位数の捻れ**があるので、★**次数 > 1 の射は 1 本も mono でない**
- ⟹ FSMI 射はちょうど `Div = 1` ⟹ ★★**鎖の長さは `Div` そのもの**
- ⟹ `Div = 1` の irreducible な **pre-step** で条件(有界性)が**成り立ってしまう**
  (`cx2_refutes_1_14_iii` / `prop_1_14_iii_mpr_false`)

★原文 p.42 の (a) の議論は「次数を大きくした prime-Frobenius `ψ` を後置」だが、
条件文は `ψ` にも FSMI を要求している(p.41「where α1, . . . , αn, ψ are FSMI-morphisms」
—— 原文で確認済み)。この圏にはその `ψ` が無い。

★**残る唯一の反証経路は「我々の `BoundedFSMIFactor` の写し方が原文と違う」ことである。**

## 2. `Proposition 1.6, (v)` の `base-trivial` / `metrically trivial` の `⟸`

分類 ② `missingMath`。不足は `Aut-ample`。`𝒞'` の同型を作ると四角形が `f` の `Base` を
1 つに指定するが、base-trivial が与えるのは**ある**同型であって底を指定した同型ではない。
原文は (vi) を「if」(片向き)としか書いておらず、著者は (v) と (vi) で向きを書き分けている。
`⟹` 向きは両方とも実装済み(`cfp_metricallyTrivial_mp` / `cfp_baseTrivial_mp`)。
★**こちらはまだ反例を構成していない**ので ③ を名乗っていない。

## 規律

`.src`(条なし = その原典項目を完全に実装したという主張)を追加仮定つきで付けることは
しない([[challenger-audit-without-context]] の Prop 1.10 取り下げで定めた規律)。
★**主張が偽と分かった項目は、そもそも付けられない。**
