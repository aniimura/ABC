---
name: mathlib-local-instances
description: mathlib には「local instance でしか付かない」インスタンスがある——Ideal.Quotient.field はその代表
metadata:
  type: feedback
---

`Ideal.Quotient.field`（`I.IsMaximal` から `Field (R ⧸ I)`）は **mathlib の大域 instance ではない**。
mathlib 自身が使う所ごとに `attribute [local instance] Ideal.Quotient.field` と書いている
（`FieldTheory/Separable.lean`、`NumberTheory/KummerDedekind.lean` など）。

**Why:** 2026-08-21、`Algebra.IsSeparable (ℤ ⧸ span {p}) (𝓞 K ⧸ P)` を出そうとして
`Field (ℤ ⧸ span {p})` が synth できず 5 往復した。`haveI : Field _ := Ideal.Quotient.field _` と
手で入れると**別の instance になって diamond**が立ち、`PerfectField.ofFinite` が繋がらない。
★正解は `attribute [local instance] Ideal.Quotient.field` をファイルの先頭に書くこと。

**How to apply:**
* 剰余体まわり（`R ⧸ p` を体として使う）では、まずこの attribute を疑う。
* 一般に「mathlib のこのインスタンスがなぜか見つからない」ときは、
  **mathlib 側で `attribute [local instance]` されていないか grep する**：
  `grep -rn "attribute \[local instance\] <name>" .lake/packages/mathlib/`
* 手で `haveI` を入れて通ったように見えても、後段で defeq が壊れることがある。
  → [[namespace-shadows-mathlib]]（同じく「見つからない」系の罠）
