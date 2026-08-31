/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Interface.GenEll.EllModuli

/-!
# 界面の測定と訂正 —— **`imageContainsSL2_of_torsionExt` は `5 ≤ l` を要る**（`Check`）

**これは原典の主張ではない**（我々の界面についての事実）ので `.src` を持たない。

## ★★★★★★★★2026-08-31 の測定と訂正（第 776）

`Interface/GenEll/EllModuli.lean` の `imageContainsSL2_of_torsionExt` は、以前

```
imageContainsSL2_of_torsionExt : ∀ (E : Curve) (l : ℕ), Nat.Prime l →
  HasMultRed (torsionExt E) → PrimeToLocalHeights (torsionExt E) l →
  ¬ HasLCyclic (torsionExt E) l → ImageContainsSL2 E l
```

であった。★★**`l ∈ {2, 3}` では witness が作れない。**

## ★★★なぜか

原文 p.20 の最終段は **`Lemma 3.1, (iv)`** を使う。その形式化
`Found/GenEll/Sl2Padic.lean` の `lemma_3_1_iv` は仮説に **`5 ≤ l`** を持つ
（原文の `Lemma 3.1, (iv)` も `l ≥ 5` を仮定している——`pro-l` 群の Frattini 型の
議論が `l ≥ 5` を要求する）。

★したがって `l = 2, 3` では「`α` が入り安定直線が無い」から `SL₂(ℤ_l) ⊆ Im` は**出ない**。

## ★★★★訂正（第 776）

欄に `5 ≤ l` を加えた。★★**`Theorem 3.8` の statement は変わらない**——
その証明の中で `5 ≤ l` が**両方の条件から出る**からである:

* 条件 (a): `23040·100·d·(ht^Falt + C·d^ϵ) ≤ l` で、`C` の取り方から括弧の中は `≥ 1`、
  `d ≥ 1` なので `l ≥ 2304000 ≥ 5`
* 条件 (b): `l` は `30` と素な素数だから `l ∉ {2, 3, 5}`、すなわち `l ≥ 7`

## ☆同じ形の測定

* `Check/GenEll/EllModuliDegInfPos.lean`（第 745）——界面は `deg∞ > 0` を強制する
* `Check/GenEll/LcyclicExcTooStrong.lean`（第 754-755）——`mem_lcyclicExc` は `l` の下界を落としていた
* 本ファイル（第 776）——`imageContainsSL2_of_torsionExt` は `5 ≤ l` を落としていた

★どれも **witness を実際に作ろうとして初めて見えた**ものである。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll

/-- ★★★★★★**訂正後に界面が主張していること**——`5 ≤ l` が仮説に入っている。 -/
theorem imageContainsSL2_needs_five (D : EllModuliData) (E : D.Curve) (l : ℕ)
    (hl : Nat.Prime l) (hl5 : 5 ≤ l)
    (hm : D.HasMultRed (D.torsionExt E)) (hp : D.PrimeToLocalHeights (D.torsionExt E) l)
    (hc : ¬ D.HasLCyclic (D.torsionExt E) l) : D.ImageContainsSL2 E l :=
  D.imageContainsSL2_of_torsionExt E l hl hl5 hm hp hc

end ABC3.Check.GenEll
