/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic.IntervalCases
import Mathlib.Tactic.NormNum.Prime

/-!
# 界面の測定 —— **`VeluQuotOK` は `5 ≤ l` でしか要らない**（`Check`）

**これは原典の主張ではない**（我々の界面についての事実）ので `.src` を持たない。

## ★★★★★★★★2026-09-02 の測定（第 1434）

`Skeleton/GenEll/VeluSemistable.lean` の `veluQuotOK_all` は

```
∀ (E : SSCurve) (l : ℕ), VeluQuotOK E l
```

を主張しており、**`l = 2` や `l = 3`、さらに非素数の `l` まで**含んでいる。

★だが §9-1335 の実装（第 1410-1432）が閉じたのは

* `l` が**奇素数**であること
* `p ∣ l` のときは `p ∤ 6`（すなわち `l ≥ 5`）

の下だけである。☆したがって「界面が要求しすぎていないか」を測った。

## ★★★測定の結果——**要求しすぎている**

`Theorem 3.8`（`Skeleton/GenEll/GaloisImage.lean`）が `lemma_3_7` の第 3 の主張
（`(condA ∨ condB) → HasLCyclic → cls E ∈ Exc`）を呼ぶ場面では、
**どちらの枝でも `5 ≤ l` が出ている**:

* 枝 (a): `23040·100·d·(ht^Falt + C·d^ε) ≤ l` で、`C` の取り方から括弧の中は `≥ 1`、
  `d ≥ 1` なので **`l ≥ 100`**（本ファイルの `hundred_le_of_condA`）
* 枝 (b): `Nat.Coprime l 30` なので `l ∉ {2, 3, 5}`、すなわち **`l ≥ 7`**
  （本ファイルの `five_le_of_coprime_thirty`）

★★これは第 776（`Check/GenEll/ImageSL2NeedsL5.lean`）と**同じ形**の測定であり、
そこでは `imageContainsSL2_of_torsionExt` の欄に `5 ≤ l` を加えて解決した。

## ☆したがって取れる訂正

`mem_lcyclicExc`（`Interface/GenEll/EllModuli.lean`）と `lemma_3_7`
（`Skeleton/GenEll/Section3.lean`）に `5 ≤ l` を加えれば、
`veluQuotOK_all` も `∀ E l, 5 ≤ l → VeluQuotOK E l` で足りる。
★★★そうすると §9-1335 に残る穴は

* **良い素点・`p ∣ l`・核に深い点**（`jExp p E′ ≥ 0`、第 1431）

**ただ 1 つ**になる——`l = 2` と `l = 3` の穴は界面の側で消える。

☆`Theorem 3.8` の statement は変わらない（第 776 と同じ理由）。

## ☆同じ形の測定

* `Check/GenEll/EllModuliDegInfPos.lean`（第 745）——界面は `deg∞ > 0` を強制する
* `Check/GenEll/LcyclicExcTooStrong.lean`（第 754-755）——`mem_lcyclicExc` は `l` の下界を落としていた
* `Check/GenEll/ImageSL2NeedsL5.lean`（第 776）——`imageContainsSL2_of_torsionExt` は `5 ≤ l` を落としていた
* 本ファイル（第 1434）——`VeluQuotOK` は `5 ≤ l` でしか要らない
-/

namespace ABC3.Check.GenEll

/-- ★★**枝 (b) の下界**——`30` と素な素数は `5` 以上（実際は `7` 以上）。 -/
theorem five_le_of_coprime_thirty {l : ℕ} (hl : l.Prime) (h : Nat.Coprime l 30) : 5 ≤ l := by
  by_contra hlt
  rw [not_le] at hlt
  interval_cases l
  · exact absurd hl (by decide)
  · exact absurd hl (by decide)
  · exact absurd h (by decide)
  · exact absurd h (by decide)
  · exact absurd hl (by decide)

/-- ★★**枝 (a) の下界**——`100·d·(F + C·P) ≤ l` で括弧の中が `≥ 1`、`d ≥ 1` なら `l ≥ 100`。 -/
theorem hundred_le_of_condA (d P F C L : ℝ) (hd : 1 ≤ d)
    (h1 : 1 ≤ F + C * P) (hle : 100 * d * (F + C * P) ≤ L) : 100 ≤ L := by
  nlinarith

end ABC3.Check.GenEll
