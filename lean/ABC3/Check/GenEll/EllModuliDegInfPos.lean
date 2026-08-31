/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Interface.GenEll.EllModuli

/-!
# 界面の測定 —— **`EllModuliData` は `deg∞ > 0` を強制する**（`Check`）

**これは原典の主張ではない**（我々の界面についての事実）ので `.src` を持たない。

## ★★★★★★★★2026-08-31 の測定（第 745）

`Interface/GenEll/EllModuli.lean` の `EllModuliData` は §4 の帳簿として

| 欄 | 内容 |
|---|---|
| `multCard_pos` | `0 < multCard E` |
| `multPrime_prime` | `multPrime E j` は素数 |
| `localHt_pos` | `0 < localHt E j` |
| `sum_localHt_eq` | `∑ h_j·log(p_j) = 23040·d·deg∞([E_L])` |

を持つ。★★この 4 本から**どの `E : Curve` に対しても `deg∞(cls E) > 0`** が出る:

* `multCard E ≥ 1` なので和は空でない
* `localHt E j ≥ 1`、`multPrime E j ≥ 2` なので各項は `≥ log 2 > 0`
* したがって左辺 `> 0`、`d > 0` なので `deg∞ > 0`

## ★★★★★これが witness の設計に効くこと

**`Curve` 欄に「至る所良還元」の楕円曲線を入れてはならない**。
（数体上には至る所良還元の楕円曲線が実在する——例えば `ℚ(√29)` 上。
`ℚ` 上には無い（Tate）が、界面は一般の数体を許している。）

★これは界面の欠陥ではなく、原文の題

> Corollary 4.3. (Full Galois Actions for **Degenerating** Elliptic Curves)

そのものである。★★`Found/GenEll/EllModuliObjects.lean` の `DegCurve`
（半安定＋乗法還元を持つ）がその形である。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll

/-- ★★★★★★★★**界面は `deg∞ > 0` を強制する**。

★§4 の帳簿の 4 本（`multCard_pos`・`multPrime_prime`・`localHt_pos`・`sum_localHt_eq`）
だけから出る。 -/
theorem degInf_pos_of_ellModuliData (D : EllModuliData) (E : D.Curve) :
    0 < D.degInf (D.cls E) := by
  have hcard := D.multCard_pos E
  have hsum := D.sum_localHt_eq E
  have hpos : (0:ℝ) < ∑ j : Fin (D.multCard E),
      (D.localHt E j : ℝ) * Real.log (D.multPrime E j) := by
    refine Finset.sum_pos (fun j _ => ?_) ⟨⟨0, hcard⟩, Finset.mem_univ _⟩
    have h1 : (1:ℝ) ≤ (D.localHt E j : ℝ) := by exact_mod_cast D.localHt_pos E j
    have h2 : (2:ℝ) ≤ (D.multPrime E j : ℝ) := by
      exact_mod_cast (D.multPrime_prime E j).two_le
    exact mul_pos (by linarith) (Real.log_pos (by linarith))
  rw [hsum] at hpos
  have hd : (0:ℝ) < (D.degOfDefinition E : ℝ) := by exact_mod_cast D.degOfDefinition_pos E
  by_contra hc
  push_neg at hc
  nlinarith

end ABC3.Check.GenEll
