import ABC3.Found.PGC.WildBreakUpperBound
import ABC3.Found.PGC.ConcreteNormedModel
import ABC3.Found.PGC.PureStepSetup

/-!
# [pGC] ★**仮説ゼロの `k = 0` 出口を具体的な体に代入する** —— `ℚ₂(√2)/ℚ₂`

`Found/PGC/WildBreakUpperBound.lean` の
`WildBreakUpper.exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_free`
(★跳びについての仮説は 0 本、残るのは構造の 6 本)を
`Found/PGC/ConcreteNormedModel.lean` の `ℚ₂(√2)/ℚ₂` に代入し、
★**仮説を 1 つも残さずに** `∃ y ∈ K, ‖x − y‖ ≤ axDecay 2 1 · ‖g₂x − x‖` を得る。

★これは「跳びの下界・上界を定理にした」ことが**空虚でない**ことの証拠である
(下界 `1 ≤ t` は `WildBreakLowerBound` §4、上界 `(p−1)t ≤ p·e` は `WildBreakUpperBound` §3)。

## ★★`hiso`(等長)の供給について —— 測って分かったこと 2 つ

持ち場は「`PureStepSetup.lean:283` の `norm_algEquiv_eq` が供給する」と示唆していた。
★**測った結果、この道は `ℚ₂(√2)` には**そのままでは**通らない。**

1. ★`norm_algEquiv_eq` は `[NormedAlgebra k M]` を要求するが、
   ★**`NormedAlgebra ℚ_[2] M2` は instance になっていない**(`infer_instance` が失敗する)。
   `ConcreteNormedModel` が置いているのは `normedFieldM2 : NormedField M2` だけである。
   mathlib の `spectralNorm.normedAlgebra` は `def` であって `instance` ではなく、
   しかも `@NormedAlgebra K L _ (seminormedRing K L)` という**別のノルム構造**に対する形をしている。
2. ★**もっと短い道が在った**: `M2` のノルムは**定義そのものが**スペクトルノルムなので

       example (z : M2) : ‖z‖ = spectralNorm ℚ_[2] M2 z := rfl      -- ★通る

   が `rfl` で通る。したがって `hiso` は mathlib の
   `spectralNorm_eq_of_equiv (σ : L ≃ₐ[K] L) (y : L) : spectralNorm K L y = spectralNorm K L (σ y)`
   ★**1 本で済む**(向きは `.symm`)。

★**名前の衝突も測った**: `norm_algEquiv_eq` は木に**2 つ**ある ——
`PureStepSetup`(`{k M}` 版)と `AdjoinIntegers.lean:1821`(`(K : PAdicLocalField p)` 版)。
`open ABC3.Found.PGC` の下では**後者が解決される**ので

    error: Invalid argument name `k` for function `norm_algEquiv_eq`
    Hint: Perhaps you meant one of the following parameter names: `p` `K` `σ` `x`

になる。★`(k := ℚ_[2])` を付けても直らない(そもそも別の定理)。
-/

namespace ABC3.Found.PGC

namespace ConcreteDegPFree

open ConcreteNormedModel GainedTowerModel.GaloisTower

/-- `M2` のノルムは**定義により**スペクトルノルム。★`rfl` で通る。 -/
theorem norm_eq_spectralNorm_M2 (z : M2) : ‖z‖ = spectralNorm ℚ_[2] M2 z := rfl

/-- ★★**`hiso`(等長)は無料** —— `spectralNorm_eq_of_equiv` 1 本。

★`PureStepSetup.norm_algEquiv_eq` は `NormedAlgebra ℚ_[2] M2` を要求し、
それは instance になっていないので使えない(冒頭 docstring)。 -/
theorem hiso_g2 (z : M2) : ‖g2 z‖ = ‖z‖ := (spectralNorm_eq_of_equiv g2 z).symm

/-- ★★★★★**仮説ゼロの出口が具体的な体で成り立つ** —— `ℚ₂(√2)/ℚ₂`。

★跳び(`t`)についての仮定も、`harith` も、`hbreak` も**一切渡していない**。 -/
theorem model_deg_p_free_exists (x : M2) :
    ∃ y : twr g2 2 0, ‖x - algebraMap (twr g2 2 0) M2 y‖ ≤ axDecay 2 1 * ‖g2 x - x‖ :=
  WildBreakUpper.exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_free
    (p := 2) (π := pi2) g2 orderOf_g2 hiso_g2 (by simpa using finrank_M2) valK
    norm_two_M2 norm_pi2_lt_one x

/-! ## `.src` と 公理 -/

def model_deg_p_free_exists.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_eq_spectralNorm_M2
#print axioms hiso_g2
#print axioms model_deg_p_free_exists

end ConcreteDegPFree

end ABC3.Found.PGC
