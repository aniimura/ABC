import ABC3.Found.GenEll.LogDiff
import ABC3.Found.GenEll.ProductFormula
import Mathlib.RingTheory.DedekindDomain.Factorization
import Mathlib.NumberTheory.NumberField.Discriminant.Different

/-!
# [GenEll] Definition 1.5, (iii) —— `log-diff` の値は判別式である(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> determines a well-deﬁned log-diﬀerent function log-diﬀX on X(Q).

## ★★何を取るか —— 定義に**値**を与える

`Found/GenEll/LogDiff.lean` で `log-diff_X(x) ≝ deg_F(δ_x)` を構成した。
★しかしそれは**式**であって**値**ではない。ここで値を与える:

> **`deg(δ_F) = log |disc F|`**、したがって **`log-diff = log|disc F| / [F:ℚ]`**

★mathlib に **`NumberField.absNorm_differentIdeal : absNorm (differentIdeal ℤ 𝒪) = |discr K|`**
があるので、あとは「算術因子の次数はイデアルノルムの対数である」を繋げばよい。

## ★機構: `deg(idealADiv I) = log(absNorm I)`

イデアルの素因子分解 `I = ∏ᶠ_v v^{count_v}`(`finprod_heightOneSpectrum_factorization`)に
`absNorm` を掛けると `absNorm I = ∏ᶠ_v (absNorm v)^{count_v}`。
対数を取れば `Σᶠ_v count_v · log(absNorm v)` で、これが `deg(idealADiv I)` の定義そのものである。

★**積公式(`ProductFormula.lean`)と同じ形**——`finprod` を `Real.log_finprod` で `finsum` にし、
`Finsupp.sum` へ落とす。道具は使い回せた。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain Ideal

variable (F : Type*) [Field F] [NumberField F]

/-! ## ★イデアルの次数はノルムの対数 -/

/-- ★**`deg(idealADiv I) = log(absNorm I)`**。 -/
theorem deg_idealADiv (I : Ideal (𝓞 F)) (hI : I ≠ 0) :
    deg (idealADiv F I) = Real.log (Ideal.absNorm I) := by
  classical
  -- 分解に `absNorm` を掛ける
  have hfact : ∏ᶠ v : FinitePlace F, v.maxPowDividing I = I :=
    finprod_heightOneSpectrum_factorization hI
  have hpos : ∀ v : FinitePlace F, (0 : ℝ) < (Ideal.absNorm (v.maxPowDividing I) : ℝ) := by
    intro v
    have : v.maxPowDividing I ≠ 0 := by
      simp only [IsDedekindDomain.HeightOneSpectrum.maxPowDividing]
      exact pow_ne_zero _ v.ne_bot
    have := (Ideal.absNorm_eq_zero_iff (I := v.maxPowDividing I)).not.2 this
    positivity
  have hprod : (Ideal.absNorm I : ℝ)
      = ∏ᶠ v : FinitePlace F, (Ideal.absNorm (v.maxPowDividing I) : ℝ) := by
    conv_lhs => rw [← hfact]
    exact ((Nat.castRingHom ℝ).toMonoidHom.comp Ideal.absNorm.toMonoidHom).map_finprod_of_preimage_one
      (by simp) _
  rw [hprod, Real.log_finprod hpos]
  -- 各項を `count_v · log q_v` に直す
  have hterm : ∀ v : FinitePlace F,
      Real.log (Ideal.absNorm (v.maxPowDividing I) : ℝ)
        = (((Associates.mk v.asIdeal).count (Associates.mk I).factors : ℤ) : ℝ)
            * Real.log (residueCard v) := by
    intro v
    simp only [IsDedekindDomain.HeightOneSpectrum.maxPowDividing, map_pow, residueCard]
    push_cast
    rw [Real.log_pow]
  rw [finsum_congr hterm]
  -- `finsum` を `Finsupp.sum` に落とす
  simp only [idealADiv, dif_neg hI, deg, ADiv.fin, ADiv.arc]
  rw [Finsupp.sum_zero_index, add_zero, Finsupp.sum]
  symm
  refine finsum_eq_finsetSum_of_support_subset _ ?_
  intro v hv
  rw [Function.mem_support] at hv
  have hc : (((Associates.mk v.asIdeal).count (Associates.mk I).factors : ℕ) : ℤ) ≠ 0 := by
    intro hcon
    apply hv
    rw [hcon]
    push_cast
    ring
  simp only [Finset.mem_coe, Finsupp.mem_support_iff]
  exact hc

/-! ## ★★`log-diff` の値 —— 判別式である -/

/-- ★**`deg(δ_F) = log |disc F|`**。

原文 (GenEll p.8):
> (iii) Let x ∈ X(F ) ⊆ X(Q), where F is a minimal field of definition of x. Then the diﬀerent ideal of F determines an eﬀective arithmetic divisor

★mathlib の `NumberField.absNorm_differentIdeal`(差積イデアルのノルムは絶対判別式)を、
上の `deg_idealADiv` に繋ぐだけである。 -/
theorem deg_differentADiv :
    deg (differentADiv F) = Real.log ((NumberField.discr F).natAbs) := by
  have hne : differentIdeal ℤ (𝓞 F) ≠ 0 := differentIdeal_ne_bot
  rw [differentADiv, deg_idealADiv F _ hne, NumberField.absNorm_differentIdeal F (𝓞 F)]

/-- ★★**[GenEll] Definition 1.5, (iii) の値**。

> **`log-diff = log|disc F| / [F:ℚ]`**

★これで `log-diff` は「型が付くだけの式」ではなく**計算できる量**になった。
`ArithDiv.lean` の下線つき `deg`(正規化次数)で割るところまで含めて、
原文の定義がそのまま古典的な公式に落ちる。 -/
theorem logDiffOfField_eq :
    logDiffOfField F
      = Real.log ((NumberField.discr F).natAbs) / (Module.finrank ℚ F : ℝ) := by
  rw [logDiffOfField, degNormalized, deg_differentADiv]

/-! ## ★出典の紐付け(`.src`) -/

def deg_idealADiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

def logDiffOfField_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5, (iii)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
