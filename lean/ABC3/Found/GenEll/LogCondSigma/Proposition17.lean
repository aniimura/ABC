/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SigmaBound
import ABC3.Found.GenEll.CartierPullback
import ABC3.Found.GenEll.LogCondSigma.Definition15

/-!
# LogCondSigma —— `[GenEll] Proposition 1.7` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

/-- ★★★★★**導手の係数は `≤ 1`** —— これが `SigmaBound` へ渡す形。 -/
theorem idealADiv_radical_coeff_le_one {F : Type*} [Field F] [NumberField F]
    (I : Ideal (𝓞 F)) (hI : I ≠ 0) (v : FinitePlace F) :
    (idealADiv F I.radical).fin v ≤ 1 := by
  classical
  have hrad0 : I.radical ≠ 0 := radical_ne_zero hI
  have hval : (idealADiv F I.radical).fin v
      = ((Associates.mk v.asIdeal).count (Associates.mk I.radical).factors : ℤ) := by
    unfold idealADiv
    rw [dif_neg hrad0]
    rfl
  rw [hval]
  exact_mod_cast count_radical_le_one I hI v

/-! ## ★★★★★★★★★`log-cond` の `Σ`-限界 -/

variable {X : AlgebraicGeometry.Scheme.{0}}

/-- ★★★★★★★★★**`Σ` の上に台を持つ `log-cond` は `Σ_{q∈Σ} log q` 以下**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★★**右辺は `Σ` だけで決まる定数**であり、点 `x` にも定義体 `F` にも依らない。
★★★これが `Remark 1.5.1`(BD-class が `(X_ℚ, D_ℚ)` だけに依る)と
`Proposition 1.7`(`Σ` 上の寄与は `≈ 0`)が共通して使う形である。 -/
theorem logCond_le_sum_log (F : Type) [Field F] [NumberField F]
    (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X)
    (hI : pullbackIdeal F D xF ≠ 0)
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ (conductorADiv F D xF).fin.support, ch v ∈ Sig)
    (hover : ∀ v ∈ (conductorADiv F D xF).fin.support,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    logCond F D xF ≤ ∑ q ∈ Sig, Real.log q := by
  rw [logCond]
  refine degNormalized_le_sum_log (conductorADiv F D xF) ?_ ?_ Sig hprime ch hmem hover
  · rw [conductorADiv]
    exact idealADiv_arc F _
  · intro v
    rw [conductorADiv]
    exact idealADiv_radical_coeff_le_one _ hI v

def logCond_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i)(Σ 上の log-cond の寄与は Σ_{q∈Σ} log q で抑えられる)",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
