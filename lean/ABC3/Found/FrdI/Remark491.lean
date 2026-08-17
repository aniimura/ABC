/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Example43
import ABC3.Found.FrdI.Prop44Gp
import ABC3.Found.FrdI.Def45

/-!
# [FrdI] Remark 4.9.1 —— `Example 4.3` の Frobenioid は rational 型でない

原文 (FrdI p.90):
> One verifies immediately that the Frobenioid of Example 4.3 is
> not of rational type.

★★**「immediate」の中身を測る。** 要点は 2 つである。

1. **`Prime(ℚ≥0)` は空でない** —— `1` が primary(`ℚ≥0` は Archimedes 的なので
   0 でない元はすべて互いに `MPrec` 同値、したがって primary)。
2. ★★**`Φ^birat(A) = 0`** —— `𝒪^×(A^birat)` の元は Frobenius 次数 1 であり、
   `Example 4.3` では**次数が決まれば零因子も決まる**(`Div φ = b − deg·a`)ので、
   代表元 `(α, φ)` の零因子が一致し、差が 0 になる。

したがって `toGp a − toGp b ∈ Φ^birat(A)` は `a = b` を強い、
`p ∈ Supp(a)` かつ `p ∉ Supp(b)` と両立しない。★これは **`ι` の取り方によらない**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

/-! ## ★1. `Prime(ℚ≥0)` は空でない -/

/-- ★`ℚ≥0` では 0 でない元はすべて `1` を `MPrec` で越える(Archimedes 性)。 -/
theorem mprec_one_nnrat {b : ℚ≥0} (hb : b ≠ 0) : MPrec (1 : ℚ≥0) b := by
  obtain ⟨n, hn⟩ := exists_nat_ge (1 / b : ℚ≥0)
  refine ⟨n + 1, Nat.succ_pos n, mle_nnrat_iff.mpr ?_⟩
  have hb' : (0 : ℚ≥0) < b := pos_iff_ne_zero.mpr hb
  have h1 : (1 : ℚ≥0) = (1 / b) * b := by
    rw [div_mul_cancel₀ _ hb]
  have h2 : (1 / b : ℚ≥0) * b ≤ ((n : ℚ≥0) + 1) * b :=
    mul_le_mul_of_nonneg_right (le_trans hn (by simp)) zero_le
  calc (1 : ℚ≥0) = (1 / b) * b := h1
    _ ≤ ((n : ℚ≥0) + 1) * b := h2
    _ = ((n + 1 : ℕ) • b) := by
        rw [nsmul_eq_mul]
        push_cast
        ring

/-- ★`1 ∈ ℚ≥0` は primary。 -/
theorem isPrimaryElt_one_nnrat : IsPrimaryElt (1 : ℚ≥0) :=
  ⟨one_ne_zero, fun _ hb _ => mprec_one_nnrat hb⟩

/-- ★`Prime(ℚ≥0)` の元。 -/
noncomputable def primeNNRat : Prime ℚ≥0 := toPrime ℚ≥0 1 isPrimaryElt_one_nnrat

/-! ## ★2. `Example 4.3` では `Φ^birat(A) = 0` -/

/-- ★次数 1 の射の零因子は始域と終域だけで決まる。 -/
theorem ex43Div_eq_of_deg_one {a b : Ex43} (f g : a ⟶ b) (hf : f.deg = 1) (hg : g.deg = 1) :
    ex43Div f = ex43Div g :=
  NNRat.coe_injective (by rw [ex43Div_deg_one f hf, ex43Div_deg_one g hg])

/-- ★★`𝒞^birat` の**次数 1** の自己射の `Φ^gp` 零因子は 0。 -/
theorem ex43_biratDivGp_eq_zero (A : BiratCat ex43P ex43_frobenioid) (δ : A ⟶ A)
    (hδ : (biratPre ex43P ex43_frobenioid).degFr δ = 1) :
    biratDivGp δ = 0 := by
  obtain ⟨Z, φ, hrep⟩ := HomBirat.exists_rep δ
  have hdeg : ex43P.degFr φ = 1 := by
    rw [← biratDeg_mk Z φ, hrep]
    exact hδ
  have hza : Z.unop.hom.hom.deg = 1 := Z.unop.hom.property.2.1
  have hdiv : ex43P.Div φ = ex43P.Div Z.unop.hom.hom :=
    ex43Div_eq_of_deg_one φ Z.unop.hom.hom hdeg hza
  rw [← hrep, biratDivGp_mk, sliceDivGpOf_eq]
  rw [hdiv, hdeg]
  rw [PNat.one_coe, one_smul, sub_self, map_zero]
  rfl

/-- ★★★`Example 4.3` では `Φ^birat(A)` は自明。 -/
theorem ex43_phiBiratAt_eq_bot (A : BiratCat ex43P ex43_frobenioid) :
    phiBiratAt ex43P ex43_frobenioid A = ⊥ := by
  refine le_antisymm ?_ bot_le
  rintro x ⟨δ, hδ, rfl⟩
  exact ex43_biratDivGp_eq_zero A (δ : A ⟶ A) hδ.1.2

/-! ## ★3. rational 型でない

原文 (FrdI p.90):
> One verifies immediately that the Frobenioid of Example 4.3 is
-/

/-- ★★どの対象も **strictly rational でない**(`ι` の取り方によらない)。 -/
theorem ex43_not_strictlyRational
    (ι : ∀ Y : Discrete PUnit, Prime (Phi43.val Y) → Pf (Phi43.val Y) → ℝ≥0) (A : Ex43) :
    ¬ IsStrictlyRational ex43P ex43_frobenioid ι A := by
  intro hsr
  obtain ⟨a0, b0, hmem, hpa, hpb⟩ := hsr primeNNRat
  obtain ⟨a, rfl⟩ : ∃ a : ℚ≥0, a = a0 := ⟨a0, rfl⟩
  obtain ⟨b, rfl⟩ : ∃ b : ℚ≥0, b = b0 := ⟨b0, rfl⟩
  rw [ex43_phiBiratAt_eq_bot] at hmem
  have hgp : toGp _ a = toGp _ b := sub_eq_zero.mp hmem
  obtain ⟨c, hc⟩ := toGp_eq_iff.mp hgp
  have hab : a = b := add_right_cancel hc
  exact hpb (hab ▸ hpa)

/-- ★★★★**[FrdI] Remark 4.9.1** —— `Example 4.3` の Frobenioid は **rational 型でない**。

原文 (FrdI p.90):
> not of rational type.
-/
theorem ex43_not_rationalType
    (ι : ∀ Y : Discrete PUnit, Prime (Phi43.val Y) → Pf (Phi43.val Y) → ℝ≥0) :
    ¬ IsOfRationalType Ex43 ex43P ex43_frobenioid ι := by
  intro hr
  obtain ⟨B, _, _, hsr⟩ := hr (Ex43.mk 0)
  exact ex43_not_strictlyRational ι B hsr

def ex43_not_rationalType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 90, item := "Remark 4.9.1",
    sectionId := "frdi-remark-4-9-1" }

end ABC3.Found.FrdI
