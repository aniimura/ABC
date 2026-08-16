import ABC3.Found.FrdI.Def31

/-!
# [FrdI] Proposition 3.3 —— Base-identity Pre-steps and Units

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.59–p.61。

原文 (FrdI p.59):
> (Base-identity Pre-steps and Units)

## ★この命題の規模(測定)

**5 条、主張は 8**:

| 条 | # | 内容 | 状態 |
|---|---|---|---|
| (i) | 1 | `𝒟` が Frobenius-slim なら `𝔽 → End(𝒞^pl-bk_A → 𝒞)^bs-iso` の `1` の像は base-identity pre-step | 未 |
| (i) | 2 | 逆に Frobenius-normalized ＋ `A` が Frobenius-trivial なら、どの base-identity pre-step もそう現れる | 未 |
| (ii) | 3 | `α₁ ≈ α₂`(unit-equivalent)⟺ `degFr`・`Div`・`Base` が一致 —— **必要性** | ★**ここで実装** |
| (ii) | 4 | 同 —— **十分性** | 未 |
| (iii) | 5 | `𝒞^istr → 𝒞^un-tr` は full かつ本質的全射、圏同値 ⟺ `𝒞^istr` が unit-trivial 型 | 未 |
| (iv) | 6 | `𝒞^un-tr → 𝔽_Φ` は忠実かつ本質的全射で、`𝒞^un-tr` に Frobenioid の構造を与える | 未 |
| (iv) | 7 | `𝒞^un-tr` の射の類型は `𝒞^istr` の射の類型から来る | 未 |
| (v) | 8 | `𝒞 → 𝔽_Φ` が圏同値 ⟺ `𝒞` が Aut-ample・unit-trivial・base-trivial 型 | 未 |

★**十分性の証明は原文どおり 4 手**(`Definition 1.3, (iv), (a)` の 3 分解、
`Definition 1.3, (ii)` の Frobenius 型射の本質的一意性、
`Definition 1.3, (i), (c)` の pull-back の圏同値、`Definition 1.3, (vi)` の
faithfulness up to units)であり、`Remark 2.7.2` の分解と同じ形をしている。
★**別途取る。**
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(ii) の必要性

原文 (FrdI p.59):
> three conditions are satisfied: (a) degFr(

★**理由は「`𝒪^×` の元は LB-invertible な base-identity linear 自己射」**
(原文どおり)——`Div = 0`、`degFr = 1`、`Base = 𝟙` なので、
3 つの量のどれにも寄与しない。
-/

/-- ★★★**[FrdI] Proposition 3.3, (ii) の必要性** —— unit-equivalent な 2 射は
`Frobenius 次数`・`Div`・`Base` がすべて一致する。

★すなわち **`𝔽_Φ` へ同じ射に写る**。 -/
theorem prop_3_3_ii_necessity {A B : C} {α₁ α₂ : A ⟶ B}
    (h : IsUnitEquivalent P α₁ α₂) :
    P.degFr α₁ = P.degFr α₂ ∧ P.Div α₁ = P.Div α₂ ∧ P.Base α₁ = P.Base α₂ := by
  obtain ⟨Cc, γ, β, δ, hδ, h₁, h₂⟩ := h
  haveI : IsIso ((δ : Cc ⟶ Cc)) := isIso_of_mem_otimes P hδ
  -- `δ` の 3 つの量
  have hdb : P.Base ((δ : Cc ⟶ Cc)) = 𝟙 _ := by
    have h : P.Base ((δ : Cc ⟶ Cc)) = P.Base (𝟙 Cc) := hδ.1.1
    rwa [P.Base_id] at h
  have hdd : P.degFr ((δ : Cc ⟶ Cc)) = 1 := hδ.1.2
  have hdv : P.Div ((δ : Cc ⟶ Cc)) = 0 := isIsometric_of_isIso P _
  -- `δ ≫ β` は `β` と 3 つの量が一致する
  have hdegβ : P.degFr (((δ : Cc ⟶ Cc)) ≫ β) = P.degFr β := by
    rw [P.degFr_comp, hdd, mul_one]
  have hbaseβ : P.Base (((δ : Cc ⟶ Cc)) ≫ β) = P.Base β := by
    rw [P.Base_comp, hdb, Category.id_comp]
  have hdivβ : P.Div (((δ : Cc ⟶ Cc)) ≫ β) = P.Div β := by
    rw [P.Div_comp, hdb, MonoidOn.map_id, hdv]
    simp
  refine ⟨?_, ?_, ?_⟩
  · rw [h₁, h₂, P.degFr_comp, P.degFr_comp, hdegβ]
  · rw [h₁, h₂, P.Div_comp, P.Div_comp, hdegβ, hdivβ]
  · rw [h₁, h₂, P.Base_comp, P.Base_comp, hbaseβ]

/-- ★**`𝔽_Φ` へ同じ射に写る**という原文の言い方そのもの。 -/
theorem prop_3_3_ii_toElem {A B : C} {α₁ α₂ : A ⟶ B} (h : IsUnitEquivalent P α₁ α₂) :
    P.toElem.map α₁ = P.toElem.map α₂ := by
  obtain ⟨hd, hv, hb⟩ := prop_3_3_ii_necessity P h
  exact ElemFrobCat.Hom.ext hb hv hd

end ABC3.Found.FrdI
