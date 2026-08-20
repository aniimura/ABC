/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop56Cos
import ABC3.Found.FrdI.Thm51Cls

/-!
# [FrdI] Proposition 5.6 —— `E` の分解と `M_p` の `shift` 条件

原文 (FrdI p.106):
> that any f ∈E admits a factorization f = f0 · f1, where f0 is an automorphism,

★`Prop56Cos.lean` は「`IsShiftStable E M` から商単系 `E_M` が出る」ところまでを
純粋な単系の代数として実装した。★**本ファイルはその仮定を Frobenioid の側で外す**。

## ★原文の 2 段

1. **分解**: `Theorem 5.1, (iii)`(isotropic 型では Frobenius-trivial 対象は Aut-ample)
   により、base-isomorphism な自己射 `f` は `f = f₀ · f₁`
   (`f₀` は自己同型、`f₁` は base-identity)に分解する。
2. **入れ替え**: `𝒞` が model 型ゆえ **Frobenius-normalized** なので
   `f₁ · m = m^d · f₁`(`d = deg_Fr(f₁)`)。したがって
   `f · m = (f₀ · m^d · f₀⁻¹) · f` となり、`m' := f₀ · m^d · f₀⁻¹`。

★★**したがって `shift` 条件は「`M` が自己同型による共役で閉じている」ことに帰着する**。
原文の `M_p`(`l ≠ p` の pro-`l` 部分が位相的に生成する閉部分群)は
**標準的に定まる**ので共役で閉じている —— その構成自体は別の葉
(鎖 `prop56` の `p56-Mp` / `p56-limit`)である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★base-isomorphism な自己射のなす部分単系(原文 p.106 の `E`)。 -/
def baseIsoEnd (P : PreFrobenioid C Φ) (A : C) : Submonoid (End A) where
  carrier := {f | IsBaseIsomorphism P ((f : A ⟶ A))}
  one_mem' := by
    show IsIso (P.Base (𝟙 A))
    rw [P.Base_id]; infer_instance
  mul_mem' {x y} hx hy := by
    show IsIso (P.Base ((y : A ⟶ A) ≫ (x : A ⟶ A)))
    haveI : IsIso (P.Base ((x : A ⟶ A))) := hx
    haveI : IsIso (P.Base ((y : A ⟶ A))) := hy
    rw [P.Base_comp]; infer_instance

/-- ★`𝒪^▷(A)` の元は base-isomorphism(base-identity だから)。 -/
theorem mem_baseIsoEnd_of_otri {A : C} {f : End A} (hf : f ∈ OTri P A) :
    f ∈ baseIsoEnd P A := by
  show IsIso (P.Base ((f : A ⟶ A)))
  have : P.Base ((f : A ⟶ A)) = P.Base (𝟙 A) := hf.1
  rw [this, P.Base_id]
  infer_instance

/-- ★★★**原文 p.106 の「`f = f₀ · f₁`」** —— base-isomorphism な自己射は
「自己同型 ＋ base-identity 自己射」に分解する。

★★`Theorem 5.1, (iii)`(isotropic 型では Frobenius-trivial 対象は Aut-ample)から。 -/
theorem baseIsoEnd_factor (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A : C} (hA : IsFrobeniusTrivial P A) {f : End A} (hf : f ∈ baseIsoEnd P A) :
    ∃ f₀ : End A, ∃ _ : IsIso ((f₀ : A ⟶ A)), ∃ f₁ : End A,
      IsBaseIdentity P f₁ ∧ f = f₀ * f₁ := by
  haveI hfi : IsIso (P.Base ((f : A ⟶ A))) := hf
  obtain ⟨f₀, hf₀iso, hf₀b⟩ := isAutAmple_of_frobTrivial G hiso hA (P.Base ((f : A ⟶ A))) hfi
  haveI := hf₀iso
  refine ⟨f₀, hf₀iso, (f : A ⟶ A) ≫ inv ((f₀ : A ⟶ A)), ?_, ?_⟩
  · show P.Base ((f : A ⟶ A) ≫ inv ((f₀ : A ⟶ A))) = P.Base (𝟙 A)
    rw [P.Base_comp, ← hf₀b, ← P.Base_comp, IsIso.hom_inv_id]
  · show (f : A ⟶ A) = ((f : A ⟶ A) ≫ inv ((f₀ : A ⟶ A))) ≫ (f₀ : A ⟶ A)
    rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id]

/-- ★★★★**原文 p.106 の「`f · m = m' · f`」の条は自動で成り立つ** ——
`M ⊆ 𝒪^×(A)` が**自己同型による共役で閉じている**なら。

★分解 `f = f₀ · f₁`(`baseIsoEnd_factor`)と Frobenius-normalized
(`f₁ · m = m^d · f₁`)から `m' = f₀ · m^d · f₀⁻¹`。 -/
theorem baseIsoEnd_shift (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A : C} (hA : IsFrobeniusTrivial P A) (hfn : IsFrobeniusNormalized P A)
    (M : Submonoid (End A))
    (hMle : ∀ m ∈ M, m ∈ OTri P A)
    (hAut : ∀ θ θi : End A, θ * θi = 1 → θi * θ = 1 → ∀ m ∈ M, θ * m * θi ∈ M) :
    ∀ m ∈ M, ∀ f ∈ baseIsoEnd P A, ∃ m' ∈ M, f * m = m' * f := by
  intro m hm f hf
  obtain ⟨f₀, hf₀iso, f₁, hf₁b, rfl⟩ := baseIsoEnd_factor G hiso hA hf
  haveI := hf₀iso
  set f₀i : End A := inv ((f₀ : A ⟶ A)) with hf₀i
  have h1 : f₀ * f₀i = 1 := by
    show (inv ((f₀ : A ⟶ A))) ≫ (f₀ : A ⟶ A) = 𝟙 A
    exact IsIso.inv_hom_id _
  have h2 : f₀i * f₀ = 1 := by
    show (f₀ : A ⟶ A) ≫ (inv ((f₀ : A ⟶ A))) = 𝟙 A
    exact IsIso.hom_inv_id _
  refine ⟨f₀ * m ^ ((P.degFr ((f₁ : A ⟶ A)) : ℕ+) : ℕ) * f₀i,
    hAut f₀ f₀i h1 h2 _ (M.pow_mem hm _), ?_⟩
  calc (f₀ * f₁) * m = f₀ * (f₁ * m) := by rw [mul_assoc]
    _ = f₀ * (m ^ ((P.degFr ((f₁ : A ⟶ A)) : ℕ+) : ℕ) * f₁) :=
        congrArg (fun t => f₀ * t) (hfn f₁ hf₁b m (hMle m hm)).symm
    _ = (f₀ * m ^ ((P.degFr ((f₁ : A ⟶ A)) : ℕ+) : ℕ) * f₀i) * (f₀ * f₁) := by
        rw [mul_assoc, mul_assoc, ← mul_assoc f₀i f₀ f₁, h2, one_mul, ← mul_assoc]

/-- ★★★★★**`IsShiftStable` が Frobenioid の側で成り立つ** ——
したがって商単系 `E_M`(`Prop56Cos.lean`)が使える。

★仮定は 2 つだけ: `M` が群であること、`M` が**自己同型による共役で閉じている**こと。
★★原文の `M_p` は標準的に定まる(`l ≠ p` の pro-`l` 部分が生成する閉部分群)ので
どちらも満たすはずである —— その構成自体は鎖 `prop56` の `p56-Mp` に残っている。 -/
theorem isShiftStable_baseIsoEnd (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A : C} (hA : IsFrobeniusTrivial P A) (hfn : IsFrobeniusNormalized P A)
    (M : Submonoid (End A))
    (hMle : ∀ m ∈ M, m ∈ OTri P A)
    (hMinv : ∀ x ∈ M, ∃ y ∈ M, x * y = 1 ∧ y * x = 1)
    (hAut : ∀ θ θi : End A, θ * θi = 1 → θi * θ = 1 → ∀ m ∈ M, θ * m * θi ∈ M) :
    IsShiftStable (baseIsoEnd P A) M where
  le := fun _ hm => mem_baseIsoEnd_of_otri (hMle _ hm)
  inv := hMinv
  shift := baseIsoEnd_shift G hiso hA hfn M hMle hAut

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.6` の証明中の 2 段(分解と入れ替え)。
★**条つき**: `M_p` の構成(pro-`l` 分解を使う)は未実装。 -/
def isShiftStable_baseIsoEnd.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 106,
    item := "Proposition 5.6 — E の分解 f = f₀·f₁ と M の shift 条件",
    sectionId := "frdi-prop-5-6" }

end ABC3.Found.FrdI
