/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop41
import ABC3.Found.FrdI.Prop110
import ABC3.Found.FrdI.Prop32Equiv
import ABC3.Found.FrdI.Def24SuppElt

/-!
# [FrdI] `Proposition 4.1` の**圏層** —— `Definition 1.3, (iii)(d)` による翻訳

原文 (FrdI p.76):
> Proof. First, we consider assertion (i). By applying the second equivalence of

原文 (FrdI p.76):
> categories of Definition 1.3, (iii), (d), to the various pre-steps over A, it follows

## ★★原文の証明は 5 条とも同じ形をしている

**`Definition 1.3, (iii), (d)` の後置の圏同値 `(𝒞^coa-pre)_A ≃ Order(Φ(A))^opp` を
「`A` の上の pre-step」に当てて、圏の条件を単系の条件に書き換える。**
そのうえで単系の側を解く。

★`Prop41.lean` が単系の側((i) の分)を持つ。★**本ファイルが圏の側**である。

## ★★★翻訳の 3 本の柱

| 柱 | 宣言 | 内容 |
|---|---|---|
| 値 | `preStepVal` | `A` の上の pre-step `φ` に `x_φ := (Base φ)⁻¹^*(Div φ) ∈ Φ(A)` を割り当てる |
| 分解 → 和 | `preStepVal_comp` | `φ = ψ ≫ χ` なら `x_φ = x_χ + (Base φ)⁻¹^*(Div ψ)` |
| 和 → 分解 | `exists_factor_of_mle'` | `a ≼ x_φ` なら `φ = ψ ≫ χ`・`x_χ = a` と分解できる(本質的全射 ＋ 充満) |
| 通す | `exists_factor_through` | `x_φ ≼ x_ζ` なら `ζ = ζ' ≫ φ`(充満性のみ) |
| 図式 | `preStepVal_of_frob_square` | `β' ≫ ζ = φ_A ≫ α_n` なら `x_ζ = n · x_{φ_A}` |

★★**「step ⟺ 値が 0 でない」**(`isStep_iff_preStepVal_ne_zero`)が isotropic 型で効く
—— 等長な pre-step が同型になるからである。

## ★(i) の形

原文 (FrdI p.75):
> (i) φ is primary if and only if, for every factorization φ = φA ◦φB, where

図式
```
B --φ_B--> B' --φ_A--> A
           |β'         |α_n
           v           v
           B'' --ζ--> A          (ζ = ζ' ≫ φ, ζ' : B'' ⟶ B は pre-step)
```
は「`x_φ ⪯ x_{φ_A}`」そのものである。★したがって (i) は
`isPrimaryElt_iff_of_perfect`(`Prop41.lean` の単系層)と繋ぐだけで閉じる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe v u w v2 u2 w'

/-! ## ★0. 単系の側の小補題 -/

section Monoid

variable {M : Type w'} {N : Type w'} [AddCommMonoid M] [AddCommMonoid N]

/-- ★`primary` は加法的全単射で移る。 -/
theorem isPrimaryElt_of_bijective (f : M →+ N) (hf : Function.Bijective f) {a : M}
    (h : IsPrimaryElt a) : IsPrimaryElt (f a) := by
  refine ⟨fun h0 => h.1 (hf.1 (h0.trans (map_zero f).symm)), fun b hb hprec => ?_⟩
  obtain ⟨b', rfl⟩ := hf.2 b
  have hb' : b' ≠ 0 := fun h0 => hb (by rw [h0, map_zero])
  obtain ⟨n, hn, c, hc⟩ := hprec
  obtain ⟨c', rfl⟩ := hf.2 c
  have hprec' : MPrec b' a := by
    refine ⟨n, hn, c', hf.1 ?_⟩
    rw [map_add, map_nsmul]
    exact hc
  obtain ⟨m, hm, d, hd⟩ := h.2 b' hb' hprec'
  exact ⟨m, hm, f d, by rw [← map_nsmul, ← map_add, hd]⟩

end Monoid

section Cat

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-! ## ★1. 「後置の値」と合成則 -/

/-- ★★**「後置の値」の合成則** —— `x_{ψ≫χ} = x_χ + (Base(ψ≫χ))⁻¹^*(Div ψ)`。

★★これが `Proposition 4.1` の圏層の**核**である ——
`Definition 1.3, (iii)(d)` の後置の圏同値で「`A` の上の pre-step」を `Φ(A)` に写すと、
**分解 `φ = ψ ≫ χ` が和の分解 `x_φ = x_χ + y` に化ける**。 -/
theorem overVal_comp (Q : PreFrobenioid C Φ) {A B B' : C} (ψ : B ⟶ B') (χ : B' ⟶ A)
    (hψ : IsIso (Q.Base ψ)) (hχ : IsIso (Q.Base χ)) (hc : IsIso (Q.Base (ψ ≫ χ)))
    (hlin : Q.degFr χ = 1) :
    Φ.map (@inv _ _ _ _ (Q.Base (ψ ≫ χ)) hc) (Q.Div (ψ ≫ χ))
      = Φ.map (@inv _ _ _ _ (Q.Base χ) hχ) (Q.Div χ)
        + Φ.map (@inv _ _ _ _ (Q.Base (ψ ≫ χ)) hc) (Q.Div ψ) := by
  have hinv : @inv _ _ _ _ (Q.Base (ψ ≫ χ)) hc
      = @inv _ _ _ _ (Q.Base χ) hχ ≫ @inv _ _ _ _ (Q.Base ψ) hψ := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [Q.Base_comp, Category.assoc, ← Category.assoc (Q.Base χ), IsIso.hom_inv_id,
      Category.id_comp]
    exact @IsIso.hom_inv_id _ _ _ _ (Q.Base ψ) hψ
  have hdiv : Q.Div (ψ ≫ χ) = Φ.map (Q.Base ψ) (Q.Div χ) + Q.Div ψ := by
    rw [Q.Div_comp, hlin]
    simp
  rw [hdiv, hinv, Φ.map_comp, map_add, MonoidOn.map_map_inv Φ (Q.Base ψ) hψ,
    map_add, ← Φ.map_comp]

variable {P : PreFrobenioid C Φ}

/-- ★★**pre-step の「後置の値」** —— `Definition 1.3, (iii)(d)` の後置の圏同値が
`A` の上の co-angular pre-step に割り当てる `Φ(A)` の元。 -/
noncomputable def preStepVal (P : PreFrobenioid C Φ) {A B : C} (φ : B ⟶ A)
    (h : IsPreStep P φ) : Φ.val (P.toElem.obj A).base :=
  Φ.map (@inv _ _ _ _ (P.Base φ) h.2) (P.Div φ)

@[simp] theorem preStepVal_eq_zero_iff {A B : C} (φ : B ⟶ A) (h : IsPreStep P φ) :
    preStepVal P φ h = 0 ↔ P.Div φ = 0 := by
  haveI := h.2
  constructor
  · intro h0
    exact Φ.map_injective (@inv _ _ _ _ (P.Base φ) h.2) (h0.trans (map_zero _).symm)
  · intro h0
    show Φ.map _ (P.Div φ) = 0
    rw [h0, map_zero]

/-- ★★**isotropic 型では「step ⟺ 後置の値が 0 でない」**。 -/
theorem isStep_iff_preStepVal_ne_zero (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (φ : B ⟶ A) (h : IsPreStep P φ) :
    IsStep P φ ↔ preStepVal P φ h ≠ 0 := by
  rw [Ne, preStepVal_eq_zero_iff]
  constructor
  · rintro ⟨-, hni⟩ h0
    exact hni (hiso B A φ h0 h)
  · intro hd
    exact ⟨h, fun _ => hd (isIsometric_of_isIso P φ)⟩

/-- ★**「後置の値」の合成則**(`preStepVal` 版)。 -/
theorem preStepVal_comp {A B B' : C} (ψ : B ⟶ B') (χ : B' ⟶ A)
    (hψ : IsPreStep P ψ) (hχ : IsPreStep P χ) (hc : IsPreStep P (ψ ≫ χ)) :
    preStepVal P (ψ ≫ χ) hc
      = preStepVal P χ hχ + Φ.map (@inv _ _ _ _ (P.Base (ψ ≫ χ)) hc.2) (P.Div ψ) :=
  overVal_comp P ψ χ hψ.2 hχ.2 hc.2 hχ.1

/-- ★**`primary` は「後置の値」と零因子で同じこと**(底に沿った全単射で移る)。 -/
theorem isPrimaryElt_preStepVal_iff {A B : C} (φ : B ⟶ A) (h : IsPreStep P φ) :
    IsPrimaryElt (preStepVal P φ h) ↔ IsPrimaryElt (P.Div φ) := by
  haveI := h.2
  have hbij : Function.Bijective (Φ.map (@inv _ _ _ _ (P.Base φ) h.2)) :=
    Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base φ) h.2).symm
  have hbij' : Function.Bijective (Φ.map (P.Base φ)) :=
    Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base φ) h.2)
  constructor
  · intro hp
    have h2 := isPrimaryElt_of_bijective (Φ.map (P.Base φ)) hbij' hp
    rwa [show Φ.map (P.Base φ) (preStepVal P φ h) = P.Div φ from
      MonoidOn.map_inv_map Φ (P.Base φ) h.2 (P.Div φ)] at h2
  · intro hp
    exact isPrimaryElt_of_bijective (Φ.map (@inv _ _ _ _ (P.Base φ) h.2)) hbij hp

/-! ## ★2. (iii)(d) の後置の圏同値を使う 2 本 -/

/-- ★★★**充満性** —— `x_φ ≼ x_ζ` なら `ζ` は `φ` を**通って**分解する。 -/
theorem exists_factor_through (G : Frobenioid P) {A B Z : C} (φ : B ⟶ A) (ζ : Z ⟶ A)
    (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    (hcζ : IsCoAngular P ζ) (hsζ : IsPreStep P ζ)
    (hle : MLe (preStepVal P φ hsφ) (preStepVal P ζ hsζ)) :
    ∃ ζ' : Z ⟶ B, IsCoAngular P ζ' ∧ IsPreStep P ζ' ∧ ζ = ζ' ≫ φ := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  let Zφ : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨φ, ⟨hcφ, hsφ⟩⟩)
  let Zζ : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨Z⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨ζ, ⟨hcζ, hsζ⟩⟩)
  obtain ⟨f, -⟩ := (coaPreOverFunctor P A).map_surjective
    (show (coaPreOverFunctor P A).obj Zζ ⟶ (coaPreOverFunctor P A).obj Zφ from
      (homOfLE (show (toOrderCat (preStepVal P φ hsφ)
          : OrderCat (Φ.val (P.toElem.obj A).base))
        ≤ toOrderCat (preStepVal P ζ hsζ) from hle)).op)
  exact ⟨f.left.hom, f.left.property.1, f.left.property.2,
    (congrArg InducedWideCategory.Hom.hom (Over.w f)).symm⟩

/-- ★★★**本質的全射 ＋ 充満** —— `a ≼ x_φ` なら `φ = ψ ≫ χ` と分解でき `x_χ = a`。 -/
theorem exists_factor_of_mle' (G : Frobenioid P) {A B : C} (φ : B ⟶ A)
    (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    {a : Φ.val (P.toElem.obj A).base} (hle : MLe a (preStepVal P φ hsφ)) :
    ∃ (B' : C) (ψ : B ⟶ B') (χ : B' ⟶ A) (_ : IsPreStep P ψ) (hχ : IsPreStep P χ),
      IsCoAngular P ψ ∧ IsCoAngular P χ ∧ φ = ψ ≫ χ ∧ preStepVal P χ hχ = a := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  let Zφ : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨φ, ⟨hcφ, hsφ⟩⟩)
  let W := (coaPreOverFunctor P A).objPreimage (Opposite.op (toOrderCat a))
  have hiso : (coaPreOverFunctor P A).obj W ≅ Opposite.op (toOrderCat a) :=
    (coaPreOverFunctor P A).objObjPreimageIso _
  haveI hWb : IsIso (P.Base W.hom.hom) := W.hom.property.2.2
  have hWa : preStepVal P W.hom.hom W.hom.property.2 = a :=
    mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2
      (leOfHom hiso.inv.unop) (leOfHom hiso.hom.unop)
  obtain ⟨f, -⟩ := (coaPreOverFunctor P A).map_surjective
    (show (coaPreOverFunctor P A).obj Zφ ⟶ (coaPreOverFunctor P A).obj W from
      (homOfLE (show (toOrderCat a : OrderCat (Φ.val (P.toElem.obj A).base))
        ≤ toOrderCat (preStepVal P φ hsφ) from hle)).op ≫ hiso.inv)
  exact ⟨W.left.obj, f.left.hom, W.hom.hom, f.left.property.2, W.hom.property.2,
    f.left.property.1, W.hom.property.1,
    (congrArg InducedWideCategory.Hom.hom (Over.w f)).symm, hWa⟩

/-! ## ★3. Frobenius の四角形を通した値 -/

/-- ★`Div-identity` 自己射の底に沿った `Φ.map` は恒等(逆射の側も)。 -/
theorem map_inv_base_of_divIdentity {A : C} (α : A ⟶ A) (hD : IsDivIdentity P α)
    (h : IsIso (P.Base α)) (x : Φ.val (P.toElem.obj A).base) :
    Φ.map (@inv _ _ _ _ (P.Base α) h) x = x := by
  have hid : Φ.map (P.Base α) = Φ.map (P.Base (𝟙 A)) := hD
  rw [P.Base_id] at hid
  have h2 : Φ.map (@inv _ _ _ _ (P.Base α) h) (Φ.map (P.Base α) x) = x := by
    rw [← Φ.map_comp, @IsIso.inv_hom_id _ _ _ _ (P.Base α) h]
    exact Φ.map_id _ x
  rw [hid, Φ.map_id] at h2
  exact h2

/-- ★★**Frobenius の四角形を通した「後置の値」** —— `x_ζ = n · x_{φ_A}`。

★★これが `Proposition 4.1, (i)` の図式を単系の `⪯` に翻訳する 1 行である。 -/
theorem preStepVal_of_frob_square {A B' Y : C} {n : ℕ+} (φA : B' ⟶ A) (hφA : IsPreStep P φA)
    (αn : A ⟶ A) (hαF : IsFrobeniusType P αn) (hαD : IsDivIdentity P αn)
    (hαn : P.degFr αn = n)
    (β' : B' ⟶ Y) (hβF : IsFrobeniusType P β') (ζ : Y ⟶ A) (hζ : IsPreStep P ζ)
    (hsq : β' ≫ ζ = φA ≫ αn) :
    preStepVal P ζ hζ = ((n : ℕ+) : ℕ) • preStepVal P φA hφA := by
  haveI hbβ : IsIso (P.Base β') := hβF.2
  haveI hbζ : IsIso (P.Base ζ) := hζ.2
  haveI hbα : IsIso (P.Base αn) := hαF.2
  haveI hbφ : IsIso (P.Base φA) := hφA.2
  have hdivL : P.Div (β' ≫ ζ) = Φ.map (P.Base β') (P.Div ζ) := by
    rw [P.Div_comp, show P.Div β' = 0 from hβF.1.2, hζ.1]
    simp
  have hdivR : P.Div (φA ≫ αn) = ((n : ℕ+) : ℕ) • P.Div φA := by
    rw [P.Div_comp, show P.Div αn = 0 from hαF.1.2, map_zero, zero_add, hαn]
  have hd : Φ.map (P.Base β') (P.Div ζ) = ((n : ℕ+) : ℕ) • P.Div φA := by
    rw [← hdivL, ← hdivR, hsq]
  have key : Φ.map (@inv _ _ _ _ (P.Base ζ) hbζ)
      (Φ.map (@inv _ _ _ _ (P.Base β') hbβ) (Φ.map (P.Base β') (P.Div ζ)))
      = preStepVal P ζ hζ := by
    rw [MonoidOn.map_map_inv Φ (P.Base β') hbβ]
    rfl
  rw [← key, hd]
  have hinvsq : @inv _ _ _ _ (P.Base ζ) hbζ ≫ @inv _ _ _ _ (P.Base β') hbβ
      = @inv _ _ _ _ (P.Base αn) hbα ≫ @inv _ _ _ _ (P.Base φA) hbφ := by
    haveI : IsIso (P.Base β' ≫ P.Base ζ) := IsIso.comp_isIso' hbβ hbζ
    have h1 : @inv _ _ _ _ (P.Base β' ≫ P.Base ζ) inferInstance
        = @inv _ _ _ _ (P.Base ζ) hbζ ≫ @inv _ _ _ _ (P.Base β') hbβ := by simp
    haveI : IsIso (P.Base φA ≫ P.Base αn) := IsIso.comp_isIso' hbφ hbα
    have h2 : @inv _ _ _ _ (P.Base φA ≫ P.Base αn) inferInstance
        = @inv _ _ _ _ (P.Base αn) hbα ≫ @inv _ _ _ _ (P.Base φA) hbφ := by simp
    rw [← h1, ← h2]
    congr 1
    rw [← P.Base_comp, ← P.Base_comp, hsq]
  have hL : Φ.map (@inv _ _ _ _ (P.Base ζ) hbζ)
      (Φ.map (@inv _ _ _ _ (P.Base β') hbβ) (((n : ℕ+) : ℕ) • P.Div φA))
      = Φ.map (@inv _ _ _ _ (P.Base αn) hbα)
        (Φ.map (@inv _ _ _ _ (P.Base φA) hbφ) (((n : ℕ+) : ℕ) • P.Div φA)) := by
    rw [← Φ.map_comp, ← Φ.map_comp, hinvsq]
  rw [hL, map_nsmul, map_inv_base_of_divIdentity αn hαD hbα]
  rfl

/-! ## ★4. ★★★★★`Proposition 4.1, (i)` -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★**[FrdI] Proposition 4.1, (i) の圏層** —— 図式の条件は単系の条件そのもの。 -/
theorem prop_4_1_i_bridge (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A : C} (α : ℕ+ → (A ⟶ A))
    (hαd : ∀ n, P.degFr (α n) = n) (hαD : ∀ n, IsDivIdentity P (α n))
    (hαF : ∀ n, IsFrobeniusType P (α n))
    {B : C} (φ : B ⟶ A) (hcφ : IsCoAngular P φ) (hst : IsStep P φ) :
    (∀ a b : Φ.val (P.toElem.obj A).base, preStepVal P φ hst.1 = a + b → a ≠ 0 → b ≠ 0 →
        MPrec (preStepVal P φ hst.1) a) ↔
      ∀ (B' : C) (φB : B ⟶ B') (φA : B' ⟶ A), IsStep P φB → IsStep P φA → φ = φB ≫ φA →
        ∃ (Y : C) (n : ℕ+) (β' : B' ⟶ Y) (ζ : Y ⟶ A) (ζ' : Y ⟶ B),
          IsFrobeniusType P β' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧
          ζ = ζ' ≫ φ ∧ β' ≫ ζ = φA ≫ α n := by
  constructor
  · intro hmono B' φB φA hsB hsA hfac
    have hval : preStepVal P φ hst.1
        = preStepVal P φA hsA.1 + Φ.map (@inv _ _ _ _ (P.Base φ) hst.1.2) (P.Div φB) := by
      subst hfac
      exact preStepVal_comp φB φA hsB.1 hsA.1 hst.1
    have ha0 : preStepVal P φA hsA.1 ≠ 0 :=
      (isStep_iff_preStepVal_ne_zero hiso φA hsA.1).mp hsA
    have hb0 : Φ.map (@inv _ _ _ _ (P.Base φ) hst.1.2) (P.Div φB) ≠ 0 := by
      intro h0
      exact ((isStep_iff_preStepVal_ne_zero hiso φB hsB.1).mp hsB)
        ((preStepVal_eq_zero_iff φB hsB.1).mpr
          (Φ.map_injective (@inv _ _ _ _ (P.Base φ) hst.1.2) (h0.trans (map_zero _).symm)))
    obtain ⟨n, hn, c, hc⟩ := hmono _ _ hval ha0 hb0
    obtain ⟨Y, β', ζ, hβ'F, -, hζs, hsq⟩ :=
      prop_1_10_ii P G.core φA hsA.1 (α ⟨n, hn⟩) (hαF ⟨n, hn⟩)
    have hζval : preStepVal P ζ hζs = ((⟨n, hn⟩ : ℕ+) : ℕ) • preStepVal P φA hsA.1 :=
      preStepVal_of_frob_square φA hsA.1 (α ⟨n, hn⟩) (hαF ⟨n, hn⟩) (hαD ⟨n, hn⟩)
        (hαd ⟨n, hn⟩) β' hβ'F ζ hζs hsq
    have hle : MLe (preStepVal P φ hst.1) (preStepVal P ζ hζs) := by
      rw [hζval]
      exact ⟨c, hc⟩
    obtain ⟨ζ', -, hζ's, hζeq⟩ := exists_factor_through G φ ζ hcφ hst.1
      (isCoAngular_of_isotropic_dom (P := P) G.core (hiso Y) ζ) hζs hle
    exact ⟨Y, ⟨n, hn⟩, β', ζ, ζ', hβ'F, hζs, hζ's, hζeq, hsq⟩
  · intro hcat a b hab ha0 hb0
    obtain ⟨B', ψ, χ, hψs, hχs, hψc, hχc, hfac, hχa⟩ :=
      exists_factor_of_mle' G φ hcφ hst.1 (⟨b, hab.symm⟩ : MLe a (preStepVal P φ hst.1))
    subst hfac
    have hval : preStepVal P (ψ ≫ χ) hst.1
        = preStepVal P χ hχs + Φ.map (@inv _ _ _ _ (P.Base (ψ ≫ χ)) hst.1.2) (P.Div ψ) :=
      preStepVal_comp ψ χ hψs hχs hst.1
    letI := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj A).base).1.1
    have hval' : preStepVal P (ψ ≫ χ) hst.1
        = a + Φ.map (@inv _ _ _ _ (P.Base (ψ ≫ χ)) hst.1.2) (P.Div ψ) := by
      rw [← hχa]; exact hval
    have hy : Φ.map (@inv _ _ _ _ (P.Base (ψ ≫ χ)) hst.1.2) (P.Div ψ) = b :=
      add_left_cancel (a := a) (hval'.symm.trans hab)
    have hχstep : IsStep P χ := (isStep_iff_preStepVal_ne_zero hiso χ hχs).mpr (hχa ▸ ha0)
    have hψstep : IsStep P ψ := by
      refine (isStep_iff_preStepVal_ne_zero hiso ψ hψs).mpr ?_
      rw [Ne, preStepVal_eq_zero_iff]
      intro h0
      exact hb0 (by rw [← hy, h0, map_zero])
    obtain ⟨Y, n, β', ζ, ζ', hβ'F, hζs, hζ's, hζeq, hsq⟩ :=
      hcat B' ψ χ hψstep hχstep rfl
    have hζval : preStepVal P ζ hζs = ((n : ℕ+) : ℕ) • preStepVal P χ hχs :=
      preStepVal_of_frob_square χ hχs (α n) (hαF n) (hαD n) (hαd n) β' hβ'F ζ hζs hsq
    have hle : MLe (preStepVal P (ψ ≫ χ) hst.1) (preStepVal P ζ hζs) := by
      subst hζeq
      exact ⟨Φ.map (@inv _ _ _ _ (P.Base (ζ' ≫ ψ ≫ χ)) hζs.2) (P.Div ζ'),
        (preStepVal_comp ζ' (ψ ≫ χ) hζ's hst.1 hζs).symm⟩
    exact ⟨((n : ℕ+) : ℕ), n.2, by rw [← hχa, ← hζval]; exact hle⟩

set_option maxHeartbeats 1000000 in
/-- ★★★★★**[FrdI] Proposition 4.1, (i)** —— **primary な step の圏論的特徴づけ**。

原文 (FrdI p.75):
> (i) φ is primary if and only if, for every factorization φ = φA ◦φB, where

★単系層(`isPrimaryElt_iff_of_perfect`、`Φ(A)` が perfect)と
圏層(`prop_4_1_i_bridge`)を繋ぐだけ。 -/
theorem prop_4_1_i (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A : C} (α : ℕ+ → (A ⟶ A))
    (hαd : ∀ n, P.degFr (α n) = n) (hαD : ∀ n, IsDivIdentity P (α n))
    (hαF : ∀ n, IsFrobeniusType P (α n))
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj A).base))
    {B : C} (φ : B ⟶ A) (hcφ : IsCoAngular P φ) (hst : IsStep P φ) :
    IsPrimaryElt (P.Div φ) ↔
      ∀ (B' : C) (φB : B ⟶ B') (φA : B' ⟶ A), IsStep P φB → IsStep P φA → φ = φB ≫ φA →
        ∃ (Y : C) (n : ℕ+) (β' : B' ⟶ Y) (ζ : Y ⟶ A) (ζ' : Y ⟶ B),
          IsFrobeniusType P β' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧
          ζ = ζ' ≫ φ ∧ β' ≫ ζ = φA ≫ α n :=
  (isPrimaryElt_preStepVal_iff φ hst.1).symm.trans
    ((isPrimaryElt_iff_of_perfect hperf
        ((isStep_iff_preStepVal_ne_zero hiso φ hst.1).mp hst)).trans
      (prop_4_1_i_bridge G hiso α hαd hαD hαF φ hcφ hst))

/-- ★★★**`Div-Frobenius-trivial` な対象は `α_n` の族を与える**。 -/
theorem prop_4_1_i_of_divFrobTrivial (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A : C} (hA : IsDivFrobeniusTrivial P A)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj A).base))
    {B : C} (φ : B ⟶ A) (hcφ : IsCoAngular P φ) (hst : IsStep P φ) :
    ∃ α : ℕ+ → (A ⟶ A), (∀ n, P.degFr (α n) = n) ∧ (∀ n, IsDivIdentity P (α n)) ∧
      (∀ n, IsFrobeniusType P (α n)) ∧
      (IsPrimaryElt (P.Div φ) ↔
        ∀ (B' : C) (φB : B ⟶ B') (φA : B' ⟶ A), IsStep P φB → IsStep P φA → φ = φB ≫ φA →
          ∃ (Y : C) (n : ℕ+) (β' : B' ⟶ Y) (ζ : Y ⟶ A) (ζ' : Y ⟶ B),
            IsFrobeniusType P β' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧
            ζ = ζ' ≫ φ ∧ β' ≫ ζ = φA ≫ α n) := by
  obtain ⟨ζ, hζd, hζp⟩ := hA
  exact ⟨fun n => ζ n, hζd, fun n => (hζp n).1, fun n => (hζp n).2,
    prop_4_1_i G hiso (fun n => ζ n) hζd (fun n => (hζp n).1) (fun n => (hζp n).2) hperf
      φ hcφ hst⟩


/-! ## ★4b. (ii) のための追加の道具 -/

/-- ★★★**零因子が等しい 2 本の co-angular pre-step は `B` の上で同型**。

★`Definition 1.3, (iii)(d)` の**前置**の圏同値(コスライス)から。
行き先が前順序圏なので射は高々 1 本であり、両向きの射が自動的に互いに逆になる。 -/
theorem exists_iso_of_div_eq (G : Frobenioid P) {B A A' : C} (φ : B ⟶ A) (φ' : B ⟶ A')
    (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    (hcφ' : IsCoAngular P φ') (hsφ' : IsPreStep P φ')
    (h : P.Div φ = P.Div φ') :
    ∃ θ : A ⟶ A', IsIso θ ∧ φ ≫ θ = φ' := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv B
  let Z : Under (⟨B⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨φ, hcφ, hsφ⟩)
  let W : Under (⟨B⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A'⟩ from ⟨φ', hcφ', hsφ'⟩)
  have hobj : (coaPreUnderFunctor P B).obj Z = (coaPreUnderFunctor P B).obj W := by
    show toOrderCat (P.Div φ) = toOrderCat (P.Div φ')
    rw [h]
  obtain ⟨f, -⟩ := (coaPreUnderFunctor P B).map_surjective (eqToHom hobj)
  obtain ⟨g, -⟩ := (coaPreUnderFunctor P B).map_surjective (eqToHom hobj.symm)
  have hfg : f ≫ g = 𝟙 Z :=
    (coaPreUnderFunctor P B).map_injective (Subsingleton.elim _ _)
  have hgf : g ≫ f = 𝟙 W :=
    (coaPreUnderFunctor P B).map_injective (Subsingleton.elim _ _)
  refine ⟨f.right.hom, ⟨g.right.hom, ?_, ?_⟩, ?_⟩
  · exact congrArg (fun t : Z ⟶ Z => t.right.hom) hfg
  · exact congrArg (fun t : W ⟶ W => t.right.hom) hgf
  · exact congrArg InducedWideCategory.Hom.hom (Under.w f)

/-- ★底の逆射の合成。 -/
theorem inv_base_comp {A B B' : C} (ψ : B ⟶ B') (χ : B' ⟶ A)
    (hψ : IsIso (P.Base ψ)) (hχ : IsIso (P.Base χ)) (hc : IsIso (P.Base (ψ ≫ χ))) :
    @inv _ _ _ _ (P.Base (ψ ≫ χ)) hc
      = @inv _ _ _ _ (P.Base χ) hχ ≫ @inv _ _ _ _ (P.Base ψ) hψ := by
  refine IsIso.inv_eq_of_hom_inv_id ?_
  rw [P.Base_comp, Category.assoc, ← Category.assoc (P.Base χ), IsIso.hom_inv_id,
    Category.id_comp]
  exact @IsIso.hom_inv_id _ _ _ _ (P.Base ψ) hψ

/-- ★★**第 1 因子の値は、`Φ(A)` の値を `ψ` で押したもの**。 -/
theorem yVal_eq {A B Cc : C} (φ : B ⟶ A) (ψ : A ⟶ Cc)
    (hsφ : IsPreStep P φ) (hsψ : IsPreStep P ψ) (hc : IsPreStep P (φ ≫ ψ)) :
    Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hc.2) (P.Div φ)
      = Φ.map (@inv _ _ _ _ (P.Base ψ) hsψ.2) (preStepVal P φ hsφ) := by
  rw [inv_base_comp φ ψ hsφ.2 hsψ.2 hc.2, Φ.map_comp]
  rfl

/-- ★第 1 因子の零因子は合成の零因子以下。 -/
theorem mle_div_comp {A B B' : C} (ψ : B ⟶ B') (χ : B' ⟶ A) (hlin : P.degFr χ = 1) :
    MLe (P.Div ψ) (P.Div (ψ ≫ χ)) := by
  refine ⟨Φ.map (P.Base ψ) (P.Div χ), ?_⟩
  rw [P.Div_comp, hlin]
  simp [add_comm]

theorem map_base_comp_apply {A B Cc : C} (φ : B ⟶ A) (ψ : A ⟶ Cc)
    (d : Φ.val (P.toElem.obj Cc).base) :
    Φ.map (P.Base φ) (Φ.map (P.Base ψ) d) = Φ.map (P.Base (φ ≫ ψ)) d := by
  rw [P.Base_comp, Φ.map_comp]

set_option maxHeartbeats 1000000 in
/-- ★★★**値 `dA` を持つ第 1 因子を取り出す**。 -/
theorem exists_first_factor_of_mle (G : Frobenioid P) {A B : C} (φ : B ⟶ A)
    (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    {dA : Φ.val (P.toElem.obj A).base} (hle : MLe dA (preStepVal P φ hsφ)) :
    ∃ (A'' : C) (φ'' : B ⟶ A'') (ζ : A'' ⟶ A), IsCoAngular P φ'' ∧ IsPreStep P φ'' ∧
      IsCoAngular P ζ ∧ IsPreStep P ζ ∧ φ = φ'' ≫ ζ ∧ P.Div φ'' = Φ.map (P.Base φ) dA := by
  letI := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj A).base).1.1
  obtain ⟨e, he⟩ := hle
  obtain ⟨A'', φ'', ζ, hsφ'', hsζ, hcφ'', hcζ, hfac, hζval⟩ :=
    exists_factor_of_mle' G φ hcφ hsφ (⟨dA, by rw [add_comm]; exact he⟩ : MLe e _)
  have key : ∀ (f : B ⟶ A) (hf : IsPreStep P f), f = φ'' ≫ ζ →
      preStepVal P f hf = e + Φ.map (@inv _ _ _ _ (P.Base f) hf.2) (P.Div φ'') := by
    rintro f hf rfl
    rw [← hζval]
    exact preStepVal_comp φ'' ζ hsφ'' hsζ hf
  have hv := key φ hsφ hfac
  have hda : Φ.map (@inv _ _ _ _ (P.Base φ) hsφ.2) (P.Div φ'') = dA := by
    refine add_left_cancel (a := e) ?_
    rw [← hv, ← he, add_comm]
  refine ⟨A'', φ'', ζ, hcφ'', hsφ'', hcζ, hsζ, hfac, ?_⟩
  rw [← hda, MonoidOn.map_inv_map Φ (P.Base φ) hsφ.2]

end Cat

/-! ## ★5. `Proposition 4.1, (ii)` の単系層 -/

/-- ★★★★**[FrdI] Proposition 4.1, (ii) の単系層**。

`q` が primary のとき、

  `p + q` が primary ⟺ どの分解 `p + q = a + b`(`a, b ≠ 0`)にも
  **`d ≼ q` かつ `d ≼ b` なる `d ≠ 0`** がある。

★★`⟹` は **perfect で `n` 割り**するだけ(`d := q/n`)——
`p+q` の primary 性から `p+q ≼ n·b` が出るので `q ≼ n·b`、
両辺を `n` で割ると `d ≼ b` になる。
★★`⟸` は ★**primary ⟺ 台が 1 点**(`isPrimaryElt_iff_suppElt_singleton`)を使う ——
台に別の素点 `P'` があれば `P'` だけを取り出す分解(`exists_split_suppElt`)が
条件を破る(`eq_zero_of_mle_of_suppElt_disjoint`)。 -/
theorem prop_4_1_ii_monoid {M : Type w'} [AddCommMonoid M] {ι : Prime M → Pf M → ℝ≥0}
    (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M)
    (hdiv : IsDivisorial M) {p q : M} (hq : IsPrimaryElt q) :
    IsPrimaryElt (p + q) ↔
      ∀ a b : M, p + q = a + b → a ≠ 0 → b ≠ 0 →
        ∃ d : M, d ≠ 0 ∧ MLe d q ∧ MLe d b := by
  have hqle : MLe q (p + q) := ⟨p, add_comm q p⟩
  constructor
  · intro hpq a b hab ha0 hb0
    obtain ⟨n, hn, c, hc⟩ :=
      hpq.2 b hb0 ⟨1, one_pos, a, by rw [one_smul, add_comm]; exact hab.symm⟩
    obtain ⟨d, hd⟩ := (hperf ⟨n, hn⟩).2 q
    have hnd : n • d = q := hd
    refine ⟨d, ?_, ?_, ?_⟩
    · intro h0
      exact hq.1 (by rw [← hnd, h0, smul_zero])
    · exact ⟨(n - 1) • d, by
        rw [← hnd, ← succ_nsmul']
        congr 1
        omega⟩
    · obtain ⟨e, he⟩ : MLe q (n • b) := mle_trans hqle ⟨c, hc⟩
      obtain ⟨e', he'⟩ := (hperf ⟨n, hn⟩).2 e
      have hne' : n • e' = e := he'
      refine ⟨e', nsmul_inj_of_divisorial hdiv hn ?_⟩
      rw [smul_add, hnd, hne']
      exact he
  · intro hcond
    obtain ⟨P, hP⟩ := suppElt_singleton_of_primary H hperf hdiv hq
    have hqsub : SuppElt ι q ⊆ SuppElt ι (p + q) := suppElt_mono_of_mle H hqle
    have hPmem : P ∈ SuppElt ι (p + q) := hqsub (by rw [hP]; rfl)
    have hsingle : SuppElt ι (p + q) = {P} := by
      refine Set.eq_singleton_iff_unique_mem.mpr ⟨hPmem, fun P' hP' => ?_⟩
      by_contra hne
      obtain ⟨y, z, hsum, hy, hz⟩ := exists_split_suppElt H hperf hdiv (p + q) {P'}ᶜ
      have hzsub : SuppElt ι z ⊆ {P'} := by
        rw [hz]
        intro x hx
        simpa using hx.1
      have hzP' : P' ∈ SuppElt ι z := by
        rw [hz]
        exact ⟨by simp, hP'⟩
      have hz0 : z ≠ 0 := by
        intro h0
        rw [h0, suppElt_zero H] at hzP'
        exact hzP'
      by_cases hy0 : y = 0
      · rw [hy0, suppElt_zero H] at hy
        have hmem : P ∈ ({P'}ᶜ : Set (Prime M)) ∩ SuppElt ι (p + q) := by
          refine ⟨?_, hPmem⟩
          simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
          exact fun h => hne h.symm
        rw [← hy] at hmem
        exact hmem
      · obtain ⟨d, hd0, hdq, hdz⟩ := hcond y z hsum hy0 hz0
        refine hd0 (eq_zero_of_mle_of_suppElt_disjoint H hdiv hdq hdz ?_)
        refine Set.eq_empty_iff_forall_notMem.mpr (fun x hx => ?_)
        have h1 : x = P := by rw [hP] at hx; exact hx.1
        have h2 : x = P' := by simpa using hzsub hx.2
        exact hne (h2.symm.trans h1)
    refine isPrimaryElt_of_suppElt_singleton H hdiv ?_ hsingle
    intro h0
    rw [h0, suppElt_zero H] at hsingle
    exact absurd hsingle.symm (by simp)

/-! ## ★6. ★★★★★`Proposition 4.1, (ii)` -/

section Cat2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

set_option maxHeartbeats 1000000 in
/-- ★★**(ii) の圏層** —— 圏の条件から単系の条件へ。 -/
theorem prop_4_1_ii_bridge_mpr (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B Cc : C} (φ : B ⟶ A) (ψ : A ⟶ Cc)
    (hstφ : IsStep P φ) (hstψ : IsStep P ψ)
    (hcat : ∀ (A' : C) (φ' : B ⟶ A') (ψ' : A' ⟶ Cc), IsStep P φ' → IsStep P ψ' →
          φ ≫ ψ = φ' ≫ ψ' →
        ∃ (A'' : C) (φ'' : B ⟶ A'') (ζ : A'' ⟶ A) (ζ' : A'' ⟶ A'),
          IsStep P φ'' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧ φ = φ'' ≫ ζ ∧ φ' = φ'' ≫ ζ')
    (a b : Φ.val (P.toElem.obj Cc).base)
    (hab : preStepVal P (φ ≫ ψ) (IsPreStep.comp P hstφ.1 hstψ.1) = a + b)
    (ha0 : a ≠ 0) (hb0 : b ≠ 0) :
    ∃ d : Φ.val (P.toElem.obj Cc).base, d ≠ 0 ∧
      MLe d (Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ))
        (IsPreStep.comp P hstφ.1 hstψ.1).2) (P.Div φ)) ∧ MLe d b := by
  have hsc : IsPreStep P (φ ≫ ψ) := IsPreStep.comp P hstφ.1 hstψ.1
  have hcc : IsCoAngular P (φ ≫ ψ) := isCoAngular_of_isotropic_dom (P := P) G.core (hiso B) _
  letI := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj Cc).base).1.1
  obtain ⟨A', u, v, hsu, hsv, hcu, hcv, hfac2, hva⟩ :=
    exists_factor_of_mle' G (φ ≫ ψ) hcc hsc (⟨b, hab.symm⟩ : MLe a _)
  have key : ∀ (f : B ⟶ Cc) (hf : IsPreStep P f), f = u ≫ v →
      preStepVal P f hf = a + Φ.map (@inv _ _ _ _ (P.Base f) hf.2) (P.Div u) := by
    rintro f hf rfl
    rw [← hva]
    exact preStepVal_comp u v hsu hsv hf
  have hval2 := key (φ ≫ ψ) hsc hfac2
  have hyu : Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (P.Div u) = b :=
    add_left_cancel (a := a) (hval2.symm.trans hab)
  have hstv : IsStep P v := (isStep_iff_preStepVal_ne_zero hiso v hsv).mpr (hva ▸ ha0)
  have hstu : IsStep P u := by
    refine (isStep_iff_preStepVal_ne_zero hiso u hsu).mpr ?_
    rw [Ne, preStepVal_eq_zero_iff]
    intro h0
    exact hb0 (by rw [← hyu, h0, map_zero])
  obtain ⟨A'', φ'', ζ, ζ', hstφ'', hζs, hζ's, hfacφ, hfacu⟩ := hcat A' u v hstu hstv hfac2
  refine ⟨Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (P.Div φ''), ?_, ?_, ?_⟩
  · intro h0
    refine ((isStep_iff_preStepVal_ne_zero hiso φ'' hstφ''.1).mp hstφ'') ?_
    rw [preStepVal_eq_zero_iff]
    exact Φ.map_injective (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (h0.trans (map_zero _).symm)
  · refine map_mle (Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2)) ?_
    rw [hfacφ]
    exact mle_div_comp φ'' ζ hζs.1
  · rw [← hyu]
    refine map_mle (Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2)) ?_
    rw [hfacu]
    exact mle_div_comp φ'' ζ' hζ's.1

set_option maxHeartbeats 1000000 in
/-- ★★**(ii) の圏層** —— 単系の条件から圏の条件へ。

★★2 つの第 1 因子 `φ''`(`φ` 側)と `φ'''`(`φ'` 側)は**零因子が等しい**ので、
`exists_iso_of_div_eq`(前置の圏同値)で同型になり、1 本にまとまる。 -/
theorem prop_4_1_ii_bridge_mp (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B Cc : C} (φ : B ⟶ A) (ψ : A ⟶ Cc)
    (hstφ : IsStep P φ) (hstψ : IsStep P ψ)
    (hmono : ∀ a b : Φ.val (P.toElem.obj Cc).base,
        preStepVal P (φ ≫ ψ) (IsPreStep.comp P hstφ.1 hstψ.1) = a + b → a ≠ 0 → b ≠ 0 →
        ∃ d : Φ.val (P.toElem.obj Cc).base, d ≠ 0 ∧
          MLe d (Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ))
            (IsPreStep.comp P hstφ.1 hstψ.1).2) (P.Div φ)) ∧ MLe d b)
    (A' : C) (φ' : B ⟶ A') (ψ' : A' ⟶ Cc) (hstφ' : IsStep P φ') (hstψ' : IsStep P ψ')
    (hfac : φ ≫ ψ = φ' ≫ ψ') :
    ∃ (A'' : C) (φ'' : B ⟶ A'') (ζ : A'' ⟶ A) (ζ' : A'' ⟶ A'),
      IsStep P φ'' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧ φ = φ'' ≫ ζ ∧ φ' = φ'' ≫ ζ' := by
  have hsc : IsPreStep P (φ ≫ ψ) := IsPreStep.comp P hstφ.1 hstψ.1
  have hcφ : IsCoAngular P φ := isCoAngular_of_isotropic_dom (P := P) G.core (hiso B) φ
  have hcφ' : IsCoAngular P φ' := isCoAngular_of_isotropic_dom (P := P) G.core (hiso B) φ'
  have key2 : ∀ (f : B ⟶ Cc) (hf : IsPreStep P f), f = φ' ≫ ψ' →
      preStepVal P f hf = preStepVal P ψ' hstψ'.1
        + Φ.map (@inv _ _ _ _ (P.Base f) hf.2) (P.Div φ') := by
    rintro f hf rfl
    exact preStepVal_comp φ' ψ' hstφ'.1 hstψ'.1 hf
  have hval2 := key2 (φ ≫ ψ) hsc hfac
  have ha0 : preStepVal P ψ' hstψ'.1 ≠ 0 :=
    (isStep_iff_preStepVal_ne_zero hiso ψ' hstψ'.1).mp hstψ'
  have hb0 : Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (P.Div φ') ≠ 0 := by
    intro h0
    refine ((isStep_iff_preStepVal_ne_zero hiso φ' hstφ'.1).mp hstφ') ?_
    rw [preStepVal_eq_zero_iff]
    exact Φ.map_injective (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (h0.trans (map_zero _).symm)
  obtain ⟨d, hd0, hdq, hdb⟩ := hmono _ _ hval2 ha0 hb0
  rw [yVal_eq φ ψ hstφ.1 hstψ.1 hsc] at hdq
  have hdA : MLe (Φ.map (P.Base ψ) d) (preStepVal P φ hstφ.1) := by
    have h := map_mle (Φ.map (P.Base ψ)) hdq
    rwa [MonoidOn.map_inv_map Φ (P.Base ψ) hstψ.1.2] at h
  obtain ⟨A'', φ'', ζ, hcφ'', hsφ'', hcζ, hsζ, hfacφ, hdivφ''⟩ :=
    exists_first_factor_of_mle G φ hcφ hstφ.1 hdA
  have key3 : ∀ (f : B ⟶ Cc) (hf : IsPreStep P f), f = φ' ≫ ψ' →
      Φ.map (@inv _ _ _ _ (P.Base f) hf.2) (P.Div φ')
        = Φ.map (@inv _ _ _ _ (P.Base ψ') hstψ'.1.2) (preStepVal P φ' hstφ'.1) := by
    rintro f hf rfl
    exact yVal_eq φ' ψ' hstφ'.1 hstψ'.1 hf
  rw [key3 (φ ≫ ψ) hsc hfac] at hdb
  have hdA' : MLe (Φ.map (P.Base ψ') d) (preStepVal P φ' hstφ'.1) := by
    have h := map_mle (Φ.map (P.Base ψ')) hdb
    rwa [MonoidOn.map_inv_map Φ (P.Base ψ') hstψ'.1.2] at h
  obtain ⟨A''', φ''', ζ', hcφ''', hsφ''', hcζ', hsζ', hfacφ', hdivφ'''⟩ :=
    exists_first_factor_of_mle G φ' hcφ' hstφ'.1 hdA'
  have hdiveq : P.Div φ'' = P.Div φ''' := by
    rw [hdivφ'', hdivφ''', map_base_comp_apply, map_base_comp_apply, hfac]
  obtain ⟨θ, hθiso, hθ⟩ := exists_iso_of_div_eq G φ'' φ''' hcφ'' hsφ'' hcφ''' hsφ''' hdiveq
  haveI := hθiso
  refine ⟨A'', φ'', ζ, θ ≫ ζ', ?_, hsζ,
    IsPreStep.comp P (isPreStep_of_isIso P θ) hsζ', hfacφ, ?_⟩
  · refine (isStep_iff_preStepVal_ne_zero hiso φ'' hsφ'').mpr ?_
    rw [Ne, preStepVal_eq_zero_iff, hdivφ'']
    intro h0
    refine hd0 ?_
    refine Φ.map_injective (P.Base ψ) ?_
    refine Φ.map_injective (P.Base φ) ?_
    rw [map_zero, map_zero]
    exact h0
  · rw [hfacφ', ← hθ, Category.assoc]

set_option maxHeartbeats 1000000 in
/-- ★★★★★**[FrdI] Proposition 4.1, (ii)**。

原文 (FrdI p.75):
> (ii) Suppose that φ is primary. Then the composite ψ ◦φ : B →C, hence

★単系層(`prop_4_1_ii_monoid`)と圏層(`prop_4_1_ii_bridge_*`)を繋ぐだけ。 -/
theorem prop_4_1_ii (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B Cc : C} (φ : B ⟶ A) (ψ : A ⟶ Cc)
    (hstφ : IsStep P φ) (hstψ : IsStep P ψ)
    {ι : Prime (Φ.val (P.toElem.obj Cc).base) → Pf (Φ.val (P.toElem.obj Cc).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj Cc).base) ι)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj Cc).base))
    (hprim : IsPrimaryElt (P.Div φ)) :
    IsPrimaryElt (P.Div (φ ≫ ψ)) ↔
      (∀ (A' : C) (φ' : B ⟶ A') (ψ' : A' ⟶ Cc), IsStep P φ' → IsStep P ψ' →
          φ ≫ ψ = φ' ≫ ψ' →
        ∃ (A'' : C) (φ'' : B ⟶ A'') (ζ : A'' ⟶ A) (ζ' : A'' ⟶ A'),
          IsStep P φ'' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧ φ = φ'' ≫ ζ ∧ φ' = φ'' ≫ ζ') := by
  have hsc : IsPreStep P (φ ≫ ψ) := IsPreStep.comp P hstφ.1 hstψ.1
  have hval : preStepVal P (φ ≫ ψ) hsc
      = preStepVal P ψ hstψ.1 + Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (P.Div φ) :=
    preStepVal_comp φ ψ hstφ.1 hstψ.1 hsc
  have hq : IsPrimaryElt (Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (P.Div φ)) := by
    rw [yVal_eq φ ψ hstφ.1 hstψ.1 hsc]
    exact isPrimaryElt_of_bijective (Φ.map (@inv _ _ _ _ (P.Base ψ) hstψ.1.2))
      (Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base ψ) hstψ.1.2).symm)
      ((isPrimaryElt_preStepVal_iff φ hstφ.1).mpr hprim)
  have hmonoiff : IsPrimaryElt (preStepVal P (φ ≫ ψ) hsc) ↔
      (∀ a b : Φ.val (P.toElem.obj Cc).base, preStepVal P (φ ≫ ψ) hsc = a + b →
        a ≠ 0 → b ≠ 0 → ∃ d : Φ.val (P.toElem.obj Cc).base, d ≠ 0 ∧
          MLe d (Φ.map (@inv _ _ _ _ (P.Base (φ ≫ ψ)) hsc.2) (P.Div φ)) ∧ MLe d b) := by
    rw [hval]
    exact prop_4_1_ii_monoid H hperf (P.divisorial _) hq
  refine (isPrimaryElt_preStepVal_iff (φ ≫ ψ) hsc).symm.trans (hmonoiff.trans ⟨?_, ?_⟩)
  · intro hm A' φ' ψ' h1 h2 h3
    exact prop_4_1_ii_bridge_mp G hiso φ ψ hstφ hstψ hm A' φ' ψ' h1 h2 h3
  · intro hc a b hab ha0 hb0
    exact prop_4_1_ii_bridge_mpr G hiso φ ψ hstφ hstψ hc a b hab ha0 hb0

end Cat2


/-! ## ★7. ★★★★★`Proposition 4.1, (iii)` -/

section Cat3

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★`ε = ε' ≫ ζ` なら `x_ζ ≼ x_ε`。 -/
theorem mle_preStepVal_of_factor {E Z F : C} (ε : E ⟶ F) (ε' : E ⟶ Z) (ζ : Z ⟶ F)
    (hsε : IsPreStep P ε) (hsε' : IsPreStep P ε') (hsζ : IsPreStep P ζ)
    (hfac : ε = ε' ≫ ζ) : MLe (preStepVal P ζ hsζ) (preStepVal P ε hsε) := by
  have key : ∀ (f : E ⟶ F) (hf : IsPreStep P f), f = ε' ≫ ζ →
      preStepVal P f hf = preStepVal P ζ hsζ
        + Φ.map (@inv _ _ _ _ (P.Base f) hf.2) (P.Div ε') := by
    rintro f hf rfl
    exact preStepVal_comp ε' ζ hsε' hsζ hf
  exact ⟨_, (key ε hsε hfac).symm⟩

/-- ★★**任意の値を持つ co-angular pre-step が取れる**(後置の圏同値の本質的全射性)。 -/
theorem exists_preStep_of_val (G : Frobenioid P) (F : C)
    (a : Φ.val (P.toElem.obj F).base) :
    ∃ (U : C) (w : U ⟶ F) (hw : IsPreStep P w), IsCoAngular P w ∧ preStepVal P w hw = a := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv F
  let W := (coaPreOverFunctor P F).objPreimage (Opposite.op (toOrderCat a))
  have hiso2 : (coaPreOverFunctor P F).obj W ≅ Opposite.op (toOrderCat a) :=
    (coaPreOverFunctor P F).objObjPreimageIso _
  haveI hWb : IsIso (P.Base W.hom.hom) := W.hom.property.2.2
  exact ⟨W.left.obj, W.hom.hom, W.hom.property.2, W.hom.property.1,
    mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2
      (leOfHom hiso2.inv.unop) (leOfHom hiso2.hom.unop)⟩

set_option maxHeartbeats 1000000 in
/-- ★★★★★**[FrdI] Proposition 4.1, (iii) の前半** ——
**台が交わらない ⟺ 共通の pre-step はすべて同型**。

原文 (FrdI p.75):
> (iii) ϵ∗(Div(ϵ)), ι∗(Div(ι)) ∈Φ(F) [where we write ϵ∗, ι∗for the respective

★単系層は `suppElt_disjoint_iff`、圏層は後置の圏同値の充満性と本質的全射性。
★「同型 ⟺ 値が 0」は isotropic 型から。 -/
theorem prop_4_1_iii (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {E I F : C} (ε : E ⟶ F) (κ : I ⟶ F) (hsε : IsPreStep P ε) (hsκ : IsPreStep P κ)
    {ι : Prime (Φ.val (P.toElem.obj F).base) → Pf (Φ.val (P.toElem.obj F).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj F).base) ι)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj F).base)) :
    SuppElt ι (preStepVal P ε hsε) ∩ SuppElt ι (preStepVal P κ hsκ) = ∅ ↔
      ∀ (Z : C) (ζ : Z ⟶ F), IsPreStep P ζ →
        (∃ ε' : E ⟶ Z, IsPreStep P ε' ∧ ε = ε' ≫ ζ) →
        (∃ κ' : I ⟶ Z, IsPreStep P κ' ∧ κ = κ' ≫ ζ) → IsIso ζ := by
  rw [suppElt_disjoint_iff H hperf (P.divisorial _)]
  constructor
  · intro hmono Z ζ hsζ ⟨ε', hsε', hfacε⟩ ⟨κ', hsκ', hfacκ⟩
    refine hiso Z F ζ ?_ hsζ
    show P.Div ζ = 0
    rw [← preStepVal_eq_zero_iff ζ hsζ]
    exact hmono _ (mle_preStepVal_of_factor ε ε' ζ hsε hsε' hsζ hfacε)
      (mle_preStepVal_of_factor κ κ' ζ hsκ hsκ' hsζ hfacκ)
  · intro hcat d hdε hdκ
    have hcε : IsCoAngular P ε := isCoAngular_of_isotropic_dom (P := P) G.core (hiso E) ε
    obtain ⟨Z, ε', ζ, hsε', hsζ, hcε', hcζ, hfacε, hζval⟩ :=
      exists_factor_of_mle' G ε hcε hsε hdε
    have hdζ : MLe (preStepVal P ζ hsζ) (preStepVal P κ hsκ) := by
      rw [hζval]; exact hdκ
    have hcκ : IsCoAngular P κ := isCoAngular_of_isotropic_dom (P := P) G.core (hiso I) κ
    obtain ⟨κ', -, hsκ', hfacκ⟩ := exists_factor_through G ζ κ hcζ hsζ hcκ hsκ hdζ
    haveI := hcat Z ζ hsζ ⟨ε', hsε', hfacε⟩ ⟨κ', hsκ', hfacκ⟩
    rw [← hζval, preStepVal_eq_zero_iff]
    exact isIsometric_of_isIso P ζ

set_option maxHeartbeats 1000000 in
/-- ★★★★**[FrdI] Proposition 4.1, (iii) の後半** —— **pre-step の圏での四角形**。

原文 (FrdI p.75):
> ι satisfying ϵ = ζ ◦ϵ′, ι = ζ ◦ι′ is, in fact, an isomorphism. In this case, we shall

★値は `x_ε + x_κ` に取る。★原文の 2 本の等式
`ε∗(ε′∗(Div ε′)) = ι∗(Div ι)`・`ι∗(ι′∗(Div ι′)) = ε∗(Div ε)` はその和の分解である。 -/
theorem prop_4_1_iii_square (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {E I F : C} (ε : E ⟶ F) (κ : I ⟶ F) (hsε : IsPreStep P ε) (hsκ : IsPreStep P κ) :
    ∃ (U : C) (ε' : U ⟶ E) (κ' : U ⟶ I) (_ : IsPreStep P ε') (_ : IsPreStep P κ')
      (_ : ε' ≫ ε = κ' ≫ κ) (hw : IsPreStep P (ε' ≫ ε)),
      Φ.map (@inv _ _ _ _ (P.Base (ε' ≫ ε)) hw.2) (P.Div ε') = preStepVal P κ hsκ ∧
      Φ.map (@inv _ _ _ _ (P.Base (ε' ≫ ε)) hw.2) (P.Div κ') = preStepVal P ε hsε := by
  letI := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj F).base).1.1
  obtain ⟨U, w, hsw, hcw, hwval⟩ :=
    exists_preStep_of_val G F (preStepVal P ε hsε + preStepVal P κ hsκ)
  have hcε : IsCoAngular P ε := isCoAngular_of_isotropic_dom (P := P) G.core (hiso E) ε
  have hcκ : IsCoAngular P κ := isCoAngular_of_isotropic_dom (P := P) G.core (hiso I) κ
  obtain ⟨ε', hcε', hsε', hfacε⟩ := exists_factor_through G ε w hcε hsε hcw hsw
    (by rw [hwval]; exact ⟨preStepVal P κ hsκ, rfl⟩)
  obtain ⟨κ', hcκ', hsκ', hfacκ⟩ := exists_factor_through G κ w hcκ hsκ hcw hsw
    (by rw [hwval]; exact ⟨preStepVal P ε hsε, add_comm _ _⟩)
  subst hfacε
  have h1 : preStepVal P (ε' ≫ ε) hsw = preStepVal P ε hsε
      + Φ.map (@inv _ _ _ _ (P.Base (ε' ≫ ε)) hsw.2) (P.Div ε') :=
    preStepVal_comp ε' ε hsε' hsε hsw
  rw [hwval] at h1
  have h2 : Φ.map (@inv _ _ _ _ (P.Base (ε' ≫ ε)) hsw.2) (P.Div ε') = preStepVal P κ hsκ :=
    (add_left_cancel (a := preStepVal P ε hsε) h1).symm
  have key : ∀ (f : U ⟶ F) (hf : IsPreStep P f), f = κ' ≫ κ →
      preStepVal P f hf = preStepVal P κ hsκ
        + Φ.map (@inv _ _ _ _ (P.Base f) hf.2) (P.Div κ') := by
    rintro f hf rfl
    exact preStepVal_comp κ' κ hsκ' hsκ hf
  have h3 := key (ε' ≫ ε) hsw hfacκ
  rw [hwval, add_comm (preStepVal P ε hsε)] at h3
  have h4 : Φ.map (@inv _ _ _ _ (P.Base (ε' ≫ ε)) hsw.2) (P.Div κ') = preStepVal P ε hsε :=
    (add_left_cancel (a := preStepVal P κ hsκ) h3).symm
  exact ⟨U, ε', κ', hsε', hsκ', hfacκ, hsw, h2, h4⟩

end Cat3


/-! ## ★出典の紐付け(条つき) -/

def prop_4_1_ii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 75, item := "Proposition 4.1, (ii)",
    sectionId := "frdi-prop-4-1" }

def prop_4_1_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 75, item := "Proposition 4.1, (i)",
    sectionId := "frdi-prop-4-1" }

end ABC3.Found.FrdI
