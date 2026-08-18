/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop41
import ABC3.Found.FrdI.Prop110
import ABC3.Found.FrdI.Prop32Equiv

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

end Cat

/-! ## ★出典の紐付け(条つき) -/

def prop_4_1_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 75, item := "Proposition 4.1, (i)",
    sectionId := "frdi-prop-4-1" }

end ABC3.Found.FrdI
