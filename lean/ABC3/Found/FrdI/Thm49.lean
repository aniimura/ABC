/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42Order

/-!
# [FrdI] Theorem 4.9 —— 除数単系の圏論性(`Ψ_Φ`)

原文 (FrdI p.88):
> Theorem 4.9. (Category-theoreticity of Divisor Monoids) For i = 1, 2, let

## ★★★★★★測って分かったこと —— **最難関の段を既に持っている**

原文の証明は「`Theorem 4.9` を示すには **`Theorem 4.2, (iii)` の
right-hand と left-hand の同型が一致する**ことを示せば十分」と言い、
そこに twin-primary steps の議論(原文 p.89-90)を費やす。

★★★**我々の `Theorem 4.2, (iii)` は最初から 1 本しかない** ——
`divEquiv` は `𝒪^▷(A)` の**自己射**の `Div` で作るので、
「`A` から出る」と「`A` へ入る」の区別が無い。
したがって **right-hand ＝ left-hand は自明**である(逸脱 (B) `hdivS` のおかげ)。

## ★★★★★残るのは**関手性**、そしてそれは 1 本の計算で出る

`f : W ⟶ A` と `u ∈ 𝒪^▷(A)`、`w ∈ 𝒪^▷(W)` が四角形
`w ≫ f = f ≫ u` を満たすとき、在庫の `div_square_frob` が
```
Φ.map (Base f) (Div u) = (degFr f) • Div w
```
を与える。★**`𝒞₂` 側でも同じ計算**をして、`divEquiv` が加法的
(したがって `ℕ` 倍と交換する)ことを使えば閉じる。
★★**場合分けが要らない** —— `div_square_frob` は `f` の型を問わない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section PsiPhi

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★`divMap` は `𝒪^▷` の `Div` を `Ψ` の像の `Div` へ送る。 -/
theorem divMap_div_otri (Ψ : C ≌ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    {A : C} {u : End A} (hu : u ∈ OTri P A) :
    divMap Ψ G hiso hdivS hOTri A (P.Div (((u : End A)) : A ⟶ A))
      = P₂.Div (Ψ.functor.map (((u : End A)) : A ⟶ A)) := by
  show P₂.Div (Ψ.functor.map ((((realizeDiv P hdivS A
    (P.Div (((u : End A)) : A ⟶ A))) : End A) : A ⟶ A))) = _
  exact div_map_eq_of_div_eq Ψ G hiso
    (realizeDiv_mem hdivS A (P.Div (((u : End A)) : A ⟶ A))) hu (by rw [realizeDiv_div])

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`Ψ_Φ` の自然性(四角形が与えられた場合)** ——
`div_square_frob` を両側で当てるだけ。★**`f` の型を問わない**。 -/
theorem divMap_naturality_of_square (Ψ : C ≌ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {W A : C} (f : W ⟶ A) {u : End A} (hu : u ∈ OTri P A) {ww : End W}
    (hw : ww ∈ OTri P W)
    (hsq : (((ww : End W)) : W ⟶ W) ≫ f = f ≫ (((u : End A)) : A ⟶ A)) :
    divMap Ψ G hiso hdivS hOTri W (Φ.map (P.Base f) (P.Div (((u : End A)) : A ⟶ A)))
      = Φ₂.map (P₂.Base (Ψ.functor.map f))
          (divMap Ψ G hiso hdivS hOTri A (P.Div (((u : End A)) : A ⟶ A))) := by
  -- ★左辺
  have hL : Φ.map (P.Base f) (P.Div (((u : End A)) : A ⟶ A))
      = ((P.degFr f : ℕ+) : ℕ) • P.Div (((ww : End W)) : W ⟶ W) :=
    div_square_frob P f hw hu hsq
  -- ★右辺の四角形(`𝒞₂` 側)
  have hsq₂ : ((Ψ.functor.map (((ww : End W)) : W ⟶ W)) : Ψ.functor.obj W ⟶ Ψ.functor.obj W)
      ≫ Ψ.functor.map f
      = Ψ.functor.map f
        ≫ ((Ψ.functor.map (((u : End A)) : A ⟶ A)) : Ψ.functor.obj A ⟶ Ψ.functor.obj A) := by
    rw [← Ψ.functor.map_comp, ← Ψ.functor.map_comp, hsq]
  have hR : Φ₂.map (P₂.Base (Ψ.functor.map f))
      (P₂.Div (Ψ.functor.map (((u : End A)) : A ⟶ A)))
      = ((P₂.degFr (Ψ.functor.map f) : ℕ+) : ℕ)
        • P₂.Div (Ψ.functor.map (((ww : End W)) : W ⟶ W)) :=
    div_square_frob P₂ (Ψ.functor.map f) (hOTri W ww hw) (hOTri A u hu) hsq₂
  rw [hL, map_nsmul, divMap_div_otri Ψ G hiso hdivS hOTri hw,
    divMap_div_otri Ψ G hiso hdivS hOTri hu, hR, hdeg f]

def divMap_naturality_of_square.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Theorem 4.9 — Ψ_Φ の自然性(四角形が与えられた場合)",
    sectionId := "frdi-thm-4-9" }

/-! ## ★2. `𝒪^▷` の引き戻し —— 3 つの射の類型 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★**Frobenius 型に沿った `𝒪^▷` の引き戻し**。

★★`Thm42DivId.lean` の `isDivIdentity_of_divIdCat` の中で組んだ構成を
**再利用できる補題として取り出したもの**である。

★手順: `Φ.map (Base γ) (Div v)` を `degFr γ` で割り(**perfect**)、
`hdivS`(逸脱 (B))で `𝒪^▷(A)` の元 `v₃` として実現し、
`Proposition 1.10, (i)` で四角形を立て、
`Base γ` が同型なので `Div u₀ = Div v` が出る。 -/
theorem otriPull_frobType (F : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    {A X : C} (γ : A ⟶ X) (hγ : IsFrobeniusType P γ)
    {v : End X} (hv : v ∈ OTri P X) :
    ∃ (v₃ : End A) (u₀ : End X), v₃ ∈ OTri P A ∧ u₀ ∈ OTri P X ∧
      P.Div (((u₀ : End X)) : X ⟶ X) = P.Div (((v : End X)) : X ⟶ X) ∧
      (((v₃ : End A)) : A ⟶ A) ≫ γ = γ ≫ (((u₀ : End X)) : X ⟶ X) := by
  haveI hγi : IsIso (P.Base γ) := hγ.2
  obtain ⟨d, hdd⟩ := (hperfM A (P.degFr γ)).2
    (Φ.map (P.Base γ) (P.Div ((((v : End X)) : X ⟶ X))))
  have hddb : ((P.degFr γ : ℕ+) : ℕ) • d
      = Φ.map (P.Base γ) (P.Div ((((v : End X)) : X ⟶ X))) := hdd
  obtain ⟨v₃, hv₃d⟩ := hdivS A d
  have hv₃ : (v₃ : End A) ∈ OTri P A := v₃.2
  obtain ⟨u₀, hsq0, -⟩ :=
    prop_1_10_i_exists_given P F ((((v₃ : End A)) : A ⟶ A)) γ hγ γ hγ rfl
  have hb0 : P.Base ((((v₃ : End A)) : A ⟶ A)) ≫ P.Base γ = P.Base γ ≫ P.Base u₀ := by
    rw [← P.Base_comp, ← P.Base_comp, hsq0]
  have hu₀b : IsBaseIdentity P u₀ := by
    show P.Base u₀ = P.Base (𝟙 X)
    rw [P.Base_id]
    rw [show P.Base ((((v₃ : End A)) : A ⟶ A)) = P.Base (𝟙 A) from hv₃.1, P.Base_id,
      Category.id_comp] at hb0
    exact ((cancel_epi (P.Base γ)).mp (by rw [Category.comp_id]; exact hb0)).symm
  have hu₀l : IsLinear P u₀ := by
    have hdd2 : P.degFr ((((v₃ : End A)) : A ⟶ A) ≫ γ) = P.degFr (γ ≫ u₀) := by rw [hsq0]
    rw [P.degFr_comp, P.degFr_comp,
      show P.degFr ((((v₃ : End A)) : A ⟶ A)) = 1 from hv₃.2, mul_one] at hdd2
    show P.degFr u₀ = 1
    exact (mul_right_cancel (b := P.degFr γ) (by rw [one_mul]; exact hdd2)).symm
  have hu₀m : u₀ ∈ OTri P X := ⟨hu₀b, hu₀l⟩
  have hdiv0 : Φ.map (P.Base γ) (P.Div u₀)
      = ((P.degFr γ : ℕ+) : ℕ) • P.Div ((((v₃ : End A)) : A ⟶ A)) :=
    div_square_frob P γ hv₃ hu₀m hsq0
  refine ⟨(v₃ : End A), u₀, hv₃, hu₀m, ?_, hsq0⟩
  refine (Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base γ) hγi)).1 ?_
  show Φ.map (P.Base γ) (P.Div u₀)
    = Φ.map (P.Base γ) (P.Div ((((v : End X)) : X ⟶ X)))
  rw [hdiv0, hv₃d]
  exact hdd

def otriPull_frobType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Theorem 4.9 — Frobenius 型に沿った 𝒪^▷ の引き戻し",
    sectionId := "frdi-thm-4-9" }

end PsiPhi

end ABC3.Found.FrdI
