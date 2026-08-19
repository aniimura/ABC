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

end PsiPhi

end ABC3.Found.FrdI
