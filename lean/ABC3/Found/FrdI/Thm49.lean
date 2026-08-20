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
theorem otriPull_frobType (F : FrobenioidCore P)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    {A X : C} (γ : A ⟶ X) (hγ : IsFrobeniusType P γ)
    {v : End X} :
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

/-! ## ★3. 自然性 —— 3 つの類型から任意の射へ -/

variable (Ψ : C ≌ C₂) (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
  (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
    ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
  (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
    ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
      ∈ OTri P₂ (Ψ.functor.obj Z))

/-- ★★`Ψ_Φ` の自然性(1 本の射について)。 -/
def DivMapNat {W A : C} (f : W ⟶ A) : Prop :=
  ∀ x : Φ.val (P.toElem.obj A).base,
    divMap Ψ G hiso hdivS hOTri W (Φ.map (P.Base f) x)
      = Φ₂.map (P₂.Base (Ψ.functor.map f)) (divMap Ψ G hiso hdivS hOTri A x)

variable {Ψ G hiso hdivS hOTri}

/-- ★自然性は合成で閉じる。 -/
theorem divMapNat_comp {W V A : C} {g : W ⟶ V} {h : V ⟶ A}
    (hg : DivMapNat Ψ G hiso hdivS hOTri g) (hh : DivMapNat Ψ G hiso hdivS hOTri h) :
    DivMapNat Ψ G hiso hdivS hOTri (g ≫ h) := by
  intro x
  have h1 : Φ.map (P.Base (g ≫ h)) x = Φ.map (P.Base g) (Φ.map (P.Base h) x) := by
    rw [P.Base_comp, Φ.map_comp]
  rw [h1, hg (Φ.map (P.Base h) x), hh x, ← Φ₂.map_comp, ← P₂.Base_comp,
    ← Ψ.functor.map_comp]

/-- ★四角形があれば自然性が出る。 -/
theorem divMapNat_of_squares
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {W A : C} (f : W ⟶ A)
    (hsq : ∀ (v : End A), v ∈ OTri P A →
      ∃ (ww : End W) (u : End A), ww ∈ OTri P W ∧ u ∈ OTri P A ∧
        P.Div (((u : End A)) : A ⟶ A) = P.Div (((v : End A)) : A ⟶ A) ∧
        (((ww : End W)) : W ⟶ W) ≫ f = f ≫ (((u : End A)) : A ⟶ A)) :
    DivMapNat Ψ G hiso hdivS hOTri f := by
  intro x
  obtain ⟨v, hvx⟩ := hdivS A x
  subst hvx
  obtain ⟨ww, u, hw, hu, hdu, hsqf⟩ := hsq (v : End A) v.2
  rw [← hdu]
  exact divMap_naturality_of_square Ψ G hiso hdivS hOTri hdeg f hu hw hsqf

/-- ★pull-back の場合。 -/
theorem divMapNat_pullBack
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {W A : C} (f : W ⟶ A) (hf : IsPullBack P f) :
    DivMapNat Ψ G hiso hdivS hOTri f :=
  divMapNat_of_squares hdeg f (fun v hv => by
    obtain ⟨ww, hw, hsq⟩ := otriPullBack P f hf hv
    exact ⟨ww, v, hw, hv, rfl, hsq.symm⟩)

/-- ★co-angular pre-step の場合。 -/
theorem divMapNat_preStep (F : FrobenioidCore P)
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {W A : C} (f : W ⟶ A) (hf : IsPreStep P f) :
    DivMapNat Ψ G hiso hdivS hOTri f :=
  divMapNat_of_squares hdeg f (fun v hv => by
    have hfc : IsCoAngular P f := prop_1_4_i P _ (fun Z _ => hiso Z)
    exact ⟨((otriPull P F f hfc hf.1 ⟨v, hv⟩ : End W)),
      v, (otriPull P F f hfc hf.1 ⟨v, hv⟩).2, hv, rfl,
      (otriPull_spec P F f hfc hf.1 ⟨v, hv⟩).symm⟩)

/-- ★Frobenius 型の場合。 -/
theorem divMapNat_frobType (F : FrobenioidCore P)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {W A : C} (f : W ⟶ A) (hf : IsFrobeniusType P f) :
    DivMapNat Ψ G hiso hdivS hOTri f :=
  divMapNat_of_squares hdeg f (fun v hv => by
    obtain ⟨v₃, u₀, hv₃, hu₀, hdu, hsq⟩ :=
      otriPull_frobType F hdivS hperfM f hf (v := v)
    exact ⟨v₃, u₀, hv₃, hu₀, hdu, hsq⟩)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Theorem 4.9 の核** —— `Ψ_Φ` は**任意の射**について自然。

★`Definition 1.3, (i), (d)` の 3 分解(`arbFactor`)で 3 つの類型に落とし、
それぞれの四角形を当てる。★合成で閉じる(`divMapNat_comp`)ので貼り合わせられる。 -/
theorem divMapNat_all (F : FrobenioidCore P)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {W A : C} (f : W ⟶ A) :
    DivMapNat Ψ G hiso hdivS hOTri f := by
  obtain ⟨Y₀, Z₀, γ, β, α₀, hfac, hγ, hβ, hα₀⟩ := F.arbFactor f
  rw [hfac]
  exact divMapNat_comp (divMapNat_frobType F hperfM hdeg γ hγ)
    (divMapNat_comp (divMapNat_preStep F hdeg β hβ)
      (divMapNat_pullBack hdeg α₀ hα₀))

def divMapNat_all.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Theorem 4.9 — Ψ_Φ の自然性",
    sectionId := "frdi-thm-4-9" }

/-! ## ★4. **`Ψ` は Div-equivalence を保つ** —— `Corollary 4.11, (ii)` が要る形 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Theorem 4.9 の系** —— `Ψ` は平行射の Div-equivalence を保つ。

★★これが `Corollary 4.11, (ii)` の要である ——
原文 (物理 p.93) が `Theorem 4.2, (ii)` に送る段。

★手筋: `divMap` は**全射**なので、`Φ₂.map (Base Ψf)` の値は
`divMap` の像の上で決まる。そこで `divMapNat_all` を両側に当てればよい。 -/
theorem divEquivalent_map (Ψ : C ≌ C₂) (G : Frobenioid P) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hdivS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
      ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (hOTri' : ∀ (Z : C) (δ : End Z),
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z) → δ ∈ OTri P Z)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {W A : C} (f g : W ⟶ A) (h : DivEquivalent P f g) :
    DivEquivalent P₂ (Ψ.functor.map f) (Ψ.functor.map g) := by
  show Φ₂.map (P₂.Base (Ψ.functor.map f)) = Φ₂.map (P₂.Base (Ψ.functor.map g))
  refine AddMonoidHom.ext (fun y => ?_)
  obtain ⟨x, hx⟩ := divMap_surjective Ψ G hiso hdivS hdivS₂ hOTri hOTri' A y
  rw [← hx]
  have hf := divMapNat_all (Ψ := Ψ) (G := G) (hiso := hiso) (hdivS := hdivS)
    (hOTri := hOTri) F hperfM hdeg f x
  have hg := divMapNat_all (Ψ := Ψ) (G := G) (hiso := hiso) (hdivS := hdivS)
    (hOTri := hOTri) F hperfM hdeg g x
  rw [← hf, ← hg]
  have hxy : Φ.map (P.Base f) x = Φ.map (P.Base g) x :=
    congrArg (fun t : Φ.val _ →+ Φ.val _ => t x) h
  rw [hxy]

def divEquivalent_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Theorem 4.9 — Ψ は Div-equivalence を保つ",
    sectionId := "frdi-thm-4-9" }

/-! ## ★5. `Theorem 4.9` 本体 —— 関手の同型 -/

variable (P) in
/-- ★`Φ` を `𝒞` の上の関手と見たもの。

原文 (FrdI p.88):
> [where we regard, for i = 1, 2, the functor Φi : Di → Mon as a functor on Ci, by -/
def phiOnC : Cᵒᵖ ⥤ AddCommMonCat.{w} := P.proj.op ⋙ Φ.functor

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Theorem 4.9** —— `Ψ_Φ : Φ₁ ≅ Φ₂`(`Ψ` の上に乗る関手の同型)。 -/
noncomputable def psiPhi (Ψ : C ≌ C₂) (G : Frobenioid P) (G₂ : Frobenioid P₂)
    (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hdivS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
      ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (hOTri' : ∀ (Z : C) (δ : End Z),
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z) → δ ∈ OTri P Z)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g) :
    phiOnC P ≅ Ψ.functor.op ⋙ phiOnC P₂ :=
  NatIso.ofComponents
    (fun X =>
      { hom := AddCommMonCat.ofHom (divMap Ψ G hiso hdivS hOTri X.unop)
        inv := AddCommMonCat.ofHom
          (divEquiv Ψ G G₂ hiso hiso₂ hdivS hdivS₂ hOTri hOTri' X.unop).symm.toAddMonoidHom
        hom_inv_id := AddCommMonCat.hom_ext (AddMonoidHom.ext (fun x =>
          (divEquiv Ψ G G₂ hiso hiso₂ hdivS hdivS₂ hOTri hOTri' X.unop).symm_apply_apply x))
        inv_hom_id := AddCommMonCat.hom_ext (AddMonoidHom.ext (fun y =>
          (divEquiv Ψ G G₂ hiso hiso₂ hdivS hdivS₂ hOTri hOTri' X.unop).apply_symm_apply y)) })
    (fun {_X _Y} f => AddCommMonCat.hom_ext (AddMonoidHom.ext (fun x =>
      divMapNat_all (Ψ := Ψ) (G := G) (hiso := hiso) (hdivS := hdivS) (hOTri := hOTri)
        F hperfM hdeg f.unop x)))

/-- ★★★★★★**[FrdI] Theorem 4.9** —— 条なしの locator。

| 主張 | 実装 |
|---|---|
| `Ψ_Φ : Φ₁ ≅ Φ₂`(関手の同型) | `psiPhi` |
| `Ψ` の上に乗る(自然性) | `divMapNat_all` |
| `Ψ_Prime`(`Theorem 4.2, (ii)`)との両立 | ★`psiPrime_eq_primeEquiv`(`Thm42Order.lean`) |

★★逸脱 (B)(`hdivS`: `Div : 𝒪^▷(A) ↠ Φ(A)`、`prop_4_4.src` で開示済)を使う。
そのおかげで原文が twin-primary steps で示す
「right-hand ＝ left-hand」が**自明**になっている。 -/
def psiPhi.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Theorem 4.9",
    sectionId := "frdi-thm-4-9" }

/-! ## ★6. `Ψ_Φ` は `Div` を `Div` へ送る —— `Corollary 4.11, (iv)` の核 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★**co-angular pre-step の場合** ——
`hdivS` で `𝒪^▷` の元に取り替え、`exists_iso_of_div_eq`(コスライスの圏同値)で
**同型でずれるだけ**にする。同型は零因子 0・次数 1 なので `Div` は変わらない。 -/
theorem divMap_div_preStep (Ψ : C ≌ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    {A B : C} (ϵ : A ⟶ B) (hs : IsPreStep P ϵ) :
    divMap Ψ G hiso hdivS hOTri A (P.Div ϵ) = P₂.Div (Ψ.functor.map ϵ) := by
  obtain ⟨u, hu⟩ := hdivS A (P.Div ϵ)
  obtain ⟨θ, hθiso, hθ⟩ := exists_iso_of_div_eq G (((u : End A)) : A ⟶ A) ϵ
    (prop_1_4_i P _ (fun Z _ => hiso Z)) (isPreStep_of_otri _ u.2)
    (prop_1_4_i P _ (fun Z _ => hiso Z)) hs hu
  haveI := hθiso
  haveI : IsIso (Ψ.functor.map θ) := Ψ.functor.map_isIso _
  have hdiv2 : P₂.Div (Ψ.functor.map ϵ)
      = P₂.Div (Ψ.functor.map (((u : End A)) : A ⟶ A)) := by
    rw [← hθ, Ψ.functor.map_comp, P₂.Div_comp,
      show P₂.Div (Ψ.functor.map θ) = 0 from isIsometric_of_isIso P₂ _,
      degFr_of_isIso P₂ (Ψ.functor.map θ), map_zero, zero_add,
      show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  rw [hdiv2, ← hu]
  exact divMap_div_otri Ψ G hiso hdivS hOTri u.2

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`Ψ_Φ` は `Div` を `Div` へ送る**(任意の射)。

★`arbFactor` の 3 分解 `φ = γ ≫ pre ≫ plb` で
**`Div γ = 0`(Frobenius 型は isometric)・`Div plb = 0`(pull-back は LB-invertible)**
なので `Div φ = Φ.map (Base γ) (Div pre)` に潰れる。
★あとは `divMapNat_all`(自然性)と pre-step の場合で降りる。 -/
theorem divMap_div_all (Ψ : C ≌ C₂) (G : Frobenioid P) (F : FrobenioidCore P)
    (F₂ : FrobenioidCore P₂)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    (hFT : ∀ {X Y : C} (g : X ⟶ Y), IsFrobeniusType P g →
      IsFrobeniusType P₂ (Ψ.functor.map g))
    (hPS : ∀ {X Y : C} (g : X ⟶ Y), IsPreStep P g → IsPreStep P₂ (Ψ.functor.map g))
    (hPB : ∀ {X Y : C} (g : X ⟶ Y), IsPullBack P g → IsPullBack P₂ (Ψ.functor.map g))
    {A B : C} (φ : A ⟶ B) :
    divMap Ψ G hiso hdivS hOTri A (P.Div φ) = P₂.Div (Ψ.functor.map φ) := by
  obtain ⟨Y₀, Z₀, γ, pre, plb, hfac, hγ, hpre, hplb, -⟩ :
      ∃ (Y₀ Z₀ : C) (γ : A ⟶ Y₀) (pre : Y₀ ⟶ Z₀) (plb : Z₀ ⟶ B),
        φ = γ ≫ pre ≫ plb ∧ IsFrobeniusType P γ ∧ IsPreStep P pre ∧ IsPullBack P plb
          ∧ True := by
    obtain ⟨Y₀, Z₀, γ, pre, plb, hfac, hγ, hpre, hplb⟩ := F.arbFactor φ
    exact ⟨Y₀, Z₀, γ, pre, plb, hfac, hγ, hpre, hplb, trivial⟩
  -- ★`Div (pre ≫ plb) = Div pre`
  have hdplb : P.Div plb = 0 := (F.pullBackLB plb hplb).1.2
  have hdlin : P.degFr plb = 1 := (F.pullBackLB plb hplb).2
  have hd1 : P.Div (pre ≫ plb) = P.Div pre := by
    rw [P.Div_comp, hdplb, hdlin, map_zero, zero_add,
      show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  -- ★`Div φ = Φ.map (Base γ) (Div pre)`
  have hdγ : P.Div γ = 0 := hγ.1.2
  have hd2 : P.Div φ = Φ.map (P.Base γ) (P.Div pre) := by
    rw [hfac, P.Div_comp, hd1, hdγ, smul_zero, add_zero]
  -- ★`𝒞₂` 側でも同じ計算
  have hdplb₂ : P₂.Div (Ψ.functor.map plb) = 0 :=
    (F₂.pullBackLB _ (hPB plb hplb)).1.2
  have hdlin₂ : P₂.degFr (Ψ.functor.map plb) = 1 := by rw [hdeg plb, hdlin]
  have hd1₂ : P₂.Div (Ψ.functor.map pre ≫ Ψ.functor.map plb)
      = P₂.Div (Ψ.functor.map pre) := by
    rw [P₂.Div_comp, hdplb₂, hdlin₂, map_zero, zero_add,
      show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  have hdγ₂ : P₂.Div (Ψ.functor.map γ) = 0 := (hFT γ hγ).1.2
  have hd2₂ : P₂.Div (Ψ.functor.map φ)
      = Φ₂.map (P₂.Base (Ψ.functor.map γ)) (P₂.Div (Ψ.functor.map pre)) := by
    rw [hfac, Ψ.functor.map_comp, Ψ.functor.map_comp, P₂.Div_comp, hd1₂, hdγ₂,
      smul_zero, add_zero]
  rw [hd2, hd2₂]
  rw [divMapNat_all (Ψ := Ψ) (G := G) (hiso := hiso) (hdivS := hdivS) (hOTri := hOTri)
    F hperfM hdeg γ (P.Div pre)]
  exact congrArg (fun t => Φ₂.map (P₂.Base (Ψ.functor.map γ)) t)
    (divMap_div_preStep Ψ G hiso hdivS hOTri pre hpre)

def divMap_div_all.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 94,
    item := "Corollary 4.11, (iv) — Ψ_Φ は Div を Div へ送る",
    sectionId := "frdi-cor-4-11" }

/-- ★`psiPhi` の成分は `divMap` そのもの(定義から)。 -/
theorem psiPhi_app_apply (Ψ : C ≌ C₂) (G : Frobenioid P) (G₂ : Frobenioid P₂)
    (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hdivS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
      ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (hOTri' : ∀ (Z : C) (δ : End Z),
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z) → δ ∈ OTri P Z)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    (A : C) (x : Φ.val (P.toElem.obj A).base) :
    ((psiPhi Ψ G G₂ F hiso hiso₂ hdivS hdivS₂ hOTri hOTri' hperfM hdeg).hom.app
        (Opposite.op A)).hom x
      = divMap Ψ G hiso hdivS hOTri A x := rfl

end PsiPhi

end ABC3.Found.FrdI
