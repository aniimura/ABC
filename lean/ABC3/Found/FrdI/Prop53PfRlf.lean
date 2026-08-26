/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfRatFn

/-!
# [FrdI] `Proposition 5.3` の図式の**右の縦の矢印** `𝒞^pf ⟶ 𝒞^rlf`

原文 (FrdI p.103):
> if C is of Frobenius-isotropic type, then there is a natural 1-commutative

★★`Proposition 5.3` の 1-可換図式は

```
     𝒞  ⟶ 𝒞^un-tr ⟶ (𝒞^un-tr)^pf ⟶ 𝒞^rlf
     │        │            │            ║
     ↓        ↓            ↓            ║
    …       …           𝒞^pf   ⟶⟶⟶⟶  𝒞^rlf
```

という形で、**右の縦の矢印は `𝒞^pf ⟶ 𝒞^rlf`** である。

## ★これまで詰まっていた理由と、詰まりが取れた理由

`𝒞^pf` の model data の零因子は `Φ^pf = Pf Φ`(**商模型**)であり、
`𝒞^rlf` の側は `Φ ⊗_ℕ ℝ≥0`(**テンソル模型**)である。零因子の側は
`pfEquivPfT`(`Pf M ≃+ ℚ≥0 ⊗_ℕ M`)と係数拡大 `ℚ≥0 → ℝ≥0` で繋がる。

★★問題は**有理関数の単系の側**であった。`𝒞^pf` の側は
`(Φ^pf)^birat`(`phiBiratOn Gpf`)であって、これが
**`ℚ·Φ^birat`(`qPhiBiratOn`)である**ことを言わないと
`ℝ·Φ^birat`(`sPhiBiratOn ℝ≥0`)へ写せない。
★これは `Prop55PfSlice.lean` の `phiBiratOn_pf_eq_qPhiBiratOn_all` で閉じた。

★★★あとは **`ℚ·Φ^birat ⟶ ℝ·Φ^birat`**(`qPhiBiratOn_map_rlf`、`Prop53QPhi.lean`)
—— 分母 `k` を払って `Φ^birat` の像に落とし、`ℝ≥0` 側ではスカラー倍で割り戻す ——
を使って `ModelDataHom` を組むだけである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section PfRlf

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. 零因子の側 —— `Φ^pf ⟶ Φ ⊗_ℕ ℝ≥0` の自然性 -/

/-- ★★`pfToRlfHom` は単系の射と可換(`Pf.map` ↔ `scMap`)。 -/
theorem pfToRlfHom_map {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (g : M →+ N)
    (x : Pf M) : pfToRlfHom N (Pf.map g x) = scMap (S := NNReal) g (pfToRlfHom M x) := by
  show scBase (NNRat.castHom NNReal) (pfToPfT (Pf.map g x)) = _
  rw [pfToPfT_natural, scBase_scMap]
  rfl

/-- ★その `Gp` 版。 -/
theorem gpMap_pfToRlfHom_map {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (g : M →+ N)
    (x : Gp (Pf M)) :
    gpMap _ (pfToRlfHom N) (gpMap _ (Pf.map g) x)
      = gpMap _ (scMap (S := NNReal) g) (gpMap _ (pfToRlfHom M) x) := by
  have hc : (pfToRlfHom N).comp (Pf.map g)
      = (scMap (S := NNReal) g).comp (pfToRlfHom M) := by
    ext y; exact pfToRlfHom_map g y
  rw [← AddMonoidHom.comp_apply, ← gpMap_comp, ← AddMonoidHom.comp_apply, ← gpMap_comp, hc]

/-! ## ★2. `ModelDataHom` -/

variable [IsConnected D] {F : FrobenioidCore P} {G : Frobenioid P}

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**`Proposition 5.3` の図式の右の縦の矢印**(model data の層)——
`𝒞^pf` の model data `(Φ^pf, (Φ^pf)^birat)` から
`𝒞^rlf` の model data `(Φ ⊗ ℝ≥0, ℝ·Φ^birat)` への射。

★零因子の側は `pfToRlfHom`(`Pf Φ ≃ ℚ≥0 ⊗ Φ ⟶ ℝ≥0 ⊗ Φ`)、
有理関数の側は `heq`(`(Φ^pf)^birat = ℚ·Φ^birat`)で書き換えたうえで
`qPhiBiratOn_map_rlf` である。

★`heq` は `phiBiratOn_pf_eq_qPhiBiratOn_all`(`Prop55PfSlice.lean`)が与える。 -/
noncomputable def pfRootToRlfHom (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hftr : ∀ X : C, IsFrobeniusTrivial P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (ζC : ∀ X : C, ℕ+ →* End X)
    (hdegC : ∀ (X : C) (m : ℕ+), P.degFr ((ζC X m : End X) : X ⟶ X) = m)
    (hpropC : ∀ (X : C) (m : ℕ+),
      IsBaseIdentity P (ζC X m) ∧ IsFrobeniusType P ((ζC X m : End X) : X ⟶ X))
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D)
    (hfnB : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A)))
    (heq : ∀ d : D, phiBiratOn Gpf d = qPhiBiratOn P G d) :
    ModelDataHom
        (pfRoot_ratFnData (F := F) hfi hiso hftr hfnC ζC hdegC hpropC hut Gpf hfsmD).model
        (scModel NNReal G hiso hfnB hcharInj' hint' hfsmD) where
  phiHom d := pfToRlfHom (Φ.val d)
  phiNat f x := pfToRlfHom_map (Φ.map f) x
  bmonHom d :=
    AddMonoidHom.codRestrict
      ((gpMap _ (pfToRlfHom (Φ.val d))).comp (phiBiratOn Gpf d).subtype) _
      (fun x => qPhiBiratOn_map_rlf G (by rw [← heq d]; exact x.2))
  bmonNat := by
    intro A B f x
    exact Subtype.ext (gpMap_pfToRlfHom_map (Φ.map f) _)
  divCompat _ _ := rfl

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**`Proposition 5.3` の図式の右の縦の矢印** `𝒞^pf ⟶ 𝒞^rlf`。

★`𝒞^pf` を model Frobenioid で表したもの(`pfRoot_ratFnData … |>.model` の対象)から
`𝒞^rlf` への関手である。 -/
noncomputable def pfRootToRlfFunctor (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hftr : ∀ X : C, IsFrobeniusTrivial P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (ζC : ∀ X : C, ℕ+ →* End X)
    (hdegC : ∀ (X : C) (m : ℕ+), P.degFr ((ζC X m : End X) : X ⟶ X) = m)
    (hpropC : ∀ (X : C) (m : ℕ+),
      IsBaseIdentity P (ζC X m) ∧ IsFrobeniusType P ((ζC X m : End X) : X ⟶ X))
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D)
    (hfnB : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A)))
    (heq : ∀ d : D, phiBiratOn Gpf d = qPhiBiratOn P G d) :
    ModelData.Obj
        (pfRoot_ratFnData (F := F) hfi hiso hftr hfnC ζC hdegC hpropC hut Gpf hfsmD).model
      ⥤ Crlf G hiso hfnB hcharInj' hint' hfsmD :=
  (pfRootToRlfHom hfi hiso hftr hfnC ζC hdegC hpropC hut Gpf hfsmD hfnB
    hcharInj' hint' heq).functor

end PfRlf

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.3` の 1-可換図式の**右の縦の矢印**
`𝒞^pf ⟶ 𝒞^rlf`。 -/
def pfRootToRlfFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 1-可換図式の右の縦の矢印 𝒞^pf ⥤ 𝒞^rlf",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
