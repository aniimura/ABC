/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfSlice
import ABC3.Found.FrdI.Prop55Std

/-!
# [FrdI] Proposition 5.5, (iv) —— `𝒞^pf` の有理関数の単系は `ℚ·Φ^birat`

原文 (FrdI p.105):
> Finally, assertion (iv) is immediate from the definitions [cf. also assertions (i), (ii);

★★`Prop55Std.lean` の `pfRoot_ratFnData` は `𝒞^pf` の有理関数の単系を
`(Φ^pf)^birat`(`biratRatFnData` 経由、`phiBiratOn Gpf` から作る)として立てているが、
原文は `Proposition 5.3` に合わせて `ℚ·Φ^birat` と書く。

★★★`Prop55PfSlice.lean` の `phiBiratOn_pf_eq_qPhiBiratOn_all` で
**両者が `𝒟` 上の部分関手として一致する**ことが示されたので、本ファイルで
`pfRoot_ratFnData` の言葉に翻訳しておく。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] [IsConnected D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P} {G : Frobenioid P}

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (iv) の単系の同定** ——
`pfRoot_ratFnData` が立てる有理関数の単系は、各 `d ∈ Ob(𝒟)` で
**`ℚ·Φ^birat(d)`**(`qPhiBiratOn`)である。

★`pfRoot_ratFnData` の `bmon` は `phiBiratOn Gpf` から作られており、
`phiBiratOn_pf_eq_qPhiBiratOn_all` がそれを `ℚ·Φ^birat` に同定する。 -/
theorem pfRoot_ratFnData_bmon_val (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hftr : ∀ X : C, IsFrobeniusTrivial P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (ζC : ∀ X : C, ℕ+ →* End X)
    (hdegC : ∀ (X : C) (m : ℕ+), P.degFr ((ζC X m : End X) : X ⟶ X) = m)
    (hpropC : ∀ (X : C) (m : ℕ+),
      IsBaseIdentity P (ζC X m) ∧ IsFrobeniusType P ((ζC X m : End X) : X ⟶ X))
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D)
    (F' : FrobenioidCore (biratPre P G))
    (hisoB : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X)
    (hfnBirat : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hAB : ∀ A : C, IsFrobeniusTrivial (biratPre P G) (biratUp P G A))
    (hfnB' : ∀ A : C, IsFrobeniusNormalized (biratPre P G)
      (rtObj (biratPre P G) F' (biratUp P G A) 1))
    (ζ : ∀ A : C, ℕ+ →* End (biratUp P G A))
    (hdeg : ∀ (A : C) (n : ℕ+), (biratPre P G).degFr
      ((ζ A n : End (biratUp P G A)) : biratUp P G A ⟶ biratUp P G A) = n)
    (hprop : ∀ (A : C) (n : ℕ+), IsBaseIdentity (biratPre P G) (ζ A n)
      ∧ IsFrobeniusType (biratPre P G) ((ζ A n : End (biratUp P G A)) : _ ⟶ _))
    (hfnPfBirat : ∀ X : BiratCat (pfRootPre P F) Gpf,
      IsFrobeniusNormalized (biratPre (pfRootPre P F) Gpf) X)
    (d : D) :
    (pfRoot_ratFnData (F := F) hfi hiso hftr hfnC ζC hdegC hpropC hut Gpf hfsmD).bmon.val d
      = ↥(qPhiBiratOn P G d) := by
  show ↥(phiBiratOn Gpf d) = ↥(qPhiBiratOn P G d)
  rw [phiBiratOn_pf_eq_qPhiBiratOn_all hfi hiso Gpf F' hisoB hfnBirat hAB hfnB'
    ζ hdeg hprop hfnPfBirat d]
  rfl

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.5, (iv)` の単系の同定を
`pfRoot_ratFnData` の言葉で述べたもの。 -/
def pfRoot_ratFnData_bmon_val.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iv) — 𝒞^pf の有理関数の単系は ℚ·Φ^birat",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
