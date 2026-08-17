import ABC3.Found.Arakelov.PicGroup

/-!
# Arakelov (B1) 第 17 ブロック —— **`Pic X` と `CommGroup`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★可逆前層の同型類

第 16 ブロックで可逆前層の群構造の**部品**が揃った。
★本ブロックはそれを**同型類の商**に降ろし、`CommGroup` を載せる。

## ★★逆の一意性が要る

`inv` を商に降ろすには「**台が同型なら逆も同型**」が要る。
★★★これはモノイダル圏の標準の議論である:

    A ≅ 𝟙 ⊗ A ≅ (M ⊗ B) ⊗ A ≅ M ⊗ (B ⊗ A) ≅ M ⊗ (A ⊗ B)
      ≅ (M ⊗ A) ⊗ B ≅ (L ⊗ A) ⊗ B ≅ 𝟙 ⊗ B ≅ B
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

namespace InvertiblePresheaf

variable {X : Scheme.{u}}

/-- ★★**逆の一意性** —— 台が同型なら逆も同型である。

★★★モノイダル圏の標準の議論(8 段)。 -/
theorem invIsoOfCarrierIso {L M : InvertiblePresheaf X}
    (e : L.carrier ≅ M.carrier) : Nonempty (L.inv ≅ M.inv) := by
  obtain ⟨eL⟩ := L.isInv
  obtain ⟨eM⟩ := M.isInv
  exact ⟨(λ_ L.inv).symm
    ≪≫ (eM.symm ⊗ᵢ Iso.refl L.inv)
    ≪≫ α_ M.carrier M.inv L.inv
    ≪≫ (Iso.refl M.carrier ⊗ᵢ β_ M.inv L.inv)
    ≪≫ (α_ M.carrier L.inv M.inv).symm
    ≪≫ ((e.symm ⊗ᵢ Iso.refl L.inv) ⊗ᵢ Iso.refl M.inv)
    ≪≫ (eL ⊗ᵢ Iso.refl M.inv)
    ≪≫ λ_ M.inv⟩

/-- ★可逆前層の同型による同値関係(台が同型であること)。 -/
def setoid (X : Scheme.{u}) : Setoid (InvertiblePresheaf X) where
  r L M := Nonempty (L.carrier ≅ M.carrier)
  iseqv :=
    { refl := fun _ => ⟨Iso.refl _⟩
      symm := fun ⟨e⟩ => ⟨e.symm⟩
      trans := fun ⟨e⟩ ⟨f⟩ => ⟨e ≪≫ f⟩ }

end InvertiblePresheaf

/-- ★★★★★**`Pic X`** —— 可逆前層の同型類。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
def PicPre (X : Scheme.{u}) : Type _ :=
  Quotient (InvertiblePresheaf.setoid X)

namespace PicPre

variable {X : Scheme.{u}}

/-- ★同型類を取る。 -/
def mk (L : InvertiblePresheaf X) : PicPre X := Quotient.mk _ L

theorem mk_eq_mk {L M : InvertiblePresheaf X} (e : L.carrier ≅ M.carrier) :
    mk L = mk M := Quotient.sound ⟨e⟩

/-- ★★★★★**`Pic X` は可換群である**。

★★★乗法・単位元・逆元・結合律・可換律・逆元の法則——
すべて第 16 ブロックの部品を商へ降ろすだけである。 -/
noncomputable instance : CommGroup (PicPre X) where
  mul := Quotient.map₂ InvertiblePresheaf.mul
    (by
      rintro L L' ⟨e⟩ M M' ⟨f⟩
      exact ⟨e ⊗ᵢ f⟩)
  one := mk (InvertiblePresheaf.one X)
  inv := Quotient.map InvertiblePresheaf.symm
    (by
      rintro L L' ⟨e⟩
      exact InvertiblePresheaf.invIsoOfCarrierIso e)
  mul_assoc := by
    rintro ⟨L⟩ ⟨M⟩ ⟨N⟩
    exact mk_eq_mk (InvertiblePresheaf.mulAssoc L M N)
  one_mul := by
    rintro ⟨L⟩
    exact mk_eq_mk (InvertiblePresheaf.oneMul L)
  mul_one := by
    rintro ⟨L⟩
    exact mk_eq_mk (ρ_ L.carrier)
  inv_mul_cancel := by
    rintro ⟨L⟩
    obtain ⟨e⟩ := L.symm.isInv
    exact mk_eq_mk e
  mul_comm := by
    rintro ⟨L⟩ ⟨M⟩
    exact mk_eq_mk (InvertiblePresheaf.mulComm L M)

end PicPre

/-! ## ★出典の紐付け(`.src`) -/

def PicPre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Pic X = 可逆前層の同型類)",
    sectionId := "genell-def-1-1-i" }

def picCommGroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Pic X が可換群であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
