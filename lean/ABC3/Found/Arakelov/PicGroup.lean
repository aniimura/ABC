import ABC3.Found.Arakelov.PicLocalTrivial

/-!
# Arakelov (B1) 第 16 ブロック —— **可逆前層の群構造**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★層化を通さない

これまで `tensorModules`(前層でテンソル → 層化)で組んできたが、
★★★**層化を通すたびに「局所自明性が保たれるか」を示さねばならない**。

★★★★**前層の段で組めば、その問いは消える**:

| 何 | 前層の段 | 層化を通す段 |
|---|---|---|
| モノイダル構造 | ★**mathlib が持っている**(第 1 ブロック) | 無い(我々が 13 ブロックかけて作った) |
| 結合律・単位律・可換律 | ★**無料** | ★局所論法が要る |
| テンソル閉性 | ★第 15 ブロックで取得 | ★層化の局所自明性が要る |

★★**局所自明な前層は自動的に層である**(層条件は局所的だから)ので、
両者は数学的に一致する。★その同一視は最後に 1 回だけ要る。

## ★本ブロックで取れるもの

`InvertiblePresheaf X`(可逆前層)と、その **群構造の全部品**:
`mul` / `one` / `symm`、および `isInv` の作り直し。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

/-- ★★★**可逆前層(直線束の前層版)**。 -/
structure InvertiblePresheaf (X : Scheme.{u}) where
  /-- 台。 -/
  carrier : X.PresheafOfModules
  /-- テンソル積についての逆。 -/
  inv : X.PresheafOfModules
  /-- ★逆であること。 -/
  isInv : Nonempty (carrier ⊗ inv ≅ 𝟙_ X.PresheafOfModules)
  /-- ★局所自明。 -/
  trivial : IsLocallyTrivial X carrier
  /-- ★逆の側も局所自明。 -/
  invTrivial : IsLocallyTrivial X inv

namespace InvertiblePresheaf

variable {X : Scheme.{u}}

/-- ★**構造層は可逆前層である**(単位元)。 -/
noncomputable def one (X : Scheme.{u}) : InvertiblePresheaf X where
  carrier := 𝟙_ X.PresheafOfModules
  inv := 𝟙_ X.PresheafOfModules
  isInv := ⟨λ_ _⟩
  trivial := fun U => ⟨⊤, (Opens.grothendieckTopology X).top_mem U,
    fun V _ _ => ⟨(restrictPresheafUnit (X := X) (U := V)).symm⟩⟩
  invTrivial := fun U => ⟨⊤, (Opens.grothendieckTopology X).top_mem U,
    fun V _ _ => ⟨(restrictPresheafUnit (X := X) (U := V)).symm⟩⟩

/-- ★★**逆を取る操作**(逆元)。 -/
noncomputable def symm (L : InvertiblePresheaf X) : InvertiblePresheaf X where
  carrier := L.inv
  inv := L.carrier
  isInv := L.isInv.map fun e => (β_ L.inv L.carrier) ≪≫ e
  trivial := L.invTrivial
  invTrivial := L.trivial

/-- ★★★**テンソル積**(乗法)。

★★★逆の存在は `tensorTensorTensorComm`(braided モノイダル圏の並べ替え)で出る:

    (A ⊗ B) ⊗ (A' ⊗ B') ≅ (A ⊗ A') ⊗ (B ⊗ B') ≅ 𝟙_ ⊗ 𝟙_ ≅ 𝟙_

★局所自明性は第 15 ブロック `IsLocallyTrivial.tensor`。 -/
noncomputable def mul (L M : InvertiblePresheaf X) : InvertiblePresheaf X where
  carrier := L.carrier ⊗ M.carrier
  inv := L.inv ⊗ M.inv
  isInv := by
    obtain ⟨eL⟩ := L.isInv
    obtain ⟨eM⟩ := M.isInv
    -- ★(A ⊗ B) ⊗ (A' ⊗ B') ≅ (A ⊗ A') ⊗ (B ⊗ B') を結合子と braiding で組む
    exact ⟨α_ L.carrier M.carrier (L.inv ⊗ M.inv)
      ≪≫ (Iso.refl L.carrier ⊗ᵢ
            ((α_ M.carrier L.inv M.inv).symm
              ≪≫ (β_ M.carrier L.inv ⊗ᵢ Iso.refl M.inv)
              ≪≫ α_ L.inv M.carrier M.inv))
      ≪≫ (α_ L.carrier L.inv (M.carrier ⊗ M.inv)).symm
      ≪≫ (eL ⊗ᵢ eM) ≪≫ λ_ _⟩
  trivial := L.trivial.tensor M.trivial
  invTrivial := L.invTrivial.tensor M.invTrivial

@[simp] theorem mul_carrier (L M : InvertiblePresheaf X) :
    (L.mul M).carrier = L.carrier ⊗ M.carrier := rfl

@[simp] theorem one_carrier (X : Scheme.{u}) :
    (one X).carrier = 𝟙_ X.PresheafOfModules := rfl

@[simp] theorem symm_carrier (L : InvertiblePresheaf X) : L.symm.carrier = L.inv := rfl

/-! ## ★★★群法則(前層のモノイダル構造から**無料**で出る) -/

/-- ★★**結合律**——前層の結合子そのもの。 -/
noncomputable def mulAssoc (L M N : InvertiblePresheaf X) :
    ((L.mul M).mul N).carrier ≅ (L.mul (M.mul N)).carrier :=
  α_ L.carrier M.carrier N.carrier

/-- ★★**可換律**——前層の braiding そのもの。 -/
noncomputable def mulComm (L M : InvertiblePresheaf X) :
    (L.mul M).carrier ≅ (M.mul L).carrier :=
  β_ L.carrier M.carrier

/-- ★★**単位律**——前層の左単位子そのもの。 -/
noncomputable def oneMul (L : InvertiblePresheaf X) :
    ((one X).mul L).carrier ≅ L.carrier :=
  λ_ L.carrier

/-- ★★**逆元の法則**——`isInv` そのもの。 -/
theorem mulSymm (L : InvertiblePresheaf X) :
    Nonempty ((L.mul L.symm).carrier ≅ (one X).carrier) :=
  L.isInv

/-! ## ★出典の紐付け(`.src`) -/

def picGroupStructure.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆前層とその群構造)",
    sectionId := "genell-def-1-1-i" }

def picGroupMul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆前層のテンソル積)",
    sectionId := "genell-def-1-1-i" }

end InvertiblePresheaf

end ABC3.Found.Arakelov
