import ABC3.Found.Arakelov.PicInvTrivial
import ABC3.Interface.Arakelov.LineBundle

/-!
# Arakelov (B1) 第 73 ブロック —— **`Found` の可逆層と `Interface` の可逆層は同じもの**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`sheafOf` 族を埋める

`Interface/Arakelov/LineBundle.lean` の `PicardData` は 14 フィールド。
★第 62 ブロックまでで **6 つ**(`Pic`・`group`・`pullback`・`_mul`・`_id`・`_comp`)が済んだ。

★★本ブロックは **`sheafOf` 族の 6 つ**を埋める:

| フィールド | 機構 |
|---|---|
| `sheafOf` | ★`Quotient.out` で代表を選ぶ |
| `sheafOf_invertible` | ★`InvSheaf` → `IsInvertibleSheaf`(そのまま) |
| `sheafOf_one` | ★`Quotient.mk_out` |
| `sheafOf_mul` | ★`Quotient.mk_out` |
| `sheafOf_pullback` | ★`Quotient.mk_out` + 第 61 の `carrier` が `rfl` |
| `sheafOf_injective` | ★`Quotient.out_eq` |
| `sheafOf_surjective` | ★★第 72 ブロック(逆も局所自明)が要る |

★★★残るのは **`equivPicRing` ただ 1 つ**である。

## ★★`Interface` は `Found` を import できない——逆向きは許される

`Found/GenEll/HeightInterface.lean` と同じ形である。
★`Interface` 側の条件は `restrictPresheafFunctor` を mathlib の
`pushforward₀OfCommRingCat` に展開して書いてあるので、`abbrev` を通じて
**定義的に一致する**。
-/

universe u

namespace ABC3.Found.Arakelov

open ABC3.Interface.Arakelov AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

variable (X : Scheme.{0})

/-! ## ★★`Found` → `Interface` -/

/-- ★**`Found` の可逆層は `Interface` の意味でも可逆層である**。 -/
theorem isInvertibleSheaf_carrier (L : InvSheaf X) : IsInvertibleSheaf L.carrier :=
  ⟨⟨L.inv, L.isInv⟩, L.trivial⟩

/-! ## ★★`Interface` → `Found`(第 72 ブロックが要る) -/

/-- ★★**`Interface` の可逆層から `Found` の可逆層を作る**。

★★★逆の側の局所自明性は第 72 ブロック(`isLocallyTrivial_of_tensor_unit`)から出る。 -/
noncomputable def invSheafOfIsInvertible (F : X.Modules) (h : IsInvertibleSheaf F) :
    InvSheaf X where
  carrier := F
  inv := h.1.choose
  isInv := h.1.choose_spec
  trivial := h.2
  invTrivial :=
    isLocallyTrivial_of_tensor_unit X F h.1.choose h.1.choose_spec.some h.2

@[simp] theorem invSheafOfIsInvertible_carrier (F : X.Modules) (h : IsInvertibleSheaf F) :
    (invSheafOfIsInvertible X F h).carrier = F := rfl

/-! ## ★★★`sheafOf` -/

/-- ★★★**`Pic X` の元の下にある可逆層**——同型類から代表を選ぶ。 -/
noncomputable def picSheafOf (L : PicSheaf X) : X.Modules := L.out.carrier

/-- ★代表は元を再現する。 -/
theorem picSheafOf_spec (L : InvSheaf X) :
    Nonempty (picSheafOf X (PicSheaf.mk L) ≅ L.carrier) :=
  Quotient.mk_out (s := InvSheaf.setoid X) L

/-- ★**下にある層は可逆層である**。 -/
theorem picSheafOf_invertible (L : PicSheaf X) : IsInvertibleSheaf (picSheafOf X L) :=
  isInvertibleSheaf_carrier X L.out

/-- ★**単位元の下にあるのは構造層**。 -/
theorem picSheafOf_one :
    Nonempty (picSheafOf X 1 ≅ SheafOfModules.unit X.ringCatSheaf) :=
  picSheafOf_spec X (InvSheaf.one X)

/-- ★★**積はテンソル積に移る**。 -/
theorem picSheafOf_mul (L M : PicSheaf X) :
    Nonempty (picSheafOf X (L * M) ≅ modTensor X (picSheafOf X L) (picSheafOf X M)) := by
  have h : L * M = PicSheaf.mk (InvSheaf.mul L.out M.out) := by
    conv_lhs => rw [← Quotient.out_eq L, ← Quotient.out_eq M]
    rfl
  rw [h]
  exact picSheafOf_spec X (InvSheaf.mul L.out M.out)

/-- ★★★★★**引き戻しは層の引き戻しである**。 -/
theorem picSheafOf_pullback {X Y : Scheme.{0}} (f : X ⟶ Y) (L : PicSheaf Y) :
    Nonempty (picSheafOf X (picPullback f L)
      ≅ (Scheme.Modules.pullback f).obj (picSheafOf Y L)) := by
  have h : picPullback f L = PicSheaf.mk (InvSheaf.pullback f L.out) := by
    conv_lhs => rw [← Quotient.out_eq L]
    rfl
  rw [h]
  exact picSheafOf_spec X (InvSheaf.pullback f L.out)

/-- ★★★**同型な層は同じ元を与える**。 -/
theorem picSheafOf_injective (L M : PicSheaf X)
    (h : Nonempty (picSheafOf X L ≅ picSheafOf X M)) : L = M := by
  conv_lhs => rw [← Quotient.out_eq L]
  conv_rhs => rw [← Quotient.out_eq M]
  exact PicSheaf.mk_eq_mk h.some

/-- ★★★**可逆層はすべて現れる**。 -/
theorem picSheafOf_surjective (F : X.Modules) (h : IsInvertibleSheaf F) :
    ∃ L : PicSheaf X, Nonempty (picSheafOf X L ≅ F) :=
  ⟨PicSheaf.mk (invSheafOfIsInvertible X F h), picSheafOf_spec X _⟩

/-! ## ★出典の紐付け(`.src`) -/

def picSheafOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Pic の元の下にある可逆層)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
