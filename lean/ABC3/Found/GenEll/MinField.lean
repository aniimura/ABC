import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.ResidueField
import Mathlib.Algebra.Field.Subfield.Basic

/-!
# [GenEll] Definition 1.5, (i) —— minimal field of definition(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> (i) Note if x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, then by considering the scheme-theoretic image of the corresponding morphism Spec(F ) → X, one obtains a well-deﬁned minimal ﬁeld of deﬁnition Fmin ⊆ F of x. In particular, it makes sense to say that F “is a minimal ﬁeld of deﬁnition of x” [i.e., that F = Fmin].

## ★何を取るか

`x ∈ X(F)` を `x_F : Spec(F) → X` として受け、

> **`F_min ≝ κ(ξ) の F の中での像`**(`ξ ≝ x_F` の像点)

を構成する。★原文は「scheme 論的像」を経由するが、
`Spec(F) → X` の scheme 論的像は像点 `ξ` の被約閉包であり、
その関数体が `κ(ξ)` である——**同じものである**。

## ★★極小性まで取る(これが本体)

「`F_min`」と名乗るには**極小であること**を言わねばならない:

> **`x_F` が部分体 `E ⊆ F` を経由するなら `F_min ⊆ E`**

★★これを証明しないまま `F_min` と名付けるのは**過剰主張**である。
(同種の過剰主張を 2026-08-17 に 1 度やって取り消している——
`deg_adivRed_le` に条なしの `.src` を付けて指標を 4→5 に動かした件。)

## ★機構 —— mathlib の `SpecToEquivOfField`

**`AlgebraicGeometry.Scheme.SpecToEquivOfField`**(2026-08-17 実測、mathlib にある):

> `(Spec K ⟶ X) ≃ Σ x : X, (κ(x) ⟶ K)`

★**全単射**なので、`F_min` の構成(`toFun`)と分解(`invFun`)が同時に手に入り、
極小性は**単射性を 1 回使うだけ**で出る。

★★またしても「正面から作る必要が無かった」——scheme 論的像の理論も、
被約閉包も、関数体も、**1 つも要らなかった**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

section MinimalField

variable (F : Type*) [Field F] {X : Scheme}

/-- ★`x_F : Spec(F) → X` に対する「像点 `ξ` と埋め込み `κ(ξ) ↪ F`」の対。

★mathlib の `SpecToEquivOfField` の像そのものである。 -/
noncomputable def minFieldData (xF : Spec (CommRingCat.of F) ⟶ X) :
    Σ x : X, X.residueField x ⟶ CommRingCat.of F :=
  Scheme.SpecToEquivOfField F X xF

/-- ★★**`F_min`** —— `κ(ξ) ↪ F` の像として得られる `F` の部分体。

原文 (GenEll p.8):
> (i) Note if x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, then by considering the scheme-theoretic image of the corresponding morphism Spec(F ) → X, one obtains a well-deﬁned minimal ﬁeld of deﬁnition Fmin ⊆ F of x.
-/
noncomputable def minField (xF : Spec (CommRingCat.of F) ⟶ X) : Subfield F :=
  (minFieldData F xF).2.hom.fieldRange

/-- ★**`x_F` は `Spec κ(ξ)` を経由する** —— `F_min` が実際に定義体であること。

★`SpecToEquivOfField` の `left_inv` そのものである。 -/
theorem specMap_minFieldData_comp (xF : Spec (CommRingCat.of F) ⟶ X) :
    Spec.map (minFieldData F xF).2
        ≫ X.fromSpecResidueField (minFieldData F xF).1 = xF :=
  Scheme.descResidueField_stalkClosedPointTo_fromSpecResidueField F X xF

/-- ★★★**極小性** —— `x_F` が部分体 `E ⊆ F` を経由するなら `F_min ⊆ E`。

★★これが `F_min` を「**minimal** field of definition」と呼べる理由である。

原文 (GenEll p.8):
> In particular, it makes sense to say that F “is a minimal ﬁeld of deﬁnition of x” [i.e., that F = Fmin].

★証明は `SpecToEquivOfField` の**単射性を 1 回**使うだけ:
`x_F = Spec(E↪F) ≫ x_E` を `invFun` の形に直すと
`Spec(ψ ≫ ι) ≫ fromSpecResidueField` になり、
単射性から `F_min` の生成元が `ι` の像に入る。 -/
theorem minField_le_of_factors (E : Subfield F)
    (xE : Spec (CommRingCat.of (E : Type _)) ⟶ X)
    (xF : Spec (CommRingCat.of F) ⟶ X)
    (h : xF = Spec.map (CommRingCat.ofHom E.subtype) ≫ xE) :
    minField F xF ≤ E := by
  set ι : CommRingCat.of (E : Type _) ⟶ CommRingCat.of F :=
    CommRingCat.ofHom E.subtype with hι
  -- `x_E` 側の対
  set dE : Σ x : X, X.residueField x ⟶ CommRingCat.of (E : Type _) :=
    Scheme.SpecToEquivOfField (E : Type _) X xE with hdE
  have hE : Spec.map dE.2 ≫ X.fromSpecResidueField dE.1 = xE :=
    Scheme.descResidueField_stalkClosedPointTo_fromSpecResidueField _ X xE
  -- `x_F` を `⟨dE.1, dE.2 ≫ ι⟩` の像として書き直す
  have hfac : Spec.map (dE.2 ≫ ι) ≫ X.fromSpecResidueField dE.1 = xF := by
    rw [Spec.map_comp, Category.assoc, hE, h]
  have hkey : minFieldData F xF = ⟨dE.1, dE.2 ≫ ι⟩ := by
    apply (Scheme.SpecToEquivOfField F X).symm.injective
    rw [minFieldData, Equiv.symm_apply_apply]
    exact hfac.symm
  -- 像が `E` に入ることを読む
  rw [minField]
  rintro y ⟨a, rfl⟩
  have : (minFieldData F xF).2 = (Scheme.residueFieldCongr (X := X)
      (congrArg Sigma.fst hkey)).hom ≫ dE.2 ≫ ι :=
    (Scheme.SpecToEquivOfField_eq_iff.mp hkey).choose_spec
  rw [this]
  exact (dE.2.hom ((Scheme.residueFieldCongr (congrArg Sigma.fst hkey)).hom.hom a)).2

/-- ★**`F` が `x` の minimal field of definition である**(原文 `F = F_min`)。 -/
def IsMinimalFieldOfDefinition (xF : Spec (CommRingCat.of F) ⟶ X) : Prop :=
  minField F xF = ⊤

/-- ★★**負の対照** —— `x_F` が真の部分体 `E < F` を経由するなら、
`F` は minimal field of definition では**ない**。

★これが出ないなら `minField` は極小性を捉えていない。 -/
theorem not_isMinimalFieldOfDefinition_of_factors (E : Subfield F) (hE : E ≠ ⊤)
    (xE : Spec (CommRingCat.of (E : Type _)) ⟶ X)
    (xF : Spec (CommRingCat.of F) ⟶ X)
    (h : xF = Spec.map (CommRingCat.ofHom E.subtype) ≫ xE) :
    ¬ IsMinimalFieldOfDefinition F xF := by
  intro hmin
  apply hE
  have := minField_le_of_factors F E xE xF h
  rw [hmin] at this
  exact top_le_iff.mp this

end MinimalField

/-! ## ★出典の紐付け(`.src`) -/

def minField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5, (i)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
