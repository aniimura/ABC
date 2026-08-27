/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricGroup
import ABC3.Meta.Claim

/-!
# **`APic(X)` —— 算術直線束の同型類がなす群**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★到達点 —— 原文の `APic(X)` が群であること

    `APicM X ≝ { (L̄, L̄⁻¹, L̄ ⊗ L̄⁻¹ ≅ Ō_X) } / 等長同型`

に **`CommGroup` のインスタンス**が入る。★積は**本物のテンソル積**であり
（`normOf_mul : |s ⊗ t| = |s| · |t|`）、単位元は `Ō_X` である。

★★2026-08-28 に `Definition 1.1` の項目全体の `.src` を下げた理由

> `TorsorMetric.base` は対象ごとの `Classical.choice` なので
> 群法則が計量のテンソル積を表さない

は、**(i) の側についてはこれで塞がった**。

## ★★★★★★★逆元は「データとして担ぐ」

本プロジェクトの `InvSheaf`（`PicSheafGroup.lean`）と同じ形である
——双対 `F^∨` を構成する代わりに、逆とその計量を**組にして持つ**。

★`mul` の逆は並べ替え

    `(A ⊗ B) ⊗ (A' ⊗ B') ≅ (A ⊗ A') ⊗ (B ⊗ B')`

で作る。★★その並べ替えは**結合律と交換律の等長版**（`AMetricGroup.lean`）を
繋いだものである。

## ★★★★★逆元は同型を除いて一意である

`L̄ ⊗ L̄⁻¹ ≅ Ō_X` かつ `L̄ ≅ M̄` なら `L̄⁻¹ ≅ M̄⁻¹`（`AInv.inv_isometric`）
——これがないと商の上で `⁻¹` が定義できない。

## ★残っている段（明示）

★★`Definition 1.1` の項目全体には (ii) の `deg_F`
（台帳 `arakelov-degF-finite-places`）も要る。
★★★`APicM` と既存の `APicOf`（捻れ集合表示の商）を繋ぐ段も別に要る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-! ## ★★同値関係の合同性 -/

theorem isometric_mul_left {L M N : AMetric X} (h : Isometric M N) :
    Isometric (L * M) (L * N) :=
  isometric_mul (isometric_refl L) h

theorem isometric_mul_right {L M N : AMetric X} (h : Isometric L M) :
    Isometric (L * N) (M * N) :=
  isometric_mul h (isometric_refl N)

/-- ★★★★**並べ替え** `(A ⊗ B) ⊗ (A' ⊗ B') ≅ (A ⊗ A') ⊗ (B ⊗ B')`。

★結合律と交換律の等長版を繋いだだけである。 -/
theorem isometric_rearrange (A B A' B' : AMetric X) :
    Isometric ((A * B) * (A' * B')) ((A * A') * (B * B')) := by
  refine isometric_trans (isometric_mul_assoc A B (A' * B')) ?_
  refine isometric_trans (isometric_mul_left (M := B * (A' * B')) (N := A' * (B * B')) ?_) ?_
  · refine isometric_trans (isometric_symm (isometric_mul_assoc B A' B')) ?_
    exact isometric_trans (isometric_mul_right (isometric_mul_comm B A'))
      (isometric_mul_assoc A' B B')
  · exact isometric_symm (isometric_mul_assoc A A' (B * B'))

/-! ## ★★★★★★逆元を担ぐ算術直線束 -/

/-- ★★★★★★**逆を持つ算術直線束**——`InvSheaf` と同じ形。 -/
structure AInv (X : Scheme.{0}) where
  /-- 台。 -/
  carrier : AMetric X
  /-- テンソル積についての逆。 -/
  inv : AMetric X
  /-- ★逆であること（等長同型として）。 -/
  isInv : Isometric (carrier * inv) 1

namespace AInv

/-- ★単位元 `Ō_X`。 -/
noncomputable def one (X : Scheme.{0}) : AInv X :=
  ⟨1, 1, isometric_mul_one 1⟩

/-- ★★逆元。 -/
noncomputable def symm (L : AInv X) : AInv X :=
  ⟨L.inv, L.carrier, isometric_trans (isometric_mul_comm L.inv L.carrier) L.isInv⟩

/-- ★★★★**積**——逆は並べ替えで出る。 -/
noncomputable def mul (L M : AInv X) : AInv X where
  carrier := L.carrier * M.carrier
  inv := L.inv * M.inv
  isInv := by
    refine isometric_trans (isometric_rearrange L.carrier M.carrier L.inv M.inv) ?_
    refine isometric_trans (isometric_mul L.isInv M.isInv) ?_
    exact isometric_mul_one 1

@[simp] theorem mul_carrier (L M : AInv X) :
    (L.mul M).carrier = L.carrier * M.carrier := rfl

@[simp] theorem one_carrier (X : Scheme.{0}) : (one X).carrier = 1 := rfl

@[simp] theorem symm_carrier (L : AInv X) : L.symm.carrier = L.inv := rfl

/-- ★★★★★★**逆元は等長同型を除いて一意である**。

★これがないと商の上で `⁻¹` が定義できない。 -/
theorem inv_isometric {L M : AInv X} (h : Isometric L.carrier M.carrier) :
    Isometric L.inv M.inv := by
  refine isometric_trans (isometric_symm (isometric_one_mul L.inv)) ?_
  refine isometric_trans (isometric_mul_right (isometric_symm M.isInv)) ?_
  refine isometric_trans (isometric_mul_right (isometric_mul_right (isometric_symm h))) ?_
  refine isometric_trans (isometric_mul_right (isometric_mul_comm L.carrier M.inv)) ?_
  refine isometric_trans (isometric_mul_assoc M.inv L.carrier L.inv) ?_
  exact isometric_trans (isometric_mul_left L.isInv) (isometric_mul_one M.inv)

/-- ★★台の等長同型による同値関係。 -/
def setoid (X : Scheme.{0}) : Setoid (AInv X) where
  r L M := Isometric L.carrier M.carrier
  iseqv := ⟨fun L => isometric_refl _, isometric_symm, isometric_trans⟩

end AInv

/-- ★★★★★★★★★★**`APic(X)`——算術直線束の同型類がなす群**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
def APicM (X : Scheme.{0}) : Type _ := Quotient (AInv.setoid X)

namespace APicM

/-- ★類を取る。 -/
def mk (L : AInv X) : APicM X := Quotient.mk _ L

theorem mk_eq_mk {L M : AInv X} (h : Isometric L.carrier M.carrier) :
    mk L = mk M := Quotient.sound h

noncomputable instance : One (APicM X) := ⟨mk (AInv.one X)⟩

noncomputable instance : Mul (APicM X) :=
  ⟨Quotient.map₂ AInv.mul (fun _ _ h₁ _ _ h₂ => isometric_mul h₁ h₂)⟩

noncomputable instance : Inv (APicM X) :=
  ⟨Quotient.map AInv.symm (fun _ _ h => AInv.inv_isometric h)⟩

@[simp] theorem mk_mul (L M : AInv X) : mk L * mk M = mk (L.mul M) := rfl

@[simp] theorem mk_inv (L : AInv X) : (mk L)⁻¹ = mk L.symm := rfl

theorem one_def (X : Scheme.{0}) : (1 : APicM X) = mk (AInv.one X) := rfl

/-- ★★★★★★★★★★**`APic(X)` は可換群である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★積は**本物のテンソル積**であり（`normOf_mul : |s ⊗ t| = |s| · |t|`）、
単位元は `Ō_X` である。 -/
noncomputable instance : CommGroup (APicM X) where
  mul_assoc := by
    refine Quotient.ind fun L => Quotient.ind fun M => Quotient.ind fun N => ?_
    exact mk_eq_mk (isometric_mul_assoc L.carrier M.carrier N.carrier)
  one_mul := by
    refine Quotient.ind fun L => ?_
    exact mk_eq_mk (isometric_one_mul L.carrier)
  mul_one := by
    refine Quotient.ind fun L => ?_
    exact mk_eq_mk (isometric_mul_one L.carrier)
  inv_mul_cancel := by
    refine Quotient.ind fun L => ?_
    exact mk_eq_mk L.symm.isInv
  mul_comm := by
    refine Quotient.ind fun L => Quotient.ind fun M => ?_
    exact mk_eq_mk (isometric_mul_comm L.carrier M.carrier)

end APicM

/-! ### ★出典の紐付け(`.src`) -/

def isometric_rearrange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)((A⊗B)⊗(A'⊗B') ≅ (A⊗A')⊗(B⊗B') の等長版)",
    sectionId := "genell-def-1-1-i" }

def AInv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(逆を持つ算術直線束)",
    sectionId := "genell-def-1-1-i" }

def AInv.inv_isometric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(逆元は等長同型を除いて一意であること)",
    sectionId := "genell-def-1-1-i" }

def APicM.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(APic(X)——算術直線束の同型類がなす可換群)",
    sectionId := "genell-def-1-1-i" }

def APicM.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isIsometry_mul / isIsometry_mul_assoc / isIsometry_mul_comm / isIsometry_mul_one"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isIsometry_mul") 3,
    .citation "[ABC3]" "LocalMetric.tensor(テンソル積の計量の構成)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.LocalMetric.tensor") 3,
    .citation "[ABC3]" "InvSheaf(逆をデータとして担ぐ形——層の側の前例)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.InvSheaf") 3,
    .implicitStep
      ("★★★★★2026-08-28 に Definition 1.1 の項目全体の .src を下げた理由" ++
       "(TorsorMetric.base が Classical.choice で群法則が計量のテンソル積を表さない)は、" ++
       "**(i) の側についてはこれで塞がった**——積は本物のテンソル積であり、" ++
       "|s ⊗ t| = |s| · |t| が成り立つ") 3,
    .implicitStep
      ("★★逆元は双対 F^∨ を構成する代わりに**データとして担いだ**" ++
       "(本プロジェクトの InvSheaf と同じ形)。" ++
       "★mul の逆は並べ替え (A⊗B)⊗(A'⊗B') ≅ (A⊗A')⊗(B⊗B') で作り、" ++
       "その並べ替えは結合律と交換律の等長版を繋いだものである") 3,
    .implicitStep
      ("★★★残っている段: Definition 1.1 の項目全体には (ii) の deg_F" ++
       "(台帳 arakelov-degF-finite-places)も要る。" ++
       "★APicM と既存の APicOf(捻れ集合表示の商)を繋ぐ段も別に要る") 3 ]

end ABC3.Found.Arakelov
