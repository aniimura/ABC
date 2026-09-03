/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateSpecialize
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.JScale
import ABC3.Found.GaloisRep.TateEquation
import ABC3.Meta.Claim
import ABC3.Found.GaloisRep.TateVelu.Lemma32

/-!
# TateVelu —— `[GenEll] Definition 3.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Found.GenEll
variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★★★★★★★★★**同じことを `veluQuotientFull` の形で**。

☆`veluQuotientFull W S = veluCurve W (veluVFull W S) (veluWFull W S)` なので、
`v`・`w` を `μ_l` の点にわたる Vélu の和に取れば上の補題がそのまま効く。 -/
theorem veluCurve_tateCurveAt_eq' [IsAdicComplete I R] (q : R) (hq : q ∈ I)
    (l : ℕ) (hql : q ^ l ∈ I) (v w : R)
    (h4 : (tateCurveAt (q ^ l) hql).a₄ = (tateCurveAt q hq).a₄ - 5 * v)
    (h6 : (tateCurveAt (q ^ l) hql).a₆ = (tateCurveAt q hq).a₆ - v - 7 * w) :
    veluCurve (tateCurveAt q hq) v w = tateCurveAt (q ^ l) hql :=
  veluCurve_tateCurveAt_eq q hq (q ^ l) hql v w h4 h6

/-- ★★★★★★**`a₄(q) = evalAdic tateA4 q`**。 -/
theorem tateCurveAt_a₄ [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₄ = evalAdic tateA4 q hq := by
  simp [tateCurveAt, tateCurve, WeierstrassCurve.map, evalAdicHom]

/-- ★★★★★★**`a₆(q) = evalAdic tateA6 q`**。 -/
theorem tateCurveAt_a₆ [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₆ = evalAdic tateA6 q hq := by
  simp [tateCurveAt, tateCurve, WeierstrassCurve.map, evalAdicHom]

/-! ## ★★★★★★★★★★★★★★★★★★★★Weierstrass の ODE の形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`(2Y + X)² = 4X³ + X² + 4a₄X + 4a₆`**——Tate 曲線の方程式の書き換え。

★★★`veluU (tateCurveAt q) X Y = (2Y+X)²` （第 810）なので、
**Vélu の `u_Q` が `X` の三次式で書ける**。

☆さらに `D = u∂_u` と置くと `DX = 2Y + X` なので、
これは Weierstrass の ODE `(DX)² = 4X³ + X² + 4a₄X + 4a₆` であり、
微分すると `D²X = 6X² + X + 2a₄` が出る——
★★これが `∑_ζ X²` を**畳み込みを経由せずに**与える道である
（`Skeleton/GenEll/SigmaConvolution.lean` の第 842 の記述）。 -/
theorem two_y_add_x_sq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    (2 * tateYpair a w q hq + tateXpair a w q hq) ^ 2
      = 4 * tateXpair a w q hq ^ 3 + tateXpair a w q hq ^ 2
        + 4 * (tateCurveAt q hq).a₄ * tateXpair a w q hq + 4 * (tateCurveAt q hq).a₆ := by
  have h := tate_equation a w q hq haw ha hw
  linear_combination 4 * h

def two_y_add_x_sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3((2Y+X)² = 4X³ + X² + 4a₄X + 4a₆。★無条件)",
    sectionId := "genell-def-3-3" }

/-! ## ★出典の紐付け(`.src`) -/

def tateCurveAt_a₄.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(a₄(q) = evalAdic tateA4 q。★無条件)",
    sectionId := "genell-def-3-3" }

def tateCurveAt_a₆.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(a₆(q) = evalAdic tateA6 q。★無条件)",
    sectionId := "genell-def-3-3" }

def tateCurveAt_b₂.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の b₂ = 1。★無条件)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
