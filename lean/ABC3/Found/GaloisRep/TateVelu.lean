/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateSpecialize
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.JScale
import ABC3.Found.GaloisRep.TateEquation
import ABC3.Meta.Claim

/-!
# Tate 曲線の Vélu の商 —— **変数変換は要らない**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★到達点

`Skeleton/GenEll/TateIsogeny.lean` の `tateModel_of_quot_mu`（`§9-1154`、第 717）が
固定した穴を、**2 つの係数の等式に還元する**。

★`tateCurveAt q` は `⟨1, 0, 0, a₄(q), a₆(q)⟩` であり、`veluCurve` は
`a₁, a₂, a₃` を変えない。★★したがって

> **`veluCurve (tateCurveAt q) v w = tateCurveAt (q^l)`**
>
> ⟺ `a₄(q^l) = a₄(q) − 5v` かつ `a₆(q^l) = a₆(q) − v − 7w`

（`b₂ = a₁² + 4a₂ = 1` だから `a₆` の補正項は `−v − 7w`）。

★★★☆**変数変換は要らない**——これで残る穴が `q`-展開の恒等式 2 本になった。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Found.GenEll

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★`tateCurveAt` の低い係数 -/

@[simp] theorem tateCurveAt_a₁ [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₁ = 1 := by
  simp [tateCurveAt, tateCurve, WeierstrassCurve.map]

@[simp] theorem tateCurveAt_a₂ [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₂ = 0 := by
  simp [tateCurveAt, tateCurve, WeierstrassCurve.map]

@[simp] theorem tateCurveAt_a₃ [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₃ = 0 := by
  simp [tateCurveAt, tateCurve, WeierstrassCurve.map]

/-- ★★★★★★**`b₂ = 1`**——`a₁ = 1`・`a₂ = 0` だから。 -/
@[simp] theorem tateCurveAt_b₂ [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).b₂ = 1 := by
  simp [WeierstrassCurve.b₂]

/-! ## ★★★★★★★★変数変換は要らない -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Tate 曲線の Vélu の商が Tate 曲線であることは係数 2 本に還元される**

    `veluCurve (tateCurveAt q) v w = tateCurveAt (q^l)`
      ⟺ `a₄(q^l) = a₄(q) − 5v` かつ `a₆(q^l) = a₆(q) − v − 7w`

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★☆**これで `Skeleton/GenEll/TateIsogeny.lean` の穴が
`q`-展開の恒等式 2 本になった**——変数変換を探す必要はない。 -/
theorem veluCurve_tateCurveAt_eq [IsAdicComplete I R] (q : R) (hq : q ∈ I)
    (q' : R) (hq' : q' ∈ I) (v w : R)
    (h4 : (tateCurveAt q' hq').a₄ = (tateCurveAt q hq).a₄ - 5 * v)
    (h6 : (tateCurveAt q' hq').a₆ = (tateCurveAt q hq).a₆ - v - 7 * w) :
    veluCurve (tateCurveAt q hq) v w = tateCurveAt q' hq' := by
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · show (tateCurveAt q hq).a₁ = (tateCurveAt q' hq').a₁
    rw [tateCurveAt_a₁, tateCurveAt_a₁]
  · show (tateCurveAt q hq).a₂ = (tateCurveAt q' hq').a₂
    rw [tateCurveAt_a₂, tateCurveAt_a₂]
  · show (tateCurveAt q hq).a₃ = (tateCurveAt q' hq').a₃
    rw [tateCurveAt_a₃, tateCurveAt_a₃]
  · show (tateCurveAt q hq).a₄ - 5 * v = (tateCurveAt q' hq').a₄
    rw [h4]
  · show (tateCurveAt q hq).a₆ - (tateCurveAt q hq).b₂ * v - 7 * w
      = (tateCurveAt q' hq').a₆
    rw [h6, tateCurveAt_b₂, one_mul]

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

/-! ## ★★★★★★★★★★★★Tate 曲線上の Vélu の量 -/

/-- ★★★★★★★★**Tate 曲線では `g^x_Q = 3x² + a₄ − y`**（`a₁ = 1`, `a₂ = 0`）。 -/
@[simp] theorem veluGx_tateCurveAt [IsAdicComplete I R] (q : R) (hq : q ∈ I) (x y : R) :
    veluGx (tateCurveAt q hq) x y = 3 * x ^ 2 + (tateCurveAt q hq).a₄ - y := by
  rw [veluGx, tateCurveAt_a₁, tateCurveAt_a₂]
  ring

/-- ★★★★★★★★**`v_Q = g^x_Q`**（`veluV2` の形）。 -/
@[simp] theorem veluV2_tateCurveAt [IsAdicComplete I R] (q : R) (hq : q ∈ I) (x y : R) :
    veluV2 (tateCurveAt q hq) x y = 3 * x ^ 2 + (tateCurveAt q hq).a₄ - y := by
  rw [veluV2, veluGx_tateCurveAt]

/-- ★★★★★★★★**Tate 曲線では `g^y_Q = −2y − x`**。 -/
@[simp] theorem veluGy_tateCurveAt [IsAdicComplete I R] (q : R) (hq : q ∈ I) (x y : R) :
    veluGy (tateCurveAt q hq) x y = -2 * y - x := by
  rw [veluGy, tateCurveAt_a₁, tateCurveAt_a₃]
  ring

/-- ★★★★★★★★**`u_Q = (2y + x)²`**。 -/
@[simp] theorem veluU_tateCurveAt [IsAdicComplete I R] (q : R) (hq : q ∈ I) (x y : R) :
    veluU (tateCurveAt q hq) x y = (2 * y + x) ^ 2 := by
  rw [veluU, veluGy_tateCurveAt]
  ring

/-! ## ★★★★★★★★★★★★★★★★★★★★`j` の一致へ -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`c₄`・`c₆` の尺度関係から `j` の一致を出す**。

★★★これが葉 1（`jExp_velu_bad`）が実際に消費する形である。
☆仮説の 2 本は `Skeleton/GenEll/TateIsogeny.lean` の
`c4_velu_tate`・`c6_velu_tate`（第 837、数値確認済み）である。 -/
theorem j_velu_tate_eq {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] (q : R) (hq : q ∈ I) (l : ℕ) (hql : q ^ l ∈ I) (v w : R)
    [(veluCurve (tateCurveAt q hq) v w).IsElliptic] [(tateCurveAt (q ^ l) hql).IsElliptic]
    (h4 : (tateCurveAt q hq).c₄ + 240 * v = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄)
    (h6 : (tateCurveAt q hq).c₆ + 504 * v + 6048 * w
      = (l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆) :
    (veluCurve (tateCurveAt q hq) v w).j = (tateCurveAt (q ^ l) hql).j := by
  refine ABC3.Found.GenEll.j_eq_of_c4_c6_scale_pos _ _ ((l : R)) ?_ ?_
  · rw [veluCurve_c₄]
    exact h4
  · rw [veluCurve_c₆, tateCurveAt_b₂]
    rw [mul_one]
    exact h6

def j_velu_tate_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₄・c₆ の尺度関係から j の一致。★無条件)",
    sectionId := "genell-lemma-3-2" }

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

def veluGx_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate 曲線では g^x_Q = 3x² + a₄ − y。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluV2_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_Q = g^x_Q の Tate 形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluGy_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate 曲線では g^y_Q = −2y − x。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluU_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(u_Q = (2y + x)² の Tate 形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluCurve_tateCurveAt_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の Vélu の商——変数変換は要らない。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateCurveAt_b₂.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の b₂ = 1。★無条件)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
