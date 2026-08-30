/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateSpecialize
import ABC3.Found.GenEll.Velu
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

/-! ## ★出典の紐付け(`.src`) -/

def veluCurve_tateCurveAt_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の Vélu の商——変数変換は要らない。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateCurveAt_b₂.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の b₂ = 1。★無条件)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
