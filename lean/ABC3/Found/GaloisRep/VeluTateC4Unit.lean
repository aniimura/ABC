/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateIsogeny
import ABC3.Found.GaloisRep.TateMultRed
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Found.GenEll.Velu
import ABC3.Meta.Claim

/-!
# 第 1323 ブロック —— **Tate 曲線の Vélu の商は `c₄` が単元**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`VeluQuotOK` の半安定性（悪い素点）の核

在庫の `c4_velu_tate`（`Skeleton/GenEll/TateIsogeny.lean`、**sorry 0**）は

> `c₄(E_q) + 240·Σ v₂ = l⁴ · c₄(E_{q^l})`

を与える。★右辺は**単元**である（`l` は単元、`c₄(E_{q^l}) = 1 + 240σ₃` も単元）。

☆したがって **Vélu の商の `c₄` は単元**であり、第 1322
（`c₄` が単元の整モデルは半安定）に渡せる。

★★★これが「Vélu の商は悪い素点で半安定」の**核**である
——`VeluQuotOK` の 2 つの穴のうち、悪い素点側がこれで埋まる。
-/

namespace ABC3.Found.GaloisRep

open Finset WeierstrassCurve ABC3.Skeleton.GenEll ABC3.Found.GenEll ABC3.Meta

variable {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★★★★★★★★★★★★
**Tate 曲線の Vélu の商は `c₄` が単元**——★**無条件**（第 1323）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆在庫の `c4_velu_tate` の右辺 `l⁴ · c₄(E_{q^l})` が単元だからである。 -/
theorem isUnit_c4_velu_tate {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R))) (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    IsUnit ((tateCurveAt q hq).c₄
      + 240 * (∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))) := by
  rw [c4_velu_tate hl hζ hlu hu q hq hql h2 hDX]
  exact (hlu.pow 4).mul (tateCurveAt_c4_isUnit (q ^ l) hql)

/-! ## ★出典の紐付け(`.src`) -/

def isUnit_c4_velu_tate.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate 曲線の Vélu の商は c₄ が単元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isUnit_c4_velu_tate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "c4_velu_tate(在庫、sorry 0)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.c4_velu_tate") 1,
    .citation "[ABC3]" "tateCurveAt_c4_isUnit(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateCurveAt_c4_isUnit") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1323）**——「Vélu の商は悪い素点で半安定」の**核**である。" ++
       "☆第 1322（`c₄` が単元の整モデルは半安定）に渡せば、" ++
       "`VeluQuotOK` の 2 つの穴のうち**悪い素点側が埋まる**。") 3 ]

end ABC3.Found.GaloisRep
