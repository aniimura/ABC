/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluTateC4Unit
import ABC3.Found.GaloisRep.LocalHeightDelta
import ABC3.Meta.Claim

/-!
# 第 1324 ブロック —— **Vélu の商の `c₄` は付値 `0`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か——第 1323 を付値の言葉に

第 1323 は「Tate 曲線の Vélu の商の `c₄` は `R` の単元」と言う。
★`R` の単元の付値は `0`（`tateDvrVal_eq_zero_of_isUnit`、在庫）なので、
分数体 `K` の付値の言葉に直せる。

☆これが第 1322（`c₄` が単元の整モデルは半安定）へ渡す一歩手前の形である。
-/

namespace ABC3.Found.GaloisRep

open Finset WeierstrassCurve ABC3.Skeleton.GenEll ABC3.Found.GenEll ABC3.Meta

variable {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★★★★★★★★★
**Vélu の商の `c₄` は付値 `0`**——★**無条件**（第 1324）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1323（`c₄` は `R` の単元）＋ `tateDvrVal_eq_zero_of_isUnit`（在庫）。 -/
theorem vAdd_c4_velu_tate_eq_zero
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R))) (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0)
    (h0 : algebraMap R K ((tateCurveAt q hq).c₄
      + 240 * (∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))) ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K ((tateCurveAt q hq).c₄
      + 240 * (∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))) h0) = 0 :=
  tateDvrVal_eq_zero_of_isUnit _ (isUnit_c4_velu_tate hl hζ hlu hu q hq hql h2 hDX) h0

/-! ## ★出典の紐付け(`.src`) -/

def vAdd_c4_velu_tate_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の c₄ は付値 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vAdd_c4_velu_tate_eq_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isUnit_c4_velu_tate(第 1323、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_c4_velu_tate") 1,
    .citation "[ABC3]" "tateDvrVal_eq_zero_of_isUnit(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDvrVal_eq_zero_of_isUnit") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1324）**——第 1322（`c₄` が単元の整モデルは半安定）へ" ++
       "渡す一歩手前の形である。") 2 ]

end ABC3.Found.GaloisRep
