/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.LocalInputsSplit
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Meta.Claim

/-!
# 第 1318 ブロック —— **悪い素点の言葉で局所の 2 入力を出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——残り (ii) の解消

第 1317 は `l ∤ v(Q)` を `mkTateSetup` の言葉で受け取っていた。
★在庫の `not_dvd_vAdd_tateParam_of_not_dvd_jExp` が
**ちょうどその形**を `jExp` の言葉から与える（`PrimeToLocalHeights` そのもの）。

☆したがって仮説は原文の言葉になる:

* `jExp p E < 0`（悪い素点）
* `¬ l ∣ jExp p E`（`PrimeToLocalHeights`）
* `ζ` が `L_v` の原始 `l` 乗根
* `E ⊗ L_v` は極小かつ分裂乗法還元
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain NumberField
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**悪い素点の言葉で局所の 2 入力を出す**——★**無条件**（第 1318）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1317 に在庫の `not_dvd_vAdd_tateParam_of_not_dvd_jExp` を渡すだけである。 -/
theorem local_inputs_of_bad_prime {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) {l : ℕ} (hl : l.Prime) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {ζ : Lvˣ} (hζ : IsPrimitiveRoot ((ζ : Lv)) l) :
    (∃ P₀ : (E.baseChange Lv).toAffine.Point, addOrderOf P₀ = l) ∧
      (∀ T : Finset ((E.baseChange Lv).toAffine.Point),
        (∀ q ∈ T, l • q = 0) → T.card ≤ l) := by
  have hnd := not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp E h hj hcop
  exact local_inputs_of_split (R := R) (E.baseChange Lv) h hl hζ hnd

/-! ## ★出典の紐付け(`.src`) -/

def local_inputs_of_bad_prime.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(悪い素点の言葉で局所の 2 入力を出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def local_inputs_of_bad_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "not_dvd_vAdd_tateParam_of_not_dvd_jExp(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp") 1,
    .citation "[ABC3]" "local_inputs_of_split(第 1317、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.local_inputs_of_split") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1318）**——残り (ii) の解消である。" ++
       "☆仮説が原文の言葉（`jExp p E < 0`・`¬ l ∣ jExp p E` ＝ `PrimeToLocalHeights`）になった。") 2 ]

end ABC3.Found.GenEll
