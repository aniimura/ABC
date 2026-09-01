/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuPairDenomFree
import ABC3.Found.GaloisRep.TateOrigin
import ABC3.Meta.Claim

/-!
# 第 1133 ブロック —— **`K` 水準の座標と分母なし版を繋ぐ橋**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か（節点 5 の消費者 B・C）

`pointCoords_tatePtPair`（`TateVeluPoints.lean:42`）は
`IsUnit (1 − a)` を仮説に取り、`μ_l` の点の座標を `algebraMap R K (tateXpair …)` と同定する。
★`p ∣ l` では `1 − ζ^i` が `R` の単元でないのでこの形は使えない。

☆しかし **`K` 水準の座標 `tateXK`・`tateYK` は最初から分母を払った形で定義されている**:

    `tateXK a w q hq = φ(tateXpairE a w q hq) · φ(1 − a)⁻²`
    `tateXpairE a w q hq = a + (1 − a)²·(尾)`

★要るのは「`R` の単元性」ではなく「`φ(1 − a) ≠ 0`」だけである。

## ★★★★★★★★★★橋の中身

`tateXpairDF`（第 1118）は同じ「尾」を `l²` 倍で持つ:

    `tateXpairDF l a w q hq = a·S(a)² + l²·(尾)`,   `S(a) = ∑_{k<l} k·a^k`

`(1 − a)·S(a) = −l`（`one_sub_mul_muS`、`IsUnit` 不要）を使うと

    `(1 − a)²·tateXpairDF = l²·tateXpairE`   （★`R` の中の恒等式、`IsUnit` 不要）

になり、`φ(1 − a)` で割れば

    `l²·tateXK = φ(tateXpairDF)`   （★`K` の中、`IsUnit` 不要）

が出る。☆`Y` 側は `(1 − a)³` と `l³` で同じ形である。
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {R K : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  [Field K] [Algebra R K]

/-! ## ★★★★★★`R` の中の恒等式 -/

/-- ★★★★★★★★**`(1 − a)²·tateXpairDF = l²·tateXpairE`**——★`IsUnit` 不要（第 1133）。

☆`(1 − a)·S(a) = −l`（`one_sub_mul_muS`）だけで出る。 -/
theorem one_sub_sq_mul_tateXpairDF {l : ℕ} {a : R} (hpow : a ^ l = 1)
    (hsum : ∑ k ∈ range l, a ^ k = 0) (w q : R) (hq : q ∈ I) :
    (1 - a) ^ 2 * tateXpairDF l a w q hq = (l : R) ^ 2 * tateXpairE a w q hq := by
  have hS : (1 - a) * muS l a = -(l : R) := one_sub_mul_muS hpow hsum
  rw [tateXpairDF, tateXpairE, tateXtermE]
  linear_combination (a * ((1 - a) * muS l a - (l : R))) * hS

/-- ★★★★★★★★**`(1 − a)³·tateYpairDF = l³·tateYpairE`**——★`IsUnit` 不要（第 1133）。 -/
theorem one_sub_cube_mul_tateYpairDF {l : ℕ} {a : R} (hpow : a ^ l = 1)
    (hsum : ∑ k ∈ range l, a ^ k = 0) (w q : R) (hq : q ∈ I) :
    (1 - a) ^ 3 * tateYpairDF l a w q hq = (l : R) ^ 3 * tateYpairE a w q hq := by
  have hS : (1 - a) * muS l a = -(l : R) := one_sub_mul_muS hpow hsum
  rw [tateYpairDF, tateYpairE, tateYtermE]
  linear_combination
    (-(a ^ 2) * (((1 - a) * muS l a) ^ 2 - (1 - a) * muS l a * (l : R) + (l : R) ^ 2)) * hS

/-! ## ★★★★★★★★★★`K` の中の橋 -/

/-- ★★★★★★★★★★★★★★★★
**`l²·tateXK = φ(tateXpairDF)`**（第 1133）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★**`IsUnit (1 − a)` を仮説に置いていない**——要るのは `φ(1 − a) ≠ 0` だけで、
`a = ζ^i`（`0 < i < l`）なら `R → K` の単射性から自動である。
☆これが `p ∣ l` で `pointCoords_tatePtPair` の代わりになる形の核である。 -/
theorem natCast_pow_mul_tateXK {l : ℕ} {a : R} (hpow : a ^ l = 1)
    (hsum : ∑ k ∈ range l, a ^ k = 0)
    (hne : algebraMap R K (1 - a) ≠ 0) (w q : R) (hq : q ∈ I) :
    (l : K) ^ 2 * (tateXK (I := I) a w q hq : K)
      = algebraMap R K (tateXpairDF l a w q hq) := by
  have h := congrArg (algebraMap R K) (one_sub_sq_mul_tateXpairDF hpow hsum w q hq)
  rw [map_mul, map_mul, map_pow, map_pow, map_natCast] at h
  rw [tateXK]
  field_simp
  linear_combination -h

/-- ★★★★★★★★★★★★★★★★**`l³·tateYK = φ(tateYpairDF)`**（第 1133）。 -/
theorem natCast_pow_mul_tateYK {l : ℕ} {a : R} (hpow : a ^ l = 1)
    (hsum : ∑ k ∈ range l, a ^ k = 0)
    (hne : algebraMap R K (1 - a) ≠ 0) (w q : R) (hq : q ∈ I) :
    (l : K) ^ 3 * (tateYK (I := I) a w q hq : K)
      = algebraMap R K (tateYpairDF l a w q hq) := by
  have h := congrArg (algebraMap R K) (one_sub_cube_mul_tateYpairDF hpow hsum w q hq)
  rw [map_mul, map_mul, map_pow, map_pow, map_natCast] at h
  rw [tateYK]
  field_simp
  linear_combination -h

/-! ## ★出典の紐付け(`.src`) -/

def natCast_pow_mul_tateXK.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(K 水準の X と分母なし版の橋。★IsUnit (1−a) を取らない)",
    sectionId := "genell-def-3-3" }

def one_sub_sq_mul_tateXpairDF.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3((1−a)²·tateXpairDF = l²·tateXpairE。★無条件)",
    sectionId := "genell-def-3-3" }

def natCast_pow_mul_tateXK.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "one_sub_mul_muS((1−η)·S(η) = −l。IsUnit 不要、第 1102、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.one_sub_mul_muS") 1,
    .citation "[ABC3]" "tateXK(K 水準の X。分母を払った形で定義されている、在庫)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateXK") 1,
    .implicitStep
      ("★★**2026-09-01（第 1133）の測定**——`tateXK` は最初から " ++
       "`φ(tateXpairE)·φ(1−a)⁻²` として定義されており、" ++
       "**`R` の単元性を要求していない**（要るのは `φ(1−a) ≠ 0` だけ）。" ++
       "☆`tateXpairDF` は同じ尾を `l²` 倍で持つので、" ++
       "`(1−a)²·tateXpairDF = l²·tateXpairE` が `R` の中で成り立ち、" ++
       "`φ(1−a)` で割れば `K` の中の橋になる。" ++
       "★これが節点 5 の消費者 B・C（`pointCoords_tatePtPair` 経由）を外す核である。") 6 ]

end ABC3.Found.GaloisRep
