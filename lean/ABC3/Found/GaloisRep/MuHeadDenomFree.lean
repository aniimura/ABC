/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuGraded
import ABC3.Found.GaloisRep.TateDSeries
import ABC3.Meta.Claim

/-!
# 第 1102 ブロック —— **分母を持たない頭項**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

第 1101 で特定した機構——在庫の `one_sub_mul_sum_nsmul`——を使う。
`η^l = 1` かつ `∑_{k<l} η^k = 0` のとき、**`IsUnit` を一切使わずに**

    `(1 − η) · S(η) = −l`,   `S(η) ≔ ∑_{k<l} k·η^k`

が成り立つ。☆すなわち `S` は `−l/(1−η)` の役をする**多項式**である。

★したがって「`l^k` を掛けた頭項」を **`S` の多項式として無条件に定義できる**。
`1 − η` が単元かどうかに依らない。

| 元の項 | 分母 | `E` 版（本ファイル） |
|---|---|---|
| `tateXterm η = η·inv(1−η)²` | 2 | `η·S²` |
| `tateYterm η = η²·inv(1−η)³` | 3 | `−η²·S³` |
| `tateDXterm η = η(1+η)·inv(1−η)³` | 3 | `−η(1+η)·S³` |
| `tateDYterm η = η²(2+η)·inv(1−η)⁴` | 4 | `η²(2+η)·S⁴` |
| `tateD2Xterm η = η(1+4η+η²)·inv(1−η)⁴` | 4 | `η(1+4η+η²)·S⁴` |
| `tateD3Xterm η = η(1+11η+11η²+η³)·inv(1−η)⁵` | 5 | `−η(1+11η+11η²+η³)·S⁵` |

☆符号は `(−1)^k`（`inv(1−η) = −S/l` の `k` 乗）。

★★本ファイルは**定義と橋**（`hu` の下で `l^k` 倍に等しいこと）を置く。
`hu` を実際に外すのは、これを使って `sum_mu_*` を書き直す次の段である。
-/

namespace ABC3.Found.GaloisRep

open Finset ABC3.Meta

variable {R : Type} [CommRing R]

/-- ★★★★**`S(η) = ∑_{k<l} k·η^k`** —— `−l/(1−η)` の役をする多項式。 -/
noncomputable def muS (l : ℕ) (η : R) : R := ∑ k ∈ range l, (k : R) * η ^ k

/-- ★★★★★★**`(1 − η)·S(η) = −l`** —— `IsUnit` を要さない。 -/
theorem one_sub_mul_muS [IsDomain R] {l : ℕ} {η : R} (hpow : η ^ l = 1)
    (hsum : ∑ k ∈ range l, η ^ k = 0) : (1 - η) * muS l η = -(l : R) :=
  one_sub_mul_sum_nsmul hpow hsum

/-! ## ★★★★分母を持たない頭項 -/

/-- ★★`l²·tateXterm η`。 -/
noncomputable def tateXtermE (l : ℕ) (η : R) : R := η * muS l η ^ 2

/-- ★★`l³·tateYterm η`。 -/
noncomputable def tateYtermE (l : ℕ) (η : R) : R := -(η ^ 2 * muS l η ^ 3)

/-- ★★`l³·tateDXterm η`。 -/
noncomputable def tateDXtermE (l : ℕ) (η : R) : R := -(η * (1 + η) * muS l η ^ 3)

/-- ★★`l⁴·tateDYterm η`。 -/
noncomputable def tateDYtermE (l : ℕ) (η : R) : R := η ^ 2 * (2 + η) * muS l η ^ 4

/-- ★★`l⁴·tateD2Xterm η`。 -/
noncomputable def tateD2XtermE (l : ℕ) (η : R) : R :=
  η * (1 + 4 * η + η ^ 2) * muS l η ^ 4

/-- ★★`l⁵·tateD3Xterm η`。 -/
noncomputable def tateD3XtermE (l : ℕ) (η : R) : R :=
  -(η * (1 + 11 * η + 11 * η ^ 2 + η ^ 3) * muS l η ^ 5)

/-! ## ★★★★★★★★橋——`hu` の下で `l^k` 倍に等しい -/

section Bridge

variable [IsDomain R] {l : ℕ} {η : R}

/-- ☆`inv(1−η) = −inv(l)·S(η)`（在庫の言い換え）。 -/
private theorem inv_eq (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    Ring.inverse (1 - η) = -(Ring.inverse ((l : R))) * muS l η :=
  ringInverse_one_sub_eq hlu hu hpow hsum

/-- ☆`l^k · inv(1−η)^k = (−1)^k · S(η)^k`。 -/
private theorem pow_bridge (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) (k : ℕ) :
    (l : R) ^ k * Ring.inverse (1 - η) ^ k = (-1) ^ k * muS l η ^ k := by
  have hc : ((l : R)) * Ring.inverse ((l : R)) = 1 := Ring.mul_inverse_cancel _ hlu
  rw [inv_eq hlu hu hpow hsum]
  have h1 : (l : R) ^ k * (-(Ring.inverse ((l : R))) * muS l η) ^ k
      = ((-1 : R)) ^ k * ((l : R) * Ring.inverse ((l : R))) ^ k * muS l η ^ k := by
    rw [mul_pow, mul_pow, neg_pow]
    ring
  rw [h1, hc, one_pow, mul_one]

/-- ★★★★★★**`X` の橋**。 -/
theorem natCast_pow_mul_tateXterm (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (l : R) ^ 2 * tateXterm η = tateXtermE l η := by
  rw [tateXterm, tateXtermE, show (l : R) ^ 2 * (η * Ring.inverse (1 - η) ^ 2)
      = η * ((l : R) ^ 2 * Ring.inverse (1 - η) ^ 2) from by ring,
    pow_bridge hlu hu hpow hsum 2]
  ring

/-- ★★★★★★**`Y` の橋**。 -/
theorem natCast_pow_mul_tateYterm (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (l : R) ^ 3 * tateYterm η = tateYtermE l η := by
  rw [tateYterm, tateYtermE, show (l : R) ^ 3 * (η ^ 2 * Ring.inverse (1 - η) ^ 3)
      = η ^ 2 * ((l : R) ^ 3 * Ring.inverse (1 - η) ^ 3) from by ring,
    pow_bridge hlu hu hpow hsum 3]
  ring

/-- ★★★★★★**`DX` の橋**。 -/
theorem natCast_pow_mul_tateDXterm (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (l : R) ^ 3 * tateDXterm η = tateDXtermE l η := by
  rw [tateDXterm, tateDXtermE,
    show (l : R) ^ 3 * (η * (1 + η) * Ring.inverse (1 - η) ^ 3)
      = η * (1 + η) * ((l : R) ^ 3 * Ring.inverse (1 - η) ^ 3) from by ring,
    pow_bridge hlu hu hpow hsum 3]
  ring

/-- ★★★★★★**`DY` の橋**。 -/
theorem natCast_pow_mul_tateDYterm (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (l : R) ^ 4 * tateDYterm η = tateDYtermE l η := by
  rw [tateDYterm, tateDYtermE,
    show (l : R) ^ 4 * (η ^ 2 * (2 + η) * Ring.inverse (1 - η) ^ 4)
      = η ^ 2 * (2 + η) * ((l : R) ^ 4 * Ring.inverse (1 - η) ^ 4) from by ring,
    pow_bridge hlu hu hpow hsum 4]
  ring

/-- ★★★★★★★★**`D²X` の橋**。 -/
theorem natCast_pow_mul_tateD2Xterm (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (l : R) ^ 4 * tateD2Xterm η = tateD2XtermE l η := by
  rw [tateD2Xterm, tateD2XtermE,
    show (l : R) ^ 4 * (η * (1 + 4 * η + η ^ 2) * Ring.inverse (1 - η) ^ 4)
      = η * (1 + 4 * η + η ^ 2) * ((l : R) ^ 4 * Ring.inverse (1 - η) ^ 4) from by ring,
    pow_bridge hlu hu hpow hsum 4]
  ring

/-- ★★★★★★★★**`D³X` の橋**。 -/
theorem natCast_pow_mul_tateD3Xterm (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - η))
    (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (l : R) ^ 5 * tateD3Xterm η = tateD3XtermE l η := by
  rw [tateD3Xterm, tateD3XtermE,
    show (l : R) ^ 5 * (η * (1 + 11 * η + 11 * η ^ 2 + η ^ 3) * Ring.inverse (1 - η) ^ 5)
      = η * (1 + 11 * η + 11 * η ^ 2 + η ^ 3)
          * ((l : R) ^ 5 * Ring.inverse (1 - η) ^ 5) from by ring,
    pow_bridge hlu hu hpow hsum 5]
  ring

end Bridge

/-! ## ★★★★★★★★★★★★第 1111 —— `E` 版は環準同型と無条件に可換

★★★★**これが `hu` を外す鍵である。**

☆在庫の `map_tateD2Xterm`（`MuDYSum.lean`）は `IsUnit (1 - t)` を要求する——
`f (Ring.inverse x) = Ring.inverse (f x)` が単元性を要するからである。
★`E` 版は **`Ring.inverse` を含まない多項式**なので、
任意の環準同型と**無条件に**可換である。 -/

section Map

variable {A B : Type} [CommRing A] [CommRing B]

/-- ★☆`muS` は環準同型で写る。 -/
theorem map_muS (f : A →+* B) (l : ℕ) (t : A) : f (muS l t) = muS l (f t) := by
  rw [muS, muS, map_sum]
  exact Finset.sum_congr rfl (fun k _ => by rw [map_mul, map_pow, map_natCast])

/-- ★★★★★★**`tateD2XtermE` は環準同型で写る**——`IsUnit` 不要。 -/
theorem map_tateD2XtermE (f : A →+* B) (l : ℕ) (t : A) :
    f (tateD2XtermE l t) = tateD2XtermE l (f t) := by
  rw [tateD2XtermE, tateD2XtermE, map_mul, map_mul, map_pow, map_muS, map_add, map_add,
    map_one, map_mul, map_pow]
  simp [map_ofNat]

/-- ★★★★★★**`tateDXtermE` は環準同型で写る**——`IsUnit` 不要。 -/
theorem map_tateDXtermE (f : A →+* B) (l : ℕ) (t : A) :
    f (tateDXtermE l t) = tateDXtermE l (f t) := by
  rw [tateDXtermE, tateDXtermE, map_neg]
  congr 1
  rw [map_mul, map_mul, map_pow, map_muS, map_add, map_one]

/-- ★★★★★★**`tateD3XtermE` は環準同型で写る**——`IsUnit` 不要。 -/
theorem map_tateD3XtermE (f : A →+* B) (l : ℕ) (t : A) :
    f (tateD3XtermE l t) = tateD3XtermE l (f t) := by
  rw [tateD3XtermE, tateD3XtermE, map_neg]
  congr 1
  rw [map_mul, map_mul, map_pow, map_muS, map_add, map_add, map_add, map_one,
    map_mul, map_mul, map_pow, map_pow]
  simp [map_ofNat]

/-- ★★★★★★**`tateXtermE` は環準同型で写る**——`IsUnit` 不要（第 1125）。 -/
theorem map_tateXtermE (f : A →+* B) (l : ℕ) (t : A) :
    f (tateXtermE l t) = tateXtermE l (f t) := by
  rw [tateXtermE, tateXtermE, map_mul, map_pow, map_muS]

/-- ★★★★★★**`tateYtermE` は環準同型で写る**——`IsUnit` 不要（第 1125）。 -/
theorem map_tateYtermE (f : A →+* B) (l : ℕ) (t : A) :
    f (tateYtermE l t) = tateYtermE l (f t) := by
  rw [tateYtermE, tateYtermE, map_neg]
  congr 1
  rw [map_mul, map_pow, map_pow, map_muS]

/-- ★★★★★★**`tateDYtermE` は環準同型で写る**——`IsUnit` 不要（第 1125）。 -/
theorem map_tateDYtermE (f : A →+* B) (l : ℕ) (t : A) :
    f (tateDYtermE l t) = tateDYtermE l (f t) := by
  rw [tateDYtermE, tateDYtermE, map_mul, map_mul, map_pow, map_pow, map_muS, map_add]
  simp [map_ofNat]

end Map

/-! ## ★出典の紐付け(`.src`) -/

def muS.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(S(η) = ∑ k·η^k は −l/(1−η) の役をする多項式)",
    sectionId := "genell-lemma-3-5" }

def tateXtermE.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(分母を持たない頭項——l^k を掛けた形を S で定義する)",
    sectionId := "genell-lemma-3-5" }

def natCast_pow_mul_tateD2Xterm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(D²X の橋——hu の下で l⁴ 倍が E 版に等しい)",
    sectionId := "genell-lemma-3-5" }

def muS.needs : List ProofObligation :=
  [ .citation "[ABC3]" "one_sub_mul_sum_nsmul((1−η)·S = −l、在庫、IsUnit 不要)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.one_sub_mul_sum_nsmul") 1,
    .implicitStep
      ("★★**2026-09-01（第 1102）**——`E` 版の定義は `IsUnit` を含まないので" ++
       "`p ∣ l` でも意味を持つ。☆残るのは `sum_mu_*` を `E` 版で述べ直し、" ++
       "その証明から `hu` を落とすことである。") 25 ]

end ABC3.Found.GaloisRep
