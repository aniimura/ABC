/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluBadPrimeAssembly
import ABC3.Meta.Claim

/-!
# 第 1435 ブロック —— **`p ∣ l` かつ `l ≥ 5` なら `p ∤ 6`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か——第 1425-1431 が要る `p ∤ 6` を作る

第 1425（`semistableAt_of_valAdd_c4_c6_scaled`）と第 1430-1431 は
`valAdd p (48) = 0`・`valAdd p (864) = 0`・`valAdd p (1728) = 0`、
すなわち **`p ∤ 6`** を要求する。

☆`p ∣ l`（＝`¬ IsUnit (l : 𝒪_p)`）で `l` が `5` 以上の素数なら、それは自動である:
`l` と `48 = 2⁴·3` は互いに素なので、Bézout `a·l + b·48 = 1` を `𝒪_p` に落とすと
**両方が非単元だと `1 ∈ 𝔪` になってしまう**。

| 定理 | 内容 |
|---|---|
| `isUnit_natCast_of_coprime` | ★★`l` が非単元で `gcd(l, n) = 1` なら `n` は単元 |
| `not_dvd_of_five_le` | ☆`l ≥ 5` の素数は `2^a·3^b` を割らない |
| `valAdd_48_eq_zero` ほか | ★★★★`p ∣ l`・`l ≥ 5` から `v_p(48) = v_p(864) = v_p(1728) = 0` |
-/

namespace ABC3.Found.GenEll

open IsDedekindDomain NumberField ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★**非単元 `l` と互いに素な `n` は単元**——Bézout を `𝒪_p` に落とすだけ。 -/
theorem isUnit_natCast_of_coprime (p : HeightOneSpectrum (𝓞 L)) {l n : ℕ}
    (hcop : Nat.Coprime l n) (hl : ¬ IsUnit ((l : primeSubring p))) :
    IsUnit ((n : primeSubring p)) := by
  by_contra hn
  have hlm : ((l : primeSubring p)) ∈ IsLocalRing.maximalIdeal (primeSubring p) :=
    (IsLocalRing.mem_maximalIdeal _).2 hl
  have hnm : ((n : primeSubring p)) ∈ IsLocalRing.maximalIdeal (primeSubring p) :=
    (IsLocalRing.mem_maximalIdeal _).2 hn
  have hb : (1 : ℤ) = l * Nat.gcdA l n + n * Nat.gcdB l n := by
    have h := Nat.gcd_eq_gcd_ab l n
    rw [hcop] at h
    exact_mod_cast h
  have h1 : (1 : primeSubring p)
      = (l : primeSubring p) * ((Nat.gcdA l n : ℤ) : primeSubring p)
        + (n : primeSubring p) * ((Nat.gcdB l n : ℤ) : primeSubring p) := by
    have hh := congrArg (fun z : ℤ => ((z : primeSubring p))) hb
    push_cast at hh ⊢
    exact hh
  have hone : (1 : primeSubring p) ∈ IsLocalRing.maximalIdeal (primeSubring p) := by
    rw [h1]
    exact Ideal.add_mem _ (Ideal.mul_mem_right _ _ hlm) (Ideal.mul_mem_right _ _ hnm)
  exact (Ideal.ne_top_iff_one _).1 (IsLocalRing.maximalIdeal.isMaximal _).ne_top hone

/-- ☆`5` 以上の素数は `2^a·3^b` を割らない。 -/
theorem not_dvd_of_five_le {l n : ℕ} (hl : l.Prime) (hl5 : 5 ≤ l) (a b : ℕ)
    (hn : n = 2 ^ a * 3 ^ b) : ¬ l ∣ n := by
  intro hd
  rw [hn] at hd
  rcases (Nat.Prime.dvd_mul hl).1 hd with h | h
  · have h2 := Nat.Prime.dvd_of_dvd_pow hl h
    have := Nat.le_of_dvd (by norm_num) h2
    omega
  · have h3 := Nat.Prime.dvd_of_dvd_pow hl h
    have := Nat.le_of_dvd (by norm_num) h3
    omega

/-- ☆`p ∣ l`・`l ≥ 5` から `v_p(n) = 0`（`n = 2^a·3^b`）。 -/
theorem valAdd_natCast_eq_zero_of_five_le (p : HeightOneSpectrum (𝓞 L)) {l : ℕ}
    (hl : l.Prime) (hl5 : 5 ≤ l) (hnu : ¬ IsUnit ((l : primeSubring p)))
    {n : ℕ} (a b : ℕ) (hn : n = 2 ^ a * 3 ^ b) (hne : ((n : ℕ) : L) ≠ 0) :
    valAdd p (Units.mk0 (((n : ℕ) : L)) hne) = 0 :=
  ABC3.Found.GenEll.valAdd_natCast_eq_zero_of_isUnit p hne
    (isUnit_natCast_of_coprime p
      ((Nat.Prime.coprime_iff_not_dvd hl).2 (not_dvd_of_five_le hl hl5 a b hn)) hnu)

/-- ★★★★**`p ∣ l`・`l ≥ 5` なら `v_p(48) = 0`**。 -/
theorem valAdd_48_eq_zero (p : HeightOneSpectrum (𝓞 L)) {l : ℕ} (hl : l.Prime)
    (hl5 : 5 ≤ l) (hnu : ¬ IsUnit ((l : primeSubring p))) :
    valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0 := by
  have hne : ((48 : ℕ) : L) ≠ 0 := by norm_num
  have h := valAdd_natCast_eq_zero_of_five_le p hl hl5 hnu 4 1 (by norm_num) hne
  have hEq : Units.mk0 (((48 : ℕ) : L)) hne = Units.mk0 (48 : L) (by norm_num) := by
    refine Units.ext ?_
    push_cast
    rfl
  rwa [hEq] at h

/-- ★★★★**`p ∣ l`・`l ≥ 5` なら `v_p(864) = 0`**。 -/
theorem valAdd_864_eq_zero (p : HeightOneSpectrum (𝓞 L)) {l : ℕ} (hl : l.Prime)
    (hl5 : 5 ≤ l) (hnu : ¬ IsUnit ((l : primeSubring p))) :
    valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0 := by
  have hne : ((864 : ℕ) : L) ≠ 0 := by norm_num
  have h := valAdd_natCast_eq_zero_of_five_le p hl hl5 hnu 5 3 (by norm_num) hne
  have hEq : Units.mk0 (((864 : ℕ) : L)) hne = Units.mk0 (864 : L) (by norm_num) := by
    refine Units.ext ?_
    push_cast
    rfl
  rwa [hEq] at h

/-- ★★★★**`p ∣ l`・`l ≥ 5` なら `v_p(1728) = 0`**。 -/
theorem valAdd_1728_eq_zero (p : HeightOneSpectrum (𝓞 L)) {l : ℕ} (hl : l.Prime)
    (hl5 : 5 ≤ l) (hnu : ¬ IsUnit ((l : primeSubring p))) :
    valAdd p (Units.mk0 (1728 : L) (by norm_num)) = 0 := by
  have hne : ((1728 : ℕ) : L) ≠ 0 := by norm_num
  have h := valAdd_natCast_eq_zero_of_five_le p hl hl5 hnu 6 3 (by norm_num) hne
  have hEq : Units.mk0 (((1728 : ℕ) : L)) hne = Units.mk0 (1728 : L) (by norm_num) := by
    refine Units.ext ?_
    push_cast
    rfl
  rwa [hEq] at h

/-! ## ★出典の紐付け(`.src`) -/

def isUnit_natCast_of_coprime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(非単元 l と互いに素な n は単元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_48_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p ∣ l・l ≥ 5 なら v_p(48) = 0)",
    sectionId := "genell-lemma-3-5" }

def valAdd_864_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p ∣ l・l ≥ 5 なら v_p(864) = 0)",
    sectionId := "genell-lemma-3-5" }

def valAdd_1728_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p ∣ l・l ≥ 5 なら v_p(1728) = 0)",
    sectionId := "genell-lemma-3-5" }

def valAdd_48_eq_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "valAdd_natCast_eq_zero_of_isUnit(第 1406、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.valAdd_natCast_eq_zero_of_isUnit") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1435）**——第 1425・1430・1431 が要る `p ∤ 6` を" ++
       "`p ∣ l`・`l ≥ 5` から作った。" ++
       "☆Bézout `a·l + b·48 = 1` を `𝒪_p` に落とすだけである" ++
       "——両方が非単元だと `1 ∈ 𝔪` になってしまう。" ++
       "★`l ≥ 5` は界面の側から取れる（`Check/GenEll/VeluQuotOKNeedsL5.lean`、第 1434）。") 17 ]

end ABC3.Found.GenEll
