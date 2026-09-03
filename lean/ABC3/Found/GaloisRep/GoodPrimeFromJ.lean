/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MinDeltaShort
import ABC3.Found.GaloisRep.PointValuation
import ABC3.Found.GaloisRep.SemistableCriterion
import ABC3.Meta.Claim

/-!
# 第 1431 ブロック —— **`j` が整なら `v(Δ) = 12S` から半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——良い素点で `p ∣ l` の残りを 1 本にする

良い素点で `p ∣ l` の場合、判別式の恒等式（第 1402）と `3 ∣ v_p(N)`（第 1396）から

    `v_p(Δ(E′)) = 12·S`（`S = −v_p(N)/3`）

は出る。★足りないのは `v_p(c₄(E′)) ≥ 4S`・`v_p(c₆(E′)) ≥ 6S` で、
そこから第 1430（`minDeltaExp_eq_zero_of_short_integral`）が `minDeltaExp = 0` を出す。

☆本ブロックは**その 2 つが `jExp p E′ ≥ 0`（`j` の整性）から出る**ことを示す:

* `j = c₄³/Δ` なので `0 ≤ jExp = 3v(c₄) − 12S` ⟹ `v(c₄) ≥ 4S`
* `1728Δ = c₄³ − c₆²` なので `v(c₆²) ≥ min(3v(c₄), v(Δ)) ≥ 12S` ⟹ `v(c₆) ≥ 6S`
  （`v(1728) = 0`、すなわち `p ∤ 6`）

★★★これで良い素点の `p ∣ l` は **`jExp p E′ ≥ 0` 1 本**になる。
☆古典的にはモジュラー多項式 `Φ_l(j, j′) = 0` の単項性から出る（mathlib にはない）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★**`jExp ≥ 0` と `v(Δ) = 12S` から `v(c₄) ≥ 4S`**（第 1431）。

☆`j = Δ⁻¹·c₄³` なので `jExp = 3·v(c₄) − v(Δ)` である。 -/
theorem valAtLeast_c4_of_jExp_nonneg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] {S : ℤ}
    (hΔ : valAdd p (Units.mk0 W.Δ W.isUnit_Δ.ne_zero) = 12 * S)
    (hj : 0 ≤ jExp p W) :
    ValAtLeast p (4 * S) W.c₄ := by
  intro hc4ne
  have hΔne : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  have hjeq : W.j = W.Δ⁻¹ * W.c₄ ^ 3 := ABC3.Found.GenEll.j_eq_inv_Delta_mul W
  have hjne : W.j ≠ 0 := by
    rw [hjeq]
    exact mul_ne_zero (inv_ne_zero hΔne) (pow_ne_zero 3 hc4ne)
  have hunit : Units.mk0 W.j hjne
      = (Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔne)⁻¹ := by
    refine Units.ext ?_
    show W.j = ((Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔne)⁻¹ : Lˣ)
    push_cast
    rw [hjeq]
    show W.Δ⁻¹ * W.c₄ ^ 3 = W.c₄ ^ 3 * (W.Δ)⁻¹
    ring
  have hJ : jExp p W = valAdd p (Units.mk0 W.j hjne) := dif_neg hjne
  rw [hunit, valAdd_mul, valAdd_pow, valAdd_inv, hΔ] at hJ
  omega

set_option maxHeartbeats 1600000 in
/-- ★★★★**`jExp ≥ 0` と `v(Δ) = 12S`・`v(1728) = 0` から `v(c₆) ≥ 6S`**（第 1431）。

☆`1728Δ = c₄³ − c₆²` と `v(c₄) ≥ 4S` から `v(c₆²) ≥ 12S`。 -/
theorem valAtLeast_c6_of_jExp_nonneg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] {S : ℤ}
    (hΔ : valAdd p (Units.mk0 W.Δ W.isUnit_Δ.ne_zero) = 12 * S)
    (hj : 0 ≤ jExp p W)
    (h1728 : valAdd p (Units.mk0 (1728 : L) (by norm_num)) = 0) :
    ValAtLeast p (6 * S) W.c₆ := by
  have hΔne : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  have h4 : ValAtLeast p (4 * S) W.c₄ := valAtLeast_c4_of_jExp_nonneg p W hΔ hj
  have h43 : ValAtLeast p (12 * S) (W.c₄ ^ 3) := by
    have h := valAtLeast_mul (valAtLeast_mul h4 h4) h4
    have hEq : W.c₄ * W.c₄ * W.c₄ = W.c₄ ^ 3 := by ring
    rw [hEq] at h
    exact valAtLeast_mono (by omega) h
  have hΔv : ValAtLeast p (12 * S) W.Δ := by
    have h := valAtLeast_unit p (Units.mk0 W.Δ hΔne)
    rw [hΔ] at h
    exact h
  have h1728v : ValAtLeast p 0 ((1728 : L)) := by
    have h := valAtLeast_unit p (Units.mk0 (1728 : L) (by norm_num))
    rw [h1728] at h
    exact h
  have hprod : ValAtLeast p (12 * S) ((1728 : L) * W.Δ) := by
    have h := valAtLeast_mul h1728v hΔv
    simpa using h
  have hc6sq : ValAtLeast p (12 * S) (W.c₆ ^ 2) := by
    have hEq : W.c₆ ^ 2 = W.c₄ ^ 3 + (-((1728 : L) * W.Δ)) := by
      have h := W.c_relation
      linear_combination h
    rw [hEq]
    exact valAtLeast_add h43 (valAtLeast_neg hprod)
  intro hc6ne
  have hsq : valAdd p (Units.mk0 (W.c₆ ^ 2) (pow_ne_zero 2 hc6ne))
      = 2 * valAdd p (Units.mk0 W.c₆ hc6ne) := by
    have h : Units.mk0 (W.c₆ ^ 2) (pow_ne_zero 2 hc6ne) = (Units.mk0 W.c₆ hc6ne) ^ 2 :=
      Units.ext (by push_cast; rfl)
    rw [h, valAdd_pow]
    norm_num
  have hge := hc6sq (pow_ne_zero 2 hc6ne)
  omega

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★
**`jExp p W ≥ 0` と `v_p(Δ) = 12·v_p(n)`・`p ∤ 6` なら `minDeltaExp p W = 0`**（第 1431）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが良い素点で `p ∣ l` の残りを **`jExp p E′ ≥ 0` 1 本**にする形である。 -/
theorem minDeltaExp_eq_zero_of_jExp_nonneg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] {n : L} (hn : n ≠ 0)
    (hΔ : valAdd p (Units.mk0 W.Δ W.isUnit_Δ.ne_zero)
      = 12 * valAdd p (Units.mk0 n hn))
    (hj : 0 ≤ jExp p W)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (h1728 : valAdd p (Units.mk0 (1728 : L) (by norm_num)) = 0) :
    minDeltaExp p W = 0 := by
  have h4 := valAtLeast_c4_of_jExp_nonneg p W hΔ hj
  have h6 := valAtLeast_c6_of_jExp_nonneg p W hΔ hj h1728
  refine minDeltaExp_eq_zero_of_short_integral p W hn ?_ ?_ hΔ
  · refine div_mem_primeSubring_of_valAdd_le p (by norm_num) h48 hn 4 (fun hc0 => ?_)
    have := h4 hc0
    push_cast
    omega
  · refine div_mem_primeSubring_of_valAdd_le p (by norm_num) h864 hn 6 (fun hc0 => ?_)
    have := h6 hc0
    push_cast
    omega

/-- ★★★★★★★★★★★★**同じ仮定から半安定**（第 1431）。 -/
theorem semistableAt_of_jExp_nonneg_of_valAdd_Delta (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] {n : L} (hn : n ≠ 0)
    (hΔ : valAdd p (Units.mk0 W.Δ W.isUnit_Δ.ne_zero)
      = 12 * valAdd p (Units.mk0 n hn))
    (hj : 0 ≤ jExp p W)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (h1728 : valAdd p (Units.mk0 (1728 : L) (by norm_num)) = 0) :
    SemistableAt p W :=
  Or.inl (minDeltaExp_eq_zero_of_jExp_nonneg p W hn hΔ hj h48 h864 h1728)

/-! ## ★出典の紐付け(`.src`) -/

def valAtLeast_c4_of_jExp_nonneg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(jExp ≥ 0 と v(Δ) = 12S から v(c₄) ≥ 4S。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAtLeast_c6_of_jExp_nonneg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(jExp ≥ 0 と v(Δ) = 12S から v(c₆) ≥ 6S。★p ∤ 6)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_zero_of_jExp_nonneg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(jExp ≥ 0 と v(Δ) = 12v(n) なら minDeltaExp = 0)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_of_jExp_nonneg_of_valAdd_Delta.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(jExp ≥ 0 と v(Δ) = 12v(n) なら半安定。★p ∤ 6)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_zero_of_jExp_nonneg.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_zero_of_short_integral(第 1430、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_zero_of_short_integral") 1,
    .citation "[mathlib]" "WeierstrassCurve.c_relation(1728Δ = c₄³ − c₆²)"
      (.inMathlib "WeierstrassCurve.c_relation") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1431）**——良い素点で `p ∣ l` の残りを" ++
       "**`jExp p E′ ≥ 0` 1 本**にした。" ++
       "☆判別式の恒等式（第 1402）と `3 ∣ v_p(N)`（第 1396）から `v_p(Δ(E′)) = 12S` は出るので、" ++
       "`j` の整性さえ言えれば `v(c₄) ≥ 4S`・`v(c₆) ≥ 6S` が従い、" ++
       "第 1430 で `minDeltaExp = 0` になる。" ++
       "★`j(E′)` の整性は古典的にはモジュラー多項式 `Φ_l(j, j′) = 0` の単項性から出る" ++
       "（mathlib にはない、2026-09-02 確認）。") 17 ]

end ABC3.Found.GaloisRep
