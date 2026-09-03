/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.KernelNormVal
import ABC3.Meta.Claim

/-!
# 第 1395 ブロック —— **深い点は倍化で深いまま**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`hdbl` を落とす

第 1394 の `three_dvd_valAdd_negYdiff` は「`2P` の `x` も整である」（`hdbl`）を
仮定していた。★`p ∤ l` の良い素点なら第 1073 がそれをくれるが、
**`p ∣ l` ではくれない**——核の点が形式群に入りうるからである。

☆本ブロックは `hdbl` を**位数の奇数性で置き換える**。要は 2 段:

| 段 | 内容 |
|---|---|
| 1 | **深い点は倍化で深いまま**（深さも変わらない） |
| 2 | 奇素数位数なら `2^{l−1}·P = P`（Fermat の小定理） |

★★★段 1 は付値の帳簿だけで出る:

    v(Ψ₂Sq(x)) = 2 v(t) = −6m,  v(Φ₂(x)) = v(x⁴) = −8m
    ⟹ v(x(2P)) = −8m + 6m = −2m = v(x)

（`v_p(2) = 0` を使う。`Φ₂` の主項 `x⁴` が他の項より厳密に深い。）

★★★段 2 と合わせると:「`x_P` が整なのに `v(t_P) > 0`」を仮定すると
`2P` が深くなり（第 1394）、段 1 で `2^j·(2P)` がすべて深く、
`j = l−2` で `2^{l−1}·P = P` が深いことになって整性に反する。

☆**したがって `v(t_P) = 0`**——`hdbl` は要らない。

## ★★★★★★★★これで `3 ∣ v_p(N)` は `p ∣ l` でも使える

`three_dvd_valAdd_negYdiff_of_addOrderOf_prime` は
「良い素点（`v_p(Δ) = 0`）・`v_p(2) = 0`・位数が奇素数」だけを仮定する。
★これが第 1393 の `dvd_minDeltaExp_of_disc_pow_eq` にそのまま渡る形である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial IsDedekindDomain NumberField ABC3.Meta
open WeierstrassCurve.Affine

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★深い点は 2-捩れでない -/

/-- ☆**深い点では `t = y − negY(x,y) ≠ 0`**——`v_p(2) = 0` のとき。

★`2y` が `a₁x + a₃` より厳密に深いので和は零にならない。 -/
theorem negYdiff_ne_zero_of_deep (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {x y : L} (hx : x ≠ 0) (hy : y ≠ 0) (heq : W.toAffine.Equation x y)
    (hneg : valAdd p (Units.mk0 x hx) < 0) :
    y ≠ W.toAffine.negY x y := by
  obtain ⟨h1m, h2m, h3m, h4m, h6m⟩ := mem_primeSubring_of_isIntegral p W
  obtain ⟨m, hm0, hxm, hym⟩ := exists_depth_of_valAdd_x_neg p W hx hy heq hneg
  have h2yne : (2 : L) * y ≠ 0 := mul_ne_zero two_ne_zero hy
  have h2yval : valAdd p (Units.mk0 ((2 : L) * y) h2yne) = valAdd p (Units.mk0 y hy) := by
    rw [valAdd_mulL p two_ne_zero hy h2yne, h2]; ring
  have hrest : ValAtLeast p (-2 * (m : ℤ)) (W.a₁ * x + W.a₃) := by
    refine valAtLeast_add ?_ ?_
    · have hxA : ValAtLeast p (-2 * (m : ℤ)) x := by
        have := valAtLeast_unit p (Units.mk0 x hx)
        rw [hxm] at this; exact this
      exact valAtLeast_mono (by omega) (valAtLeast_mul (valAtLeast_of_mem h1m) hxA)
    · exact valAtLeast_mono (by omega) (valAtLeast_of_mem h3m)
  have hlt : valAdd p (Units.mk0 ((2 : L) * y) h2yne) < -2 * (m : ℤ) := by
    rw [h2yval, hym]; omega
  have hsumne : (2 : L) * y + (W.a₁ * x + W.a₃) ≠ 0 :=
    add_ne_zero_of_valAdd_lt h2yne hrest hlt
  intro hc
  apply hsumne
  rw [WeierstrassCurve.Affine.negY] at hc
  linear_combination hc

/-! ## ★★★★★★★★★★★★★★★★段 1 —— 倍化で深さは変わらない -/

/-- ★★★★★★★★★★★★★★★★
**深い点は倍化しても同じ深さ**——★**`v_p(2) = 0`**（第 1395）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v(Ψ₂Sq(x)) = −6m`・`v(Φ₂(x)) = −8m` なので `v(x(2P)) = −2m = v(x)` である。

★★★これが形式群の「`[2]` は母数を `2` 倍する」の付値版である。 -/
theorem valAdd_x_two_smul_of_deep (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {x y : L} (h : W.toAffine.Nonsingular x y) (hx : x ≠ 0)
    (hneg : valAdd p (Units.mk0 x hx) < 0) :
    ∃ (x₂ y₂ : L) (h₂ : W.toAffine.Nonsingular x₂ y₂) (hx₂ : x₂ ≠ 0),
      (2 : ℕ) • (Point.some x y h) = Point.some x₂ y₂ h₂ ∧
      valAdd p (Units.mk0 x₂ hx₂) = valAdd p (Units.mk0 x hx) := by
  obtain ⟨hb2, hb4, hb6, hb8⟩ := mem_primeSubring_b p W
  have hy : y ≠ 0 := y_ne_zero_of_valAdd_x_neg p W hx h.left hneg
  obtain ⟨m, hm0, hxm, hym⟩ := exists_depth_of_valAdd_x_neg p W hx hy h.left hneg
  have hty : y ≠ W.toAffine.negY x y := negYdiff_ne_zero_of_deep p W h2 hx hy h.left hneg
  have htne : y - W.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hty
  have htval : valAdd p (Units.mk0 (y - W.toAffine.negY x y) htne) = -3 * (m : ℤ) := by
    rw [valAdd_negYdiff_eq_of_deep p W h2 hx hy h.left hneg hty]; exact hym
  have hsq : W.Ψ₂Sq.eval x = (y - W.toAffine.negY x y) ^ 2 := psi2Sq_eval W h.left
  have hsqne : (y - W.toAffine.negY x y) ^ 2 ≠ 0 := pow_ne_zero 2 htne
  have hpsine : W.Ψ₂Sq.eval x ≠ 0 := by rw [hsq]; exact hsqne
  have hpsival : valAdd p (Units.mk0 (W.Ψ₂Sq.eval x) hpsine) = -6 * (m : ℤ) := by
    rw [valAdd_congr p hpsine hsqne hsq, valAdd_powL p htne 2 hsqne, htval]; push_cast; ring
  have hxA : ValAtLeast p (-2 * (m : ℤ)) x := by
    have := valAtLeast_unit p (Units.mk0 x hx); rw [hxm] at this; exact this
  have hx4ne : x ^ 4 ≠ 0 := pow_ne_zero 4 hx
  have hx4val : valAdd p (Units.mk0 (x ^ 4) hx4ne) = -8 * (m : ℤ) := by
    rw [valAdd_powL p hx 4 hx4ne, hxm]; push_cast; ring
  have h2mem : ((2 : L)) ∈ primeSubring p := by simp
  have hrest : ValAtLeast p (-4 * (m : ℤ))
      (-(W.b₄ * x ^ 2) + -(2 * W.b₆ * x) + -W.b₈) := by
    have hx2A : ValAtLeast p (-4 * (m : ℤ)) (x ^ 2) := by
      have hmm := valAtLeast_mul hxA hxA
      have hsq2 : x * x = x ^ 2 := by ring
      rw [hsq2] at hmm
      exact valAtLeast_mono (by omega) hmm
    refine valAtLeast_add (valAtLeast_add ?_ ?_) ?_
    · exact valAtLeast_neg
        (valAtLeast_mono (by omega) (valAtLeast_mul (valAtLeast_of_mem hb4) hx2A))
    · refine valAtLeast_neg ?_
      have h2b6 : ValAtLeast p 0 ((2 : L) * W.b₆) := valAtLeast_of_mem (mul_mem h2mem hb6)
      exact valAtLeast_mono (by omega) (valAtLeast_mul h2b6 hxA)
    · exact valAtLeast_neg (valAtLeast_mono (by omega) (valAtLeast_of_mem hb8))
  have hsplitP : (W.Φ 2).eval x
      = x ^ 4 + (-(W.b₄ * x ^ 2) + -(2 * W.b₆ * x) + -W.b₈) := by
    rw [Phi2_eval_expand]; ring
  have hltP : valAdd p (Units.mk0 (x ^ 4) hx4ne) < -4 * (m : ℤ) := by rw [hx4val]; omega
  have hPsum : x ^ 4 + (-(W.b₄ * x ^ 2) + -(2 * W.b₆ * x) + -W.b₈) ≠ 0 :=
    add_ne_zero_of_valAdd_lt hx4ne hrest hltP
  have hΦne : (W.Φ 2).eval x ≠ 0 := by rw [hsplitP]; exact hPsum
  have hΦval : valAdd p (Units.mk0 ((W.Φ 2).eval x) hΦne) = -8 * (m : ℤ) := by
    rw [valAdd_congr p hΦne hPsum hsplitP,
      valAdd_add_eq_of_lt hx4ne hPsum hrest hltP, hx4val]
  have h2ne : (W.ΨSq ((2 : ℕ) : ℤ)).eval x ≠ 0 := by
    rw [show ((2 : ℕ) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two]; exact hpsine
  obtain ⟨x₂, y₂, h₂, hP2, hx2eq, -⟩ := mulOK_two W h h2ne
  rw [show ((2 : ℕ) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two] at hx2eq
  have hx₂ne : x₂ ≠ 0 := by
    intro hc; apply hΦne; rw [← hx2eq, hc, zero_mul]
  refine ⟨x₂, y₂, h₂, hx₂ne, hP2, ?_⟩
  have hprod : x₂ * W.Ψ₂Sq.eval x ≠ 0 := by rw [hx2eq]; exact hΦne
  have hval := valAdd_mulL p hx₂ne hpsine hprod
  rw [valAdd_congr p hprod hΦne hx2eq, hΦval, hpsival] at hval
  rw [hxm]
  omega

/-- ★★★★★★★★**`2^j` 倍しても深いまま**——★（第 1395）。 -/
theorem exists_deep_two_pow_smul (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {x y : L} (h : W.toAffine.Nonsingular x y) (hx : x ≠ 0)
    (hneg : valAdd p (Units.mk0 x hx) < 0) (j : ℕ) :
    ∃ (x' y' : L) (h' : W.toAffine.Nonsingular x' y') (hx' : x' ≠ 0),
      (2 ^ j : ℕ) • (Point.some x y h) = Point.some x' y' h' ∧
      valAdd p (Units.mk0 x' hx') < 0 := by
  induction j with
  | zero => exact ⟨x, y, h, hx, by simp, hneg⟩
  | succ j ih =>
      obtain ⟨x', y', h', hx', heq', hneg'⟩ := ih
      obtain ⟨x₂, y₂, h₂, hx₂, heq₂, hval₂⟩ :=
        valAdd_x_two_smul_of_deep p W h2 h' hx' hneg'
      refine ⟨x₂, y₂, h₂, hx₂, ?_, by omega⟩
      have hpow : (2 ^ (j + 1) : ℕ) = 2 ^ j * 2 := by ring
      rw [hpow, mul_nsmul (Point.some x y h) (2 ^ j) 2, heq', heq₂]

/-! ## ★★★★★★★★段 2 —— Fermat の小定理で自分に戻る -/

omit [NumberField L] in
/-- ☆**奇素数位数の点は 2-捩れでない**。 -/
theorem negYdiff_ne_zero_of_addOrderOf_prime (W : WeierstrassCurve L)
    {x y : L} (h : W.toAffine.Nonsingular x y) {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hQ : addOrderOf (Point.some x y h) = l) :
    y ≠ W.toAffine.negY x y := by
  intro hc
  have hz : (2 : ℕ) • (Point.some x y h) = 0 :=
    (two_smul_eq_zero_iff W h).2 ((psi2Sq_eval_eq_zero_iff W h.left).2 hc)
  have hdvd : addOrderOf (Point.some x y h) ∣ 2 := addOrderOf_dvd_of_nsmul_eq_zero hz
  rw [hQ] at hdvd
  rcases (Nat.Prime.eq_one_or_self_of_dvd Nat.prime_two l hdvd) with h1 | h1
  · exact hl.one_lt.ne' h1
  · exact hodd h1

omit [NumberField L] in
/-- ★★★★★★**奇素数位数の点は `2^{l−1}` 倍で自分に戻る**——Fermat の小定理。 -/
theorem two_pow_smul_self (W : WeierstrassCurve L)
    {x y : L} (h : W.toAffine.Nonsingular x y) {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hQ : addOrderOf (Point.some x y h) = l) :
    ((2 : ℕ) ^ (l - 1)) • (Point.some x y h) = Point.some x y h := by
  have hcop : Nat.Coprime 2 l := (Nat.coprime_primes Nat.prime_two hl).2 (Ne.symm hodd)
  have hf := Nat.ModEq.pow_totient hcop
  rw [Nat.totient_prime hl] at hf
  have h1le : 1 ≤ (2 : ℕ) ^ (l - 1) := Nat.one_le_two_pow
  have hdvd : l ∣ (2 : ℕ) ^ (l - 1) - 1 := (Nat.modEq_iff_dvd' h1le).1 hf.symm
  obtain ⟨k, hk⟩ := hdvd
  have hsplit : (2 : ℕ) ^ (l - 1) = 1 + l * k := by omega
  rw [hsplit, add_nsmul, one_nsmul, mul_nsmul (Point.some x y h) l k,
    ← hQ, addOrderOf_nsmul_eq_zero, smul_zero, add_zero]

/-! ## ★★★★★★★★★★★★★★★★★★★★結論——`hdbl` の要らない形 -/

/-- ★★★★★★★★★★★★★★★★
**整な奇素数位数の点では `v_p(t_P) = 0`**——★（第 1395）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v_p(t_P) > 0` なら `2P` が深くなり（第 1394）、
深さは倍化で保たれるので `2^{l−1}·P = P` も深い——`x_P` の整性に反する。

★★★これで第 1394 の `hdbl`（`2P` の整性）は要らなくなる。 -/
theorem valAdd_negYdiff_eq_zero_of_addOrderOf_prime
    (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] (hΔ : W.Δ ≠ 0)
    (hΔ0 : valAdd p (Units.mk0 W.Δ hΔ) = 0)
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {x y : L} (h : W.toAffine.Nonsingular x y) (hx : x ∈ primeSubring p)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hQ : addOrderOf (Point.some x y h) = l)
    (hty : y ≠ W.toAffine.negY x y) :
    valAdd p (Units.mk0 (y - W.toAffine.negY x y) (sub_ne_zero.mpr hty)) = 0 := by
  have hy := mem_primeSubring_y_of_mem_x p W h.left hx
  have htne : y - W.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hty
  have hge : 0 ≤ valAdd p (Units.mk0 (y - W.toAffine.negY x y) htne) :=
    valAtLeast_of_mem (mem_primeSubring_negYdiff p W hx hy) htne
  by_contra hne
  have hpos : 0 < valAdd p (Units.mk0 (y - W.toAffine.negY x y) htne) := by omega
  obtain ⟨x₂, y₂, h₂, hx₂ne, hP2, hval⟩ :=
    valAdd_x_two_smul_eq p W hΔ hΔ0 h hx hty hpos
  have hneg2 : valAdd p (Units.mk0 x₂ hx₂ne) < 0 := by omega
  obtain ⟨x', y', h', hx', heq', hneg'⟩ :=
    exists_deep_two_pow_smul p W h2 h₂ hx₂ne hneg2 (l - 2)
  have hl3 : 3 ≤ l := by
    have := hl.two_le
    omega
  have hchain : ((2 : ℕ) ^ (l - 2)) • ((2 : ℕ) • Point.some x y h) = Point.some x y h := by
    rw [← mul_nsmul' (Point.some x y h) ((2 : ℕ) ^ (l - 2)) 2]
    have hpe : (2 : ℕ) ^ (l - 2) * 2 = 2 ^ (l - 1) := by
      rw [← pow_succ]; congr 1; omega
    rw [hpe, two_pow_smul_self W h hl hodd hQ]
  rw [hP2, heq', Point.some.injEq] at hchain
  obtain ⟨hxx, -⟩ := hchain
  have hx0 : x ≠ 0 := by rw [← hxx]; exact hx'
  have hxneg : valAdd p (Units.mk0 x hx0) < 0 := by
    rw [← valAdd_congr p hx' hx0 hxx]; exact hneg'
  have := valAtLeast_of_mem hx hx0
  omega

/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] 核のノルムの因子の付値は 3 の倍数**——★**良い素点・`v_p(2) = 0`**（第 1395）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆整な点では `0`、深い点では `−3m`。★**`p ∣ l` でも通る形**である。

★★★これが第 1393 の `dvd_minDeltaExp_of_disc_pow_eq` に渡る `3 ∣ v_p(N)` の
**点ごとの段の最終形**である。 -/
theorem three_dvd_valAdd_negYdiff_of_addOrderOf_prime
    (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] (hΔ : W.Δ ≠ 0)
    (hΔ0 : valAdd p (Units.mk0 W.Δ hΔ) = 0)
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {x y : L} (h : W.toAffine.Nonsingular x y) {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hQ : addOrderOf (Point.some x y h) = l)
    (hty : y ≠ W.toAffine.negY x y) :
    (3 : ℤ) ∣ valAdd p (Units.mk0 (y - W.toAffine.negY x y) (sub_ne_zero.mpr hty)) := by
  by_cases hx : x ∈ primeSubring p
  · rw [valAdd_negYdiff_eq_zero_of_addOrderOf_prime p W hΔ hΔ0 h2 h hx hl hodd hQ hty]
    exact dvd_zero 3
  · have hx0 : x ≠ 0 := fun h0 => hx (by rw [h0]; exact zero_mem _)
    have hneg : valAdd p (Units.mk0 x hx0) < 0 := by
      by_contra hge
      rw [not_lt] at hge
      exact hx ((mem_primeSubring_iff p x).2 ((valAdd_nonneg_iff p _).1 hge))
    have hy0 : y ≠ 0 := y_ne_zero_of_valAdd_x_neg p W hx0 h.left hneg
    rw [valAdd_negYdiff_eq_of_deep p W h2 hx0 hy0 h.left hneg hty]
    obtain ⟨m, hm0, hxm, hym⟩ := exists_depth_of_valAdd_x_neg p W hx0 hy0 h.left hneg
    exact ⟨-(m : ℤ), by omega⟩

/-! ## ★出典の紐付け(`.src`) -/

def negYdiff_ne_zero_of_deep.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い点では t = y − negY(x,y) ≠ 0。★v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def valAdd_x_two_smul_of_deep.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い点は倍化しても同じ深さ。★v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def valAdd_x_two_smul_of_deep.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1395）**——形式群の「`[2]` は母数を `2` 倍する」の付値版。" ++
       "☆`v(Ψ₂Sq(x)) = −6m`・`v(Φ₂(x)) = v(x⁴) = −8m` なので " ++
       "`v(x(2P)) = −2m = v(x)` である。`Φ₂` の主項 `x⁴` が他より厳密に深いのが要。") 17 ]

def exists_deep_two_pow_smul.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2^j 倍しても深いまま。★v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def negYdiff_ne_zero_of_addOrderOf_prime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(奇素数位数の点は 2-捩れでない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def two_pow_smul_self.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(奇素数位数の点は 2^{l−1} 倍で自分に戻る——Fermat の小定理。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_negYdiff_eq_zero_of_addOrderOf_prime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(整な奇素数位数の点では v_p(t_P) = 0。★良い素点・v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def three_dvd_valAdd_negYdiff_of_addOrderOf_prime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(核のノルムの因子の付値は 3 の倍数——hdbl 不要の形。★良い素点)",
    sectionId := "genell-lemma-3-5" }

def three_dvd_valAdd_negYdiff_of_addOrderOf_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "three_dvd_valAdd_negYdiff(第 1394、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.three_dvd_valAdd_negYdiff") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1395）**——第 1394 の `hdbl`（`2P` の整性）を" ++
       "位数の奇数性で置き換えた。☆`p ∤ l` では第 1073 が `hdbl` をくれるが、" ++
       "`p ∣ l` ではくれない——核の点が形式群に入りうるからである。" ++
       "★深さが倍化で保たれること（第 1395 段 1）と Fermat の小定理で置き換わる。") 17 ]

end ABC3.Found.GaloisRep
