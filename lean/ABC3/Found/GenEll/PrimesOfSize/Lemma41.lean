/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.Chebyshev

/-!
# PrimesOfSize —— `[GenEll] Lemma 4.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open Real Finset

/-! ## ★準備 —— `θ` の被和関数 -/

/-- `Chebyshev.theta` の被和関数を `Finset.filter` から外に出したもの。 -/
noncomputable def logPrime (p : ℕ) : ℝ := if p.Prime then Real.log p else 0

theorem logPrime_nonneg (p : ℕ) : 0 ≤ logPrime p := by
  unfold logPrime
  split_ifs with hp
  · exact Real.log_natCast_nonneg p
  · exact le_rfl

theorem logPrime_of_prime {p : ℕ} (hp : p.Prime) : logPrime p = Real.log p := if_pos hp

theorem theta_eq_sum_logPrime (x : ℝ) :
    Chebyshev.theta x = ∑ p ∈ Finset.Ioc 0 ⌊x⌋₊, logPrime p := by
  rw [Chebyshev.theta, Finset.sum_filter]
  rfl

/-- `θ(y) − θ(z) = Σ_{⌊z⌋ < p ≤ ⌊y⌋} logPrime p`。 -/
theorem theta_sub_theta (z y : ℝ) (hzy : ⌊z⌋₊ ≤ ⌊y⌋₊) :
    Chebyshev.theta y - Chebyshev.theta z
      = ∑ p ∈ Finset.Ioc ⌊z⌋₊ ⌊y⌋₊, logPrime p := by
  rw [theta_eq_sum_logPrime, theta_eq_sum_logPrime, sub_eq_iff_eq_add']
  exact (Finset.sum_Ioc_consecutive _ (Nat.zero_le _) hzy).symm

/-! ## ★原文の WLOG を消す 2 本の補題 -/

/-- ★**条件 (i) を `ℝ′_{>0}` から `ℝ_{≥0}` 全体へ延ばす。**

原文 (GenEll p.20):
> Lemma 4.1. (The Existence of Primes of Prescribed Size) Write

原文の条件 (i) は `x ∈ ℝ′_{>0}`(素数でない正実数)についてしか与えられていない。
★**`θ(x) = θ(⌊x⌋₊)` を使えば、`x` を同じ床を持つ非自然数へずらせる。**
上からの評価なので、ずらし幅を `0` に近づけて `<` を `≤` に落とす。 -/
theorem theta_le_of_cond_i {Ceps : ℝ} (_hCeps : 0 < Ceps)
    (hi1 : ∀ x : ℝ, 0 < x → (∀ p : ℕ, p.Prime → (p : ℝ) ≠ x) →
      Chebyshev.theta x < 5 / 4 * x + Ceps) :
    ∀ x : ℝ, 0 ≤ x → Chebyshev.theta x ≤ 5 / 4 * x + Ceps := by
  intro x hx
  by_contra hcon
  push_neg at hcon
  set d : ℝ := Chebyshev.theta x - (5 / 4 * x + Ceps) with hddef
  have hd : 0 < d := by simp only [hddef]; linarith
  set n : ℕ := ⌊x⌋₊ with hndef
  have hnx : (n : ℝ) ≤ x := Nat.floor_le hx
  have hxn : x < (n : ℝ) + 1 := Nat.lt_floor_add_one x
  -- ずらし幅
  set η : ℝ := min (((n : ℝ) + 1 - x) / 2) (d * 2 / 5) with hηdef
  have hη0 : 0 < η := lt_min (by linarith) (by positivity)
  have hηa : η ≤ ((n : ℝ) + 1 - x) / 2 := min_le_left _ _
  have hηb : η ≤ d * 2 / 5 := min_le_right _ _
  -- `x + η` は `n` と同じ床を持ち、自然数ではない
  have hfl : ⌊x + η⌋₊ = n := by
    rw [hndef, Nat.floor_eq_iff (by positivity)]
    constructor
    · linarith
    · push_cast; linarith
  have hnp : ∀ p : ℕ, p.Prime → (p : ℝ) ≠ x + η := by
    intro p _ hpe
    have h1 : (n : ℝ) < (p : ℝ) := by rw [hpe]; linarith
    have h2 : (p : ℝ) < (n : ℝ) + 1 := by rw [hpe]; linarith
    have h1' : n < p := by exact_mod_cast h1
    have h2' : p < n + 1 := by exact_mod_cast h2
    omega
  have hθ : Chebyshev.theta (x + η) = Chebyshev.theta x := by
    rw [Chebyshev.theta_eq_theta_coe_floor (x + η), hfl,
      ← Chebyshev.theta_eq_theta_coe_floor]
  have := hi1 (x + η) (by linarith) hnp
  rw [hθ] at this
  simp only [hddef] at hd ⊢
  nlinarith

/-- ★**条件 (ii) を `ℝ′_{>0}` から `ℝ_{>0}` 全体へ延ばす。**

下からの評価なので、`x` を同じ床を持つ**わずかに大きい**非自然数へずらすだけでよく、
極限を取る必要がない(`<` が保たれる)。 -/
theorem theta_gt_of_cond_ii {eps xeps : ℝ} (heps1 : eps < 1)
    (hi2 : ∀ x : ℝ, 0 < x → (∀ p : ℕ, p.Prime → (p : ℝ) ≠ x) → xeps ≤ x →
      (1 - eps) * x < Chebyshev.theta x) :
    ∀ x : ℝ, 0 < x → xeps ≤ x → (1 - eps) * x < Chebyshev.theta x := by
  intro x hx hxe
  set n : ℕ := ⌊x⌋₊ with hndef
  have hnx : (n : ℝ) ≤ x := Nat.floor_le hx.le
  have hxn : x < (n : ℝ) + 1 := Nat.lt_floor_add_one x
  set η : ℝ := ((n : ℝ) + 1 - x) / 2 with hηdef
  have hη0 : 0 < η := by simp only [hηdef]; linarith
  have hfl : ⌊x + η⌋₊ = n := by
    rw [hndef, Nat.floor_eq_iff (by positivity)]
    constructor
    · simp only [hηdef]; linarith
    · push_cast; simp only [hηdef]; linarith
  have hnp : ∀ p : ℕ, p.Prime → (p : ℝ) ≠ x + η := by
    intro p _ hpe
    have h1 : (n : ℝ) < (p : ℝ) := by rw [hpe]; linarith
    have h2 : (p : ℝ) < (n : ℝ) + 1 := by rw [hpe]; simp only [hηdef]; linarith
    have h1' : n < p := by exact_mod_cast h1
    have h2' : p < n + 1 := by exact_mod_cast h2
    omega
  have hθ : Chebyshev.theta (x + η) = Chebyshev.theta x := by
    rw [Chebyshev.theta_eq_theta_coe_floor (x + η), hfl,
      ← Chebyshev.theta_eq_theta_coe_floor]
  have := hi2 (x + η) (by linarith) hnp (by linarith)
  rw [hθ] at this
  nlinarith

/-! ## Lemma 4.1 本体 -/

/-- **[GenEll] Lemma 4.1**(The Existence of Primes of Prescribed Size)。

原文 (GenEll p.20):
> Lemma 4.1. (The Existence of Primes of Prescribed Size) Write

`ℝ′_{>0}` を `ℝ_{>0}` から素数を除いたもの、`θ(x) ≝ Σ_{p<x} log(p)` とし、
`M` 正整数、`ϵ, x_ϵ, C_ϵ ∈ ℝ_{>0}` が `0 < ϵ < 1/4`、`ϵ·x_ϵ > C_ϵ` と条件 (i)(ii) を
満たすとする。このとき任意の非負 `h` と `x_A > x_ϵ` なる任意の有限素数集合 `A` に対し、
**`M` 個の相異なる素数 `p_1,…,p_M` が存在して `p_j ∉ A` かつ
`h ≤ p_j ≤ (1+6ϵ)·x_A + 8h`**。

## ★原文の証明(p.21)の道筋をそのまま追った

結論が偽と仮定する。`δ ≝ 6ϵ`、`y_A ≝ (1+δ)·x_A + 8h` と置くと、
`h ≤ p ≤ y_A` なる素数はすべて——`M−1` 個の例外を除いて——`A` に属する。
すると `x_A ≥ −M·log(y_A) − θ((1+δ)h) + θ(y_A)` であり、これを

- `θ(y_A) > (1−ϵ)y_A`(条件 (i) の後半)
- `θ((1+δ)h) ≤ (5/4)(1+δ)h + C_ϵ`(条件 (i) の前半)
- `M·log(y_A) ≤ ϵ·y_A`(条件 (ii))

で評価すると `x_A > x_A` になり矛盾する。

★最後の矛盾は `ϵ·x_A > ϵ·x_ϵ > C_ϵ` と `5(1+δ)/4 < 4` から出る。

## ★★原文の WLOG は使っていない

原文は「`y_A, (1+δ)h ∈ ℝ′_{>0}` を仮定してよい」と書くが、
**我々は `θ` の側を `ℝ_{>0}` 全体へ延ばした**(`theta_le_of_cond_i` /
`theta_gt_of_cond_ii`)。★これにより `δ, h` を動かさずに済む。 -/
theorem lemma_4_1 (M : ℕ) (hM : 0 < M)
    (eps xeps Ceps : ℝ) (heps : 0 < eps) (hxeps : 0 < xeps) (hCeps : 0 < Ceps)
    (heps4 : eps < 1 / 4) (hxC : Ceps < eps * xeps)
    (hi1 : ∀ x : ℝ, 0 < x → (∀ p : ℕ, p.Prime → (p : ℝ) ≠ x) →
      Chebyshev.theta x < 5 / 4 * x + Ceps)
    (hi2 : ∀ x : ℝ, 0 < x → (∀ p : ℕ, p.Prime → (p : ℝ) ≠ x) → xeps ≤ x →
      (1 - eps) * x < Chebyshev.theta x)
    (hii : ∀ x : ℝ, 0 < x → xeps ≤ x → (M : ℝ) * Real.log x ≤ eps * x)
    (h : ℝ) (hh : 0 ≤ h) (A : Finset ℕ) (hA : ∀ p ∈ A, p.Prime)
    (hAx : xeps < ∑ p ∈ A, Real.log p) :
    ∃ P : Finset ℕ, P.card = M ∧ (∀ p ∈ P, p.Prime) ∧ (∀ p ∈ P, p ∉ A) ∧
      ∀ p ∈ P, h ≤ (p : ℝ) ∧
        (p : ℝ) ≤ (1 + 6 * eps) * (∑ q ∈ A, Real.log q) + 8 * h := by
  classical
  set xA : ℝ := ∑ p ∈ A, Real.log p with hxAdef
  set Y : ℝ := (1 + 6 * eps) * xA + 8 * h with hYdef
  set H : ℝ := (1 + 6 * eps) * h with hHdef
  have hxA0 : 0 < xA := lt_trans hxeps hAx
  have hY0 : 0 < Y := by simp only [hYdef]; nlinarith
  have hH0 : 0 ≤ H := by simp only [hHdef]; nlinarith
  have hHY : H ≤ Y := by simp only [hHdef, hYdef]; nlinarith
  have hYxe : xeps ≤ Y := by simp only [hYdef]; nlinarith
  -- ★候補集合 `S` —— `h ≤ p ≤ Y` なる `A` の外の素数
  set S : Finset ℕ :=
    (Finset.range (⌊Y⌋₊ + 1)).filter
      (fun p => p.Prime ∧ p ∉ A ∧ h ≤ (p : ℝ) ∧ (p : ℝ) ≤ Y) with hSdef
  -- ★★`Y ≥ 2` は**仮定から強制される**(置いた仮定ではない)。
  -- `Y < 2` なら `θ(Y) = 0` だが、延長した条件 (ii) は `(1−ϵ)Y < θ(Y)` を要求する。
  have hthY : (1 - eps) * Y < Chebyshev.theta Y :=
    theta_gt_of_cond_ii (by linarith) hi2 Y hY0 hYxe
  have hY2 : (2 : ℝ) ≤ Y := by
    by_contra hlt
    push_neg at hlt
    rw [Chebyshev.theta_eq_zero_of_lt_two hlt] at hthY
    nlinarith
  have hlogY0 : 0 ≤ Real.log Y := Real.log_nonneg (by linarith)
  by_contra hcon
  -- ★`S.card < M` —— そうでなければ `M` 元部分集合が取れて結論が成り立つ
  have hScard : S.card < M := by
    by_contra hge
    push_neg at hge
    obtain ⟨P, hPS, hPcard⟩ := Finset.exists_subset_card_eq hge
    refine hcon ⟨P, hPcard, ?_, ?_, ?_⟩
    · intro p hp
      have := hPS hp
      simp only [hSdef, Finset.mem_filter] at this
      exact this.2.1
    · intro p hp
      have := hPS hp
      simp only [hSdef, Finset.mem_filter] at this
      exact this.2.2.1
    · intro p hp
      have := hPS hp
      simp only [hSdef, Finset.mem_filter] at this
      exact ⟨this.2.2.2.1, this.2.2.2.2⟩
  -- ★`θ(Y) − θ(H)` を `Ioc ⌊H⌋ ⌊Y⌋` 上の和として書く
  have hfloor : ⌊H⌋₊ ≤ ⌊Y⌋₊ := Nat.floor_le_floor hHY
  have hsplit := theta_sub_theta H Y hfloor
  -- ★`Ioc ⌊H⌋ ⌊Y⌋` の素数は `A` か `S` に入る
  have hkey : ∑ p ∈ Finset.Ioc ⌊H⌋₊ ⌊Y⌋₊, logPrime p
      ≤ xA + (M - 1 : ℕ) * Real.log Y := by
    set T : Finset ℕ := Finset.Ioc ⌊H⌋₊ ⌊Y⌋₊ with hTdef
    have hTA : ∑ p ∈ T ∩ A, logPrime p ≤ xA := by
      calc ∑ p ∈ T ∩ A, logPrime p ≤ ∑ p ∈ A, logPrime p :=
            Finset.sum_le_sum_of_subset_of_nonneg Finset.inter_subset_right
              (fun p _ _ => logPrime_nonneg p)
        _ = xA := by
            rw [hxAdef]
            exact Finset.sum_congr rfl fun p hp => logPrime_of_prime (hA p hp)
    have hTS : ∑ p ∈ T \ A, logPrime p ≤ (M - 1 : ℕ) * Real.log Y := by
      have hsub : (T \ A).filter (fun p => p.Prime) ⊆ S := by
        intro p hp
        simp only [Finset.mem_filter, Finset.mem_sdiff, hTdef, Finset.mem_Ioc] at hp
        obtain ⟨⟨⟨hlo, hhi⟩, hnA⟩, hpp⟩ := hp
        simp only [hSdef, Finset.mem_filter, Finset.mem_range]
        have hpY : (p : ℝ) ≤ Y := by
          calc (p : ℝ) ≤ (⌊Y⌋₊ : ℝ) := by exact_mod_cast hhi
            _ ≤ Y := Nat.floor_le hY0.le
        have hpH : H < (p : ℝ) := by
          have : (⌊H⌋₊ : ℝ) + 1 ≤ (p : ℝ) := by exact_mod_cast hlo
          have hHfl : H < (⌊H⌋₊ : ℝ) + 1 := Nat.lt_floor_add_one H
          linarith
        have hph : h ≤ (p : ℝ) := by
          have : h ≤ H := by simp only [hHdef]; nlinarith
          linarith
        exact ⟨by omega, hpp, hnA, hph, hpY⟩
      have hzero : ∀ p ∈ T \ A, ¬ p.Prime → logPrime p = 0 := by
        intro p _ hnp; unfold logPrime; exact if_neg hnp
      calc ∑ p ∈ T \ A, logPrime p
          = ∑ p ∈ (T \ A).filter (fun p => p.Prime), logPrime p := by
            rw [Finset.sum_filter_of_ne]
            intro p hp hne
            by_contra hnp
            exact hne (hzero p hp hnp)
        _ ≤ ∑ p ∈ S, logPrime p :=
            Finset.sum_le_sum_of_subset_of_nonneg hsub (fun p _ _ => logPrime_nonneg p)
        _ ≤ ∑ _p ∈ S, Real.log Y := by
            refine Finset.sum_le_sum fun p hp => ?_
            simp only [hSdef, Finset.mem_filter] at hp
            rw [logPrime_of_prime hp.2.1]
            exact Real.log_le_log (by exact_mod_cast hp.2.1.pos) hp.2.2.2.2
        _ = (S.card : ℝ) * Real.log Y := by rw [Finset.sum_const, nsmul_eq_mul]
        _ ≤ (M - 1 : ℕ) * Real.log Y := by
            have hcard : (S.card : ℝ) ≤ ((M - 1 : ℕ) : ℝ) := by
              have : S.card ≤ M - 1 := by omega
              exact_mod_cast this
            exact mul_le_mul_of_nonneg_right hcard hlogY0
    calc ∑ p ∈ T, logPrime p
        = ∑ p ∈ T ∩ A, logPrime p + ∑ p ∈ T \ A, logPrime p := by
          rw [Finset.sum_inter_add_sum_sdiff]
      _ ≤ xA + (M - 1 : ℕ) * Real.log Y := add_le_add hTA hTS
  -- ★残る 2 つの評価
  have hthH : Chebyshev.theta H ≤ 5 / 4 * H + Ceps :=
    theta_le_of_cond_i hCeps hi1 H hH0
  have hMlog : (M : ℝ) * Real.log Y ≤ eps * Y := hii Y hY0 hYxe
  have hMsub : ((M - 1 : ℕ) : ℝ) ≤ (M : ℝ) := by
    have : (M - 1 : ℕ) ≤ M := Nat.sub_le _ _
    exact_mod_cast this
  -- ★★原文 p.21 の不等式の連鎖
  have hchain : xA ≥ Chebyshev.theta Y - Chebyshev.theta H - (M : ℝ) * Real.log Y := by
    have h1 : Chebyshev.theta Y - Chebyshev.theta H ≤ xA + ((M - 1 : ℕ) : ℝ) * Real.log Y := by
      rw [hsplit]; exact hkey
    have h2 : ((M - 1 : ℕ) : ℝ) * Real.log Y ≤ (M : ℝ) * Real.log Y :=
      mul_le_mul_of_nonneg_right hMsub hlogY0
    linarith
  -- `x_A > (1−2ϵ)·y_A − (5/4)·(1+δ)h − C_ϵ`
  have hstep : xA > (1 - 2 * eps) * Y - 5 / 4 * H - Ceps := by linarith
  rw [hYdef, hHdef] at hstep
  -- ★`(1−2ϵ)(1+6ϵ) ≥ 1 + ϵ`(`ϵ < 1/4` が効く唯一の箇所)
  have hcoef1 : (1 - 2 * eps) * (1 + 6 * eps) ≥ 1 + eps := by nlinarith
  -- ★`8(1−2ϵ) − (5/4)(1+6ϵ) ≥ 0` —— 原文の「`5(1+δ)/4 < 4`」に対応
  have hcoef2 : 8 * (1 - 2 * eps) - 5 / 4 * (1 + 6 * eps) ≥ 0 := by nlinarith
  have hfinal : xA > (1 + eps) * xA - Ceps := by nlinarith
  -- ★`ϵ·x_A > ϵ·x_ϵ > C_ϵ` との矛盾
  nlinarith [mul_lt_mul_of_pos_left hAx heps]

def lemma_4_1.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 20, item := "Lemma 4.1",
    sectionId := "genell-lemma-4-1" }

def theta_le_of_cond_i.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 20, item := "Lemma 4.1",
    sectionId := "genell-lemma-4-1" }

def theta_gt_of_cond_ii.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 20, item := "Lemma 4.1",
    sectionId := "genell-lemma-4-1" }

end ABC3.Found.GenEll
