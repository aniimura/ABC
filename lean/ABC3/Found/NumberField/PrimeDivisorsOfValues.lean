/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Tactic.Linarith

/-!
# Schur の補題 —— 多項式の値を割る素数は無限個ある

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何のためか

`Theorem 6.4, (iv)` は Tchebotarev の密度定理を **3 箇所**で使う。
そのうち (a)「`[L_i:ℚ]` は `deg(L_i, v_i)` の最大値」は、
**「Galois 拡大で完全分解する素数は無限個ある」**で足りる
(依存グラフ `ResearchPaper/frdi-decomposition.json` の鎖 `cheb` の `cheb-spl-infinite`)。

★★そしてそれは **Chebotarev 全体を使わずに出る** ——
`L = K(θ)`、`f` を `θ` の最小多項式とすると、`𝔭 ∤ disc(f)` について
「`𝔭` が完全分解 ⟺ `f` が `mod 𝔭` で根を持つ」(`L/K` が Galois のとき、Dedekind–Kummer)。
そして「`f` の値を割る素数は無限個」が本ファイルの **Schur の補題**である。

## ★証明(初等的、極限を使わない)

素数の集合 `S` が有限だとする。`c := f(n₀) ≠ 0` を取り(`f` は定数でないので根は有限個)、
`P := ∏_{p ∈ S} p`、`N := c·P` と置く。

* `x = n₀ + N·m` に対し `N ∣ f(x) − c`(`Polynomial.sub_dvd_eval_sub`)。
* よって `f(x) = c·u`、`u = 1 + P·j`。★**`u ≡ 1 (mod P)`** が要点である。
* `u ∈ {0, ±1}` は `f(x) ∈ {0, ±c}` と同値で、それは**有限個の `x`** でしか起きない
  (`f − C a` は非零多項式で根は有限個)。★**極限を使わずに済むのはここである。**
* そこで `|u| ≥ 2` なる `x` を取り、`u` の素因子 `q` を取ると
  `q ∣ f(x)` から `q ∈ S`、ゆえに `q ∣ P`、ゆえに `q ∣ u − P·j = 1` で矛盾。

★mathlib に同値の主張は見当たらなかった(`Schur` は Schur 補題・Schur 積などの別物、
`infinite_setOf_prime` 系は算術級数のもの)。
-/

namespace ABC3.Found.NumberField

open Polynomial

/-- ★`f` の値が `a` になる整数は有限個(`f` が定数でないとき)。 -/
theorem finite_eval_eq {f : ℤ[X]} (hf : 0 < f.natDegree) (a : ℤ) :
    {x : ℤ | f.eval x = a}.Finite := by
  have hne : f - C a ≠ 0 := by
    intro h
    have : f = C a := sub_eq_zero.mp h
    rw [this] at hf
    simp at hf
  refine (Polynomial.finite_setOf_isRoot hne).subset (fun x hx => ?_)
  show (f - C a).eval x = 0
  simp only [eval_sub, eval_C]
  rw [show f.eval x = a from hx]
  ring

/-- ★★★★**Schur の補題** —— 定数でない整数係数多項式の値を割る素数は**無限個**ある。

★★`Theorem 6.4, (iv)` が Tchebotarev を使う 3 箇所のうち (a) を、
**密度を使わずに**出すための鍵である(鎖 `cheb` の `cheb-spl-infinite`)。 -/
theorem infinite_primes_dvd_eval (f : ℤ[X]) (hf : 0 < f.natDegree) :
    {p : ℕ | p.Prime ∧ ∃ n : ℤ, (p : ℤ) ∣ f.eval n}.Infinite := by
  classical
  by_contra hfin
  rw [Set.not_infinite] at hfin
  -- ★1. 根でない `n₀` を取る(根は有限個だから)
  obtain ⟨n₀, hn₀⟩ : ∃ n₀ : ℤ, f.eval n₀ ≠ 0 := by
    by_contra hall
    push_neg at hall
    exact (Set.infinite_univ (α := ℤ)).mono (Set.subset_univ _)
      |>.elim (by simpa [hall] using (finite_eval_eq hf 0))
  set c : ℤ := f.eval n₀ with hc
  -- ★2. `P = ∏_{p ∈ S} p`、`N = c·P`
  set P : ℤ := ∏ p ∈ hfin.toFinset, (p : ℤ) with hP
  have hPpos : 0 < P := by
    refine Finset.prod_pos (fun p hp => ?_)
    have : p ∈ {p : ℕ | p.Prime ∧ ∃ n : ℤ, (p : ℤ) ∣ f.eval n} := hfin.mem_toFinset.mp hp
    exact_mod_cast this.1.pos
  set N : ℤ := c * P with hN
  have hNne : N ≠ 0 := mul_ne_zero hn₀ hPpos.ne'
  -- ★3. `f(x) ∈ {0, ±c}` となる `x` は有限個、`n₀ + N·ℤ` は無限個
  have hbad : ({x : ℤ | f.eval x = c} ∪ {x : ℤ | f.eval x = -c} ∪
      {x : ℤ | f.eval x = 0}).Finite :=
    ((finite_eval_eq hf c).union (finite_eval_eq hf (-c))).union (finite_eval_eq hf 0)
  have hAinf : (Set.range (fun m : ℤ => n₀ + N * m)).Infinite := by
    refine Set.infinite_range_of_injective (fun a b hab => ?_)
    have : N * a = N * b := by simpa using hab
    exact mul_left_cancel₀ hNne this
  obtain ⟨x, hxA, hxB⟩ := (hAinf.sdiff hbad).nonempty
  obtain ⟨m, rfl⟩ := hxA
  simp only [Set.mem_union, Set.mem_setOf_eq, not_or] at hxB
  obtain ⟨⟨hv1, hv2⟩, hv0⟩ := hxB
  set x : ℤ := n₀ + N * m with hx
  -- ★4. `f(x) = c·u`、`u = 1 + P·j`
  have hdvd : N ∣ f.eval x - c := by
    have h1 : (x - n₀) ∣ f.eval x - f.eval n₀ := Polynomial.sub_dvd_eval_sub x n₀ f
    have h2 : x - n₀ = N * m := by rw [hx]; ring
    rw [h2] at h1
    exact dvd_trans (Dvd.intro m rfl) h1
  obtain ⟨j, hj⟩ := hdvd
  set u : ℤ := 1 + P * j with hu
  have hval : f.eval x = c * u := by
    rw [hN] at hj
    rw [hu]
    linarith [hj]
  -- ★5. `u` は 0 でも ±1 でもない
  have hu0 : u ≠ 0 := fun h => hv0 (by rw [hval, h, mul_zero])
  have hu1 : u ≠ 1 := fun h => hv1 (by rw [hval, h, mul_one])
  have hum1 : u ≠ -1 := fun h => hv2 (by rw [hval, h]; ring)
  have hnat : u.natAbs ≠ 1 := by
    intro h
    rcases Int.natAbs_eq u with he | he
    · exact hu1 (by rw [he, h]; rfl)
    · exact hum1 (by rw [he, h]; rfl)
  -- ★6. `u` の素因子 `q` は `S` に入り、かつ `q ∣ 1` になる
  obtain ⟨q, hqp, hqu0⟩ := Nat.exists_prime_and_dvd hnat
  have hqu : (q : ℤ) ∣ u := Int.dvd_natAbs.mp (Int.natCast_dvd_natCast.mpr hqu0)
  have hqv : (q : ℤ) ∣ f.eval x := hval ▸ Dvd.dvd.mul_left hqu c
  have hqS : q ∈ {p : ℕ | p.Prime ∧ ∃ n : ℤ, (p : ℤ) ∣ f.eval n} := ⟨hqp, x, hqv⟩
  have hqP : (q : ℤ) ∣ P :=
    Finset.dvd_prod_of_mem (fun p : ℕ => (p : ℤ)) (hfin.mem_toFinset.mpr hqS)
  have hq1 : (q : ℤ) ∣ 1 := by
    have h1 : (q : ℤ) ∣ P * j := hqP.mul_right j
    have h2 : (q : ℤ) ∣ u - P * j := hqu.sub h1
    rw [hu] at h2
    simpa using h2
  have hq1' : q ∣ 1 := by exact_mod_cast hq1
  exact hqp.one_lt.ne' (Nat.eq_one_of_dvd_one hq1')

/-! ### ★出典の紐付け -/

/-- ★locator —— `Theorem 6.4, (iv)` が使う Tchebotarev の (a) を密度なしで出す鍵。 -/
def infinite_primes_dvd_eval.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 多項式の値を割る素数は無限個(Schur)",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NumberField
