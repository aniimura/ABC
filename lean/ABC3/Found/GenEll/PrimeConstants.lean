import ABC3.Found.GenEll.PrimeNumberTheorem

/-!
# [GenEll] §4 の定数 —— `Lemma 4.1` の仮説 (i)(ii) を具体化する(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.20-p.23。

原文 (GenEll p.21):
> which satisfy condition (ii) of Lemma 4.1 is entirely elementary. On the other hand, with regard to condition (i), the fact that

## ★★★★★なぜこれが要るのか

`Lemma 4.1` は『条件 (i)(ii) を満たす `ϵ, x_ϵ, C_ϵ` が**与えられたとき**』を主張する。
★`Corollary 4.3` / `Corollary 4.4` はそれを**適用する**側なので、
ここで初めて「そのような `ϵ, x_ϵ, C_ϵ` が**存在する**」ことが要る。

★★`ϵ := 1/6` に固定する——原文 p.23 が『"1 + 6ϵ" to be 2』と取るからである。

## ★★★★定数の順番が鍵である

★`C_ϵ < ϵ·x_ϵ` という制約があるので、**`C_ϵ` を先に固めてから `x_ϵ` を大きく取る**。

    θ(x)/x → 1 から X₀ を取る(x ≥ X₀ で 5/6·x < θ(x) < 5/4·x)
    C_ϵ := θ(X₀) + 1        ← 0 < x < X₀ では θ の単調性で押さえる
    x_ϵ := max X₀ (6(C_ϵ+1))  ← こう取れば C_ϵ < (1/6)·x_ϵ

★★`x_ϵ` を大きくしても (i) の第 2 式と (ii) は壊れない——どちらも
「`x ≥ x_ϵ` なら」の形だからである。

## ★本ファイルで取れるもの

| 定理 | 内容 |
|---|---|
| `theta_mono` | ★`θ` は非減少 |
| `exists_cond_i_ii` | ★★★★**(i)(ii) を満たす `x_ϵ, C_ϵ` が `ϵ = 1/6` で存在する** |
| `exists_finset_primes_sum_log_gt` | ★★★**対数和がいくらでも大きい有限素数集合がある** |

★最後の 1 本が、原文の『by enlarging `S` … we may always assume that `x_ϵ ≤ x_S`』の中身である。
-/

namespace ABC3.Found.GenEll

open Real Finset

/-! ## ★`θ` の単調性 -/

/-- ★**`θ` は非減少**である——被和関数が非負だから。 -/
theorem theta_mono {x y : ℝ} (h : x ≤ y) : Chebyshev.theta x ≤ Chebyshev.theta y := by
  have hf : ⌊x⌋₊ ≤ ⌊y⌋₊ := Nat.floor_le_floor h
  have hsub := theta_sub_theta x y hf
  have hnn : (0:ℝ) ≤ ∑ p ∈ Finset.Ioc ⌊x⌋₊ ⌊y⌋₊, logPrime p :=
    Finset.sum_nonneg (fun p _ => logPrime_nonneg p)
  linarith

theorem theta_nonneg (x : ℝ) : 0 ≤ Chebyshev.theta x := by
  rw [theta_eq_sum_logPrime]
  exact Finset.sum_nonneg (fun p _ => logPrime_nonneg p)

/-! ## ★★★★条件 (i)(ii) を満たす定数の存在 -/

/-- ★★★★**`ϵ = 1/6` に対して `Lemma 4.1` の条件 (i)(ii) を満たす `x_ϵ, C_ϵ` が存在する**。

原文 (GenEll p.21):
> which satisfy condition (ii) of Lemma 4.1 is entirely elementary. On the other hand, with regard to condition (i), the fact that

★(i) の 2 式はここでは**素数でない `x`** という制限なしに取れるので、
`lemma_4_1` の仮説にそのまま渡せる。 -/
theorem exists_cond_i_ii (M : ℕ) :
    ∃ xeps Ceps : ℝ, 0 < xeps ∧ 0 < Ceps ∧ Ceps < (1/6 : ℝ) * xeps
      ∧ (∀ x : ℝ, 0 < x → Chebyshev.theta x < 5 / 4 * x + Ceps)
      ∧ (∀ x : ℝ, 0 < x → xeps ≤ x → (1 - (1/6 : ℝ)) * x < Chebyshev.theta x)
      ∧ (∀ x : ℝ, 0 < x → xeps ≤ x → (M : ℝ) * Real.log x ≤ (1/6 : ℝ) * x) := by
  have htend := theta_div_tendsto_one
  have hev : ∀ᶠ x : ℝ in Filter.atTop, |Chebyshev.theta x / x - 1| < 1/6 := by
    have := Metric.tendsto_nhds.1 htend (1/6) (by norm_num)
    simpa [Real.dist_eq] using this
  obtain ⟨X₁, hX₁⟩ := Filter.eventually_atTop.mp hev
  set X₀ : ℝ := max X₁ 1 with hX₀def
  have hX₀pos : (0:ℝ) < X₀ := lt_of_lt_of_le zero_lt_one (le_max_right _ _)
  have hkey : ∀ x : ℝ, X₀ ≤ x → 0 < x →
      (5/6 : ℝ) * x < Chebyshev.theta x ∧ Chebyshev.theta x < 5/4 * x := by
    intro x hx hxpos
    have h := hX₁ x (le_trans (le_max_left _ _) hx)
    have h2 := abs_lt.1 h
    constructor
    · have hlt : (1 - 1/6 : ℝ) < Chebyshev.theta x / x := by linarith [h2.1]
      rw [lt_div_iff₀ hxpos] at hlt
      linarith
    · have hlt : Chebyshev.theta x / x < 1 + 1/6 := by linarith [h2.2]
      rw [div_lt_iff₀ hxpos] at hlt
      linarith
  set Ceps : ℝ := Chebyshev.theta X₀ + 1 with hCdef
  have hCpos : 0 < Ceps := by
    have := theta_nonneg X₀
    simp only [hCdef]
    linarith
  obtain ⟨y₀, hy₀pos, hy₀⟩ := exists_xeps_cond_ii M (1/6 : ℝ) (by norm_num)
  refine ⟨max (max X₀ y₀) (6 * (Ceps + 1)), Ceps, ?_, hCpos, ?_, ?_, ?_, ?_⟩
  · exact lt_of_lt_of_le hX₀pos (le_trans (le_max_left _ _) (le_max_left _ _))
  · have h6 : 6 * (Ceps + 1) ≤ max (max X₀ y₀) (6 * (Ceps + 1)) := le_max_right _ _
    linarith
  · intro x hxpos
    rcases le_or_gt X₀ x with hx | hx
    · have hup := (hkey x hx hxpos).2
      linarith
    · have hle : Chebyshev.theta x ≤ Chebyshev.theta X₀ := theta_mono hx.le
      have hq : (0:ℝ) < 5 / 4 * x := by linarith
      simp only [hCdef]
      linarith
  · intro x hxpos hx
    have hx0 : X₀ ≤ x := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hx
    have hlow := (hkey x hx0 hxpos).1
    linarith
  · intro x _ hx
    exact hy₀ x (le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hx)

/-! ## ★★★対数和の大きい有限素数集合 -/

/-- ★★★**いくらでも大きい対数和を持つ有限素数集合がある**。

原文 (GenEll p.23):
> Corollary 4.3, except that instead of applying condition (a) of Theorem 3.8, we

★原文は『by enlarging `S` [and possibly increasing the "C" of condition (c)],
we may always assume that `x_ϵ ≤ x_S`』と 1 文で済ませている。
★★その中身は `θ(x) → ∞` であり、`θ(x) > (5/6)·x` から直ちに出る。 -/
theorem exists_finset_primes_sum_log_gt (y : ℝ) :
    ∃ A : Finset ℕ, (∀ p ∈ A, p.Prime) ∧ y < ∑ p ∈ A, Real.log p := by
  obtain ⟨xeps, Ceps, hxpos, _, _, _, hlow, _⟩ := exists_cond_i_ii 1
  set x : ℝ := max xeps (max (2 * y + 1) 1) with hxdef
  have hxpos' : (0:ℝ) < x :=
    lt_of_lt_of_le zero_lt_one (le_trans (le_max_right _ _) (le_max_right _ _))
  have hx : xeps ≤ x := le_max_left _ _
  have h := hlow x hxpos' hx
  have hy : 2 * y + 1 ≤ x := le_trans (le_max_left _ _) (le_max_right _ _)
  refine ⟨(Finset.Ioc 0 ⌊x⌋₊).filter Nat.Prime, fun p hp => (Finset.mem_filter.1 hp).2, ?_⟩
  have hsum : ∑ p ∈ (Finset.Ioc 0 ⌊x⌋₊).filter Nat.Prime, Real.log p
      = Chebyshev.theta x := by
    rw [theta_eq_sum_logPrime, Finset.sum_filter]
    exact Finset.sum_congr rfl (fun p _ => by
      by_cases hp : p.Prime
      · simp [hp, logPrime_of_prime hp]
      · simp [hp, logPrime])
  rw [hsum]
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def exists_cond_i_ii.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Remark 4.1.1(条件 (i)(ii) を満たす定数の存在)",
    sectionId := "genell-rem-4-1-1" }

def exists_finset_primes_sum_log_gt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 23,
    item := "Corollary 4.3 の証明(by enlarging S … we may always assume that x_ϵ ≤ x_S)",
    sectionId := "genell-cor-4-3" }

end ABC3.Found.GenEll
