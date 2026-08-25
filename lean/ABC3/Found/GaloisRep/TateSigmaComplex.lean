import ABC3.Found.GaloisRep.TateDoubleSum

/-!
# Galois (G6) 第 229 ブロック —— **★★★★★★★`s₁` の q 展開が ℂ の上で取れた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★§9-540 で書き留めた穴を塞ぐ

    ∑_{m≥1} f(qᵐ) = ∑_{N≥1} σ₁(N) qᴺ        (ℂ の上、`‖q‖ < 1`)

★形式側では第 111 ブロック(`tateXtail_one`)で取れているが、`evalAdic` の議論は
`I` 進完備環でしか使えないので、**ℂ の上では別に要る**。

## ★★★★繊維で括り直す

第 228 で `∑_{m≥1}f(qᵐ) = ∑_{(m,d)} (d+1)q^{(m+1)(d+1)}` まで開いた。
これを `N = (m+1)(d+1) − 1` の**繊維**で括り直す:

| 段 | 道具 |
|---|---|
| `ℕ×ℕ ≃ Σ N, 繊維 N` | `Equiv.sigmaFiberEquiv` |
| 二重和を繰り返し和に | `Summable.tsum_sigma` |
| 繊維 ≃ 約数対の集合 | ★★★★本ブロックの `fiberEquiv` |
| 有限型の `tsum` を `Finset` の和に | `Finset.tsum_subtype` |
| `∑_{(a,b):ab=n} b = σ₁(n)` | `Nat.sum_divisorsAntidiagonal` + `Nat.sum_div_divisors` |

★★`fiberEquiv N : {p // (p.1+1)(p.2+1)−1 = N} ≃ ↥((N+1).divisorsAntidiagonal)` を
**非依存な型の間の同値**として作るのが要点である。`Σ` の中で直接作ろうとすると
第二成分が依存型になり、`Sigma.mk.injEq` の書き換えで motive が壊れる。
★`Equiv.sigmaFiberEquiv` に繊維分解は任せ、**繊維ごとの同値だけを手で作る**。

## ★★★`ℕ` の引き算に注意

繊維の側は `(m+1, d+1)`、約数対の側は `(a, b)` で `a, b ≥ 1`。
行き来は `+1` と `−1` だが、`ℕ` の切り捨て引き算なので `1 ≤ a` を出してから
`omega` に渡す(`pos_of_mem_divAntidiag`)。
★★`omega` は `(m+1, d+1).1 * (m+1, d+1).2` と `(m+1)*(d+1)` を**別の原子**として扱う。
`show (m+1)*(d+1) = …` で形を揃えてから渡すこと。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pos_of_mem_divAntidiag` | ★約数対の成分は `1` 以上 |
| `fiberEquiv` | ★★★★**繊維 ≃ 約数対の集合** |
| `sum_divAntidiag_snd` | ★★`∑_{(a,b):ab=n} b = σ₁(n)` |
| `tsum_fiber_card` / `tsum_fiber_eq` | ★★★★★繊維の和 |
| `tsum_nat_tateXterm_eq_sigma` | ★★★★★★★**`s₁` の q 展開(ℂ の上)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★約数対 -/

theorem one_le_mul_succ (m d : ℕ) : 1 ≤ (m + 1) * (d + 1) :=
  Nat.one_le_iff_ne_zero.2 (by positivity)

theorem pos_of_mem_divAntidiag {N : ℕ} {x : ℕ × ℕ}
    (hx : x ∈ Nat.divisorsAntidiagonal (N + 1)) : 1 ≤ x.1 ∧ 1 ≤ x.2 := by
  rw [Nat.mem_divisorsAntidiagonal] at hx
  obtain ⟨hprod, _⟩ := hx
  constructor
  · rcases Nat.eq_zero_or_pos x.1 with h | h
    · rw [h, zero_mul] at hprod
      omega
    · exact h
  · rcases Nat.eq_zero_or_pos x.2 with h | h
    · rw [h, mul_zero] at hprod
      omega
    · exact h

/-- ★★★★**繊維は約数対の集合と同じもの**。

★`Σ` の中で直接作ると第二成分が依存型になって壊れる。
繊維分解は `Equiv.sigmaFiberEquiv` に任せ、**繊維ごとの同値だけ**を手で作る。 -/
def fiberEquiv (N : ℕ) :
    {p : ℕ × ℕ // (p.1 + 1) * (p.2 + 1) - 1 = N}
      ≃ ((N + 1).divisorsAntidiagonal : Finset (ℕ × ℕ)) where
  toFun p := ⟨(p.val.1 + 1, p.val.2 + 1), by
    rw [Nat.mem_divisorsAntidiagonal]
    have h1 := one_le_mul_succ p.val.1 p.val.2
    have h2 := p.property
    refine ⟨?_, ?_⟩
    · show (p.val.1 + 1) * (p.val.2 + 1) = N + 1
      omega
    · omega⟩
  invFun x := ⟨(x.val.1 - 1, x.val.2 - 1), by
    obtain ⟨ha, hb⟩ := pos_of_mem_divAntidiag x.property
    have hprod := (Nat.mem_divisorsAntidiagonal.1 x.property).1
    have e1 : x.val.1 - 1 + 1 = x.val.1 := by omega
    have e2 : x.val.2 - 1 + 1 = x.val.2 := by omega
    show (x.val.1 - 1 + 1) * (x.val.2 - 1 + 1) - 1 = N
    rw [e1, e2, hprod]
    omega⟩
  left_inv p := by
    apply Subtype.ext
    simp
  right_inv x := by
    apply Subtype.ext
    obtain ⟨ha, hb⟩ := pos_of_mem_divAntidiag x.property
    have e1 : x.val.1 - 1 + 1 = x.val.1 := by omega
    have e2 : x.val.2 - 1 + 1 = x.val.2 := by omega
    simp only [e1, e2]

/-- ★★**`∑_{(a,b) : ab = n} b = σ₁(n)`**。 -/
theorem sum_divAntidiag_snd (n : ℕ) :
    ∑ x ∈ Nat.divisorsAntidiagonal n, x.2 = ArithmeticFunction.sigma 1 n := by
  rw [Nat.sum_divisorsAntidiagonal (fun _ b => b)]
  rw [show (∑ i ∈ n.divisors, n / i) = ∑ i ∈ n.divisors, id (n / i) from rfl,
    Nat.sum_div_divisors n id, ArithmeticFunction.sigma_apply]
  simp

/-! ## ★★★★★繊維の和 -/

theorem tsum_fiber_card (N : ℕ) :
    (∑' x : {p : ℕ × ℕ // (p.1 + 1) * (p.2 + 1) - 1 = N}, ((x.val.2 + 1 : ℕ) : ℂ))
      = ((ArithmeticFunction.sigma 1 (N + 1) : ℕ) : ℂ) := by
  have hEq : (∑' x : {p : ℕ × ℕ // (p.1 + 1) * (p.2 + 1) - 1 = N},
      ((x.val.2 + 1 : ℕ) : ℂ))
      = ∑' b : ((N + 1).divisorsAntidiagonal : Finset (ℕ × ℕ)), ((b.val.2 : ℕ) : ℂ) := by
    rw [← (fiberEquiv N).tsum_eq
      (fun b : ((N + 1).divisorsAntidiagonal : Finset (ℕ × ℕ)) => ((b.val.2 : ℕ) : ℂ))]
    exact tsum_congr fun c => rfl
  rw [hEq, Finset.tsum_subtype ((N + 1).divisorsAntidiagonal)
    (fun x : ℕ × ℕ => ((x.2 : ℕ) : ℂ)), ← Nat.cast_sum, sum_divAntidiag_snd]

theorem tsum_fiber_eq (q : ℂ) (N : ℕ) :
    (∑' c : {p : ℕ × ℕ // (p.1 + 1) * (p.2 + 1) - 1 = N},
        ((c.val.2 + 1 : ℕ) : ℂ) * q ^ ((c.val.1 + 1) * (c.val.2 + 1)))
      = ((ArithmeticFunction.sigma 1 (N + 1) : ℕ) : ℂ) * q ^ (N + 1) := by
  have hstep : ∀ c : {p : ℕ × ℕ // (p.1 + 1) * (p.2 + 1) - 1 = N},
      ((c.val.2 + 1 : ℕ) : ℂ) * q ^ ((c.val.1 + 1) * (c.val.2 + 1))
        = ((c.val.2 + 1 : ℕ) : ℂ) * q ^ (N + 1) := by
    intro c
    have h1 := one_le_mul_succ c.val.1 c.val.2
    have h2 := c.property
    rw [show (c.val.1 + 1) * (c.val.2 + 1) = N + 1 by omega]
  rw [tsum_congr hstep, tsum_mul_right, tsum_fiber_card]

/-! ## ★★★★★★★`s₁` の q 展開 -/

/-- ★★★★★★★**`∑_{m≥1} f(qᵐ) = ∑_{N≥1} σ₁(N) qᴺ`(ℂ の上)**。

★形式側の第 111 ブロック(`tateXtail_one`)に対応する解析側の主張である。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tsum_nat_tateXterm_eq_sigma (q : ℂ) (hq : ‖q‖ < 1) :
    (∑' m : ℕ, tateXterm (q ^ (m + 1)))
      = ∑' N : ℕ, ((ArithmeticFunction.sigma 1 (N + 1) : ℕ) : ℂ) * q ^ (N + 1) := by
  have hsum : Summable (fun p : ℕ × ℕ => ((p.2 + 1 : ℕ) : ℂ) * q ^ ((p.1 + 1) * (p.2 + 1))) :=
    summable_double_family q hq
  have hcomp : Summable ((fun p : ℕ × ℕ => ((p.2 + 1 : ℕ) : ℂ) * q ^ ((p.1 + 1) * (p.2 + 1)))
      ∘ (Equiv.sigmaFiberEquiv (fun x : ℕ × ℕ => (x.1 + 1) * (x.2 + 1) - 1))) :=
    (Equiv.summable_iff _).2 hsum
  rw [← tsum_double_family_eq q hq,
    ← (Equiv.sigmaFiberEquiv (fun x : ℕ × ℕ => (x.1 + 1) * (x.2 + 1) - 1)).tsum_eq
      (fun p : ℕ × ℕ => ((p.2 + 1 : ℕ) : ℂ) * q ^ ((p.1 + 1) * (p.2 + 1)))]
  refine Eq.trans ?_ (tsum_congr fun N => tsum_fiber_eq q N)
  exact hcomp.tsum_sigma

/-! ## ★出典の紐付け(`.src`) -/

def fiberEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——繊維と約数対の同値)",
    sectionId := "genell-def-3-3" }

def tsum_nat_tateXterm_eq_sigma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——s1 の q 展開、ℂ の上)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
