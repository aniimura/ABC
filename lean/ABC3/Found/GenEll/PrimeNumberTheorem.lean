import ABC3.Found.GenEll.PrimesOfSize
import PrimeNumberTheoremAnd.Consequences

/-!
# [GenEll] Remark 4.1.1 の第 2 項 —— `θ(x)/x → 1`(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.21。

原文 (GenEll p.21):
> which satisfy condition (ii) of Lemma 4.1 is entirely elementary. On the other hand, with regard to condition (i), the fact that

## ★★これは我々が証明したものではない

原文は素数定理を証明せず、外部文献([Edw], p. 76)を指す。
★我々も同じく**外部から借りる**——`PrimeNumberTheoremAnd` の
`chebyshev_asymptotic : θ ~[atTop] id` である。

★★**`θ` は mathlib の `Chebyshev.theta` そのもの**である
(`PrimeNumberTheoremAnd/Consequences.lean` は `open scoped Chebyshev` して
`Chebyshev.theta_le_psi` / `Chebyshev.theta_eq_sum_Icc` を使う。2026-08-17 実測)。
ゆえに橋は要らず、漸近同値を商の極限に直すだけである。

## ★橋

mathlib の `Asymptotics.isEquivalent_iff_tendsto_one`:

> `(∀ᶠ x in l, v x ≠ 0) → (u ~[l] v ↔ Tendsto (u / v) l (𝓝 1))`

`v = id` は `atTop` で最終的に非零なので、そのまま適用できる。
-/

namespace ABC3.Found.GenEll

open Filter Asymptotics

/-- ★★**`θ(x)/x → 1`** —— 素数定理の帰結。

★★**我々の証明ではない。** `PrimeNumberTheoremAnd.chebyshev_asymptotic` を
商の極限の形に直しただけである。原文も同じく外部文献を指している。 -/
theorem theta_div_tendsto_one :
    Tendsto (fun x : ℝ => Chebyshev.theta x / x) atTop (nhds 1) := by
  have hz : ∀ᶠ x : ℝ in atTop, id x ≠ 0 := by
    filter_upwards [eventually_gt_atTop (0 : ℝ)] with x hx
    exact ne_of_gt hx
  have h := (isEquivalent_iff_tendsto_one hz).1 chebyshev_asymptotic
  have heq : (Chebyshev.theta / id : ℝ → ℝ) = fun x : ℝ => Chebyshev.theta x / x := by
    funext x; rfl
  rwa [heq] at h

/-! ## ★★[GenEll] Remark 4.1.1 —— 2 項をまとめる -/

/-- ★★★**[GenEll] Remark 4.1.1**。

原文 (GenEll p.21):
> which satisfy condition (ii) of Lemma 4.1 is entirely elementary. On the other hand, with regard to condition (i), the fact that

原文は 2 つのことを言う:

1. `Lemma 4.1` の条件 (ii) を満たす `ϵ, x_ϵ, C_ϵ` を見つけるのは **entirely elementary**
2. `θ(x)/x → 1` は **素数定理のよく知られた帰結**である(cf. [Edw], p. 76)

## ★★出どころを分けて書く

- **第 1 項は我々が証明した**(`PrimesOfSize.lean` の `exists_xeps_cond_ii`)。
  `Real.isLittleO_log_id_atTop` である。
- ★★**第 2 項は外部から借りた**(`PrimeNumberTheoremAnd.chebyshev_asymptotic`)。
  **我々は素数定理を証明していない。原文も証明していない**——[Edw] を指すだけである。
  ★借用先は `#print axioms` で標準 3 公理のみを実測した(2026-08-17)。
  同 repo の `Wiener.lean` に sorry 2 件があるが**経路外**であることも
  同じ `#print axioms` が示している。 -/
theorem remark_4_1_1 :
    (∀ (M : ℕ) (eps : ℝ), 0 < eps →
        ∃ xeps : ℝ, 0 < xeps ∧ ∀ x : ℝ, xeps ≤ x → (M : ℝ) * Real.log x ≤ eps * x)
  ∧ Tendsto (fun x : ℝ => Chebyshev.theta x / x) atTop (nhds 1) :=
  ⟨exists_xeps_cond_ii, theta_div_tendsto_one⟩

/-! ## ★出典の紐付け(`.src`) -/

def theta_div_tendsto_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21,
    item := "Remark 4.1.1(第 2 項——素数定理は外部から借りた)",
    sectionId := "genell-rem-4-1-1" }

/-- ★★**条なしの `.src`** —— `Remark 4.1.1` の 2 項が揃った。

★★ただし**第 2 項は借り物である**ことを上の docstring に明記してある。
原文も外部文献を指しているので、原文の主張としては完全に写せている。 -/
def remark_4_1_1.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Remark 4.1.1",
    sectionId := "genell-rem-4-1-1" }

end ABC3.Found.GenEll
