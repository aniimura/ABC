/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SigmaBound

/-!
# [GenEll] Proposition 1.4, (iii) / Proposition 1.7 —— **係数が有界な因子の一様上界**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★`SigmaBound.lean` の一般化

`SigmaBound.lean` の `degNormalized_le_sum_log` / `degNormalized_sub_le_sum_log` は
**係数が `≤ 1`** であることを使う——導手が `(−)_red` で被約だからである。

★★しかし `Proposition 1.4, (iii)` が要るのは**垂直因子**の寄与の限界である:

原文 (GenEll p.6):
> line bundle LQ on XQ. In particular, it makes sense to write [htL] or [htLQ] for the

`D_ℚ ∼ E_ℚ`（`X_ℚ` の上で線形同値）なら `D − E = div(f) + V` で `V` は
**有限個のファイバーに台をもつ垂直因子**である。★`V` の係数は `1` 以下とは限らない
——固定の整数 `n_i` である。

★★★そこで**係数の絶対値が `N` で抑えられる場合**に一般化する。上界は `N` 倍になるだけで、
**`F` にも点にも依らない**という肝心の性質は保たれる。

## ★★★★到達点

| 主張 | 宣言 |
|---|---|
| 係数が `N` 以下の因子の次数 | ★`abs_degNormalized_le_of_bounded` |
| 2 つの因子の差（高さの差の形） | ★★`abs_sub_le_of_bounded_diff` |

★`SigmaBound.lean` の `sum_log_residueCard_le_of_over`（台が `Σ` の上にある集まりの
`log N v` の総和の限界）をそのまま使う——**そこが再利用できる形に切り出されていた**。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★**「垂直因子である」ことそのものは仮定として受ける**——
`a.fin − b.fin` の台が `Σ` の上にあり係数が `N` 以下、という形で。
★★`D_ℚ ∼ E_ℚ ⟹ D − E = div(f) + V` という分解（`Div(X) → Div(X_ℚ)` の全射性と
その核が垂直因子であること）は**別の仕事**である。
-/

namespace ABC3.Found.GenEll

open NumberField

variable {F : Type} [Field F] [NumberField F]

/-- ★★★★★★**係数が `N` で抑えられる因子の次数の一様上界**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★`SigmaBound.lean` の `degNormalized_le_sum_log` は係数 `≤ 1` を使う。
本定理は `|係数| ≤ N` に一般化したもので、上界が `N` 倍になるだけである。

★★★**`F` にも点にも依らない**——それが「一様」の中身である。 -/
theorem abs_degNormalized_le_of_bounded
    (N : ℕ) (a : ADiv F) (harc : a.arc = 0)
    (hN : ∀ v, |a.fin v| ≤ (N : ℤ))
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ a.fin.support, ch v ∈ Sig)
    (hover : ∀ v ∈ a.fin.support, (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    |degNormalized a| ≤ (N : ℝ) * ∑ q ∈ Sig, Real.log q := by
  classical
  have hF : (0:ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have harc0 : a.arc.sum (fun _ r => r) = 0 := by rw [harc]; simp
  have hbound : |deg a| ≤ (N : ℝ) * ((Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q) := by
    rw [deg, harc0, add_zero, Finsupp.sum]
    calc |∑ v ∈ a.fin.support, ((a.fin v : ℤ) : ℝ) * Real.log (residueCard v)|
        ≤ ∑ v ∈ a.fin.support, |((a.fin v : ℤ) : ℝ) * Real.log (residueCard v)| :=
          Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ v ∈ a.fin.support, (N : ℝ) * Real.log (residueCard v) := by
          refine Finset.sum_le_sum (fun v hv => ?_)
          rw [abs_mul, abs_of_nonneg (Real.log_natCast_nonneg _)]
          refine mul_le_mul_of_nonneg_right ?_ (Real.log_natCast_nonneg _)
          rw [← Int.cast_abs]
          exact_mod_cast hN v
      _ = (N : ℝ) * ∑ v ∈ a.fin.support, Real.log (residueCard v) := by rw [Finset.mul_sum]
      _ ≤ (N : ℝ) * ((Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q) := by
          refine mul_le_mul_of_nonneg_left ?_ (Nat.cast_nonneg N)
          exact sum_log_residueCard_le_of_over _ Sig hprime ch hmem hover
  rw [degNormalized, abs_div, abs_of_pos hF, div_le_iff₀ hF]
  nlinarith [hbound]

/-- ★★★★★★★**2 つの高さの差** —— 引き戻しの差が `Σ` 上で係数 `N` 以下なら BD-同値。

原文 (GenEll p.6):
> line bundle LQ on XQ. In particular, it makes sense to write [htL] or [htLQ] for the

★★★これが `Proposition 1.4, (iii)` の**垂直因子の側**である
——`D_ℚ ∼ E_ℚ` なら `D − E` は主因子と垂直因子の和で、垂直因子の寄与が本定理で抑えられる。

★**上界 `N · ∑_{q∈Σ} log q` は `F` にも点にも依らない**ので BD-同値が出る。 -/
theorem abs_sub_le_of_bounded_diff
    (N : ℕ) (a b : ADiv F) (harc : a.arc = b.arc)
    (hN : ∀ v, |a.fin v - b.fin v| ≤ (N : ℤ))
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ (a.fin - b.fin).support, ch v ∈ Sig)
    (hover : ∀ v ∈ (a.fin - b.fin).support,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    |degNormalized a - degNormalized b| ≤ (N : ℝ) * ∑ q ∈ Sig, Real.log q := by
  have hsub : degNormalized a - degNormalized b = degNormalized (a - b) := by
    have hadd := degNormalized_add (a - b) b
    have hab : a - b + b = a := by abel
    rw [hab] at hadd
    linarith
  rw [hsub]
  refine abs_degNormalized_le_of_bounded N (a - b) ?_ ?_ Sig hprime ch ?_ ?_
  · show a.arc - b.arc = 0
    rw [harc]; abel
  · intro v; exact hN v
  · intro v hv; exact hmem v hv
  · intro v hv; exact hover v hv

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Proposition 1.4, (iii)` には
`D_ℚ ∼ E_ℚ ⟹ D − E = div(f) + V`（`Div(X) → Div(X_ℚ)` の全射性と核が垂直因子であること）
が要り、本ファイルは **`V` の寄与の限界**だけである。 -/

def abs_degNormalized_le_of_bounded.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.4, (iii)(係数が N で抑えられる因子の次数の一様上界)",
    sectionId := "genell-prop-1-4" }

def abs_sub_le_of_bounded_diff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.4, (iii)(垂直因子の側——分解 D − E = div(f) + V は含まない)",
    sectionId := "genell-prop-1-4" }

def abs_sub_le_of_bounded_diff.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "sum_log_residueCard_le_of_over(台が Σ の上にある集まりの限界)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sum_log_residueCard_le_of_over") 9,
    .citation "[ABC3]" "abs_degNormalized_le_of_bounded(係数 N 版の一様上界)"
      (.inProject "ABC3" "ABC3.Found.GenEll.abs_degNormalized_le_of_bounded") 9,
    .implicitStep
      ("★★SigmaBound.lean の degNormalized_le_sum_log は係数 ≤ 1(導手が (−)_red で" ++
       "被約だから)を使う。Proposition 1.4, (iii) が要る垂直因子の係数は 1 以下とは" ++
       "限らないので、|係数| ≤ N に一般化した。上界が N 倍になるだけで" ++
       "「F にも点にも依らない」性質は保たれる") 9,
    .implicitStep
      ("★逸脱: 「垂直因子である」ことそのものは仮定として受ける(台が Σ の上にあり" ++
       "係数が N 以下、という形)。D_ℚ ∼ E_ℚ ⟹ D − E = div(f) + V という分解" ++
       "(Div(X) → Div(X_ℚ) の全射性とその核が垂直因子であること)は別の仕事である") 9 ]

end ABC3.Found.GenEll
