/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SigmaBoundN

/-!
# [GenEll] Lemma 3.7 の最初の観察 —— **次数は 1 素点の寄与を超える**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Then if (a) is satisfied, then l is > the local heights of EL at all the primes of

## ★★★原文が「observe」で畳んだ 1 行

`Lemma 3.7` の証明の 1 行目である。★中身は **剰余体が少なくとも 2 元をもつ**ことだけ:

    deg_F(a) = ∑_w a_w · log #(𝓞_F/w) ≥ a_v · log 2

有効因子（`0 ≤ a_w`）なら、和のうち `v` の項以外はすべて非負なので落とせる。
★★`log #(𝓞_F/w) ≥ log 2` は `1 < absNorm w`（mathlib の
`NumberField.HeightOneSpectrum.one_lt_absNorm`）から。

## ★★★★これが `Lemma 3.7` のどこで効くか

原文はこの観察と `Proposition 3.4` を合わせて

> l ≥ 100d · deg([EL]) / 14 ≥ 100 log(2) / 14 · v > v

を出す——★**「`l` が局所高さより大きい」**という結論である。
`Proposition 3.4` の側は `[Silv2]` に依るので未だが、**本観察は独立に取れる**。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★原文の `deg([E_L])` は `M_ell` 上の高さ（Faltings 高さの無限素点での寄与）であり、
本ファイルの `deg` は算術因子の次数である。
★★**両者を繋ぐのは `Definition 3.3` の `deg_∞(E) = v_K(q_E)·log #(𝓞_K/𝔪)`** であり、
その接続は本ファイルには入っていない——ここで取るのは**不等式の骨**だけである。
-/

namespace ABC3.Found.GenEll

open NumberField

/-- ★★**剰余体は少なくとも 2 元をもつ** —— `log #(𝓞_F/v) ≥ log 2`。

原文 (GenEll p.18):
> Then if (a) is satisfied, then l is > the local heights of EL at all the primes of

★mathlib の `NumberField.HeightOneSpectrum.one_lt_absNorm` から。 -/
theorem log_two_le_log_residueCard {F : Type} [Field F] [NumberField F]
    (v : FinitePlace F) :
    Real.log 2 ≤ Real.log (residueCard v) := by
  have h1 : 1 < residueCard v := NumberField.HeightOneSpectrum.one_lt_absNorm v
  have h2 : (2 : ℝ) ≤ (residueCard v : ℝ) := by exact_mod_cast h1
  exact Real.log_le_log (by norm_num) h2

/-- ★★★★★**[GenEll] Lemma 3.7 の最初の観察** ——
有効因子の次数は、どの 1 素点の寄与よりも大きい。

原文 (GenEll p.18):
> Then if (a) is satisfied, then l is > the local heights of EL at all the primes of

★★★**`deg_F(a) ≥ a_v · log 2`**。中身は「和の他の項が非負」と
「剰余体が 2 元以上」だけである。

★`v` が台に無い場合も込みで扱う（そのときは左辺が `0`）。 -/
theorem coeff_mul_log_two_le_deg {F : Type} [Field F] [NumberField F]
    (a : ADiv F) (harc : a.arc = 0) (heff : ∀ w, 0 ≤ a.fin w)
    (v : FinitePlace F) :
    ((a.fin v : ℤ) : ℝ) * Real.log 2 ≤ deg a := by
  classical
  have harc0 : a.arc.sum (fun _ r => r) = 0 := by rw [harc]; simp
  rw [deg, harc0, add_zero, Finsupp.sum]
  by_cases hv : v ∈ a.fin.support
  · have hterm : ((a.fin v : ℤ) : ℝ) * Real.log 2
        ≤ ((a.fin v : ℤ) : ℝ) * Real.log (residueCard v) := by
      have h0 : (0:ℝ) ≤ ((a.fin v : ℤ) : ℝ) := by exact_mod_cast heff v
      exact mul_le_mul_of_nonneg_left (log_two_le_log_residueCard v) h0
    refine le_trans hterm (Finset.single_le_sum (f := fun w =>
      ((a.fin w : ℤ) : ℝ) * Real.log (residueCard w)) (fun w _ => ?_) hv)
    have h0 : (0:ℝ) ≤ ((a.fin w : ℤ) : ℝ) := by exact_mod_cast heff w
    exact mul_nonneg h0 (Real.log_natCast_nonneg _)
  · have hz : a.fin v = 0 := Finsupp.notMem_support_iff.1 hv
    rw [hz]
    simp only [Int.cast_zero, zero_mul]
    refine Finset.sum_nonneg (fun w _ => ?_)
    have h0 : (0:ℝ) ≤ ((a.fin w : ℤ) : ℝ) := by exact_mod_cast heff w
    exact mul_nonneg h0 (Real.log_natCast_nonneg _)

/-- ★★★**正規化した形** —— `[F:ℚ]·deg_F(a) ≥ a_v · log 2`。

★原文の `d · deg([E_L]) ≥ v · log(2)` の `d` はまさに `[L:ℚ]` である。 -/
theorem coeff_mul_log_two_le_finrank_mul_degNormalized {F : Type} [Field F] [NumberField F]
    (a : ADiv F) (harc : a.arc = 0) (heff : ∀ w, 0 ≤ a.fin w)
    (v : FinitePlace F) :
    ((a.fin v : ℤ) : ℝ) * Real.log 2
      ≤ (Module.finrank ℚ F : ℝ) * degNormalized a := by
  have hF : (0:ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have h := coeff_mul_log_two_le_deg a harc heff v
  rw [degNormalized, mul_div_cancel₀ _ (ne_of_gt hF)]
  exact h

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Lemma 3.7` の残りは
`Proposition 3.4`（`[Silv2]` に依る）と、例外集合 `Exc` の Galois-有限性である。 -/

def log_two_le_log_residueCard.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(剰余体は少なくとも 2 元をもつ)",
    sectionId := "genell-lemma-3-7" }

def coeff_mul_log_two_le_deg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(最初の観察——次数は 1 素点の寄与を超える)",
    sectionId := "genell-lemma-3-7" }

def coeff_mul_log_two_le_deg.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "NumberField.HeightOneSpectrum.one_lt_absNorm(剰余体は 2 元以上)"
      (.inMathlib "NumberField.HeightOneSpectrum.one_lt_absNorm") 18,
    .implicitStep
      ("★原文が「observe」で畳んだ 1 行である。中身は「和の他の項が非負」と" ++
       "「剰余体が 2 元以上」だけ") 18,
    .implicitStep
      ("★逸脱: 原文の deg([E_L]) は M_ell 上の高さ(Faltings 高さの無限素点での寄与)であり、" ++
       "本ファイルの deg は算術因子の次数である。両者を繋ぐのは Definition 3.3 の " ++
       "deg_∞(E) = v_K(q_E)·log #(𝓞_K/𝔪) だが、その接続は本ファイルには入っていない" ++
       "——ここで取るのは不等式の骨だけである") 18,
    .implicitStep
      ("★★Lemma 3.7 の残りは Proposition 3.4([Silv2] に依る)と、" ++
       "例外集合 Exc の Galois-有限性である") 18 ]

def coeff_mul_log_two_le_finrank_mul_degNormalized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(正規化した形——原文の d·deg ≥ v·log 2)",
    sectionId := "genell-lemma-3-7" }

end ABC3.Found.GenEll
