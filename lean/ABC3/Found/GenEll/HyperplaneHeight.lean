/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightArithDegree
import ABC3.Found.GenEll.HeightMetric
import ABC3.Found.GenEll.LogDiffValue
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★超平面因子の `htArith` は素朴高さである（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★これは何か —— 段 C2c の到達点

`§9-866`（段 C2c-4）は**実数の等式**

    `log N((x₀/x_{i₀})) + Σ_v mult(v)·log( sup_i v(x_i)/v(x₀) ) = log H(x)`

を取った。★本ファイルはそれを**プロジェクトの高さ `htArith` の言葉**に載せる:

    `htArith F D̄ xF = log H(x) / [F : ℚ]`

（`D̄ = (超平面因子, Fubini–Study 型の Green 関数)`）。

## ★★★機構 —— 在庫の分解に §9-866 を差し込むだけ

★`§9-Heights` の `htArith_eq_add`:

    `htArith F D xF = degNormalized (idealADiv F (pullbackIdeal F D.divisor xF))
                       + (archADiv F D.green xF).sum / [F:ℚ]`

★★有限側は `§9-865`（引き戻しは `(x₀/x_{i₀})`）と `deg_idealADiv`
（`deg (idealADiv I) = log N(I)`）で `log N((x₀/x_{i₀}))` になる。
★★★アルキメデス側は `archADiv_apply` で `Σ_v mult(v)·g(archPoint xF v)` になり、
`g` が Fubini–Study なら `§9-866` の第 2 項そのものである。

## ★★★★★測定の記録 —— 層 `O(1)` は要らなかった

★台帳は当初「`O(1)` の層が mathlib に無い」ことを段 C2c の障害としていたが、
**因子表示 `ArithCartier`（イデアル層 ＋ Green 関数）で足りる**。
★★`GreenFn X` は `complexPoints X → ℝ` という**ただの実関数**なので、
Fubini–Study を載せるのに層も計量の構造も要らない（2026-08-28 実測）。

## ★残っている段（明示）

★★★★★仮定 `hgreen`（`g (archPoint xF v) = log( sup_i v(x_i)/v(x₀) )`）が残る。
これは「`v` が定める複素点で座標を読むと、その絶対値が `v(x_i)` である」という
**素点と複素点の対応**であり、`archPoint` の定義（`archRingHom v` を前合成）と
`projPointCoord`（`§9-C2b`）を突き合わせれば出る。★新しい数学は要らない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★アルキメデス側の総和 -/

/-- ★**`archADiv` の総和は無限素点全体の和である**。

★`Finsupp.sum` は台の上の和だが、台の外の項は `0` なので `univ` の和に等しい。 -/
theorem archADiv_sum_eq {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) :
    (archADiv F g xF).sum (fun _ r => r)
      = ∑ v : InfinitePlace F, (InfinitePlace.mult v : ℝ) * g (archPoint xF v) := by
  classical
  rw [Finsupp.sum]
  have hsub : (archADiv F g xF).support ⊆ (Finset.univ : Finset (InfinitePlace F)) :=
    fun v _ => Finset.mem_univ v
  rw [Finset.sum_subset hsub (fun v _ hv => Finsupp.notMem_support_iff.1 hv)]
  exact Finset.sum_congr rfl (fun v _ => archADiv_apply F g xF v)

/-! ## ★★★★★★★★★★★★本体 -/

/-- ★★★★★★★★★★★★**超平面因子の `htArith` は素朴高さである**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `htArith F D̄ xF = log H(x) / [F : ℚ]`

★仮定は 4 つで、いずれも既に取れている／取れる形である:

| 仮定 | 出どころ |
|---|---|
| `hI`（点がチャート `D₊(x_{i₀})` を通る） | `§9-C2b` ＋ `§9-854` |
| `hz`・`hj`（`x₀ = z·x_{i₀}`、点は超平面上にない） | 同上 |
| `hpull`（引き戻しは `(z)`） | ★`§9-865`（段 C2c-1） |
| `hgreen`（Green 関数は Fubini–Study） | ★★残る段（素点と複素点の対応） |

★★★機構は `htArith_eq_add` の 2 項に `deg_idealADiv` と `§9-866` を差し込むだけである。 -/
theorem htArith_eq_log_mulHeight {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) {ι : Type} [Finite ι]
    (x : ι → 𝓞 F) (hx : x ≠ 0) (i₀ j : ι) (z : 𝓞 F)
    (hz : x j = z * x i₀) (hj : x j ≠ 0)
    (hI : Ideal.span (Set.range x) = Ideal.span {x i₀})
    (hpull : pullbackIdeal F D.divisor xF = Ideal.span {z})
    (hgreen : ∀ v : InfinitePlace F,
      D.green (archPoint xF v) = Real.log ((⨆ i, v ((x i : F))) / v ((x j : F)))) :
    htArith F D xF
      = Real.log (Height.mulHeight (fun i => (x i : F))) / (Module.finrank ℚ F : ℝ) := by
  have hzne : z ≠ 0 := by rintro rfl; exact hj (by rw [hz, zero_mul])
  have hspan : Ideal.span {z} ≠ (0 : Ideal (𝓞 F)) := by
    intro h0
    exact hzne (by simpa using (Ideal.span_singleton_eq_bot.1 h0))
  rw [htArith_eq_add F D xF, hpull, degNormalized, deg_idealADiv F _ hspan,
    archADiv_sum_eq F D.green xF]
  have hg : ∑ v : InfinitePlace F, (InfinitePlace.mult v : ℝ) * D.green (archPoint xF v)
      = ∑ v : InfinitePlace F,
          (InfinitePlace.mult v : ℝ) * Real.log ((⨆ i, v ((x i : F))) / v ((x j : F))) :=
    Finset.sum_congr rfl (fun v _ => by rw [hgreen v])
  rw [hg, ← add_div,
    log_degFin_add_degArch_eq_log_mulHeight F x hx i₀ j z hz hj hI]

/-- ★★**指数を取った形** —— `Northcott` の消費側（`hcmp`）が読む向き。

★`Height.mulHeight` は**相対高さ**（`H_K`）なので `[F:ℚ]` 倍が付く。 -/
theorem mulHeight_eq_exp_htArith {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) {ι : Type} [Finite ι]
    (x : ι → 𝓞 F) (hx : x ≠ 0) (i₀ j : ι) (z : 𝓞 F)
    (hz : x j = z * x i₀) (hj : x j ≠ 0)
    (hI : Ideal.span (Set.range x) = Ideal.span {x i₀})
    (hpull : pullbackIdeal F D.divisor xF = Ideal.span {z})
    (hgreen : ∀ v : InfinitePlace F,
      D.green (archPoint xF v) = Real.log ((⨆ i, v ((x i : F))) / v ((x j : F)))) :
    Height.mulHeight (fun i => (x i : F))
      = Real.exp ((Module.finrank ℚ F : ℝ) * htArith F D xF) := by
  have hpos : (0:ℝ) < Height.mulHeight (fun i => (x i : F)) :=
    lt_of_lt_of_le zero_lt_one (Height.one_le_mulHeight _)
  have hfr : ((Module.finrank ℚ F : ℕ) : ℝ) ≠ 0 := by
    have : (0:ℝ) < (Module.finrank ℚ F : ℝ) := by exact_mod_cast Module.finrank_pos
    exact this.ne'
  rw [htArith_eq_log_mulHeight F D xF x hx i₀ j z hz hj hI hpull hgreen,
    mul_div_cancel₀ _ hfr, Real.exp_log hpos]

/-! ## ★出典の紐付け(`.src`) -/

def archADiv_sum_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(archADiv の総和は無限素点全体の和である)",
    sectionId := "genell-prop-1-4" }

def htArith_eq_log_mulHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(超平面因子の htArith は素朴高さである)",
    sectionId := "genell-prop-1-4" }

def mulHeight_eq_exp_htArith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(素朴高さは htArith の指数である——Northcott が読む向き)",
    sectionId := "genell-prop-1-4" }

def htArith_eq_log_mulHeight.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "log_degFin_add_degArch_eq_log_mulHeight(段 C2c-4、§9-866)"
      (.inProject "ABC3" "ABC3.Found.GenEll.log_degFin_add_degArch_eq_log_mulHeight") 2,
    .citation "[ABC3]" "pullbackIdeal_hyperplane_point(引き戻しは (x₀/x_i)、段 C2c-1、§9-865)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdeal_hyperplane_point") 2,
    .citation "[ABC3]" "htArith_eq_add(高さ = 有限側 + アルキメデス側)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_eq_add") 2,
    .citation "[ABC3]" "deg_idealADiv(deg (idealADiv I) = log N(I))"
      (.inProject "ABC3" "ABC3.Found.GenEll.deg_idealADiv") 2,
    .implicitStep
      ("★★★★★測定: 台帳は当初「O(1) の層が mathlib に無い」ことを段 C2c の障害と" ++
       "していたが、**因子表示 ArithCartier(イデアル層 ＋ Green 関数)で足りる**。" ++
       "GreenFn X は complexPoints X → ℝ という**ただの実関数**なので、" ++
       "Fubini–Study を載せるのに層も計量の構造も要らない(2026-08-28 実測)") 4,
    .implicitStep
      ("★★仮定 hgreen(g (archPoint xF v) = log( sup_i v(x_i)/v(x₀) ))が残る。" ++
       "これは「v が定める複素点で座標を読むと、その絶対値が v(x_i) である」という" ++
       "**素点と複素点の対応**であり、archPoint の定義(archRingHom v を前合成)と " ++
       "projPointCoord(§9-C2b)を突き合わせれば出る。★新しい数学は要らない") 4 ]

end ABC3.Found.GenEll
