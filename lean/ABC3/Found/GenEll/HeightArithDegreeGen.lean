/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightArithDegree
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★`deg_fin + deg_arch = log H` —— チャートの仮定を外す（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★★★これは何か —— `chart` の仮定の**算術の側**を外す

`§9-866` の `log_degFin_add_degArch_eq_log_mulHeight` は仮定に

    `hI : Ideal.span (Set.range x) = Ideal.span {x i₀}`

——**座標の生成するイデアルが主イデアルで、しかも `x_{i₀}` が生成する**——を持っていた。
★これは `Spec 𝓞_F ⟶ ℙᴺ` が**チャート `D₊(x_{i₀})` を通る**ことの代数的な形であり、
`§9-928` で測ったとおり**イデアル類が自明でないと成り立たない**。

★★★★**本ファイルはそれを外す**。必要なのは**ノルムの等式**だけである:

    `N(I) · N(𝔞) = N((x_j))`   （`𝔞 ≔ (x_0,…,x_N)`）

——これは `I = x_j·𝔞⁻¹`（＝古典的な「超平面の引き戻し」）が満たす等式そのものである。

## ★★★機構 —— もともと `hI` は 1 行しか使っていなかった

`§9-866` の証明は

* `log_mulHeight_add_log_absNorm`（`log H + log N(𝔞) = log ∏_v sup^mult`）——**一般**
* `prod_infinitePlace_eq_absNorm`（`∏_v v(x_j)^mult = N((x_j))`）——**一般**

の 2 本でできており、`hI` は `N(𝔞)` を `N((x_{i₀}))` に**書き換えるためだけ**に使われていた。
★★そこを「`N(I)·N(𝔞) = N((x_j))`」で受ければ、`hI` は要らない。

## ★これで何が残ったか

★★★`Proposition 1.4, (iv)` の `chart` 仮定のうち、**算術の側は外れた**。
残るのは**幾何の側**——

    `pullbackIdeal F (超平面) xF` のノルムが `N((x_0))/N(𝔞)` であること

——すなわち `pullbackIdeal = x_0·𝔞⁻¹` である。
★道筋は `§9-928` に記録したとおり、`Spec 𝓞_F` を `{xF⁻¹(X_{s_i})}` の中の
基本開集合で細分して貼ることである。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★ノルムの等式だけで足りる -/

/-- ★★★★★★★★★★★★**`hI`（チャート）を外した形**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★仮定は**ノルムの等式** `N(I)·N(𝔞) = N((x_j))` だけである
（`𝔞 ≔ (x_0,…,x_N)`）。
★★`§9-866` では `hI : 𝔞 = (x_{i₀})` を課していたが、
それは `N(𝔞)` を `N((x_{i₀}))` に書き換えるためだけに使われていた。 -/
theorem log_absNorm_add_archHeight' (K : Type) [Field K] [NumberField K] {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (hx : x ≠ 0) (j : ι) (I : Ideal (𝓞 K)) (hj : x j ≠ 0)
    (hI : I.absNorm * (Ideal.span (Set.range x)).absNorm = (Ideal.span {x j}).absNorm) :
    Real.log ((I.absNorm : ℝ))
        + (Real.log (∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult)
          - Real.log (∏ v : InfinitePlace K, v ((x j : K)) ^ v.mult))
      = Real.log (Height.mulHeight (fun i => (x i : K))) := by
  have hane : ((Ideal.span (Set.range x)).absNorm : ℝ) ≠ 0 := absNorm_span_range_ne_zero K x hx
  have hjne : ((Ideal.span {x j}).absNorm : ℝ) ≠ 0 := absNorm_span_singleton_ne_zero K (x j) hj
  have hIne : ((I.absNorm : ℝ)) ≠ 0 := by
    intro hc
    apply hjne
    rw [← hI]
    push_cast
    rw [hc, zero_mul]
  have h3 := log_mulHeight_add_log_absNorm K x hx
  have h1 : Real.log (∏ v : InfinitePlace K, v ((x j : K)) ^ v.mult)
      = Real.log (((Ideal.span {x j}).absNorm : ℝ)) := by
    rw [prod_infinitePlace_eq_absNorm]
  have h2 : (((Ideal.span {x j}).absNorm : ℝ))
      = ((I.absNorm : ℝ)) * ((Ideal.span (Set.range x)).absNorm : ℝ) := by
    rw [← hI]; push_cast; ring
  rw [h1, h2, Real.log_mul hIne hane]
  linarith [h3]

/-- ★★★★★★★★★★★★★★**段 C2c 後半（チャートの仮定なし）** ——
`deg_fin + deg_arch = log H(x)`。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★これが `Proposition 1.4, (iv)` の高さの計算の**一般形**である
——点がチャートを通らなくてもよい。 -/
theorem log_degFin_add_degArch_eq_log_mulHeight' (K : Type) [Field K] [NumberField K]
    {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (hx : x ≠ 0) (j : ι) (I : Ideal (𝓞 K)) (hj : x j ≠ 0)
    (hI : I.absNorm * (Ideal.span (Set.range x)).absNorm = (Ideal.span {x j}).absNorm) :
    Real.log ((I.absNorm : ℝ))
        + ∑ v : InfinitePlace K,
            (v.mult : ℝ) * Real.log ((⨆ i, v ((x i : K))) / v ((x j : K)))
      = Real.log (Height.mulHeight (fun i => (x i : K))) := by
  rw [← log_arch_ratio K x j hj]
  exact log_absNorm_add_archHeight' K x hx j I hj hI

/-! ## ★★★`§9-866` の形は特別な場合である -/

/-- ★★★**`§9-866` の `hI` からノルムの等式が出る**——一般化であることの確認。 -/
theorem absNorm_mul_of_span_eq (K : Type) [Field K] [NumberField K] {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (i₀ j : ι) (z : 𝓞 K) (hz : x j = z * x i₀)
    (hI : Ideal.span (Set.range x) = Ideal.span {x i₀}) :
    (Ideal.span {z}).absNorm * (Ideal.span (Set.range x)).absNorm
      = (Ideal.span {x j}).absNorm := by
  rw [hI, hz, ← Ideal.span_singleton_mul_span_singleton, Ideal.absNorm.map_mul]

/-! ## ★出典の紐付け(`.src`) -/

def log_absNorm_add_archHeight'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(ノルムの等式だけで deg_fin + deg_arch = log H が出る)",
    sectionId := "genell-prop-1-4" }

def log_degFin_add_degArch_eq_log_mulHeight'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(段 C2c 後半——チャートの仮定なし)",
    sectionId := "genell-prop-1-4" }

def absNorm_mul_of_span_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(チャートの仮定からノルムの等式が出る——一般化であることの確認)",
    sectionId := "genell-prop-1-4" }

def log_degFin_add_degArch_eq_log_mulHeight'.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "log_mulHeight_add_log_absNorm(log H + log N(𝔞) = log ∏_v sup^mult)"
      (.inProject "ABC3" "ABC3.Found.GenEll.log_mulHeight_add_log_absNorm") 2,
    .citation "[ABC3]" "prod_infinitePlace_eq_absNorm(積公式)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prod_infinitePlace_eq_absNorm") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-866 の hI(座標の生成するイデアルが (x_{i₀}))は " ++
       "**N(𝔞) を N((x_{i₀})) に書き換えるためだけ**に使われていた。" ++
       "そこを「N(I)·N(𝔞) = N((x_j))」で受ければ hI は要らない" ++
       "——これは I = x_j·𝔞⁻¹(＝古典的な超平面の引き戻し)が満たす等式そのものである") 4,
    .implicitStep
      ("★★★これで Proposition 1.4, (iv) の chart 仮定のうち**算術の側は外れた**。" ++
       "残るのは幾何の側——pullbackIdeal F (超平面) xF のノルムが N((x_0))/N(𝔞) であること" ++
       "(すなわち pullbackIdeal = x_0·𝔞⁻¹)。" ++
       "★道筋は Spec 𝓞_F を {xF⁻¹(X_{s_i})} の中の基本開集合で細分して貼ることである") 5 ]

end ABC3.Found.GenEll
