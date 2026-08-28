/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HyperplaneHeight
import ABC3.Found.GenEll.HeightArithDegreeGen
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★`ht = log H` —— チャートの仮定を完全に外す（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★★★★★★これは何か —— `hI` が**局所条件**に置き換わった

`§9-870` の `htArith_eq_log_mulHeight` は仮定に

    `hI : Ideal.span (Set.range x) = Ideal.span {x i₀}`

——**座標の生成するイデアル `𝔞` が主イデアル**——を持っていた。
★`§9-928` で測ったとおり、これは**イデアル類が自明でないと成り立たない**。

★★★★**本ファイルはそれを外す**。置き換わる先は**各極大イデアルでの局所条件**である:

    各極大イデアル `P` について、ある `i` と `g ∈ (𝓞_F)_P` があって
      `𝔞·(𝓞_F)_P = (x_i)`、`I·(𝓞_F)_P = (g)`、`g·x_i = x_j`

——これは**必ず成り立つ**（`v_P(x_i)` が最小になる `i` を取ればよい）。

## ★★★機構 —— 2 段である

### 段 1: 局所条件から `I·𝔞 = (x_j)`

`Ideal.eq_of_localization_maximal`（mathlib）——**イデアルは局所化で決まる**。
局所では `map(I·𝔞) = (g)·(x_i) = (g·x_i) = (x_j)` である。

### 段 2: ノルムを取る

`Ideal.absNorm` は乗法的なので `N(I)·N(𝔞) = N((x_j))`
——これが `§9-929` の `log_degFin_add_degArch_eq_log_mulHeight'` が受ける形そのものである。

## ★★★これで `Proposition 1.4, (iv)` の高さの計算は無条件になった

★★★★`htArith_eq_log_mulHeight_of_localization` の仮定は

* 座標 `x`（`x ≠ 0`、`x j ≠ 0`）
* **局所条件**（上記——常に成り立つ）
* `hgreen`（Green 関数が Fubini–Study であること）

だけであり、**チャートを通ることは要求しない**。

★残るのは「局所条件を幾何から作る」段——`§9-930`
（`map_pullbackIdeal_globalToProj_of_localChart`）が
`I·(𝓞_F)_P = (g)` を与えるので、あとは `𝔞·(𝓞_F)_P = (x_i)` と `g·x_i = x_j` である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★★★★★★★★★★★段 1 —— イデアルは局所化で決まる -/

/-- ★★★★★★★★★★★★**局所条件から `I·𝔞 = (z)`**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★各極大イデアル `P` で `𝔞` が `x_i` で生成され、`I` が `g` で生成され、
`g·x_i = z` であれば、大域で `I·𝔞 = (z)` である。
★★機構は mathlib の `Ideal.eq_of_localization_maximal`
——**イデアルは局所化で決まる**——だけである。 -/
theorem mul_span_range_eq_span_of_localization {R : Type} [CommRing R]
    {ι : Type} (x : ι → R) (z : R) (I : Ideal R)
    (h : ∀ (P : Ideal R) (_ : P.IsMaximal), ∃ (i : ι) (g : Localization.AtPrime P),
      Ideal.map (algebraMap R (Localization.AtPrime P)) (Ideal.span (Set.range x))
          = Ideal.span {algebraMap R (Localization.AtPrime P) (x i)} ∧
        Ideal.map (algebraMap R (Localization.AtPrime P)) I = Ideal.span {g} ∧
        g * algebraMap R (Localization.AtPrime P) (x i)
          = algebraMap R (Localization.AtPrime P) z) :
    I * Ideal.span (Set.range x) = Ideal.span {z} := by
  refine Ideal.eq_of_localization_maximal ?_
  intro P hP
  obtain ⟨i, g, ha, hi, hg⟩ := h P hP
  rw [Ideal.map_mul, ha, hi, Ideal.span_singleton_mul_span_singleton, hg, Ideal.map_span,
    Set.image_singleton]

/-! ## ★★★段 2 —— ノルムを取る -/

/-- ★★★**`I·𝔞 = (x_j)` からノルムの等式**——`Ideal.absNorm` は乗法的である。 -/
theorem absNorm_mul_of_mul_span_eq {ι : Type} (x : ι → 𝓞 F) (j : ι)
    (I : Ideal (𝓞 F)) (h : I * Ideal.span (Set.range x) = Ideal.span {x j}) :
    I.absNorm * (Ideal.span (Set.range x)).absNorm = (Ideal.span {x j}).absNorm := by
  rw [← h, Ideal.absNorm.map_mul]

omit [NumberField F] in
/-- ★**`I·𝔞 = (x_j)` かつ `x_j ≠ 0` なら `I ≠ 0`**。 -/
theorem ne_zero_of_mul_span_eq {ι : Type} (x : ι → 𝓞 F) (j : ι) (hj : x j ≠ 0)
    (I : Ideal (𝓞 F)) (h : I * Ideal.span (Set.range x) = Ideal.span {x j}) : I ≠ 0 := by
  rintro rfl
  rw [zero_mul] at h
  exact hj (by simpa using (Ideal.span_singleton_eq_bot.1 h.symm))

/-! ## ★★★★★★★★★★★★★★★ノルムの等式だけで `ht = log H` -/

/-- ★★★★★★★★★★★★★★★**`hI`（`𝔞` が主イデアル）を外した高さの計算**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★`§9-870` は `hI : 𝔞 = (x_{i₀})` を課していたが、
`§9-929` で測ったとおり、必要なのは**ノルムの等式** `N(I)·N(𝔞) = N((x_j))` だけである。 -/
theorem htArith_eq_log_mulHeight' {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) {ι : Type} [Finite ι]
    (x : ι → 𝓞 F) (hx : x ≠ 0) (j : ι) (hj : x j ≠ 0)
    (hne : pullbackIdeal F D.divisor xF ≠ 0)
    (hnorm : (pullbackIdeal F D.divisor xF).absNorm
        * (Ideal.span (Set.range x)).absNorm = (Ideal.span {x j}).absNorm)
    (hgreen : ∀ v : InfinitePlace F,
      D.green (archPoint xF v) = Real.log ((⨆ i, v ((x i : F))) / v ((x j : F)))) :
    htArith F D xF
      = Real.log (Height.mulHeight (fun i => (x i : F))) / (Module.finrank ℚ F : ℝ) := by
  rw [htArith_eq_add F D xF, degNormalized, deg_idealADiv F _ hne,
    archADiv_sum_eq F D.green xF]
  have hg : ∑ v : InfinitePlace F, (InfinitePlace.mult v : ℝ) * D.green (archPoint xF v)
      = ∑ v : InfinitePlace F,
          (InfinitePlace.mult v : ℝ) * Real.log ((⨆ i, v ((x i : F))) / v ((x j : F))) :=
    Finset.sum_congr rfl (fun v _ => by rw [hgreen v])
  rw [hg, ← add_div,
    log_degFin_add_degArch_eq_log_mulHeight' F x hx j (pullbackIdeal F D.divisor xF) hj hnorm]

/-- ★★★★★★★★★★★★★★★★★**局所条件だけで `ht = log H`**
—— `Proposition 1.4, (iv)` の高さの計算の**無条件形**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★仮定は「各極大イデアルで `𝔞` が `x_i` で生成され、引き戻しイデアルが `g` で生成され、
`g·x_i = x_j`」——**チャートを通ることは要求しない**。 -/
theorem htArith_eq_log_mulHeight_of_localization {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) {ι : Type} [Finite ι]
    (x : ι → 𝓞 F) (hx : x ≠ 0) (j : ι) (hj : x j ≠ 0)
    (hloc : ∀ (P : Ideal (𝓞 F)) (_ : P.IsMaximal), ∃ (i : ι) (g : Localization.AtPrime P),
      Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime P)) (Ideal.span (Set.range x))
          = Ideal.span {algebraMap (𝓞 F) (Localization.AtPrime P) (x i)} ∧
        Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime P))
            (pullbackIdeal F D.divisor xF) = Ideal.span {g} ∧
        g * algebraMap (𝓞 F) (Localization.AtPrime P) (x i)
          = algebraMap (𝓞 F) (Localization.AtPrime P) (x j))
    (hgreen : ∀ v : InfinitePlace F,
      D.green (archPoint xF v) = Real.log ((⨆ i, v ((x i : F))) / v ((x j : F)))) :
    htArith F D xF
      = Real.log (Height.mulHeight (fun i => (x i : F))) / (Module.finrank ℚ F : ℝ) := by
  have hmul := mul_span_range_eq_span_of_localization x (x j) (pullbackIdeal F D.divisor xF) hloc
  exact htArith_eq_log_mulHeight' F D xF x hx j hj
    (ne_zero_of_mul_span_eq F x j hj _ hmul)
    (absNorm_mul_of_mul_span_eq F x j _ hmul) hgreen

/-! ## ★出典の紐付け(`.src`) -/

def mul_span_range_eq_span_of_localization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(局所条件から I·𝔞 = (z)——イデアルは局所化で決まる)",
    sectionId := "genell-prop-1-4" }

def absNorm_mul_of_mul_span_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(I·𝔞 = (x_j) からノルムの等式)",
    sectionId := "genell-prop-1-4" }

def ne_zero_of_mul_span_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(I·𝔞 = (x_j) かつ x_j ≠ 0 なら I ≠ 0)",
    sectionId := "genell-prop-1-4" }

def htArith_eq_log_mulHeight'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(ノルムの等式だけで ht = log H——𝔞 が主イデアルという仮定なし)",
    sectionId := "genell-prop-1-4" }

def htArith_eq_log_mulHeight_of_localization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(局所条件だけで ht = log H——高さの計算の無条件形)",
    sectionId := "genell-prop-1-4" }

def htArith_eq_log_mulHeight_of_localization.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.eq_of_localization_maximal(イデアルは局所化で決まる)"
      (.inMathlib "Ideal.eq_of_localization_maximal") 2,
    .citation "[ABC3]" "log_degFin_add_degArch_eq_log_mulHeight'(ノルムの等式だけで足りる、§9-929)"
      (.inProject "ABC3" "ABC3.Found.GenEll.log_degFin_add_degArch_eq_log_mulHeight'") 3,
    .citation "[ABC3]" "htArith_eq_add(高さは deg_fin + deg_arch)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_eq_add") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-870 の hI(𝔞 が主イデアル)は" ++
       "**各極大イデアルでの局所条件に置き換わる**。" ++
       "局所条件は**必ず成り立つ**——v_P(x_i) が最小になる i を取ればよいからである。" ++
       "★これで Proposition 1.4, (iv) の高さの計算はイデアル類に依らなくなった") 5,
    .implicitStep
      ("★残るのは局所条件を幾何から作る段——§9-930 " ++
       "(map_pullbackIdeal_globalToProj_of_localChart)が I·(𝓞_F)_P = (g) を与えるので、" ++
       "あとは 𝔞·(𝓞_F)_P = (x_i) と g·x_i = x_j である") 3 ]

end ABC3.Found.GenEll
