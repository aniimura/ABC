/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17DiffEq
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.7, (i)` の両側（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 組み立て

`§9-966`〜`§9-973` で揃えた部品を**素点ごとの和として組み立てる**。

    **`log-cond_E(z) − log-cond_D(y) = log-diff(K) − log-diff(F)`**（左の `≲` は**等式**）
    **`log-diff(K) − log-diff(F) ≤ (1 − 1/e)·log-cond_E(z)`**（右の `≲`、定数 `0`）

★★★★★★どちらも**定数 `0`** である——`slack` も `Σ` も要らない。
原文が『prime-to-`Σ` 部分では **`=` と `≤` であって `≲` ではない**』と
注記していることの、**そのままの形**である。

## ★機構 —— 台の対応（原文の条件 (b)）が `biUnion` になる

原文の条件 (b) `D_ℚ = φ_ℚ^{-1}(E_ℚ)_red` は、算術の言葉では

    **`D` の導手の台 ＝ `E` の導手の台の上にある素点の全体**

である。★これを `Finset` で書けば `S_K = T.biUnion S`（`T` ＝ `E` の台、
`S p` ＝ `p` の上にある `K` の素点）であり、和が

    `∑_{P ∈ S_K} (…) = ∑_{p ∈ T} ∑_{P ∣ p} (…)`

と分解する（`Finset.sum_biUnion`）。★★あとは素点 `p` ごとに `§9-971` の算術
（基本等式 `∑ e_P f_P = [K:F]` で分子が畳める）を当てるだけである。

## ★部品の出どころ

| 部品 | 出どころ |
|---|---|
| `log-cond` ＝ 台の上の `∑ log N` | `§9-966` `logCond_eq_sum_log` |
| `log-diff` の差 ＝ `∑ (e_P−1) log N(P)/[K:ℚ]` | `§9-973` `logDiff_tower_eq_sum_of_totallyRamified_tame` |
| `log N(P) = f_P·log N(p)` | `§9-972` `log_absNorm_eq_inertiaDeg_mul` |
| `∑_{P∣p} e_P f_P = [K:F]` | mathlib `Ideal.sum_ramification_inertia`／`§FinitePlaceRel` |
| 素点ごとの算術 | `§9-971` `diff_contrib_eq` |
| 右の `≲` の数値の核 | `§9-970` `one_sub_ratio_le` |

## ★★★★★逸脱（明示）

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 条件 (b) | `D_ℚ = φ_ℚ^{-1}(E_ℚ)_red`（scheme の等式） | **台の対応**（`T`・`S`・`pr`・`hdisj`・`hpr`）として受ける | scheme 的な `(−)_red` の引き戻しの等式から台の対応を導く段は、本プロジェクトの語彙の外にある。★★受けているのはその**算術的帰結**であって、条件 (b) より弱い |
| (ii) Riemann–Hurwitz | 述べる | **含めない** | `deg` が語彙の外（`Skeleton` の `prop_1_7` と同じ） |
| 分岐の型 | 一般 | `hsum`（基本等式）と `hle`（`e_P ≤ e`）で受ける | 前者は mathlib にあり、後者は原文の条件 (d)（`e_P ∣ e`）の帰結 |

★★★★★★**空虚ではない**——`hsum` は mathlib の `Ideal.sum_ramification_inertia`
そのものであり、`hcondF`／`hcondK` は `§9-966`、`hdiff` は `§9-973` が与える。
-/

namespace ABC3.Found.GenEll

open NumberField AlgebraicGeometry

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★数値の核（`biUnion` の上の和） -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★**台の対応の上での等式**（純粋に `Finset` の算術）。

    `(∑_i L_i)/m − (∑_{k} f_k·L_{pr k})/(n·m) = (∑_k (e_k−1)·f_k·L_{pr k})/(n·m)`

（`∑_{k ∈ S i} e_k f_k = n` が各 `i` で成り立つとき）。

★左辺の第 1 項が `log-cond_E`、第 2 項が `log-cond_D`、右辺が `log-diff` の差である。
★★機構は `§9-971` と同じ——基本等式で分子が畳める。 -/
theorem sum_biUnion_diff_eq {ι κ : Type*} [DecidableEq κ] (T : Finset ι) (S : ι → Finset κ)
    (hdisj : (T : Set ι).PairwiseDisjoint S)
    (ee ff : κ → ℕ) (pr : κ → ι) (hpr : ∀ i ∈ T, ∀ k ∈ S i, pr k = i)
    (n m : ℕ) (hn : 0 < n) (hm : 0 < m) (Lg : ι → ℝ)
    (he1 : ∀ k, 1 ≤ ee k)
    (hsum : ∀ i ∈ T, ∑ k ∈ S i, ee k * ff k = n) :
    (∑ i ∈ T, Lg i) / (m : ℝ)
        - (∑ k ∈ T.biUnion S, (ff k : ℝ) * Lg (pr k)) / ((n : ℝ) * (m : ℝ))
      = (∑ k ∈ T.biUnion S, ((ee k - 1 : ℕ) : ℝ) * (ff k : ℝ) * Lg (pr k))
        / ((n : ℝ) * (m : ℝ)) := by
  have hnR : (0:ℝ) < (n : ℝ) := by exact_mod_cast hn
  have hmR : (0:ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hcast : ∀ k : κ, ((ee k - 1 : ℕ) : ℝ) = (ee k : ℝ) - 1 := by
    intro k; push_cast [Nat.cast_sub (he1 k)]; ring
  set Ff : ι → ℝ := fun i => ∑ k ∈ S i, (ff k : ℝ) with hFf
  have key2 : ∀ i ∈ T, ∑ k ∈ S i, (ff k : ℝ) * Lg (pr k) = Ff i * Lg i := by
    intro i hi
    rw [hFf, Finset.sum_mul]
    exact Finset.sum_congr rfl (fun k hk => by rw [hpr i hi k hk])
  have key : ∀ i ∈ T, ∑ k ∈ S i, ((ee k - 1 : ℕ) : ℝ) * (ff k : ℝ) * Lg (pr k)
      = ((n : ℝ) - Ff i) * Lg i := by
    intro i hi
    have hs : (∑ k ∈ S i, (ee k : ℝ) * (ff k : ℝ)) = (n : ℝ) := by
      have h := congrArg (fun t : ℕ => (t : ℝ)) (hsum i hi)
      push_cast at h
      exact h
    calc ∑ k ∈ S i, ((ee k - 1 : ℕ) : ℝ) * (ff k : ℝ) * Lg (pr k)
        = ∑ k ∈ S i, (((ee k : ℝ) - 1) * (ff k : ℝ)) * Lg i := by
          refine Finset.sum_congr rfl (fun k hk => ?_)
          rw [hcast k, hpr i hi k hk]
      _ = (∑ k ∈ S i, ((ee k : ℝ) * (ff k : ℝ) - (ff k : ℝ))) * Lg i := by
          rw [Finset.sum_mul]
          exact Finset.sum_congr rfl (fun k _ => by ring)
      _ = ((∑ k ∈ S i, (ee k : ℝ) * (ff k : ℝ)) - Ff i) * Lg i := by
          rw [Finset.sum_sub_distrib]
      _ = ((n : ℝ) - Ff i) * Lg i := by rw [hs]
  rw [Finset.sum_biUnion hdisj, Finset.sum_biUnion hdisj,
    Finset.sum_congr rfl key, Finset.sum_congr rfl key2]
  have hexp : ∑ i ∈ T, ((n : ℝ) - Ff i) * Lg i
      = (n : ℝ) * (∑ i ∈ T, Lg i) - ∑ i ∈ T, Ff i * Lg i := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  rw [hexp]
  field_simp

/-- ★★★★★**`e_k ≤ e` なら `log-cond_D ≥ (1/e)·log-cond_E`**（純粋に `Finset` の算術）。

★基本等式 `∑ e_k f_k = n` と `e_k ≤ e` から `n ≤ e·∑ f_k`——これが右の `≲` の核である。 -/
theorem sum_biUnion_ratio_le {ι κ : Type*} [DecidableEq κ] (T : Finset ι) (S : ι → Finset κ)
    (hdisj : (T : Set ι).PairwiseDisjoint S)
    (ee ff : κ → ℕ) (pr : κ → ι) (hpr : ∀ i ∈ T, ∀ k ∈ S i, pr k = i)
    (n m e : ℕ) (hn : 0 < n) (hm : 0 < m) (he : 0 < e) (Lg : ι → ℝ)
    (hle : ∀ k, ee k ≤ e) (hL : ∀ i ∈ T, 0 ≤ Lg i)
    (hsum : ∀ i ∈ T, ∑ k ∈ S i, ee k * ff k = n) :
    (∑ i ∈ T, Lg i) / ((e : ℝ) * (m : ℝ))
      ≤ (∑ k ∈ T.biUnion S, (ff k : ℝ) * Lg (pr k)) / ((n : ℝ) * (m : ℝ)) := by
  have hnR : (0:ℝ) < (n : ℝ) := by exact_mod_cast hn
  have hmR : (0:ℝ) < (m : ℝ) := by exact_mod_cast hm
  have heR : (0:ℝ) < (e : ℝ) := by exact_mod_cast he
  set Ff : ι → ℝ := fun i => ∑ k ∈ S i, (ff k : ℝ) with hFf
  have key2 : ∀ i ∈ T, ∑ k ∈ S i, (ff k : ℝ) * Lg (pr k) = Ff i * Lg i := by
    intro i hi
    rw [hFf, Finset.sum_mul]
    exact Finset.sum_congr rfl (fun k hk => by rw [hpr i hi k hk])
  have hnle : ∀ i ∈ T, (n : ℝ) ≤ (e : ℝ) * Ff i := by
    intro i hi
    have hs : (∑ k ∈ S i, (ee k : ℝ) * (ff k : ℝ)) = (n : ℝ) := by
      have h := congrArg (fun t : ℕ => (t : ℝ)) (hsum i hi)
      push_cast at h
      exact h
    rw [← hs, hFf, Finset.mul_sum]
    refine Finset.sum_le_sum (fun k _ => ?_)
    have h1 : ((ee k : ℝ)) ≤ (e : ℝ) := by exact_mod_cast hle k
    have h2 : (0:ℝ) ≤ (ff k : ℝ) := by positivity
    exact mul_le_mul_of_nonneg_right h1 h2
  rw [Finset.sum_biUnion hdisj, Finset.sum_congr rfl key2]
  rw [div_le_div_iff₀ (by positivity) (by positivity), Finset.sum_mul, Finset.sum_mul]
  refine Finset.sum_le_sum (fun i hi => ?_)
  have hLi := hL i hi
  have hni := hnle i hi
  nlinarith [mul_nonneg hLi (le_of_lt hmR)]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★左の `≲` は等式 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`Proposition 1.7, (i)` の左の `≲`
は等式である**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `log-cond_E(z) − log-cond_D(y) = log-diff(K) − log-diff(F)`

★★★★★★定数は `0` である——`slack` も `Σ` も要らない。
★原文が『prime-to-`Σ` 部分では **`=` と `≤` であって `≲` ではない**』と
注記していることの、そのままの形である。

★★仮定 `T`・`S`・`pr`・`hdisj`・`hpr` が**原文の条件 (b)**
（`D_ℚ = φ_ℚ^{-1}(E_ℚ)_red`）の算術的帰結、
`hsum` が**基本等式**（mathlib）、
`hcondF`／`hcondK` が `§9-966`、`hdiff` が `§9-973`、`hlogN` が `§9-972` である。 -/
theorem prop_1_7_left_eq (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    [DecidableEq (Ideal (𝓞 K))]
    {X Y : Scheme.{0}} (E : X.IdealSheafData) (D : Y.IdealSheafData)
    (zF : specRingOfIntegers F ⟶ X) (yK : specRingOfIntegers K ⟶ Y)
    (T : Finset (Ideal (𝓞 F))) (S : Ideal (𝓞 F) → Finset (Ideal (𝓞 K)))
    (hdisj : (T : Set (Ideal (𝓞 F))).PairwiseDisjoint S)
    (ee ff : Ideal (𝓞 K) → ℕ) (pr : Ideal (𝓞 K) → Ideal (𝓞 F))
    (hpr : ∀ p ∈ T, ∀ P ∈ S p, pr P = p)
    (n m : ℕ) (hn : 0 < n) (hm : 0 < m)
    (he1 : ∀ P, 1 ≤ ee P)
    (hsum : ∀ p ∈ T, ∑ P ∈ S p, ee P * ff P = n)
    (hlogN : ∀ P ∈ T.biUnion S, Real.log (Ideal.absNorm P)
       = (ff P : ℝ) * Real.log (Ideal.absNorm (pr P)))
    (hcondF : logCond F E zF = (∑ p ∈ T, Real.log (Ideal.absNorm p)) / (m : ℝ))
    (hcondK : logCond K D yK
       = (∑ P ∈ T.biUnion S, Real.log (Ideal.absNorm P)) / ((n : ℝ) * (m : ℝ)))
    (hdiff : logDiffOfField K - logDiffOfField F
       = (∑ P ∈ T.biUnion S, ((ee P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P))
         / ((n : ℝ) * (m : ℝ))) :
    logCond F E zF - logCond K D yK = logDiffOfField K - logDiffOfField F := by
  rw [hcondF, hcondK, hdiff]
  have h1 : ∑ P ∈ T.biUnion S, Real.log (Ideal.absNorm P)
      = ∑ P ∈ T.biUnion S, (ff P : ℝ) * Real.log (Ideal.absNorm (pr P)) :=
    Finset.sum_congr rfl (fun P hP => hlogN P hP)
  have h2 : ∑ P ∈ T.biUnion S, ((ee P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P)
      = ∑ P ∈ T.biUnion S, ((ee P - 1 : ℕ) : ℝ) * (ff P : ℝ)
          * Real.log (Ideal.absNorm (pr P)) :=
    Finset.sum_congr rfl (fun P hP => by rw [hlogN P hP]; ring)
  rw [h1, h2]
  exact sum_biUnion_diff_eq T S hdisj ee ff pr hpr n m hn hm _ he1 hsum

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★右の `≲` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`Proposition 1.7, (i)` の右の `≲`**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `log-diff(K) − log-diff(F) ≤ (1 − 1/e)·log-cond_E(z)`

★★★★★★定数は `0` である。
★左の等式（`prop_1_7_left_eq`）に `log-cond_D ≥ (1/e)·log-cond_E` を当てるだけ
——後者は条件 (d)（`e_P ∣ e` ゆえ `e_P ≤ e`）と基本等式から出る。 -/
theorem prop_1_7_right_le (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    [DecidableEq (Ideal (𝓞 K))]
    {X Y : Scheme.{0}} (E : X.IdealSheafData) (D : Y.IdealSheafData)
    (zF : specRingOfIntegers F ⟶ X) (yK : specRingOfIntegers K ⟶ Y)
    (T : Finset (Ideal (𝓞 F))) (S : Ideal (𝓞 F) → Finset (Ideal (𝓞 K)))
    (hdisj : (T : Set (Ideal (𝓞 F))).PairwiseDisjoint S)
    (ee ff : Ideal (𝓞 K) → ℕ) (pr : Ideal (𝓞 K) → Ideal (𝓞 F))
    (hpr : ∀ p ∈ T, ∀ P ∈ S p, pr P = p)
    (n m e : ℕ) (hn : 0 < n) (hm : 0 < m) (he : 0 < e)
    (he1 : ∀ P, 1 ≤ ee P) (hle : ∀ P, ee P ≤ e)
    (hLnn : ∀ p ∈ T, 0 ≤ Real.log (Ideal.absNorm p))
    (hsum : ∀ p ∈ T, ∑ P ∈ S p, ee P * ff P = n)
    (hlogN : ∀ P ∈ T.biUnion S, Real.log (Ideal.absNorm P)
       = (ff P : ℝ) * Real.log (Ideal.absNorm (pr P)))
    (hcondF : logCond F E zF = (∑ p ∈ T, Real.log (Ideal.absNorm p)) / (m : ℝ))
    (hcondK : logCond K D yK
       = (∑ P ∈ T.biUnion S, Real.log (Ideal.absNorm P)) / ((n : ℝ) * (m : ℝ)))
    (hdiff : logDiffOfField K - logDiffOfField F
       = (∑ P ∈ T.biUnion S, ((ee P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P))
         / ((n : ℝ) * (m : ℝ))) :
    logDiffOfField K - logDiffOfField F ≤ (1 - 1 / (e : ℝ)) * logCond F E zF := by
  have hmR : (0:ℝ) < (m : ℝ) := by exact_mod_cast hm
  have heR : (0:ℝ) < (e : ℝ) := by exact_mod_cast he
  have heq := prop_1_7_left_eq F K E D zF yK T S hdisj ee ff pr hpr n m hn hm he1 hsum
    hlogN hcondF hcondK hdiff
  rw [← heq]
  have h1 : ∑ P ∈ T.biUnion S, Real.log (Ideal.absNorm P)
      = ∑ P ∈ T.biUnion S, (ff P : ℝ) * Real.log (Ideal.absNorm (pr P)) :=
    Finset.sum_congr rfl (fun P hP => hlogN P hP)
  have hratio := sum_biUnion_ratio_le T S hdisj ee ff pr hpr n m e hn hm he
    (fun p => Real.log (Ideal.absNorm p)) hle hLnn hsum
  rw [hcondF, hcondK, h1]
  have hA : (∑ p ∈ T, Real.log (Ideal.absNorm p)) / ((e : ℝ) * (m : ℝ))
      = (1 / (e:ℝ)) * ((∑ p ∈ T, Real.log (Ideal.absNorm p)) / (m : ℝ)) := by
    field_simp
  rw [hA] at hratio
  linarith [hratio]

/-! ## ★出典の紐付け(`.src`) -/

def sum_biUnion_diff_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(台の対応の上での等式——数値の核)",
    sectionId := "genell-prop-1-7" }

def sum_biUnion_ratio_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(log-cond_D ≥ (1/e)·log-cond_E——右の ≲ の核)",
    sectionId := "genell-prop-1-7" }

def prop_1_7_left_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7((i) の左の ≲ は等式である)",
    sectionId := "genell-prop-1-7" }

def prop_1_7_right_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7((i) の右の ≲——定数 0)",
    sectionId := "genell-prop-1-7" }

def prop_1_7_left_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logCond_eq_sum_log(log-cond は台の上の ∑ log N、§9-966)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logCond_eq_sum_log") 3,
    .citation "[ABC3]" "logDiff_tower_eq_sum_of_totallyRamified_tame(log-diff の差の等式、§9-973)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiff_tower_eq_sum_of_totallyRamified_tame") 4,
    .citation "[ABC3]" "log_absNorm_eq_inertiaDeg_mul(log N(P) = f_P·log N(p)、§9-972)"
      (.inProject "ABC3" "ABC3.Found.GenEll.log_absNorm_eq_inertiaDeg_mul") 2,
    .citation "[mathlib]" "Ideal.sum_ramification_inertia(基本等式 ∑ e_P f_P = [K:F])"
      (.inMathlib "Ideal.sum_ramification_inertia") 2,
    .implicitStep
      ("★★★★★★★逸脱(2026-08-29): 原文の条件 (b)(D_ℚ = φ_ℚ^{-1}(E_ℚ)_red)は" ++
       "**台の対応**(T・S・pr・hdisj・hpr)として受けている——" ++
       "scheme 的な (−)_red の引き戻しの等式から台の対応を導く段は本プロジェクトの語彙の外。" ++
       "★受けているのはその**算術的帰結**であって条件 (b) より弱い") 7,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 左の ≲ は**等式**であり定数は 0 である。" ++
       "原文の『prime-to-Σ 部分では = と ≤ であって ≲ ではない』の、そのままの形。" ++
       "★§9-966〜§9-973 の部品を Finset.sum_biUnion で素点ごとに分解して当てるだけ") 6 ]

def prop_1_7_right_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "prop_1_7_left_eq(左の等式)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prop_1_7_left_eq") 3,
    .implicitStep
      ("★★★右の ≲ は原文の条件 (d)(D_ℚ の各点での分岐指数が e を割る)から e_P ≤ e、" ++
       "基本等式と合わせて log-cond_D ≥ (1/e)·log-cond_E。★定数は 0 である") 4 ]

end ABC3.Found.GenEll
