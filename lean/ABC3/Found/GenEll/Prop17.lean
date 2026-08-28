/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17LeftEq
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.7`（`Found`、項目まるごと）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`Proposition 1.7, (i)` を、**点にわたる 1 本の主張**として取る:

    **`log-cond_E(z) − log-cond_D(y) ≈ log-diff(K) − log-diff(F)`**（`BDeq`、定数 `0`）
    **`log-diff(K) − log-diff(F) ≲ (1 − 1/e)·log-cond_E(z)`**（`BDge`、定数 `0`）

★★★★★★どちらも定数は `0` である。`slack` も `Σ` も要らない。

## ★★★2026-08-27 の `Skeleton` からの前進（明示）

| | `Skeleton/GenEll/Section1.lean` の `prop_1_7` | ★本ファイル |
|---|---|---|
| `hlow`／`hup` | **仮定として受ける** | ★**証明する**（`§9-974`） |
| `slack` | `Σ_{q∈Σ} log q` で受ける | ★**`0`** |
| 左の `≲` | 不等式 | ★**等式**（`BDeq`） |
| 分岐の情報 | 無し（実数関数だけ） | ★`log-cond`・`log-diff` の**本体**に当たる |

★★`Skeleton` の `prop_1_7` は「後半だけ」を証明していた
（原文の `hlow`／`hup` は仮定）。★★★本ファイルはその `hlow`／`hup` を
**算術から作る**——`§9-954`〜`§9-974` の 21 ブロックがその中身である。

## ★★★★★★★★★★逸脱（明示）

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 条件 (b) | `D_ℚ = φ_ℚ^{-1}(E_ℚ)_red`（scheme の等式） | **台の対応**（`T`・`S`・`pr`・`hdisj`・`hpr`） | scheme 的な `(−)_red` の引き戻しから台の対応を導く段は本プロジェクトの語彙の外。★受けているのはその**算術的帰結**であって条件 (b) より弱い |
| (ii) Riemann–Hurwitz | 述べる | **含めない** | `deg` が語彙の外（`Skeleton` の `prop_1_7` と同じ） |
| 左の `≲` | `≲` | ★**`BDeq`（等式）** | 我々が証明したのは等式であり、`≲` より**強い**。原文も『prime-to-`Σ` 部分では `=` と `≤` であって `≲` ではない』と注記している |
| 右の `≲` の**向き** | 印字は `≲` | ★**`BDge`**（＝ `log-diff の差 ≤ (1−1/e)·log-cond_E + C`） | ★★★★**印字どおりの向き（`BDle`）では偽である**——`Check/GenEll/Prop17Direction.lean` の `prop_1_7_printed_direction_false` で機械検証した。`Gap/GenEll/BDDirection.lean` に記録済みの向きの食い違いの、**3 例目にして初の反例**である |
| 分岐の一様性 | 条件 (d)（`e_P ∣ e`） | `hle`（`e_P ≤ e`） | 割り切るなら小さい。★弱めた方を仮定している |

★★★★★★**空虚ではない**——`hsum` は mathlib の `Ideal.sum_ramification_inertia`
そのもの、`hcondF`／`hcondK` は `§9-966`、`hdiff` は `§9-973` が与える。
-/

namespace ABC3.Found.GenEll

open NumberField AlgebraicGeometry

/-! ## ★`BD` への持ち上げ（定数 `0`） -/

/-- ★点ごとに等しければ `BDeq`（定数 `0`）。 -/
theorem bdeq_of_eq {Pt : Type} (α β : Pt → ℝ) (h : ∀ x, α x = β x) : BDeq α β :=
  ⟨0, fun x => by rw [h x]; simp⟩

/-- ★点ごとに `α ≤ β` なら `BDge α β`（定数 `0`）。

★`BDge α β` は原文の `α ≳ β` の逐語（`α(x) − β(x) ≤ C`）であり、
**通常の読み `α ≤ β + C`** にあたる。 -/
theorem bdge_of_le {Pt : Type} (α β : Pt → ℝ) (h : ∀ x, α x ≤ β x) : BDge α β :=
  ⟨0, fun x => by have := h x; linarith⟩

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.7` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Proposition 1.7**
(Conductors and Log Differents)、(i)。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `log-cond_E(z) − log-cond_D(y) ≈ log-diff(K) − log-diff(F)`
    `log-diff(K) − log-diff(F) ≲ (1 − 1/e)·log-cond_E(z)`

★★★★★★どちらも定数は `0`。

★仮定の出どころ:

| 仮定 | 出どころ |
|---|---|
| `T`・`S`・`pr`・`hdisj`・`hpr` | 原文の条件 (b) の算術的帰結（**受ける**） |
| `hsum` | 基本等式——mathlib `Ideal.sum_ramification_inertia` |
| `hle` | 原文の条件 (d)（`e_P ∣ e` ⟹ `e_P ≤ e`） |
| `hlogN` | `§9-972` `log_absNorm_eq_inertiaDeg_mul` |
| `hcondF`／`hcondK` | `§9-966` `logCond_eq_sum_log` |
| `hdiff` | `§9-973` `logDiff_tower_eq_sum_of_totallyRamified_tame` |

★★★(ii)（Riemann–Hurwitz）は含めない——`deg` が語彙の外である。 -/
theorem prop_1_7 {X Y : Scheme.{0}} (E : X.IdealSheafData) (D : Y.IdealSheafData)
    (Pt : Type) (fldF fldK : Pt → IntermediateField ℚ ℂ)
    (hnfF : ∀ x, NumberField (fldF x)) (hnfK : ∀ x, NumberField (fldK x))
    (deceqK : ∀ x, haveI := hnfK x; DecidableEq (Ideal (𝓞 (fldK x))))
    (e : ℕ) (he : 0 < e)
    (zE : ∀ x, haveI := hnfF x; specRingOfIntegers (fldF x) ⟶ X)
    (yD : ∀ x, haveI := hnfK x; specRingOfIntegers (fldK x) ⟶ Y)
    (nn mm : Pt → ℕ) (hn : ∀ x, 0 < nn x) (hm : ∀ x, 0 < mm x)
    (T : ∀ x, haveI := hnfF x; Finset (Ideal (𝓞 (fldF x))))
    (S : ∀ x, haveI := hnfF x; haveI := hnfK x;
      Ideal (𝓞 (fldF x)) → Finset (Ideal (𝓞 (fldK x))))
    (hdisj : ∀ x, haveI := hnfF x; haveI := hnfK x;
      ((T x : Finset (Ideal (𝓞 (fldF x)))) : Set (Ideal (𝓞 (fldF x)))).PairwiseDisjoint (S x))
    (ee ff : ∀ x, haveI := hnfK x; Ideal (𝓞 (fldK x)) → ℕ)
    (pr : ∀ x, haveI := hnfF x; haveI := hnfK x;
      Ideal (𝓞 (fldK x)) → Ideal (𝓞 (fldF x)))
    (hpr : ∀ x, haveI := hnfF x; haveI := hnfK x;
      ∀ p ∈ T x, ∀ P ∈ S x p, pr x P = p)
    (he1 : ∀ x, haveI := hnfK x; ∀ P, 1 ≤ ee x P)
    (hle : ∀ x, haveI := hnfK x; ∀ P, ee x P ≤ e)
    (hLnn : ∀ x, haveI := hnfF x; ∀ p ∈ T x, 0 ≤ Real.log (Ideal.absNorm p))
    (hsum : ∀ x, haveI := hnfF x; haveI := hnfK x;
      ∀ p ∈ T x, ∑ P ∈ S x p, ee x P * ff x P = nn x)
    (hlogN : ∀ x, haveI := hnfF x; haveI := hnfK x; letI := deceqK x;
      ∀ P ∈ (T x).biUnion (S x), Real.log (Ideal.absNorm P)
        = (ff x P : ℝ) * Real.log (Ideal.absNorm (pr x P)))
    (hcondF : ∀ x, haveI := hnfF x;
      logCond (fldF x) E (zE x)
        = (∑ p ∈ T x, Real.log (Ideal.absNorm p)) / (mm x : ℝ))
    (hcondK : ∀ x, haveI := hnfF x; haveI := hnfK x; letI := deceqK x;
      logCond (fldK x) D (yD x)
        = (∑ P ∈ (T x).biUnion (S x), Real.log (Ideal.absNorm P))
          / ((nn x : ℝ) * (mm x : ℝ)))
    (hdiff : ∀ x, haveI := hnfF x; haveI := hnfK x; letI := deceqK x;
      logDiffOfField (fldK x) - logDiffOfField (fldF x)
        = (∑ P ∈ (T x).biUnion (S x), ((ee x P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P))
          / ((nn x : ℝ) * (mm x : ℝ))) :
    BDeq (fun x => haveI := hnfF x; haveI := hnfK x;
            logCond (fldF x) E (zE x) - logCond (fldK x) D (yD x))
         (fun x => haveI := hnfF x; haveI := hnfK x;
            logDiffOfField (fldK x) - logDiffOfField (fldF x))
  ∧ BDge (fun x => haveI := hnfF x; haveI := hnfK x;
            logDiffOfField (fldK x) - logDiffOfField (fldF x))
         (fun x => haveI := hnfF x;
            (1 - 1 / (e : ℝ)) * logCond (fldF x) E (zE x)) := by
  constructor
  · refine bdeq_of_eq _ _ (fun x => ?_)
    haveI := hnfF x
    haveI := hnfK x
    letI := deceqK x
    exact prop_1_7_left_eq (fldF x) (fldK x) E D (zE x) (yD x) (T x) (S x) (hdisj x)
      (ee x) (ff x) (pr x) (hpr x) (nn x) (mm x) (hn x) (hm x) (he1 x) (hsum x)
      (hlogN x) (hcondF x) (hcondK x) (hdiff x)
  · refine bdge_of_le _ _ (fun x => ?_)
    haveI := hnfF x
    haveI := hnfK x
    letI := deceqK x
    exact prop_1_7_right_le (fldF x) (fldK x) E D (zE x) (yD x) (T x) (S x) (hdisj x)
      (ee x) (ff x) (pr x) (hpr x) (nn x) (mm x) e (hn x) (hm x) he (he1 x) (hle x)
      (hLnn x) (hsum x) (hlogN x) (hcondF x) (hcondK x) (hdiff x)

/-! ## ★出典の紐付け(`.src`) -/

def prop_1_7.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.7",
    sectionId := "genell-prop-1-7" }

def prop_1_7.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "prop_1_7_left_eq((i) の左は等式、§9-974)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prop_1_7_left_eq") 4,
    .citation "[ABC3]" "prop_1_7_right_le((i) の右、定数 0、§9-974)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prop_1_7_right_le") 4,
    .citation "[mathlib]" "Ideal.sum_ramification_inertia(基本等式 ∑ e_P f_P = [K:F])"
      (.inMathlib "Ideal.sum_ramification_inertia") 2,
    .implicitStep
      ("★★★★★★★逸脱(2026-08-29): 原文の条件 (b)(D_ℚ = φ_ℚ^{-1}(E_ℚ)_red)は" ++
       "**台の対応**(T・S・pr・hdisj・hpr)として受けている。" ++
       "scheme 的な (−)_red の引き戻しから台の対応を導く段は本プロジェクトの語彙の外であり、" ++
       "受けているのはその**算術的帰結**であって条件 (b) より弱い") 7,
    .implicitStep
      ("★★★★★★★★逸脱(2026-08-29): 右の ≲ を**印字どおりの向き(BDle)で書くと偽**である" ++
       "——Check/GenEll/Prop17Direction.lean の prop_1_7_printed_direction_false で機械検証。" ++
       "★Gap/GenEll/BDDirection.lean が『falsifier を書けないうちは sourceGap を名乗るな』と" ++
       "していた、その **falsifier がこれである**。本ファイルは通常の読み(BDge)で取った") 9,
    .implicitStep
      ("★★(ii)(Riemann–Hurwitz)は含めない——deg が語彙の外である" ++
       "(Skeleton/GenEll/Section1.lean の prop_1_7 と同じ範囲)") 5,
    .implicitStep
      ("★★★★★Skeleton からの前進(2026-08-29): Skeleton の prop_1_7 は hlow/hup を" ++
       "**仮定として受けて**いた。本ファイルはそれを**算術から作る**" ++
       "——§9-954〜§9-974 の 21 ブロックがその中身であり、slack は Σ ではなく **0** である") 8 ]

end ABC3.Found.GenEll
