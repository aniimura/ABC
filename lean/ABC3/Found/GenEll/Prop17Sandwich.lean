/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentTameGlobal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★`log-diff` が `log-cond` で両側から挟まれた（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 挟み込み

`§9-967` は `hup`（右側の `≲`）を原文の `log-cond` で書いた。
`§9-968` は `hlow`（左側の `≲`）を大域の定理にした。
★★★★★★本ファイルは後者も **`log-cond` の言葉**にして、**両側から挟む**:

    `(1 − 1/e)·log-cond_E(x)  ≲  log-diff(K_x) − log-diff(F_x)  ≲  (e − 1)·log-cond_E(x)`

★★どちらの `≲` も**定数 `0`** である（`slack` も `Σ` も要らない）。

## ★★★機構

* 上から: `§9-968`（完全分岐・馴なら `𝔡` の指数は `e_P − 1` 以下）＋
  `e_P ≤ e` から `(e_P − 1) ≤ (e − 1)`
* 導手への読み替え: `§9-966`（`log-cond` は台の上の `∑ log N(P)`）

## ★これが原文の何にあたるか

原文 `Proposition 1.7, (i)` は

    `log-cond_E − log-cond_D  ≲  log-diff_Y − log-diff_Z  ≲  (1 − 1/e)·log-cond_E`

である。★★本ファイルの挟み込みは**同じ形の両側評価**であり、
係数と `log-cond_D` の項が違う:

| | 原文 | 本ファイル |
|---|---|---|
| 下から | `log-cond_E − log-cond_D` | `(1 − 1/e)·log-cond_E` |
| 上から | `(1 − 1/e)·log-cond_E` | `(e − 1)·log-cond_E` |

★★★**`log-cond_D`（底の側の導手）の項が本プロジェクトにまだ無い**
——それが `Proposition 1.7` に残る最後の差である。
★分岐の側からは `(1 − 1/e)` と `(e − 1)` の 2 つの係数しか出ず、
原文の形にするには底の側の導手が要る。
-/

namespace ABC3.Found.GenEll

open NumberField AlgebraicGeometry

/-! ## ★★★★係数の比較 -/

/-- ★★★★**`e_P ≤ e` なら `∑ (e_P − 1)·log N(P) ≤ (e − 1)·∑ log N(P)`**。 -/
theorem sum_sub_one_le_mul (K : Type) [Field K] [NumberField K]
    (e : ℕ) (he : 0 < e) (S : Finset (Ideal (𝓞 K))) (hne : ∀ P ∈ S, P ≠ ⊥)
    (ram : Ideal (𝓞 K) → ℕ) (hle : ∀ P ∈ S, ram P ≤ e) :
    (∑ P ∈ S, ((ram P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P))
      ≤ ((e : ℝ) - 1) * ∑ P ∈ S, Real.log (Ideal.absNorm P) := by
  have hcast : ((e - 1 : ℕ) : ℝ) = (e : ℝ) - 1 := by
    have h1 : 1 ≤ e := he
    push_cast [Nat.cast_sub h1]
    ring
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum (fun P hP => ?_)
  have hNpos : 1 ≤ (Ideal.absNorm P : ℝ) := by
    have h : Ideal.absNorm P ≠ 0 := by rw [Ne, Ideal.absNorm_eq_zero_iff]; exact hne P hP
    exact_mod_cast Nat.one_le_iff_ne_zero.mpr h
  have hlog : 0 ≤ Real.log (Ideal.absNorm P) := Real.log_nonneg hNpos
  refine mul_le_mul_of_nonneg_right ?_ hlog
  have h1 : ((ram P - 1 : ℕ) : ℝ) ≤ ((e - 1 : ℕ) : ℝ) := by
    have := hle P hP
    exact_mod_cast Nat.sub_le_sub_right this 1
  rwa [hcast] at h1

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★上からの `≲`（`log-cond` の言葉） -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**`log-diff` の差は `(e − 1)·log-cond` 以下**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★`§9-968`（`hlow` の大域形）を `§9-966`（導手の等式）で読み替え、
`e_P ≤ e` で係数を揃えただけである。★★定数は `0`。 -/
theorem hlow_logCond_pointwise {X : Scheme.{0}} (E : X.IdealSheafData)
    (P : Type) (fldF fldK : P → IntermediateField ℚ ℂ)
    (hnfF : ∀ x, NumberField (fldF x)) (hnfK : ∀ x, NumberField (fldK x))
    (alg : ∀ x, Algebra (fldF x) (fldK x))
    (deceq : ∀ x, haveI := hnfK x; DecidableEq (Ideal (𝓞 (fldK x))))
    (e : ℕ) (he : 0 < e)
    (xE : ∀ x, haveI := hnfK x; specRingOfIntegers (fldK x) ⟶ X)
    (S : ∀ x, haveI := hnfK x; Finset (Ideal (𝓞 (fldK x))))
    (hmax : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp.IsMaximal)
    (hne : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp ≠ ⊥)
    (hI : ∀ x, haveI := hnfK x; pullbackIdeal (fldK x) E (xE x) ≠ 0)
    (hmem : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp ∣ pullbackIdeal (fldK x) E (xE x))
    (hsupp : ∀ x, haveI := hnfK x; ∀ Q : Ideal (𝓞 (fldK x)), Irreducible Q →
      Q ∣ pullbackIdeal (fldK x) E (xE x) → Q ∈ S x)
    (pr : ∀ x, haveI := hnfK x; haveI := hnfF x;
      Ideal (𝓞 (fldK x)) → Ideal (𝓞 (fldF x)))
    (ram : ∀ x, haveI := hnfK x; Ideal (𝓞 (fldK x)) → ℕ)
    (hram1 : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, 1 ≤ ram x Pp)
    (hramle : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, ram x Pp ≤ e)
    (hfac : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x; ∀ Pp ∈ S x,
      Pp ^ (ram x Pp) = Ideal.map (algebraMap (𝓞 (fldF x)) (𝓞 (fldK x))) (pr x Pp))
    (htame : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x; ∀ Pp ∈ S x,
      (Algebra.intTrace (𝓞 (fldF x)) (𝓞 (fldK x))) 1 ∉ pr x Pp)
    (hD : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x;
      differentIdeal (𝓞 (fldF x)) (𝓞 (fldK x)) ≠ 0)
    (hprodne : ∀ x, haveI := hnfK x; (∏ Pp ∈ S x, Pp ^ (ram x Pp - 1)) ≠ 0)
    (hramsupp : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x;
      ∀ Q : Ideal (𝓞 (fldK x)), ∀ _ : Q.IsPrime,
        ¬ Algebra.IsUnramifiedAt (𝓞 (fldF x)) Q → Q ∈ S x) :
    BDle
      (fun x => haveI := hnfK x; ((e : ℝ) - 1) * logCond (fldK x) E (xE x))
      (fun x => haveI := hnfK x; haveI := hnfF x;
        logDiffOfField (fldK x) - logDiffOfField (fldF x)) := by
  refine bdle_of_le _ _ (fun x => ?_)
  haveI := hnfF x
  haveI := hnfK x
  haveI := deceq x
  letI := alg x
  have hfr : (0:ℝ) < (Module.finrank ℚ (fldK x) : ℝ) := by exact_mod_cast Module.finrank_pos
  have h1 := logDiff_tower_le_sum_of_totallyRamified_tame (fldF x) (fldK x) (S x) (hne x)
    (pr x) (ram x) (hram1 x) (hfac x) (htame x) (hD x) (hprodne x) (hramsupp x)
  refine le_trans h1 ?_
  rw [logCond_eq_sum_log (fldK x) E (xE x) (hI x) (S x) (hmax x) (hne x) (hmem x) (hsupp x),
    ← mul_div_assoc]
  exact (div_le_div_iff_of_pos_right hfr).mpr
    (sum_sub_one_le_mul (fldK x) e he (S x) (hne x) (ram x) (hramle x))

/-! ## ★出典の紐付け(`.src`) -/

def sum_sub_one_le_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(e_P ≤ e なら係数は e − 1 で揃う)",
    sectionId := "genell-prop-1-7" }

def hlow_logCond_pointwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(log-diff の差は (e − 1)·log-cond 以下)",
    sectionId := "genell-prop-1-7" }

def hlow_logCond_pointwise.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logDiff_tower_le_sum_of_totallyRamified_tame(hlow の大域形、§9-968)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiff_tower_le_sum_of_totallyRamified_tame") 4,
    .citation "[ABC3]" "logCond_eq_sum_log(導手の等式、§9-966)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logCond_eq_sum_log") 3,
    .implicitStep
      ("★★★★★★測定(2026-08-29): §9-967(hup)と本ファイル(hlow)で" ++
       "**log-diff の差が log-cond で両側から挟まれた**: " ++
       "(1 − 1/e)·log-cond ≲ log-diff(K) − log-diff(F) ≲ (e − 1)·log-cond。" ++
       "★どちらの ≲ も定数 0 である(slack も Σ も要らない)") 6,
    .implicitStep
      ("★★★原文 Proposition 1.7, (i) との差は**log-cond_D(底の側の導手)の項**である" ++
       "——分岐の側からは (1 − 1/e) と (e − 1) の 2 つの係数しか出ず、" ++
       "原文の形(log-cond_E − log-cond_D ≲ … ≲ (1 − 1/e)·log-cond_E)にするには" ++
       "底の側の導手が要る。★これが Proposition 1.7 に残る最後の差である") 5 ]

end ABC3.Found.GenEll
