/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HyperplaneBDClass
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★正規化の橋 —— `htArith`（絶対）と `mulHeight`（相対）（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★これは何か —— 測定と、その橋

★★★★★**測定（2026-08-28）**: `northcott_of_projModel` の仮定 `hcmp` は

    `Height.mulHeight (crd p) ≤ exp(ht p + const)`

だが、`Height.mulHeight` は**相対高さ** `H_K`（`[K:ℚ]` 乗根を取らない）であり、
`htArith` は `[F:ℚ]` で割った**絶対高さ**である。
★したがって `log H_K = [F:ℚ]·ht` であって、
**`ht` をそのまま入れると次数の分だけ足りない**。

★★橋を架けるには**次数の上界 `d`** を使う:

    `log H_K ≤ d·(ht + C)`   （`ht + C ≥ 0` は `log H_K ≥ 0` から自動）

★★★したがって `hcmp` は `ht' := d·ht` の形で満たされる。
★★★★`{p | ht p ≤ C'} ⊆ {p | ht' p ≤ d·C'}` なので、Northcott の結論は変わらない。

## ★★★機構 —— 実数の不等式 2 本だけ

★`|ht − log H/d_K| ≤ C` から `log H/d_K ≤ ht + C`、
`log H ≥ 0` から `ht + C ≥ 0`、あとは `d_K ≤ d` を掛けるだけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★実数の橋 -/

/-- ★**次数の上界で相対高さの対数を抑える**。

★`ht + C ≥ 0` は `H ≥ 0` から**自動で出る**（仮定に置く必要がない）。 -/
theorem log_le_of_abs_sub_le (H t C dk d : ℝ) (hdk : 0 < dk) (hdkd : dk ≤ d)
    (hH : 0 ≤ H) (h : |t - H / dk| ≤ C) : H ≤ d * t + d * C := by
  have hHdk : 0 ≤ H / dk := div_nonneg hH hdk.le
  have h1 : H / dk ≤ t + C := by have := (abs_le.1 h).1; linarith
  have hlow : 0 ≤ t + C := le_trans hHdk h1
  have h2 : H ≤ dk * (t + C) := by rw [div_le_iff₀ hdk] at h1; linarith
  have h3 : dk * (t + C) ≤ d * (t + C) := mul_le_mul_of_nonneg_right hdkd hlow
  linarith

/-- ★★**相対高さは `exp(d·ht + d·C)` で抑えられる** —— `hcmp` の形。 -/
theorem mulHeight_le_exp_of_abs_sub_le (F : Type) [Field F] [NumberField F]
    {ι : Type} [Finite ι] (y : ι → F) (t C : ℝ) (d : ℕ)
    (hd : (Module.finrank ℚ F : ℝ) ≤ (d : ℝ))
    (h : |t - Real.log (Height.mulHeight y) / (Module.finrank ℚ F : ℝ)| ≤ C) :
    Height.mulHeight y ≤ Real.exp ((d : ℝ) * t + (d : ℝ) * C) := by
  have hpos : (0:ℝ) < Height.mulHeight y :=
    lt_of_lt_of_le zero_lt_one (Height.one_le_mulHeight _)
  have hlog : 0 ≤ Real.log (Height.mulHeight y) :=
    Real.log_nonneg (Height.one_le_mulHeight _)
  have hdk : (0:ℝ) < (Module.finrank ℚ F : ℝ) := by exact_mod_cast Module.finrank_pos
  rw [← Real.log_le_iff_le_exp hpos]
  exact log_le_of_abs_sub_le _ t C _ _ hdk hd hlog h

/-! ## ★★★★★★★★★★★`hcmp` の形 -/

/-- ★★★★★★★★★★★**`northcott_of_projModel` の `hcmp` が満たされる**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `Height.mulHeight (crd) ≤ exp( d·ht_E(x) + const )`

★`§9-876`（段 C2c）に正規化の橋を掛けただけである。
★★`ht' := d·ht` と取れば `northcott_of_projModel` の `hcmp` そのものになる。 -/
theorem mulHeight_le_exp_htArith (N : ℕ) (d : ℕ)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p)) :
    ∃ const : ℝ, 0 ≤ const ∧
      ∀ (F : Type) [Field F] [NumberField F], (Module.finrank ℚ F : ℝ) ≤ (d : ℝ) →
        ∀ (i₀ : Fin (N+1))
          (ψ : CommRingCat.of (HomogeneousLocalization.Away
              (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
            ⟶ CommRingCat.of (NumberField.RingOfIntegers F))
          (x : Fin (N+1) → NumberField.RingOfIntegers F), x ≠ 0 →
          (∀ k, x k = ψ.hom (projCoord N ℤ i₀ k) * x i₀) → x 0 ≠ 0 →
          Height.mulHeight (fun k => ((x k : F)))
            ≤ Real.exp ((d : ℝ) * htArith F E (Spec.map ψ ≫ chartA N ℤ i₀) + const) := by
  obtain ⟨C, hC, h⟩ := abs_htArith_sub_log_mulHeight_le N E hdiv hcont
  refine ⟨(d : ℝ) * C, mul_nonneg (Nat.cast_nonneg d) hC, ?_⟩
  intro F _ _ hd i₀ ψ x hx hw h0
  exact mulHeight_le_exp_of_abs_sub_le F _ _ _ d hd (h F i₀ ψ x hx hw h0)

/-! ## ★出典の紐付け(`.src`) -/

def log_le_of_abs_sub_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(次数の上界で相対高さの対数を抑える)",
    sectionId := "genell-prop-1-4" }

def mulHeight_le_exp_of_abs_sub_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(相対高さは exp(d·ht + d·C) で抑えられる)",
    sectionId := "genell-prop-1-4" }

def mulHeight_le_exp_htArith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(northcott_of_projModel の hcmp が満たされる)",
    sectionId := "genell-prop-1-4" }

def mulHeight_le_exp_htArith.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "abs_htArith_sub_log_mulHeight_le(段 C2c、§9-876)"
      (.inProject "ABC3" "ABC3.Found.GenEll.abs_htArith_sub_log_mulHeight_le") 2,
    .citation "[mathlib]" "Height.one_le_mulHeight(相対高さは 1 以上)"
      (.inMathlib "Height.one_le_mulHeight") 2,
    .implicitStep
      ("★★★★★測定(2026-08-28): northcott_of_projModel の仮定 hcmp は " ++
       "Height.mulHeight (crd p) ≤ exp(ht p + const) だが、" ++
       "Height.mulHeight は**相対高さ** H_K([K:ℚ] 乗根を取らない)であり、" ++
       "htArith は [F:ℚ] で割った**絶対高さ**である。" ++
       "★log H_K = [F:ℚ]·ht なので **ht をそのまま入れると次数の分だけ足りない**") 4,
    .implicitStep
      ("★★橋を架けるには**次数の上界 d** を使う: log H_K ≤ d·(ht + C)。" ++
       "ht + C ≥ 0 は log H_K ≥ 0 から**自動で出る**(仮定に置く必要がない)。" ++
       "★したがって hcmp は ht' := d·ht の形で満たされ、" ++
       "{p | ht p ≤ C'} ⊆ {p | ht' p ≤ d·C'} なので Northcott の結論は変わらない") 3 ]

end ABC3.Found.GenEll
