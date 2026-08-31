/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjArcModel
import ABC3.Found.GenEll.HeightMetricDiff
import ABC3.Found.GenEll.GreenGlobal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★段 C2c —— **超平面因子の高さは素朴高さと BD-同値**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.6):
> The BD-class of htL depends only on the isomorphism class of the line bundle LQ on XQ.

## ★★★★★★★★★★★★★★★これは何か —— 段 C2c の到達点

★`§9-871`: `htArith F (hyperplaneArith N) xF = log H(x)/[F:ℚ]`（Fubini–Study のとき）
★★`§9-872`: 同じ因子で**差が連続**なら高さの差は一様に有界
★★★`§9-875`: `ℙᴺ_ℤ` の射影モデル `projArcModel`（コンパクト性の出どころ）

★★★★これを繋ぐと、**任意の**（Fubini–Study との差が連続な）計量について

    `|ht_E(x) − log H(x)/[F:ℚ]| ≤ C`   （`C` は `F` にも点にも依らない）

が出る。★これが `northcott_of_projModel` の `hcmp` が要求する形である。

## ★★★仮定の読み方

★`hcont`（`E.green − greenFS` が連続）は、原典 `Definition 1.1` の
「計量は連続である」の**正しい形**である——`§9-872` の測定のとおり、
Green 関数**単独**は台の上で発散するので連続ではありえず、
連続なのは同じ因子の 2 つの計量の**比**（＝差）である。

★★したがって本補題の仮定は「`E` の計量が Fubini–Study に対して連続」であり、
これは**与えられた計量についての条件**である（証明すべきものではない）。
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★`ℙᴺ` は複素点を持つ -/

/-- ★**`ℙᴺ_ℤ` は複素点を持つ**（`[1 : 1 : … : 1]`）。 -/
instance projComplexPointsNonempty (N : ℕ) :
    Nonempty (complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :=
  ⟨projPointOfCoords N (fun _ => (1 : ℂ)) 0 one_ne_zero⟩

/-! ## ★★★★★★★★★★★★★★超平面因子の BD-類 -/

/-- ★★★★★★★★★★★★★**Fubini–Study と任意の計量の高さは一様に近い**。

原文 (GenEll p.6):
> The BD-class of htL depends only on the isomorphism class of the line bundle LQ on XQ.

★`§9-872`（差の連続性版）に `§9-875`（`ℙᴺ` の射影モデル）を渡すだけである。 -/
theorem htArith_hyperplane_sub_abs_le (N : ℕ)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p)) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F]
        (xF : specRingOfIntegers F ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)),
        |htArith F E xF - htArith F (hyperplaneArith N) xF| ≤ C :=
  htArith_sub_abs_le_of_diff (projArcModel N) E (hyperplaneArith N) hdiv hcont

/-- ★★★★★★★★★★★★★★★**段 C2c** —— 超平面因子の高さは素朴高さと一様に近い。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `|ht_E(x) − log H(x)/[F:ℚ]| ≤ C`   （`C` は `F` にも点にも依らない）

★これが `northcott_of_projModel` の `hcmp` が要求する形である。
★★`§9-871`（Fubini–Study のときの等式）と `§9-872`＋`§9-875`（計量の取り替え）の合成。 -/
theorem abs_htArith_sub_log_mulHeight_le (N : ℕ)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p)) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (i₀ : Fin (N+1))
        (ψ : CommRingCat.of (HomogeneousLocalization.Away
            (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
          ⟶ CommRingCat.of (NumberField.RingOfIntegers F))
        (x : Fin (N+1) → NumberField.RingOfIntegers F), x ≠ 0 →
        (∀ k, x k = ψ.hom (projCoord N ℤ i₀ k) * x i₀) → x 0 ≠ 0 →
        |htArith F E (Spec.map ψ ≫ chartA N ℤ i₀)
          - Real.log (Height.mulHeight (fun k => ((x k : F))))
            / (Module.finrank ℚ F : ℝ)| ≤ C := by
  obtain ⟨C, hC, h⟩ := htArith_sub_abs_le_of_diff (projArcModel N) E (hyperplaneArith N)
    hdiv hcont
  refine ⟨C, hC, ?_⟩
  intro F _ _ i₀ ψ x hx hw h0
  rw [← htArith_hyperplaneArith F N i₀ ψ x hx hw h0]
  exact h F _

/-! ## ★出典の紐付け(`.src`) -/

def projComplexPointsNonempty.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(ℙᴺ_ℤ は複素点を持つ)",
    sectionId := "genell-prop-1-4" }

def htArith_hyperplane_sub_abs_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(Fubini–Study と任意の計量の高さは一様に近い)",
    sectionId := "genell-prop-1-4" }

def abs_htArith_sub_log_mulHeight_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(段 C2c——超平面因子の高さは素朴高さと一様に近い)",
    sectionId := "genell-prop-1-4" }

def abs_htArith_sub_log_mulHeight_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_hyperplaneArith(Fubini–Study のときの等式、§9-871)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_hyperplaneArith") 2,
    .citation "[ABC3]" "htArith_sub_abs_le_of_diff(差の連続性版、§9-872)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_sub_abs_le_of_diff") 2,
    .citation "[ABC3]" "projArcModel(ℙᴺ_ℤ の射影モデル、§9-875)"
      (.inProject "ABC3" "ABC3.Found.GenEll.projArcModel") 2,
    .implicitStep
      ("★仮定 hcont(E.green − greenFS が連続)は、原典 Definition 1.1 の" ++
       "「計量は連続である」の**正しい形**である——§9-872 の測定のとおり、" ++
       "Green 関数**単独**は台の上で発散するので連続ではありえず、" ++
       "連続なのは同じ因子の 2 つの計量の**比**(＝差)である") 3,
    .implicitStep
      ("★★したがって本補題の仮定は「E の計量が Fubini–Study に対して連続」であり、" ++
       "これは**与えられた計量についての条件**である(証明すべきものではない)") 3,
    .implicitStep
      ("★★★これが northcott_of_projModel の hcmp が要求する形である。" ++
       "残るのは hcmp の指数の形(mulHeight ≤ exp(ht + const))へ整える段だけで、" ++
       "mulHeight_eq_exp_htArith(§9-867)がその橋である") 4 ]

end ABC3.Found.GenEll
