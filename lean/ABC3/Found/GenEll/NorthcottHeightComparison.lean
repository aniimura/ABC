/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottGlobalToProj
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★与えられた因子への移送 —— 高さの比較だけでよい（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★これは何か —— 仮定を**点ごとの等式**まで弱める

`§9-882` の `northcott_of_veryAmple` は仮定に

    `hDn : D.npow n = E.comap ψ`（**算術因子の等式**）

を持っていた。★これは段 E0（切断の零点因子 `divisorOfSection`）を要求する重い仮定である。

★★★★**しかし Northcott が実際に使うのは高さの等式だけ**である:

    `ht_{ψ^*Ē}(x) = n · ht_D(x)`（**点ごと**）

★本ファイルはその形に弱めた版を取る。`hDn` からは
`htArith_npow`（`§9-880`）で直ちに出るので、**一般化であって弱体化ではない**。

## ★これで何が残ったか

★★★`Proposition 1.4, (iv)` に残るのは

1. 点ごとの高さの等式（＝段 E0、`D^{⊗n}` と `ψ^*Ē` の一致）
2. 点がどれかのチャートを通ること（`X` が固有なら自動のはず）
3. 座標が点を分けること（`hinj`）——`ψ` が埋め込み（`§9-920`）であることの点版

の 3 つだけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★高さが定数倍なら Northcott は移る -/

/-- ★★★**高さが定数倍なら Northcott は移る**（純粋に順序の話）。 -/
theorem northcott_le_of_mul {P : Type} (h1 h2 : P → ℝ) (n : ℕ)
    (hcmp : ∀ p, h1 p = (n : ℝ) * h2 p)
    (hfin : ∀ C : ℝ, {p : P | h1 p ≤ C}.Finite) (C : ℝ) :
    {p : P | h2 p ≤ C}.Finite := by
  refine Set.Finite.subset (hfin ((n : ℝ) * C)) (fun p hp => ?_)
  have hp' : h2 p ≤ C := hp
  show h1 p ≤ (n : ℝ) * C
  rw [hcmp p]
  exact mul_le_mul_of_nonneg_left hp' (Nat.cast_nonneg n)

/-- ★★★★**因子の等式から点ごとの高さの等式が出る**——`§9-882` の仮定はこれより強い。 -/
theorem htArith_comap_eq_nsmul {Y : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (D : ArithCartier X) (E : ArithCartier Y) (ψ : X ⟶ Y) (n : ℕ)
    (hDn : D.npow n = E.comap ψ)
    (xF : specRingOfIntegers F ⟶ X)
    (hDne : pullbackIdeal F D.divisor xF ≠ 0) :
    htArith F (E.comap ψ) xF = (n : ℝ) * htArith F D xF := by
  rw [← hDn, htArith_npow F D xF hDne n]

/-! ## ★★★★★★★★★★★★★★★★★与えられた因子についての Northcott -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★**与えられた算術因子 `D` についての Northcott**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★仮定 `hht` は**点ごとの高さの等式**であり、
`§9-882` の因子の等式 `D.npow n = E.comap ψ` より弱い
（`htArith_comap_eq_nsmul` で後者から出る）。

★★これで `Proposition 1.4, (iv)` に残るのは
(1) `hht`（＝段 E0）、(2) 点がチャートを通ること、(3) `hinj` の 3 つだけである。 -/
theorem northcott_globalToProj_of_height (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ} (d : ℕ)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p))
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (chart : P → Fin (N + 1))
    (xF : ∀ p, haveI := hnf p;
      specRingOfIntegers (fld p) ⟶ (nonVanishing M (s (chart p))).toScheme)
    (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
    (hx : ∀ p, x p ≠ 0)
    (hw : ∀ p, haveI := hnf p; ∀ k, x p k =
      (Spec.preimage (xF p ≫ globalChartMorphism M hM φ s (chart p))).hom
        (projCoord N ℤ (chart p) k) * x p (chart p))
    (h0 : ∀ p, x p 0 ≠ 0)
    (idx : Fin (N + 1))
    (hinj : Function.Injective (fun (p : P) (k : Fin (N + 1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (D : ArithCartier X) (n : ℕ)
    (hht : ∀ p, haveI := hnf p;
      htArith (fld p) (E.comap (globalToProj M hM φ s hcov))
          (xF p ≫ (nonVanishing M (s (chart p))).ι)
        = (n : ℝ) * htArith (fld p) D (xF p ≫ (nonVanishing M (s (chart p))).ι))
    (C : ℝ) :
    {p : P | haveI := hnf p;
      htArith (fld p) D (xF p ≫ (nonVanishing M (s (chart p))).ι) ≤ C}.Finite := by
  refine northcott_le_of_mul
    (fun p => haveI := hnf p; htArith (fld p) (E.comap (globalToProj M hM φ s hcov))
      (xF p ≫ (nonVanishing M (s (chart p))).ι))
    (fun p => haveI := hnf p; htArith (fld p) D (xF p ≫ (nonVanishing M (s (chart p))).ι))
    n hht (fun C' => ?_) C
  exact northcott_globalToProj M hM d φ s hcov E hdiv hcont fld hnf hdeg chart xF x hx hw h0
    idx hinj C'

/-! ## ★出典の紐付け(`.src`) -/

def northcott_le_of_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(高さが定数倍なら Northcott は移る)",
    sectionId := "genell-prop-1-4" }

def htArith_comap_eq_nsmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(因子の等式から点ごとの高さの等式が出る)",
    sectionId := "genell-prop-1-4" }

def northcott_globalToProj_of_height.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(与えられた算術因子についての Northcott)",
    sectionId := "genell-prop-1-4" }

def northcott_globalToProj_of_height.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_globalToProj(貼った射に沿った Northcott、§9-921)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_globalToProj") 3,
    .citation "[ABC3]" "htArith_npow(ht_{D^n} = n·ht_D、§9-880)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_npow") 1,
    .implicitStep
      ("★★仮定 hht は**点ごとの高さの等式**であり、" ++
       "§9-882 の因子の等式 D.npow n = E.comap ψ より弱い" ++
       "(htArith_comap_eq_nsmul で後者から出る)。" ++
       "★段 E0(切断の零点因子)を通さずに済むので、Northcott の機構と E0 が**分離した**") 3,
    .implicitStep
      ("★★★これで Proposition 1.4, (iv) に残るのは " ++
       "(1) hht(＝段 E0、D^{⊗n} と ψ^*Ē の一致)、" ++
       "(2) 点がどれかのチャートを通ること(X が固有なら自動のはず)、" ++
       "(3) 座標が点を分けること(hinj——ψ が埋め込みであることの点版)" ++
       "の 3 つだけである") 5 ]

end ABC3.Found.GenEll
