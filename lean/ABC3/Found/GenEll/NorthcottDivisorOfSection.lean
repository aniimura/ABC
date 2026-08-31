/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DivisorOfSectionComap
import ABC3.Found.GenEll.NorthcottGlobalToProj
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★切断の零点因子についての Northcott（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★これは何か —— `Proposition 1.4, (iv)` の到達形

★`§9-926` で `div(s₀) = ψ^*超平面` が出たので、
**与えられた算術因子 `D` が「ある切断の零点因子＋引き戻した計量」であれば
その高さについて Northcott が成り立つ**と言える:

    `D.divisor = div(s₀)` かつ `D.green = Ē.green ∘ ψ`  ⟹  `{p | ht_D(x_p) ≤ C}` は有限

★★これが原文の
「`L_ℚ` が豊富 ⟹ `ht_L` で有界な `X(ℚ̄)^{≤d}` の点は有限」の**構成の側での形**である。

## ★★★これまでの鎖（本セッション 25 ブロック）

| 段 | 内容 |
|---|---|
| `§9-907`〜`§9-911` | **段 E1d** ——`ψ : X ⟶ ℙᴺ` が貼れた |
| `§9-912`〜`§9-913` | `ψ⁻¹(D₊(x_i)) = X_{s_i}`・**`ψ` は埋め込み** |
| `§9-914`〜`§9-916` | 被覆・アフィン性・部分族 |
| `§9-917`〜`§9-920` | **段 E3** ——チャート写像は全射／**ample なら埋め込める** |
| `§9-921`〜`§9-922` | 高さの側 ——`ht_{ψ^*Ē}(x) = log H(座標)/[F:ℚ]` |
| `§9-923`〜`§9-926` | **段 E0** ——`div(s₀) = ψ^*超平面` |
| `§9-927` | **本ファイル ——到達** |

## ★残っている仮定（明示）

★★★★**幾何の側は「有限アフィン自明化被覆」に集約されている**
（`U`・`hU`・`eU`——`X` が分離的・準コンパクトで `M` が可逆なら出るはず）。

★★点の側に残るのは 2 つだけ:

* `chart`・`xF`——**点がどれかのチャートを通ること**（`X` が固有なら自動のはず）
* `hinj`——**座標が点を分けること**（`ψ` が埋め込み（`§9-920`）であることの点版）

★原文もこの 2 つを（`X(ℚ̄)^{≤d}` という書き方の中に）暗黙に含んでいる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★算術因子の外延性 -/

/-- ★**算術因子は因子と Green 関数で決まる**。 -/
theorem ArithCartier.ext' {D D' : ArithCartier X}
    (h1 : D.divisor = D'.divisor) (h2 : D.green = D'.green) : D = D' := by
  cases D; cases D'; simp_all

/-! ## ★★★★★★★★★★★★★★★★★★★★到達 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★**切断の零点因子についての Northcott**
—— `Proposition 1.4, (iv)` の到達形。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `D.divisor = div(s₀)` かつ `D.green = Ē.green ∘ ψ`
      ⟹ `{p | ht_D(x_p) ≤ C}` は有限

★★`§9-926`（`div(s₀) = ψ^*超平面`）で `D = Ē.comap ψ` になり、
`§9-921`（貼った射に沿った Northcott）がそのまま当たる。 -/
theorem northcott_of_divisorOfSection (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ} (d : ℕ)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    {ι : Type} (U : ι → X.Opens) (hcovU : (⨆ j, U j) = ⊤)
    (hU : ∀ j, IsAffineOpen (U j))
    (eU : ∀ j, (restrictPresheafFunctor X (U j)).obj M ≅ 𝟙_ (PresheafModulesOn X (U j)))
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p))
    (D : ArithCartier X)
    (hDdiv : D.divisor = divisorOfSection M (s 0))
    (hDgreen : D.green = fun p => E.green (p ≫ globalToProj M hM φ s hcov))
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
    (C : ℝ) :
    {p : P | haveI := hnf p;
      htArith (fld p) D (xF p ≫ (nonVanishing M (s (chart p))).ι) ≤ C}.Finite := by
  have hD : D = E.comap (globalToProj M hM φ s hcov) := by
    refine ArithCartier.ext' ?_ hDgreen
    rw [hDdiv, ArithCartier.comap_divisor, hdiv]
    exact divisorOfSection_eq_comap_globalToProj M hM φ s hcov haff U hcovU hU eU
  rw [hD]
  exact northcott_globalToProj M hM d φ s hcov E hdiv hcont fld hnf hdeg chart xF x hx hw h0
    idx hinj C

/-! ## ★出典の紐付け(`.src`) -/

def ArithCartier.ext'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(算術因子は因子と Green 関数で決まる)",
    sectionId := "genell-def-1-2" }

def northcott_of_divisorOfSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(切断の零点因子についての Northcott——到達形)",
    sectionId := "genell-prop-1-4" }

def northcott_of_divisorOfSection.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "divisorOfSection_eq_comap_globalToProj(段 E0、§9-926)"
      (.inProject "ABC3" "ABC3.Found.GenEll.divisorOfSection_eq_comap_globalToProj") 5,
    .citation "[ABC3]" "northcott_globalToProj(貼った射に沿った Northcott、§9-921)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_globalToProj") 4,
    .citation "[ABC3]" "globalToProj(貼った射、§9-911)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalToProj") 3,
    .implicitStep
      ("★★★★幾何の側は「有限アフィン自明化被覆」(U・hU・eU)に集約されている" ++
       "——X が分離的・準コンパクトで M が可逆なら出るはずのものである") 4,
    .implicitStep
      ("★★点の側に残るのは 2 つだけ: " ++
       "(1) chart・xF——点がどれかのチャートを通ること(X が固有なら自動のはず)、" ++
       "(2) hinj——座標が点を分けること(ψ が埋め込み(§9-920)であることの点版)。" ++
       "★原文もこの 2 つを X(ℚ̄)^{≤d} という書き方の中に暗黙に含んでいる") 5 ]

end ABC3.Found.GenEll
