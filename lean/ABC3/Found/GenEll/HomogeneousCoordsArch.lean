/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HomogeneousCoords
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★無限素点の整合は構成から出る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★これは何か —— `§9-940` の 2 本目が消える

`§9-940` は点ごとに**無限素点の整合**

    `∀ p v, ∃ i ρ, (複素点が D₊(x_i) を通る) ∧ (σ_v(x_k) = ρ(x_k/x_i)·σ_v(x_i)) ∧ (x_i ≠ 0)`

を仮定として受けていた。★★★★本ファイルはそれを **`§9-941` の同次座標から導く**。

## ★★★機構 —— 無限素点は生成点を経由する

`archRingHom F v = v.embedding ∘ (𝓞_F ↪ F)` だから、`Spec` の反変性で

    `archPoint(x) ≫ ψ = Spec(v.embedding) ≫ (生成点 ≫ ψ)`

★したがって**生成点のチャートがそのまま無限素点のチャートになる**
——`ρ_v ≔ (生成点のチャート射) ≫ v.embedding` と取ればよい。
★★`x_i ≠ 0` は `§9-941` が与えている（`x_i = b`、分母払いの分母）。

## ★これで残ったのは有限素点の整合だけ

★★★`§9-940` の 2 本のうち **1 本が消えた**。残るのは

    `∀ p Q, ∃ i y, (局所化した点が X_{s_i} を通る) ∧ (𝔞_Q = (x_i)) ∧ ((s_0/s_i)(y)·x_i = x_0)`

であり、★これは `v_Q(x_j)` が最小の `j` を取る段
——`𝓞_F` が Dedekind であることを使う——である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial NumberField ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★無限素点は生成点を経由する -/

/-- ★★★**`Spec ℂ ⟶ Spec 𝓞_F` は生成点 `Spec F` を経由する**。

★`archRingHom F v = v.embedding ∘ (𝓞_F ↪ F)` と `Spec` の反変性だけである。 -/
theorem archSpecMap_factor (F : Type) [Field F] [NumberField F] (v : InfinitePlace F) :
    archSpecMap F v = Spec.map (CommRingCat.ofHom (v.embedding : F →+* ℂ))
      ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) := by
  rw [archSpecMap, archRingHom, ← Spec.map_comp]
  rfl

/-! ## ★★★★★★★★★★★★★★★★★無限素点のチャートは生成点のチャートである -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★**生成点のチャートがそのまま無限素点のチャートになる**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`ρ_v ≔ (生成点のチャート射) ≫ v.embedding` と取ればよい。
★★これが `§9-940` の**無限素点の整合**そのものである。 -/
theorem archChart_of_coords (N : ℕ) (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X)
    (ψ : X ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (i : Fin (N + 1))
    (hx : Set.range (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ ψ).base
      ⊆ Set.range (chartA N ℤ i).base)
    (x : Fin (N + 1) → 𝓞 F)
    (hrel : ∀ k, ((x k : F)) = projPointCoord N ℤ F
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ ψ) i hx k * ((x i : F)))
    (v : InfinitePlace F) :
    ∃ (ρ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ),
      archPoint xF v ≫ ψ = Spec.map ρ ≫ chartA N ℤ i ∧
      ∀ k, archRingHom F v (x k)
        = ρ.hom (projCoord N ℤ i k) * archRingHom F v (x i) := by
  set q := Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ ψ with hq
  refine ⟨projChartHom N ℤ F q i hx ≫ CommRingCat.ofHom (v.embedding : F →+* ℂ), ?_, ?_⟩
  · rw [Spec.map_comp, Category.assoc, specMap_projChartHom N ℤ F q i hx, hq,
      archPoint, archSpecMap_factor]
    simp only [Category.assoc]
  · intro k
    have h1 : (projChartHom N ℤ F q i hx ≫ CommRingCat.ofHom (v.embedding : F →+* ℂ)).hom
        (projCoord N ℤ i k) = v.embedding (projPointCoord N ℤ F q i hx k) := by
      rw [CommRingCat.hom_comp, CommRingCat.hom_ofHom]
      rfl
    rw [h1]
    show v.embedding (algebraMap (𝓞 F) F (x k)) = _
    show _ = _ * v.embedding (algebraMap (𝓞 F) F (x i))
    rw [show algebraMap (𝓞 F) F (x k) = ((x k : F)) from rfl,
      show algebraMap (𝓞 F) F (x i) = ((x i : F)) from rfl, hrel k, map_mul]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★**`§9-940` の「無限素点の整合」はもう仮定ではない**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-941` の同次座標を与えれば、`§9-940` の 2 本目の存在命題は**定理として出る**。 -/
theorem harch_of_homogeneous_coords (N : ℕ) (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X)
    (ψ : X ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (i : Fin (N + 1))
    (hx : Set.range (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ ψ).base
      ⊆ Set.range (chartA N ℤ i).base)
    (x : Fin (N + 1) → 𝓞 F) (hxi : x i ≠ 0)
    (hrel : ∀ k, ((x k : F)) = projPointCoord N ℤ F
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ ψ) i hx k * ((x i : F)))
    (v : InfinitePlace F) :
    ∃ (j : Fin (N + 1)) (ρ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X j))
      ⟶ CommRingCat.of ℂ),
      archPoint xF v ≫ ψ = Spec.map ρ ≫ chartA N ℤ j ∧
      (∀ k, archRingHom F v (x k)
        = ρ.hom (projCoord N ℤ j k) * archRingHom F v (x j)) ∧
      x j ≠ 0 := by
  obtain ⟨ρ, hfac, hcv⟩ := archChart_of_coords N F xF ψ i hx x hrel v
  exact ⟨i, ρ, hfac, hcv, hxi⟩

/-! ## ★出典の紐付け(`.src`) -/

def archSpecMap_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(Spec ℂ ⟶ Spec 𝓞_F は生成点を経由する)",
    sectionId := "genell-prop-1-4" }

def archChart_of_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(生成点のチャートがそのまま無限素点のチャートになる)",
    sectionId := "genell-prop-1-4" }

def harch_of_homogeneous_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(無限素点の整合はもう仮定ではない)",
    sectionId := "genell-prop-1-4" }

def harch_of_homogeneous_coords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_homogeneous_coords(同次座標の構成、§9-941)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_homogeneous_coords") 3,
    .citation "[ABC3]" "specMap_projChartHom(チャート射で点が復元される)"
      (.inProject "ABC3" "ABC3.Found.GenEll.specMap_projChartHom") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-940 の 2 本の存在命題のうち" ++
       "**無限素点の整合は定理である**。" ++
       "archRingHom F v = v.embedding ∘ (𝓞_F ↪ F) だから Spec の反変性で" ++
       "archPoint(x) ≫ ψ = Spec(v.embedding) ≫ (生成点 ≫ ψ) となり、" ++
       "生成点のチャートがそのまま無限素点のチャートになる") 5,
    .implicitStep
      ("★残るのは有限素点の整合だけである——" ++
       "∀ p Q, ∃ i y, (局所化した点が X_{s_i} を通る) ∧ (𝔞_Q = (x_i)) " ++
       "∧ ((s_0/s_i)(y)·x_i = x_0)。" ++
       "★これは v_Q(x_j) が最小の j を取る段(𝓞_F が Dedekind であることを使う)である") 4 ]

end ABC3.Found.GenEll
