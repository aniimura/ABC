/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HyperplaneChartSection
import ABC3.Found.GenEll.PullbackLocalization
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★段 C2c-1 —— 点に沿った超平面の引き戻しは `(x₀/x_i)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★これは何か —— 段 C2c-1 の到達点

段 C2c-1 の文言は「点 `xF : Spec 𝓞_F ⟶ ℙᴺ_ℤ` に沿った超平面因子の引き戻しが
**座標の生成するイデアル**になること」であった。★本ファイルがそれを取る:

    `pullbackIdeal F (hyperplaneIdeal N ℤ) xF = (ψ(x₀/x_i))`

（`xF` がチャート `D₊(x_i)` を通る、すなわち `xF = Spec ψ ≫ chartA i` のとき）。

## ★★★積み上がった段

| 段 | 内容 | 場所 |
|---|---|---|
| (a) | 可換な四角 `Φ ∘ Ψ = Ψ′ ∘ σ` | `§9-862` |
| (b) | `ker σ = (x₀)` | `§9-861` |
| (c) | `ker Φ = (x₀/x_i)`（代数の側） | `§9-863` |
| C2c-1d | 切断環 `Γ(Proj 𝒜, D₊(x_i))` への翻訳 | `§9-864` |
| **C2c-1e** | **点に沿った引き戻しへの接続** | ★★**本ファイル** |

## ★★★★★機構 —— 同型の逆と同型が打ち消し合う

★`pullbackIdealOf B D x` は定義上 `(D.comap x).ideal ⊤` を `ΓSpecIso` で読んだものである。
★★開埋め込みに沿った `comap` は mathlib の `ideal_comap_of_isOpenImmersion` で
**像での値の引き戻し**になるので、`§9-864` がそのまま入る。
★★★残るのは

    `(ΓSpecIso⁻¹ ≫ appIso⁻¹) ≫ (chartA.app ≫ eqToHom ≫ ΓSpecIso) = 𝟙`

という**同型の打ち消し**だけで、`Scheme.Hom.appIso_hom` を逆向きに使えば出る。

## ★測定の記録

★`rw [← Scheme.Hom.appIso_hom]` は**そのままでは通らない**——`≫` が右結合なので
`app ≫ (map ≫ ΓSpecIso)` となっており、`app ≫ map` が**部分項ではない**。
先に `rw [← Category.assoc …]` で括り直す必要がある（2026-08-28 実測）。
★★このとき出る "not type-correct under the `instances` transparency level" は
**副次的な症状**であって、原因は結合の向きである。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★同型の打ち消し -/

set_option maxHeartbeats 2000000 in
/-- ★★**`Γ(Proj 𝒜, D₊(x_i)) ≅ A⁰_{x_i}` の逆と表が打ち消し合う**。

★`rw [← Scheme.Hom.appIso_hom]` を通すには先に `← Category.assoc` で
`app ≫ (map ≫ ΓSpecIso)` を `(app ≫ map) ≫ ΓSpecIso` に括り直す必要がある。 -/
theorem chartIso_inv_comp (N : ℕ) (R : Type) [CommRing R] (i : Fin (N+1)) :
    ((Scheme.ΓSpecIso (CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i)))).inv
        ≫ (Scheme.Hom.appIso (chartA N R i) ⊤).inv) ≫
      (Scheme.Hom.app (chartA N R i) ((chartA N R i) ''ᵁ ⊤) ≫
        (Spec (CommRingCat.of (HomogeneousLocalization.Away
          (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
          (MvPolynomial.X i)))).presheaf.map
            (eqToHom (Scheme.Hom.preimage_image_eq (chartA N R i) ⊤).symm).op ≫
          (Scheme.ΓSpecIso (CommRingCat.of (HomogeneousLocalization.Away
            (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i)))).hom)
      = 𝟙 _ := by
  rw [← Category.assoc (Scheme.Hom.app (chartA N R i) ((chartA N R i) ''ᵁ ⊤)),
    ← Scheme.Hom.appIso_hom (chartA N R i) ⊤]
  simp

/-! ## ★★★★★★★★★★チャートに沿った引き戻し -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**チャートに沿った超平面の引き戻しは `(x₀/x_i)` である**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-864`（切断環で読んだ形）に mathlib の `ideal_comap_of_isOpenImmersion` を当て、
同型の打ち消し（`chartIso_inv_comp`）で `comap` が消える。 -/
theorem pullbackIdealOf_hyperplane_chart (N : ℕ) (R : Type) [CommRing R] (i : Fin (N+1)) :
    pullbackIdealOf (CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i)))
      (hyperplaneIdeal N R) (chartA N R i)
      = Ideal.span {projCoord N R i 0} := by
  rw [pullbackIdealOf, Scheme.IdealSheafData.equivOfIsAffine_apply,
    Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion (hyperplaneIdeal N R)
      (chartA N R i) ⟨⊤, isAffineOpen_top _⟩,
    hyperplaneIdeal_apply N R ⟨(chartA N R i) ''ᵁ ⊤,
      (isAffineOpen_top _).image_of_isOpenImmersion (chartA N R i)⟩,
    ker_app_hyperplaneι N R i, Ideal.comap_comap, Ideal.comap_comap,
    ← CommRingCat.hom_comp, ← CommRingCat.hom_comp, chartIso_inv_comp,
    CommRingCat.hom_id, Ideal.comap_id]

/-! ## ★★★★★★★★★★★点に沿った引き戻し —— 段 C2c-1 -/

/-- ★★★★★★★★★★★**点に沿った超平面の引き戻し** —— 段 C2c-1。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `pullbackIdealOf B (hyperplaneIdeal) (Spec ψ ≫ chartA i) = (ψ(x₀/x_i))`

★`pullbackIdealOf_specMap`（任意の環準同型に沿って拡大イデアルになる）に
`pullbackIdealOf_hyperplane_chart` を渡すだけである。 -/
theorem pullbackIdealOf_hyperplane_point (N : ℕ) (R : Type) [CommRing R] (i : Fin (N+1))
    (B : CommRingCat.{0})
    (ψ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i)) ⟶ B) :
    pullbackIdealOf B (hyperplaneIdeal N R) (Spec.map ψ ≫ chartA N R i)
      = Ideal.span {ψ.hom (projCoord N R i 0)} := by
  rw [pullbackIdealOf_specMap, pullbackIdealOf_hyperplane_chart, Ideal.map_span,
    Set.image_singleton]

/-- ★★**数体の点の場合** —— `pullbackIdeal F` の形。

★これが `northcott_of_projModel` の `hcmp` を作るときに読む形である
——`degFin` は `log N(pullbackIdeal)` だから。 -/
theorem pullbackIdeal_hyperplane_point (F : Type) [Field F] [NumberField F]
    (N : ℕ) (i : Fin (N+1))
    (ψ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of (NumberField.RingOfIntegers F))
    (xF : specRingOfIntegers F ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (hx : xF = Spec.map ψ ≫ chartA N ℤ i) :
    pullbackIdeal F (hyperplaneIdeal N ℤ) xF
      = Ideal.span {ψ.hom (projCoord N ℤ i 0)} := by
  rw [← pullbackIdealOf_eq_pullbackIdeal, hx, pullbackIdealOf_hyperplane_point]

/-! ## ★出典の紐付け(`.src`) -/

def chartIso_inv_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートの切断環の同型が打ち消し合うこと)",
    sectionId := "genell-prop-1-4" }

def pullbackIdealOf_hyperplane_chart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートに沿った超平面の引き戻しは (x₀/x_i))",
    sectionId := "genell-prop-1-4" }

def pullbackIdealOf_hyperplane_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-1——点に沿った超平面の引き戻し)",
    sectionId := "genell-prop-1-4" }

def pullbackIdeal_hyperplane_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-1——数体の点の場合)",
    sectionId := "genell-prop-1-4" }

def pullbackIdealOf_hyperplane_point.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ker_app_hyperplaneι(切断環で読んだ形、段 C2c-1d、§9-864)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_app_hyperplaneι") 2,
    .citation "[ABC3]" "pullbackIdealOf_specMap(環準同型に沿って拡大イデアルになる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdealOf_specMap") 2,
    .citation "[mathlib]" "IdealSheafData.ideal_comap_of_isOpenImmersion"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion") 2,
    .implicitStep
      ("★測定: rw [← Scheme.Hom.appIso_hom] は**そのままでは通らない**" ++
       "——≫ が右結合なので app ≫ (map ≫ ΓSpecIso) となっており、" ++
       "app ≫ map が**部分項ではない**。先に rw [← Category.assoc …] で括り直す。" ++
       "★★このとき出る not type-correct under the instances transparency level は" ++
       "**副次的な症状**であって、原因は結合の向きである(2026-08-28 実測)") 3,
    .implicitStep
      ("★★仮定 hx(点がチャート D₊(x_i) を通ること)は §9-C2b(exists_X_notMem——" ++
       "体値の点はどれかのチャートに入る)から作る。" ++
       "★残るのは degFin(= log N(pullbackIdeal))を素朴高さに繋ぐ段 C2c 後半である") 4 ]

end ABC3.Found.GenEll
