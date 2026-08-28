/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HyperplanePullbackChart
import ABC3.Found.GenEll.GlobalChartToProj
import ABC3.Found.GenEll.HyperplanePullback
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★点に沿った超平面の引き戻しは比の切断である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★これは何か —— 幾何側と数論側の接続

★`§9-885`（段 E3e）は**環の水準**で
`(globalAwayHom i)_* ( ker (Away.map hyperplaneHom (x_i)) ) = ( s_0/s_i )` を取った。
★★`§9-865`（段 C2c-1）は**点に沿った引き戻し**を計算する道具である。

★★★本ファイルはその 2 つを繋ぐ:

    `pullbackIdeal F (超平面) (x ≫ globalChartToProj i) = ( (s_0/s_i)(x) )`

★★★★すなわち「点 `x` での超平面の引き戻しは、**比の切断 `s_0/s_i` の値**が生成する」。

## ★★★機構 —— チャート射は `Spec` の射に潰れる

★`globalChartToProj i = globalChartMorphism i ≫ awayι (x_i)` であり、
`globalChartMorphism i` の**終域はアフィン**（`Spec A⁰_{x_i}`）なので、
点 `x : Spec 𝓞_F ⟶ X_{s_i}` を合成すると `Spec.map ρ` の形になる（`Spec.map_preimage`）。
★★あとは `§9-865` をそのまま当てるだけである。

## ★★★★★これで何が繋がったか

| 側 | 到達点 |
|---|---|
| 幾何 | `ample ⟹ 切断で生成 ⟹ チャート射 ⟹ 超平面の引き戻し = 比の切断` |
| 数論 | `超平面の引き戻し ⟹ deg_fin ⟹ ht = log H ⟹ Northcott` |

★**両側が同じ対象（点に沿った超平面の引き戻しイデアル）で出会った**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★`Spec.preimage` は後合成と可換である -/

/-- ★**`Spec.preimage (a ≫ Spec.map g) = g ≫ Spec.preimage a`**。 -/
theorem specPreimage_comp_specMap {A B C : CommRingCat.{0}}
    (a : Spec C ⟶ Spec B) (g : A ⟶ B) :
    Spec.preimage (a ≫ Spec.map g) = g ≫ Spec.preimage a := by
  apply Spec.map_injective
  rw [Spec.map_preimage, Spec.map_comp, Spec.map_preimage]

/-! ## ★★★★★★★★★★★★★点に沿った引き戻し -/

/-- ★★★★★★★★★★★★**点に沿った超平面の引き戻し**（チャート射に沿って）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★`globalChartMorphism` の終域がアフィンなので、点を合成すると `Spec.map ρ` の形になり、
`§9-865`（段 C2c-1）がそのまま当たる。 -/
theorem pullbackIdeal_hyperplane_globalChart (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme) :
    pullbackIdeal F (hyperplaneIdeal N ℤ) (xF ≫ globalChartToProj M hM φ s i)
      = Ideal.span {(Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)).hom
          (projCoord N ℤ i 0)} := by
  refine pullbackIdeal_hyperplane_point F N i
    (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)) _ ?_
  rw [Spec.map_preimage, globalChartToProj, Category.assoc]
  rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★**生成元は比の切断の値である**。

    `ρ(x_j/x_i) = (s_j/s_i)(x)`

★これで幾何側（`§9-885`）と数論側（`§9-865`）が**同じ対象で出会う**。 -/
theorem preimage_globalChartMorphism_projCoord (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme) :
    (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)).hom (projCoord N ℤ i j)
      = (Spec.preimage (xF ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
          (((nonVanishing M (s i)).topIso.inv.hom) (globalRatio M hM (s j) (s i))) := by
  rw [globalChartMorphism, ← Category.assoc, specPreimage_comp_specMap]
  show (Spec.preimage (xF ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
    (((nonVanishing M (s i)).topIso.inv.hom)
      (globalAwayHom M hM φ s i (projCoord N ℤ i j))) = _
  rw [globalAwayHom_projCoord]

/-- ★★★★★★★★★★★★★★**まとめ** —— 点に沿った超平面の引き戻しは `(s_0/s_i)(x)` が生成する。 -/
theorem pullbackIdeal_hyperplane_globalRatio (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme) :
    pullbackIdeal F (hyperplaneIdeal N ℤ) (xF ≫ globalChartToProj M hM φ s i)
      = Ideal.span {(Spec.preimage (xF ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
          (((nonVanishing M (s i)).topIso.inv.hom) (globalRatio M hM (s 0) (s i)))} := by
  rw [pullbackIdeal_hyperplane_globalChart, preimage_globalChartMorphism_projCoord]

/-! ## ★出典の紐付け(`.src`) -/

def specPreimage_comp_specMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(Spec.preimage は後合成と可換)",
    sectionId := "genell-prop-1-4" }

def pullbackIdeal_hyperplane_globalChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(点に沿った超平面の引き戻し——チャート射に沿って)",
    sectionId := "genell-prop-1-4" }

def pullbackIdeal_hyperplane_globalRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(点に沿った超平面の引き戻しは (s_0/s_i)(x) が生成する)",
    sectionId := "genell-prop-1-4" }

def pullbackIdeal_hyperplane_globalRatio.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullbackIdeal_hyperplane_point(段 C2c-1、§9-865)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdeal_hyperplane_point") 2,
    .citation "[ABC3]" "globalAwayHom_projCoord(x_j/x_i ↦ s_j/s_i、§9-842)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalAwayHom_projCoord") 2,
    .implicitStep
      ("★globalChartMorphism の**終域はアフィン**(Spec A⁰_{x_i})なので、" ++
       "点 x : Spec 𝓞_F ⟶ X_{s_i} を合成すると Spec.map ρ の形になる(Spec.map_preimage)。" ++
       "★★あとは §9-865 をそのまま当てるだけである") 2,
    .implicitStep
      ("★★★これで幾何側(ample ⟹ 切断で生成 ⟹ チャート射 ⟹ 超平面の引き戻し = 比の切断)と " ++
       "数論側(超平面の引き戻し ⟹ deg_fin ⟹ ht = log H ⟹ Northcott)が" ++
       "**同じ対象(点に沿った超平面の引き戻しイデアル)で出会った**") 3 ]

end ABC3.Found.GenEll
