/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalAwayHom
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★点の水準での段 E3d —— 座標が点を分ける（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★測定 —— 段 E1d は Prop 1.4 (iv) の**必須路ではない**

`Found/GenEll/NorthcottCoord.lean` の `northcott_of_projModel` が受け取るのは

| 受けるもの | 中身 |
|---|---|
| `crd p : ι → fld p` | ★**同次座標** |
| `idx : ι` | 割る成分 |
| `hinj` | 正規化座標 `crd p j / crd p idx` が**単射** |

であって、**スキームの射 `ψ : X ⟶ ℙᴺ_R` そのものではない**（2026-08-28 実測）。

★したがって段 E1d（チャートの射を貼って `ψ` を作る）は
**Prop 1.4 (iv) の消費側から見れば迂回できる**——要るのは「座標が点を分ける」ことだけである。

★★本ファイルはその「点を分ける」を型にする:

> `A⁰_{x_i} → Γ(X, X_{s_i})` が全射なら（`§9-847`）、
> `X_{s_i}` の 2 つの `Spec F`-点は **`A⁰_{x_i}` の像での値が一致すれば等しい**。

## ★★★機構 —— 2 本

| 道具 | 役割 |
|---|---|
| `§9-847` の全射性 | `Γ(X_{s_i}, ⊤)` は `A⁰_{x_i}` の像で尽きる |
| `ext_of_isAffine`（mathlib） | アフィンなら射は `Γ` の準同型で決まる |

## ★残っている段（明示）

★★★**`A⁰_{x_i}` が `R` 上 `x_j/x_i` で生成されること**が要る
——本ファイルは「`A⁰_{x_i}` **全体**での一致」を仮定にしているが、
`northcott_of_projModel` が渡すのは**座標 `x_j/x_i` での一致**だけだからである。
★これは `Proj` の標準チャート `D₊(x_i) ≅ 𝔸ᴺ` という古典的事実だが、
mathlib に見当たらない（2026-08-28 実測）。★★純粋に代数の主張なので、
段 E1d（射の貼り合わせ）より**小さい穴**である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★チャートの座標環から `Γ(X_{s_i}, ⊤)` へ -/

/-- ★**チャートの座標環から `Γ(X_{s_i}, ⊤)` への環準同型**。

★`§9-842` の `globalAwayHom` に `topIso`（開部分スキームの `Γ(⊤)` と `Γ(X, U)` の同型）を繋いだもの。 -/
noncomputable def chartRingHom (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)
      →+* (Γ((nonVanishing M (s i)).toScheme, (⊤ : (nonVanishing M (s i)).toScheme.Opens))
        : Type) :=
  ((nonVanishing M (s i)).topIso.inv.hom).comp (globalAwayHom M hM φ s i)

/-- ★★**`§9-847` の全射性はそのまま移る**（`topIso` は同型だから）。 -/
theorem surjective_chartRingHom (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (hsurj : Function.Surjective (globalAwayHom M hM φ s i)) :
    Function.Surjective (chartRingHom M hM φ s i) := by
  have hiso : Function.Surjective ((nonVanishing M (s i)).topIso.inv.hom) :=
    (ConcreteCategory.bijective_of_isIso (nonVanishing M (s i)).topIso.inv).2
  show Function.Surjective (((nonVanishing M (s i)).topIso.inv.hom).comp
    (globalAwayHom M hM φ s i))
  rw [RingHom.coe_comp]
  exact hiso.comp hsurj

/-! ## ★★★★★★★★★座標が点を分ける -/

/-- ★★★★★★★★★**座標が点を分ける** —— 点の水準での段 E3d。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`X_{s_i}` がアフィンで `A⁰_{x_i} → Γ(X, X_{s_i})` が全射（`§9-847`）なら、
★`X_{s_i}` の 2 つの `Spec F`-点は **`A⁰_{x_i}` の像での値が一致すれば等しい**。

★★これが `northcott_of_projModel`（`Found/GenEll/NorthcottCoord.lean`）が要求する
`hinj`（正規化座標が単射）の中身である
——**スキームの射 `ψ : X ⟶ ℙᴺ_R` を作らなくてよい**。

★★★機構は `§9-847` の全射性と `ext_of_isAffine`（mathlib）だけである。 -/
theorem ext_of_chart (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i)))
    (hsurj : Function.Surjective (globalAwayHom M hM φ s i))
    {F : Type} [CommRing F]
    (p q : Spec (CommRingCat.of F) ⟶ (nonVanishing M (s i)).toScheme)
    (h : ∀ z, p.appTop.hom (chartRingHom M hM φ s i z)
      = q.appTop.hom (chartRingHom M hM φ s i z)) : p = q := by
  haveI : IsAffine (nonVanishing M (s i)).toScheme := haff
  refine ext_of_isAffine ?_
  refine CommRingCat.hom_ext (RingHom.ext (fun y => ?_))
  obtain ⟨z, rfl⟩ := surjective_chartRingHom M hM φ s i hsurj y
  exact h z

/-! ## ★出典の紐付け(`.src`) -/

def chartRingHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの座標環から Γ(X_{s_i}, ⊤) への環準同型)",
    sectionId := "genell-prop-1-4" }

def surjective_chartRingHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(全射性は topIso で移る)",
    sectionId := "genell-prop-1-4" }

def ext_of_chart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(座標が点を分ける——点の水準での段 E3d)",
    sectionId := "genell-prop-1-4" }

def ext_of_chart.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_coords_surjective(段 E3c、§9-847)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_coords_surjective") 2,
    .citation "[mathlib]" "ext_of_isAffine(アフィンなら射は Γ の準同型で決まる)"
      (.inMathlib "AlgebraicGeometry.ext_of_isAffine") 2,
    .citation "[ABC3]" "northcott_of_projModel(消費側——座標と単射性しか受け取らない)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_projModel") 2,
    .implicitStep
      ("★★★★★**測定**: northcott_of_projModel が受け取るのは同次座標と単射性であって、" ++
       "**スキームの射 ψ : X ⟶ ℙᴺ_R そのものではない**(2026-08-28 実測)。" ++
       "★したがって段 E1d(チャートの射を貼って ψ を作る)は" ++
       "**Prop 1.4 (iv) の消費側から見れば迂回できる**") 8,
    .implicitStep
      ("★★★**A⁰_{x_i} が R 上 x_j/x_i で生成されること**が残る" ++
       "——本ファイルは A⁰_{x_i} 全体での一致を仮定にしているが、" ++
       "northcott_of_projModel が渡すのは座標 x_j/x_i での一致だけだからである。" ++
       "★これは Proj の標準チャート D₊(x_i) ≅ 𝔸ᴺ という古典的事実だが mathlib に見当たらない" ++
       "(2026-08-28 実測)。★★純粋に代数の主張なので、段 E1d より**小さい穴**である") 6 ]

end ABC3.Found.GenEll
