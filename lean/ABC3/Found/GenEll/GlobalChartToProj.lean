/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalAwayHom
import ABC3.Found.GenEll.ProjectiveSpace
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★チャートから `ℙᴺ_R` への射（`⊓ V` なし）—— 段 E1d（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か

`§9-816` の `chartToProj` は `X_{s_i} **⊓ V** ⟶ ℙᴺ_R` という形だった
——`§9-841` で測ったとおり、そこが段 E3d と噛み合わない。

★本ファイルは `§9-842` の `globalAwayHom`（`⊓ V` なし）から同じ射を作り直す:

    `X_{s_i} ⟶ Spec Γ(X, X_{s_i}) ⟶ Spec A⁰_{x_i} ⟶ ℙᴺ_R`

★★これで**チャートは `X_{s_i}` そのもの**になり、`§9-847` の全射性がそのまま
「チャートで閉埋め込み」に化ける（段 E3d、`§9-834`）。

## ★★★機構 —— mathlib の 3 本を繋ぐだけ

| 道具 | 役割 |
|---|---|
| `Scheme.toSpecΓ` | `X_{s_i} ⟶ Spec Γ(X_{s_i})` |
| `Scheme.Opens.topIso` | `Γ(X_{s_i}, ⊤) ≅ Γ(X, X_{s_i})` |
| `Proj.awayι` | `Spec A⁰_{x_i} ⟶ Proj 𝒜`（開埋め込み、像はちょうど `D₊(x_i)`） |

## ★残っている段（明示）

★★**`chartToProjG i` が重なりで一致すること**が段 E1d の残りである。
`§9-836` の `glueOpens` に渡せば `ψ : X ⟶ ℙᴺ_R` が出る。
★★★機構は `sectionRatio_agree`（段 D3）と同じ「比の比は比」だが、
`Spec` への射の一致は**環準同型の一致**に落とすのが定石である
——`Spec A⁰_{x_i}` と `Spec A⁰_{x_j}` は別のチャートなので、
`D₊(x_i) ⊓ D₊(x_j)` へ両方を落として比べる段が要る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★チャートから `Spec A⁰_{x_i}` への射 -/

/-- ★★★★★★★★**チャートの射** `X_{s_i} ⟶ Spec A⁰_{x_i}`（`⊓ V` なし）。

★`§9-842` の `globalAwayHom` を `Spec` へ渡し、`Scheme.toSpecΓ` と繋ぐだけである。 -/
noncomputable def globalChartMorphism (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    (nonVanishing M (s i)).toScheme
      ⟶ Spec (CommRingCat.of (HomogeneousLocalization.Away
          (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i))) :=
  (nonVanishing M (s i)).toScheme.toSpecΓ
    ≫ Spec.map (CommRingCat.ofHom
        (((nonVanishing M (s i)).topIso.inv.hom).comp (globalAwayHom M hM φ s i)))

/-- ★★★★★★★★★★**チャートから `ℙᴺ_R` への射**（`⊓ V` なし）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★`Proj.awayι`（mathlib、`Spec A⁰_f ⟶ Proj 𝒜` の開埋め込み）と繋ぐだけである。
★★`§9-816` と違い**チャートは `X_{s_i}` そのもの**である。 -/
noncomputable def globalChartToProj (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    (nonVanishing M (s i)).toScheme ⟶ projSpace N R :=
  globalChartMorphism M hM φ s i
    ≫ Proj.awayι _ (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos

/-- ★★★★★**チャートの射の像は `D₊(x_i)` に入る**。

★`Proj.opensRange_awayι`（mathlib）そのものである。 -/
theorem range_globalChartToProj_le (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    Set.range (globalChartToProj M hM φ s i).base
      ⊆ (↑(Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
          (MvPolynomial.X i)) :
          Set (Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))) := by
  rintro y ⟨z, rfl⟩
  have h : (globalChartToProj M hM φ s i).base z
      = (Proj.awayι _ (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos).base
        ((globalChartMorphism M hM φ s i).base z) := rfl
  rw [h, ← Proj.opensRange_awayι _ (MvPolynomial.X i)
    (MvPolynomial.isHomogeneous_X R i) one_pos]
  exact ⟨_, rfl⟩

/-! ## ★★★★★★★★段 E3c から段 E3d への橋 -/

/-- ★★★★★★★★**チャートの射は閉埋め込みである**（座標が全射なら）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`X_{s_i}` がアフィンで `A⁰_{x_i} → Γ(X, X_{s_i})` が全射（`§9-847`）なら、
★`X_{s_i} ⟶ Spec A⁰_{x_i}` は**閉埋め込み**である。

★★機構は「アフィンなら `toSpecΓ` は同型」＋
`IsClosedImmersion.spec_of_surjective`（mathlib）だけである。

## ★測定の記録

★★★`haveI` で置いた `IsClosedImmersion (Spec.map …)` を
instance 解決が**拾わなかった**（2026-08-28 実測）。
`@IsClosedImmersion.comp … inferInstance hspec` と**明示的に渡す**と通る。 -/
theorem isClosedImmersion_globalChartMorphism (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i)))
    (hsurj : Function.Surjective (globalAwayHom M hM φ s i)) :
    IsClosedImmersion (globalChartMorphism M hM φ s i) := by
  haveI : IsAffine (nonVanishing M (s i)).toScheme := haff
  have hiso : Function.Surjective ((nonVanishing M (s i)).topIso.inv.hom) :=
    (ConcreteCategory.bijective_of_isIso (nonVanishing M (s i)).topIso.inv).2
  have hspec : IsClosedImmersion (Spec.map (CommRingCat.ofHom
      (((nonVanishing M (s i)).topIso.inv.hom).comp (globalAwayHom M hM φ s i)))) := by
    refine IsClosedImmersion.spec_of_surjective _ ?_
    show Function.Surjective (((nonVanishing M (s i)).topIso.inv.hom).comp
      (globalAwayHom M hM φ s i))
    rw [RingHom.coe_comp]
    exact hiso.comp hsurj
  exact @IsClosedImmersion.comp _ _ _ (nonVanishing M (s i)).toScheme.toSpecΓ
    (Spec.map (CommRingCat.ofHom
      (((nonVanishing M (s i)).topIso.inv.hom).comp (globalAwayHom M hM φ s i))))
    inferInstance hspec

/-- ★★★★★★★★★**したがってチャートの射は埋め込みである** —— 原典の「embedding」。

★閉埋め込み（本ファイル）のあとに開埋め込み `Proj.awayι` を合成するので、
全体は **immersion** になる。 -/
theorem isImmersion_globalChartToProj (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i)))
    (hsurj : Function.Surjective (globalAwayHom M hM φ s i)) :
    IsImmersion (globalChartToProj M hM φ s i) := by
  haveI := isClosedImmersion_globalChartMorphism M hM φ s i haff hsurj
  exact inferInstanceAs (IsImmersion (globalChartMorphism M hM φ s i ≫
    Proj.awayι _ (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos))

/-! ## ★出典の紐付け(`.src`) -/

def globalChartMorphism.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの射 X_{s_i} ⟶ Spec A⁰_{x_i}——⊓ V なし)",
    sectionId := "genell-prop-1-4" }

def globalChartToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートから ℙᴺ_R への射——⊓ V なし)",
    sectionId := "genell-prop-1-4" }

def range_globalChartToProj_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの射の像は D₊(x_i) に入る——⊓ V なし)",
    sectionId := "genell-prop-1-4" }

def isClosedImmersion_globalChartMorphism.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの射は閉埋め込みである)",
    sectionId := "genell-prop-1-4" }

def isImmersion_globalChartToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの射は埋め込みである)",
    sectionId := "genell-prop-1-4" }

def globalChartToProj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "globalAwayHom(大域のチャート写像、§9-842)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalAwayHom") 2,
    .citation "[mathlib]" "Scheme.toSpecΓ / Scheme.Opens.topIso / Proj.awayι"
      (.inMathlib "AlgebraicGeometry.Proj.awayι") 2,
    .implicitStep
      ("★§9-816 の chartToProj は X_{s_i} **⊓ V** ⟶ ℙᴺ_R という形で、" ++
       "§9-841 で測ったとおり段 E3d と噛み合わなかった。" ++
       "★★本ファイルは §9-842 の globalAwayHom から同じ射を作り直したもので、" ++
       "**チャートは X_{s_i} そのもの**になる") 4,
    .implicitStep
      ("★★★**chartToProjG i が重なりで一致すること**が段 E1d の残りである。" ++
       "§9-836 の glueOpens に渡せば ψ : X ⟶ ℙᴺ_R が出る。" ++
       "★機構は sectionRatio_agree(段 D3)と同じ「比の比は比」だが、" ++
       "Spec A⁰_{x_i} と Spec A⁰_{x_j} は別のチャートなので " ++
       "D₊(x_i) ⊓ D₊(x_j) へ両方を落として比べる段が要る") 7 ]

end ABC3.Found.GenEll
