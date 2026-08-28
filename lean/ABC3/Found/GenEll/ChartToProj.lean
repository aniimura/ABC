/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayHom
import ABC3.Found.GenEll.ProjectiveSpace
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★チャートから `ℙᴺ_R` への射 —— 段 E1d の前半（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か

`§9-815` で環準同型

    `A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)`

が入った。★本ファイルはそれを**射**に直す:

    `X_{s_i} ⊓ V ⟶ Spec A⁰_{x_i} ⟶ ℙᴺ_R`

## ★★★★★機構 —— 2 本の在庫を繋ぐだけ

| 段 | 道具 |
|---|---|
| 環準同型 ⟶ `Spec` への射 | `Y.toSpecΓ ≫ Spec.map (…)`（mathlib） |
| `Γ(U.toScheme, ⊤) ≅ Γ(X, U)` | `Scheme.Opens.topIso`（mathlib） |
| `Spec A⁰_{x_i} ⟶ Proj 𝒜` | `Proj.awayι`（mathlib、開埋め込み） |

★`toSpecΓ` の行き先は `Spec Γ(Y, ⊤)` なので、`topIso.inv` で
`Γ(X, X_{s_i} ⊓ V)` から `Γ(Y, ⊤)` へ移してから `Spec.map` する。

## ★残っている段（明示）

★★**貼り合わせは本ファイルに無い**。要るのは 2 つ:

1. **重なりでの一致**——`X_{s_i} ⊓ X_{s_j}` の上で `chartToProj i` と `chartToProj j` が一致すること。
   ★機構は `sectionRatio_agree`（段 D3）と同じ「比の比は比」である。
2. **被覆**——`⋃_i X_{s_i} = X`（段 E2 の前半、`IsAmple` の準コンパクト性から）。

★★★この 2 つが済めば `X ⟶ ℙᴺ_R` が貼り合い、段 E1 が閉じる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-- ★★★★★★★★★**チャートの射** `X_{s_i} ⊓ V ⟶ Spec A⁰_{x_i}`。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★`toSpecΓ` の行き先は `Spec Γ(Y, ⊤)` なので、`topIso.inv` で
`Γ(X, X_{s_i} ⊓ V)` から `Γ(Y, ⊤)` へ移してから `Spec.map` する。 -/
noncomputable def chartMorphism (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    (nonVanishing M (s i) ⊓ V).toScheme
      ⟶ Spec (CommRingCat.of (HomogeneousLocalization.Away
          (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i))) :=
  (nonVanishing M (s i) ⊓ V).toScheme.toSpecΓ
    ≫ Spec.map (CommRingCat.ofHom
        (((nonVanishing M (s i) ⊓ V).topIso.inv.hom).comp (awayToSectionHom M V e φ s i)))

/-- ★★★★★★★★★★**チャートから `ℙᴺ_R` への射**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★`Proj.awayι`（mathlib、`Spec A⁰_f ⟶ Proj 𝒜` の開埋め込み）と繋ぐだけである。 -/
noncomputable def chartToProj (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    (nonVanishing M (s i) ⊓ V).toScheme ⟶ projSpace N R :=
  chartMorphism M V e φ s i
    ≫ Proj.awayι _ (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos

/-- ★★★★★**チャートの射の像は `D₊(x_i)` に入る**。

★`Proj.opensRange_awayι`（mathlib、`awayι` の像はちょうど `D₊(f)`）そのものである。
★★貼り合わせ（段 E1d の残り）で「どのチャートがどこを覆うか」を言うのに要る。

★★★型の記録: `projSpace N R` は `def`（semireducible）なので
`(projSpace N R).Opens` では `Proj.basicOpen` と**構文的に合わない**。
`Proj (…)` の綴りで書くと通る（2026-08-28 実測）。 -/
theorem range_chartToProj_le (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    Set.range (chartToProj M V e φ s i).base
      ⊆ (↑(Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
          (MvPolynomial.X i)) :
          Set (Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))) := by
  rintro y ⟨z, rfl⟩
  have h : (chartToProj M V e φ s i).base z
      = (Proj.awayι _ (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos).base
        ((chartMorphism M V e φ s i).base z) := rfl
  rw [h, ← Proj.opensRange_awayι _ (MvPolynomial.X i)
    (MvPolynomial.isHomogeneous_X R i) one_pos]
  exact ⟨_, rfl⟩

/-! ## ★出典の紐付け(`.src`) -/

def range_chartToProj_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの射の像は D₊(x_i) に入る)",
    sectionId := "genell-prop-1-4" }

def chartMorphism.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの射 X_{s_i} ⊓ V ⟶ Spec A⁰_{x_i})",
    sectionId := "genell-prop-1-4" }

def chartToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートから ℙᴺ_R への射。貼り合わせは含まない)",
    sectionId := "genell-prop-1-4" }

def chartToProj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "awayToSectionHom(環準同型 A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)、§9-815)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayToSectionHom") 7,
    .citation "[mathlib]" "Scheme.toSpecΓ / Scheme.Opens.topIso"
      (.inMathlib "AlgebraicGeometry.Scheme.toSpecΓ") 7,
    .citation "[mathlib]" "Proj.awayι(Spec A⁰_f ⟶ Proj 𝒜 の開埋め込み)"
      (.inMathlib "AlgebraicGeometry.Proj.awayι") 7,
    .implicitStep
      ("★**貼り合わせは本ファイルに無い**。要るのは (1) 重なり X_{s_i} ⊓ X_{s_j} での一致" ++
       "(機構は sectionRatio_agree(段 D3)と同じ「比の比は比」)と、" ++
       "(2) 被覆 ⋃_i X_{s_i} = X(段 E2 の前半、IsAmple の準コンパクト性から)である。" ++
       "★★この 2 つが済めば X ⟶ ℙᴺ_R が貼り合い、段 E1 が閉じる") 7 ]

end ABC3.Found.GenEll
