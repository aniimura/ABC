/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalToProj
import ABC3.Found.GenEll.BasicOpenRatio
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★`ψ : X ⟶ ℙᴺ_R` は埋め込みである（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★★★★これは何か —— 原文の「embedding」

`§9-911` で `ψ : X ⟶ ℙᴺ_R` が貼れた。★本ファイルはそれが**埋め込みである**ことを示す。

★★機構は「`IsImmersion` は target について局所的である」
（mathlib の `instance : IsZariskiLocalAtTarget @IsImmersion`）ことである:

| 段 | 内容 |
|---|---|
| 1 | `ψ⁻¹(D₊(x_i)) = X_{s_i}` |
| 2 | `ψ ∣_ D₊(x_i) ≫ ι = X_{s_i}.ι ≫ ψ = globalChartToProj i`（`§9-848` で埋め込み） |
| 3 | `IsImmersion.of_comp` で `ψ ∣_ D₊(x_i)` が埋め込み |
| 4 | `IsZariskiLocalAtTarget.of_iSup_eq_top` で `ψ` が埋め込み |

## ★★★段 1 の計算

`globalChartToProj j ⁻¹ᵁ D₊(x_i)` を 3 つの引き戻しに分ける:

* `Proj.awayι_preimage_basicOpen`: `awayι(x_j)⁻¹ D₊(x_i) = D(x_i/x_j)`
* `SpecMap_preimage_basicOpen`: `Spec.map φ ⁻¹ D(r) = D(φ r)`
* `Scheme.toSpecΓ_preimage_basicOpen`: `toSpecΓ ⁻¹ D(a) = basicOpen a`

★`globalAwayHom … (x_i/x_j) = s_i/s_j`（`§9-842`）なので、
結果は **`X_{s_j}` の中の `s_i/s_j` の非零点**になる。
★★`§9-912` の `basicOpen_globalRatio`（`X.basicOpen(s/t) = X_s ⊓ X_t`）で
それが `X_{s_i} ⊓ X_{s_j}` だと分かる。
★★★`⨆_j X_{s_j} = ⊤` で貼れば `ψ⁻¹(D₊(x_i)) = X_{s_i}` である。

## ★残っている仮定（明示）

★★`haff`（チャートがアフィン）と `hsurj`（チャート写像が全射）は
`§9-848` の `isImmersion_globalChartToProj` がそのまま要求するものであり、
**段 E3 の内容**である（`§9-833`・`§9-847` が道具を与えている）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★段 1 —— チャートごとの引き戻し -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★**`globalChartToProj j ⁻¹ᵁ D₊(x_i)` は `s_i/s_j` の非零点である**。

★引き戻しを `awayι`・`Spec.map`・`toSpecΓ` の 3 段に分けるだけである。 -/
theorem globalChartToProj_preimage_basicOpen (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) :
    globalChartToProj M hM φ s j ⁻¹ᵁ
        (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i))
      = (nonVanishing M (s j)).toScheme.basicOpen
          ((nonVanishing M (s j)).topIso.inv (globalRatio M hM (s i) (s j))) := by
  rw [globalChartToProj, globalChartMorphism, Scheme.Hom.comp_preimage,
    Scheme.Hom.comp_preimage,
    Proj.awayι_preimage_basicOpen _ _ _ (MvPolynomial.isHomogeneous_X R i) one_pos,
    isLocalizationElem_eq_projCoord', SpecMap_preimage_basicOpen,
    Scheme.toSpecΓ_preimage_basicOpen]
  show (nonVanishing M (s j)).toScheme.basicOpen
    ((nonVanishing M (s j)).topIso.inv (globalAwayHom M hM φ s j (projCoord N R j i))) = _
  rw [globalAwayHom_projCoord]

/-! ## ★★★★★★★★★★★★段 1 —— 大域の引き戻し -/

/-- ★★★★★★★★★★★★**`ψ⁻¹(D₊(x_i)) = X_{s_i}`**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★これが「`ψ` は `s_i` が消えない所をちょうど `i` 番目のチャートに送る」ことである。
★★`§9-912` の `basicOpen_globalRatio` と `⨆_j X_{s_j} = ⊤` で貼る。 -/
theorem globalToProj_preimage_basicOpen (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1)) :
    globalToProj M hM φ s hcov ⁻¹ᵁ
        (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i))
      = nonVanishing M (s i) := by
  set A : X.Opens := globalToProj M hM φ s hcov ⁻¹ᵁ
    (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)) with hA
  have key : ∀ j, A ⊓ nonVanishing M (s j)
      = nonVanishing M (s i) ⊓ nonVanishing M (s j) := by
    intro j
    have h2 : (nonVanishing M (s j)).ι ⁻¹ᵁ A
        = globalChartToProj M hM φ s j ⁻¹ᵁ
          (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
            (MvPolynomial.X i)) := by
      rw [hA, ← Scheme.Hom.comp_preimage, ι_globalToProj]
    calc A ⊓ nonVanishing M (s j)
        = (nonVanishing M (s j)).ι ''ᵁ ((nonVanishing M (s j)).ι ⁻¹ᵁ A) := by
          rw [Scheme.Hom.image_preimage_eq_opensRange_inf, Scheme.Opens.opensRange_ι, inf_comm]
      _ = (nonVanishing M (s j)).ι ''ᵁ ((nonVanishing M (s j)).toScheme.basicOpen
            ((nonVanishing M (s j)).topIso.inv (globalRatio M hM (s i) (s j)))) := by
          rw [h2, globalChartToProj_preimage_basicOpen]
      _ = nonVanishing M (s i) ⊓ nonVanishing M (s j) := image_basicOpen_globalRatio M hM _ _
  calc A = A ⊓ ⊤ := (inf_top_eq A).symm
    _ = A ⊓ ⨆ k, nonVanishing M (s k) := by rw [hcov]
    _ = ⨆ k, A ⊓ nonVanishing M (s k) := inf_iSup_eq _ _
    _ = ⨆ k, nonVanishing M (s i) ⊓ nonVanishing M (s k) := iSup_congr key
    _ = nonVanishing M (s i) ⊓ ⨆ k, nonVanishing M (s k) := (inf_iSup_eq _ _).symm
    _ = nonVanishing M (s i) ⊓ ⊤ := by rw [hcov]
    _ = nonVanishing M (s i) := inf_top_eq _

/-! ## ★★★★★★★★★★★★★★★★★段 2-4 —— 埋め込みであること -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★**`ψ` は埋め込みである** —— 原文の「embedding」。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★`IsImmersion` は target について局所的（mathlib の
`instance : IsZariskiLocalAtTarget @IsImmersion`）なので、
`D₊(x_i)` の上で見れば `globalChartToProj i`（`§9-848` で埋め込み）に戻る。

★★仮定 `haff`（チャートがアフィン）・`hsurj`（チャート写像が全射）は
`§9-848` がそのまま要求するもので、**段 E3 の内容**である。 -/
theorem isImmersion_globalToProj (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    (hsurj : ∀ i, Function.Surjective (globalAwayHom M hM φ s i)) :
    IsImmersion (globalToProj M hM φ s hcov) := by
  refine IsZariskiLocalAtTarget.of_iSup_eq_top
    (fun i : Fin (N + 1) => Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (MvPolynomial.X i)) (iSup_basicOpen_X_eq_top (Fin (N + 1)) R) ?_
  intro i
  haveI : IsImmersion (globalToProj M hM φ s hcov ∣_
      (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i))
      ≫ (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i)).ι) := by
    rw [morphismRestrict_ι, globalToProj_preimage_basicOpen M hM φ s hcov i, ι_globalToProj]
    exact isImmersion_globalChartToProj M hM φ s i (haff i) (hsurj i)
  exact IsImmersion.of_comp _
    (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)).ι

/-! ## ★出典の紐付け(`.src`) -/

def globalChartToProj_preimage_basicOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャート射に沿った D₊(x_i) の引き戻し)",
    sectionId := "genell-prop-1-4" }

def globalToProj_preimage_basicOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(ψ⁻¹(D₊(x_i)) = X_{s_i})",
    sectionId := "genell-prop-1-4" }

def isImmersion_globalToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(ψ : X ⟶ ℙᴺ_R は埋め込みである)",
    sectionId := "genell-prop-1-4" }

def isImmersion_globalToProj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsZariskiLocalAtTarget @IsImmersion(埋め込みは target について局所的)"
      (.inMathlib "AlgebraicGeometry.IsZariskiLocalAtTarget") 1,
    .citation "[mathlib]" "IsImmersion.of_comp・morphismRestrict_ι"
      (.inMathlib "AlgebraicGeometry.IsImmersion.of_comp") 1,
    .citation "[ABC3]" "isImmersion_globalChartToProj(チャートごとの埋め込み、§9-848)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isImmersion_globalChartToProj") 2,
    .citation "[ABC3]" "basicOpen_globalRatio(比が消えない所は X_s ∩ X_t、§9-912)"
      (.inProject "ABC3" "ABC3.Found.GenEll.basicOpen_globalRatio") 2,
    .citation "[ABC3]" "globalToProj(貼った射、§9-911)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalToProj") 3,
    .implicitStep
      ("★★仮定 haff(チャートがアフィン)と hsurj(チャート写像が全射)は " ++
       "§9-848 の isImmersion_globalChartToProj がそのまま要求するものであり、" ++
       "**段 E3 の内容**である(§9-833・§9-847 が道具を与えている)") 4,
    .implicitStep
      ("★★★これで原文の「[some positive tensor power of] the ample line bundle L_ℚ " ++
       "yields an embedding」の**射と埋め込み性**が揃った。" ++
       "残るのは (1) 被覆条件 ⨆_i X_{s_i} = ⊤ を ample から出すこと(段 E2)、" ++
       "(2) haff・hsurj を ample から出すこと(段 E3)、" ++
       "(3) 高さの側で D^{⊗n} = ψ^*(超平面) を言うこと") 6 ]

end ABC3.Found.GenEll
