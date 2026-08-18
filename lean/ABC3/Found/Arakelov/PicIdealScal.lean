import ABC3.Found.Arakelov.PicIdealAway

/-!
# Arakelov (B2) 第 157 ブロック —— **係数を `Γ(X, D(f))` へ広げる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★第 125 の `idealSections` 版

    idealAwayEquivScalar : I_f ≃ₗ[Γ(X, D f)] idealSections D (D f)

★`tilde` 版(第 125)では `Module` と `IsScalarTower` を**別の経路**から取る必要があったが、
イデアル層では素直に通った——台が `Γ(X, D f)` の部分加群だからである。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `modOnLocalizedXInst` | ★instance 化 |
| `towerOnLocalizedXInst` | ★足場(`of_algebraMap_smul` + `AlgEquiv.commutes`) |
| `idealAwayEquivScalar` | ★★★係数を広げた同型 |

★★★**一発で通った**——第 125 が 20 回以上詰まったのと対照的である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData) (A : X.affineOpens)

noncomputable instance modOnLocalizedXInst (f : (Γ(X, A.1) : Type u)) :
    Module (Γ(X, X.basicOpen f) : Type u)
      (LocalizedModule (Submonoid.powers f) (D.ideal A)) :=
  modOnLocalizedX D A f

instance towerOnLocalizedXInst (f : (Γ(X, A.1) : Type u)) :
    letI := resAlg A f
    IsScalarTower (Γ(X, A.1) : Type u) (Γ(X, X.basicOpen f) : Type u)
      (LocalizedModule (Submonoid.powers f) (D.ideal A)) := by
  letI := resAlg A f
  refine IsScalarTower.of_algebraMap_smul fun r x => ?_
  show ((awayRingEquivX A f) (algebraMap (Γ(X, A.1) : Type u)
    (Γ(X, X.basicOpen f) : Type u) r)) • x = r • x
  rw [AlgEquiv.commutes, IsScalarTower.algebraMap_smul]

noncomputable def idealAwayEquivScalar (f : (Γ(X, A.1) : Type u)) :
    LocalizedModule (Submonoid.powers f) (D.ideal A)
      ≃ₗ[(Γ(X, X.basicOpen f) : Type u)] (idealSections D (X.basicOpen f)) :=
  letI := resAlg A f
  haveI := isLocalizationAway_bo A f
  haveI := towerOnLocalizedXInst D A f
  (idealAwayEquiv D A f).extendScalarsOfIsLocalization (Submonoid.powers f)
    (Γ(X, X.basicOpen f) : Type u)


/-! ## ★出典の紐付け(`.src`) -/

def idealAwayEquivScalar.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——係数を Γ(X,D f) へ広げる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
