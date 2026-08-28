/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottCoords
import Mathlib.AlgebraicGeometry.IdealSheaf.Subscheme
import Mathlib.AlgebraicGeometry.ProjectiveSpectrum.Proper
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★射影埋め込みから `ℤ`-モデルを作る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★★★★★★★★★★★★★★★★これは何か —— `Remark 1.4.1` 第 2 文の道具

`Remark141.lean` の末尾に記録したとおり、`Remark 1.4.1` に残っていたのは**第 2 文**

> this theory may be applied to any normal, projective scheme `Y` over `ℚ`

であり、それは**射影埋め込みから `ℤ`-モデルを作る**段を要求する。

★★★★本ファイルはその段を取る:

    `e : Y ⟶ ℙᴺ_ℤ`  ⟹  **スキーム論的像** `e.image ⊆ ℙᴺ_ℤ` は
    `Spec ℤ` 上**固有**であり、分離的かつ準コンパクトである

★★`Y` が `ℚ` 上の射影的スキームなら `Y ⟶ ℙᴺ_ℚ ⟶ ℙᴺ_ℤ` の像がその `ℤ`-モデルである。

## ★★★機構 —— mathlib のスキーム論的像

* `Scheme.Hom.image = f.ker.subscheme`（`f` の核イデアル層が切る閉部分スキーム）
* `Scheme.Hom.imageι` は**閉埋め込み**（インスタンス）
* `Proj.toSpecZero` は**固有**（mathlib、`Proj` の固有性）
* ★閉埋め込みは固有だから、合成も固有

★★準コンパクト性は `QuasiCompact.compactSpace_of_compactSpace` で
`Spec ℤ` のコンパクト性から降りる。

## ★これで何が繋がるか

★★★`§9-920`（`ample` なら `ψ : X ⟶ ℙᴺ` が埋め込み）と合わせると

    **豊富な直線束を持つスキームは `Spec ℤ` 上固有な `ℤ`-モデルを持つ**

が出る。★これは `Proposition 1.6`（`Prop16Proper.lean`）が仮定として受けている
固有性を**供給する側**でもある。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★★★★★★★★★★★★★スキーム論的像としての `ℤ`-モデル -/

/-- ★★★★★★★★★★★★★★**射影埋め込みの `ℤ`-モデル** —— スキーム論的像。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

★`Y ⟶ ℙᴺ_ℤ` の像は `ℙᴺ_ℤ` の閉部分スキームであり、これが `ℤ`-モデルである。 -/
noncomputable abbrev zModel {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) : Scheme.{0} :=
  e.image

/-- ★**`ℤ`-モデルの `ℙᴺ_ℤ` への閉埋め込み**。 -/
noncomputable abbrev zModelι {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    zModel e ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) :=
  e.imageι

/-- ★**もとのスキームから `ℤ`-モデルへの射**。 -/
noncomputable abbrev toZModel {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) : Y ⟶ zModel e :=
  e.toImage

/-- ★**分解する** —— `Y ⟶ ℤ-モデル ⟶ ℙᴺ_ℤ`。 -/
theorem toZModel_zModelι {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    toZModel e ≫ zModelι e = e :=
  Scheme.Hom.toImage_imageι e

/-- ★★**`ℤ`-モデルは `ℙᴺ_ℤ` の閉部分スキームである**。 -/
theorem isClosedImmersion_zModelι {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    IsClosedImmersion (zModelι e) :=
  inferInstance

/-! ## ★★★★★★★★★★★★★★★★★★固有性 -/

/-- ★★★★★★★★★★★★★★★★★★**`ℤ`-モデルは `Spec ℤ` 上固有である**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

★`Proj` の固有性（mathlib）＋「閉埋め込みは固有」＋「固有射の合成は固有」。
★★これが `Remark 1.4.1` 第 2 文（`ℚ` 上の射影的スキームに理論が適用できる）の
**`ℤ`-モデルの側**である。 -/
theorem isProper_zModel {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    IsProper (zModelι e ≫
      Proj.toSpecZero (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :=
  inferInstance

/-- ★★★**`ℤ`-モデルは分離的である**。

★閉埋め込みは分離的で、`ℙᴺ_ℤ` も分離的（mathlib）だから、合成が分離的である。 -/
theorem isSeparated_zModel {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    Scheme.IsSeparated (zModel e) := by
  refine ⟨?_⟩
  have h : terminal.from (zModel e) = zModelι e ≫ terminal.from
      (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) := Subsingleton.elim _ _
  rw [h]
  infer_instance

/-- ★★★**`ℤ`-モデルの台空間はコンパクトである**。

★`Spec ℤ` がコンパクトで、`ℤ`-モデル ⟶ `Spec ℤ` が準コンパクト（固有性の一部）だから。
★★これが `Proposition 1.6`（`Prop16Proper.lean`）が仮定として受けている
`[CompactSpace X]` を**供給する側**である。 -/
theorem compactSpace_zModel {Y : Scheme.{0}} {N : ℕ}
    (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    CompactSpace (zModel e : Type) := by
  haveI : QuasiCompact (zModelι e ≫ Proj.toSpecZero
    (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) := inferInstance
  exact QuasiCompact.compactSpace_of_compactSpace
    (zModelι e ≫ Proj.toSpecZero (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))

/-! ## ★出典の紐付け(`.src`) -/

def zModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(射影埋め込みの ℤ-モデル——スキーム論的像)",
    sectionId := "genell-rem-1-4-1" }

def toZModel_zModelι.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(Y ⟶ ℤ-モデル ⟶ ℙᴺ_ℤ と分解する)",
    sectionId := "genell-rem-1-4-1" }

def isClosedImmersion_zModelι.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℤ-モデルは ℙᴺ_ℤ の閉部分スキームである)",
    sectionId := "genell-rem-1-4-1" }

def isProper_zModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℤ-モデルは Spec ℤ 上固有である)",
    sectionId := "genell-rem-1-4-1" }

def isSeparated_zModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℤ-モデルは分離的である)",
    sectionId := "genell-rem-1-4-1" }

def compactSpace_zModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℤ-モデルの台空間はコンパクトである)",
    sectionId := "genell-rem-1-4-1" }

def isProper_zModel.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.image(スキーム論的像)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.image") 2,
    .citation "[mathlib]" "Proj.toSpecZero の固有性"
      (.inMathlib "AlgebraicGeometry.Proj.toSpecZero") 3,
    .citation "[mathlib]" "QuasiCompact.compactSpace_of_compactSpace"
      (.inMathlib "AlgebraicGeometry.QuasiCompact.compactSpace_of_compactSpace") 1,
    .implicitStep
      ("★★★★測定(2026-08-29): Remark141.lean の末尾に「残るのは原文第 2 文の側だけ" ++
       "——射影埋め込みから ℤ-モデルを作る段」と記録していたが、" ++
       "**mathlib のスキーム論的像でそのまま取れる**。" ++
       "e.imageι は閉埋め込み(インスタンス)、Proj は Spec ℤ 上固有、" ++
       "閉埋め込みは固有だから合成も固有である") 5,
    .implicitStep
      ("★★§9-920(ample なら ψ : X ⟶ ℙᴺ が埋め込み)と合わせると" ++
       "「豊富な直線束を持つスキームは Spec ℤ 上固有な ℤ-モデルを持つ」が出る。" ++
       "★これは Proposition 1.6(Prop16Proper.lean)が仮定として受けている" ++
       "固有性・コンパクト性を供給する側でもある") 4 ]

end ABC3.Found.GenEll
