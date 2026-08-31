/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Hyperplane
import ABC3.Found.GenEll.AmpleDef
import Mathlib.AlgebraicGeometry.IdealSheaf.Functorial
import Mathlib.AlgebraicGeometry.Morphisms.Immersion
import ABC3.Meta.Claim

/-!
# ★★★★★★★**very ample な因子**の定義（因子表示）（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か

原文 `Proposition 1.4, (iv)` の証明は「ample な `L_ℚ` のある正冪が**埋め込み**
`ϵ_ℚ : X_ℚ ⊆ ℙ_ℚ` を与える」で始まる。★その「埋め込み」を型にするのが本ファイルである。

★★**因子表示**（本プロジェクトが §1 を通して採っている表示）では、
`O(1)` の層を作らずに済む——`O(1)` は**超平面因子**だからである
（`Hyperplane.lean` の `hyperplaneIdeal`）。

    `D` が very ample  ⟺  `∃ N` と immersion `i : X ⟶ ℙᴺ_R` で `D = i^*(超平面因子)`

★★★`Scheme.IdealSheafData.comap`（mathlib）が引き戻しを与える。

## ★★★★★台帳との対応

`ResearchPaper/mathlib-gap.json` の `ample-and-projective-embedding` の段 **C2c** は

> 残るのは「超平面因子の高さが素朴高さであること」と、
> very ample（`n·D` が閉埋め込みによる超平面の引き戻し）の**定義**である。

と書いていた。★**本ファイルはその後半（定義）である。**

## ★★★★★★残っているもの（明示）

| 段 | 内容 | 状態 |
|---|---|---|
| C2c 前半 | 超平面因子の高さ = 素朴高さ `log max|x_i|` | ★未着手 |
| **E** | **Serre の定理**（ample のある冪が very ample） | ★★★未着手（[Stacks] Lemma 29.40.4） |
| —— | `IsAmple`（加群層、`AmpleDef.lean`）と因子 `D` を繋ぐ橋 | ★未着手 |

★★**本ファイルは `IsAmple` との同値を主張しない**——それが段 E である。
★★★`IsVeryAmpleDiv` を「ample」の**定義に採ってはならない**
（[Stacks] 29.40.4 の (1) ⟺ (3) は定理であって定義ではない）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

/-- ★★★★★★★**very ample な因子**（因子表示、[EGA2] Définition 4.4.2 ／[Stacks] 29.38）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`D` が `ℙᴺ_R` への immersion `i` による**超平面因子の引き戻し**であること。

★★`O(1)` の層を作らずに済むのが因子表示の利点である。 -/
def IsVeryAmpleDiv {X : Scheme.{0}} (R : Type) [CommRing R] (D : X.IdealSheafData) : Prop :=
  ∃ (N : ℕ) (i : X ⟶ projSpace N R), IsImmersion i
    ∧ D = (hyperplaneIdeal N R).comap i

/-- ★★**空虚封じ**——`ℙᴺ_R` の超平面因子そのものは very ample である。

★`i = 𝟙` と取ればよい（恒等射は open immersion、したがって immersion）。 -/
theorem isVeryAmpleDiv_hyperplaneIdeal (N : ℕ) (R : Type) [CommRing R] :
    IsVeryAmpleDiv R (hyperplaneIdeal N R) := by
  refine ⟨N, 𝟙 _, ?_, (Scheme.IdealSheafData.comap_id _).symm⟩
  exact IsImmersion.instOfIsOpenImmersion (𝟙 (projSpace N R))

/-- ★★★**very ample 性は immersion で引き戻せる**。

★原文の証明が `Y → X`（吹き上げ）と `Y → ℙ` を並べる段で要る形である。 -/
theorem IsVeryAmpleDiv.comap {X Y : Scheme.{0}} (R : Type) [CommRing R]
    {D : Y.IdealSheafData} (hD : IsVeryAmpleDiv R D) (g : X ⟶ Y) [IsImmersion g] :
    IsVeryAmpleDiv R (D.comap g) := by
  obtain ⟨N, i, hi, hEq⟩ := hD
  refine ⟨N, g ≫ i, IsImmersion.comp g i, ?_⟩
  rw [hEq, ← Scheme.IdealSheafData.comap_comp]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def IsVeryAmpleDiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(語彙——very ample な因子。Serre の定理は含まない)",
    sectionId := "genell-prop-1-4" }

def isVeryAmpleDiv_hyperplaneIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(空虚封じ——超平面因子は very ample)",
    sectionId := "genell-prop-1-4" }

def IsVeryAmpleDiv.comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(very ample 性は immersion で引き戻せる)",
    sectionId := "genell-prop-1-4" }

def IsVeryAmpleDiv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "hyperplaneIdeal(超平面因子——因子表示での O(-1))"
      (.inProject "ABC3" "ABC3.Found.GenEll.hyperplaneIdeal") 6,
    .citation "[mathlib]" "Scheme.IdealSheafData.comap(イデアル層の引き戻し)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.comap") 7,
    .citation "[mathlib]" "AlgebraicGeometry.IsImmersion"
      (.inMathlib "AlgebraicGeometry.IsImmersion") 7,
    .implicitStep
      ("★★★本ファイルは `IsAmple`(加群層、AmpleDef.lean)との**同値を主張しない**" ++
       "——それが [Stacks] Lemma 29.40.4(Serre の定理)であり、" ++
       "mathlib-gap.json の ample-and-projective-embedding の**段 E**(未着手)である。" ++
       "★消費側 northcott_of_projModel が要求する「正規化座標の単射性」は" ++
       "**閉埋め込みからしか来ない**ので、段 E を迂回して " ++
       "Proposition 1.4, (iv) を取ることはできない(2026-08-28 の測定)") 7 ]

end ABC3.Found.GenEll
