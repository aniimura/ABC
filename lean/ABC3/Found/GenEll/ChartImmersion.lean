/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjSpaceCover
import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.Morphisms.ClosedImmersion
import Mathlib.AlgebraicGeometry.Morphisms.Immersion

/-!
# ★★★★★★★段 E3d —— チャートごとの閉埋め込みは大域の閉埋め込みになる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★これは何か —— 段 E3d

段 E3c は「各チャートで `A⁰_{x_i} → Γ(X, X_{s_i})` が全射」を言う。
★そこから「`ψ : X ⟶ ℙᴺ_R` が**埋め込み**である」へ渡すのが段 E3d である。

★★橋は 2 本しかない:

| 段 | 道具 |
|---|---|
| 全射 ⟹ チャートで閉埋め込み | `IsClosedImmersion.of_surjective_of_isAffine`（mathlib） |
| チャートで閉埋め込み ⟹ 大域で閉埋め込み | `IsZariskiLocalAtTarget @IsClosedImmersion`（mathlib） |

★★★あとは「`D₊(x_i)` が `ℙᴺ` を覆う」ことだけであり、それは
`§9-Cover` の `exists_X_notMem`（無関係イデアルは変数で生成される）から出る。

## ★測定の記録

★**閉埋め込みが得られる**（immersion より強い）——`D₊(x_i)` は `ℙᴺ` **全体**を覆うからである。
★★原典が「embedding」としか書かないのは、`X_ℚ` の側で使うときに
生成ファイバーへ落とすと閉性が落ちうるためであろう（段 F2c の話）。
本ファイルは**強い方（閉埋め込み）を出し、そこから `IsImmersion` を導く**。

## ★残っている段（明示）

★★★段 E3c の残り（「`a/s^n` が `A⁰_{x_i}` の像に入る」）と、
`§9-816` の `chartToProj` を貼って `ψ` を作る段（段 E1d）は本ファイルに無い。
★本ファイルは**その 2 つが揃ったときに何が言えるか**を型で固定したものである。
-/

namespace ABC3.Found.GenEll

open MvPolynomial ProjectiveSpectrum AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★★★`D₊(x_i)` は `ℙᴺ` を覆う（`Proj.Opens` の側） -/

/-- ★★★★**標準チャート `D₊(x_i)` は `ℙᴺ` を覆う**（`(Proj 𝒜).Opens` の側）。

★`§9-Cover` の `iSup_basicOpen_X_eq_top` は `ProjectiveSpectrum` の開の側だった。
本補題は `Proj.basicOpen`（スキームの開）で言い直したもので、
★★`IsZariskiLocalAtTarget` の被覆条件がこの形を要求する。 -/
theorem iSup_projBasicOpen_X_eq_top (N : ℕ) (R : Type) [CommRing R] :
    (⨆ i : Fin (N + 1),
      Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i)) = ⊤ := by
  refine eq_top_iff.2 (fun z _ => TopologicalSpace.Opens.mem_iSup.2 ?_)
  obtain ⟨i, hi⟩ := exists_X_notMem (Fin (N + 1)) R z
  exact ⟨i, hi⟩

/-! ## ★★★★★★★段 E3d -/

/-- ★★★★★★**閉埋め込みは対象で局所的である**（mathlib の言い換え）。

★被覆の各員の上で閉埋め込みなら、全体で閉埋め込みである。 -/
theorem isClosedImmersion_of_charts {X Y : Scheme.{0}} (ψ : X ⟶ Y) {ι : Type}
    (V : ι → Y.Opens) (hcov : (⨆ i, V i) = ⊤)
    (h : ∀ i, IsClosedImmersion (ψ ∣_ V i)) : IsClosedImmersion ψ :=
  (IsZariskiLocalAtTarget.iff_of_iSup_eq_top V hcov).2 h

/-- ★★★★★★★**段 E3d** —— 各標準チャートで閉埋め込みなら `ψ : X ⟶ ℙᴺ_R` は閉埋め込みである。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★機構は「`D₊(x_i)` が `ℙᴺ` を覆う」（本ファイル）＋
「閉埋め込みは対象で局所的」（mathlib）だけである。 -/
theorem isClosedImmersion_of_projCharts {X : Scheme.{0}} (N : ℕ) (R : Type) [CommRing R]
    (ψ : X ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))
    (h : ∀ i : Fin (N + 1), IsClosedImmersion
      (ψ ∣_ Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i))) :
    IsClosedImmersion ψ :=
  isClosedImmersion_of_charts ψ _ (iSup_projBasicOpen_X_eq_top N R) h

/-- ★★★★★★★**したがって埋め込みである** —— 原典の「yields an embedding」の形。

★閉埋め込みは埋め込みより強い。`D₊(x_i)` が `ℙᴺ` **全体**を覆うので強い方が出る。 -/
theorem isImmersion_of_projCharts {X : Scheme.{0}} (N : ℕ) (R : Type) [CommRing R]
    (ψ : X ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))
    (h : ∀ i : Fin (N + 1), IsClosedImmersion
      (ψ ∣_ Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i))) :
    IsImmersion ψ :=
  haveI := isClosedImmersion_of_projCharts N R ψ h
  inferInstance

/-- ★★★★★**段 E3c から段 E3d への橋** —— アフィンなチャートで `Γ` が全射なら閉埋め込み。

★これが段 E3c（`A⁰_{x_i} → Γ(X, X_{s_i})` の全射性）を消費する口である。
★★機構は mathlib の `IsClosedImmersion.of_surjective_of_isAffine` そのものである。 -/
theorem isClosedImmersion_restrict_of_surjective {X Y : Scheme.{0}} (ψ : X ⟶ Y) (V : Y.Opens)
    [IsAffine (ψ ⁻¹ᵁ V : X.Opens)] [IsAffine (V : Y.Opens)]
    (h : Function.Surjective (ψ ∣_ V).appTop.hom) :
    IsClosedImmersion (ψ ∣_ V) :=
  IsClosedImmersion.of_surjective_of_isAffine _ h

/-! ## ★出典の紐付け(`.src`) -/

def iSup_projBasicOpen_X_eq_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(標準チャートは ℙᴺ を覆う——Proj.Opens の側)",
    sectionId := "genell-prop-1-4" }

def isClosedImmersion_of_projCharts.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3d——チャートごとの閉埋め込みは大域の閉埋め込み)",
    sectionId := "genell-prop-1-4" }

def isImmersion_of_projCharts.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(したがって埋め込みである)",
    sectionId := "genell-prop-1-4" }

def isClosedImmersion_restrict_of_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c から段 E3d への橋——Γ が全射なら閉埋め込み)",
    sectionId := "genell-prop-1-4" }

def isClosedImmersion_of_projCharts.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsZariskiLocalAtTarget @IsClosedImmersion(閉埋め込みは対象で局所的)"
      (.inMathlib "AlgebraicGeometry.IsClosedImmersion.isZariskiLocalAtTarget") 2,
    .citation "[mathlib]" "IsClosedImmersion.of_surjective_of_isAffine(アフィン＋全射)"
      (.inMathlib "AlgebraicGeometry.IsClosedImmersion.of_surjective_of_isAffine") 2,
    .citation "[ABC3]" "exists_X_notMem(無関係イデアルは変数で生成される、段 C2a)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_X_notMem") 2,
    .implicitStep
      ("★**閉埋め込みが得られる**(immersion より強い)——D₊(x_i) が ℙᴺ **全体**を覆うからである。" ++
       "★★原典が「embedding」としか書かないのは、X_ℚ の側で使うときに" ++
       "生成ファイバーへ落とすと閉性が落ちうるためであろう(段 F2c の話)") 6,
    .implicitStep
      ("★★★段 E3c の残り(「a/s^n が A⁰_{x_i} の像に入る」)と、" ++
       "§9-816 の chartToProj を貼って ψ を作る段(段 E1d)は本ファイルに無い。" ++
       "★本ファイルは**その 2 つが揃ったときに何が言えるか**を型で固定したものである") 7 ]

end ABC3.Found.GenEll
