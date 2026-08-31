/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.Morphisms.Proper
import Mathlib.AlgebraicGeometry.Morphisms.Immersion
import ABC3.Meta.Claim

/-!
# [Stacks] Lemma 29.42.7 —— **固有な immersion は閉埋め込み**（`Found`）

原典: The Stacks Project, Chapter 29 Morphisms of Schemes, §29.42 Proper morphisms、
物理 p.2520（tag `01W6`）。

原文 (Stacks p.2520):
> In particular, in both cases the image of X in Y is closed.

## ★★★★★★★★なぜこれを作ったのか

`[GenEll] Proposition 1.4, (iv)`（Northcott 性）の証明は、原文 p.7 で

    「[some positive tensor power of] the ample line bundle L_ℚ yields
     an embedding ϵ_ℚ : X_ℚ ⊆ P_ℚ」

と 1 行で畳んでいる。★中身は **2 段**である:

| 段 | 原典 | mathlib |
|---|---|---|
| 1 | `[Stacks] Lemma 29.40.4`: `L_ℚ` が ample なら `ℙⁿ_ℚ` への **immersion** が取れる | ★**無い**——可逆層そのものが無い |
| 2 | `[Stacks] Lemma 29.42.7`: `X_ℚ` が固有なら、その immersion は **閉**埋め込み | ★★**書ける**——本ファイル |

★★★段 2 は `ample` の語彙を一切使わない。**だから今のうちに閉じられる**。

## ★★★★★★段 2 が Northcott にどう効くか

閉埋め込みは mono（`IsClosedImmersion → IsImmersion → IsPreimmersion → Mono`）なので、
**任意の `T` について `T`-点が単射**になる:

    (T ⟶ X_ℚ)  ↪  (T ⟶ ℙⁿ_ℚ)

★`T = Spec ℚ̄` と取れば、これが `Found/GenEll/NorthcottCoord.lean` の
`northcott_of_projModel` が要求する `hinj`（座標比が点を分ける）の**幾何の半分**である。
★★残る半分は「`ℙⁿ(K)` = 座標の組」——`Proj` の点の関手であり、これも mathlib に無い。

## ★★★★mathlib で使った道具（2026-08-27 実測）

| 使ったもの | 場所 |
|---|---|
| `IsProper.of_comp` | `Morphisms/Proper.lean` ——★**本補題 (2) そのもの** |
| `UniversallyClosed.of_comp_of_isSeparated` | 同上 ——★本補題 (1) そのもの |
| `Scheme.Hom.isClosedMap` | `Morphisms/UniversallyClosed.lean` |
| `IsClosedImmersion.of_isPreimmersion` | `Morphisms/ClosedImmersion.lean` |
| `IsPreimmersion → Mono` | `Morphisms/Preimmersion.lean`（インスタンス） |

★★**mathlib は Stacks 29.42.7 を「(1) と (2)」の形で既に持っていた**——
足りなかったのは「immersion ＋ 普遍閉 ⟹ 閉埋め込み」の 1 行だけである。
-/

namespace ABC3.Found.Stacks

open CategoryTheory AlgebraicGeometry

universe u

variable {X Y S : Scheme.{u}}

/-! ## ★★★★★★immersion ＋ 普遍閉 ⟹ 閉埋め込み -/

/-- ★★★★★★**immersion であって普遍閉なら閉埋め込み**。

原文 (Stacks p.2520):
> In particular, in both cases the image of X in Y is closed.

★機構は 1 行——普遍閉なら底空間の写像が閉写像なので、
像 `f '' univ = range f` が閉。あとは `IsClosedImmersion.of_isPreimmersion`。 -/
theorem isClosedImmersion_of_universallyClosed (i : X ⟶ Y)
    [IsImmersion i] [UniversallyClosed i] : IsClosedImmersion i :=
  .of_isPreimmersion i <| by
    simpa using i.isClosedMap _ isClosed_univ

/-! ## ★★★★★★★★[Stacks] Lemma 29.42.7 の形 -/

/-- ★★★★★★★★**[Stacks] Lemma 29.42.7, (1) の帰結** ——
可換三角形 `X ⟶ Y ⟶ S` で `Y` が `S` 上分離的、`X ⟶ S` が普遍閉、
`i : X ⟶ Y` が immersion なら、`i` は**閉埋め込み**。

原文 (Stacks p.2520):
> In particular, in both cases the image of X in Y is closed.

★`UniversallyClosed.of_comp_of_isSeparated` が「`X ⟶ Y` も普遍閉」を出し、
そこへ上の補題を流す。 -/
theorem isClosedImmersion_of_universallyClosed_comp (i : X ⟶ Y) (g : Y ⟶ S)
    [IsImmersion i] [IsSeparated g] [UniversallyClosed (i ≫ g)] :
    IsClosedImmersion i :=
  haveI : UniversallyClosed i := UniversallyClosed.of_comp_of_isSeparated i g
  isClosedImmersion_of_universallyClosed i

/-- ★★★★★★★★**[Stacks] Lemma 29.42.7, (2) の帰結** ——
`X` が `S` 上**固有**で `Y` が `S` 上分離的なら、
`X` から `Y` への immersion は**閉埋め込み**。

原文 (Stacks p.2520):
> In particular, in both cases the image of X in Y is closed.

★★★これが `[GenEll] Proposition 1.4, (iv)` の
「`ϵ_ℚ : X_ℚ ⊆ P_ℚ`」の**閉**性の根拠である——
`X_ℚ` は `ℚ` 上固有、`ℙⁿ_ℚ` は `ℚ` 上分離的。 -/
theorem isClosedImmersion_of_isProper_comp (i : X ⟶ Y) (g : Y ⟶ S)
    [IsImmersion i] [IsSeparated g] [IsProper (i ≫ g)] :
    IsClosedImmersion i :=
  haveI : IsProper i := IsProper.of_comp i g
  isClosedImmersion_of_universallyClosed i

/-! ## ★★★★★★点の単射性 —— Northcott の `hinj` の幾何の半分 -/

/-- ★★**mono なら `T`-点が単射**。 -/
theorem points_injective_of_mono (i : X ⟶ Y) [Mono i] (T : Scheme.{u}) :
    Function.Injective (fun a : T ⟶ X => a ≫ i) :=
  fun _ _ h => (cancel_mono i).1 h

/-- ★★★★★★★★**閉埋め込みなら `T`-点が単射**。

★これが `Found/GenEll/NorthcottCoord.lean` の `northcott_of_projModel` が要求する
`hinj`（座標の比が点を分ける）の**幾何の半分**である。
★★残る半分は「`ℙⁿ(K)` = 座標の組」——`Proj` の点の関手であり、mathlib に無い。 -/
theorem points_injective_of_isClosedImmersion (i : X ⟶ Y) [IsClosedImmersion i]
    (T : Scheme.{u}) :
    Function.Injective (fun a : T ⟶ X => a ≫ i) :=
  points_injective_of_mono i T

/-- ★★★★★★★★**到達点** —— `S` 上固有な `X` から `S` 上分離的な `Y` への immersion は、
すべての `T`-点で単射。

原文 (Stacks p.2520):
> In particular, in both cases the image of X in Y is closed.

★★★`[GenEll] Proposition 1.4, (iv)` で `X = X_ℚ`・`Y = ℙⁿ_ℚ`・`S = Spec ℚ`・
`T = Spec ℚ̄` と取ったものが、原文の
「`ϵ_ℚ : X_ℚ ⊆ P_ℚ` により `X(ℚ̄)` が `P(ℚ̄)` に埋まる」である。 -/
theorem points_injective_of_isProper_comp (i : X ⟶ Y) (g : Y ⟶ S)
    [IsImmersion i] [IsSeparated g] [IsProper (i ≫ g)] (T : Scheme.{u}) :
    Function.Injective (fun a : T ⟶ X => a ≫ i) :=
  haveI := isClosedImmersion_of_isProper_comp i g
  points_injective_of_isClosedImmersion i T

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`)

★★**項目全体の `.src` は Lemma 29.42.7 に対して置く。** 本ファイルは
(1)(2) の両方と「像が閉」の帰結を含んでおり、原文の主張を覆っている。 -/

def isClosedImmersion_of_universallyClosed_comp.src : ABC3.Meta.Source :=
  { paper := "Stacks", pdfPage := 2520,
    item := "Lemma 29.42.7, (1)(X → S が普遍閉 ⟹ X → Y が普遍閉、像は閉)",
    sectionId := "stacks-lemma-29-42-7" }

def isClosedImmersion_of_isProper_comp.src : ABC3.Meta.Source :=
  { paper := "Stacks", pdfPage := 2520,
    item := "Lemma 29.42.7(固有な immersion は閉埋め込み)",
    sectionId := "stacks-lemma-29-42-7" }

/-- ★原文 p.2520 を 150 dpi で目視して数えた（2026-08-27）。 -/
def isClosedImmersion_of_isProper_comp.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsProper.of_comp(本補題 (2) そのもの)"
      (.inMathlib "AlgebraicGeometry.IsProper.of_comp") 2520,
    .citation "[mathlib]" "UniversallyClosed.of_comp_of_isSeparated(本補題 (1) そのもの)"
      (.inMathlib "AlgebraicGeometry.UniversallyClosed.of_comp_of_isSeparated") 2520,
    .citation "[mathlib]" "IsClosedImmersion.of_isPreimmersion(像が閉なら閉埋め込み)"
      (.inMathlib "AlgebraicGeometry.IsClosedImmersion.of_isPreimmersion") 2520,
    .citation "[mathlib]" "Scheme.Hom.isClosedMap(普遍閉なら底空間で閉写像)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.isClosedMap") 2520,
    .implicitStep
      ("★原文は X → X ×_S Y → Y と分解して、第 1 射がグラフゆえ閉埋め込み" ++
       "(Schemes, Lemma 26.21.10)、第 2 射が底変換(Lemma 29.42.5)と論じる。" ++
       "★★形式化では mathlib の of_comp が既にその分解を内蔵しているので、" ++
       "残ったのは「immersion ＋ 普遍閉 ⟹ 閉埋め込み」の 1 行だけだった") 2520,
    .otherPaper "[Stacks]"
      ("Lemma 29.40.4(ample ⟹ ℙⁿ_S への immersion)——★これが段 1 で、" ++
       "本ファイルは段 2 である。★★段 1 は mathlib に無い" ++
       "(可逆層そのものが無い)ので、ResearchPaper/mathlib-gap.json の " ++
       "ample-and-projective-embedding に残っている") 2516,
    .otherPaper "[GenEll]"
      ("Proposition 1.4, (iv)——原文 p.7「[some positive tensor power of] the " ++
       "ample line bundle LQ yields an embedding ϵQ : XQ ⊆ PQ」が本ファイルの消費者。" ++
       "★Northcott 側(northcott_of_projModel)は既に閉じている") 7 ]

def points_injective_of_isProper_comp.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsPreimmersion ⟹ Mono(インスタンス)"
      (.inMathlib "AlgebraicGeometry.IsPreimmersion") 2520,
    .citation "[ABC3]" "northcott_of_projModel の hinj がこれを消費する"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_projModel") 7,
    .implicitStep
      ("★★残る半分は「ℙⁿ(K) = 座標の組」——Proj の点の関手である。" ++
       "mathlib の ProjectiveSpectrum/Functor.lean は Proj の射の関手性しか持たず、" ++
       "点の関手(座標表示)は無い(2026-08-27 実測)") 2516 ]

end ABC3.Found.Stacks
