/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PicEquivRing
import Mathlib.RingTheory.PicardGroup
import Mathlib.NumberTheory.NumberField.ClassNumber
import ABC3.Meta.Claim

/-!
# `Pic(Spec R) ≃* Cl(R)` —— **類群との橋**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★★★★★台帳の D2（`APic(Spec 𝓞_F) ≅ ADiv(F)/APrc(F)`）の**有限素点側**

`Definition 1.2` の `ht_M̄(x) ≝ deg_F(x_F^* M̄)` を線束表示で書くには、
`Spec 𝓞_F` 上の算術直線束の群を `ADiv(F)/APrc(F)` と繋ぐ必要がある。
★その**有限素点側**は 2 つの同型を継ぐだけで出る:

    `PicSheaf (Spec R) ≃* CommRing.Pic R`   —— `Found/Arakelov/PicEquivRing.lean`
    `ClassGroup R ≃* CommRing.Pic R`        —— ★mathlib の `ClassGroup.equivPic`

★★`ClassGroup.equivPic` は mathlib に**ある**（2026-08-28 実測、`RingTheory/PicardGroup.lean`）。

## ★★★★★★★★★測定 —— `deg_F` が捻れ集合表示で書けない理由（設計の問題）

`APicOf X = (Pic X × Multiplicative (共役不変な連続関数)) / Γ(X,𝓞^×)` は
**原文どおりの同型類の群である**（`APicQuot.lean`）。
★`deg_F` をそこに載せようとすると

    `deg(L, g) = deg(L₀, base_{[L]}) + (アルキメデス側の g の寄与)`

となるが、`base_{[L]}` は `TorsorMetric.base` の **`Classical.choice`** である。

★★★**すると `deg` の加法性が落ちる**——`base_{[L·M]}` と
`base_{[L]} ⊗ base_{[M]}` は一致しないからである。

★★★★したがって **`Classical.choice` の基準計量では `deg_F` は作れない**。
道は「基準計量を**整合的に**選ぶ」ことであり、それは
**分数イデアル `𝔞 ⊂ F` を標準の計量つきで代表に取る**こと、すなわち
`ADiv(F)` の側から作ることに他ならない。★本ファイルはその第 1 歩である。

## ★残っている段（明示）

1. ★`ADiv(F)` の有限素点部分 ≅ 分数イデアル（Dedekind の分解）
2. ★★`ADiv(F)/APrc(F) ≃* APicOf (Spec 𝓞_F)`——代表を `𝔞` に取る
3. ★★★`deg_F` の転送と、`degNormalizedAPic`（在庫）との一致
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite NumberField

/-- ★★★★★★**`Pic(Spec R) ≃* Cl(R)`**（`R` は整域）。

★2 つの同型を継ぐだけである——
`equivPicRingSheaf`（本プロジェクト）と `ClassGroup.equivPic`（mathlib）。 -/
noncomputable def picSheafEquivClassGroup (R : CommRingCat.{0}) [IsDomain R] :
    PicSheaf (Spec R) ≃* ClassGroup (R : Type) :=
  (equivPicRingSheaf R).trans (ClassGroup.equivPic (R : Type)).symm

/-- ★★★★★★★**`Pic(Spec 𝓞_F) ≃* Cl(F)`**——数体の場合。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x. -/
noncomputable def picSheafEquivClassGroupOF (F : Type) [Field F] [NumberField F] :
    PicSheaf (Spec (CommRingCat.of (𝓞 F))) ≃* ClassGroup (𝓞 F) :=
  (equivPicRingSheaf (CommRingCat.of (𝓞 F))).trans (ClassGroup.equivPic (𝓞 F)).symm

/-! ### ★出典の紐付け(`.src`) -/

def picSheafEquivClassGroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(層 B——Pic(Spec R) ≃* Cl(R))",
    sectionId := "genell-def-1-2-i" }

def picSheafEquivClassGroupOF.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(層 D——Pic(Spec 𝓞_F) ≃* Cl(F)。deg_F の橋の有限素点側)",
    sectionId := "genell-def-1-2-i" }

def picSheafEquivClassGroup.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "equivPicRingSheaf(Pic(Spec R) ≃* CommRing.Pic R)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.equivPicRingSheaf") 4,
    .citation "[mathlib]" "ClassGroup.equivPic(類群 ≃* Picard 群)"
      (.inMathlib "ClassGroup.equivPic") 4,
    .implicitStep
      ("★★★測定(2026-08-28): deg_F は捻れ集合表示 (L, g) では**作れない**。" ++
       "deg(L, g) = deg(L₀, base_{[L]}) + (g の寄与) となるが " ++
       "base は TorsorMetric.base の Classical.choice なので、" ++
       "base_{[L·M]} ≠ base_{[L]} ⊗ base_{[M]} となり**加法性が落ちる**") 4,
    .implicitStep
      ("★★★★したがって道は「基準計量を整合的に選ぶ」ことであり、" ++
       "それは分数イデアル 𝔞 ⊂ F を標準の計量つきで代表に取ること、" ++
       "すなわち ADiv(F) の側から作ることに他ならない。" ++
       "★残る段: (1) ADiv(F) の有限部分 ≅ 分数イデアル、" ++
       "(2) ADiv(F)/APrc(F) ≃* APicOf (Spec 𝓞_F)、(3) deg_F の転送") 4 ]

end ABC3.Found.Arakelov
