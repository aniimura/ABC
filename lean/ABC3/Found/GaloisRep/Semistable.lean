/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★半安定還元を mathlib の言葉で（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か

`Interface/GenEll/EllModuli.lean` の `SemiStable`（原文の
「`E_L` … with **semi-stable reduction** at all the finite primes of `L`」）は
**posit** されている。★★★2026-08-29 の再実測で、mathlib に

    `Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean`

があり `IsMinimal`・`HasGoodReduction`・`HasMultiplicativeReduction`・
`HasAdditiveReduction` と**三分律**まで揃っていることが分かった
（2026-08-16 の測定は古かった）。

★★★★**本ファイルは半安定還元をその上で定義し、基本的な性質を取る。**

    `HasSemistableReduction R W  ≔  HasGoodReduction R W ∨ HasMultiplicativeReduction R W`

* `hasSemistableReduction_or_hasAdditiveReduction`——三分律の言い換え
* `hasSemistableReduction_iff`——★**`半安定 ⟺ ¬加法的`**
  （`Δ` と `c₄` の付値で排他になる）

## ★★残るのは判定定理

☆「`E[m]` が有理的な `m ≥ 3` があれば半安定」（Néron モデル／
Néron–Ogg–Shafarevich）——★これが `Theorem 3.8` の `torsionExt` に要るものであり、
mathlib には**まだ無い**（2026-08-29 に `#check` で確認）。
★★本ファイルはその**述語の側**を mathlib native にした
——判定定理を書くときの受け皿である。

★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve

variable (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-! ## ★★★★★★★★半安定還元 -/

/-- ★★★★★★★★**半安定還元 ＝ 良還元 または 乗法還元**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `WeierstrassCurve.HasGoodReduction`／`HasMultiplicativeReduction`
（`AlgebraicGeometry/EllipticCurve/Reduction.lean`）の上で書いた。 -/
def HasSemistableReduction (W : WeierstrassCurve K) : Prop :=
  W.HasGoodReduction R ∨ W.HasMultiplicativeReduction R

/-- ★★★★★**三分律の言い換え**——極小方程式なら半安定か加法的である。 -/
theorem hasSemistableReduction_or_hasAdditiveReduction (W : WeierstrassCurve K)
    [WeierstrassCurve.IsMinimal R W] :
    HasSemistableReduction R W ∨ W.HasAdditiveReduction R := by
  rcases hasGoodReduction_or_hasMultiplicativeReduction_or_hasAdditiveReduction R
    (W := W) with h | h | h
  · exact Or.inl (Or.inl h)
  · exact Or.inl (Or.inr h)
  · exact Or.inr h

/-- ★★★★★★**半安定なら加法的でない**——`Δ` と `c₄` の付値で排他になる。 -/
theorem not_hasAdditiveReduction_of_semistable (W : WeierstrassCurve K)
    (h : HasSemistableReduction R W) : ¬ W.HasAdditiveReduction R := by
  intro hadd
  rcases h with hgood | hmult
  · exact absurd hgood.goodReduction (ne_of_lt hadd.badReduction)
  · exact absurd hmult.multiplicativeReduction (ne_of_lt hadd.additiveReduction)

/-- ★★★★★★★★**半安定 ⟺ 加法的でない**（極小方程式のもとで）。

★原文が「semi-stable reduction」と言うのは、まさに「加法的還元を持たない」ことである。 -/
theorem hasSemistableReduction_iff (W : WeierstrassCurve K)
    [WeierstrassCurve.IsMinimal R W] :
    HasSemistableReduction R W ↔ ¬ W.HasAdditiveReduction R := by
  refine ⟨not_hasAdditiveReduction_of_semistable R W, fun h => ?_⟩
  rcases hasSemistableReduction_or_hasAdditiveReduction R W with hs | hadd
  · exact hs
  · exact absurd hadd h

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def HasSemistableReduction.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(半安定還元の定義——mathlib の Reduction.lean の上で)",
    sectionId := "genell-thm-3-8" }

def hasSemistableReduction_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(半安定 ⟺ 加法的でない)",
    sectionId := "genell-thm-3-8" }

def hasSemistableReduction_iff.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      ("WeierstrassCurve.hasGoodReduction_or_hasMultiplicativeReduction_or_hasAdditiveReduction" ++
       "(還元の三分律)")
      (.inMathlib "WeierstrassCurve.hasGoodReduction_or_hasMultiplicativeReduction_or_hasAdditiveReduction") 2,
    .folklore
      ("半安定還元の**判定定理**: 『E[m] が有理的な m ≥ 3 があれば半安定』" ++
       "(Néron モデル／Néron–Ogg–Shafarevich)。" ++
       "★これが Theorem 3.8 の torsionExt に要るものであり、" ++
       "mathlib には**まだ無い**(2026-08-29 に #check で確認)。**残る**") 8,
    .implicitStep
      ("★★★★★測定(2026-08-29): Interface/GenEll/EllModuli.lean の SemiStable は" ++
       "posit されているが、mathlib の AlgebraicGeometry/EllipticCurve/Reduction.lean に" ++
       "IsMinimal・HasGoodReduction・HasMultiplicativeReduction・HasAdditiveReduction と" ++
       "**三分律**が揃っている(2026-08-16 の測定は古かった)。" ++
       "★本ファイルは半安定還元をその上で定義し、" ++
       "**半安定 ⟺ 加法的でない**まで取った——判定定理を書くときの受け皿である") 6 ]

end ABC3.Found.GaloisRep
