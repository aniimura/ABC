/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.ResidueFieldFinite
import ABC3.Found.GenEll.SplitDichotomy

/-!
# 第 992 ブロック —— **★★★★★★★★★★★★★★★★★★★★完備化での分裂の二者択一**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か

第 982（分裂性の二者択一）は `[Fintype k]` を要求していたので、
**数体の素点の完備化には当てられなかった**（第 983 で実測）。
☆第 989 でその有限性を内製したので、**いま当てられる**。

★仮説は 3 つだけである:

* `W.a₁ = 0`・`W.a₃ = 0`（体の側の正規化。第 981 が整モデルへ運ぶ）
* 完備化で乗法還元をもつこと（第 976 が与える。第 983 が `c₄ ≠ 0` に直す）

☆結論は「整モデルの 2 次式が分裂する」か「ある捻りで分裂する」かの二者択一。
★捻り `d` は完備化の整数環の元だが、第 990 で `𝓞 L` の元に取り替えられる。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll IsDedekindDomain NumberField WeierstrassCurve Polynomial
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★**完備化での分裂の二者択一**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 992）**——第 989（剰余体の有限性）が入ったので、
第 982 の二者択一が数体の素点の完備化にそのまま当たる。 -/
theorem splits_or_twist_at_completion (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (ha1 : W.a₁ = 0) (ha3 : W.a₃ = 0)
    [WeierstrassCurve.HasMultiplicativeReduction (p.adicCompletionIntegers L)
      (W.baseChange (p.adicCompletion L))] :
    Splits (Polynomial.map
        (algebraMap _ (IsLocalRing.ResidueField (p.adicCompletionIntegers L)))
        (C (integralModel (p.adicCompletionIntegers L)
              (W.baseChange (p.adicCompletion L))).c₄ * X ^ 2
          + C ((integralModel (p.adicCompletionIntegers L)
                (W.baseChange (p.adicCompletion L))).a₁
              * (integralModel (p.adicCompletionIntegers L)
                (W.baseChange (p.adicCompletion L))).c₄) * X
          - C (54 * (integralModel (p.adicCompletionIntegers L)
                  (W.baseChange (p.adicCompletion L))).b₆
              - 3 * (integralModel (p.adicCompletionIntegers L)
                    (W.baseChange (p.adicCompletion L))).b₂
                  * (integralModel (p.adicCompletionIntegers L)
                    (W.baseChange (p.adicCompletion L))).b₄
              + (integralModel (p.adicCompletionIntegers L)
                    (W.baseChange (p.adicCompletion L))).a₂
                  * (integralModel (p.adicCompletionIntegers L)
                    (W.baseChange (p.adicCompletion L))).c₄)))
      ∨ ∃ d : (p.adicCompletionIntegers L),
        (algebraMap _ (IsLocalRing.ResidueField (p.adicCompletionIntegers L))) d ≠ 0
        ∧ Splits (Polynomial.map
          (algebraMap _ (IsLocalRing.ResidueField (p.adicCompletionIntegers L)))
          (C (quadTwist (integralModel (p.adicCompletionIntegers L)
                (W.baseChange (p.adicCompletion L))) d).c₄ * X ^ 2
            + C ((quadTwist (integralModel (p.adicCompletionIntegers L)
                  (W.baseChange (p.adicCompletion L))) d).a₁
                * (quadTwist (integralModel (p.adicCompletionIntegers L)
                  (W.baseChange (p.adicCompletion L))) d).c₄) * X
            - C (54 * (quadTwist (integralModel (p.adicCompletionIntegers L)
                    (W.baseChange (p.adicCompletion L))) d).b₆
                - 3 * (quadTwist (integralModel (p.adicCompletionIntegers L)
                      (W.baseChange (p.adicCompletion L))) d).b₂
                    * (quadTwist (integralModel (p.adicCompletionIntegers L)
                      (W.baseChange (p.adicCompletion L))) d).b₄
                + (quadTwist (integralModel (p.adicCompletionIntegers L)
                      (W.baseChange (p.adicCompletion L))) d).a₂
                    * (quadTwist (integralModel (p.adicCompletionIntegers L)
                      (W.baseChange (p.adicCompletion L))) d).c₄))) := by
  haveI : Fintype (IsLocalRing.ResidueField (p.adicCompletionIntegers L)) :=
    Fintype.ofFinite _
  have ha1' : (W.baseChange (p.adicCompletion L)).a₁ = 0 := by
    show (W.map (algebraMap L (p.adicCompletion L))).a₁ = 0
    rw [WeierstrassCurve.map_a₁, ha1, map_zero]
  have ha3' : (W.baseChange (p.adicCompletion L)).a₃ = 0 := by
    show (W.map (algebraMap L (p.adicCompletion L))).a₃ = 0
    rw [WeierstrassCurve.map_a₃, ha3, map_zero]
  haveI := isCharNeTwoNF_integralModel (R := p.adicCompletionIntegers L)
    (W.baseChange (p.adicCompletion L)) ha1' ha3'
  exact splits_or_exists_twist_splits''
    (algebraMap (p.adicCompletionIntegers L)
      (IsLocalRing.ResidueField (p.adicCompletionIntegers L)))
    Ideal.Quotient.mk_surjective
    (integralModel (p.adicCompletionIntegers L) (W.baseChange (p.adicCompletion L)))
    (residue_c4_ne_zero_of_multiplicativeReduction _)

/-! ## ★出典の紐付け(`.src`) -/

def splits_or_twist_at_completion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(完備化での分裂の二者択一。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
