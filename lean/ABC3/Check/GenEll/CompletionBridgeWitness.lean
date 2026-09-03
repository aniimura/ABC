/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GaloisRep.AdicCompleteIntegers

/-!
# 界面の測定 —— **★★★★★★★★★★完備化の橋は非空虚である**（`Check`）

**これは原典の主張ではない**（我々のモデルについての事実）ので `.src` を持たない。

## ★★★★★★★★★★ 2026-08-31 の測定（第 898）

`Found/GaloisRep/DegInfLocal.lean`・`DegInfTateParam.lean` の局所の定理群は
抽象な `R`・`Lv` の上で述べられている。
★**それらが数体の素点の完備化で実際に合成できること**をここで確かめる。

    `R := p.adicCompletionIntegers L`   `Lv := p.adicCompletion L`

☆第 896 の測定では `IsAdicComplete` が無くてこのファイルは書けなかった。
★第 897 でそれを作ったので、今は**通る**。

## ☆これが意味すること

`Lemma 3.5` の残る道は「大域の `E/L` を各悪い素点で完備化に落とす」ことであるが、
**その落とし先の型はもうある**。
☆残るのは内容（分裂乗法還元、Tate 一意化、`H` の像）だけである。

## ★配管の注意（実測）

インスタンスの存在を `have i : C X := inferInstance` を**並べて**確かめてはいけない。
★`have i1 : Field Lv := inferInstance` を置いた瞬間、局所文脈に **2 つ目の
`Field Lv`** が入るので、続く `Algebra L Lv` の合成が菱形で落ちる。
☆各インスタンスは**別々の `example`** で確かめること。
-/

namespace ABC3.Check.GenEll

open IsDedekindDomain NumberField IsLocalRing WeierstrassCurve ABC3.Found.GaloisRep

/-! ## ★局所の台は揃っている -/

noncomputable example (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) :
    Algebra L (p.adicCompletion L) := inferInstance

example (L : Type) [Field L] [NumberField L] (p : HeightOneSpectrum (𝓞 L)) :
    IsDomain (p.adicCompletionIntegers L) := inferInstance

example (L : Type) [Field L] [NumberField L] (p : HeightOneSpectrum (𝓞 L)) :
    IsDiscreteValuationRing (p.adicCompletionIntegers L) := inferInstance

noncomputable example (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) :
    Algebra (p.adicCompletionIntegers L) (p.adicCompletion L) := inferInstance

example (L : Type) [Field L] [NumberField L] (p : HeightOneSpectrum (𝓞 L)) :
    IsFractionRing (p.adicCompletionIntegers L) (p.adicCompletion L) := inferInstance

/-- ★★★★★★★★★★第 897 で作ったもの。 -/
example (L : Type) [Field L] [NumberField L] (p : HeightOneSpectrum (𝓞 L)) :
    IsAdicComplete (maximalIdeal (p.adicCompletionIntegers L))
      (p.adicCompletionIntegers L) := inferInstance

/-! ## ★★★★★★★★★★第 892 の定理が実際に当てられる -/

/-- ★★★★★★★★★★**`minDeltaExp_eq_mul_of_tateParamR`（第 892）は
数体の素点の完備化で実際に使える**。

☆型が通ること自体が測定である——仮説を全部受けたうえで結論を出す。 -/
theorem minDeltaExp_eq_mul_at_completion
    (L : Type) [Field L] [NumberField L] (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) (l : ℕ)
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (h' : (E'.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    (hq : tateParamR (E'.baseChange (p.adicCompletion L)) h'
        = (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) :
    minDeltaExp p E' = l * minDeltaExp p E :=
  minDeltaExp_eq_mul_of_tateParamR (R := p.adicCompletionIntegers L) E E' l h h' p hp
    C hC hc4ne hc4 C' hC' hc4ne' hc4' hq

end ABC3.Check.GenEll
