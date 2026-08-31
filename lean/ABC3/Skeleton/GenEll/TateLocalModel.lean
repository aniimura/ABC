/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.HtFinJ
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace

/-!
# `Lemma 3.5` に残る 2 つ —— **局所モデルの完備化への移行**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★これは何か

2026-08-31（第 899-906）の時点で、`Lemma 3.5` の連鎖は

    `stableLine_is_mu_of_coprime`（第 906）
      → `tateParam_quot_velu_dvr`（第 900）
      → `minDeltaExp_eq_mul_of_veluMu`（第 904）
      → `lemma_3_5_velu_bad_delta`（第 903）

と全部繋がっている。★残っているのは**局所モデルを完備化に移す 2 段**だけである。
本ファイルはその 2 つを**依存グラフの節点として型で固定する**。

| 節点 | 内容 |
|---|---|
| `isMinimal_baseChange` | 極小モデルは完備化でも極小 |
| `hasSplitMultiplicativeReduction_baseChange` | 半安定＋悪いなら分裂乗法還元 |
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

/-- **[GenEll] 残る 2 つの (1)**——**極小モデルは完備化でも極小**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆極小判別式の付値は完備化で変わらない（Tate のアルゴリズムは剰余体と付値しか見ない）。
★`DegInfLocal.lean` の定理群はこれを**インスタンス仮説として受けている**。 -/
theorem isMinimal_baseChange {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (hmin : IsMinimal (primeSubring p) W) :
    (W.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L) := by
  sorry

def isMinimal_baseChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小モデルは完備化でも極小)",
    sectionId := "genell-lemma-3-5" }

def isMinimal_baseChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMinimal_baseChange_of_c4(★★★★Lemma 3.5 が要る場合は埋まった、第 908)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isMinimal_baseChange_of_c4") 1,
    .implicitStep
      ("★★★★**後日談（2026-08-31、第 908）——本節点は `Lemma 3.5` には要らない**。" ++
       "`Lemma 3.5` が使うのは**乗法還元の場合だけ**で、" ++
       "そこでは `v_p(c₄) = 0` なので `isMinimal_of_c4_vAdd_eq_zero` が**直接**極小性を与える" ++
       "（`Found/GaloisRep/DegInfTateParam.lean` の `isMinimal_baseChange_of_c4`）。" ++
       "☆本節点（一般の極小性の移行）は弱近似を要するが、臨界の道には乗っていない") 1,
    .citation "[ABC3]" "isIntegral_baseChange_of_isIntegral(整性の移行、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isIntegral_baseChange_of_isIntegral") 1,
    .implicitStep
      ("★極小性の移行: `R` 上の整モデルの判別式の付値の下限が " ++
       "`primeSubring p` 上のそれと一致すること。" ++
       "☆剰余体も付値群も完備化で変わらないので、Tate のアルゴリズムの出力は同じである。" ++
       "★★mathlib には Tate のアルゴリズムは無いので、" ++
       "`v(Δ_min)` の定義（整モデルの下限）から直接見るのが早い——" ++
       "`≤` は底変換ですぐ出る。`≥` は `R` の整モデルを " ++
       "`primeSubring p` に近似する段（弱近似）を要する") 6 ]

/-- ★★★★★★★★★★★★**[GenEll] 残る 1 つ**——**分裂乗法還元は
「剰余体で 2 次式が分裂する」だけで得られる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-08-31（第 910）**——第 907 では `hsplit : True` という**空欄**だった。
第 909 で乗法還元の側（`v(Δ) < 1`・`v(c₄) = 1`）が埋まったので、
★**空欄を mathlib の構造体のフィールそのもの（`Splits`）に置き換え、定理を証明した**。

☆逆差の記録（第 907 から継続）——原文は「semi-stable reduction」としか言わないので
乗法還元は**非分裂**でもよい。非分裂の場合は不分岐 2 次拡大で分裂になり、
局所高さ `v(q_E)` は変わらない——その降下は別節点である。 -/
theorem hasSplitMultiplicativeReduction_baseChange {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0)
    (hΔpos : 0 < valAdd p (Units.mk0 W.Δ hΔ))
    [(W.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (hsplits : Polynomial.Splits (Polynomial.map
      (algebraMap (p.adicCompletionIntegers L)
        (IsLocalRing.ResidueField (p.adicCompletionIntegers L)))
      (Polynomial.C ((W.baseChange (p.adicCompletion L)).integralModel
            (p.adicCompletionIntegers L)).c₄ * Polynomial.X ^ 2
        + Polynomial.C (((W.baseChange (p.adicCompletion L)).integralModel
              (p.adicCompletionIntegers L)).a₁
            * ((W.baseChange (p.adicCompletion L)).integralModel
              (p.adicCompletionIntegers L)).c₄) * Polynomial.X
        - Polynomial.C (54 * ((W.baseChange (p.adicCompletion L)).integralModel
              (p.adicCompletionIntegers L)).b₆
            - 3 * ((W.baseChange (p.adicCompletion L)).integralModel
              (p.adicCompletionIntegers L)).b₂
              * ((W.baseChange (p.adicCompletion L)).integralModel
              (p.adicCompletionIntegers L)).b₄
            + ((W.baseChange (p.adicCompletion L)).integralModel
              (p.adicCompletionIntegers L)).a₂
              * ((W.baseChange (p.adicCompletion L)).integralModel
              (p.adicCompletionIntegers L)).c₄)))) :
    (W.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L) := by
  haveI := hasMultiplicativeReduction_baseChange
    (R := p.adicCompletionIntegers L) p hp W hΔ hc4ne hc4 hΔpos
  exact { splitMultiplicativeReduction := hsplits }

def hasSplitMultiplicativeReduction_baseChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半安定で悪ければ分裂乗法還元)",
    sectionId := "genell-lemma-3-5" }

def hasSplitMultiplicativeReduction_baseChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "hasMultiplicativeReduction_baseChange(乗法還元の側、第 909、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hasMultiplicativeReduction_baseChange") 1,
    .implicitStep
      ("★★★★**2026-08-31（第 910）——本節点は証明済みになった**。" ++
       "第 907 の `hsplit : True` を mathlib の構造体のフィール（`Splits`）そのものに置き換えた。" ++
       "☆したがって `Lemma 3.5` に残るのは" ++
       "「非分裂の場合を不分岐 2 次拡大へ上げて降りる」段だけである" ++
       "（局所高さ v(q_E) は変わらない）") 5 ]

end ABC3.Skeleton.GenEll
