/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.HtFinJ
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GenEll.VeluImage
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

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
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

/-- **[GenEll] 残る 1 つ**——**非分裂乗法還元の降下**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 918）の測定**——原文は「semi-stable reduction」としか言わないので、
乗法還元は**非分裂**でもよい。しかし Tate 一意化（`tateParamR`）は
**分裂**を要求する。

☆古典的な事実: 非分裂の場合、**不分岐 2 次拡大**で分裂になり、
局所高さ `v(q_E)` も極小判別式の付値 `v_p(Δ_min)` も**変わらない**
（不分岐なので付値群が変わらず、極小モデルも変わらない）。

★したがって `Lemma 3.5` の結論（`deg∞` と `ht^Falt` の不等式）は
分裂の場合から降りる。本節点がその降下である。

☆要るもの: 不分岐 2 次拡大 `Lv′/Lv` とその整数環、
`IsAdicComplete`（第 897 と同じ道）、そして `minDeltaExp` の不変性。 -/
theorem minDeltaExp_descend_of_nonsplit {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ)
    (hss : SemistableAt p E) (hbad : jExp p E < 0)
    (hdescend : True) :
    (l : ℤ) * minDeltaExp p E ≤ minDeltaExp p E' := by
  sorry

def minDeltaExp_descend_of_nonsplit.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(非分裂乗法還元の降下)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_descend_of_nonsplit.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_of_twist(★★★★降下自体は埋まった、第 919)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_twist") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 919）——降下の部分は埋まった**。" ++
       "半安定なら `v_p(Δ_min) = max(0, −v_p(j))` なので " ++
       "`minDeltaExp` は `j` だけで決まり、捧りで変わらない。" ++
       "☆したがって残るのは**`j` が同じ分裂乗法還元の対を 1 つ作る**ことだけである" ++
       "（不分岐な 2 次の捧り。`c₄ ↦ d²c₄`・`c₆ ↦ d³c₆` なので `j` は変わらない）") 4,
    .implicitStep
      ("☆分裂の場合はすでに閉じている——" ++
       "stableLine_is_mu_of_coprime(906) → tateParam_quot_velu_j_dvr(914) " ++
       "→ minDeltaExp_eq_mul_of_veluMu(904) → lemma_3_5_velu_bad_delta(903)") 1,
    .citation "[ABC3]" "hasSplitMultiplicativeReduction_baseChange(分裂の場合、第 910、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.hasSplitMultiplicativeReduction_baseChange") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★`Lemma 3.5` に残る**ただ 1 つ**の節点 -/

/-- ★**原文の「H は乗法還元の各素点で `μ_l` に対応する」の帰結**を型にしたもの。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`lemma_3_5_velu_bad_delta`（第 903、証明済み）が**そのまま消費する形**である。 -/
def IsMuAtBadPrimes {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic] (l : ℕ) : Prop :=
  ∀ p : HeightOneSpectrum (𝓞 L),
    jExp p E < 0 → minDeltaExp p E' = l * minDeltaExp p E

/-- **[GenEll] `Lemma 3.5` に残るただ 1 つの節点**——
`l` が局所高さと互いに素なら、Vélu の商は悪い素点で `Δ_min` を `l` 倍する。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 941）の測定**——**部品はすべて Lean にある**。
本節点はそれらを並べる（長いが機械的な）組み立てである。

☆分裂の場合の連鎖:

    `stableLine_is_mu_of_coprime`(906) → `tateParam_quot_veluCurve_dvr`(927)
      → `minDeltaExp_eq_mul_of_veluMu`(904)

☆非分裂の場合:

    `splits_quadTwist_of_not_isSquare`(938) で捧りを分裂にし、
    `quadTwist_veluQuotientFull`(925) ➕ `veluVFull_quadTwist`(940) で `v`・`w` を送り、
    `j_veluCurve_variableChange`(926) で Tate モデルへ移し、
    上の分裂の連鎖を当ててから `minDeltaExp_eq_mul_of_nonsplit`(929) で降りる。

☆完備化の台: `isAdicComplete_adicCompletionIntegers`(897)、
`isMinimal_baseChange_of_c4`(908)、`hasMultiplicativeReduction_baseChange`(909)、
`dvrTatePhiAddEquiv`(899)。

★これが閉じれば `lemma_3_5_velu_bad_delta`(903) に流せて `Lemma 3.5` が終わる。 -/
theorem isMuAtBadPrimes_of_veluQuotient {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((Finset.range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hssE : ∀ p, SemistableAt p E) (hssE' : ∀ p, SemistableAt p E')
    (hcop : True) :
    IsMuAtBadPrimes E E' l := by
  sorry

def IsMuAtBadPrimes.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(H が μ_l に対応することの帰結を型にしたもの)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残るただ 1 つの節点——部品はすべて揃っている)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_of_veluMu(分裂の連鎖の終点、第 904、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_veluMu") 1,
    .citation "[ABC3]" "stableLine_is_mu_of_coprime(H ↦ μ_l、第 906、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.stableLine_is_mu_of_coprime") 1,
    .citation "[ABC3]" "splits_quadTwist_of_not_isSquare(非分裂を捧りで分裂に、第 938、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.splits_quadTwist_of_not_isSquare") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_nonsplit(捧った対から降りる、第 929、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_nonsplit") 1,
    .citation "[ABC3]" "isAdicComplete_adicCompletionIntegers(完備化の台、第 897、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isAdicComplete_adicCompletionIntegers") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 941-942）の測定**——" ++
       "局所の数学の部品はすべて揃っているが、" ++
       "組み立てには**もう 1 層**が要る。" ++
       "☆`lemma_3_2_i_tate_all`（第 906 が使う）は `K ⊆ L`（`IsGalois K L`）を受けている" ++
       "——`μ_l` や `l`-捉れ点は完備化 `Lv` の中にとは限らないからである。" ++
       "★したがって `Lv` の有限 Galois 拡大 `L′` とその整数環を立て、" ++
       "そこで `IsAdicComplete`（第 897 と同じ 3 段）を取り、" ++
       "`Δ_min` が拡大でどう変わるか（分岐指数）を扱う層が別に要る") 8,
    .implicitStep
      ("☆hcop : True は「l が局所高さと互いに素」（原文の仮定）を型に書く段") 3,
    .implicitStep
      ("☆逆差の記録: 本プロジェクトの `lemma_3_5_velu` 系は " ++
       "`Q : E.toAffine.Point`（`L`-有理な点）で `H` を与えているが、" ++
       "原文の `H_F` は部分群スキームであり点が有理とは限らない。" ++
       "★Vélu の商の式は対称式なので係数は `L` に落ちるが、" ++
       "その降下を型で書く段が別に要る") 5 ]

end ABC3.Skeleton.GenEll
