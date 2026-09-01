/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.HtFinJ
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GenEll.VeluImage
import ABC3.Found.GaloisRep.SplitAtCompletion
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import ABC3.Skeleton.GenEll.TateIsogeny

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
    (hcop : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) :
    IsMuAtBadPrimes E E' l := by
  sorry

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★第 1005 —— `IsMuAtBadPrimes` まで組む

★第 1004 で悪い素点 1 個あたりの局所入力は `hmin`・`h`・`hlu` の 3 本になった。
☆`hmin` は第 954（半安定性から `C`・`hC`・`hc4ne`・`hc4`）＋第 973 で**自前で作れる**。

★★したがって `IsMuAtBadPrimes` は **残り 2 本**を仮説に置けば閉じる:

* `hlu`——`l` が `p` の完備化の整数環で単元（すなわち `p ∤ l`）
* `hsplit`——極小化した `C • E` が完備化で**分裂**乗法還元をもつこと

☆後者は第 976＋993 が「分裂または捻りで分裂」まで詰めており、
**`p ∣ 2` で非分裂の場合だけ**が不分岐 2 次拡大（943＋944）待ちである。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] `IsMuAtBadPrimes`——
残る仮説は `p ∤ l` と分裂性の 2 本**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1005）**——第 944 から始まった連鎖の到達点である。
☆`hsplit` と `hlu` を仮説に置けば `Lemma 3.5` の残る節点はこれで閉じる。 -/
theorem isMuAtBadPrimes_of_veluQuotient_of_split {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hssE : ∀ p, SemistableAt p E) (hssE' : ∀ p, SemistableAt p E')
    (hcop : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E))
    (hlu : ∀ p : HeightOneSpectrum (𝓞 L), IsUnit ((l : (p.adicCompletionIntegers L))))
    (hsplit : ∀ (p : HeightOneSpectrum (𝓞 L)), jExp p E < 0 →
      ∀ (C : WeierstrassCurve.VariableChange L)
        (_hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
          ((C • E).baseChange (p.adicCompletion L))),
        ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
          (p.adicCompletionIntegers L)) :
    IsMuAtBadPrimes E E' l := by
  intro p hbad
  obtain ⟨C, hC, hc4ne, hc4⟩ :=
    ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg p E (hssE p) hbad
  haveI hCE : ((C • E).baseChange (p.adicCompletion L)).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, isUnit_iff_ne_zero]
    show ((C • E).map (algebraMap L (p.adicCompletion L))).Δ ≠ 0
    rw [WeierstrassCurve.map_Δ]
    exact (map_ne_zero_iff _ (algebraMap L (p.adicCompletion L)).injective).2
      (variableChange_Delta_ne_zero E E.isUnit_Δ.ne_zero C)
  haveI hCE' : ((C • E').baseChange (p.adicCompletion L)).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, isUnit_iff_ne_zero]
    show ((C • E').map (algebraMap L (p.adicCompletion L))).Δ ≠ 0
    rw [WeierstrassCurve.map_Δ]
    exact (map_ne_zero_iff _ (algebraMap L (p.adicCompletion L)).injective).2
      (variableChange_Delta_ne_zero E' E'.isUnit_Δ.ne_zero C)
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  have hmin := ABC3.Found.GaloisRep.isMinimal_baseChange_at_bad_prime p hp E C hC hc4ne hc4
  exact minDeltaExp_eq_mul_at_bad_prime_full p E E' (hssE p) (hssE' p) hbad hl hodd
    (hcop p hbad) C hmin (hsplit p hbad C hmin) (hlu p) hQ hE'

def isMuAtBadPrimes_of_veluQuotient_of_split.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsMuAtBadPrimes——残る仮説は p ∤ l と分裂性の 2 本)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient_of_split.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_full(第 1004、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_full") 1,
    .citation "[ABC3]" "exists_minimal_c4_unit_of_jExp_neg(C・hC・hc4、第 954、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg") 1,
    .citation "[ABC3]" "isMinimal_baseChange_at_bad_prime(hmin、第 973、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isMinimal_baseChange_at_bad_prime") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1005）の測定**——残るのは 2 本だけである。" ++
       "☆`hlu` は `p ∤ l`。★`hsplit` は第 976＋993 が" ++
       "「分裂または捻りで分裂」まで詰めており、" ++
       "`p ∣ 2` で非分裂の場合だけが不分岐 2 次拡大（943＋944）待ちである") 2 ]

def IsMuAtBadPrimes.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(H が μ_l に対応することの帰結を型にしたもの)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残るただ 1 つの節点——部品はすべて揃っている)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_of_globalVelu(大域の Vélu の商で受ける形、第 972、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_globalVelu") 1,
    .citation "[ABC3]" "exists_minimal_c4_unit_of_jExp_neg(C・hC・hc4、第 954、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg") 1,
    .citation "[ABC3]" "hasMultiplicativeReduction_at_bad_prime(極小性と乗法還元、第 976、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hasMultiplicativeReduction_at_bad_prime") 1,
    .citation "[ABC3]" "valuation_algebraMap_adicCompletion(hp、第 964、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion") 1,
    .citation "[ABC3]" "isUnit_natCast_of_valAdd_eq_zero(hlu、第 953、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_natCast_of_valAdd_eq_zero") 1,
    .citation "[ABC3]" "not_dvd_vAdd_tateParam_of_not_dvd_jExp(hcop、第 978、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp") 1,
    .citation "[ABC3]" "tateModel_map_Delta_ne_zero(hΔ、第 977、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateModel_map_Delta_ne_zero") 1,
    .citation "[ABC3]" "exists_veluW_mu(hvw の w、第 961、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_veluW_mu") 1,
    .citation "[ABC3]" "isElliptic_veluCurve_tate_map(hvw の楼円性、第 962、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isElliptic_veluCurve_tate_map") 1,
    .citation "[ABC3]" "splits_or_exists_twist_splits''(分裂性の二者択一、第 982、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.splits_or_exists_twist_splits''") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_nonsplit(非分裂の降下、第 929、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_nonsplit") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 944-985）の測定**——" ++
       "「大域データ（Q・hQ・hE′）＋局所データ → 結論」は" ++
       "第 972 で 1 本の証明済み定理になった。" ++
       "☆その局所データ 9 項目のうち 8 項目に供給元が付いている" ++
       "（954・964・953・978・977・961・962・970+947+971）。") 2,
    .implicitStep
      ("★★**残る 1 項目は分裂性であり、その唯一の障害は mathlib の穴である**: " ++
       "`Finite (IsLocalRing.ResidueField (p.adicCompletionIntegers L))` が無い" ++
       "（第 983・985 で実測）。" ++
       "☆第 982 が `[Fintype k]` を要求するのでこれが要る。" ++
       "★道筋: `L` の稠密性（`denseRange_algebraMap`、mathlib にあり）で " ++
       "`𝓞 L → ResidueField R` の全射性を示し、" ++
       "`Ideal.finiteQuotientOfFreeOfNeBot`（mathlib にあり）に帰着させる。" ++
       "☆第 897（IsAdicComplete の内製）と同型の作業である。") 5,
    .implicitStep
      ("☆残りの配管: (i) `C • E` を `a₁ = a₃ = 0` に正規化して `p` での整性を保つ段" ++
       "（第 981 により体の側でやればよい。`p ∣ 2` は第 979 の測定）、" ++
       "(ii) 捧り `d` を `Lv` の整数環から `L` に降ろす段、" ++
       "(iii) 972 に並べて `minDeltaExp_variableChange` で `E`・`E′` に戻す段。") 3 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 993 —— 悪い素点での分裂／捻りの二者択一 -/

open Polynomial in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点では、完備化で分裂乗法還元をもつか、
ある捻りで 2 次式が分裂するかのどちらかである**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 993）**——分裂性の連鎖の終点である:

    976（乗法還元）→ 983（剰余体で c₄ ≠ 0）→ 981（整モデルの IsCharNeTwoNF）
    → 989（剰余体の有限性・内製）→ 982（二者択一）→ 992（完備化で当てる）
    → **993（HasSplitMultiplicativeReduction に流す）**

☆捻りの側は 990 で `d` を `𝓞 L` に降ろし、929 で `Δ_min` を降ろす。 -/
theorem hasSplit_or_twist_at_bad_prime {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (ha1 : W.a₁ = 0) (ha3 : W.a₃ = 0)
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0)
    (hΔpos : 0 < valAdd p (Units.mk0 W.Δ hΔ))
    [WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      (W.baseChange (p.adicCompletion L))]
    [WeierstrassCurve.HasMultiplicativeReduction (p.adicCompletionIntegers L)
      (W.baseChange (p.adicCompletion L))] :
    WeierstrassCurve.HasSplitMultiplicativeReduction (p.adicCompletionIntegers L)
        (W.baseChange (p.adicCompletion L))
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
  rcases splits_or_twist_at_completion L p W ha1 ha3 with hs | ⟨d, hd, hs⟩
  · exact Or.inl (hasSplitMultiplicativeReduction_baseChange p W hp hΔ hc4ne hc4 hΔpos hs)
  · exact Or.inr ⟨d, hd, hs⟩

def hasSplit_or_twist_at_bad_prime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での分裂／捻りの二者択一)",
    sectionId := "genell-lemma-3-5" }

def hasSplit_or_twist_at_bad_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "splits_or_twist_at_completion(完備化での二者択一、第 992、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.splits_or_twist_at_completion") 1,
    .citation "[ABC3]" "hasSplitMultiplicativeReduction_baseChange(Splits から分裂還元へ、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.hasSplitMultiplicativeReduction_baseChange") 1 ]

end ABC3.Skeleton.GenEll
