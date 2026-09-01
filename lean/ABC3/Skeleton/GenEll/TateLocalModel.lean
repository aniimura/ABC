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
import ABC3.Found.GaloisRep.Lemma35Ineq
import ABC3.Found.GaloisRep.UnramQuad
import ABC3.Found.GaloisRep.TorsionIntegralGood

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

と全部繋がっている。

★★★★**2026-09-01（第 1123）の整理**——本ファイルにあった 3 つの `sorry` 節点は
すべて**後継の証明済み定理に置き換いった**ので削除した:

| 削除した節点 | 代わりに立っているもの |
|---|---|
| `isMinimal_baseChange` | `isMinimal_baseChange_of_c4`（第 908、`Found/GaloisRep/DegInfTateParam.lean`） |
| `minDeltaExp_descend_of_nonsplit` | `minDeltaExp_eq_mul_of_nonsplit`（第 929、`Found/GaloisRep`） |
| `isMuAtBadPrimes_of_veluQuotient` | `isMuAtBadPrimes_of_veluQuotient_of_split`（第 1005、本ファイル） |

☆これで本ファイルの `sorry` は 0 になった。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

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

/-! ## ★★★★★★★★★★★★★★★★★★★★`Lemma 3.5` が消費する型——`IsMuAtBadPrimes` -/

/-- ★**原文の「H は乗法還元の各素点で `μ_l` に対応する」の帰結**を型にしたもの。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`lemma_3_5_velu_bad_delta`（第 903、証明済み）が**そのまま消費する形**である。 -/
def IsMuAtBadPrimes {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic] (l : ℕ) : Prop :=
  ∀ p : HeightOneSpectrum (𝓞 L),
    jExp p E < 0 → minDeltaExp p E' = l * minDeltaExp p E

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

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★第 1006 —— `Lemma 3.5` の不等式まで通す

★第 1005（`IsMuAtBadPrimes`）を第 903（`lemma_3_5_velu_bad_delta`）に流す。
☆これで `Lemma 3.5` の不等式

    `(1/(12(1+ϵ)))·l·deg∞(E) ≤ ht^Falt(E) + 2log(l) + C`

が、**残る仮説の全リスト**のもとで通る。★その全リストが本ブロックの型である。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Lemma 3.5 ——
Vélu の商で受ける形（残る仮説の全リスト）**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1006）**——残る仮説は次の 3 群である:

| 群 | 仮説 | 状況 |
|---|---|---|
| 分裂性 | `hsplit` | 第 976＋993 が「分裂または捻りで分裂」まで。`p ∣ 2` の非分裂が残件 |
| `l` の位置 | `hlu` | `p ∤ l` |
| 素点ごとの整性・アルキメデス | `hmin0`・`hint`・`P`・`Cv` | 第 903 がもともと受けていたもの |

☆`hcop`（`l` は局所高さと素）と半安定性は**原文自身の仮定**である。 -/
theorem lemma_3_5_velu_mu (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime → l ≠ 2 →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (∀ p : HeightOneSpectrum (𝓞 L), IsUnit ((l : (p.adicCompletionIntegers L)))) →
      (∀ (p : HeightOneSpectrum (𝓞 L)), jExp p E < 0 →
        ∀ (C : WeierstrassCurve.VariableChange L)
          (_hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
            ((C • E).baseChange (p.adicCompletion L))),
          ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
            (p.adicCompletionIntegers L)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := ABC3.Found.GaloisRep.lemma_3_5_velu_bad_delta eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl hodd Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint
    hssE hssE' hcop hlu hsplit => ?_⟩
  exact hC L E E' l hl.pos Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint hssE hssE'
    (isMuAtBadPrimes_of_veluQuotient_of_split E E' hl hodd Q hQ hE' hssE hssE'
      hcop hlu hsplit)

def lemma_3_5_velu_mu.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商で受ける形——残る仮説は分裂性と p ∤ l だけ)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_mu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_of_split(第 1005、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_of_split") 1,
    .citation "[ABC3]" "lemma_3_5_velu_bad_delta(不等式の側、第 903、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_velu_bad_delta") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1006）の測定**——`Lemma 3.5` の不等式は" ++
       "**残る仮説の全リスト**のもとで通った。☆残るのは 2 群だけである: " ++
       "(1) `hsplit`（`p ∣ 2` の非分裂＝不分岐 2 次拡大待ち）、(2) `hlu`（`p ∤ l`）。" ++
       "★その他（`P`・`Cv`・`neronExp = 0`・`IsIntegral`）は第 903 が" ++
       "もともと受けていたものである") 2 ]

/-! ## ★★★★★★★★★★★★★★★★第 1009 —— 非分裂を含めた節点

★第 1004／1005 は `hsplit`（完備化で**分裂**乗法還元）を仮説に置いていた。
☆本節点はそれを**外した形**であり、残る仮説は `hlu`（`p ∤ l`）だけである。

★★分裂の場合は**すでに証明済み**（第 1004）。残るのは非分裂の場合で、
古典的には**不分岐 2 次拡大**に上げて分裂に帰着する:

| 段 | 状況 |
|---|---|
| 捻り `d`（完備化の整数環の単元）を取る | ★第 992／993 で**済** |
| `d` が平方でなければ `K` でも平方でない | ★第 1007 で**済** |
| `X² − d` が既約 → `K(√d)` は 2 次拡大 | ★第 1008 で**済** |
| `R[√d]` が完備離散付値環（`e = 1`） | ☆残 |
| 付値の橋 `hp′` が `L` 上でそのまま | ☆残（不分岐なので `e = 1`） |
| そこで分裂するので第 997 が当たる | ☆残 |

☆第 997（`jExp_eq_mul_of_j_tate_pow`）は `Lv`・`R` について**一般**なので、
上の 3 段が済めばそのまま当たる。 -/

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係
（分裂性を仮定しない形）**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1009）**——分裂の場合は第 1004 で証明済み。
☆非分裂の場合が不分岐 2 次拡大待ちである（葉は第 1007／1008 で 2 つ済んだ）。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_any {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletion L p
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletionIntegers L p
  obtain ⟨C, hC, hc4ne, hc4⟩ :=
    ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg p E hssE hjneg
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hmin :=
    ABC3.Found.GaloisRep.isMinimal_baseChange_at_bad_prime p hp E C hC hc4ne hc4
  haveI hmult :=
    ABC3.Found.GaloisRep.hasMultiplicativeReduction_at_bad_prime p hp E C hC hc4ne hc4 hjneg
  have hA := ABC3.Found.GaloisRep.integralModel_c4_isUnit
    (R := p.adicCompletionIntegers L) ((C • E).baseChange (p.adicCompletion L))
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L))
  · exact minDeltaExp_eq_mul_at_bad_prime_full p E E' hssE hssE' hjneg hl hodd hcop
      C hmin hs hlu hQ hE'
  · -- ★非分裂——不分岐 2 次拡大を通す（第 1034）
    have hirr := ABC3.Found.GaloisRep.irreducible_map_residue_of_not_splits
      (ABC3.Found.GaloisRep.monic_splitQuadPoly _ hA)
      (ABC3.Found.GaloisRep.natDegree_splitQuadPoly _ hA)
      (fun hsp => hs (ABC3.Found.GaloisRep.hasSplit_of_splits_splitQuadPoly _ hA hsp))
    exact minDeltaExp_eq_mul_at_bad_prime_nonsplit p hp E E' C hC hc4ne hc4 hA hirr
      hssE hssE' hjneg hl hodd hcop hlu hQ hE'

def minDeltaExp_eq_mul_at_bad_prime_any.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——分裂性を仮定しない形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_any.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_full(分裂の場合、第 1004、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_full") 1,
    .citation "[ABC3]" "splits_or_twist_at_completion(捻り d を取る、第 992、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.splits_or_twist_at_completion") 1,
    .citation "[ABC3]" "not_isSquare_in_fractionField(d は K でも平方でない、第 1007、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_isSquare_in_fractionField") 1,
    .citation "[ABC3]"
      "irreducible_X_sq_sub_C_fractionField(X² − d は既約、第 1008、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.irreducible_X_sq_sub_C_fractionField") 1,
    .citation "[ABC3]" "jExp_eq_mul_of_j_tate_pow(Lv・R について一般、第 997、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.jExp_eq_mul_of_j_tate_pow") 1,
    .implicitStep
      ("☆残る 3 段: (1) `R[√d]` が完備離散付値環であること" ++
       "（`d` は単元なので不分岐、`𝔪R[√d]` が極大で主。" ++
       "mathlib は `integralClosure.isDedekindDomain` を持つ）、" ++
       "(2) 付値の橋 `hp′`（不分岐なので `e = 1`）、" ++
       "(3) そこで分裂すること（剰余体の 2 次拡大で 2 次式が分裂する）") 3 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] `IsMuAtBadPrimes`——
残る仮説は `p ∤ l` の 1 本だけ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1009）**——第 1005 から `hsplit` を外した形である。 -/
theorem isMuAtBadPrimes_of_veluQuotient_of_coprime {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hssE : ∀ p, SemistableAt p E) (hssE' : ∀ p, SemistableAt p E')
    (hcop : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E))
    (hlu : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → IsUnit ((l : (p.adicCompletionIntegers L)))) :
    IsMuAtBadPrimes E E' l := fun p hbad =>
  minDeltaExp_eq_mul_at_bad_prime_any p E E' (hssE p) (hssE' p) hbad hl hodd
    (hcop p hbad) (hlu p hbad) hQ hE'

def isMuAtBadPrimes_of_veluQuotient_of_coprime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsMuAtBadPrimes——残る仮説は p ∤ l の 1 本)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient_of_coprime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_any(第 1009)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_any") 1 ]

/-! ## ★★★★★★★★★★★★★★★★第 1019 —— `hlu`（`p ∤ l`）の出どころ

★第 1009 に残った 2 本目の仮説は `hlu : IsUnit ((l : R_p))` である。
☆`R_p` は局所環なので、これは **`l` が `p` の剰余標数と異なる**ことに等しい
（`l` は素数なので `char k_p ∣ l` ⟺ `char k_p = l`）。

★★**なぜ落ちるか**——悪い素点で `l` 捉れ点があると、第 946／947 により
`μ_l ⊂ Lv` である（`hcop` が効いて `x^l = q^k` の `k` が `l` で割れる）。
☆ところが `l` が剰余標数なら `μ_l ⊂ Lv` は分岐指数 `e ≥ l − 1` を要求する
（`ℚ_l(ζ_l)/ℚ_l` は完全分岐で次数 `l − 1`）。
★したがって `[L:ℚ] ≥ [Lv:ℚ_l] ≥ l − 1` となり、`l` が `[L:ℚ] + 1` より大きければ矛盾する。

☆`Lemma 3.7`／`Theorem 3.8` は `l ≥ 100d·(…)`（`d = [L:ℚ]`）を仮定するので、
**この条件は原文の仮定から出る**。 -/

/-- ★★★★**局所環では `l` が単元 ⟺ `l` が極大イデアルに入らない**。 -/
theorem isUnit_natCast_iff_notMem {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (l : ℕ) :
    IsUnit ((l : (p.adicCompletionIntegers L)))
      ↔ (l : (p.adicCompletionIntegers L)) ∉
        IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) :=
  (IsLocalRing.notMem_maximalIdeal).symm

/-! ## ★★★★★★★★第 1038 —— 残るのは分岐指数の上界だけ

★第 1037 で `l ∈ 𝔪^{l−1}`（＝ `e_p ≥ l − 1`）まで来た。
☆あとは **`e_p ≤ [L:ℚ]`** を言えば `l ≤ [L:ℚ] + 1` となり、
`hd : [L:ℚ] + 1 < l` に矛盾する。

★mathlib は `Ideal.sum_ramification_inertia`（`∑_{P|p} e_P f_P = [L:K]`）を持つ。
☆`Ideal.ramificationIdx` は `sSup {n | map p ≤ P^n}` なので、
`l ∈ p^{l−1}` から `l − 1 ≤ e_p` を出すには**上に有界**であることが要る。 -/

/-! ## ★★★★★★★★★★★★第 1040 —— 完備化から `𝓞 L` へ降ろす

★第 1037 は `l ∈ 𝔪_R^{l−1}`（`R` は完備化の整数環）を与える。
☆第 1039 が要るのは `l ∈ p^{l−1}`（`𝓞 L` の側）である。
★両者は付値で繋がる——`intValuation` は `hp`（第 964）で一致する。 -/

open IsDedekindDomain NumberField in
/-- ★★★★★★★★**自然数の `intValuation` は完備化で変わらない**（第 1040）。 -/
theorem intValuation_natCast_completion {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (m : ℕ) :
    (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)).intValuation
        ((m : ℕ) : (p.adicCompletionIntegers L))
      = p.intValuation ((m : ℕ) : 𝓞 L) := by
  rw [← HeightOneSpectrum.valuation_of_algebraMap (K := p.adicCompletion L),
    ← HeightOneSpectrum.valuation_of_algebraMap (K := L)]
  rw [← ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
    ((algebraMap (𝓞 L) L) ((m : ℕ) : 𝓞 L))]
  congr 1
  push_cast
  ring

open IsDedekindDomain NumberField in
/-- ★★★★★★★★★★★★**`𝔪_R^k` の元は `p^k` の元**（自然数について、第 1040）。 -/
theorem natCast_mem_pow_of_mem_pow_completion {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) {m k : ℕ}
    (h : ((m : ℕ) : (p.adicCompletionIntegers L))
      ∈ (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)).asIdeal ^ k) :
    ((m : ℕ) : 𝓞 L) ∈ p.asIdeal ^ k := by
  rw [← HeightOneSpectrum.intValuation_le_pow_iff_mem] at h ⊢
  rwa [intValuation_natCast_completion] at h

def intValuation_natCast_completion.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(自然数の intValuation は完備化で変わらない。★無条件)",
    sectionId := "genell-lemma-3-7" }

def intValuation_natCast_completion.needs : List ProofObligation := []

def natCast_mem_pow_of_mem_pow_completion.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(𝔪_R^k の元は p^k の元——自然数について。★無条件)",
    sectionId := "genell-lemma-3-7" }

def natCast_mem_pow_of_mem_pow_completion.needs : List ProofObligation := []

/-- ☆**節点**: `p^{l−1}` が `l` を含むなら `l − 1 ≤ [L:ℚ]`。

★中身は `∑_{P ∣ l} e_P f_P = [L:ℚ]`（mathlib `Ideal.sum_ramification_inertia`）と
`e_p ≥ l − 1`（第 1037）である。 -/
theorem le_finrank_of_natCast_mem_pow {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) {l : ℕ} (hl : l.Prime)
    (hmem : ((l : ℕ) : 𝓞 L) ∈ p.asIdeal ^ (l - 1)) :
    l - 1 ≤ Module.finrank ℚ L := by
  have hnorm : Ideal.absNorm (Ideal.span {((l : ℕ) : 𝓞 L)}) = l ^ Module.finrank ℚ L := by
    rw [Ideal.absNorm_span_singleton]
    have h1 : ((l : ℕ) : 𝓞 L) = algebraMap ℤ (𝓞 L) (l : ℤ) := by push_cast; rfl
    rw [h1, Algebra.norm_algebraMap, NumberField.RingOfIntegers.rank]
    simp
  have hle : Ideal.span {((l : ℕ) : 𝓞 L)} ≤ p.asIdeal ^ (l - 1) := by
    rw [Ideal.span_le, Set.singleton_subset_iff]; exact hmem
  have hle2 : Ideal.span {((l : ℕ) : 𝓞 L)} ≤ p.asIdeal :=
    le_trans hle (Ideal.pow_le_self (by have := hl.two_le; omega))
  have hdvd2 : Ideal.absNorm p.asIdeal ∣ l ^ Module.finrank ℚ L := by
    rw [← hnorm]; exact Ideal.absNorm_dvd_absNorm_of_le hle2
  have hnetop : p.asIdeal ≠ ⊤ := p.isPrime.ne_top
  have hne1 : Ideal.absNorm p.asIdeal ≠ 1 := by
    rw [Ne, Ideal.absNorm_eq_one_iff]; exact hnetop
  have hldvd : l ∣ Ideal.absNorm p.asIdeal := by
    obtain ⟨k, hk, hkk⟩ := (Nat.dvd_prime_pow hl).1 hdvd2
    rcases Nat.eq_zero_or_pos k with rfl | hkpos
    · rw [pow_zero] at hkk; exact absurd hkk hne1
    · rw [hkk]; exact dvd_pow_self l hkpos.ne'
  have hdvd1 : Ideal.absNorm (p.asIdeal ^ (l - 1)) ∣ l ^ Module.finrank ℚ L := by
    rw [← hnorm]; exact Ideal.absNorm_dvd_absNorm_of_le hle
  rw [map_pow] at hdvd1
  exact (Nat.pow_dvd_pow_iff_le_right hl.one_lt).1
    (dvd_trans (pow_dvd_pow_of_dvd hldvd (l - 1)) hdvd1)

def le_finrank_of_natCast_mem_pow.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(p^{l−1} が l を含むなら l − 1 ≤ [L:ℚ])",
    sectionId := "genell-lemma-3-7" }

def le_finrank_of_natCast_mem_pow.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Ideal.absNorm_dvd_absNorm_of_le(イデアルノルムは包含で割る)"
      (.inMathlib "Ideal.absNorm_dvd_absNorm_of_le") 1,
    .citation "[mathlib]" "Ideal.absNorm_span_singleton(単元生成イデアルのノルム)"
      (.inMathlib "Ideal.absNorm_span_singleton") 1 ]

/-! ## ★★★★★★★★★★★★★★★★第 1042 —— 分裂の場合の `hlu`

★第 1037（`l ∈ 𝔪_R^{l−1}`）→ 第 1040（`l ∈ p^{l−1}`）→ 第 1039（`l − 1 ≤ [L:ℚ]`）
の連鎖に、**入口**（`ζ` の生成）を繋ぐ。

☆入口は第 1041（Tate モデルの上の位数 `l` の点）＋第 947（点から原始 `l` 乗根へ）である。 -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**分裂の場合の `hlu`**（第 1042）。

☆`l` が単元でないとすると `μ_l ⊂ R` から `l ∈ 𝔪^{l−1}`（第 1037）、
`𝓞 L` に降ろして（第 1040）`l − 1 ≤ [L:ℚ]`（第 1039）——`hd` に矛盾する。 -/
theorem isUnit_natCast_of_split {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    {l : ℕ} (hl : l.Prime) (hd : Module.finrank ℚ L + 1 < l)
    (hbad : jExp p E < 0) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l) :
    IsUnit ((l : (p.adicCompletionIntegers L))) := by
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletion L p
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletionIntegers L p
  by_contra hnu
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  have hΔ := ABC3.Found.GaloisRep.tateModel_map_Delta_ne_zero
    (E.baseChange (p.adicCompletion L)) h
  have hcop' :=
    ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp E h hbad hcop
  obtain ⟨P, hP, hP0⟩ := ABC3.Found.GaloisRep.exists_point_tateModel E h hl hQ
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hcop' P hP hP0
  have hmem :=
    ABC3.Found.GaloisRep.natCast_mem_maximalIdeal_pow_of_not_isUnit hl hζ hnu
  have hmem2 := natCast_mem_pow_of_mem_pow_completion p (m := l) (k := l - 1) hmem
  have hle := le_finrank_of_natCast_mem_pow p hl hmem2
  omega

def isUnit_natCast_of_split.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(分裂の場合の hlu——l が定義体の次数より大きければ単元)",
    sectionId := "genell-lemma-3-7" }

def isUnit_natCast_of_split.needs : List ProofObligation :=
  [ .citation "[ABC3]" "natCast_mem_maximalIdeal_pow_of_not_isUnit(第 1037、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_mem_maximalIdeal_pow_of_not_isUnit") 1,
    .citation "[ABC3]" "natCast_mem_pow_of_mem_pow_completion(第 1040、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.natCast_mem_pow_of_mem_pow_completion") 1,
    .citation "[ABC3]" "le_finrank_of_natCast_mem_pow(第 1039、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.le_finrank_of_natCast_mem_pow") 1 ]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**[GenEll] 悪い素点で `l` は単元**——
`l` が定義体の次数より十分大きければ。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★★★**2026-09-01（第 1019）の測定**——`hlu` の出どころを型で固定した。
☆中身は「`l` 捉れ点 ⇒ `μ_l ⊂ Lv`（第 946／947）」と
「`μ_l ⊂ Lv` かつ `l` が剰余標数 ⇒ `e ≥ l − 1`」の 2 段である。 -/
theorem isUnit_natCast_at_bad_prime {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hd : Module.finrank ℚ L + 1 < l)
    (hssE : SemistableAt p E)
    (hbad : jExp p E < 0) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    IsUnit ((l : (p.adicCompletionIntegers L))) := by
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletion L p
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletionIntegers L p
  obtain ⟨C, hC, hc4ne, hc4⟩ :=
    ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg p E hssE hbad
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hmin :=
    ABC3.Found.GaloisRep.isMinimal_baseChange_at_bad_prime p hp E C hC hc4ne hc4
  haveI hmult :=
    ABC3.Found.GaloisRep.hasMultiplicativeReduction_at_bad_prime p hp E C hC hc4ne hc4 hbad
  have hjC : jExp p (C • E) < 0 := by
    rw [ABC3.Found.GaloisRep.jExp_variableChange p E C]; exact hbad
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • E)) := by
    rw [ABC3.Found.GaloisRep.jExp_variableChange p E C]; exact hcop
  have hQ' : addOrderOf (ABC3.Found.GenEll.vcPoint C E Q) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint C E Q]; exact hQ
  have hA := ABC3.Found.GaloisRep.integralModel_c4_isUnit
    (R := p.adicCompletionIntegers L) ((C • E).baseChange (p.adicCompletion L))
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L))
  · exact isUnit_natCast_of_split p (C • E) hs hl hd hjC hcopC hQ'
  · -- ★非分裂——不分岐 2 次拡大で ζ を作り、付値で降ろす
    by_contra hnu
    have hfm := ABC3.Found.GaloisRep.monic_splitQuadPoly _ hA
    have hfd := ABC3.Found.GaloisRep.natDegree_splitQuadPoly _ hA
    have hirr := ABC3.Found.GaloisRep.irreducible_map_residue_of_not_splits hfm hfd
      (fun hsp => hs (ABC3.Found.GaloisRep.hasSplit_of_splits_splitQuadPoly _ hA hsp))
    haveI hdom := ABC3.Found.GaloisRep.isDomain_adjoinRoot hfm hirr
    haveI hlocR := ABC3.Found.GaloisRep.isLocalRing_adjoinRoot hfm hfd hirr
    haveI hdvr := ABC3.Found.GaloisRep.isDiscreteValuationRing_adjoinRoot hfm hfd hirr
    haveI hcomp0 : IsAdicComplete
        (IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
        (AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
          (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
            ((C • E).baseChange (p.adicCompletion L))) hA)) :=
      ABC3.Found.GaloisRep.isAdicComplete_adjoinRoot _ hfm hfd
    haveI hcomp : IsAdicComplete
        (IsLocalRing.maximalIdeal (AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
          (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
            ((C • E).baseChange (p.adicCompletion L))) hA)))
        (AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
          (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
            ((C • E).baseChange (p.adicCompletion L))) hA)) := by
      rw [ABC3.Found.GaloisRep.maximalIdeal_adjoinRoot_eq_map hfm hfd hirr]
      exact ABC3.Found.GaloisRep.isAdicComplete_map_algebraMap _
    haveI hchar0 : CharZero (AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
        (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
          ((C • E).baseChange (p.adicCompletion L))) hA)) :=
      charZero_of_injective_algebraMap
        (ABC3.Found.GaloisRep.algebraMap_adjoinRoot_injective hfm hfd)
    haveI hchar0K : CharZero (FractionRing (AdjoinRoot
        (ABC3.Found.GaloisRep.splitQuadPoly
          (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
            ((C • E).baseChange (p.adicCompletion L))) hA))) :=
      IsFractionRing.charZero _
    letI : Algebra L (FractionRing (AdjoinRoot
        (ABC3.Found.GaloisRep.splitQuadPoly (WeierstrassCurve.integralModel
          (p.adicCompletionIntegers L)
          ((C • E).baseChange (p.adicCompletion L))) hA))) :=
      ((ABC3.Found.GaloisRep.quadFieldHom (K := p.adicCompletion L) hfm hfd).comp
        (algebraMap L (p.adicCompletion L))).toAlgebra
    have halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot
        (ABC3.Found.GaloisRep.splitQuadPoly (WeierstrassCurve.integralModel
          (p.adicCompletionIntegers L)
          ((C • E).baseChange (p.adicCompletion L))) hA))) x
        = ABC3.Found.GaloisRep.quadFieldHom hfm hfd
          (algebraMap L (p.adicCompletion L) x) := fun _ => rfl
    haveI hmin' :=
      ABC3.Found.GaloisRep.isMinimal_baseChange_ext p hp hfm hfd hirr halg E C hC hc4ne hc4
    haveI hmult' := ABC3.Found.GaloisRep.hasMultiplicativeReduction_ext p hp hfm hfd hirr
      halg E C hC hc4ne hc4 hbad
    have h' := ABC3.Found.GaloisRep.hasSplitMultiplicativeReduction_ext (C • E) hA rfl halg
    have hp' := ABC3.Found.GaloisRep.valuation_algebraMap_ext p hp hfm hfd hirr halg
    have hΔ' := ABC3.Found.GaloisRep.tateModel_map_Delta_ne_zero
      ((C • E).baseChange (FractionRing (AdjoinRoot _))) h'
    letI : DecidableEq (FractionRing (AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
      (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
        ((C • E).baseChange (p.adicCompletion L))) hA))) := Classical.decEq _
    have hcop'' := ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp'
      (C • E) h' hjC hcopC
    obtain ⟨P, hP, hP0⟩ := ABC3.Found.GaloisRep.exists_point_tateModel (C • E) h' hl hQ'
    obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
      ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point
        (tateParamR ((C • E).baseChange (FractionRing (AdjoinRoot _))) h')
        (tateParamR_mem ((C • E).baseChange (FractionRing (AdjoinRoot _))) h')
        (tateParamR_ne_zero ((C • E).baseChange (FractionRing (AdjoinRoot _))) h')
        hΔ' hl hcop'' P hP hP0
    have hlm : ((l : ℕ) : (p.adicCompletionIntegers L))
        ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) := by
      by_contra hc
      exact hnu (IsLocalRing.notMem_maximalIdeal.1 hc)
    have hnu' : ¬ IsUnit ((l : AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
        (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
          ((C • E).baseChange (p.adicCompletion L))) hA))) := by
      intro hu
      refine IsLocalRing.notMem_maximalIdeal.2 hu ?_
      have hcast : ((l : ℕ) : AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
        (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
          ((C • E).baseChange (p.adicCompletion L))) hA))
          = algebraMap (p.adicCompletionIntegers L)
            (AdjoinRoot (ABC3.Found.GaloisRep.splitQuadPoly
        (WeierstrassCurve.integralModel (p.adicCompletionIntegers L)
          ((C • E).baseChange (p.adicCompletion L))) hA))
            ((l : ℕ) : (p.adicCompletionIntegers L)) := (map_natCast _ l).symm
      rw [ABC3.Found.GaloisRep.maximalIdeal_adjoinRoot_eq_map hfm hfd hirr, hcast]
      exact Ideal.mem_map_of_mem _ hlm
    have hmemR' :=
      ABC3.Found.GaloisRep.natCast_mem_maximalIdeal_pow_of_not_isUnit hl hζ hnu'
    have hmemR := ABC3.Found.GaloisRep.natCast_mem_pow_of_mem_pow_quadExt
      (K := p.adicCompletion L) hfm hfd hirr hmemR'
    have hmem2 := natCast_mem_pow_of_mem_pow_completion p (m := l) (k := l - 1) hmemR
    have hle := le_finrank_of_natCast_mem_pow p hl hmem2
    omega

def isUnit_natCast_iff_notMem.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(局所環では l が単元 ⟺ l が極大イデアルに入らない。★無条件)",
    sectionId := "genell-lemma-3-7" }

def isUnit_natCast_iff_notMem.needs : List ProofObligation := []

def isUnit_natCast_at_bad_prime.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(悪い素点で l は単元——l が定義体の次数より大きければ)",
    sectionId := "genell-lemma-3-7" }

def isUnit_natCast_at_bad_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_mu_point_dvr(l 捉れ点は μ_l に対応する、第 946、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_mu_point_dvr") 1,
    .implicitStep
      ("☆`μ_l ⊂ Lv` かつ `l` が剰余標数なら分岐指数は `e ≥ l − 1`" ++
       "（`ℚ_l(ζ_l)/ℚ_l` は完全分岐で次数 `l − 1`）。" ++
       "★したがって `[L:ℚ] ≥ [Lv:ℚ_l] ≥ l − 1` で `hd` に矛盾する。" ++
       "☆mathlib に `ℚ_p(ζ_p)` の分岐は無い（2026-09-01 実測）") 3 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] `IsMuAtBadPrimes`——
残る仮説は原文自身の仮定だけ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1020）**——第 1009 の `hlu` を第 1019 で供給した形。
☆`hd : [L:ℚ] + 1 < l` は `Lemma 3.7`／`Theorem 3.8` が仮定する
`l ≥ 100d·(…)` から出るので、**原文にない仮説は無くなった**。 -/
theorem isMuAtBadPrimes_of_veluQuotient_of_large {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hd : Module.finrank ℚ L + 1 < l)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hssE : ∀ p, SemistableAt p E) (hssE' : ∀ p, SemistableAt p E')
    (hcop : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) :
    IsMuAtBadPrimes E E' l :=
  isMuAtBadPrimes_of_veluQuotient_of_coprime E E' hl hodd Q hQ hE' hssE hssE' hcop
    (fun p hbad => isUnit_natCast_at_bad_prime p E hl hd (hssE p) hbad (hcop p hbad) Q hQ)

def isMuAtBadPrimes_of_veluQuotient_of_large.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsMuAtBadPrimes——残る仮説は原文自身の仮定だけ)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient_of_large.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_of_coprime(第 1009)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_of_coprime") 1,
    .citation "[ABC3]" "isUnit_natCast_at_bad_prime(hlu の供給、第 1019)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isUnit_natCast_at_bad_prime") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★第 1045 —— `Lemma 3.5` の不等式（到達形）

★第 1044 で `IsMuAtBadPrimes` が `sorry`-free になったので、
第 1006 の `hsplit`・`hlu` は**消える**。

☆残る仮説はすべて**原文自身の仮定**か、第 903 がもともと受けていた
アルキメデスの周期対である。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Lemma 3.5 の不等式**
（`IsMuAtBadPrimes` は証明済み、第 1045）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1045）**——`hd : [L:ℚ] + 1 < l` は
`Lemma 3.7`／`Theorem 3.8` の `l ≥ 100d·(…)` から出る。
☆したがって残るのは**アルキメデスの周期対**（`P`・`Cv`）と
素点ごとの整性（`neronExp = 0`・`E′.IsIntegral`）だけである。 -/
theorem lemma_3_5_velu_large (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime → l ≠ 2 →
      Module.finrank ℚ L + 1 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := ABC3.Found.GaloisRep.lemma_3_5_velu_bad_delta eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl hodd hd Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint
    hssE hssE' hcop => ?_⟩
  exact hC L E E' l hl.pos Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint hssE hssE'
    (isMuAtBadPrimes_of_veluQuotient_of_large E E' hl hodd hd Q hQ hE' hssE hssE' hcop)

def lemma_3_5_velu_large.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(不等式——IsMuAtBadPrimes は証明済み、残るはアルキメデスの周期対)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_large.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_of_large(第 1044、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_of_large") 1,
    .citation "[ABC3]" "lemma_3_5_velu_bad_delta(不等式の側、第 903、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_velu_bad_delta") 1 ]

/-! ## ★★★★★★★★★★★★★★★★第 1046 —— アルキメデスの周期対は自前で作れる

★第 1045 は `P`・`Cv`・`hΔ`・`hPC`・`hell1`・`hell2` を仮説で受けていた。
☆しかしこれらは第 348（`exists_periodPair_of_isElliptic`、**無条件**）から
`choose` で作れる——複素楕円曲線はつねに格子で一意化されるからである。 -/

open ABC3.Found.GenEll in
/-- ★★★★★★★★★★★★★★★★**アルキメデスの周期対はつねに取れる**（第 1046）。 -/
theorem exists_periodPair_family {L : Type} [Field L] [NumberField L]
    (E : WeierstrassCurve L) [E.IsElliptic] :
    ∃ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → WeierstrassCurve.VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0)
      ∧ (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ))
      ∧ (∀ σ : L →+* ℂ, (E.map σ).IsElliptic)
      ∧ (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) := by
  have hell1 : ∀ σ : L →+* ℂ, (E.map σ).IsElliptic := by
    intro σ
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.map_Δ, isUnit_iff_ne_zero]
    exact (map_ne_zero_iff _ σ.injective).2 E.isUnit_Δ.ne_zero
  choose P Cv hPC using fun σ : L →+* ℂ => exists_periodPair_of_isElliptic (E.map σ) (hell1 σ)
  refine ⟨P, Cv, ?_, hPC, hell1, ?_⟩
  · intro σ
    haveI := hell1 σ
    haveI : (latticeCurve (P σ)).IsElliptic := by
      rw [← hPC σ, WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
      exact (((Cv σ).u⁻¹).isUnit.pow 12).mul (E.map σ).isUnit_Δ
    exact latticeDisc_ne_zero_of_isElliptic (P σ)
  · intro σ
    haveI := hell1 σ
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact (((Cv σ).u⁻¹).isUnit.pow 12).mul (E.map σ).isUnit_Δ

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Lemma 3.5 の不等式**
（アルキメデスの周期対も自前、第 1046）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1046）**——残る仮説は
`neronExp p E = 0` と `E′.IsIntegral (primeSubring p)` の 2 本だけになった。 -/
theorem lemma_3_5_velu_arch (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime → l ≠ 2 →
      Module.finrank ℚ L + 1 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_large eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl hodd hd Q hQ hE' hmin0 hint hssE hssE' hcop => ?_⟩
  obtain ⟨P, Cv, hΔ, hPC, hell1, hell2⟩ := exists_periodPair_family E
  exact hC L E E' l hl hodd hd Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint hssE hssE' hcop

def exists_periodPair_family.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(アルキメデスの周期対はつねに取れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_periodPair_family.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_periodPair_of_isElliptic(第 348、無条件、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_periodPair_of_isElliptic") 1 ]

def lemma_3_5_velu_arch.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(不等式——残るは neronExp = 0 と E′ の整性の 2 本)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_arch.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu_large(第 1045、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.lemma_3_5_velu_large") 1,
    .citation "[ABC3]" "exists_periodPair_family(第 1046、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_periodPair_family") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★第 1050 —— `Lemma 3.5` の最終形

★第 1049 で有限素点側が **1 本の不等式 `hfin`** になった。
☆連鎖は浅い（`lemma_3_5_of_isogeny_estimate_le` が `hfalt` を直接受ける）ので、
一気に組める。

★★これで `Lemma 3.5` の仮説は
**原文自身の仮定 ＋ `hfin` の 1 本**だけになる。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Lemma 3.5**
——残る仮説は `hfin` の 1 本（第 1050）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `(1/(12(1+ϵ)))·l·deg_∞(E) ≤ ht^Falt(E) + 2·log(l) + C`

★★★★**2026-09-01（第 1050）**——`IsMuAtBadPrimes`（第 1044）と
アルキメデスの周期対（第 1046）と有限素点側の 1 本化（第 1049）を合わせた形。
☆`hfin` は「`E′` のモデルは `E` のモデルより極小から遠くない」であり、
第 704 の `hmin`・`hint` はその十分条件にすぎなかった。 -/
theorem lemma_3_5_velu_defect (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime → l ≠ 2 →
      Module.finrank ℚ L + 1 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ((∑ᶠ p : HeightOneSpectrum (𝓞 L),
            (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
          - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
            (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (3 / 2) * (Module.finrank ℚ L : ℝ) * Real.log l) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := ABC3.Found.GaloisRep.lemma_3_5_of_isogeny_estimate_le eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl hodd hd Q hQ hE' hfin hssE hssE' hcop => ?_⟩
  obtain ⟨P, Cv, hΔ, hPC, hell1, hell2⟩ := exists_periodPair_family E
  refine hC L E E' l hssE' ?_
    (ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le_of_defect E E' l hl.pos Q hQ hE'
      P Cv hΔ hPC hell1 hell2 hfin)
  refine ABC3.Found.GaloisRep.degInfOf_ge_of_local E E' l (fun p => ?_)
  exact ABC3.Found.GaloisRep.minDeltaExp_le_of_bad_delta p E E' (hssE p) l
    (fun hb => isMuAtBadPrimes_of_veluQuotient_of_large E E' hl hodd hd Q hQ hE'
      hssE hssE' hcop p hb)

def lemma_3_5_velu_defect.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る仮説は hfin の 1 本——E′ のモデルは E より極小から遠くない)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_defect.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_of_large(第 1044、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_of_large") 1,
    .citation "[ABC3]" "htFalt_veluQuotientFull_le_of_defect(第 1049、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le_of_defect") 1,
    .citation "[ABC3]" "exists_periodPair_family(アルキメデス、第 1046、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_periodPair_family") 1 ]

/-! ## ★★★★★★★★第 1052 —— `hfin` の残りを 2 つの節点に割る

★第 1051 の分解により、`hfin` は次の 2 つから出る:

| 節点 | 内容 |
|---|---|
| `hmin` | `E` のモデルが各素点で極小（`∀ p, neronExp p E = 0`）——モデルの取り方 |
| `hint` | `E′`（Vélu の商）が各素点で整 |

☆`hint` は古典的には「剰余標数と素な位数の捻れ点は極小モデルで整」から出る。
★そして `l ≠ char(k_p)` は**第 1044 で証明済み**（`isUnit_natCast_at_bad_prime`）である。 -/

open Finset in
/-- ☆**節点**: `l` が剰余標数と異なれば Vélu の商は整である。

★中身は「位数が剰余標数と素な捻れ点は極小モデルで整座標をもつ」であり、
Vélu の係数（`veluVFull`・`veluWFull`）はその多項式だから整になる。
☆`l ≠ char(k_p)` は第 1044（`isUnit_natCast_at_bad_prime`）が与える。 -/
theorem isIntegral_veluQuotientFull_of_coprime {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hlu : IsUnit ((l : (primeSubring p))))
    (hEint : E.IsIntegral (primeSubring p))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))) :
    E'.IsIntegral (primeSubring p) := by
  subst hE'
  haveI := hEint
  exact ABC3.Found.GaloisRep.isIntegral_veluQuotientFull_of_addOrderOf_prime p E hl hodd hlu Q hQ


/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★第 1084 —— `hfin` を外した形 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Lemma 3.5 —— 曲線の水準で無条件**（第 1084）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1084）**——第 1083 で `hfin` が無条件に出たので、
第 1050 の最後の仮説が消えた。
☆残る仮説は「`l` が奇素数で `d+1 < l`」「`Q` の位数が `l`」
「`E′` が Vélu の商」「半安定」「`l ∤ v_p(j)`」だけである。 -/
theorem lemma_3_5_velu (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime → l ≠ 2 →
      Module.finrank ℚ L + 1 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_defect eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl hodd hd Q hQ hE' hssE hssE' hcop => ?_⟩
  exact hC L E E' l hl hodd hd Q hQ hE'
    (ABC3.Found.GaloisRep.hfin_of_veluQuotientFull E E' hl hodd Q hQ hE') hssE hssE' hcop

def lemma_3_5_velu.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(曲線の水準で無条件——hfin を外した形)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu_defect(第 1050、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.lemma_3_5_velu_defect") 1,
    .citation "[ABC3]" "hfin_of_veluQuotientFull(hfin そのもの、第 1083、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hfin_of_veluQuotientFull") 1 ]

def isIntegral_veluQuotientFull_of_coprime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l が剰余標数と異なれば Vélu の商は整)",
    sectionId := "genell-lemma-3-5" }

def isIntegral_veluQuotientFull_of_coprime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isUnit_natCast_at_bad_prime(l ≠ char(k_p)、第 1044、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isUnit_natCast_at_bad_prime") 1,
    .citation "[ABC3]" "exists_vw_tate_mu(悪い素点側はこれで出る、第 1003、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_vw_tate_mu") 1,
    .implicitStep
      ("★★**2026-09-01（第 1064）の測定**——本節点は 2 つに割れる。" ++
       "☆**悪い素点**ではすでに出る——第 1003 が Vélu の係数 `v`・`w` を" ++
       "**`R` の元として**与えるので `veluCurve (tate q) v w` は `R`-曲線であり、" ++
       "第 1058 の輸送と `C₀ : VariableChange R`（`u` は単元）で `E′ ⊗ Lv` に移る。" ++
       "★残るのは**良い素点**だけで、そこでは Tate の機械が使えないので" ++
       "「位数が剰余標数と素な捧れ点は極小モデルで整」（形式群）が要る") 2 ]

/-! ## ★★★★★★★★★★★★★★★★第 1053 —— `hfin` はモデルに依らない

★★★★**測定（2026-09-01、第 1053）**——`hfin` の左辺
`∑ᶠ p, neronExp p E · log N(p) − ∑ᶠ p, neronExp p E′ · log N(p)` は
**`E` と `E′` を同じ変数変換で動かしても変わらない**。

☆`neronExp p (C • W) = neronExp p W − valAdd p C.u`（`neronExp_variableChange`）で、
`E` と `E′` の両方から同じ `valAdd p C.u` が引かれるからである。

★★そして第 969（`veluQuotientFull_vcPoint_eq`）により
`C • E′ = veluQuotientFull (C • E) (vcPoint C E Q の集合)` なので、
**`E` を各素点で極小なモデルに取り替えてよい**。
☆したがって第 1050 の `hfin` に `hmin`（大域極小モデルの存在——類数の障害がある）は要らない。 -/

/-- ★★★★★★★★★★★★★★★★**`neronExp` の差は変数変換で不変**（第 1053）。

☆これが「`hfin` はモデルに依らない」ことの中身である。 -/
theorem neronExp_sub_variableChange {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (W W' : WeierstrassCurve L)
    (hΔ : W.Δ ≠ 0) (hΔ' : W'.Δ ≠ 0) (C : WeierstrassCurve.VariableChange L) :
    neronExp p (C • W) - neronExp p (C • W') = neronExp p W - neronExp p W' := by
  rw [ABC3.Found.GaloisRep.neronExp_variableChange p W hΔ C,
    ABC3.Found.GaloisRep.neronExp_variableChange p W' hΔ' C]
  ring

def neronExp_sub_variableChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(neronExp の差は変数変換で不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def neronExp_sub_variableChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "neronExp_variableChange(第 319 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.neronExp_variableChange") 1 ]

def IsMuAtBadPrimes.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(H が μ_l に対応することの帰結を型にしたもの)",
    sectionId := "genell-lemma-3-5" }

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
