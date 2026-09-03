/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.TateParamJ
import ABC3.Found.GaloisRep.TateDXNeZero
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GaloisRep.TateSetupDvr
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GenEll.MuPrimitiveRoot
import ABC3.Found.GenEll.CyclotomicUnits
import ABC3.Found.GaloisRep.TateModelPoint
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.CompletionValuationBridge
import ABC3.Found.GaloisRep.UnramQuad
import ABC3.Found.GaloisRep.TateMuInvolution
import ABC3.Found.GaloisRep.VeluTateDelta
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import ABC3.Found.GaloisRep.VeluMuSum
import ABC3.Found.GenEll.JScale
import ABC3.Meta.Claim
import ABC3.Skeleton.GenEll.TateODE
import ABC3.Skeleton.GenEll.TateIsogeny.GlobalVelu

/-!
# TateIsogeny —— 極小化した対に当てる・`hvw` を作る・内側で作る

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1001-1002 —— 極小化した対に当てる

★第 973／976 が与えるのは `E` ではなく **`C • E`** についての極小性・乗法還元である
（`C` は `p` で極小にする変数変換、第 954 が `SemistableAt` から取り出す）。
☆したがって第 1000 は `C • E` に当てることになる。

★そのとき Vélu の商も `C • E′` に移す必要があるが、
第 969（`veluQuotientFull_vcPoint_eq`）が**曲線の等式として**それを与える:

    `C • E′ = veluQuotientFull (C • E) (vcPoint C E Q の生成する集合)`

☆`Δ_min` も `jExp` も変数変換で不変（`minDeltaExp_variableChange`・`jExp_variableChange`）、
`SemistableAt` も不変（第 1001）、点の位数も不変（`addOrderOf_vcPoint`）。
★★したがって**そのまま `E`・`E′` の主張に戻る**。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
極小化した対に当てる形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1002）**——第 1000 を `C • E` に当て、
第 969 で Vélu の商を `C • E′` に移し、変数変換の不変性で `E`・`E′` に戻す。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_vc {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (C : WeierstrassCurve.VariableChange L)
    [((C • E).baseChange (p.adicCompletion L)).IsElliptic]
    [((C • E').baseChange (p.adicCompletion L)).IsElliptic]
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L)))
    (h : ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
            (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hCE' : (C • E').IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E'.isUnit_Δ
  have h2L : (2 : L) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (ABC3.Found.GenEll.vcPoint C E Q) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint C E Q]; exact hQ
  have hEq := ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq C E E' hQ h2L hE'
  have hjC : jExp p (C • E) < 0 := by rw [jExp_variableChange p E C]; exact hjneg
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • E)) := by
    rw [jExp_variableChange p E C]; exact hcop
  have hkey := minDeltaExp_eq_mul_at_bad_prime p (C • E) (C • E')
    (semistableAt_variableChange p E C hssE) (semistableAt_variableChange p E' C hssE')
    hjC hl hodd hcopC hmin h hlu hvw hQ' hEq
  rwa [minDeltaExp_variableChange p E' C, minDeltaExp_variableChange p E C] at hkey

def minDeltaExp_eq_mul_at_bad_prime_vc.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——極小化した対に当てる形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_vc.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime(第 1000、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime") 1,
    .citation "[ABC3]" "veluQuotientFull_vcPoint_eq(Vélu の商を変数変換で移す、第 969、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq") 1,
    .citation "[ABC3]" "semistableAt_variableChange(半安定性の不変性、第 1001、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_variableChange") 1 ]

/-! ## ★★★★★★★★★★★★★★★★第 1003 —— `hvw` を作る

★第 1000／1002 の `hvw` は「原始 `l` 乗根 `ζ` ごとに Vélu の係数 `v`・`w` が `R` に取れ、
その `veluCurve` が楕円である」ことを要求する。

☆`v` は和そのものだから自明、`w` は第 961（`exists_veluW_mu`）——
和が **`i ↦ l−i` の対で偶数になる**（第 956-960）ことから来る。
★楕円性は第 962（`isElliptic_veluCurve_tate_map`）と `c₄`・`c₆` の関係式である。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★**[GenEll] `μ_l` に対する Vélu の係数は `R` に取れて、
その商は楕円である**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 1003）**——第 1000／1002 の `hvw` の中身である。 -/
theorem exists_vw_tate_mu {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : R)))
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R) (h2 : (2 : R) ≠ 0)
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (ζ : R) (hζ : IsPrimitiveRoot ζ l) :
    ∃ v w : R,
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt q hq)
              (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      ∧ ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic
      ∧ (tateCurveAt q hq).c₄ + 240 * v = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄
      ∧ (tateCurveAt q hq).c₆ + 504 * v + 6048 * w
        = (l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆ := by
  have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
  -- ★`ζ` を `Kˣ` に上げる
  have hζ0 : ζ ≠ 0 := by
    intro hc
    have := hζ.pow_eq_one
    rw [hc, zero_pow hl.pos.ne'] at this
    exact zero_ne_one this
  have hne : algebraMap R K ζ ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hζ0
  have hζu : algebraMap R K ζ = ((Units.mk0 (algebraMap R K ζ) hne : Kˣ) : K) := rfl
  have hζl : (Units.mk0 (algebraMap R K ζ) hne) ^ l = 1 := by
    ext
    rw [Units.val_pow_eq_pow_val, Units.val_mk0, ← map_pow, hζ.pow_eq_one, map_one,
      Units.val_one]
  have hord : ∀ n : ℕ, 0 < n → n < l →
      (Units.mk0 (algebraMap R K ζ) hne) ^ n ≠ 1 := by
    intro n hn hnl hcon
    have h1 : (algebraMap R K ζ) ^ n = 1 := by
      rw [← Units.val_mk0 hne, ← Units.val_pow_eq_pow_val, hcon, Units.val_one]
    have hz : ζ ^ n = 1 := by
      refine IsFractionRing.injective R K ?_
      rw [map_pow, h1, map_one]
    exact absurd (Nat.le_of_dvd hn ((hζ.pow_eq_one_iff_dvd n).1 hz)) (by omega)
  -- ☆`tateDXpair ≠ 0`
  have hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0 := fun i hi =>
    tateDXpair_ne_zero_of_mu (mkTateSetup q hq hq0) hΔ
      (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl hodd ζ _ hζu hζl hord i
      (Finset.mem_erase.1 hi).1 (Finset.mem_range.1 (Finset.mem_erase.1 hi).2) (hu i hi)
  have h4 := c4_velu_tate hl hζ hlu hu q hq hql h2 hDX
  have h6 := c6_velu_tate hl hζ hu q hq hql h2 hDX
  -- ★`l` は奇素数なので `l = 2m+1`
  obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
  obtain ⟨w, hw0⟩ := ABC3.Found.GaloisRep.exists_veluW_mu (mkTateSetup q hq hq0) hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) m ζ _ hζu hζl hord hu
  -- ☆`(mkTateSetup q hq hq0).q` を `q` に言い直す
  have hw : 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0,
      (veluU (tateCurveAt q hq)
          (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)) := hw0
  have h6' : (tateCurveAt q hq).c₆
      + 504 * (∑ i ∈ (range (2 * m + 1)).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq))
      + 6048 * w
      = ((2 * m + 1 : ℕ) : R) ^ 6 * (tateCurveAt (q ^ (2 * m + 1)) hql).c₆ := by
    rw [← h6]
    linear_combination (3024 : R) * hw
  refine ⟨_, w, rfl, hw, ?_, h4, h6'⟩
  exact ABC3.Found.GaloisRep.isElliptic_veluCurve_tate_map q hq (2 * m + 1) hql _ w hlu
    h4 h6' inferInstance

def exists_vw_tate_mu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l に対する Vélu の係数は R に取れ、商は楕円である)",
    sectionId := "genell-lemma-3-2" }

def exists_vw_tate_mu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_veluW_mu(w が取れる、第 961、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_veluW_mu") 1,
    .citation "[ABC3]" "isElliptic_veluCurve_tate_map(商は楕円、第 962、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isElliptic_veluCurve_tate_map") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1004 —— `hvw` を内側で作る

★第 1003 が `hvw` の中身を与えるので、第 1002 から `hvw` を落とせる。
☆残る局所入力は **`hmin`・`h`・`hlu` の 3 本だけ**である。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
残る入力は極小性・分裂性・`l` の単元性の 3 本**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1004）**——第 1003 で `hvw` を内側から供給する。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_full {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (C : WeierstrassCurve.VariableChange L)
    [((C • E).baseChange (p.adicCompletion L)).IsElliptic]
    [((C • E').baseChange (p.adicCompletion L)).IsElliptic]
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L)))
    (h : ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := hmin
  have hql : (tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) :=
    pow_mem_of_mem_ideal
      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h) hl.pos
  -- ★`E_{q^l} ⊗ Lv` の楕円性（第 1000 と同じ作り方）
  have hqlne : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2
      (pow_ne_zero l (tateParamR_ne_zero ((C • E).baseChange (p.adicCompletion L)) h))
  have hc4T' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (tateCurveAt ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l)
        hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit
        ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql
  have hev' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (evalAdic tateJinvSeries
        ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne
      ((hu'.map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero)
  haveI : ((tateCurveAt ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l)
      hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  exact minDeltaExp_eq_mul_at_bad_prime_vc p E E' hssE hssE' hjneg hl hodd hcop C hmin h hlu
    (fun ζ hζ => by
      obtain ⟨v, w, hv, hw, hell, _, _⟩ := exists_vw_tate_mu _
        (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero ((C • E).baseChange (p.adicCompletion L)) h)
        (tateModel_map_Delta_ne_zero ((C • E).baseChange (p.adicCompletion L)) h)
        hl hodd hlu hql (two_ne_zero_adicCompletionIntegers L p) ζ hζ
      exact ⟨v, w, hv, hw, hell⟩)
    hQ hE'

def minDeltaExp_eq_mul_at_bad_prime_full.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——残る入力は極小性・分裂性・l の単元性)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_full.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_vc(第 1002、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_vc") 1,
    .citation "[ABC3]" "exists_vw_tate_mu(hvw の中身、第 1003、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_vw_tate_mu") 1 ]

end ABC3.Skeleton.GenEll
