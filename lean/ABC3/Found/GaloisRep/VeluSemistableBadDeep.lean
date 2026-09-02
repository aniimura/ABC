/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluSemistableBadFree
import ABC3.Found.GaloisRep.TateDeepCoordMem
import ABC3.Found.GaloisRep.TateCoordDescend
import ABC3.Found.GaloisRep.VeluIntegralFromCoords
import ABC3.Meta.Claim

/-!
# 第 1424 ブロック —— **`p ∣ l` でも「半安定か、核が `μ_l` 型か」**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`hlu` を完全に外した悪い素点

第 1416 の `semistableAt_veluQuotient_bad_ram_free` は `hcop` を落としたが、
まだ `hlu : IsUnit (l : R)`（すなわち `p ∤ l`）を使っていた。
☆使いどころは 2 つだけである（第 1419 の測定）:

1. `isIntegral_veluQuotientFull_of_addOrderOf_prime`（第 1074）——商の整性
2. `exists_vw_tate_mu`（第 1003）——`μ_l` 側で `1 − ζ^i` が単元であること

★★**深い核の側では `hlu` を一度も使わない**。しかも座標が `𝔪` に入るので (1) も外せる
（第 1420-1422 の 3 枚）。

☆本ブロックはそれを組み立てて、**`hlu` なし**で

    `SemistableAt p E′`  または  「核が `μ_l` 型（`R` に原始 `l` 乗根がある）」

を出す。★★★`p ∤ l` なら左が第 1416 で無条件に出るので、
**残るのは `p ∣ l` かつ核が `μ_l` 型の場合だけ**である
（`blocked-leaves.json` の `pDivLMuGapNamed2026_09_02`——Kraus か弱近似）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta

open scoped Classical

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では「Vélu の商が半安定」か「核が `μ_l` 型」**——★**`hlu` なし**（第 1424）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1416 から **`hlu : IsUnit (l : R)`（＝`p ∤ l`）が落ちた**形である
——その代わり結論が二者択一になった。

★★★深い核の側は完全に無条件:

* 商の整性——第 1420（座標が整なら商も整）＋第 1421（座標を `p` へ降ろす）
  ＋第 1422（深い核では座標が `𝔪`）
* `c₄(veluCurve) = c₄(E_q) + 240v` が単元——第 1412（`v ∈ 𝔪`）
* Vélu の `w`——第 1415（対分けは核の形に依らない）

☆`μ_l` 側は右の選択肢に落とす。`p ∤ l` ならそちらも第 1416 で閉じている。 -/
theorem semistableAt_veluQuotient_bad_deepOrMu {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [WeierstrassCurve.IsIntegral (primeSubring p) E]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hΔL : E'.Δ ≠ 0) :
    SemistableAt p E'
      ∨ ∃ (ζ : R) (uζ : Lvˣ), IsPrimitiveRoot ζ l
          ∧ algebraMap R Lv ζ = (uζ : Lv) ∧ uζ ^ l = 1
          ∧ ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1 := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE, hcoords⟩ :=
    exists_variableChange_veluQuotient_tateModel_coords E E' h hl hQ h2K hE'
  have hu0 := vAdd_tateModel_u_eq_zero (E.baseChange Lv) h C₀ hCE hq0
  rcases ABC3.Found.GenEll.primitiveRoot_or_deep_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl P hP hP0 with
    ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ | ⟨y, hPz, hdeep⟩
  · -- ☆`μ_l` 型——右の選択肢
    exact Or.inr ⟨ζ, uζ, hζ, hζu, hζl, hord⟩
  · -- ★★★深い核——無条件に半安定
    left
    obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
    -- ★段 1: 整性（第 1420-1422）
    have hTmem := pointCoords_tatePhi_mem_of_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      y hdeep
    rw [← hPz] at hTmem
    have hmem := pointCoords_mem_primeSubring_of_image_mem he p hpe C₀ E Q _ hcoords hTmem
    haveI hint : (veluQuotientFull E (((range (2 * m + 1)).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))).IsIntegral (primeSubring p) :=
      isIntegral_veluQuotientFull_of_pointCoords_mem p E hl Q hQ hmem
    haveI hintE' : WeierstrassCurve.IsIntegral (primeSubring p) E' := by rw [hE']; exact hint
    -- ★段 2: Vélu の `v`・`w`（第 1412・1415）
    have hyl : (2 * m + 1) • tatePhi
        (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
        (QuotientGroup.mk y) = 0 := by rw [← hPz]; exact hP
    obtain ⟨w, hw⟩ := exists_veluW_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      m y hyl hdeep
    have hquot := veluQuotientFull_tate_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      y hdeep _ w h2K rfl hw
    rw [hPz] at hcurve
    -- ★段 3: `c₄` が単元（第 1412）
    exact semistableAt_velu_of_veluCurve_eq_ram he p hpe E' hΔL
      (tateParamR (E.baseChange Lv) h) hq _ w
      (isUnit_c4_add_240_deep (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0)
        y hdeep _ rfl)
      C₀ (hcurve.symm.trans hquot) hu0

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuotient_bad_deepOrMu.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では商が半安定か核が μ_l 型か。★p ∤ l 不要)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotient_bad_deepOrMu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isIntegral_veluQuotientFull_of_pointCoords_mem(第 1420、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isIntegral_veluQuotientFull_of_pointCoords_mem") 1,
    .citation "[ABC3]" "pointCoords_mem_primeSubring_of_image_mem(第 1421、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.pointCoords_mem_primeSubring_of_image_mem") 1,
    .citation "[ABC3]" "pointCoords_tatePhi_mem_of_deep(第 1422、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.pointCoords_tatePhi_mem_of_deep") 1,
    .citation "[ABC3]" "isUnit_c4_add_240_deep(第 1412、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_c4_add_240_deep") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1424）**——第 1416 から `hlu`（`p ∤ l`）が落ちた。" ++
       "☆深い核の側は完全に無条件である（整性は第 1420-1422、" ++
       "`c₄` の単元性は第 1412、`w` は第 1415）。" ++
       "★★★残るのは **`p ∣ l` かつ核が `μ_l` 型**の場合だけであり、" ++
       "そこでは `E′ ⊗ L_v ≅ E_{q^l}` が言える一方で " ++
       "`L` の上の変数変換で `v_p(c₄) = 0` の整モデルを作るのに" ++
       "Kraus か弱近似が要る（`pDivLMuGapNamed2026_09_02`）。") 17 ]

end ABC3.Found.GaloisRep
