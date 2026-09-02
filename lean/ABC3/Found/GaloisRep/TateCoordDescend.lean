/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateModelPoint
import ABC3.Found.GaloisRep.RamifiedLocalData
import ABC3.Found.GaloisRep.VeluIntegralFromCoords
import ABC3.Meta.Claim

/-!
# 第 1421 ブロック —— **Tate モデルの座標を `p` の側へ降ろす**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1420 の仮説 `hmem` を作る道具

第 1420 は「核の座標が `primeSubring p` に入れば Vélu の商も整」を与えた。
★本ブロックはその仮説を **Tate モデルの側から**作るための 4 本を並べる。

| 定理 | 内容 |
|---|---|
| `mem_primeSubring_of_exists_preimage_ram` | ★★`φ x` が `R` から来るなら `x ∈ primeSubring p` |
| `vc_preimage` | ★★変数変換の逆——`vcX`・`vcY` が `R` なら元も `R` |
| `exists_variableChange_veluQuotient_tateModel_coords` | ★★★★第 1341 に**座標の対応**を足した形 |
| `pointCoords_mem_primeSubring_of_image_mem` | ★★★★★★★★**核の座標が `p` で整** |

☆最後の 1 本が第 1420 の `hmem` そのものである。

★★★これが効くのは **`p ∣ l` の深い核**（第 1410-1417）である
——そこでは Tate の座標が `𝔪 ⊂ R` に入るので、`hlu`（`p ∤ l`）なしで整性が出る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta

open scoped Classical

/-! ## ★★付値の橋（`R` から来るなら `p` で整） -/

section Bridge

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★★**`φ x` が `R` から来るなら `x ∈ primeSubring p`**——分岐版（第 1421）。

☆`exists_preimage_of_mem_primeSubring_ram`（第 1375）の**逆向き**である。
★`v_{Lv}(φ x) ≤ 1` と `hpe` から `(v_p x)^e ≤ 1`、`e ≥ 1` なので `v_p x ≤ 1`。 -/
theorem mem_primeSubring_of_exists_preimage_ram {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (x : L) (hx : ∃ r : R, algebraMap R Lv r = algebraMap L Lv x) :
    x ∈ primeSubring p := by
  obtain ⟨r, hr⟩ := hx
  have h1 : (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
      (algebraMap L Lv x) ≤ 1 := by
    rw [← hr]
    exact HeightOneSpectrum.valuation_le_one _ r
  rw [hpe x] at h1
  refine (mem_primeSubring_iff p x).2 ?_
  by_contra hgt
  rw [not_le] at hgt
  exact absurd h1 (not_le.2 (one_lt_pow₀ hgt (by omega : e ≠ 0)))

end Bridge

/-! ## ★★変数変換の逆 -/

/-- ★★**`vcX`・`vcY` の値が `R` から来るなら、もとの座標も `R` から来る**（第 1421）。

☆`x = u²·vcX + r`、`y = u³·vcY + s·(x − r) + t` で、`u, r, s, t` はすべて `R` の元である。 -/
theorem vc_preimage {R Lv : Type} [CommRing R] [Field Lv] [Algebra R Lv]
    (C : WeierstrassCurve.VariableChange R) (x y : Lv) (a b : R)
    (ha : algebraMap R Lv a = ABC3.Found.GenEll.vcX (C.map (algebraMap R Lv)) x)
    (hb : algebraMap R Lv b = ABC3.Found.GenEll.vcY (C.map (algebraMap R Lv)) x y) :
    x = algebraMap R Lv ((C.u : R) ^ 2 * a + C.r)
      ∧ y = algebraMap R Lv ((C.u : R) ^ 3 * b + C.s * ((C.u : R) ^ 2 * a) + C.t) := by
  have hu : ((C.map (algebraMap R Lv)).u : Lv)
      * ((((C.map (algebraMap R Lv)).u)⁻¹ : Lvˣ) : Lv) = 1 := by
    rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
  have hur : algebraMap R Lv (C.u : R) = ((C.map (algebraMap R Lv)).u : Lv) := rfl
  have hrr : algebraMap R Lv C.r = (C.map (algebraMap R Lv)).r := rfl
  have hss : algebraMap R Lv C.s = (C.map (algebraMap R Lv)).s := rfl
  have htt : algebraMap R Lv C.t = (C.map (algebraMap R Lv)).t := rfl
  rw [ABC3.Found.GenEll.vcX] at ha
  rw [ABC3.Found.GenEll.vcY] at hb
  have hx : x = algebraMap R Lv ((C.u : R) ^ 2 * a + C.r) := by
    rw [map_add, map_mul, map_pow, hur, ha, hrr]
    have h2 : ((C.map (algebraMap R Lv)).u : Lv) ^ 2
        * (((((C.map (algebraMap R Lv)).u)⁻¹ : Lvˣ) : Lv) ^ 2
            * (x - (C.map (algebraMap R Lv)).r))
        = (x - (C.map (algebraMap R Lv)).r) * (((C.map (algebraMap R Lv)).u : Lv)
            * ((((C.map (algebraMap R Lv)).u)⁻¹ : Lvˣ) : Lv)) ^ 2 := by ring
    rw [h2, hu, one_pow, mul_one]
    ring
  refine ⟨hx, ?_⟩
  have hxr : x - (C.map (algebraMap R Lv)).r
      = algebraMap R Lv ((C.u : R) ^ 2 * a) := by
    rw [hx, map_add, hrr]
    ring
  rw [map_add, map_add, map_mul, map_pow, map_mul, hur, hb, hss, htt, ← hxr]
  have h3 : ((C.map (algebraMap R Lv)).u : Lv) ^ 3
      * (((((C.map (algebraMap R Lv)).u)⁻¹ : Lvˣ) : Lv) ^ 3
          * (y - (C.map (algebraMap R Lv)).s * (x - (C.map (algebraMap R Lv)).r)
              - (C.map (algebraMap R Lv)).t))
      = (y - (C.map (algebraMap R Lv)).s * (x - (C.map (algebraMap R Lv)).r)
          - (C.map (algebraMap R Lv)).t)
        * (((C.map (algebraMap R Lv)).u : Lv)
            * ((((C.map (algebraMap R Lv)).u)⁻¹ : Lvˣ) : Lv)) ^ 3 := by ring
  rw [h3, hu, one_pow, mul_one]
  ring

/-! ## ★★★★Tate モデルへの運び方に座標の対応を足す -/

set_option maxHeartbeats 1600000 in
/-- ★★★★**Vélu の商を Tate モデルへ運ぶ——★座標の対応つき**（第 1421）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1341 の `exists_variableChange_veluQuotient_tateModel` に
**核の座標集合の対応**（`vcX`・`vcY` と `algebraMap L Lv` の合成）を足しただけである。
★これが `pointCoords_mem_primeSubring_of_image_mem` の入力になる。 -/
theorem exists_variableChange_veluQuotient_tateModel_coords {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (h2K : (2 : Lv) ≠ 0)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    ∃ (C₀ : WeierstrassCurve.VariableChange R)
      (P : ((tateCurveAt (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv)).toAffine.Point),
      l • P = 0 ∧ P ≠ 0 ∧
      veluQuotientFull
          ((tateCurveAt (tateParamR (E.baseChange Lv) h)
            (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv))
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
        = (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      ∧ C₀ • WeierstrassCurve.integralModel R (E.baseChange Lv)
        = tateCurveAt (tateParamR (E.baseChange Lv) h)
          (tateParamR_mem (E.baseChange Lv) h)
      ∧ (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
        = (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))).image
            (fun z : L × L =>
              ((ABC3.Found.GenEll.vcX (C₀.map (algebraMap R Lv)) (algebraMap L Lv z.1),
                ABC3.Found.GenEll.vcY (C₀.map (algebraMap R Lv)) (algebraMap L Lv z.1)
                  (algebraMap L Lv z.2)) : Lv × Lv)) := by
  set φL : L →+* Lv := algebraMap L Lv with hφL
  set φR : R →+* Lv := algebraMap R Lv with hφR
  obtain ⟨hq, C₀, hne, hCE⟩ := tateParamR_spec (E.baseChange Lv) h
  have hbase : (tateCurveAt (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)).map φR
      = (C₀.map φR) • (E.map φL) :=
    tateModel_baseChange (E.baseChange Lv) h hCE
  haveI hell1 : (E.map φL).IsElliptic := inferInstanceAs ((E.baseChange Lv).IsElliptic)
  haveI hell2 : ((C₀.map φR) • (E.map φL)).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact (((C₀.map φR).u⁻¹).isUnit.pow 12).mul (E.map φL).isUnit_Δ
  haveI hellT : ((tateCurveAt (tateParamR (E.baseChange Lv) h)
      (tateParamR_mem (E.baseChange Lv) h)).map φR).IsElliptic := by
    rw [hbase]; exact hell2
  have hQ1 : addOrderOf (rhPoint φL E Q) = l := by rw [addOrderOf_rhPoint φL E Q, hQ]
  have hQ2 : addOrderOf (ABC3.Found.GenEll.vcPoint (C₀.map φR) (E.map φL)
      (rhPoint φL E Q)) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint (C₀.map φR) (E.map φL) (rhPoint φL E Q), hQ1]
  obtain ⟨P, hP, hP0, himg⟩ := exists_point_image_eq hbase hl _ hQ2
  refine ⟨C₀, P, hP, hP0, ?_, hCE, ?_⟩
  · have hE'map : E'.map φL = veluQuotientFull (E.map φL)
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • rhPoint φL E Q))) := by
      rw [hE', veluQuotientFull_map, image_pointCoords_rhPoint_nsmul φL E hQ]
    have h969 := veluQuotientFull_vcPoint_eq (C₀.map φR) (E.map φL)
      (E'.map φL) hQ1 h2K hE'map
    rw [himg, hbase]
    exact h969.symm
  · rw [himg, image_pointCoords_vcPoint_nsmul (C₀.map φR) (E.map φL) hQ1,
      image_pointCoords_rhPoint_nsmul φL E hQ, Finset.image_image]
    rfl

/-! ## ★★★★★★★★核の座標が `p` で整 -/

/-- ★★★★★★★★**Tate の側の座標が `R` に入れば、核の座標は `p` で整**（第 1421）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆これが第 1420 の `isIntegral_veluQuotientFull_of_pointCoords_mem` の仮説 `hmem` である。
★★★**`p ∤ l` を使わない**——深い核（第 1410-1417）では Tate の座標が `𝔪` に入るので
そのまま適用できる。 -/
theorem pointCoords_mem_primeSubring_of_image_mem {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv]
    {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (C₀ : WeierstrassCurve.VariableChange R) {l : ℕ}
    (E : WeierstrassCurve L) [E.IsElliptic] (Q : E.toAffine.Point)
    (T : Finset (Lv × Lv))
    (himg : T = (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))).image
      (fun z : L × L =>
        ((ABC3.Found.GenEll.vcX (C₀.map (algebraMap R Lv)) (algebraMap L Lv z.1),
          ABC3.Found.GenEll.vcY (C₀.map (algebraMap R Lv)) (algebraMap L Lv z.1)
            (algebraMap L Lv z.2)) : Lv × Lv)))
    (hTmem : ∀ z ∈ T, (∃ a : R, algebraMap R Lv a = z.1)
      ∧ (∃ b : R, algebraMap R Lv b = z.2)) :
    ∀ k ∈ (range l).erase 0,
      (pointCoords (k • Q)).1 ∈ primeSubring p
        ∧ (pointCoords (k • Q)).2 ∈ primeSubring p := by
  intro k hk
  set x := (pointCoords (k • Q)).1 with hxdef
  set y := (pointCoords (k • Q)).2 with hydef
  have hz : ((ABC3.Found.GenEll.vcX (C₀.map (algebraMap R Lv)) (algebraMap L Lv x),
      ABC3.Found.GenEll.vcY (C₀.map (algebraMap R Lv)) (algebraMap L Lv x)
        (algebraMap L Lv y)) : Lv × Lv) ∈ T := by
    rw [himg]
    exact Finset.mem_image.2 ⟨(x, y), Finset.mem_image.2 ⟨k, hk, rfl⟩, rfl⟩
  obtain ⟨⟨a, ha⟩, ⟨b, hb⟩⟩ := hTmem _ hz
  obtain ⟨hx, hy⟩ := vc_preimage C₀ (algebraMap L Lv x) (algebraMap L Lv y) a b ha hb
  exact ⟨mem_primeSubring_of_exists_preimage_ram he p hpe x ⟨_, hx.symm⟩,
    mem_primeSubring_of_exists_preimage_ram he p hpe y ⟨_, hy.symm⟩⟩

/-! ## ★出典の紐付け(`.src`) -/

def mem_primeSubring_of_exists_preimage_ram.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(φ x が R から来るなら x は p で整——分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vc_preimage.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(vcX・vcY が R なら元の座標も R。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_variableChange_veluQuotient_tateModel_coords.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商を Tate モデルへ運ぶ——座標の対応つき。★無条件)",
    sectionId := "genell-lemma-3-5" }

def pointCoords_mem_primeSubring_of_image_mem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate の座標が R なら核の座標は p で整。★p ∤ l 不要)",
    sectionId := "genell-lemma-3-5" }

def pointCoords_mem_primeSubring_of_image_mem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_variableChange_veluQuotient_tateModel(第 1341、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_variableChange_veluQuotient_tateModel") 1,
    .citation "[ABC3]" "image_pointCoords_vcPoint_nsmul(第 915 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.image_pointCoords_vcPoint_nsmul") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1421）**——第 1420 の仮説 `hmem` を" ++
       "**Tate モデルの側から**作る道具を並べた。" ++
       "☆核の座標 `x` は `x = u²·vcX + r`（`u, r ∈ R`）で復元でき、" ++
       "`vcX` が `R` に入れば `φ x` も `R` から来るので `v_p(x) ≥ 0` である。" ++
       "★★★これは `p ∤ l` を一切使わない——`p ∣ l` の深い核で効く。") 17 ]

end ABC3.Found.GaloisRep
