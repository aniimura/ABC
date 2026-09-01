/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GaloisRep.VeluSumKDenomFree
import ABC3.Meta.Claim

/-!
# 第 1135 ブロック —— **`E_q/μ_l` の組み立ての `K` 水準版**（`Found`、節点 5 の配管）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★これは何か

`veluQuotientFull_tate_mu`（第 890）は

* `hu : ∀ i, IsUnit (1 − ζ^i)`（`p ∤ l` から来る）
* Vélu の係数 `v`・`w` が **`R` に取れる**こと

を要求する。★`p ∣ l` ではどちらも成り立たない。

☆本ブロックは同じ主張を **`K` の座標で受ける形**に書き直す:

* `hu` → **`φ(1 − ζ^i) ≠ 0`**（`R → K` の単射性から自動、第 1134 の `zeta_pow_sub_ne_zero_K`）
* `v w : R` → **`v w : K`**（第 1134 が `l⁶·∑ = φ(SV)` を与えるので実際に取れる）
* 座標は `tateXK`・`tateYK`（第 1133 で見たとおり**最初から分母を払った形**）

★これで `IsUnit` が完全に落ちる。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GenEll Finset QuotientGroup
open scoped Classical

/-! ## ★★★★体の上の Vélu の商 -/

section Field

variable {F : Type} [Field F]

/-- ★★★★★★**座標の対の像で取った Vélu の商**（体の上、底変換なし）。 -/
theorem veluQuotientFull_image_eq_field (W : WeierstrassCurve F)
    (s : Finset ℕ) (X Y : ℕ → F)
    (hinj : ∀ i ∈ s, ∀ j ∈ s, ((X i, Y i) : F × F) = ((X j, Y j) : F × F) → i = j)
    (v w : F) (h2 : (2 : F) ≠ 0)
    (hv : v = ∑ i ∈ s, veluV2 W (X i) (Y i))
    (hw : 2 * w = ∑ i ∈ s, (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i))) :
    veluQuotientFull W (s.image (fun i : ℕ => ((X i, Y i) : F × F))) = veluCurve W v w := by
  have h := veluQuotientFull_image_eq (R := F) (K := F) W s X Y
    (by simpa using hinj) v w h2 hv hw
  simpa using h

/-- ★★★★★★★★**点の族が単射で原点を避けるなら Vélu の商は一致する**（体の上）。 -/
theorem veluQuotientFull_points_eq_field (W : WeierstrassCurve F)
    (s : Finset ℕ) (X Y : ℕ → F) (P : ℕ → W.toAffine.Point)
    (hP : ∀ i ∈ s, pointCoords (P i) = ((X i, Y i) : F × F))
    (hPne : ∀ i ∈ s, P i ≠ 0)
    (hPinj : ∀ i ∈ s, ∀ j ∈ s, P i = P j → i = j)
    (v w : F) (h2 : (2 : F) ≠ 0)
    (hv : v = ∑ i ∈ s, veluV2 W (X i) (Y i))
    (hw : 2 * w = ∑ i ∈ s, (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i))) :
    veluQuotientFull W (s.image (fun i : ℕ => pointCoords (P i))) = veluCurve W v w := by
  have himg : s.image (fun i : ℕ => pointCoords (P i))
      = s.image (fun i : ℕ => ((X i, Y i) : F × F)) :=
    Finset.image_congr (fun i hi => hP i hi)
  rw [himg]
  refine veluQuotientFull_image_eq_field W s X Y (fun i hi j hj hij => ?_) v w h2 hv hw
  exact hPinj i hi j hj
    (pointCoords_injective (hPne i hi) (hPne j hj) (by rw [hP i hi, hP j hj]; exact hij))

end Field

/-! ## ★★★★`ha` を要らない座標 -/

section Coord

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★**Tate の点の座標は `tateXK`・`tateYK` そのもの**——★`IsUnit (1 − a)` 不要。

☆`tatePtPair` は `Point.some` なので、座標を読むだけである。
★`pointCoords_tatePtPair`（第 886）はさらに `tateXK_eq` で `R` の元の像に直すが、
そこに `IsUnit` が要る。**その一歩手前で止めれば `IsUnit` は要らない**。 -/
theorem pointCoords_tatePtPair_K (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hw : IsUnit (1 - w)) (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    pointCoords (tatePtPair a w q hq haw hw hne hΔ)
      = ((tateXK (I := I) a w q hq, tateYK (I := I) a w q hq) : K × K) := by
  rw [tatePtPair, pointCoords_some]

end Coord

/-! ## ★★★★★★★★★★★★`E_q/μ_l` の `K` 水準の組み立て -/

section Mu

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `⟨Φ(ζ)⟩` による Vélu の商は `veluCurve (E_q ⊗ K) v w`**（第 1135）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★**`IsUnit (1 − ζ^i)` を仮説に置いていない**——要るのは `φ(1 − ζ^i) ≠ 0` だけである。
☆`v`・`w` は `R` ではなく `K` に取る。第 1134 が
`l⁶·∑ veluV2 = φ(SV)`・`l⁸·∑(…) = φ(SW)` を与えるので実際に取れる。 -/
theorem veluQuotientFull_tate_mu_K (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : 0 < l) (ζ : R) (uζ : Kˣ)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hne : ∀ i ∈ (range l).erase 0, algebraMap R K (1 - ζ ^ i) ≠ 0)
    (v w : K) (h2 : (2 : K) ≠ 0)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU ((tateCurveAt S.q S.hq).map (algebraMap R K))
          (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
          (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        + 2 * (veluV2 ((tateCurveAt S.q S.hq).map (algebraMap R K))
                (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
                (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
              * tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq))) :
    veluQuotientFull ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • tatePhi S hΔ (QuotientGroup.mk uζ))))
      = veluCurve ((tateCurveAt S.q S.hq).map (algebraMap R K)) v w := by
  have hζlR : ζ ^ l = 1 := pow_eq_one_of_map S ζ uζ hζu hζl
  have hmap : ∀ i : ℕ, algebraMap R K (ζ ^ i) = ((uζ ^ i : Kˣ) : K) := by
    intro i
    rw [map_pow, hζu, Units.val_pow_eq_pow_val]
  have hpowl : ∀ i : ℕ, (uζ ^ i) ^ l = 1 := by
    intro i
    rw [← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have hawR : ∀ i : ℕ, (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1)) = S.q := by
    intro i
    have h1 : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      have hl' : l - 1 + 1 = l := Nat.succ_pred_eq_of_pos hl
      calc (ζ ^ i) * (ζ ^ i) ^ (l - 1) = (ζ ^ i) ^ (l - 1 + 1) := by ring
        _ = (ζ ^ i) ^ l := by rw [hl']
        _ = (ζ ^ l) ^ i := by rw [← pow_mul, mul_comm, pow_mul]
        _ = 1 := by rw [hζlR, one_pow]
    calc (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1))
        = S.q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = S.q := by rw [h1, mul_one]
  have hwuR : ∀ i : ℕ, IsUnit (1 - S.q * (ζ ^ i) ^ (l - 1)) := fun i =>
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ S.hq)
  have hmk0 : (QuotientGroup.mk (uζ ^ 0) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
    rw [pow_zero]
    rfl
  have hcne : ∀ i ∈ (range l).erase 0,
      (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 := by
    intro i hi hc
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    exact hi0 (mk_pow_injOn S.v S.Q S.hQ hl uζ hζl hord hil hl (by rw [hc, hmk0]))
  refine veluQuotientFull_points_eq_field
    ((tateCurveAt S.q S.hq).map (algebraMap R K)) ((range l).erase 0)
    (fun i => tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
    (fun i => tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
    (fun k => k • tatePhi S hΔ (QuotientGroup.mk uζ)) ?_ ?_ ?_ v w h2 hv hw
  · intro i hi
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ i,
      tatePhi_of_pow_eq_one S hΔ hl (ζ ^ i) (uζ ^ i) (hmap i) (hpowl i) (hcne i hi)
        (hawR i) (hwuR i) (hne i hi),
      pointCoords_tatePtPair_K _ _ _ S.hq _ _ _ hΔ]
  · intro i hi
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ i]
    exact tatePhi_ne_zero S hΔ (hcne i hi)
  · intro i hi j hj hij
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ i, nsmul_tatePhi S hΔ Φ hΦ uζ j, ← hΦ, ← hΦ] at hij
    have hcl : (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q)
        = QuotientGroup.mk (uζ ^ j) := Φ.injective hij
    exact mk_pow_injOn S.v S.Q S.hQ hl uζ hζl hord
      (Finset.mem_range.1 (Finset.mem_erase.1 hi).2)
      (Finset.mem_range.1 (Finset.mem_erase.1 hj).2) hcl

/-- ★★★★★★★★★★★★**商の楕円性の `K` 水準版**（第 1135）。

☆`veluQuotientFull_tate_mu_K` で書き換えるだけである。
★これが節点 5 の消費者 C（`isElliptic_veluQuotient_tate_mu`）の `hlu` なし版である。 -/
theorem isElliptic_veluQuotientFull_tate_mu_K (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : 0 < l) (ζ : R) (uζ : Kˣ)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hne : ∀ i ∈ (range l).erase 0, algebraMap R K (1 - ζ ^ i) ≠ 0)
    (v w : K) (h2 : (2 : K) ≠ 0)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU ((tateCurveAt S.q S.hq).map (algebraMap R K))
          (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
          (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        + 2 * (veluV2 ((tateCurveAt S.q S.hq).map (algebraMap R K))
                (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
                (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
              * tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)))
    (hell : (veluCurve ((tateCurveAt S.q S.hq).map (algebraMap R K)) v w).IsElliptic) :
    (veluQuotientFull ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • tatePhi S hΔ (QuotientGroup.mk uζ))))).IsElliptic := by
  rw [veluQuotientFull_tate_mu_K S hΔ Φ hΦ hl ζ uζ hζu hζl hord hne v w h2 hv hw]
  exact hell

end Mu

/-! ## ★出典の紐付け(`.src`) -/

def veluQuotientFull_tate_mu_K.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(E_q/μ_l の組み立ての K 水準版。★IsUnit (1−ζ^i) を取らない)",
    sectionId := "genell-lemma-3-5" }

def pointCoords_tatePtPair_K.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate の点の座標は tateXK・tateYK そのもの。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_tate_mu_K.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_image_eq(座標の対の像での商、第 885、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_image_eq") 1,
    .citation "[ABC3]" "tatePhi_of_pow_eq_one(Φ(ζ^i) = tatePtPair、第 884、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhi_of_pow_eq_one") 1,
    .citation "[ABC3]" "zeta_pow_sub_ne_zero_K(φ(1−ζ^i) ≠ 0、第 1134、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.zeta_pow_sub_ne_zero_K") 1,
    .implicitStep
      ("★★**2026-09-01（第 1135）**——第 890 の `hu`（`IsUnit (1 − ζ^i)`）は " ++
       "`pointCoords_tatePtPair` が `tateXK_eq` で `R` の元の像に直すために要っていた。" ++
       "☆**その一歩手前で止めれば要らない**——`tatePtPair` は `Point.some` なので" ++
       "座標は `tateXK`・`tateYK` そのものである。" ++
       "★あとは `veluQuotientFull_image_eq` を `R = K` で使えば体の上の形が出る。") 6 ]

end ABC3.Found.GaloisRep
