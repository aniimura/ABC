/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluBadPrimeAll
import ABC3.Found.GenEll.VeluJExpAssemble
import ABC3.Found.GaloisRep.VeluSemistableBadMu
import ABC3.Meta.Claim

/-!
# 第 1446 ブロック —— **深い核でも乗法還元は `l`-同種で保たれる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か

`Skeleton/GenEll/VeluSemistable.lean` に残っていた `sorry`
（`JExpNegVeluStable`、`Found/GenEll/VeluJExpAssemble.lean` の `def`）を埋める。

☆主張は **`jExp p E < 0` なら `jExp p (E/⟨Q⟩) < 0`**（局所、`l` は奇素数）で、
第 1428（`semistableAt_veluQuotient_bad_ram_all`）と**同じ二者択一**を辿る:

| 核 | `v_p(c₄)` | `v_p(Δ)` | `jExp` |
|---|---|---|---|
| `μ_l` 型 | `4·v(l)` | `12·v(l) + l·v(q)` | `−l·v(q) < 0` |
| 深い核 | `0` | `> 0` | `−v(Δ) < 0` |

## ★★★★★★★★鍵——`u` は消えるので変数変換を追わなくてよい

☆`jExp` は `j = c₄³/Δ` の付値なので、変数変換の `u` は

    `3·(v(c₄) + 4v(u)) − (v(Δ) + 12v(u)) = 3·v(c₄) − v(Δ)`

と**打ち消し合う**。★したがって `jExp_lt_of_local`（本ファイル）は
`vAdd (C₀.u) = 0` を仮定しない——第 1404 の `semistableAt_velu_of_veluCurve_eq_ram` が
`hu` を要求したのは `c₄` だけを見ていたからである。

## ★★★★★★★★深い核で `Δ ∈ 𝔪` を出す道

☆第 1412 の `veluV_deep_mem` は **`v ∈ 𝔪`** を与える。
★同じ議論が `w` にも通る（`exists_veluW_deep_mem`、本ファイル）——
`exists_veluW_of_inv`（第 960）の証人 `w = Σ_{i=1}^{m} c_i` は
`veluU`・`veluV2` の和で、Tate 曲線では `a₃ = 0`・`a₄ ∈ 𝔪` なので
座標が深ければ各項が `𝔪` に入る。☆**`2` の可逆性は要らない**。

★あとは `veluCurve W v w` が `v, w` を `𝔪` に落とすと `W` に戻ること
（`veluCurve_zero`）から、剰余環で `Δ(veluCurve) = Δ(E_q) = 0`、
すなわち `Δ(veluCurve) ∈ 𝔪` である。

## ★逸脱の記録

☆**逸脱なし**。★むしろ仮説が減った——`JExpNegVeluStable` が持っていた
`p ∤ 6`（`v_p(48) = v_p(864) = 0`）と `l ∣ jExp p Y` は
本ファイルの証明では**一度も使わない**ので、
`jExp_neg_veluQuot_badPrime_all` は `l` が奇素数であることだけを要求する。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-! ## ★★★★★★★★配管——`jExp` を局所の付値で読む -/

/-- ☆**`v_p(j) = 3·v_p(c₄) − v_p(Δ)`**——`j = c₄³/Δ` の付値をとるだけ。

☆`valAdd_Delta_eq_neg_jExp`（在庫）の `v_p(c₄) = 0` を外した形である。 -/
theorem jExp_eq_valAdd {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) [W.IsElliptic]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0) :
    jExp p W = 3 * valAdd p (Units.mk0 W.c₄ hc4ne) - valAdd p (Units.mk0 W.Δ hΔ) := by
  have hjeq : W.j = W.Δ⁻¹ * W.c₄ ^ 3 := ABC3.Found.GenEll.j_eq_inv_Delta_mul W
  have hj : W.j ≠ 0 := by
    rw [hjeq]
    exact mul_ne_zero (inv_ne_zero hΔ) (pow_ne_zero 3 hc4ne)
  have hunit : Units.mk0 W.j hj
      = (Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ := by
    ext
    show W.j = ((Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ : Lˣ)
    push_cast
    rw [hjeq]
    show W.Δ⁻¹ * W.c₄ ^ 3 = W.c₄ ^ 3 * (W.Δ)⁻¹
    ring
  have hJ : jExp p W = valAdd p (Units.mk0 W.j hj) := dif_neg hj
  rw [hJ, hunit, valAdd_mul, valAdd_pow, valAdd_inv]
  omega

/-- ★★★★★★★★★★★★
**`e·jExp p E′ = 3·v(c₄ V) − v(Δ V)`**——`V` は `E′ ⊗ Lv` と変数変換で結ばれた曲線。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★変数変換の `u` は `3·4 − 12 = 0` で**打ち消し合う**ので、
`vAdd (C.u) = 0` を仮定しなくてよい。 -/
theorem jExp_mul_eq_local {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E' : WeierstrassCurve L) [E'.IsElliptic]
    (V : WeierstrassCurve Lv) (C : WeierstrassCurve.VariableChange Lv)
    (hEq : C • (E'.baseChange Lv) = V) (hVc4 : V.c₄ ≠ 0) (hVΔ : V.Δ ≠ 0) :
    (e : ℤ) * jExp p E'
      = 3 * vAdd (tateDvrVal R Lv) (Units.mk0 V.c₄ hVc4)
        - vAdd (tateDvrVal R Lv) (Units.mk0 V.Δ hVΔ) := by
  subst hEq
  have hinj : Function.Injective (algebraMap L Lv) := (algebraMap L Lv).injective
  have hΔ : E'.Δ ≠ 0 := E'.isUnit_Δ.ne_zero
  have hc4map : (E'.baseChange Lv).c₄ = algebraMap L Lv E'.c₄ := WeierstrassCurve.map_c₄ _ _
  have hΔmap : (E'.baseChange Lv).Δ = algebraMap L Lv E'.Δ := WeierstrassCurve.map_Δ _ _
  have hΔLv : (E'.baseChange Lv).Δ ≠ 0 := by
    rw [hΔmap]; exact (map_ne_zero_iff _ hinj).2 hΔ
  have hc4Lv : (E'.baseChange Lv).c₄ ≠ 0 := by
    intro h0
    refine hVc4 ?_
    rw [WeierstrassCurve.variableChange_c₄, h0, mul_zero]
  have hc4ne : E'.c₄ ≠ 0 := by
    intro h0; rw [hc4map, h0, map_zero] at hc4Lv; exact hc4Lv rfl
  have hc4Lv' : algebraMap L Lv E'.c₄ ≠ 0 := by rw [← hc4map]; exact hc4Lv
  have hΔLv' : algebraMap L Lv E'.Δ ≠ 0 := by rw [← hΔmap]; exact hΔLv
  have h1 := vAdd_c4_variableChange (R := R) (E'.baseChange Lv) hc4Lv C hVc4
  have h2 := vAdd_Delta_variableChange (R := R) (E'.baseChange Lv) hΔLv C
  have h2' : vAdd (tateDvrVal R Lv) (Units.mk0 ((C • E'.baseChange Lv).Δ) hVΔ)
      = vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).Δ) hΔLv)
        - 12 * vAdd (tateDvrVal R Lv) C.u := h2
  have h3 : vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).c₄) hc4Lv)
      = (e : ℤ) * valAdd p (Units.mk0 E'.c₄ hc4ne) := by
    rw [vAdd_mk0_congr (R := R) hc4Lv hc4Lv' hc4map]
    exact vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe E'.c₄ hc4ne hc4Lv'
  have h4 : vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).Δ) hΔLv)
      = (e : ℤ) * valAdd p (Units.mk0 E'.Δ hΔ) := by
    rw [vAdd_mk0_congr (R := R) hΔLv hΔLv' hΔmap]
    exact vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe E'.Δ hΔ hΔLv'
  rw [h1, h2', h3, h4, jExp_eq_valAdd p E' hΔ hc4ne]
  ring

/-- ★★★★★★★★★★★★★★★★
**`3·v(c₄ V) < v(Δ V)` なら `jExp p E′ < 0`**（分岐版）。 -/
theorem jExp_lt_of_local {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv]
    {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E' : WeierstrassCurve L) [E'.IsElliptic]
    (V : WeierstrassCurve Lv) (C : WeierstrassCurve.VariableChange Lv)
    (hEq : C • (E'.baseChange Lv) = V) (hVc4 : V.c₄ ≠ 0) (hVΔ : V.Δ ≠ 0)
    (hlt : 3 * vAdd (tateDvrVal R Lv) (Units.mk0 V.c₄ hVc4)
      < vAdd (tateDvrVal R Lv) (Units.mk0 V.Δ hVΔ)) :
    jExp p E' < 0 := by
  have hmul := jExp_mul_eq_local (R := R) p hpe E' V C hEq hVc4 hVΔ
  have hepos : (0 : ℤ) < (e : ℤ) := by exact_mod_cast he
  nlinarith [hmul, hlt, hepos]

/-- ☆`x = n^k·z` の付値の帳簿。 -/
theorem vAdd_mk0_pow_mul {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    {n z x : K} (hn : n ≠ 0) (hz : z ≠ 0) (hx : x ≠ 0) (k : ℕ) (hxe : x = n ^ k * z) :
    vAdd (tateDvrVal R K) (Units.mk0 x hx)
      = (k : ℤ) * vAdd (tateDvrVal R K) (Units.mk0 n hn)
        + vAdd (tateDvrVal R K) (Units.mk0 z hz) := by
  have h : Units.mk0 x hx = (Units.mk0 n hn) ^ k * Units.mk0 z hz := by
    ext
    rw [Units.val_mk0, Units.val_mul, Units.val_pow_eq_pow_val, Units.val_mk0, Units.val_mk0]
    exact hxe
  rw [h, vAdd_mul, vAdd_pow]

/-- ☆`c₄ = n⁴c₄′`・`c₆ = n⁶c₆′` なら `Δ = n¹²Δ′`（体の版）。

☆`1728Δ = c₄³ − c₆²` を使うだけで、`Delta_velu_tate_eq`（第 962）の体版である。 -/
theorem Delta_scaled_of_c4_c6 {K : Type} [Field K] [CharZero K]
    (V T : WeierstrassCurve K) (n : K)
    (h4 : V.c₄ = n ^ 4 * T.c₄) (h6 : V.c₆ = n ^ 6 * T.c₆) :
    V.Δ = n ^ 12 * T.Δ := by
  have e1 := WeierstrassCurve.c_relation V
  have e2 := WeierstrassCurve.c_relation T
  have h1728 : (1728 : K) ≠ 0 := by norm_num
  refine mul_left_cancel₀ h1728 ?_
  calc (1728 : K) * V.Δ = V.c₄ ^ 3 - V.c₆ ^ 2 := e1
    _ = (n ^ 4 * T.c₄) ^ 3 - (n ^ 6 * T.c₆) ^ 2 := by rw [h4, h6]
    _ = n ^ 12 * (T.c₄ ^ 3 - T.c₆ ^ 2) := by ring
    _ = n ^ 12 * ((1728 : K) * T.Δ) := by rw [← e2]
    _ = (1728 : K) * (n ^ 12 * T.Δ) := by ring

/-! ## ★★★★★★★★深い核——`v`・`w` はともに `𝔪` に入る -/

/-- ☆**`v, w ∈ J` かつ `Δ(W) ∈ J` なら `Δ(veluCurve W v w) ∈ J`**。

☆剰余環 `A ⧸ J` へ落とすと `veluCurve` の `v`・`w` が `0` になり
`veluCurve_zero` で `W` に戻るからである。★`2` や `1728` の可逆性は要らない。 -/
theorem veluCurve_Delta_mem {A : Type} [CommRing A] (J : Ideal A)
    (W : WeierstrassCurve A) (v w : A) (hv : v ∈ J) (hw : w ∈ J) (hΔ : W.Δ ∈ J) :
    (veluCurve W v w).Δ ∈ J := by
  have hq : (Ideal.Quotient.mk J) ((veluCurve W v w).Δ) = (Ideal.Quotient.mk J) W.Δ := by
    rw [← WeierstrassCurve.map_Δ, ← WeierstrassCurve.map_Δ, veluCurve_map,
      (Ideal.Quotient.eq_zero_iff_mem.2 hv), (Ideal.Quotient.eq_zero_iff_mem.2 hw),
      veluCurve_zero]
  rw [← Ideal.Quotient.eq_zero_iff_mem] at hΔ ⊢
  rw [hq, hΔ]

/-- ★★★★★★★★**Vélu の `w` は `J` の中で取れる**——
`exists_veluW_of_inv`（第 960）に所属を足した形。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆証人 `w = Σ_{i=1}^{m}(u_i + x_i(v_i + v_i′))` の各項が `J` に入ることを見るだけである。
★`a₃ ∈ J`・`a₄ ∈ J` が要る——Tate 曲線では `a₃ = 0`・`a₄ = −5s₃(q) ∈ 𝔪` で成り立つ。 -/
theorem exists_veluW_of_inv_mem {A : Type} [CommRing A] (J : Ideal A)
    (W : WeierstrassCurve A) (m : ℕ) (X Y : ℕ → A)
    (hX : ∀ i ∈ Icc 1 m, X (2 * m + 1 - i) = X i)
    (hY : ∀ i ∈ Icc 1 m, Y (2 * m + 1 - i) = W.toAffine.negY (X i) (Y i))
    (ha3 : W.a₃ ∈ J) (ha4 : W.a₄ ∈ J)
    (hXJ : ∀ i ∈ Icc 1 m, X i ∈ J) (hYJ : ∀ i ∈ Icc 1 m, Y i ∈ J) :
    ∃ w : A, w ∈ J ∧ 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0,
      (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i)) := by
  have hgy : ∀ x y : A, x ∈ J → y ∈ J → veluGy W x y ∈ J := by
    intro x y hx hy
    rw [veluGy]
    exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.mul_mem_left _ _ hy)
      (Ideal.mul_mem_left _ _ hx)) ha3
  have hnegY : ∀ x y : A, x ∈ J → y ∈ J → W.toAffine.negY x y ∈ J := by
    intro x y hx hy
    rw [WeierstrassCurve.Affine.negY]
    exact Ideal.sub_mem _ (Ideal.sub_mem _ (Submodule.neg_mem _ hy)
      (Ideal.mul_mem_left _ _ hx)) ha3
  have hv2 : ∀ x y : A, x ∈ J → y ∈ J → veluV2 W x y ∈ J := by
    intro x y hx hy
    rw [veluV2, veluGx]
    refine Ideal.sub_mem _ (Ideal.add_mem _ (Ideal.add_mem _ ?_ (Ideal.mul_mem_left _ _ hx))
      ha4) (Ideal.mul_mem_left _ _ hy)
    exact Ideal.mul_mem_left _ _ (by rw [sq]; exact Ideal.mul_mem_left _ _ hx)
  refine ⟨∑ i ∈ Icc 1 m, (veluU W (X i) (Y i)
      + X i * (veluV2 W (X i) (Y i) + veluV2 W (X i) (W.toAffine.negY (X i) (Y i)))), ?_, ?_⟩
  · refine Ideal.sum_mem _ (fun i hi => ?_)
    refine Ideal.add_mem _ ?_ (Ideal.mul_mem_left _ _ (Ideal.add_mem _
      (hv2 _ _ (hXJ i hi) (hYJ i hi))
      (hv2 _ _ (hXJ i hi) (hnegY _ _ (hXJ i hi) (hYJ i hi)))))
    rw [veluU, sq]
    exact Ideal.mul_mem_left _ _ (hgy _ _ (hXJ i hi) (hYJ i hi))
  · refine two_mul_sum_eq_of_pair_even m _ _ ?_
    intro i hi
    rw [hX i hi, hY i hi]
    exact veluTerm_pair_even W (X i) (Y i)

/-- ★★★★★★★★★★★★**深い核では Vélu の `w` は `𝔪` の中で取れる**（第 1446）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`exists_veluW_deep`（第 1415）に所属を足した形である。
★これが `veluV_deep_mem`（第 1412、`v ∈ 𝔪`）の `w` 版で、
両方あって初めて `Δ(veluCurve) ∈ 𝔪` が出る。 -/
theorem exists_veluW_deep_mem {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {K : Type} [Field K] [Algebra R K] (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (m : ℕ) (y : Kˣ)
    (hyl : (2 * m + 1) • tatePhi S hΔ (QuotientGroup.mk y) = 0)
    (hdeep : ∀ i : ℕ, 0 < i → i < 2 * m + 1 →
      ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y)) :
    ∃ w : R, w ∈ I ∧ 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0,
      (veluU (tateCurveAt S.q S.hq)
          (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
            (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
          (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
            (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
        + 2 * (veluV2 (tateCurveAt S.q S.hq)
                (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
                  (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
                (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
                  (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
              * tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
                  (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)) := by
  have hmem : ∀ i ∈ Finset.Icc 1 m, i ∈ (range (2 * m + 1)).erase 0 := by
    intro i hi
    rw [Finset.mem_Icc] at hi
    rw [Finset.mem_erase, Finset.mem_range]
    omega
  have hAm : ∀ i ∈ Finset.Icc 1 m, tateAOf S (QuotientGroup.mk (y ^ i)) ∈ I := by
    intro i hi
    have hi' := hmem i hi
    exact tateAOf_pow_mem S y i (hdeep i (Nat.pos_of_ne_zero (Finset.mem_erase.1 hi').1)
      (Finset.mem_range.1 (Finset.mem_erase.1 hi').2))
  refine exists_veluW_of_inv_mem I (tateCurveAt S.q S.hq) m
    (fun i => tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
      (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
    (fun i => tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
      (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq) ?_ ?_ ?_ ?_ ?_ ?_
  · intro i hi
    exact (tateXpair_deep_inv S hΔ Φ hΦ y hyl hdeep i (hmem i hi)).1
  · intro i hi
    exact (tateXpair_deep_inv S hΔ Φ hΦ y hyl hdeep i (hmem i hi)).2
  · rw [tateCurveAt_a₃]; exact Ideal.zero_mem I
  · exact tateCurveAt_a₄_mem S.q S.hq
  · intro i hi
    exact tateXpair_mem _ _ _ S.hq (hAm i hi) (tateWOf_mem S _)
  · intro i hi
    exact tateYpair_mem _ _ _ S.hq (hAm i hi) (tateWOf_mem S _)

/-! ## ★★★★★★★★★★★★深い核の枚 -/

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★
**深い核の Vélu の商は `jExp < 0`**（第 1446）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v, w ∈ 𝔪` から `Δ(veluCurve) ∈ 𝔪`（`veluCurve_Delta_mem`）、
`c₄(veluCurve) = c₄(E_q) + 240v` は単元（第 1412）なので
`3·v(c₄) − v(Δ) = −v(Δ) < 0` である。 -/
theorem jExp_neg_of_veluCurve_deep {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E' : WeierstrassCurve L) [E'.IsElliptic]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (v w : R)
    (hvmem : v ∈ IsLocalRing.maximalIdeal R) (hwmem : w ∈ IsLocalRing.maximalIdeal R)
    (hunit : IsUnit ((tateCurveAt q hq).c₄ + 240 * v))
    (C₀ : WeierstrassCurve.VariableChange R)
    (hEq : (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      = (veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv)) :
    jExp p E' < 0 := by
  have hinj : Function.Injective (algebraMap R Lv) := IsFractionRing.injective R Lv
  have hEΔ : (E'.baseChange Lv).Δ ≠ 0 := by
    have hmap : (E'.baseChange Lv).Δ = algebraMap L Lv E'.Δ := WeierstrassCurve.map_Δ _ _
    rw [hmap]
    exact (map_ne_zero_iff _ (algebraMap L Lv).injective).2 E'.isUnit_Δ.ne_zero
  have hc4unit : IsUnit ((veluCurve (tateCurveAt q hq) v w).c₄) := by
    rw [veluCurve_c₄]; exact hunit
  have hVc4map : ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv)).c₄
      = algebraMap R Lv ((veluCurve (tateCurveAt q hq) v w).c₄) :=
    WeierstrassCurve.map_c₄ _ _
  have hVΔmap : ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv)).Δ
      = algebraMap R Lv ((veluCurve (tateCurveAt q hq) v w).Δ) :=
    WeierstrassCurve.map_Δ _ _
  have hVc4ne : ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv)).c₄ ≠ 0 := by
    rw [hVc4map]; exact (hc4unit.map (algebraMap R Lv)).ne_zero
  have hVΔne : ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv)).Δ ≠ 0 := by
    rw [← hEq]
    exact variableChange_Delta_ne_zero (E'.baseChange Lv) hEΔ (C₀.map (algebraMap R Lv))
  have hΔq : (tateCurveAt q hq).Δ ∈ IsLocalRing.maximalIdeal R := by
    obtain ⟨uu, huu, hΔeq⟩ := tateCurveAt_Delta_eq_mul_unit
      (I := IsLocalRing.maximalIdeal R) q hq
    rw [hΔeq]; exact Ideal.mul_mem_right _ _ hq
  have hΔmem := veluCurve_Delta_mem (IsLocalRing.maximalIdeal R)
    (tateCurveAt q hq) v w hvmem hwmem hΔq
  have hne2 : algebraMap R Lv ((veluCurve (tateCurveAt q hq) v w).c₄) ≠ 0 := by
    rw [← hVc4map]; exact hVc4ne
  have hne3 : algebraMap R Lv ((veluCurve (tateCurveAt q hq) v w).Δ) ≠ 0 := by
    rw [← hVΔmap]; exact hVΔne
  have hzero : vAdd (tateDvrVal R Lv)
      (Units.mk0 (((veluCurve (tateCurveAt q hq) v w).map
        (algebraMap R Lv)).c₄) hVc4ne) = 0 := by
    rw [vAdd_mk0_congr (R := R) hVc4ne hne2 hVc4map]
    exact tateDvrVal_eq_zero_of_isUnit _ hc4unit hne2
  have hpos : 0 < vAdd (tateDvrVal R Lv)
      (Units.mk0 (((veluCurve (tateCurveAt q hq) v w).map
        (algebraMap R Lv)).Δ) hVΔne) := by
    rw [vAdd_mk0_congr (R := R) hVΔne hne3 hVΔmap]
    exact tateDvrVal_pos_of_mem_max _ ⟨_, hΔmem, rfl⟩
  refine jExp_lt_of_local he p hpe E' _ (C₀.map (algebraMap R Lv)) hEq hVc4ne hVΔne ?_
  rw [hzero]
  simpa using hpos

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★
**深い核の枚**——Tate モデルの上の位数 `2m+1` の深い点で割った商は `jExp < 0`。 -/
theorem jExp_neg_veluQuotient_bad_deep {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {m : ℕ} (h2K : (2 : Lv) ≠ 0)
    (C₀ : WeierstrassCurve.VariableChange R)
    (P : ((tateCurveAt (tateParamR (E.baseChange Lv) h)
      (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv)).toAffine.Point)
    (hP : (2 * m + 1) • P = 0)
    (hcurve : veluQuotientFull
        ((tateCurveAt (tateParamR (E.baseChange Lv) h)
          (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv))
        (((range (2 * m + 1)).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv))
    (y : Lvˣ)
    (hPz : P = tatePhi (mkTateSetup (tateParamR (E.baseChange Lv) h)
      (tateParamR_mem (E.baseChange Lv) h)
      (tateParamR_ne_zero (E.baseChange Lv) h))
      (tateModel_map_Delta_ne_zero (E.baseChange Lv) h) (QuotientGroup.mk y))
    (hdeep : ∀ i : ℕ, 0 < i → i < 2 * m + 1 →
      ¬ (vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
            (tateParamR_mem (E.baseChange Lv) h)
            (tateParamR_ne_zero (E.baseChange Lv) h)).v
          (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
            (tateParamR_mem (E.baseChange Lv) h)
            (tateParamR_ne_zero (E.baseChange Lv) h)).Q
        ∣ (i : ℤ) * vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
            (tateParamR_mem (E.baseChange Lv) h)
            (tateParamR_ne_zero (E.baseChange Lv) h)).v y)) :
    jExp p E' < 0 := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hyl : (2 * m + 1) • tatePhi
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (QuotientGroup.mk y) = 0 := by rw [← hPz]; exact hP
  obtain ⟨w, hwmem, hw⟩ := exists_veluW_deep_mem
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
    m y hyl hdeep
  have hquot := veluQuotientFull_tate_deep
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
    y hdeep _ w h2K rfl hw
  rw [hPz] at hcurve
  have hvmem := veluV_deep_mem
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) y hdeep _ rfl
  have hunit := isUnit_c4_add_240_deep
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) y hdeep _ rfl
  exact jExp_neg_of_veluCurve_deep he p hpe E' _ hq _ w hvmem hwmem hunit C₀
    (hcurve.symm.trans hquot)

/-! ## ★★★★★★★★★★★★`μ_l` 型の枚 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★
**`μ_l` 型の核の Vélu の商は `jExp < 0`**（第 1446）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`c₄ = l⁴·c₄(E_{q^l})`・`c₆ = l⁶·c₆(E_{q^l})`（第 1136 系）から
`Δ = l¹²·Δ(E_{q^l})`、したがって

    `3·v(c₄) − v(Δ) = 12·v(l) − (12·v(l) + l·v(q)) = −l·v(q) < 0`

である。★`p ∣ l` でも `v(l)` が打ち消し合うので条件は要らない。 -/
theorem jExp_neg_veluQuotient_bad_mu {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (h2K : (2 : Lv) ≠ 0)
    (C₀ : WeierstrassCurve.VariableChange R)
    (P : ((tateCurveAt (tateParamR (E.baseChange Lv) h)
      (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv)).toAffine.Point)
    (hcurve : veluQuotientFull
        ((tateCurveAt (tateParamR (E.baseChange Lv) h)
          (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv))
    (ζ : R) (uζ : Lvˣ) (hζ : IsPrimitiveRoot ζ l)
    (hζu : algebraMap R Lv ζ = (uζ : Lv)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hPz : P = tatePhi (mkTateSetup (tateParamR (E.baseChange Lv) h)
      (tateParamR_mem (E.baseChange Lv) h)
      (tateParamR_ne_zero (E.baseChange Lv) h))
      (tateModel_map_Delta_ne_zero (E.baseChange Lv) h) (QuotientGroup.mk uζ)) :
    jExp p E' < 0 := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hinj : Function.Injective (algebraMap R Lv) := IsFractionRing.injective R Lv
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    pow_mem_of_mem_ideal hq hl.pos
  have hqne : algebraMap R Lv (tateParamR (E.baseChange Lv) h) ≠ 0 :=
    (map_ne_zero_iff _ hinj).2 hq0
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ hinj).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit _ hql).map (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu'u, hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'u.map (algebraMap R Lv)).ne_zero)
  haveI hellQl : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  have hΔl : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).Δ ≠ 0 :=
    ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).isUnit_Δ.ne_zero
  obtain ⟨v, w, hv, hw, hell⟩ := exists_vw_tate_mu_K hl hζ
    (tateParamR (E.baseChange Lv) h) hq hql hΔl
  haveI := hell
  have hne : ∀ i ∈ (range l).erase 0, algebraMap R Lv (1 - ζ ^ i) ≠ 0 :=
    fun i hi => zeta_pow_sub_ne_zero_K hinj hζ hi
  have hquot := veluQuotientFull_tate_mu_K
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
    hl.pos ζ uζ hζu hζl hord hne v w h2K hv hw
  have hquot' : veluQuotientFull ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
        (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • tatePhi
        (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ (QuotientGroup.mk uζ))))
      = veluCurve ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
        v w := hquot
  have hvS : (l : Lv) ^ 6 * v
      = algebraMap R Lv (veluSV l ζ (tateParamR (E.baseChange Lv) h) hq) := by
    rw [hv]; exact sum_natCast_pow_mul_veluV2_K hinj hζ _ _ hq
  have hwS : (l : Lv) ^ 8 * (2 * w)
      = algebraMap R Lv (veluSW l ζ (tateParamR (E.baseChange Lv) h) hq) := by
    rw [hw]; exact sum_natCast_pow_mul_veluW_K hinj hζ _ _ hq
  obtain ⟨h4, h6⟩ := c4_c6_veluCurve_tate_field hl hζ
    (tateParamR (E.baseChange Lv) h) hq hql v w hvS hwS
  rw [hPz] at hcurve
  have hEq := hcurve.symm.trans hquot'
  have hlLv : ((l : Lv)) ≠ 0 := by
    simpa using (Nat.cast_ne_zero (R := Lv)).2 hl.ne_zero
  have hzc4 : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).c₄ = algebraMap R Lv
        (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ :=
    WeierstrassCurve.map_c₄ _ _
  have hzc4ne : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).c₄ ≠ 0 := by rw [hzc4]; exact hc4T'
  have hz4 : vAdd (tateDvrVal R Lv) (Units.mk0 (((tateCurveAt
      ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).c₄) hzc4ne) = 0 := by
    rw [vAdd_mk0_congr (R := R) hzc4ne hc4T' hzc4]
    exact tateDvrVal_eq_zero_of_isUnit _ (tateCurveAt_c4_isUnit _ hql) hc4T'
  have hVc4ne : (veluCurve ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)) v w).c₄ ≠ 0 := by
    rw [h4]; exact mul_ne_zero (pow_ne_zero 4 hlLv) hzc4ne
  have hVΔne : (veluCurve ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)) v w).Δ ≠ 0 :=
    (veluCurve ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)) v w).isUnit_Δ.ne_zero
  have hΔV := Delta_scaled_of_c4_c6 _ _ ((l : Lv)) h4 h6
  have hc4val := vAdd_mk0_pow_mul (R := R) hlLv hzc4ne hVc4ne 4 h4
  have hΔval := vAdd_mk0_pow_mul (R := R) hlLv hΔl hVΔne 12 hΔV
  have hΔmap : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).Δ = algebraMap R Lv
        (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).Δ :=
    WeierstrassCurve.map_Δ _ _
  have hΔne2 : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).Δ ≠ 0 := by
    rw [← hΔmap]; exact hΔl
  have hTΔ : vAdd (tateDvrVal R Lv) (Units.mk0 (((tateCurveAt
      ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).Δ) hΔl)
      = (l : ℤ) * vAdd (tateDvrVal R Lv)
        (Units.mk0 (algebraMap R Lv (tateParamR (E.baseChange Lv) h)) hqne) := by
    rw [vAdd_mk0_congr (R := R) hΔl hΔne2 hΔmap]
    exact vAdd_Delta_tateCurveAt_pow _ hql hΔne2 hqne
  have hqpos : 0 < vAdd (tateDvrVal R Lv)
      (Units.mk0 (algebraMap R Lv (tateParamR (E.baseChange Lv) h)) hqne) :=
    tateDvrVal_pos_of_mem_max _ ⟨_, hq, rfl⟩
  have hlpos : (0 : ℤ) < (l : ℤ) := by exact_mod_cast hl.pos
  refine jExp_lt_of_local he p hpe E' _ (C₀.map (algebraMap R Lv)) hEq hVc4ne hVΔne ?_
  rw [hc4val, hz4, hΔval, hTΔ]
  push_cast
  nlinarith [hqpos, hlpos]

/-! ## ★★★★★★★★★★★★★★★★二枚を合わせる -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**分裂乗法還元の Vélu の商は `jExp < 0`**（分岐版、第 1446）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1428 の `semistableAt_veluQuotient_bad_ram_all` と**同じ二者択一**
（`primitiveRoot_or_deep_of_torsion_point`）を辿る。
★`p ∤ 6` も核の座標の整性も要らない——`jExp` は `c₄` と `Δ` の比だけで決まるからである。 -/
theorem jExp_neg_veluQuotient_bad_ram_all {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    jExp p E' < 0 := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE⟩ :=
    exists_variableChange_veluQuotient_tateModel E E' h hl hQ h2K hE'
  rcases ABC3.Found.GenEll.primitiveRoot_or_deep_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl P hP hP0 with
    ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ | ⟨y, hPz, hdeep⟩
  · exact jExp_neg_veluQuotient_bad_mu he p hpe E E' h hl h2K C₀ P hcurve
      ζ uζ hζ hζu hζl hord hPz
  · obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
    exact jExp_neg_veluQuotient_bad_deep he p hpe E E' h h2K C₀ P hP hcurve y hPz hdeep

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**乗法還元の Vélu の商は `jExp < 0`**——★**分裂性も `p ∤ l` も仮定しない**（第 1446）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1428 の `semistableAt_veluQuot_multRed_local_all` と同じ形——
非分裂なら不分岐 2 次拡大へ上げるだけである。 -/
theorem jExp_neg_veluQuot_multRed_local_all
    (L : Type) [Field L] [NumberField L] [inst : DecidableEq L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (W : WeierstrassCurve L) [W.IsElliptic] (C : WeierstrassCurve.VariableChange L)
    [(C • W).IsElliptic]
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • W))
    (hc4ne : (C • W).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • W).c₄) hc4ne) = 0)
    [((C • W).baseChange Lv).IsElliptic]
    [WeierstrassCurve.HasMultiplicativeReduction R ((C • W).baseChange Lv)]
    (hj : jExp p W < 0) {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (h2K : (2 : Lv) ≠ 0)
    (Q : (C • W).toAffine.Point) (hQ : addOrderOf Q = l)
    [hVell : (veluQuotientFull (C • W)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic] :
    jExp p (veluQuotientFull (C • W)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) < 0 := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  haveI hCint : WeierstrassCurve.IsIntegral (primeSubring p) (C • W) := by
    haveI := hC; infer_instance
  have hA := integralModel_c4_isUnit (R := R) ((C • W).baseChange Lv)
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction R ((C • W).baseChange Lv)
  · exact jExp_neg_veluQuotient_bad_ram_all he p hpe (C • W) _ hs hl hodd h2K hQ rfl
  · set F := splitQuadPoly (WeierstrassCurve.integralModel R ((C • W).baseChange Lv)) hA
      with hFdef
    have hfm : F.Monic := monic_splitQuadPoly _ hA
    have hfd : F.natDegree = 2 := natDegree_splitQuadPoly _ hA
    have hirr := irreducible_map_residue_of_not_splits hfm hfd
      (fun hsp => hs (hasSplit_of_splits_splitQuadPoly _ hA hsp))
    haveI hdom := isDomain_adjoinRoot hfm hirr
    haveI hlocR := isLocalRing_adjoinRoot hfm hfd hirr
    haveI hdvr := isDiscreteValuationRing_adjoinRoot hfm hfd hirr
    haveI hac := isAdicComplete_maximalIdeal_adjoinRoot hfm hfd hirr
    haveI hcz : CharZero (AdjoinRoot F) :=
      charZero_of_injective_algebraMap (algebraMap_adjoinRoot_injective (R := R) hfm hfd)
    haveI hcz2 : CharZero (FractionRing (AdjoinRoot F)) := by
      refine charZero_of_injective_algebraMap (R := AdjoinRoot F) ?_
      exact IsFractionRing.injective (AdjoinRoot F) (FractionRing (AdjoinRoot F))
    letI : Algebra L (FractionRing (AdjoinRoot F)) :=
      ((quadFieldHom (K := Lv) hfm hfd).comp (algebraMap L Lv)).toAlgebra
    have halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot F)) x
        = quadFieldHom hfm hfd (algebraMap L Lv x) := fun _ => rfl
    haveI hmin' := isMinimal_baseChange_ext_ram p hpe hfm hfd hirr halg W C hC hc4ne hc4
    haveI hmult' :=
      hasMultiplicativeReduction_ext_ram he p hpe hfm hfd hirr halg W C hC hc4ne hc4 hj
    have h' := hasSplitMultiplicativeReduction_ext (C • W) hA hFdef halg
    haveI hell' : ((C • W).baseChange (FractionRing (AdjoinRoot F))).IsElliptic :=
      isElliptic_baseChange (C • W)
    have h2K' : (2 : FractionRing (AdjoinRoot F)) ≠ 0 := two_ne_zero
    exact jExp_neg_veluQuotient_bad_ram_all
      (Lv := FractionRing (AdjoinRoot F)) (R := AdjoinRoot F) he p
      (valuation_algebraMap_ext_ram p hpe hfm hfd hirr halg) (C • W) _ h' hl hodd h2K' hQ rfl

/-! ## ★★★★★★★★★★★★★★★★★★★★大域の形 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商も `jExp < 0`**——★**無条件**（第 1446）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1436 の `semistableAt_veluQuot_badPrime_all` の `jExp` 版である。
★仮説は `SemistableAt p E`・`jExp p E < 0`・`l` が奇素数・`addOrderOf Q = l` だけ
——**`p ∤ 6` も `l ∤ v_p(j)` も要らない**。 -/
theorem jExp_neg_veluQuot_badPrime_all {L : Type} [Field L] [NumberField L]
    [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E) (hj : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    [hVell : (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic] :
    jExp p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) < 0 := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  obtain ⟨C, hC, hc4ne, hc4⟩ := exists_minimal_c4_unit_of_jExp_neg p E hss hj
  obtain ⟨e, he1, hle, hIe⟩ := exists_locCyc_package p hl
  have hpe := valuation_algebraMap_locCycField p hl he1 hIe
  haveI hCell : (C • E).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (variableChange_Delta_ne_zero E hΔ C)⟩
  haveI hbcell : ((C • E).baseChange (locCycField p hl)).IsElliptic :=
    isElliptic_baseChange (C • E)
  haveI hmin := isMinimal_baseChange_at_bad_prime_ram (Lv := locCycField p hl)
    (R := locCycRing p hl) p hpe E C hC hc4ne hc4
  haveI hmult := hasMultiplicativeReduction_at_bad_prime_ram (Lv := locCycField p hl)
    (R := locCycRing p hl) he1 p hpe E C hC hc4ne hc4 hj
  have h2K : (2 : locCycField p hl) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (vcPoint C E Q) = l := by rw [addOrderOf_vcPoint, hQ]
  haveI hVell' := isElliptic_veluQuotientFull_nsmul_nf' L (C • E) hQ'
  have heq := veluQuotientFull_vcPoint_eq C E _ hQ two_ne_zero rfl
  have hres := jExp_neg_veluQuot_multRed_local_all L he1 p hpe E C hC hc4ne hc4 hj hl hodd
    h2K (vcPoint C E Q) hQ'
  have h1 : jExp p (C • veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
      = jExp p (veluQuotientFull (C • E)
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))) :=
    jExp_congr_j p _ _ (j_congr_curve heq)
  have h2 : jExp p (C • veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
      = jExp p (veluQuotientFull E
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :=
    jExp_variableChange p _ C
  omega

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★節点を埋める -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`JExpNegVeluStable` は成り立つ**（第 1446）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`Skeleton/GenEll/VeluSemistable.lean` に残っていた最後の `sorry` である。
★`l ∣ v_P(j)`（深い核）も `l ∤ v_P(j)`（`μ_l` 型）も
`jExp_neg_veluQuot_badPrime_all` が**一本で**扱う。 -/
theorem jExpNegVeluStable_holds : JExpNegVeluStable := by
  intro K _ _ P Y _ hss hj l hl hl5 _h48 _h864 _hdvd Q hQ _
  exact jExp_neg_veluQuot_badPrime_all P Y hss hj hl (by omega) Q hQ

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def jExp_lt_of_local.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(3·v(c₄) < v(Δ) なら jExp < 0——分岐版。★変数変換の u は打ち消し合う)",
    sectionId := "genell-lemma-3-5" }

def exists_veluW_deep_mem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核では Vélu の w は 𝔪 の中で取れる。★2 の可逆性は要らない)",
    sectionId := "genell-lemma-3-5" }

def jExp_neg_veluQuotient_bad_deep.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核の Vélu の商は jExp < 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def jExp_neg_veluQuotient_bad_mu.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(μ_l 型の核の Vélu の商は jExp < 0。★p ∣ l を許す)",
    sectionId := "genell-lemma-3-5" }

def jExp_neg_veluQuot_badPrime_all.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商も jExp < 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def jExpNegVeluStable_holds.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(乗法還元は l-同種で保たれる——深い核の場合も含めて閉じた)",
    sectionId := "genell-lemma-3-5" }

def jExpNegVeluStable_holds.needs : List ProofObligation :=
  [ .citation "[ABC3]" "primitiveRoot_or_deep_of_torsion_point(第 1410、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.primitiveRoot_or_deep_of_torsion_point") 1,
    .citation "[ABC3]" "isUnit_c4_add_240_deep(第 1412、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_c4_add_240_deep") 1,
    .citation "[ABC3]" "veluQuotientFull_tate_deep(第 1412、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_deep") 1,
    .citation "[ABC3]" "veluQuotientFull_tate_mu_K(第 1135、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu_K") 1,
    .citation "[ABC3]" "c4_c6_veluCurve_tate_field(第 1136 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.c4_c6_veluCurve_tate_field") 1,
    .citation "[ABC3]" "vAdd_Delta_tateCurveAt_pow(第 1057、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_Delta_tateCurveAt_pow") 1,
    .implicitStep
      ("★★★★★**2026-09-06(第 1446)**——`Skeleton/GenEll/VeluSemistable.lean` の" ++
       "最後の `sorry` を埋めた。☆鍵は 2 つである: " ++
       "(1) `jExp` は `j = c₄³/Δ` の付値なので**変数変換の `u` が打ち消し合う**" ++
       "（`3·4 − 12 = 0`）——だから第 1404 が要求した `vAdd (C₀.u) = 0` も" ++
       "核の座標の整性も要らない。" ++
       "(2) 深い核では `v ∈ 𝔪`（第 1412）に加えて **`w ∈ 𝔪`** も出る" ++
       "（`exists_veluW_of_inv` の証人は `veluU`・`veluV2` の和で、" ++
       "Tate 曲線は `a₃ = 0`・`a₄ ∈ 𝔪` だから）。" ++
       "☆すると `veluCurve` は剰余環で `E_q` に戻り `Δ(veluCurve) ∈ 𝔪`、" ++
       "`c₄` は単元なので `jExp = −v(Δ) < 0` である。" ++
       "★結果として仮説は**減った**——`p ∤ 6` も `l ∣ v_P(j)` も使っていない。") 17 ]

end ABC3.Found.GenEll
