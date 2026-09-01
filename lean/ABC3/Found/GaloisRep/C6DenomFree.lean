/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.C4DenomFree
import ABC3.Meta.Claim

/-!
# 第 1130 ブロック —— **`c₆` の恒等式の分母なし版**（`Found`、§3 の枠の節点 4 の後半）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`c6_velu_tate`（`Skeleton/GenEll/TateIsogeny.lean:192`）は `hu : ∀ i, IsUnit (1 − ζ^i)`
と `hDX` を仮説に取る。`p ∣ l` の悪い素点では前者が成り立たない。

★`c₆` の式には `veluU` と `veluV2·X` が現れる。`l` の重みは

| 量 | 分母を払う `l` の冪 |
|---|---|
| `x` | `l²` |
| `y` | `l³` |
| `veluV2 = 3x² + 2a₂x + a₄ − a₁y` | `l⁶` |
| `veluU = (−2y − a₁x − a₃)²` | `l⁶` |
| `veluV2·x` | `l⁸` |

なので**両辺に `l⁸` を掛ける**と全部が分母なしになる:

    `l⁸·c₆ + 504·l²·∑ veluV2DF + 3024·∑ (l²·veluUDF + 2·veluV2DF·tateXpairDF)
       = l¹⁴·c₆(E_{q^l})`

☆道は第 1128-1129 とまったく同じ 3 段である。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll ABC3.Skeleton.GenEll PowerSeries Finset

/-! ## ★★★★`veluU` の分母なし版 -/

section VeluU

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★**Vélu の `g^y` の分母なし版**——`X`・`Y` には `l²x`・`l³y` を入れる。 -/
noncomputable def veluGyDF (l : ℕ) (W : WeierstrassCurve R) (X Y : R) : R :=
  -2 * Y - W.a₁ * (l : R) * X - W.a₃ * (l : R) ^ 3

/-- ★★★★**Vélu の `u` の分母なし版**。 -/
noncomputable def veluUDF (l : ℕ) (W : WeierstrassCurve R) (X Y : R) : R :=
  (veluGyDF l W X Y) ^ 2

/-- ★★★★★★★★**無条件の代数恒等式**——`IsUnit` を一切使わない。 -/
theorem natCast_pow_mul_veluU (l : ℕ) (W : WeierstrassCurve R) (x y : R) :
    (l : R) ^ 6 * veluU W x y = veluUDF l W ((l : R) ^ 2 * x) ((l : R) ^ 3 * y) := by
  rw [veluU, veluGy, veluUDF, veluGyDF]
  ring

theorem map_veluGyDF {R' : Type} [CommRing R'] (φ : R →+* R') (l : ℕ)
    (W : WeierstrassCurve R) (X Y : R) :
    φ (veluGyDF l W X Y) = veluGyDF l (W.map φ) (φ X) (φ Y) := by
  rw [veluGyDF, veluGyDF, map_sub, map_sub, map_mul, map_mul, map_mul, map_mul,
    map_neg, map_ofNat, map_pow, map_natCast]
  rfl

/-- ★★★★★★**`veluUDF` は環準同型で移る**——★無条件（多項式だけ）。 -/
theorem map_veluUDF {R' : Type} [CommRing R'] (φ : R →+* R') (l : ℕ)
    (W : WeierstrassCurve R) (X Y : R) :
    φ (veluUDF l W X Y) = veluUDF l (W.map φ) (φ X) (φ Y) := by
  rw [veluUDF, veluUDF, map_pow, map_veluGyDF]

/-- ★★★★★★★★**Tate の座標での形**（`hu` の下で）。 -/
theorem natCast_pow_mul_veluU_tate [IsAdicComplete I R] [IsDomain R] {l : ℕ} {a : R}
    (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - a))
    (hpow : a ^ l = 1) (hsum : ∑ k ∈ range l, a ^ k = 0)
    (W : WeierstrassCurve R) (w q : R) (hq : q ∈ I) :
    (l : R) ^ 6 * veluU W (tateXpair a w q hq) (tateYpair a w q hq)
      = veluUDF l W (tateXpairDF l a w q hq) (tateYpairDF l a w q hq) := by
  rw [natCast_pow_mul_veluU, natCast_pow_mul_tateXpair hlu hu hpow hsum,
    natCast_pow_mul_tateYpair hlu hu hpow hsum]

end VeluU

/-! ## ★★★★★★`c₆` の分母なし左辺——定義を切る

☆式が長いので左辺に名前を付ける。★これで「環準同型で移る」が 1 本の補題になる。 -/

section Lhs

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★★**`c₆` の恒等式の分母なし左辺**（`l⁸` を掛けた形）。 -/
noncomputable def c6DFlhs (l : ℕ) (a q : R) (hq : q ∈ I) : R :=
  (l : R) ^ 8 * (tateCurveAt q hq).c₆
    + 504 * (l : R) ^ 2 * ∑ i ∈ (range l).erase 0,
        veluV2DF l (tateCurveAt q hq)
          (tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
          (tateYpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
    + 3024 * ∑ i ∈ (range l).erase 0,
        ((l : R) ^ 2 * veluUDF l (tateCurveAt q hq)
              (tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
              (tateYpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
          + 2 * (veluV2DF l (tateCurveAt q hq)
                  (tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
                  (tateYpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
                * tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq))

/-- ★★★★★★★★★★**`c6DFlhs` は環準同型で移る**——★`a` 側の単元性は要らない。 -/
theorem map_c6DFlhs {R' : Type} [CommRing R'] {J : Ideal R'} [IsAdicComplete J R']
    (φ : R →+* R') (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (l : ℕ) (a q : R) (hq : q ∈ I)
    (hwu : ∀ i : ℕ, IsUnit (1 - q * (a ^ i) ^ (l - 1))) :
    φ (c6DFlhs l a q hq)
      = c6DFlhs l (φ a) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  have hq' : φ q ∈ J := by simpa using hφ 1 q (by simpa using hq)
  have hcurve : (tateCurveAt q hq).map φ
      = tateCurveAt (φ q) (by simpa using hφ 1 q (by simpa using hq)) :=
    map_tateCurveAt φ hφ q hq
  have hmX : ∀ i : ℕ,
      φ (tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
        = tateXpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q)
            (by simpa using hφ 1 q (by simpa using hq)) := by
    intro i
    rw [map_tateXpairDF φ hφ _ _ _ _ _ (hwu i)]
    simp only [map_pow, map_mul]
  have hmY : ∀ i : ℕ,
      φ (tateYpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
        = tateYpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q)
            (by simpa using hφ 1 q (by simpa using hq)) := by
    intro i
    rw [map_tateYpairDF φ hφ _ _ _ _ _ (hwu i)]
    simp only [map_pow, map_mul]
  rw [c6DFlhs, c6DFlhs, map_add, map_add, map_mul, map_mul, map_mul, map_mul,
    map_sum, map_sum]
  rw [show φ ((l : R) ^ 8) = (l : R') ^ 8 from by rw [map_pow, map_natCast],
    show φ ((l : R) ^ 2) = (l : R') ^ 2 from by rw [map_pow, map_natCast],
    show φ (504 : R) = (504 : R') from by rw [map_ofNat],
    show φ (3024 : R) = (3024 : R') from by rw [map_ofNat],
    ← WeierstrassCurve.map_c₆, hcurve]
  have hs1 : ∑ i ∈ (range l).erase 0,
        φ (veluV2DF l (tateCurveAt q hq)
            (tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
            (tateYpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq))
      = ∑ i ∈ (range l).erase 0,
          veluV2DF l (tateCurveAt (φ q) hq')
            (tateXpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q) hq')
            (tateYpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q) hq') :=
    Finset.sum_congr rfl (fun i _ => by rw [map_veluV2DF, hcurve, hmX i, hmY i])
  have hs2 : ∑ i ∈ (range l).erase 0,
        φ ((l : R) ^ 2 * veluUDF l (tateCurveAt q hq)
              (tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
              (tateYpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
          + 2 * (veluV2DF l (tateCurveAt q hq)
                  (tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
                  (tateYpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq)
                * tateXpairDF l (a ^ i) (q * (a ^ i) ^ (l - 1)) q hq))
      = ∑ i ∈ (range l).erase 0,
          ((l : R') ^ 2 * veluUDF l (tateCurveAt (φ q) hq')
                (tateXpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q) hq')
                (tateYpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q) hq')
            + 2 * (veluV2DF l (tateCurveAt (φ q) hq')
                    (tateXpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q) hq')
                    (tateYpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q) hq')
                  * tateXpairDF l ((φ a) ^ i) (φ q * ((φ a) ^ i) ^ (l - 1)) (φ q) hq')) :=
    Finset.sum_congr rfl (fun i _ => by
      rw [map_add, map_mul, map_mul, map_mul, map_veluUDF, map_veluV2DF, hcurve, hmX i, hmY i,
        show φ ((l : R) ^ 2) = (l : R') ^ 2 from by rw [map_pow, map_natCast],
        show φ (2 : R) = (2 : R') from by rw [map_ofNat]])
  rw [hs1, hs2]

end Lhs

/-! ## ★★★★★★段 1 —— 体の上 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**段 1**——体 `K` の上の `PowerSeries K` での `c₆` の分母なし恒等式。 -/
theorem c6_velu_tateDF_field {K : Type} [Field K] [CharZero K] {l : ℕ} (hl : l.Prime)
    (hodd : l ≠ 2) {ζ : K} (hζ : IsPrimitiveRoot ζ l) :
    (l : PowerSeries K) ^ 8
        * (tateCurveAt (PowerSeries.X : PowerSeries K)
            (Ideal.mem_span_singleton_self _)).c₆
        + 504 * (l : PowerSeries K) ^ 2 * ∑ i ∈ (range l).erase 0,
            veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries K)
                (Ideal.mem_span_singleton_self _))
              (tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                PowerSeries.X (Ideal.mem_span_singleton_self _))
              (tateYpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                PowerSeries.X (Ideal.mem_span_singleton_self _))
        + 3024 * ∑ i ∈ (range l).erase 0,
            ((l : PowerSeries K) ^ 2
                * veluUDF l (tateCurveAt (PowerSeries.X : PowerSeries K)
                    (Ideal.mem_span_singleton_self _))
                  (tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                    (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                    PowerSeries.X (Ideal.mem_span_singleton_self _))
                  (tateYpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                    (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                    PowerSeries.X (Ideal.mem_span_singleton_self _))
              + 2 * (veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries K)
                      (Ideal.mem_span_singleton_self _))
                    (tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                      PowerSeries.X (Ideal.mem_span_singleton_self _))
                    (tateYpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                      PowerSeries.X (Ideal.mem_span_singleton_self _))
                  * tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                      PowerSeries.X (Ideal.mem_span_singleton_self _)))
      = (l : PowerSeries K) ^ 14
          * (tateCurveAt ((PowerSeries.X : PowerSeries K) ^ l)
              (Ideal.pow_mem_of_mem _ (Ideal.mem_span_singleton_self _) l hl.pos)).c₆ := by
  haveI : CharZero (PowerSeries K) := algebraRat.charZero _
  have hX : (PowerSeries.X : PowerSeries K)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries K)} := Ideal.mem_span_singleton_self _
  have hXl : ((PowerSeries.X : PowerSeries K) ^ l)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries K)} :=
    Ideal.pow_mem_of_mem _ hX l hl.pos
  have hζPS : IsPrimitiveRoot (PowerSeries.C ζ : PowerSeries K) l :=
    hζ.map_of_injective PowerSeries.C_injective
  have hlu : IsUnit ((l : PowerSeries K)) :=
    isUnit_natCast_powerSeries (by exact_mod_cast hl.pos.ne')
  have hζ0 : ζ ≠ 0 := by
    intro h
    have hh := hζ.pow_eq_one
    rw [h, zero_pow hl.pos.ne'] at hh
    exact zero_ne_one hh
  have hne : ∀ i ∈ (range l).erase 0, ζ ^ i ≠ 1 := by
    intro i hi h
    exact (Finset.mem_erase.1 hi).1
      (Nat.eq_zero_of_dvd_of_lt (hζ.dvd_of_pow_eq_one i h)
        (Finset.mem_range.1 (Finset.mem_erase.1 hi).2))
  have hu : ∀ i ∈ (range l).erase 0,
      IsUnit (1 - (PowerSeries.C ζ : PowerSeries K) ^ i) := by
    intro i hi
    rw [← map_pow]
    exact isUnit_one_sub_C (hne i hi)
  have hpowK : ∀ i : ℕ, (ζ ^ i) ^ l = 1 := fun i => by
    rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
  have hpowi : ∀ i : ℕ, ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ l = 1 := by
    intro i
    rw [← map_pow, ← map_pow, hpowK i, map_one]
  have hsumi : ∀ i ∈ (range l).erase 0,
      ∑ k ∈ range l, ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ k = 0 := by
    intro i hi
    have hcong : ∑ k ∈ range l, ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ k
        = PowerSeries.C (∑ k ∈ range l, (ζ ^ i) ^ k) := by
      rw [map_sum]
      exact Finset.sum_congr rfl (fun k _ => by rw [← map_pow, ← map_pow])
    rw [hcong, sum_pow_eq_zero_of_ne_one' (hpowK i) (hne i hi), map_zero]
  have h2 : (2 : PowerSeries K) ≠ 0 := by
    intro h
    have hcc := congrArg (PowerSeries.constantCoeff (R := K)) h
    rw [map_ofNat, map_zero] at hcc
    exact two_ne_zero hcc
  have hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
        (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
        PowerSeries.X hX ≠ 0 := by
    intro i hi
    rw [← map_pow]
    refine tateDXpair_C_ne_zero (pow_ne_zero _ hζ0) ?_ (hne i hi) ?_
    · intro h
      have hm1 : ζ ^ i = -1 := by linear_combination h
      have hpow : ((-1 : K)) ^ l = 1 := by rw [← hm1, hpowK i]
      rw [(hl.odd_of_ne_two hodd).neg_one_pow] at hpow
      exact absurd hpow (by norm_num)
    · exact Ideal.mul_mem_right _ _ hX
  have hc6 := c6_velu_tate hl hζPS hu PowerSeries.X hX hXl h2 hDX
  -- ★`veluV2` の和
  have hkeyV : ∑ i ∈ (range l).erase 0,
      veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
        (tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
          (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1)) PowerSeries.X hX)
        (tateYpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
          (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1)) PowerSeries.X hX)
      = (l : PowerSeries K) ^ 6 * ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
            (tateXpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
              (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
              PowerSeries.X hX)
            (tateYpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
              (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
              PowerSeries.X hX) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun i hi =>
      (natCast_pow_mul_veluV2_tate hlu (hu i hi) (hpowi i) (hsumi i hi) _ _ _ hX).symm)
  -- ★`veluU + 2·veluV2·X` の和
  have hkeyU : ∑ i ∈ (range l).erase 0,
      ((l : PowerSeries K) ^ 2
          * veluUDF l (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
            (tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
              (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
              PowerSeries.X hX)
            (tateYpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
              (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
              PowerSeries.X hX)
        + 2 * (veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
                (tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                  (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                  PowerSeries.X hX)
                (tateYpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                  (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                  PowerSeries.X hX)
              * tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                  (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                  PowerSeries.X hX))
      = (l : PowerSeries K) ^ 8 * ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
              (tateXpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
                (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                PowerSeries.X hX)
              (tateYpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
                (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                PowerSeries.X hX)
            + 2 * (veluV2 (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
                    (tateXpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
                      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                      PowerSeries.X hX)
                    (tateYpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
                      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                      PowerSeries.X hX)
                  * tateXpair ((PowerSeries.C ζ : PowerSeries K) ^ i)
                      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                      PowerSeries.X hX)) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    have hU := natCast_pow_mul_veluU_tate hlu (hu i hi) (hpowi i) (hsumi i hi)
      (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1)) PowerSeries.X hX
    have hV := natCast_pow_mul_veluV2_tate hlu (hu i hi) (hpowi i) (hsumi i hi)
      (tateCurveAt (PowerSeries.X : PowerSeries K) hX)
      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1)) PowerSeries.X hX
    have hXp := natCast_pow_mul_tateXpair hlu (hu i hi) (hpowi i) (hsumi i hi)
      (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1)) PowerSeries.X hX
    rw [← hU, ← hV, ← hXp]
    ring
  rw [hkeyV, hkeyU]
  linear_combination (l : PowerSeries K) ^ 8 * hc6

/-- ☆段 1 を `c6DFlhs` の形に畳んだもの。 -/
theorem c6DFlhs_field {K : Type} [Field K] [CharZero K] {l : ℕ} (hl : l.Prime)
    (hodd : l ≠ 2) {ζ : K} (hζ : IsPrimitiveRoot ζ l) :
    c6DFlhs l (PowerSeries.C ζ : PowerSeries K) PowerSeries.X
        (Ideal.mem_span_singleton_self _)
      = (l : PowerSeries K) ^ 14
          * (tateCurveAt ((PowerSeries.X : PowerSeries K) ^ l)
              (Ideal.pow_mem_of_mem _ (Ideal.mem_span_singleton_self _) l hl.pos)).c₆ := by
  rw [c6DFlhs]
  exact c6_velu_tateDF_field hl hodd hζ

/-! ## ★★★★★★★★段 2 —— 整域の `PowerSeries` へ降ろす -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**段 2**——`A` が標数 0 の整域なら `PowerSeries A` でも成り立つ。 -/
theorem c6DFlhs_powerSeries {A : Type} [CommRing A] [IsDomain A] [CharZero A] {l : ℕ}
    (hl : l.Prime) (hodd : l ≠ 2) {ζ : A} (hζ : IsPrimitiveRoot ζ l) :
    c6DFlhs l (PowerSeries.C ζ : PowerSeries A) PowerSeries.X
        (Ideal.mem_span_singleton_self _)
      = (l : PowerSeries A) ^ 14
          * (tateCurveAt ((PowerSeries.X : PowerSeries A) ^ l)
              (Ideal.pow_mem_of_mem _ (Ideal.mem_span_singleton_self _) l hl.pos)).c₆ := by
  have hψ : Function.Injective (algebraMap A (FractionRing A)) :=
    IsFractionRing.injective A (FractionRing A)
  refine PowerSeries.map_injective (algebraMap A (FractionRing A)) hψ ?_
  have hcont : ∀ (n : ℕ) (f : PowerSeries A),
      f ∈ (Ideal.span {(PowerSeries.X : PowerSeries A)}) ^ n →
      PowerSeries.map (algebraMap A (FractionRing A)) f
        ∈ (Ideal.span {(PowerSeries.X : PowerSeries (FractionRing A))}) ^ n :=
    fun n f hf => map_mem_span_X_pow _ n hf
  have hX : (PowerSeries.X : PowerSeries A)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries A)} := Ideal.mem_span_singleton_self _
  have hwu : ∀ i : ℕ,
      IsUnit (1 - PowerSeries.X * ((PowerSeries.C ζ : PowerSeries A) ^ i) ^ (l - 1)) :=
    fun i => isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries A)})
      (Ideal.mul_mem_right _ _ hX)
  rw [map_c6DFlhs _ hcont _ _ _ _ hwu, map_mul, map_pow, map_natCast,
    ← WeierstrassCurve.map_c₆, map_tateCurveAt _ hcont]
  simp only [PowerSeries.map_C, PowerSeries.map_X, map_pow]
  exact c6DFlhs_field hl hodd (hζ.map_of_injective hψ)

/-! ## ★★★★★★★★★★★★段 3 —— 完備環 `R` へ特殊化する -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `c₆` の恒等式の分母なし版**（第 1130、§3 の枠の節点 4 の後半）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★**`IsUnit (1 − ζ^i)` も `hDX` も仮説に置いていない**——
`p ∣ l` の悪い素点でもそのまま意味を持ち、そこで成り立つ。 -/
theorem c6_velu_tateDF {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) {ζ : R}
    (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    c6DFlhs l ζ q hq = (l : R) ^ 14 * (tateCurveAt (q ^ l) hql).c₆ := by
  have hX : (PowerSeries.X : PowerSeries R)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries R)} := Ideal.mem_span_singleton_self _
  have hwu : ∀ i : ℕ,
      IsUnit (1 - PowerSeries.X * ((PowerSeries.C ζ : PowerSeries R) ^ i) ^ (l - 1)) :=
    fun i => isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries R)})
      (Ideal.mul_mem_right _ _ hX)
  have hev : ∀ (n : ℕ) (f : PowerSeries R),
      f ∈ (Ideal.span {(PowerSeries.X : PowerSeries R)}) ^ n →
      evalAdicMapHom (RingHom.id R) q hq f ∈ I ^ n :=
    fun n f hf => evalAdicMapHom_mem_pow (RingHom.id R) q hq n f hf
  have hps := c6DFlhs_powerSeries (A := R) hl hodd hζ
  have hmain := congrArg (evalAdicMapHom (RingHom.id R) q hq) hps
  rw [map_c6DFlhs _ hev _ _ _ _ hwu, map_mul, map_pow, map_natCast,
    ← WeierstrassCurve.map_c₆, map_tateCurveAt _ hev] at hmain
  simp only [evalAdicMapHom_apply, evalAdicMap_C, evalAdicMap_X, RingHom.id_apply,
    map_pow] at hmain
  exact hmain

/-! ## ★出典の紐付け(`.src`) -/

def veluUDF.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の u の分母なし版——l⁶ を掛けて x・y の次数差を吸収する)",
    sectionId := "genell-lemma-3-5" }

def c6_velu_tateDF_field.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₆ の恒等式の E 版——体の PowerSeries の段)",
    sectionId := "genell-lemma-3-2" }

def c6_velu_tateDF.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₆ の恒等式の E 版。hu も hDX も取らない)",
    sectionId := "genell-lemma-3-2" }

def c6_velu_tateDF.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "c6_velu_tate(hu つきの原型、第 867、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.c6_velu_tate") 1,
    .citation "[ABC3]" "map_c6DFlhs(左辺は環準同型で移る、第 1130、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.map_c6DFlhs") 1,
    .citation "[ABC3]" "tateDXpair_C_ne_zero(側条件は定数項で出る、第 1128、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDXpair_C_ne_zero") 1,
    .citation "[mathlib]" "WeierstrassCurve.map_c₆(c₆ は環準同型で移る)"
      (.inMathlib "WeierstrassCurve.map_c₆") 1,
    .implicitStep
      ("★★**2026-09-01（第 1130）**——`c₆` の式には `veluU` と `veluV2·x` が現れる。" ++
       "☆`x` は `l²`、`y` は `l³`、`veluV2`・`veluU` は `l⁶`、`veluV2·x` は `l⁸` で" ++
       "分母が払えるので、**両辺に `l⁸` を掛ける**と全部が分母なしになる。" ++
       "★式が長いので左辺に `c6DFlhs` と名前を付け、" ++
       "「環準同型で移る」を 1 本の補題（`map_c6DFlhs`）にした。") 5 ]

def natCast_pow_mul_veluU.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★**2026-09-01（第 1130）**——`l⁶·veluU` を `l²x`・`l³y` だけで書いた。" ++
       "☆`veluU = (−2y − a₁x − a₃)²` なので `l³` を括り出して 2 乗すればよい。" ++
       "★これは**無条件の代数恒等式**で `ring` で閉じる。") 3 ]

end ABC3.Found.GaloisRep
