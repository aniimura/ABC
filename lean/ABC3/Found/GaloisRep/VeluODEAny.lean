/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.GenericParam
import ABC3.Found.GaloisRep.VeluDYDenomFree
import ABC3.Skeleton.GenEll.TateODE
import ABC3.Meta.Claim

/-!
# 第 1144 ブロック —— **`DY = veluV2` を `DX ≠ 0` なしで**（`Found`、`l = 2` の枝の節点 1）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★これは何か

`veluV2_eq_tateDYpair`（`TateODE.lean:89`）は `hDX : tateDXpair ≠ 0` を仮説に取る
——`tate_ode_mul` が与えるのは `DX·(DY − veluV2) = 0` だからである。

★`l = 2` では点が **2-捩れ**なので `DX = 2Y + X = 0` になり、この形は使えない。

☆本ブロックは **`hDX` を取らない版**を証明する。道は第 1143 の一般の径数である:

| 段 | 場所 | 何をするか |
|---|---|---|
| 1 | `PowerSeries U`（`U = ℤ[T]` を `T(1−T)` で局所化） | `DX` の定数項 `T(1+T)(1−T)⁻³ ≠ 0` |
| 2 | 在庫の `veluV2_eq_tateDYpair` を当てる | `hu`・`hDX` とも成り立つ |
| 3 | `evalAdicMapHom (genericSpec a) q hq` で特殊化 | `T ↦ a`、`X ↦ q` |

★`a` が単元で `1 − a` も単元なら `T ↦ a` は環準同型になる（第 1143）。
☆`w` は `a·w = q` で決まるので、特殊化の像が `w` に一致する。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll ABC3.Skeleton.GenEll PowerSeries Finset

/-! ## ★整域版の `DX ≠ 0` -/

/-- ★★★★★★**整域でも `DX ≠ 0` は定数項で出る**（第 1144）。

☆第 1128 の体版を整域に広げたもの。`Ring.inverse (1−α)` は単元だから `0` でない。 -/
theorem tateDXpair_C_ne_zero_domain {A : Type} [CommRing A] [IsDomain A] {α : A}
    (hα0 : α ≠ 0) (hαneg : 1 + α ≠ 0) (hu : IsUnit (1 - α))
    {w : PowerSeries A} (hw : w ∈ Ideal.span {(PowerSeries.X : PowerSeries A)}) :
    tateDXpair (PowerSeries.C α) w PowerSeries.X (Ideal.mem_span_singleton_self _) ≠ 0 := by
  intro h0
  have hcc := congrArg (PowerSeries.constantCoeff (R := A)) h0
  rw [tateDXpair, map_sub, map_add, map_add] at hcc
  have h1 : PowerSeries.constantCoeff (tateDXterm (PowerSeries.C α : PowerSeries A))
      = tateDXterm α := by
    rw [← map_tateDXterm (PowerSeries.C (R := A)) hu]; simp
  have h2 : PowerSeries.constantCoeff
      (tateDXtail (PowerSeries.C α : PowerSeries A) PowerSeries.X
        (Ideal.mem_span_singleton_self _)) = 0 :=
    constantCoeff_of_mem_span_X (tateDXtail_mem _ _ _)
  have h3mem : tateDXterm w ∈ Ideal.span {(PowerSeries.X : PowerSeries A)} := by
    have hh := tateDXterm_mem_pow (I := Ideal.span {(PowerSeries.X : PowerSeries A)})
      (k := 1) (t := w) (by simpa using hw)
    simpa using hh
  have h3 : PowerSeries.constantCoeff (tateDXterm w) = 0 := constantCoeff_of_mem_span_X h3mem
  have h4 : PowerSeries.constantCoeff
      (tateDXtail w PowerSeries.X (Ideal.mem_span_singleton_self _)) = 0 :=
    constantCoeff_of_mem_span_X (tateDXtail_mem _ _ _)
  rw [h1, h2, h3, h4] at hcc
  simp only [add_zero, sub_zero, map_zero] at hcc
  rw [tateDXterm] at hcc
  have hinv : Ring.inverse (1 - α) ≠ 0 := Ring.isUnit_iff_inverse_ne_zero.mp hu
  exact (mul_ne_zero (mul_ne_zero hα0 hαneg) (pow_ne_zero _ hinv)) hcc

/-! ## ★環準同型で移る 2 本（第 221 に無かった分） -/

/-- ★★★★**`veluV2` は環準同型で移る**——★無条件。 -/
theorem map_veluV2 {R R' : Type} [CommRing R] [CommRing R'] (φ : R →+* R')
    (W : WeierstrassCurve R) (x y : R) :
    φ (veluV2 W x y) = veluV2 (W.map φ) (φ x) (φ y) := by
  rw [veluV2, veluV2, veluGx, veluGx, map_sub, map_add, map_add, map_mul, map_mul,
    map_mul, map_mul, map_pow]
  simp [WeierstrassCurve.map, map_ofNat]

/-- ★★★★**`DY` の対は環準同型で移る**（単元の下で）。 -/
theorem map_tateDYpair {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {J : Ideal R'}
    [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hwu : IsUnit (1 - w)) :
    φ (tateDYpair a w q hq)
      = tateDYpair (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateDYpair, tateDYpair, map_add, map_add, map_add, map_add, map_add,
    map_tateDYterm φ ha, map_tateDYtail φ hφ a q hq,
    map_tateDXterm φ hwu, map_tateDXtail φ hφ w q hq,
    map_tateDYterm φ hwu, map_tateDYtail φ hφ w q hq]

/-! ## ★★★★径数の逆元 -/

/-- ☆`T` の逆元。 -/
noncomputable def genericTinv : GenericParamRing := ↑(isUnit_genericT.unit⁻¹)

theorem genericT_mul_inv : genericT * genericTinv = 1 := by
  have h : (isUnit_genericT.unit : GenericParamRing) * genericTinv = 1 := by
    rw [genericTinv]
    exact_mod_cast (isUnit_genericT.unit).mul_inv
  rwa [isUnit_genericT.unit_spec] at h

theorem genericT_ne_zero : genericT ≠ 0 := by
  intro h
  have hm := genericT_mul_inv
  rw [h, zero_mul] at hm
  exact zero_ne_one hm

/-! ## ★★★★★★★★★★段 1-2 —— 一般の径数での恒等式 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**一般の径数での `DY = veluV2`**（第 1144）。

☆`PowerSeries U` では `DX` の定数項が `T(1+T)(1−T)⁻³ ≠ 0` なので、
在庫の `veluV2_eq_tateDYpair` がそのまま効く。 -/
theorem veluV2_eq_tateDYpair_generic :
    veluV2 (tateCurveAt (PowerSeries.X : PowerSeries GenericParamRing)
        (Ideal.mem_span_singleton_self _))
      (tateXpair (PowerSeries.C genericT)
        (PowerSeries.C genericTinv * PowerSeries.X) PowerSeries.X
        (Ideal.mem_span_singleton_self _))
      (tateYpair (PowerSeries.C genericT)
        (PowerSeries.C genericTinv * PowerSeries.X) PowerSeries.X
        (Ideal.mem_span_singleton_self _))
      = tateDYpair (PowerSeries.C genericT)
        (PowerSeries.C genericTinv * PowerSeries.X) PowerSeries.X
        (Ideal.mem_span_singleton_self _) := by
  have hX : (PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mem_span_singleton_self _
  have hwmem : (PowerSeries.C genericTinv * PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mul_mem_left _ _ hX
  have haw : (PowerSeries.C genericT : PowerSeries GenericParamRing)
      * (PowerSeries.C genericTinv * PowerSeries.X) = PowerSeries.X := by
    rw [← mul_assoc, ← map_mul, genericT_mul_inv, map_one, one_mul]
  have ha : IsUnit (1 - (PowerSeries.C genericT : PowerSeries GenericParamRing)) := by
    have h : (1 : PowerSeries GenericParamRing) - PowerSeries.C genericT
        = PowerSeries.C (1 - genericT) := by rw [map_sub, map_one]
    rw [h]
    exact isUnit_C_of_isUnit isUnit_one_sub_genericT
  have hwu : IsUnit (1 - (PowerSeries.C genericTinv * PowerSeries.X
      : PowerSeries GenericParamRing)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)}) hwmem
  have hDX := tateDXpair_C_ne_zero_domain (A := GenericParamRing)
    genericT_ne_zero one_add_genericT_ne_zero isUnit_one_sub_genericT hwmem
  exact veluV2_eq_tateDYpair _ _ PowerSeries.X hX haw ha hwu hDX

/-! ## ★★★★★★★★★★★★段 3 —— 特殊化（★節点 1 の本体） -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `veluV2 = DY`——`DX ≠ 0` を取らない**（第 1144、`l = 2` の枝の節点 1）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★**`hDX` を仮説に置いていない**——`l = 2` の 2-捩れ点（`DX = 0`）でも成り立つ。
☆代わりに `IsUnit a` を置く（`a = ζ^i` は 1 の冪根なので自動）。 -/
theorem veluV2_eq_tateDYpair_any {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hau : IsUnit a) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    veluV2 (tateCurveAt q hq) (tateXpair a w q hq) (tateYpair a w q hq)
      = tateDYpair a w q hq := by
  have haa : IsUnit (a * (1 - a)) := hau.mul ha
  have hX : (PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mem_span_singleton_self _
  have hwmem : (PowerSeries.C genericTinv * PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mul_mem_left _ _ hX
  have haU : IsUnit (1 - (PowerSeries.C genericT : PowerSeries GenericParamRing)) := by
    have h : (1 : PowerSeries GenericParamRing) - PowerSeries.C genericT
        = PowerSeries.C (1 - genericT) := by rw [map_sub, map_one]
    rw [h]
    exact isUnit_C_of_isUnit isUnit_one_sub_genericT
  have hwU : IsUnit (1 - (PowerSeries.C genericTinv * PowerSeries.X
      : PowerSeries GenericParamRing)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)}) hwmem
  have hev := evalAdicMapHom_mem_pow (genericSpec a haa) q hq
  have h := congrArg (evalAdicMapHom (genericSpec a haa) q hq) veluV2_eq_tateDYpair_generic
  rw [map_veluV2, map_tateXpair _ hev _ _ _ hX haU hwU,
    map_tateYpair _ hev _ _ _ hX haU hwU, map_tateCurveAt _ hev,
    map_tateDYpair _ hev _ _ _ hX haU hwU] at h
  simp only [evalAdicMapHom_apply, evalAdicMap_C, evalAdicMap_X, genericSpec_T, map_mul] at h
  have hwe : genericSpec a haa genericTinv * q = w := by
    have hone : a * genericSpec a haa genericTinv = 1 := by
      have h1 : genericSpec a haa genericT * genericSpec a haa genericTinv = 1 := by
        rw [← map_mul, genericT_mul_inv, map_one]
      rwa [genericSpec_T] at h1
    refine hau.mul_left_cancel ?_
    rw [← mul_assoc, hone, one_mul, haw]
  rw [hwe] at h
  exact h

/-! ## ★★★★★★★★第 1145 —— 高階の ODE も `hDX` なしで

☆`tate_d2x`・`tate_d3x`・`tate_d4x` も同じく `DX ≠ 0` で割っている。
★同じ径数の道でどれも `hDX` が外れる。要るのは各部品の `map_*` だけである。 -/

section HigherODE

variable {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {J : Ideal R'}

theorem map_tateD2Xterm' (φ : R →+* R') {t : R} (ht : IsUnit (1 - t)) :
    φ (tateD2Xterm t) = tateD2Xterm (φ t) := by
  rw [tateD2Xterm, tateD2Xterm, map_mul, map_mul, map_pow, map_ring_inverse φ ht,
    map_sub, map_one]
  congr 1
  simp [map_ofNat]

theorem map_tateD3Xterm' (φ : R →+* R') {t : R} (ht : IsUnit (1 - t)) :
    φ (tateD3Xterm t) = tateD3Xterm (φ t) := by
  rw [tateD3Xterm, tateD3Xterm, map_mul, map_mul, map_pow, map_ring_inverse φ ht,
    map_sub, map_one]
  congr 1
  simp [map_ofNat]

theorem map_tateD4Xterm' (φ : R →+* R') {t : R} (ht : IsUnit (1 - t)) :
    φ (tateD4Xterm t) = tateD4Xterm (φ t) := by
  rw [tateD4Xterm, tateD4Xterm, map_mul, map_mul, map_pow, map_ring_inverse φ ht,
    map_sub, map_one]
  congr 1
  simp [map_ofNat]

theorem map_tateD2Xtail [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n) (u q : R) (hq : q ∈ I) :
    φ (tateD2Xtail u q hq)
      = tateD2Xtail (φ u) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateD2Xtail, map_adicSum φ hφ, tateD2Xtail]
  refine adicSum_congr _ _ fun n => ?_
  rw [map_tateD2Xterm' φ (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), map_mul, map_pow]

theorem map_tateD3Xtail [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n) (u q : R) (hq : q ∈ I) :
    φ (tateD3Xtail u q hq)
      = tateD3Xtail (φ u) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateD3Xtail, map_adicSum φ hφ, tateD3Xtail]
  refine adicSum_congr _ _ fun n => ?_
  rw [map_tateD3Xterm' φ (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), map_mul, map_pow]

theorem map_tateD4Xtail [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n) (u q : R) (hq : q ∈ I) :
    φ (tateD4Xtail u q hq)
      = tateD4Xtail (φ u) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateD4Xtail, map_adicSum φ hφ, tateD4Xtail]
  refine adicSum_congr _ _ fun n => ?_
  rw [map_tateD4Xterm' φ (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), map_mul, map_pow]

theorem map_tateDXpair [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hwu : IsUnit (1 - w)) :
    φ (tateDXpair a w q hq)
      = tateDXpair (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateDXpair, tateDXpair, map_sub, map_add, map_add,
    map_tateDXterm φ ha, map_tateDXtail φ hφ a q hq,
    map_tateDXterm φ hwu, map_tateDXtail φ hφ w q hq]

theorem map_tateD2Xpair [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hwu : IsUnit (1 - w)) :
    φ (tateD2Xpair a w q hq)
      = tateD2Xpair (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateD2Xpair, tateD2Xpair, map_add, map_add, map_add,
    map_tateD2Xterm' φ ha, map_tateD2Xtail φ hφ a q hq,
    map_tateD2Xterm' φ hwu, map_tateD2Xtail φ hφ w q hq]

theorem map_tateD3Xpair [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hwu : IsUnit (1 - w)) :
    φ (tateD3Xpair a w q hq)
      = tateD3Xpair (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateD3Xpair, tateD3Xpair, map_sub, map_add, map_add,
    map_tateD3Xterm' φ ha, map_tateD3Xtail φ hφ a q hq,
    map_tateD3Xterm' φ hwu, map_tateD3Xtail φ hφ w q hq]

theorem map_tateD4Xpair [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hwu : IsUnit (1 - w)) :
    φ (tateD4Xpair a w q hq)
      = tateD4Xpair (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateD4Xpair, tateD4Xpair, map_add, map_add, map_add,
    map_tateD4Xterm' φ ha, map_tateD4Xtail φ hφ a q hq,
    map_tateD4Xterm' φ hwu, map_tateD4Xtail φ hφ w q hq]

end HigherODE

/-! ## ★★★★★★★★★★★★高階の ODE の `hDX` なし版（第 1145） -/

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★**`D²X = 6X² + X + 2a₄`——`DX ≠ 0` を取らない**（第 1145）。 -/
theorem tate_d2x_any {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hau : IsUnit a) (ha : IsUnit (1 - a)) :
    tateD2Xpair a w q hq
      = 6 * tateXpair a w q hq ^ 2 + tateXpair a w q hq
        + 2 * (tateCurveAt q hq).a₄ := by
  have haa : IsUnit (a * (1 - a)) := hau.mul ha
  have hX : (PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mem_span_singleton_self _
  have hwmem : (PowerSeries.C genericTinv * PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mul_mem_left _ _ hX
  have hawU : (PowerSeries.C genericT : PowerSeries GenericParamRing)
      * (PowerSeries.C genericTinv * PowerSeries.X) = PowerSeries.X := by
    rw [← mul_assoc, ← map_mul, genericT_mul_inv, map_one, one_mul]
  have haU : IsUnit (1 - (PowerSeries.C genericT : PowerSeries GenericParamRing)) := by
    have h : (1 : PowerSeries GenericParamRing) - PowerSeries.C genericT
        = PowerSeries.C (1 - genericT) := by rw [map_sub, map_one]
    rw [h]; exact isUnit_C_of_isUnit isUnit_one_sub_genericT
  have hwU : IsUnit (1 - (PowerSeries.C genericTinv * PowerSeries.X
      : PowerSeries GenericParamRing)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)}) hwmem
  have hDXU := tateDXpair_C_ne_zero_domain (A := GenericParamRing)
    genericT_ne_zero one_add_genericT_ne_zero isUnit_one_sub_genericT hwmem
  have hev := evalAdicMapHom_mem_pow (genericSpec a haa) q hq
  have h := congrArg (evalAdicMapHom (genericSpec a haa) q hq)
    (tate_d2x (PowerSeries.C genericT)
      (PowerSeries.C genericTinv * PowerSeries.X) PowerSeries.X hX hawU haU hwU hDXU)
  simp only [map_tateD2Xpair _ hev _ _ _ hX haU hwU,
    map_tateD3Xpair _ hev _ _ _ hX haU hwU, map_tateD4Xpair _ hev _ _ _ hX haU hwU,
    map_tateDXpair _ hev _ _ _ hX haU hwU, map_tateXpair _ hev _ _ _ hX haU hwU,
    map_tateCurveAt_a4 _ hev, map_add, map_mul, map_pow, map_ofNat,
    evalAdicMapHom_apply, evalAdicMap_C, evalAdicMap_X, genericSpec_T] at h
  have hwe : genericSpec a haa genericTinv * q = w := by
    have hone : a * genericSpec a haa genericTinv = 1 := by
      have h1 : genericSpec a haa genericT * genericSpec a haa genericTinv = 1 := by
        rw [← map_mul, genericT_mul_inv, map_one]
      rwa [genericSpec_T] at h1
    refine hau.mul_left_cancel ?_
    rw [← mul_assoc, hone, one_mul, haw]
  rw [hwe] at h
  exact h

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★**`D³X = 12X·DX + DX`——`DX ≠ 0` を取らない**（第 1145）。 -/
theorem tate_d3x_any {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hau : IsUnit a) (ha : IsUnit (1 - a)) :
    tateD3Xpair a w q hq
      = 12 * tateXpair a w q hq * tateDXpair a w q hq + tateDXpair a w q hq := by
  have haa : IsUnit (a * (1 - a)) := hau.mul ha
  have hX : (PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mem_span_singleton_self _
  have hwmem : (PowerSeries.C genericTinv * PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mul_mem_left _ _ hX
  have hawU : (PowerSeries.C genericT : PowerSeries GenericParamRing)
      * (PowerSeries.C genericTinv * PowerSeries.X) = PowerSeries.X := by
    rw [← mul_assoc, ← map_mul, genericT_mul_inv, map_one, one_mul]
  have haU : IsUnit (1 - (PowerSeries.C genericT : PowerSeries GenericParamRing)) := by
    have h : (1 : PowerSeries GenericParamRing) - PowerSeries.C genericT
        = PowerSeries.C (1 - genericT) := by rw [map_sub, map_one]
    rw [h]; exact isUnit_C_of_isUnit isUnit_one_sub_genericT
  have hwU : IsUnit (1 - (PowerSeries.C genericTinv * PowerSeries.X
      : PowerSeries GenericParamRing)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)}) hwmem
  have hDXU := tateDXpair_C_ne_zero_domain (A := GenericParamRing)
    genericT_ne_zero one_add_genericT_ne_zero isUnit_one_sub_genericT hwmem
  have hev := evalAdicMapHom_mem_pow (genericSpec a haa) q hq
  have h := congrArg (evalAdicMapHom (genericSpec a haa) q hq)
    (tate_d3x (PowerSeries.C genericT)
      (PowerSeries.C genericTinv * PowerSeries.X) PowerSeries.X hX hawU haU hwU hDXU)
  simp only [map_tateD2Xpair _ hev _ _ _ hX haU hwU,
    map_tateD3Xpair _ hev _ _ _ hX haU hwU, map_tateD4Xpair _ hev _ _ _ hX haU hwU,
    map_tateDXpair _ hev _ _ _ hX haU hwU, map_tateXpair _ hev _ _ _ hX haU hwU,
    map_tateCurveAt_a4 _ hev, map_add, map_mul, map_pow, map_ofNat,
    evalAdicMapHom_apply, evalAdicMap_C, evalAdicMap_X, genericSpec_T] at h
  have hwe : genericSpec a haa genericTinv * q = w := by
    have hone : a * genericSpec a haa genericTinv = 1 := by
      have h1 : genericSpec a haa genericT * genericSpec a haa genericTinv = 1 := by
        rw [← map_mul, genericT_mul_inv, map_one]
      rwa [genericSpec_T] at h1
    refine hau.mul_left_cancel ?_
    rw [← mul_assoc, hone, one_mul, haw]
  rw [hwe] at h
  exact h

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★**`D⁴X = 12X·D²X + 12(DX)² + D²X`——`DX ≠ 0` を取らない**（第 1145）。 -/
theorem tate_d4x_any {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hau : IsUnit a) (ha : IsUnit (1 - a)) :
    tateD4Xpair a w q hq
      = 12 * tateXpair a w q hq * tateD2Xpair a w q hq
        + 12 * tateDXpair a w q hq ^ 2 + tateD2Xpair a w q hq := by
  have haa : IsUnit (a * (1 - a)) := hau.mul ha
  have hX : (PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mem_span_singleton_self _
  have hwmem : (PowerSeries.C genericTinv * PowerSeries.X : PowerSeries GenericParamRing)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)} :=
    Ideal.mul_mem_left _ _ hX
  have hawU : (PowerSeries.C genericT : PowerSeries GenericParamRing)
      * (PowerSeries.C genericTinv * PowerSeries.X) = PowerSeries.X := by
    rw [← mul_assoc, ← map_mul, genericT_mul_inv, map_one, one_mul]
  have haU : IsUnit (1 - (PowerSeries.C genericT : PowerSeries GenericParamRing)) := by
    have h : (1 : PowerSeries GenericParamRing) - PowerSeries.C genericT
        = PowerSeries.C (1 - genericT) := by rw [map_sub, map_one]
    rw [h]; exact isUnit_C_of_isUnit isUnit_one_sub_genericT
  have hwU : IsUnit (1 - (PowerSeries.C genericTinv * PowerSeries.X
      : PowerSeries GenericParamRing)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries GenericParamRing)}) hwmem
  have hDXU := tateDXpair_C_ne_zero_domain (A := GenericParamRing)
    genericT_ne_zero one_add_genericT_ne_zero isUnit_one_sub_genericT hwmem
  have hev := evalAdicMapHom_mem_pow (genericSpec a haa) q hq
  have h := congrArg (evalAdicMapHom (genericSpec a haa) q hq)
    (tate_d4x (PowerSeries.C genericT)
      (PowerSeries.C genericTinv * PowerSeries.X) PowerSeries.X hX hawU haU hwU hDXU)
  simp only [map_tateD2Xpair _ hev _ _ _ hX haU hwU,
    map_tateD3Xpair _ hev _ _ _ hX haU hwU, map_tateD4Xpair _ hev _ _ _ hX haU hwU,
    map_tateDXpair _ hev _ _ _ hX haU hwU, map_tateXpair _ hev _ _ _ hX haU hwU,
    map_tateCurveAt_a4 _ hev, map_add, map_mul, map_pow, map_ofNat,
    evalAdicMapHom_apply, evalAdicMap_C, evalAdicMap_X, genericSpec_T] at h
  have hwe : genericSpec a haa genericTinv * q = w := by
    have hone : a * genericSpec a haa genericTinv = 1 := by
      have h1 : genericSpec a haa genericT * genericSpec a haa genericTinv = 1 := by
        rw [← map_mul, genericT_mul_inv, map_one]
      rwa [genericSpec_T] at h1
    refine hau.mul_left_cancel ?_
    rw [← mul_assoc, hone, one_mul, haw]
  rw [hwe] at h
  exact h

/-! ## ★出典の紐付け(`.src`) -/

def veluV2_eq_tateDYpair_any.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(veluV2 = DY——★DX ≠ 0 を取らない。l = 2 でも成り立つ)",
    sectionId := "genell-lemma-3-2" }

def tateDXpair_C_ne_zero_domain.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(整域でも DX ≠ 0 は定数項で出ること)",
    sectionId := "genell-def-3-3" }

def veluV2_eq_tateDYpair_any.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "veluV2_eq_tateDYpair(hDX つきの原型、第 851、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.veluV2_eq_tateDYpair") 1,
    .citation "[ABC3]" "genericSpec(径数の特殊化、第 1143、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.genericSpec") 1,
    .citation "[ABC3]" "evalAdicMapHom_mem_pow(特殊化の連続性、第 1127、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdicMapHom_mem_pow") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1144）**——`tate_ode_mul` が与えるのは " ++
       "`DX·(DY − veluV2) = 0` なので、`DX = 0` の点（`l = 2` の 2-捩れ点）では" ++
       "情報が出ない。☆しかし `a` を一般の径数 `T` に取れば " ++
       "`DX` の定数項は `T(1+T)(1−T)⁻³ ≠ 0` になるので在庫がそのまま効き、" ++
       "`T ↦ a` で特殊化すればよい。★これが `l = 2` の枝の節点 1 である。") 10 ]

end ABC3.Found.GaloisRep
