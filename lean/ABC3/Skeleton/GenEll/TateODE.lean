/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateDSeries
import ABC3.Found.GaloisRep.MuDYSum
import ABC3.Found.GaloisRep.DualTate
import ABC3.Found.GaloisRep.MuD2XSum
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Meta.Claim

/-!
# `Skeleton` —— **★★★★★★★★★★Tate 曲線の ODE と `v = ∑_ζ DY`**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★葉 1 は 2 つの節点に落ちた

第 846 で `v = ∑_ζ DY` が見えた（数値確認済み）。その根拠は ODE 一つである:

    `DY = 3X² − Y + a₄`        …… `tate_ode`（本ファイル）

これを `veluV2 = 3x² + a₄ − y`（`veluV2_tateCurveAt`、証明済み）に入れると

    `veluV2(X_i, Y_i) = DY(ζ^i)`

がそのまま出る（`veluV2_eq_tateDYpair`、★本ファイルで**証明済み**）。
★★★★**`X²` は消え、畳み込み（Besge の恒等式）は要らなくなった**。

## ★★★★★★`tate_ode` はどう出るか

`tate_equation`（証明済み）は `(2Y+X)² = 4X³ + X² + 4a₄X + 4a₆` であり、
`tateDXpair_eq`（第 846、証明済み）は `DX = 2Y + X` である。よって

    `(DX)² = 4X³ + X² + 4a₄X + 4a₆`

万有な環 `TateUniv` の上でこれを `D`（第 845）で微分すると

    `2·DX·D²X = (12X² + 2X + 4a₄)·DX`

★`TateUniv` は整域（UFD `ℤ[A,W]` の局所化）なので `2DX` で割れて

    `D²X = 6X² + X + 2a₄`,  すなわち  `DY = 3X² − Y + a₄`
    （`D²X = D(2Y+X) = 2DY + DX = 2DY + 2Y + X` を使う）

☆級数は `TateUniv` に無く完備化の中にあるので、実際には `q^N` で切り詰めて
`TateUniv` の中で整除性を見る——`tate_equation` を通した道（第 222–229）と同じ形である。

## ★★★★★定数項の指標和

`q = 0` では `DY(u,0) = u²(2+u)/(1−u)⁴` なので、`v` の定数項は

    `240 · ∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1`     …… ★★**第 849 で証明済み**

★`l = 2, 3, 5, 7, 11, 13` で厳密一致（2026-08-31）。
☆第 835 が測った `v(0) = (l⁴−1)/240` と同じ値である。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep ABC3.Found.GenEll Finset

/-- ★★★★★★★★★★**Tate 曲線の ODE**——第 850-851 で**証明済み**になった。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★`tate_ode_mul`（`Found/GaloisRep/DualTate.lean`）は双対数 `R[ε]` の中で
`tate_equation` を `(a+εa, w−εw, q)` に適用し、`ε` 成分を取って

    `DX · (DY − (3X² − Y + a₄)) = 0`

を与える。☆あとは `DX ≠ 0` で割るだけである（`R` は整域）。 -/
theorem tate_ode {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    (a w q : R) (hq : q ∈ I) (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hDX : tateDXpair a w q hq ≠ 0) :
    tateDYpair a w q hq
      = 3 * tateXpair a w q hq ^ 2 - tateYpair a w q hq + (tateCurveAt q hq).a₄ := by
  have h := tate_ode_mul a w q hq haw ha hw
  rcases mul_eq_zero.1 h with h1 | h2
  · exact absurd h1 hDX
  · exact sub_eq_zero.1 h2

/-- ★★★★★★★★★★**`veluV2` は `DY` そのものである**——★ODE から**直ちに**出る。

これが第 846 の発見「`v = ∑_ζ DY`」の中身であり、`X²` が消える理由である。 -/
theorem veluV2_eq_tateDYpair {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) (hDX : tateDXpair a w q hq ≠ 0) :
    veluV2 (tateCurveAt q hq) (tateXpair a w q hq) (tateYpair a w q hq)
      = tateDYpair a w q hq := by
  rw [veluV2_tateCurveAt, tate_ode a w q hq haw ha hw hDX]
  ring

/-- ★★★★★★★★**`v = ∑_{i≠0} DY(ζ^i)`**——`c4_velu_tate` の左辺の `v` を置き換える形。 -/
theorem sum_veluV2_eq_sum_tateDYpair {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : 0 < l) {ζ : R} (hζl : ζ ^ l = 1)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    (∑ i ∈ (range l).erase 0,
        veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = ∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
    rw [← pow_succ']
    rw [Nat.sub_add_cancel hl, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have haw : (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hwu : IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) :=
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  exact veluV2_eq_tateDYpair _ _ q hq haw (hu i hi) hwu (hDX i hi)

/-- ★★★★★★★★★★**葉 1 の定数項は閉じた**——
`240·∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1`。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★第 847-849 で**証明済み**になった（`Found/GaloisRep/MuDYSum.lean`）。
道は `muPow`/`muM` の機械だけであり、`p₄` に新しい種は要らなかった。

☆`l = 2, 3, 5, 7, 11, 13, 17` で厳密一致（2026-08-31、第 846）。 -/
theorem sum_mu_dyterm {F : Type} [Field F] [CharZero F] {l : ℕ} (hl : l.Prime) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) :
    240 * (∑ i ∈ (range l).erase 0, tateDYterm (ζ ^ i)) = (l : F) ^ 4 - 1 :=
  sum_mu_dyterm_field hl hζ

/-! ## ★★★★★★★★`∑_ζ D²X` への還元 -/

/-- **[GenEll] 葉 1 の残り (1)**——`∑_{i≠0} DX(ζ^i) = 0`。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★理由: `∑_ζ X(ζ)` は `q` だけの関数であり、`D` （`u∂_u`）で消える。
☆定数項は `sum_mu_dxterm_field`（第 852、証明済み）であり、
`q^N` 係数は `∑_{d∣N}d²∑_{i≠0}(ζ^{id} − ζ^{-id}) = 0`（`i ↦ l−i` の対称性）。 -/
theorem sum_mu_dxpair_zero {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (h2 : (2 : R) ≠ 0) :
    ∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq = 0 := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hpow : ∀ i ∈ (range l).erase 0, (ζ ^ i) ^ (l - 1) = ζ ^ (l - i) := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hl1 : 1 ≤ l := hl.one_lt.le
    have hi1 : 1 ≤ i := Nat.one_le_iff_ne_zero.2 hi0
    have e1 : i * (l - 1) + i = i * l := by
      calc i * (l - 1) + i = i * ((l - 1) + 1) := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hl1]
    have e2 : (i - 1) * l + l = i * l := by
      calc (i - 1) * l + l = ((i - 1) + 1) * l := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hi1]
    have hidx : i * (l - 1) = (i - 1) * l + (l - i) := by omega
    rw [← pow_mul, hidx, pow_add, mul_comm (i - 1) l, pow_mul, hζl, one_pow, one_mul]
  have hterm : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
        = tateDXterm (ζ ^ i) + tateDXtail (ζ ^ i) q hq
          - tateDXtail (ζ ^ (l - i)) q hq := by
    intro i hi
    rw [tateDXpair, hpow i hi, ← tateDXtail_rec]
  have hswap : ∑ i ∈ (range l).erase 0, tateDXtail (ζ ^ (l - i)) q hq
      = ∑ i ∈ (range l).erase 0, tateDXtail (ζ ^ i) q hq :=
    sum_erase_reflect (fun j => tateDXtail (ζ ^ j) q hq)
  rw [Finset.sum_congr rfl hterm, Finset.sum_sub_distrib, Finset.sum_add_distrib, hswap]
  have hS : ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ i) = 0 := by
    have hanti : ∀ i ∈ (range l).erase 0,
        tateDXterm (ζ ^ (l - i)) = - tateDXterm (ζ ^ i) := by
      intro i hi
      have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
      have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
      have hmem' : l - i ∈ (range l).erase 0 :=
        Finset.mem_erase.2 ⟨by omega, Finset.mem_range.2 (by omega)⟩
      have huv : ζ ^ i * ζ ^ (l - i) = 1 := by
        rw [← pow_add, show i + (l - i) = l by omega, hζl]
      exact tateDXterm_inv huv (hu i hi) (hu (l - i) hmem')
    have hr1 : ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ (l - i))
        = ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ i) :=
      sum_erase_reflect (fun j => tateDXterm (ζ ^ j))
    have hr2 : ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ (l - i))
        = - ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ i) := by
      rw [Finset.sum_congr rfl hanti, Finset.sum_neg_distrib]
    have hzero : (2 : R) * ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ i) = 0 := by
      linear_combination hr2 - hr1
    exact (mul_eq_zero.1 hzero).resolve_left h2
  rw [hS]
  ring

/-- **[GenEll] 葉 1 の残り (2)**——`∑_{i≠0} D²X(ζ^i)` の閉じた式。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★`D²X` の `q^N` 係数は `∑_{d∣N}d³(u^d + u^{-d})` なので

    `∑_{i≠0} D²X|_{q^N} = ∑_{d∣N} d³·2(l[l∣d] − 1) = 2(l⁴σ₃(N/l)[l∣N] − σ₃(N))`

☆定数項は `sum_mu_d2xterm_field`（第 852、証明済み）で `120·∑ = l⁴ − 1`。 -/
theorem sum_mu_d2xpair {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    120 * ∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = ((l : R) ^ 4 - 1)
        + 240 * ((l : R) ^ 4 * evalAdic (sigmaSeries 3) (q ^ l) hql
            - evalAdic (sigmaSeries 3) q hq) := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hpow : ∀ i ∈ (range l).erase 0, (ζ ^ i) ^ (l - 1) = ζ ^ (l - i) := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hl1 : 1 ≤ l := hl.one_lt.le
    have hi1 : 1 ≤ i := Nat.one_le_iff_ne_zero.2 hi0
    have e1 : i * (l - 1) + i = i * l := by
      calc i * (l - 1) + i = i * ((l - 1) + 1) := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hl1]
    have e2 : (i - 1) * l + l = i * l := by
      calc (i - 1) * l + l = ((i - 1) + 1) * l := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hi1]
    have hidx : i * (l - 1) = (i - 1) * l + (l - i) := by omega
    rw [← pow_mul, hidx, pow_add, mul_comm (i - 1) l, pow_mul, hζl, one_pow, one_mul]
  have hterm : ∀ i ∈ (range l).erase 0,
      tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
        = tateD2Xterm (ζ ^ i) + tateD2Xtail (ζ ^ i) q hq
          + tateD2Xtail (ζ ^ (l - i)) q hq := by
    intro i hi
    rw [tateD2Xpair, hpow i hi, ← tateD2Xtail_rec]
  have hswap : ∑ i ∈ (range l).erase 0, tateD2Xtail (ζ ^ (l - i)) q hq
      = ∑ i ∈ (range l).erase 0, tateD2Xtail (ζ ^ i) q hq :=
    sum_erase_reflect (fun j => tateD2Xtail (ζ ^ j) q hq)
  rw [Finset.sum_congr rfl hterm, Finset.sum_add_distrib, Finset.sum_add_distrib, hswap,
    sum_mu_d2xtail_sigma hl hζ q hq hql]
  have hconst := sum_mu_d2xterm hl hζ hu
  linear_combination hconst

/-! ## ★★★★★★★★`D²X` と `D³X`——c₆ 側へ -/

/-- ★★★★★★**`D²X = 6X² + X + 2a₄`**——ODE の書き換え。 -/
theorem tate_d2x {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    (a w q : R) (hq : q ∈ I) (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hDX : tateDXpair a w q hq ≠ 0) :
    tateD2Xpair a w q hq
      = 6 * tateXpair a w q hq ^ 2 + tateXpair a w q hq
        + 2 * (tateCurveAt q hq).a₄ := by
  have h := tate_d2x_mul a w q hq haw ha hw
  rcases mul_eq_zero.1 h with h1 | h2
  · exact absurd h1 hDX
  · exact sub_eq_zero.1 h2

/-- ★★★★★★★★**`D³X = 12X·DX + DX = 24XY + 12X² + 2Y + X`**。 -/
theorem tate_d3x {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    (a w q : R) (hq : q ∈ I) (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hDX : tateDXpair a w q hq ≠ 0) :
    tateD3Xpair a w q hq
      = 12 * tateXpair a w q hq * tateDXpair a w q hq + tateDXpair a w q hq := by
  have h := tate_d3x_mul a w q hq haw ha hw
  have h2 := tate_d2x a w q hq haw ha hw hDX
  have hzero : tateD2Xpair a w q hq
      * (tateD2Xpair a w q hq - 6 * tateXpair a w q hq ^ 2 - tateXpair a w q hq
        - 2 * (tateCurveAt q hq).a₄) = 0 := by
    rw [h2]
    ring
  rw [hzero, add_zero] at h
  rcases mul_eq_zero.1 h with h1 | h3
  · exact absurd h1 hDX
  · linear_combination h3

/-- ★★★★★★**`∑_{i≠0} D³X(ζ^i) = 0`**——`DX` と同じ理由（奇数回の微分）。 -/
theorem sum_mu_d3xpair_zero {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (h2 : (2 : R) ≠ 0) :
    ∑ i ∈ (range l).erase 0, tateD3Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq = 0 := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hpow : ∀ i ∈ (range l).erase 0, (ζ ^ i) ^ (l - 1) = ζ ^ (l - i) := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hl1 : 1 ≤ l := hl.one_lt.le
    have hi1 : 1 ≤ i := Nat.one_le_iff_ne_zero.2 hi0
    have e1 : i * (l - 1) + i = i * l := by
      calc i * (l - 1) + i = i * ((l - 1) + 1) := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hl1]
    have e2 : (i - 1) * l + l = i * l := by
      calc (i - 1) * l + l = ((i - 1) + 1) * l := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hi1]
    have hidx : i * (l - 1) = (i - 1) * l + (l - i) := by omega
    rw [← pow_mul, hidx, pow_add, mul_comm (i - 1) l, pow_mul, hζl, one_pow, one_mul]
  have hterm : ∀ i ∈ (range l).erase 0,
      tateD3Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
        = tateD3Xterm (ζ ^ i) + tateD3Xtail (ζ ^ i) q hq
          - tateD3Xtail (ζ ^ (l - i)) q hq := by
    intro i hi
    rw [tateD3Xpair, hpow i hi, ← tateD3Xtail_rec]
  have hswap : ∑ i ∈ (range l).erase 0, tateD3Xtail (ζ ^ (l - i)) q hq
      = ∑ i ∈ (range l).erase 0, tateD3Xtail (ζ ^ i) q hq :=
    sum_erase_reflect (fun j => tateD3Xtail (ζ ^ j) q hq)
  rw [Finset.sum_congr rfl hterm, Finset.sum_sub_distrib, Finset.sum_add_distrib, hswap]
  have hS : ∑ i ∈ (range l).erase 0, tateD3Xterm (ζ ^ i) = 0 := by
    have hanti : ∀ i ∈ (range l).erase 0,
        tateD3Xterm (ζ ^ (l - i)) = - tateD3Xterm (ζ ^ i) := by
      intro i hi
      have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
      have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
      have hmem' : l - i ∈ (range l).erase 0 :=
        Finset.mem_erase.2 ⟨by omega, Finset.mem_range.2 (by omega)⟩
      have huv : ζ ^ i * ζ ^ (l - i) = 1 := by
        rw [← pow_add, show i + (l - i) = l by omega, hζl]
      exact tateD3Xterm_inv huv (hu i hi) (hu (l - i) hmem')
    have hr1 : ∑ i ∈ (range l).erase 0, tateD3Xterm (ζ ^ (l - i))
        = ∑ i ∈ (range l).erase 0, tateD3Xterm (ζ ^ i) :=
      sum_erase_reflect (fun j => tateD3Xterm (ζ ^ j))
    have hr2 : ∑ i ∈ (range l).erase 0, tateD3Xterm (ζ ^ (l - i))
        = - ∑ i ∈ (range l).erase 0, tateD3Xterm (ζ ^ i) := by
      rw [Finset.sum_congr rfl hanti, Finset.sum_neg_distrib]
    have hzero : (2 : R) * ∑ i ∈ (range l).erase 0, tateD3Xterm (ζ ^ i) = 0 := by
      linear_combination hr2 - hr1
    exact (mul_eq_zero.1 hzero).resolve_left h2
  rw [hS]
  ring

/-! ## ★★★★★★★★`D⁴X`——σ₅ が出る段 -/

/-- ★★★★★★★★**`D⁴X = 12X·D²X + 12(DX)² + D²X`**。 -/
theorem tate_d4x {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    (a w q : R) (hq : q ∈ I) (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hDX : tateDXpair a w q hq ≠ 0) :
    tateD4Xpair a w q hq
      = 12 * tateXpair a w q hq * tateD2Xpair a w q hq
        + 12 * tateDXpair a w q hq ^ 2 + tateD2Xpair a w q hq := by
  have h := tate_d4x_mul a w q hq haw ha hw
  have h2 := tate_d2x a w q hq haw ha hw hDX
  have h3 := tate_d3x a w q hq haw ha hw hDX
  have hz2 : tateD3Xpair a w q hq
      * (tateD2Xpair a w q hq - 6 * tateXpair a w q hq ^ 2 - tateXpair a w q hq
        - 2 * (tateCurveAt q hq).a₄) = 0 := by
    rw [h2]
    ring
  have hz3 : 2 * tateD2Xpair a w q hq
      * (tateD3Xpair a w q hq - 12 * tateXpair a w q hq * tateDXpair a w q hq
        - tateDXpair a w q hq) = 0 := by
    rw [h3]
    ring
  rw [hz2, hz3, add_zero, add_zero] at h
  rcases mul_eq_zero.1 h with h1 | h4
  · exact absurd h1 hDX
  · linear_combination h4

/-- ★★★★★★★★**`252·∑_{i≠0} D⁴X(ζ^i) = (1 − l⁶) + 504(l⁶s₅(q^l) − s₅(q))`**。 -/
theorem sum_mu_d4xpair {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    252 * ∑ i ∈ (range l).erase 0, tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = (1 - (l : R) ^ 6)
        + 504 * ((l : R) ^ 6 * evalAdic (sigmaSeries 5) (q ^ l) hql
            - evalAdic (sigmaSeries 5) q hq) := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hpow : ∀ i ∈ (range l).erase 0, (ζ ^ i) ^ (l - 1) = ζ ^ (l - i) := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hl1 : 1 ≤ l := hl.one_lt.le
    have hi1 : 1 ≤ i := Nat.one_le_iff_ne_zero.2 hi0
    have e1 : i * (l - 1) + i = i * l := by
      calc i * (l - 1) + i = i * ((l - 1) + 1) := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hl1]
    have e2 : (i - 1) * l + l = i * l := by
      calc (i - 1) * l + l = ((i - 1) + 1) * l := by ring
        _ = i * l := by rw [Nat.sub_add_cancel hi1]
    have hidx : i * (l - 1) = (i - 1) * l + (l - i) := by omega
    rw [← pow_mul, hidx, pow_add, mul_comm (i - 1) l, pow_mul, hζl, one_pow, one_mul]
  have hterm : ∀ i ∈ (range l).erase 0,
      tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
        = tateD4Xterm (ζ ^ i) + tateD4Xtail (ζ ^ i) q hq
          + tateD4Xtail (ζ ^ (l - i)) q hq := by
    intro i hi
    rw [tateD4Xpair, hpow i hi, ← tateD4Xtail_rec]
  have hswap : ∑ i ∈ (range l).erase 0, tateD4Xtail (ζ ^ (l - i)) q hq
      = ∑ i ∈ (range l).erase 0, tateD4Xtail (ζ ^ i) q hq :=
    sum_erase_reflect (fun j => tateD4Xtail (ζ ^ j) q hq)
  rw [Finset.sum_congr rfl hterm, Finset.sum_add_distrib, Finset.sum_add_distrib, hswap,
    sum_mu_d4xtail_sigma hl hζ q hq hql]
  have hconst := sum_mu_d4xterm hl hζ hu
  linear_combination hconst

/-! ## ★★★★★★★★★★★★★★★★★★★★`∑_ζ X²` の閉じた式——Besge の壁は無い -/

/-- ★★★★★★★★★★★★★★★★★★★★**`∑_{i≠0} X(ζ^i)²` の閉じた式**。

★★★これが第 822 以来「**Besge の畳み込み恒等式が要る**」と思われていた量である。
`tate_d2x`（`D²X = 6X² + X + 2a₄`）を足すだけで出る:

    `6·∑X² = ∑D²X − ∑X − 2(l−1)a₄`

であり、`120·∑D²X = (l⁴−1) + 240(l⁴s₃(q^l) − s₃(q))`（第 859）を入れれば良い。
☆**畳み込みは一切現れない**。 -/
theorem sum_mu_X_sq {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    720 * ∑ i ∈ (range l).erase 0, tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 2
      = ((l : R) ^ 4 - 1)
        + 240 * ((l : R) ^ 4 * evalAdic (sigmaSeries 3) (q ^ l) hql
            - evalAdic (sigmaSeries 3) q hq)
        - 120 * (∑ i ∈ (range l).erase 0,
            tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - 240 * ((l : R) - 1) * (tateCurveAt q hq).a₄ := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hawi : ∀ i ∈ (range l).erase 0,
      (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    intro i hi
    have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      rw [← pow_succ']
      rw [Nat.sub_add_cancel hl.pos, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hwu : ∀ i : ℕ, IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) := fun i =>
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  have hterm : ∀ i ∈ (range l).erase 0,
      6 * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 2
        = tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - 2 * (tateCurveAt q hq).a₄ := by
    intro i hi
    have hd2 := tate_d2x (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hu i hi)
      (hwu i) (hDX i hi)
    linear_combination -hd2
  have hsum : 6 * ∑ i ∈ (range l).erase 0,
      tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 2
      = (∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - (∑ i ∈ (range l).erase 0, tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - ((l : R) - 1) * (2 * (tateCurveAt q hq).a₄) := by
    rw [Finset.mul_sum, Finset.sum_congr rfl hterm, Finset.sum_sub_distrib,
      Finset.sum_sub_distrib, Finset.sum_const, Finset.card_erase_of_mem
        (Finset.mem_range.2 hl.pos), Finset.card_range, nsmul_eq_mul]
    have hcard : ((l - 1 : ℕ) : R) = (l : R) - 1 := by
      rw [Nat.cast_sub hl.one_lt.le, Nat.cast_one]
    rw [hcard]
  have hd2sum := sum_mu_d2xpair hl hζ hu q hq hql
  linear_combination (120 : R) * hsum + hd2sum

/-- ★★★★★★★★**`∑_{i≠0} X(ζ^i)·Y(ζ^i)`**——`∑D³X = 0` から出る。

`D³X = 24XY + 12X² + 2Y + X` と `∑D³X = 0` で

    `24·∑XY = −12·∑X² − 2·∑Y − ∑X`

☆**`∑D³X` の値を知らなくても出る**のがここの味味である。 -/
theorem sum_mu_XY {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    24 * ∑ i ∈ (range l).erase 0,
        (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          * tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      = -(12 * ∑ i ∈ (range l).erase 0,
            tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 2)
        - 2 * (∑ i ∈ (range l).erase 0,
            tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - (∑ i ∈ (range l).erase 0,
            tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq) := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hawi : ∀ i ∈ (range l).erase 0,
      (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    intro i hi
    have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      rw [← pow_succ']
      rw [Nat.sub_add_cancel hl.pos, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hwu : ∀ i : ℕ, IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) := fun i =>
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  have hterm : ∀ i ∈ (range l).erase 0,
      24 * (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
              * tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        = tateD3Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - 12 * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 2
          - 2 * tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    intro i hi
    have hd3 := tate_d3x (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hu i hi)
      (hwu i) (hDX i hi)
    have hdx := tateDXpair_eq (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hu i hi) (hwu i)
    rw [hd3, hdx]
    ring
  have hz3 := sum_mu_d3xpair_zero hl hζ hu q hq h2
  rw [Finset.mul_sum, Finset.sum_congr rfl hterm, Finset.sum_sub_distrib,
    Finset.sum_sub_distrib, Finset.sum_sub_distrib, hz3, ← Finset.mul_sum,
    ← Finset.mul_sum]
  ring

/-- ★★★★★★★★**`∑_{i≠0} X(ζ^i)³`**——`D⁴X` から出る。

`D⁴X = 120X³ + 30X² + 72a₄X + X + 48a₆ + 2a₄`（`tate_d4x` と `tate_equation` から）
なので

    `120·∑X³ = ∑D⁴X − 30∑X² − (72a₄ + 1)∑X − (l−1)(48a₆ + 2a₄)`

☆`∑D⁴X` は `sum_mu_d4xpair`（第 866）で閉じた式になっている。 -/
theorem sum_mu_X_cube {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    120 * ∑ i ∈ (range l).erase 0,
        tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 3
      = (∑ i ∈ (range l).erase 0, tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - 30 * (∑ i ∈ (range l).erase 0,
            tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 2)
        - (72 * (tateCurveAt q hq).a₄ + 1)
          * (∑ i ∈ (range l).erase 0, tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - ((l : R) - 1) * (48 * (tateCurveAt q hq).a₆ + 2 * (tateCurveAt q hq).a₄) := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hawi : ∀ i ∈ (range l).erase 0,
      (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    intro i hi
    have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      rw [← pow_succ']
      rw [Nat.sub_add_cancel hl.pos, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hwu : ∀ i : ℕ, IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) := fun i =>
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  have hterm : ∀ i ∈ (range l).erase 0,
      120 * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 3
        = tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - 30 * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ^ 2
          - (72 * (tateCurveAt q hq).a₄ + 1)
            * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - (48 * (tateCurveAt q hq).a₆ + 2 * (tateCurveAt q hq).a₄) := by
    intro i hi
    have hd4 := tate_d4x (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hu i hi)
      (hwu i) (hDX i hi)
    have hd2 := tate_d2x (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hu i hi)
      (hwu i) (hDX i hi)
    have hdx := tateDXpair_eq (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hu i hi) (hwu i)
    have heq := tate_equation (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi)
      (hu i hi) (hwu i)
    rw [hd4, hd2, hdx]
    linear_combination (-48 : R) * heq
  rw [Finset.mul_sum, Finset.sum_congr rfl hterm, Finset.sum_sub_distrib,
    Finset.sum_sub_distrib, Finset.sum_sub_distrib, Finset.sum_const,
    Finset.card_erase_of_mem (Finset.mem_range.2 hl.pos), Finset.card_range,
    nsmul_eq_mul, ← Finset.mul_sum, ← Finset.mul_sum]
  have hcard : ((l - 1 : ℕ) : R) = (l : R) - 1 := by
    rw [Nat.cast_sub hl.one_lt.le, Nat.cast_one]
  rw [hcard]

/-! ## ★出典の紐付け(`.src`)と証明義務(`.needs`) -/

def tate_ode.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の ODE——DY = 3X² − Y + a₄)",
    sectionId := "genell-lemma-3-2" }

def tate_ode.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_ode_mul(双対数による ODE、第 851、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_ode_mul") 1,
    .implicitStep "☆DX ≠ 0 で割る段(R は整域)" 1 ]

def veluV2_eq_tateDYpair.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(veluV2 は DY そのもの。★ODE から直ちに)",
    sectionId := "genell-lemma-3-2" }

def veluV2_eq_tateDYpair.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluV2_tateCurveAt(Tate 曲線では v_Q = 3x² + a₄ − y、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluV2_tateCurveAt") 1,
    .citation "[ABC3]" "tate_ode(本ファイルの ODE)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tate_ode") 1 ]

def sum_veluV2_eq_sum_tateDYpair.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(v = ∑_{i≠0} DY(ζ^i))",
    sectionId := "genell-lemma-3-2" }

def sum_veluV2_eq_sum_tateDYpair.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluV2_eq_tateDYpair(各項で veluV2 = DY)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.veluV2_eq_tateDYpair") 1,
    .implicitStep "☆ζ^i · (q·(ζ^i)^{l-1}) = q を ζ^l = 1 から出す段" 1 ]

def sum_mu_dyterm.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(定数項の指標和——240·∑ ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_dyterm.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_mu_dyterm_field(第 849、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_dyterm_field") 1 ]

def sum_mu_dxpair_zero.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{i≠0} DX(ζ^i) = 0)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_dxpair_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateDXterm_inv(Df(1/t) = −Df(t)、第 854、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDXterm_inv") 1,
    .citation "[ABC3]" "tateDXtail_rec(尾の漸化式、第 854、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDXtail_rec") 1,
    .citation "[ABC3]" "sum_erase_reflect(i ↦ l−i の全単射性、第 854、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_erase_reflect") 1 ]

def sum_mu_d2xpair.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{i≠0} D²X(ζ^i) の閉じた式)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_d2xpair.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_mu_d2xterm(定数項 120·∑ = l⁴−1、第 855、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_d2xterm") 1,
    .citation "[ABC3]" "sum_mu_d2xtail_sigma(尾 = l⁴s₃(q^l) − s₃(q)、第 858、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_d2xtail_sigma") 1,
    .citation "[ABC3]" "tateD2Xtail_rec(尾の漸化式、第 855、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateD2Xtail_rec") 1 ]

def tate_d2x.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D²X = 6X² + X + 2a₄)",
    sectionId := "genell-lemma-3-2" }

def tate_d2x.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d2x_mul(掛けた形、第 861、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_d2x_mul") 1 ]

def tate_d3x.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D³X = 12X·DX + DX)",
    sectionId := "genell-lemma-3-2" }

def tate_d3x.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d3x_mul(双対数で tate_d2x_mul を微分したもの、第 861)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_d3x_mul") 1 ]

def sum_mu_d3xpair_zero.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{i≠0} D³X(ζ^i) = 0)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_d3xpair_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateD3Xterm_inv(D³f(1/t) = −D³f(t)、第 860、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateD3Xterm_inv") 1 ]

def tate_d4x.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D⁴X = 12X·D²X + 12(DX)² + D²X)",
    sectionId := "genell-lemma-3-2" }

def tate_d4x.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d4x_mul(双対数で tate_d3x_mul を微分したもの、第 863)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_d4x_mul") 1 ]

def sum_mu_d4xpair.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(252·∑ D⁴X = (1−l⁶) + 504(l⁶s₅(q^l) − s₅(q)))",
    sectionId := "genell-lemma-3-2" }

def sum_mu_d4xpair.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_mu_d4xterm(定数項 252·∑ = 1−l⁶、第 865、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_d4xterm") 1,
    .citation "[ABC3]" "sum_mu_d4xtail_sigma(尾 = l⁶s₅(q^l) − s₅(q)、第 864、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_d4xtail_sigma") 1 ]

def sum_mu_X_sq.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{i≠0} X(ζ^i)² の閉じた式——Besge は要らない)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_X_sq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d2x(D²X = 6X² + X + 2a₄、第 861、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tate_d2x") 1,
    .citation "[ABC3]" "sum_mu_d2xpair(120·∑D²X の閉じた式、第 859、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_mu_d2xpair") 1 ]

def sum_mu_XY.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(24·∑XY = −12∑X² − 2∑Y − ∑X)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_XY.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d3x(D³X = 12X·DX + DX、第 861、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tate_d3x") 1,
    .citation "[ABC3]" "sum_mu_d3xpair_zero(∑D³X = 0、第 861、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_mu_d3xpair_zero") 1 ]

def sum_mu_X_cube.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(120·∑X³ の閉じた式)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_X_cube.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d4x(D⁴X = 12X·D²X + 12(DX)² + D²X、第 866、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tate_d4x") 1,
    .citation "[ABC3]" "tate_equation(Y² + XY = X³ + a₄X + a₆、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_equation") 1 ]

end ABC3.Skeleton.GenEll
