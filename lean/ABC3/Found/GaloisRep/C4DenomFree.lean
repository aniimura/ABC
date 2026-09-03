/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluDYDenomFree
import ABC3.Skeleton.GenEll.TateIsogenyAny
import ABC3.Meta.Claim

/-!
# 第 1129 ブロック —— **`c₄` の恒等式の分母なし版**（`Found`、§3 の枠の節点 4 の前半）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`c4_velu_tate`（`Skeleton/GenEll/TateIsogeny.lean:154`）は

* `hlu : IsUnit ((l : R))`
* `hu : ∀ i, IsUnit (1 − ζ^i)`
* `hDX : ∀ i, tateDXpair ... ≠ 0`

の 3 つを仮説に取る。`p ∣ l` の悪い素点では最初の 2 つが成り立たない。

★本ブロックは両辺に `l⁶` を掛けた形

    `l⁶·c₄(E_q) + 240·∑_i veluV2DF = l¹⁰·c₄(E_{q^l})`

を、**どれも仮説に取らずに**証明する。☆道は第 1128 とまったく同じ 3 段である。

| 段 | 場所 | 何をするか |
|---|---|---|
| 1 | `PowerSeries K`（`K` は体） | 3 つとも成り立つので `c4_velu_tate` がそのまま効く |
| 2 | `PowerSeries A → PowerSeries (FractionRing A)` | `PowerSeries.map` の単射性で降ろす |
| 3 | `PowerSeries R → R`（`X ↦ q`） | `evalAdicMapHom` で特殊化する |

☆`hDX` は第 1128 の `tateDXpair_C_ne_zero` で出る（定数項 `α(1+α)/(1−α)³`）。
★そこで `α ≠ −1` を使うので **`l ≠ 2`** が要る——`Lemma 3.5` はすでにそれを置いている。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll ABC3.Skeleton.GenEll PowerSeries Finset

/-- ☆体では `η^l = 1`・`η ≠ 1` から `∑_{k<l} η^k = 0`。 -/
theorem sum_pow_eq_zero_of_ne_one' {F : Type} [Field F] {l : ℕ} {η : F}
    (hpow : η ^ l = 1) (hne : η ≠ 1) : ∑ k ∈ range l, η ^ k = 0 := by
  have h := geom_sum_mul η l
  rw [hpow, sub_self] at h
  rcases mul_eq_zero.1 h with h1 | h2
  · exact h1
  · exact absurd (by linear_combination h2) hne

/-! ## ★★★★★★段 1 —— 体の上 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**段 1**——体 `K` の上の `PowerSeries K` での `c₄` の分母なし恒等式。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆`1 − C ζ^i` も `(l)` も単元なので、`c4_velu_tate` がそのまま効く。 -/
theorem c4_velu_tateDF_field {K : Type} [Field K] [CharZero K] {l : ℕ} (hl : l.Prime) {ζ : K} (hζ : IsPrimitiveRoot ζ l) :
    (l : PowerSeries K) ^ 6
        * (tateCurveAt (PowerSeries.X : PowerSeries K)
            (Ideal.mem_span_singleton_self _)).c₄
        + 240 * ∑ i ∈ (range l).erase 0,
            veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries K)
                (Ideal.mem_span_singleton_self _))
              (tateXpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                PowerSeries.X (Ideal.mem_span_singleton_self _))
              (tateYpairDF l ((PowerSeries.C ζ : PowerSeries K) ^ i)
                (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries K) ^ i) ^ (l - 1))
                PowerSeries.X (Ideal.mem_span_singleton_self _))
      = (l : PowerSeries K) ^ 10
          * (tateCurveAt ((PowerSeries.X : PowerSeries K) ^ l)
              (Ideal.pow_mem_of_mem _ (Ideal.mem_span_singleton_self _) l hl.pos)).c₄ := by
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
  have hc4 := c4_velu_tate_any hl hζPS hlu hu PowerSeries.X hX hXl h2
  have hkey : ∑ i ∈ (range l).erase 0,
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
  rw [hkey]
  linear_combination (l : PowerSeries K) ^ 6 * hc4

/-! ## ★★★★★★★★段 2 —— 整域の `PowerSeries` へ降ろす -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**段 2**——`A` が標数 0 の整域なら `PowerSeries A` でも成り立つ。 -/
theorem c4_velu_tateDF_powerSeries {A : Type} [CommRing A] [IsDomain A] [CharZero A] {l : ℕ}
    (hl : l.Prime) {ζ : A} (hζ : IsPrimitiveRoot ζ l) :
    (l : PowerSeries A) ^ 6 * (tateCurveAt (PowerSeries.X : PowerSeries A)
          (Ideal.mem_span_singleton_self _)).c₄
      + 240 * ∑ i ∈ (range l).erase 0,
          veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries A)
              (Ideal.mem_span_singleton_self _))
            (tateXpairDF l ((PowerSeries.C ζ : PowerSeries A) ^ i)
              (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries A) ^ i) ^ (l - 1))
              PowerSeries.X (Ideal.mem_span_singleton_self _))
            (tateYpairDF l ((PowerSeries.C ζ : PowerSeries A) ^ i)
              (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries A) ^ i) ^ (l - 1))
              PowerSeries.X (Ideal.mem_span_singleton_self _))
    = (l : PowerSeries A) ^ 10 * (tateCurveAt ((PowerSeries.X : PowerSeries A) ^ l)
            (Ideal.pow_mem_of_mem _ (Ideal.mem_span_singleton_self _) l hl.pos)).c₄ := by
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
      IsUnit (1 - (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries A) ^ i) ^ (l - 1))) :=
    fun i => isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries A)})
      (Ideal.mul_mem_right _ _ hX)
  rw [map_add, map_mul, map_mul, map_mul, map_sum,
    show ∀ (f : PowerSeries A →+* PowerSeries (FractionRing A)),
      f ((l : PowerSeries A) ^ 6) = (l : PowerSeries (FractionRing A)) ^ 6 from
      fun f => by rw [map_pow, map_natCast],
    show ∀ (f : PowerSeries A →+* PowerSeries (FractionRing A)),
      f ((l : PowerSeries A) ^ 10) = (l : PowerSeries (FractionRing A)) ^ 10 from
      fun f => by rw [map_pow, map_natCast],
    show ∀ (f : PowerSeries A →+* PowerSeries (FractionRing A)),
      f (240 : PowerSeries A) = (240 : PowerSeries (FractionRing A)) from
      fun f => by rw [map_ofNat],
    ← WeierstrassCurve.map_c₄, ← WeierstrassCurve.map_c₄,
    map_tateCurveAt _ hcont, map_tateCurveAt _ hcont]
  have hterm : ∀ i ∈ (range l).erase 0,
      PowerSeries.map (algebraMap A (FractionRing A))
        (veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries A) hX)
          (tateXpairDF l ((PowerSeries.C ζ : PowerSeries A) ^ i)
            (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries A) ^ i) ^ (l - 1))
            PowerSeries.X hX)
          (tateYpairDF l ((PowerSeries.C ζ : PowerSeries A) ^ i)
            (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries A) ^ i) ^ (l - 1))
            PowerSeries.X hX))
      = veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries (FractionRing A))
            (Ideal.mem_span_singleton_self _))
          (tateXpairDF l ((PowerSeries.C (algebraMap A (FractionRing A) ζ)) ^ i)
            (PowerSeries.X * ((PowerSeries.C (algebraMap A (FractionRing A) ζ)) ^ i) ^ (l - 1))
            PowerSeries.X (Ideal.mem_span_singleton_self _))
          (tateYpairDF l ((PowerSeries.C (algebraMap A (FractionRing A) ζ)) ^ i)
            (PowerSeries.X * ((PowerSeries.C (algebraMap A (FractionRing A) ζ)) ^ i) ^ (l - 1))
            PowerSeries.X (Ideal.mem_span_singleton_self _)) := by
    intro i _
    rw [map_veluV2DF, map_tateXpairDF _ hcont _ _ _ _ _ (hwu i),
      map_tateYpairDF _ hcont _ _ _ _ _ (hwu i), map_tateCurveAt _ hcont]
    simp only [PowerSeries.map_C, PowerSeries.map_X, map_mul, map_pow]
  rw [Finset.sum_congr rfl hterm]
  simp only [PowerSeries.map_X, map_pow]
  exact c4_velu_tateDF_field hl (hζ.map_of_injective hψ)

/-! ## ★★★★★★★★★★★★段 3 —— 完備環 `R` へ特殊化する -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `c₄` の恒等式の分母なし版**（第 1129、§3 の枠の節点 4 の前半）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★**`IsUnit (l)` も `IsUnit (1 − ζ^i)` も `hDX` も仮説に置いていない**——
`p ∣ l` の悪い素点でもそのまま意味を持ち、そこで成り立つ。 -/
theorem c4_velu_tateDF {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R}
    (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    (l : R) ^ 6 * (tateCurveAt q hq).c₄
        + 240 * ∑ i ∈ (range l).erase 0,
            veluV2DF l (tateCurveAt q hq)
              (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      = (l : R) ^ 10 * (tateCurveAt (q ^ l) hql).c₄ := by
  have hX : (PowerSeries.X : PowerSeries R)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries R)} := Ideal.mem_span_singleton_self _
  have hwu : ∀ i : ℕ,
      IsUnit (1 - (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries R) ^ i) ^ (l - 1))) :=
    fun i => isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries R)})
      (Ideal.mul_mem_right _ _ hX)
  have hev : ∀ (n : ℕ) (f : PowerSeries R),
      f ∈ (Ideal.span {(PowerSeries.X : PowerSeries R)}) ^ n →
      evalAdicMapHom (RingHom.id R) q hq f ∈ I ^ n :=
    fun n f hf => evalAdicMapHom_mem_pow (RingHom.id R) q hq n f hf
  have hps := c4_velu_tateDF_powerSeries (A := R) hl hζ
  have hmain := congrArg (evalAdicMapHom (RingHom.id R) q hq) hps
  rw [map_add, map_mul, map_mul, map_mul, map_sum,
    show ∀ (f : PowerSeries R →+* R),
      f ((l : PowerSeries R) ^ 6) = (l : R) ^ 6 from
      fun f => by rw [map_pow, map_natCast],
    show ∀ (f : PowerSeries R →+* R),
      f ((l : PowerSeries R) ^ 10) = (l : R) ^ 10 from
      fun f => by rw [map_pow, map_natCast],
    show ∀ (f : PowerSeries R →+* R),
      f (240 : PowerSeries R) = (240 : R) from fun f => by rw [map_ofNat],
    ← WeierstrassCurve.map_c₄, ← WeierstrassCurve.map_c₄,
    map_tateCurveAt _ hev, map_tateCurveAt _ hev] at hmain
  have hterm : ∀ i ∈ (range l).erase 0,
      evalAdicMapHom (RingHom.id R) q hq
        (veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries R) hX)
          (tateXpairDF l ((PowerSeries.C ζ : PowerSeries R) ^ i)
            (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries R) ^ i) ^ (l - 1))
            PowerSeries.X hX)
          (tateYpairDF l ((PowerSeries.C ζ : PowerSeries R) ^ i)
            (PowerSeries.X * ((PowerSeries.C ζ : PowerSeries R) ^ i) ^ (l - 1))
            PowerSeries.X hX))
      = veluV2DF l (tateCurveAt q hq)
          (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq) := by
    intro i _
    rw [map_veluV2DF, map_tateXpairDF _ hev _ _ _ _ _ (hwu i),
      map_tateYpairDF _ hev _ _ _ _ _ (hwu i), map_tateCurveAt _ hev]
    simp only [evalAdicMapHom_apply, evalAdicMap_C, evalAdicMap_X, RingHom.id_apply,
      map_mul, map_pow]
  rw [Finset.sum_congr rfl hterm] at hmain
  simp only [evalAdicMapHom_apply, evalAdicMap_X, map_pow] at hmain
  exact hmain

/-! ## ★出典の紐付け(`.src`) -/

def c4_velu_tateDF.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₄ の恒等式の E 版。hlu も hu も hDX も取らない)",
    sectionId := "genell-lemma-3-2" }

def c4_velu_tateDF.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "c4_velu_tate(hu つきの原型、第 866、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.c4_velu_tate") 1,
    .citation "[ABC3]" "tateDXpair_C_ne_zero(側条件は定数項で出る、第 1128、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDXpair_C_ne_zero") 1,
    .citation "[ABC3]" "evalAdicMapHom_mem_pow(特殊化の連続性、第 1127、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdicMapHom_mem_pow") 1,
    .citation "[mathlib]" "WeierstrassCurve.map_c₄(c₄ は環準同型で移る)"
      (.inMathlib "WeierstrassCurve.map_c₄") 1,
    .implicitStep
      ("★★**2026-09-01（第 1129）**——第 1128 とまったく同じ 3 段（体 → 整域 → 完備環）。" ++
       "☆`hDX` は `PowerSeries K` の定数項 `α(1+α)/(1−α)³` で出る。" ++
       "★そこで `α ≠ −1` を使うので `l ≠ 2` が要る。") 5 ]

end ABC3.Found.GaloisRep
