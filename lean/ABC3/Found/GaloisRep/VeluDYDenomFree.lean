/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.PowerSeriesUniversal
import ABC3.Found.GaloisRep.MuPairFunctorial
import ABC3.Skeleton.GenEll.TateODE
import ABC3.Meta.Claim

/-!
# 第 1128 ブロック —— **`veluV2 = DY` の分母なし版**（`Found`、§3 の枠の節点 3）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か

`veluV2_eq_tateDYpair`（`Skeleton/GenEll/TateODE.lean:89`）は
**`IsUnit (1 − a)` を要求する**——`X`・`Y` が `Ring.inverse (1 − a)` で書かれているからである。
`p ∣ l` の悪い素点では `1 − ζ^i` が単元でないので、この形は使えない。

★本ブロックは**分母を払った形**

    `veluV2DF l E_q (l²X) (l³Y) = l² · (l⁴·DY)`

を、`IsUnit (1 − a)` も `IsUnit (l)` も**要求せずに**証明する。

## ★★★★★★★★道（第 1127 の設計どおり）

| 段 | 場所 | 何をするか |
|---|---|---|
| 1 | `A₁ = PowerSeries K`（`K` は体） | `1 − C α` も `(l)` も単元なので、**在庫の `hu` つきの補題がそのまま使える** |
| 2 | `A₀ = PowerSeries A → A₁ = PowerSeries (FractionRing A)` | `PowerSeries.map` の単射性で `A₀` に降ろす（第 1125 の `map_*DF` が通り道） |
| 3 | `A₀ → R`（`X ↦ q`） | `evalAdicMapHom`（第 1126）で特殊化する |

☆`PowerSeries` を経由するのは、**係数環を取り替えても `(X)`-adic 完備性が壊れない**からである
（第 1091 で「商体に移すと `q` の収束が壊れる」と測った行き止まりの抜け道）。

## ★★★★★★側条件 `hDX ≠ 0` はどこへ行ったか

`veluV2_eq_tateDYpair` は `tateDXpair ≠ 0` を仮説に取る。
★`A₁` では**定数項を見るだけで出る**（`tateDXpair_C_ne_zero`）——
定数項は `α(1+α)/(1−α)³` で、`α ≠ 0`・`α ≠ −1`・`α ≠ 1` なら `0` でない。
☆尾はすべて `(X)` に入るので定数項に効かない。
★`α ≠ −1` は `l` が奇素数であることから出る（`Lemma 3.5` はすでに `l ≠ 2` を置いている）。
-/

namespace ABC3.Found.GaloisRep

open PowerSeries Finset

/-! ## ★定数項 -/

/-- ☆`(X)` の元の定数項は `0`（`constantCoeff` の形）。 -/
theorem constantCoeff_of_mem_span_X {A : Type} [CommRing A] {f : PowerSeries A}
    (hf : f ∈ Ideal.span {(PowerSeries.X : PowerSeries A)}) :
    PowerSeries.constantCoeff f = 0 := by
  simpa using coeff_zero_of_mem_span_X hf

/-! ## ★★★★★★段 1 —— 体の上（`A₁`）では側条件が定数項で出る -/

/-- ★★★★★★★★**`DX ≠ 0` は定数項で出る**（第 1128）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

☆`tateDXpair (C α) w X` の定数項は `α(1+α)·(1−α)⁻³` である——
尾（`tateDXtail`）も `w` 側（`w ∈ (X)`）もすべて `(X)` に入るので効かない。
★`TateDXNeZero.lean` の `tateDXpair_ne_zero_of_mu` は `TateSetup` と Tate 一意化を要するが、
`A₁` では**その重装備は要らない**。 -/
theorem tateDXpair_C_ne_zero {K : Type} [Field K] {α : K}
    (hα0 : α ≠ 0) (hαneg : α + 1 ≠ 0) (hα1 : α ≠ 1)
    {w : PowerSeries K} (hw : w ∈ Ideal.span {(PowerSeries.X : PowerSeries K)}) :
    tateDXpair (PowerSeries.C α) w PowerSeries.X (Ideal.mem_span_singleton_self _) ≠ 0 := by
  have hu : IsUnit (1 - α) := isUnit_iff_ne_zero.2 (sub_ne_zero.2 (Ne.symm hα1))
  intro h0
  have hcc := congrArg (PowerSeries.constantCoeff (R := K)) h0
  rw [tateDXpair, map_sub, map_add, map_add] at hcc
  have h1 : PowerSeries.constantCoeff (tateDXterm (PowerSeries.C α : PowerSeries K))
      = tateDXterm α := by
    rw [← map_tateDXterm (PowerSeries.C (R := K)) hu]; simp
  have h2 : PowerSeries.constantCoeff
      (tateDXtail (PowerSeries.C α : PowerSeries K) PowerSeries.X
        (Ideal.mem_span_singleton_self _)) = 0 :=
    constantCoeff_of_mem_span_X (tateDXtail_mem _ _ _)
  have h3mem : tateDXterm w ∈ Ideal.span {(PowerSeries.X : PowerSeries K)} := by
    have hh := tateDXterm_mem_pow (I := Ideal.span {(PowerSeries.X : PowerSeries K)})
      (k := 1) (t := w) (by simpa using hw)
    simpa using hh
  have h3 : PowerSeries.constantCoeff (tateDXterm w) = 0 := constantCoeff_of_mem_span_X h3mem
  have h4 : PowerSeries.constantCoeff
      (tateDXtail w PowerSeries.X (Ideal.mem_span_singleton_self _)) = 0 :=
    constantCoeff_of_mem_span_X (tateDXtail_mem _ _ _)
  rw [h1, h2, h3, h4] at hcc
  simp only [add_zero, sub_zero, map_zero] at hcc
  have hinv : Ring.inverse (1 - α) ≠ 0 := by
    rw [Ring.inverse_eq_inv]
    exact inv_ne_zero (sub_ne_zero.2 (Ne.symm hα1))
  rw [tateDXterm] at hcc
  have hne : α * (1 + α) * Ring.inverse (1 - α) ^ 3 ≠ 0 := by
    refine mul_ne_zero (mul_ne_zero hα0 ?_) (pow_ne_zero _ hinv)
    rw [add_comm]; exact hαneg
  exact hne hcc

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★**段 1**——体 `K` の上の `PowerSeries K` での分母なし恒等式。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆ここでは `1 − C α` も `(l)` も単元なので、在庫の `veluV2_eq_tateDYpair` と
第 1118-1119 の橋（`natCast_pow_mul_veluV2_tate`・`natCast_pow_mul_tateDYpair`）が
そのまま使える。★結論の形には `Ring.inverse` が現れない。 -/
theorem veluV2DF_eq_tateDYpairDF_field {K : Type} [Field K] {l : ℕ} (hl0 : (l : K) ≠ 0)
    {α : K} (hα1 : α ≠ 1) (hα0 : α ≠ 0) (hαneg : α + 1 ≠ 0)
    (hpow : α ^ l = 1) (hsum : ∑ k ∈ range l, α ^ k = 0) :
    veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries K) (Ideal.mem_span_singleton_self _))
        (tateXpairDF l (PowerSeries.C α)
          ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X) PowerSeries.X
          (Ideal.mem_span_singleton_self _))
        (tateYpairDF l (PowerSeries.C α)
          ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X) PowerSeries.X
          (Ideal.mem_span_singleton_self _))
      = (l : PowerSeries K) ^ 2 * tateDYpairDF l (PowerSeries.C α)
          ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X) PowerSeries.X
          (Ideal.mem_span_singleton_self _) := by
  have hq : (PowerSeries.X : PowerSeries K)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries K)} := Ideal.mem_span_singleton_self _
  have hw : ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries K)} := Ideal.mul_mem_left _ _ hq
  have hlpos : 0 < l := Nat.pos_of_ne_zero (fun h => hl0 (by simp [h]))
  have hpowPS : (PowerSeries.C α : PowerSeries K) ^ l = 1 := by
    rw [← map_pow, hpow, map_one]
  have haw : (PowerSeries.C α : PowerSeries K)
      * ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X) = PowerSeries.X := by
    rw [← mul_assoc, ← pow_succ', Nat.sub_add_cancel hlpos, hpowPS, one_mul]
  have hu : IsUnit (1 - (PowerSeries.C α : PowerSeries K)) := isUnit_one_sub_C hα1
  have hwu : IsUnit (1 - ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries K)}) hw
  have hlu : IsUnit ((l : PowerSeries K)) := isUnit_natCast_powerSeries hl0
  have hsumPS : ∑ k ∈ range l, (PowerSeries.C α : PowerSeries K) ^ k = 0 := by
    have hcong : ∑ k ∈ range l, (PowerSeries.C α : PowerSeries K) ^ k
        = PowerSeries.C (∑ k ∈ range l, α ^ k) := by
      rw [map_sum]
      exact Finset.sum_congr rfl (fun k _ => (map_pow _ _ _).symm)
    rw [hcong, hsum, map_zero]
  have hDX := tateDXpair_C_ne_zero hα0 hαneg hα1 hw
  have hode := ABC3.Skeleton.GenEll.veluV2_eq_tateDYpair
    (PowerSeries.C α) ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X)
    PowerSeries.X hq haw hu hwu hDX
  have h6 := natCast_pow_mul_veluV2_tate hlu hu hpowPS hsumPS
    (tateCurveAt (PowerSeries.X : PowerSeries K) hq)
    ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X) PowerSeries.X hq
  have h4 := natCast_pow_mul_tateDYpair hlu hu hpowPS hsumPS
    ((PowerSeries.C α : PowerSeries K) ^ (l - 1) * PowerSeries.X) PowerSeries.X hq
  rw [← h6, hode, ← h4]
  ring

/-! ## ★★★★★★★★段 2 —— 整域の `PowerSeries` へ降ろす -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★**段 2**——`A` が整域なら `PowerSeries A` でも成り立つ。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`A → FractionRing A` は単射だから `PowerSeries.map` も単射（`PowerSeries.map_injective`）。
★第 1125 の `map_*DF` は `a` 側の単元性を要らないので、**この降下が通る**。 -/
theorem veluV2DF_eq_tateDYpairDF_powerSeries {A : Type} [CommRing A] [IsDomain A] {l : ℕ}
    (hl0 : (l : A) ≠ 0) {α : A} (hα1 : α ≠ 1) (hα0 : α ≠ 0) (hαneg : α + 1 ≠ 0)
    (hpow : α ^ l = 1) (hsum : ∑ k ∈ range l, α ^ k = 0) :
    veluV2DF l (tateCurveAt (PowerSeries.X : PowerSeries A) (Ideal.mem_span_singleton_self _))
        (tateXpairDF l (PowerSeries.C α)
          ((PowerSeries.C α : PowerSeries A) ^ (l - 1) * PowerSeries.X) PowerSeries.X
          (Ideal.mem_span_singleton_self _))
        (tateYpairDF l (PowerSeries.C α)
          ((PowerSeries.C α : PowerSeries A) ^ (l - 1) * PowerSeries.X) PowerSeries.X
          (Ideal.mem_span_singleton_self _))
      = (l : PowerSeries A) ^ 2 * tateDYpairDF l (PowerSeries.C α)
          ((PowerSeries.C α : PowerSeries A) ^ (l - 1) * PowerSeries.X) PowerSeries.X
          (Ideal.mem_span_singleton_self _) := by
  have hψ : Function.Injective (algebraMap A (FractionRing A)) :=
    IsFractionRing.injective A (FractionRing A)
  refine PowerSeries.map_injective (algebraMap A (FractionRing A)) hψ ?_
  have hcont : ∀ (n : ℕ) (f : PowerSeries A),
      f ∈ (Ideal.span {(PowerSeries.X : PowerSeries A)}) ^ n →
      PowerSeries.map (algebraMap A (FractionRing A)) f
        ∈ (Ideal.span {(PowerSeries.X : PowerSeries (FractionRing A))}) ^ n :=
    fun n f hf => map_mem_span_X_pow _ n hf
  have hq : (PowerSeries.X : PowerSeries A)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries A)} := Ideal.mem_span_singleton_self _
  have hw : ((PowerSeries.C α : PowerSeries A) ^ (l - 1) * PowerSeries.X)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries A)} := Ideal.mul_mem_left _ _ hq
  have hwu : IsUnit (1 - ((PowerSeries.C α : PowerSeries A) ^ (l - 1) * PowerSeries.X)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries A)}) hw
  rw [map_veluV2DF, map_tateXpairDF _ hcont _ _ _ _ _ hwu,
    map_tateYpairDF _ hcont _ _ _ _ _ hwu, map_tateCurveAt _ hcont]
  conv_rhs => rw [map_mul, map_pow, map_natCast, map_tateDYpairDF _ hcont _ _ _ _ _ hwu]
  simp only [PowerSeries.map_C, PowerSeries.map_X, map_mul, map_pow]
  -- ★体の側の仮説は、すべて単射性 `hψ` で `A` の側から出る
  have hψ0 : ∀ x : A, algebraMap A (FractionRing A) x = 0 → x = 0 :=
    fun x h => hψ (by simpa using h)
  have hl0' : ((l : ℕ) : FractionRing A) ≠ 0 := fun h =>
    hl0 (hψ0 _ (by rw [map_natCast]; exact h))
  have hα1' : algebraMap A (FractionRing A) α ≠ 1 := fun h =>
    hα1 (hψ (by rw [h, map_one]))
  have hα0' : algebraMap A (FractionRing A) α ≠ 0 := fun h => hα0 (hψ0 _ h)
  have hαneg' : algebraMap A (FractionRing A) α + 1 ≠ 0 := fun h =>
    hαneg (hψ0 _ (by rw [map_add, map_one]; exact h))
  have hpow' : (algebraMap A (FractionRing A) α) ^ l = 1 := by
    rw [← map_pow, hpow, map_one]
  have hsum' : ∑ k ∈ range l, (algebraMap A (FractionRing A) α) ^ k = 0 := by
    have hcong : ∑ k ∈ range l, (algebraMap A (FractionRing A)) α ^ k
        = algebraMap A (FractionRing A) (∑ k ∈ range l, α ^ k) := by
      rw [map_sum]
      exact Finset.sum_congr rfl (fun k _ => (map_pow _ _ _).symm)
    rw [hcong, hsum, map_zero]
  exact veluV2DF_eq_tateDYpairDF_field hl0' hα1' hα0' hαneg' hpow' hsum'

/-! ## ★★★★★★★★★★★★段 3 —— 完備環 `R` へ特殊化する（★節点 3 の核） -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `veluV2 = DY` の分母なし版**（第 1128、§3 の枠の節点 3 の核）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★**`IsUnit (1 − α)` も `IsUnit (l)` も仮説に置いていない**——
`p ∣ l` の悪い素点でもそのまま意味を持ち、そこで成り立つ。
☆置いているのは `α ≠ 0`・`α ≠ 1`・`α ≠ −1` と `α^l = 1`・`∑ α^k = 0` だけで、
どれも `α = ζ^i`（`ζ` は原始 `l` 乗根、`l` は奇素数、`0 < i < l`）から出る。 -/
theorem veluV2DF_eq_tateDYpairDF {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl0 : (l : R) ≠ 0) {α : R} (hα1 : α ≠ 1) (hα0 : α ≠ 0)
    (hαneg : α + 1 ≠ 0) (hpow : α ^ l = 1) (hsum : ∑ k ∈ range l, α ^ k = 0)
    (q : R) (hq : q ∈ I) :
    veluV2DF l (tateCurveAt q hq) (tateXpairDF l α (α ^ (l - 1) * q) q hq)
        (tateYpairDF l α (α ^ (l - 1) * q) q hq)
      = (l : R) ^ 2 * tateDYpairDF l α (α ^ (l - 1) * q) q hq := by
  have hX : (PowerSeries.X : PowerSeries R)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries R)} := Ideal.mem_span_singleton_self _
  have hwX : ((PowerSeries.C α : PowerSeries R) ^ (l - 1) * PowerSeries.X)
      ∈ Ideal.span {(PowerSeries.X : PowerSeries R)} := Ideal.mul_mem_left _ _ hX
  have hwu : IsUnit (1 - ((PowerSeries.C α : PowerSeries R) ^ (l - 1) * PowerSeries.X)) :=
    isUnit_one_sub (I := Ideal.span {(PowerSeries.X : PowerSeries R)}) hwX
  have hev : ∀ (n : ℕ) (f : PowerSeries R),
      f ∈ (Ideal.span {(PowerSeries.X : PowerSeries R)}) ^ n →
      evalAdicMapHom (RingHom.id R) q hq f ∈ I ^ n :=
    fun n f hf => evalAdicMapHom_mem_pow (RingHom.id R) q hq n f hf
  have hps := veluV2DF_eq_tateDYpairDF_powerSeries (A := R) hl0 hα1 hα0 hαneg hpow hsum
  have hmain := congrArg (evalAdicMapHom (RingHom.id R) q hq) hps
  rw [map_veluV2DF, map_tateXpairDF _ hev _ _ _ _ _ hwu,
    map_tateYpairDF _ hev _ _ _ _ _ hwu, map_tateCurveAt _ hev] at hmain
  conv_rhs at hmain => rw [map_mul, map_pow, map_natCast,
    map_tateDYpairDF _ hev _ _ _ _ _ hwu]
  simp only [evalAdicMapHom_apply, evalAdicMap_C, evalAdicMap_X, RingHom.id_apply,
    map_mul, map_pow] at hmain
  exact hmain

/-! ## ★出典の紐付け(`.src`) -/

def veluV2DF_eq_tateDYpairDF.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(節点 3——veluV2 = DY の E 版。hu も IsUnit(l) も取らない)",
    sectionId := "genell-lemma-3-5" }

def tateDXpair_C_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(PowerSeries の上では DX ≠ 0 が定数項で出ること)",
    sectionId := "genell-def-3-3" }

def veluV2DF_eq_tateDYpairDF.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "veluV2_eq_tateDYpair(hu つきの原型、第 851、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.veluV2_eq_tateDYpair") 1,
    .citation "[ABC3]" "map_tateXpairDF(DF 形は環準同型を通り抜ける、第 1125、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.map_tateXpairDF") 1,
    .citation "[ABC3]" "evalAdicMapHom(S⟦q⟧ →+* R、第 1126、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdicMapHom") 1,
    .citation "[ABC3]" "evalAdicMapHom_mem_pow(特殊化の連続性、第 1127、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdicMapHom_mem_pow") 1,
    .citation "[mathlib]" "PowerSeries.map_injective(単射な係数写像は冪級数でも単射)"
      (.inMathlib "PowerSeries.map_injective") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1128）**——`hu`（`IsUnit (1 − a)`）と " ++
       "`hlu`（`IsUnit (l)`）の両方が結論から消えた。" ++
       "☆道は 3 段: (1) `PowerSeries K`（体）で在庫の `hu` つきの補題を当てる、" ++
       "(2) `PowerSeries.map` の単射性で `PowerSeries A`（整域）に降ろす、" ++
       "(3) `evalAdicMapHom` で完備環 `R` に特殊化する。" ++
       "★側条件 `hDX ≠ 0` は `A₁` の定数項 `α(1+α)/(1−α)³` を見るだけで出た。") 6 ]

end ABC3.Found.GaloisRep
