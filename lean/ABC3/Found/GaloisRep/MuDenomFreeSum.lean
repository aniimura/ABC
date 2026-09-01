/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuHeadDenomFree
import ABC3.Found.GaloisRep.MuDYSum
import ABC3.Found.GaloisRep.MuPowerSum
import ABC3.Found.GaloisRep.MuD2XSum
import ABC3.Meta.Claim

/-!
# 第 1112 ブロック —— **頭項の和を `hu` なしで**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か

`MuDYSum.lean` は既に「**分数体に埋め込む → 体版を使う → 単射性で戻す**」という
転送を実装している（`sum_mu_d2xterm`）。★そこで `hu` が要るのは
`map_tateD2Xterm`（`f (Ring.inverse x) = Ring.inverse (f x)` に単元性が要る）
**ただ 1 箇所**だった。

☆第 1111 で `E` 版が**環準同型と無条件に可換**であることを示したので、
そこを差し替えれば **`hu` が消える**。

## ★★★★これで何が出るか

`sum_mu_d2xtermE`（本ファイル）は `hu` を取らない。
★`p ∣ l`（`1 − ζ^i` が単元でない）でも成り立つ。
☆これが `Lemma 3.5` から `d + 1 < l` を外す道の**要**である。
-/

namespace ABC3.Found.GaloisRep

open Finset ABC3.Meta

section Field

variable {F : Type} [Field F] [CharZero F]

/-- ☆`ζ^i ≠ 1` なら `∑_{k<l} (ζ^i)^k = 0`。 -/
theorem sum_pow_eq_zero_of_ne_one {l : ℕ} {η : F} (hpow : η ^ l = 1) (hne : η ≠ 1) :
    ∑ k ∈ range l, η ^ k = 0 := by
  have h := geom_sum_mul η l
  rw [hpow, sub_self] at h
  rcases mul_eq_zero.1 h with h1 | h2
  · exact h1
  · exact absurd (sub_eq_zero.1 h2) hne

/-- ★★★★★★★★★★**体版**: `120·∑ D²X^E = l⁴(l⁴ − 1)`。 -/
theorem sum_mu_d2xtermE_field {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    120 * (∑ i ∈ (range l).erase 0, tateD2XtermE l (ζ ^ i))
      = (l : F) ^ 4 * ((l : F) ^ 4 - 1) := by
  have hlne : ((l : F)) ≠ 0 := by
    simpa using (Nat.cast_ne_zero (R := F)).2 hl.ne_zero
  have hlu : IsUnit ((l : F)) := hlne.isUnit
  have hbridge : ∀ i ∈ (range l).erase 0,
      tateD2XtermE l (ζ ^ i) = (l : F) ^ 4 * tateD2Xterm (ζ ^ i) := by
    intro i hi
    have hne : (1 : F) - ζ ^ i ≠ 0 := one_sub_zeta_ne_zero hζ hi
    have hu : IsUnit ((1 : F) - ζ ^ i) := hne.isUnit
    have hpow : (ζ ^ i) ^ l = 1 := by
      rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
    have hzne : ζ ^ i ≠ 1 := by
      intro hc
      exact hne (by rw [hc, sub_self])
    exact (natCast_pow_mul_tateD2Xterm hlu hu hpow
      (sum_pow_eq_zero_of_ne_one hpow hzne)).symm
  rw [Finset.sum_congr rfl hbridge, ← Finset.mul_sum]
  rw [show (120 : F) * ((l : F) ^ 4 * ∑ i ∈ (range l).erase 0, tateD2Xterm (ζ ^ i))
      = (l : F) ^ 4 * (120 * ∑ i ∈ (range l).erase 0, tateD2Xterm (ζ ^ i)) from by ring,
    sum_mu_d2xterm_field hl hζ]

/-- ★★★★★★★★**体版**: `∑ DX^E = 0`。 -/
theorem sum_mu_dxtermE_field {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    ∑ i ∈ (range l).erase 0, tateDXtermE l (ζ ^ i) = 0 := by
  have hlne : ((l : F)) ≠ 0 := by
    simpa using (Nat.cast_ne_zero (R := F)).2 hl.ne_zero
  have hlu : IsUnit ((l : F)) := hlne.isUnit
  have hbridge : ∀ i ∈ (range l).erase 0,
      tateDXtermE l (ζ ^ i) = (l : F) ^ 3 * tateDXterm (ζ ^ i) := by
    intro i hi
    have hne : (1 : F) - ζ ^ i ≠ 0 := one_sub_zeta_ne_zero hζ hi
    have hu : IsUnit ((1 : F) - ζ ^ i) := hne.isUnit
    have hpow : (ζ ^ i) ^ l = 1 := by
      rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
    have hzne : ζ ^ i ≠ 1 := fun hc => hne (by rw [hc, sub_self])
    exact (natCast_pow_mul_tateDXterm hlu hu hpow
      (sum_pow_eq_zero_of_ne_one hpow hzne)).symm
  rw [Finset.sum_congr rfl hbridge, ← Finset.mul_sum, sum_mu_dxterm_field hl hζ, mul_zero]

end Field

section Domain

variable {A : Type} [CommRing A] [IsDomain A] [CharZero A]

/-- ★★★★★★★★★★★★★★★★★★★★
**整域版（`hu` を取らない）**: `120·∑ D²X^E = l⁴(l⁴ − 1)`（第 1112）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1112）**——`E` 版が環準同型と無条件に可換（第 1111）なので、
在庫の転送（分数体 → 体版 → 単射性）から **`hu` が落ちる**。
☆`p ∣ l` でも成り立つ。 -/
theorem sum_mu_d2xtermE {l : ℕ} (hl : l.Prime) {ζ : A} (hζ : IsPrimitiveRoot ζ l) :
    120 * (∑ i ∈ (range l).erase 0, tateD2XtermE l (ζ ^ i))
      = (l : A) ^ 4 * ((l : A) ^ 4 - 1) := by
  have hinj : Function.Injective (algebraMap A (FractionRing A)) :=
    IsFractionRing.injective A (FractionRing A)
  haveI : CharZero (FractionRing A) := charZero_of_injective_algebraMap hinj
  refine hinj ?_
  have hmap : (algebraMap A (FractionRing A))
      (120 * ∑ i ∈ (range l).erase 0, tateD2XtermE l (ζ ^ i))
      = 120 * ∑ i ∈ (range l).erase 0,
          tateD2XtermE l ((algebraMap A (FractionRing A)) ζ ^ i) := by
    rw [map_mul, map_sum]
    congr 1
    · exact map_ofNat _ 120
    · refine Finset.sum_congr rfl (fun i _ => ?_)
      rw [map_tateD2XtermE, map_pow]
  rw [hmap, map_mul, map_sub, map_pow, map_natCast, map_one]
  exact sum_mu_d2xtermE_field hl (hζ.map_of_injective hinj)

/-- ★★★★★★★★★★★★★★★★**整域版（`hu` を取らない）**: `∑ DX^E = 0`（第 1113）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`sum_mu_dxpair_zero` の頭項側である。★`p ∣ l` でも成り立つ。 -/
theorem sum_mu_dxtermE {l : ℕ} (hl : l.Prime) {ζ : A} (hζ : IsPrimitiveRoot ζ l) :
    ∑ i ∈ (range l).erase 0, tateDXtermE l (ζ ^ i) = 0 := by
  have hinj : Function.Injective (algebraMap A (FractionRing A)) :=
    IsFractionRing.injective A (FractionRing A)
  haveI : CharZero (FractionRing A) := charZero_of_injective_algebraMap hinj
  refine hinj ?_
  have hmap : (algebraMap A (FractionRing A))
      (∑ i ∈ (range l).erase 0, tateDXtermE l (ζ ^ i))
      = ∑ i ∈ (range l).erase 0,
          tateDXtermE l ((algebraMap A (FractionRing A)) ζ ^ i) := by
    rw [map_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [map_tateDXtermE, map_pow]
  rw [hmap, map_zero]
  exact sum_mu_dxtermE_field hl (hζ.map_of_injective hinj)

/-! ## ★★★★★★★★★★第 1115 —— 対の `E` 版

☆第 1114 で測ったとおり、`sum_mu_dxpair_zero` の証明で `hu` を使うのは
頭項の和（`hS`）だけである。★そこを第 1113 の `sum_mu_dxtermE` に差し替えれば
`hu` が消える。☆項ごとの書き換え（`hterm`）で
`tateDXpair (ζ^i) (q·(ζ^i)^{l-1}) = tateDXterm(ζ^i) + tateDXtail(ζ^i) − tateDXtail(ζ^{l-i})`
となるので、`l³` を掛けた形を直接述べる。 -/

variable {I : Ideal A} [IsAdicComplete I A]

/-- ★★★★★★★★★★★★★★★★★★★★
**対の `E` 版（`hu` を取らない）**（第 1115）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`sum_mu_dxpair_zero` の分母なし版である。☆`p ∣ l` でも成り立つ。 -/
theorem sum_mu_dxpairE_zero {l : ℕ} (hl : l.Prime) {ζ : A} (hζ : IsPrimitiveRoot ζ l)
    (q : A) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0,
        (tateDXtermE l (ζ ^ i)
          + (l : A) ^ 3 * (tateDXtail (ζ ^ i) q hq - tateDXtail (ζ ^ (l - i)) q hq)) = 0 := by
  rw [Finset.sum_add_distrib, sum_mu_dxtermE hl hζ, zero_add, ← Finset.mul_sum,
    Finset.sum_sub_distrib, sum_erase_reflect (fun j => tateDXtail (ζ ^ j) q hq),
    sub_self, mul_zero]

def sum_mu_dxpairE_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(å¯¾ã® E çââhu ãåããªããp â£ l ã§ãæãç«ã¤)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_dxpairE_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_mu_dxtermE(é ­é ã®åãç¬¬ 1113ãè¨¼ææ¸ã¿)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_dxtermE") 1,
    .citation "[ABC3]" "sum_erase_reflect(i â¦ lâi ã®å¥ãæ¿ããå¨åº«)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_erase_reflect") 1 ]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`D²X` の対の `E` 版（`hu` を取らない）**（第 1116）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`sum_mu_d2xpair` の分母なし版である。☆`c4_velu_tate` の本体。
★頭項は第 1112（`hu` なし）、尾は在庫の `sum_mu_d2xtail_sigma`（これも `hu` なし）。 -/
theorem sum_mu_d2xpairE {l : ℕ} (hl : l.Prime) {ζ : A} (hζ : IsPrimitiveRoot ζ l)
    (q : A) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    120 * ∑ i ∈ (range l).erase 0,
        (tateD2XtermE l (ζ ^ i)
          + (l : A) ^ 4 * (tateD2Xtail (ζ ^ i) q hq + tateD2Xtail (ζ ^ (l - i)) q hq))
      = (l : A) ^ 4 * (((l : A) ^ 4 - 1)
          + 240 * ((l : A) ^ 4 * evalAdic (sigmaSeries 3) (q ^ l) hql
              - evalAdic (sigmaSeries 3) q hq)) := by
  have hswap : ∑ i ∈ (range l).erase 0, tateD2Xtail (ζ ^ (l - i)) q hq
      = ∑ i ∈ (range l).erase 0, tateD2Xtail (ζ ^ i) q hq :=
    sum_erase_reflect (fun j => tateD2Xtail (ζ ^ j) q hq)
  have htail := sum_mu_d2xtail_sigma hl hζ q hq hql
  have hhead := sum_mu_d2xtermE hl hζ
  rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_add_distrib, hswap, htail]
  linear_combination hhead

def sum_mu_d2xpairE.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(DÂ²X ã®å¯¾ã® E çââhu ãåããªããc4_velu_tate ã®æ¬ä½)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_d2xpairE.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_mu_d2xtermE(é ­é ã®åãç¬¬ 1112ãè¨¼ææ¸ã¿)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_d2xtermE") 1,
    .citation "[ABC3]" "sum_mu_d2xtail_sigma(å°¾ã®åãå¨åº«ãhu ä¸è¦)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_d2xtail_sigma") 1 ]

end Domain

/-! ## ★出典の紐付け(`.src`) -/

def sum_mu_d2xtermE.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(D²X の頭項の和を hu なしで——p ∣ l でも成り立つ)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_d2xtermE.needs : List ProofObligation :=
  [ .citation "[ABC3]" "map_tateD2XtermE(E 版は環準同型と無条件に可換、第 1111、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.map_tateD2XtermE") 1,
    .citation "[ABC3]" "sum_mu_d2xterm_field(体版、在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_d2xterm_field") 1,
    .implicitStep
      ("★★**2026-09-01（第 1112）**——これが `Lemma 3.5` から `d + 1 < l` を" ++
       "外す道の要である。☆同じ形で `DX`・`D³X` の版も出る。" ++
       "★そのうえで `sum_mu_d2xpair` を `E` 版で述べ直せば `hu` が消える。") 8 ]

end ABC3.Found.GaloisRep
