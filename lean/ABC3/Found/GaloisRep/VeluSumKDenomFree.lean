/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateKDenomFree
import ABC3.Found.GaloisRep.VeluTateMuField
import ABC3.Meta.Claim

/-!
# 第 1134 ブロック —— **`K` 水準の Vélu の和は分母なし版の像**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か（節点 5 の消費者 C）

第 1131 の `exists_vw_tate_mu_field` は `v ≔ φ(SV)/l⁶`・`w ≔ φ(SW)/(2l⁸)` を
**天下り**に取っていた。★本ブロックはそれが**本当に `K` 水準の Vélu の和である**ことを示す:

    `l⁶ · ∑_i veluV2 (E_q ⊗ K) (tateXK ζ^i …) (tateYK ζ^i …) = φ(SV)`
    `l⁸ · ∑_i (veluU … + 2·veluV2 … · tateXK …) = φ(SW)`

☆機構は 3 つの無条件恒等式の合成だけである:

| 段 | 使うもの | 第 |
|---|---|---|
| `l⁶·veluV2 W x y = veluV2DF l W (l²x) (l³y)` | `natCast_pow_mul_veluV2` | 1119 |
| `l²·tateXK = φ(tateXpairDF)`・`l³·tateYK = φ(tateYpairDF)` | 第 1133 | 1133 |
| `φ(veluV2DF l W X Y) = veluV2DF l (W.map φ) (φX) (φY)` | `map_veluV2DF` | 1125 |

★どれも `IsUnit (1 − ζ^i)` も `IsUnit (l)` も要求しない。
☆要るのは `φ(1 − ζ^i) ≠ 0` だけで、`R → K` の単射性から出る。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll Finset

section Pointwise

variable {R K : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  [Field K] [Algebra R K]

/-- ★★★★★★★★**`l⁶·veluV2`（`K` の座標）は `veluV2DF` の像**——★`IsUnit` 不要。 -/
theorem natCast_pow_mul_veluV2_K {l : ℕ} {a : R} (hpow : a ^ l = 1)
    (hsum : ∑ k ∈ range l, a ^ k = 0)
    (hne : algebraMap R K (1 - a) ≠ 0) (W : WeierstrassCurve R) (w q : R) (hq : q ∈ I) :
    (l : K) ^ 6 * veluV2 (W.map (algebraMap R K))
        (tateXK (I := I) a w q hq) (tateYK (I := I) a w q hq)
      = algebraMap R K (veluV2DF l W (tateXpairDF l a w q hq) (tateYpairDF l a w q hq)) := by
  rw [natCast_pow_mul_veluV2, natCast_pow_mul_tateXK hpow hsum hne w q hq,
    natCast_pow_mul_tateYK hpow hsum hne w q hq, map_veluV2DF]

/-- ★★★★★★★★**`l⁶·veluU`（`K` の座標）は `veluUDF` の像**——★`IsUnit` 不要。 -/
theorem natCast_pow_mul_veluU_K {l : ℕ} {a : R} (hpow : a ^ l = 1)
    (hsum : ∑ k ∈ range l, a ^ k = 0)
    (hne : algebraMap R K (1 - a) ≠ 0) (W : WeierstrassCurve R) (w q : R) (hq : q ∈ I) :
    (l : K) ^ 6 * veluU (W.map (algebraMap R K))
        (tateXK (I := I) a w q hq) (tateYK (I := I) a w q hq)
      = algebraMap R K (veluUDF l W (tateXpairDF l a w q hq) (tateYpairDF l a w q hq)) := by
  rw [natCast_pow_mul_veluU, natCast_pow_mul_tateXK hpow hsum hne w q hq,
    natCast_pow_mul_tateYK hpow hsum hne w q hq, map_veluUDF]

/-- ★★★★★★★★★★**`l⁸·(veluU + 2·veluV2·X)` は `SW` の被加数の像**——★`IsUnit` 不要。 -/
theorem natCast_pow_mul_veluW_K {l : ℕ} {a : R} (hpow : a ^ l = 1)
    (hsum : ∑ k ∈ range l, a ^ k = 0)
    (hne : algebraMap R K (1 - a) ≠ 0) (W : WeierstrassCurve R) (w q : R) (hq : q ∈ I) :
    (l : K) ^ 8 * (veluU (W.map (algebraMap R K))
          (tateXK (I := I) a w q hq) (tateYK (I := I) a w q hq)
        + 2 * (veluV2 (W.map (algebraMap R K))
                (tateXK (I := I) a w q hq) (tateYK (I := I) a w q hq)
              * tateXK (I := I) a w q hq))
      = algebraMap R K ((l : R) ^ 2
            * veluUDF l W (tateXpairDF l a w q hq) (tateYpairDF l a w q hq)
          + 2 * (veluV2DF l W (tateXpairDF l a w q hq) (tateYpairDF l a w q hq)
                * tateXpairDF l a w q hq)) := by
  have hU := natCast_pow_mul_veluU_K (K := K) hpow hsum hne W w q hq
  have hV := natCast_pow_mul_veluV2_K (K := K) hpow hsum hne W w q hq
  have hX := natCast_pow_mul_tateXK (K := K) hpow hsum hne w q hq
  rw [map_add, map_mul, map_mul, map_mul, map_pow, map_natCast, map_ofNat,
    ← hU, ← hV, ← hX]
  ring

end Pointwise

/-! ## ★★★★★★★★★★和の水準 -/

section Sum

variable {R K : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  [Field K] [Algebra R K]

/-- ☆整域では `η^l = 1`・`η ≠ 1` から `∑_{k<l} η^k = 0`。 -/
theorem sum_pow_eq_zero_of_ne_one_domain {l : ℕ} {η : R} (hpow : η ^ l = 1) (hne : η ≠ 1) :
    ∑ k ∈ range l, η ^ k = 0 := by
  have h := geom_sum_mul η l
  rw [hpow, sub_self] at h
  rcases mul_eq_zero.1 h with h1 | h2
  · exact h1
  · exact absurd (by linear_combination h2) hne

omit [IsDomain R] in
/-- ☆`0 < i < l` なら `ζ^i ≠ 1`。 -/
theorem zeta_pow_ne_one {l : ℕ} {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    {i : ℕ} (hi : i ∈ (range l).erase 0) : ζ ^ i ≠ 1 := fun h =>
  (Finset.mem_erase.1 hi).1
    (Nat.eq_zero_of_dvd_of_lt (hζ.dvd_of_pow_eq_one i h)
      (Finset.mem_range.1 (Finset.mem_erase.1 hi).2))

omit [IsDomain R] in
/-- ☆`(ζ^i)^l = 1`。 -/
theorem zeta_pow_pow {l : ℕ} {ζ : R} (hζ : IsPrimitiveRoot ζ l) (i : ℕ) :
    (ζ ^ i) ^ l = 1 := by
  rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]

/-- ☆`R → K` が単射なら `φ(1 − ζ^i) ≠ 0`。 -/
theorem zeta_pow_sub_ne_zero_K (hinj : Function.Injective (algebraMap R K))
    {l : ℕ} {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    {i : ℕ} (hi : i ∈ (range l).erase 0) : algebraMap R K (1 - ζ ^ i) ≠ 0 := by
  intro h
  have h0 : (1 : R) - ζ ^ i = 0 := hinj (by rw [h, map_zero])
  exact zeta_pow_ne_one hζ hi (by linear_combination -h0)

/-- ★★★★★★★★★★★★★★★★
**`K` 水準の `∑ veluV2` は `SV` の像**（第 1134）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★**`IsUnit (1 − ζ^i)` も `IsUnit (l)` も仮説に置いていない**。
☆これで第 1131 の `v = φ(SV)/l⁶` が**本当に Vélu の和である**ことが確かめられる。 -/
theorem sum_natCast_pow_mul_veluV2_K (hinj : Function.Injective (algebraMap R K))
    {l : ℕ} {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (W : WeierstrassCurve R) (q : R) (hq : q ∈ I) :
    (l : K) ^ 6 * ∑ i ∈ (range l).erase 0,
        veluV2 (W.map (algebraMap R K))
          (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      = algebraMap R K (∑ i ∈ (range l).erase 0,
          veluV2DF l W (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)) := by
  rw [Finset.mul_sum, map_sum]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  exact natCast_pow_mul_veluV2_K (zeta_pow_pow hζ i)
    (sum_pow_eq_zero_of_ne_one_domain (zeta_pow_pow hζ i) (zeta_pow_ne_one hζ hi))
    (zeta_pow_sub_ne_zero_K hinj hζ hi) W _ q hq

/-- ★★★★★★★★★★★★★★★★
**`K` 水準の `∑ (veluU + 2·veluV2·X)` は `SW` の像**（第 1134）。

☆これで第 1131 の `w = φ(SW)/(2l⁸)` が `2w = ∑(…)` を満たすことが確かめられる。 -/
theorem sum_natCast_pow_mul_veluW_K (hinj : Function.Injective (algebraMap R K))
    {l : ℕ} {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (W : WeierstrassCurve R) (q : R) (hq : q ∈ I) :
    (l : K) ^ 8 * ∑ i ∈ (range l).erase 0,
        (veluU (W.map (algebraMap R K))
            (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          + 2 * (veluV2 (W.map (algebraMap R K))
                  (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                * tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = algebraMap R K (∑ i ∈ (range l).erase 0,
          ((l : R) ^ 2 * veluUDF l W
                (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2DF l W
                    (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))) := by
  rw [Finset.mul_sum, map_sum]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  exact natCast_pow_mul_veluW_K (zeta_pow_pow hζ i)
    (sum_pow_eq_zero_of_ne_one_domain (zeta_pow_pow hζ i) (zeta_pow_ne_one hζ hi))
    (zeta_pow_sub_ne_zero_K hinj hζ hi) W _ q hq

end Sum

/-! ## ★出典の紐付け(`.src`) -/

def sum_natCast_pow_mul_veluV2_K.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(K 水準の ∑ veluV2 は SV の像。★IsUnit を取らない)",
    sectionId := "genell-lemma-3-5" }

def sum_natCast_pow_mul_veluW_K.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(K 水準の ∑(veluU + 2·veluV2·X) は SW の像。★IsUnit を取らない)",
    sectionId := "genell-lemma-3-5" }

def sum_natCast_pow_mul_veluV2_K.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "natCast_pow_mul_veluV2(l⁶·veluV2 の分母なし版、第 1119、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_pow_mul_veluV2") 1,
    .citation "[ABC3]" "natCast_pow_mul_tateXK(K 水準の座標の橋、第 1133、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_pow_mul_tateXK") 1,
    .citation "[ABC3]" "map_veluV2DF(veluV2DF は環準同型で移る、第 1125、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.map_veluV2DF") 1,
    .implicitStep
      ("★★**2026-09-01（第 1134）**——第 1131 の `v = φ(SV)/l⁶` は天下りに取っていたが、" ++
       "本ブロックでそれが**本当に `K` 水準の Vélu の和である**ことが確かめられた。" ++
       "☆機構は 3 つの無条件恒等式（第 1119・1133・1125）の合成だけである。") 4 ]

end ABC3.Found.GaloisRep
