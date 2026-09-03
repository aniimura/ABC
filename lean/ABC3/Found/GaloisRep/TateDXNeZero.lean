/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluPoints
import ABC3.Found.GaloisRep.TateDSeries

/-!
# Galois (G6) 第 893 ブロック —— **★★★★★★★★`DX = 0` は `P = −P` と同じ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か

`c4_velu_tate`・`c6_velu_tate`（第 853・867）は侧条件

    `hDX : tateDXpair (ζⁱ) (q(ζⁱ)^{l-1}) q hq ≠ 0`

を受けている。これは ODE `DX·(DY − (3X² − Y + a₄)) = 0` から `DX` を
約すためのものである。

★しかし `tateDXpair = 2Y + X` であり、Tate 曲線は `a₁ = 1`・`a₃ = 0` なので
`2Y + X = 2y + a₁x + a₃` —— つまり**`DX = 0` は `P = −P`（2-捉れ）と同じ**である。

☆だから `μ_l`（`l` は奇素数）の点では `hDX` は**仮説ではなく定理**になる。
本ブロックはその前半（座標の側）である。

| 定理 | 内容 |
|---|---|
| `tateCurveAtMap_a₁` / `_a₃` | ★`a₁ = 1`・`a₃ = 0` |
| `eq_neg_of_tateDXpair_eq_zero` | ★★★★★★★★**`DX = 0` なら `P = −P`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★Tate 曲線の `a₁` は `1`。 -/
@[simp] theorem tateCurveAtMap_a₁ (q : R) (hq : q ∈ I) :
    ((tateCurveAt q hq).map (algebraMap R K)).a₁ = 1 := by
  rw [tateCurveAt, WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₁]
  show algebraMap R K (evalAdicHom q hq (1 : PowerSeries ℤ)) = 1
  rw [map_one, map_one]

/-- ★Tate 曲線の `a₃` は `0`。 -/
@[simp] theorem tateCurveAtMap_a₃ (q : R) (hq : q ∈ I) :
    ((tateCurveAt q hq).map (algebraMap R K)).a₃ = 0 := by
  rw [tateCurveAt, WeierstrassCurve.map_a₃, WeierstrassCurve.map_a₃]
  show algebraMap R K (evalAdicHom q hq (0 : PowerSeries ℤ)) = 0
  rw [map_zero, map_zero]

/-- ★★★★★★★★**`DX = 0` なら `P = −P`**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

☆`tateDXpair = 2Y + X` で、Tate 曲線では `negY x y = −y − x` だから。 -/
theorem eq_neg_of_tateDXpair_eq_zero (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hw : IsUnit (1 - w)) (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (ha : IsUnit (1 - a)) (h0 : tateDXpair a w q hq = 0) :
    tatePtPair a w q hq haw hw hne hΔ = -(tatePtPair a w q hq haw hw hne hΔ) := by
  have hRel : 2 * tateYpair a w q hq + tateXpair a w q hq = 0 := by
    rw [← tateDXpair_eq a w q hq ha hw]
    exact h0
  have hK : (2 : K) * (tateYK a w q hq : K) + (tateXK a w q hq : K) = 0 := by
    rw [tateXK_eq a w q hq ha, tateYK_eq a w q hq ha, ← map_ofNat (algebraMap R K) 2,
      ← map_mul, ← map_add, hRel, map_zero]
  rw [tatePtPair, Point.neg_some, Point.some.injEq]
  refine ⟨rfl, ?_⟩
  rw [WeierstrassCurve.Affine.negY]
  show (tateYK a w q hq : K)
    = -(tateYK a w q hq : K)
      - ((tateCurveAt q hq).map (algebraMap R K)).a₁ * (tateXK a w q hq : K)
      - ((tateCurveAt q hq).map (algebraMap R K)).a₃
  rw [tateCurveAtMap_a₁, tateCurveAtMap_a₃]
  linear_combination hK

/-! ## ★★★★★★★★★★`μ_l` の点では `hDX` は定理になる -/

section Mu

variable [IsDomain R]
open scoped Classical

/-- ★`ζ` の位数がちょうど `l` なら、`ζⁿ = 1` は `l ∣ n` を意味する。 -/
theorem dvd_of_pow_eq_one {l : ℕ} (hl : 0 < l) (uζ : Kˣ) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (n : ℕ) (hn : uζ ^ n = 1) : l ∣ n := by
  by_contra hnd
  refine hord (n % l) (Nat.pos_of_ne_zero (fun h => hnd (Nat.dvd_of_mod_eq_zero h)))
    (Nat.mod_lt _ hl) ?_
  have hsplit : uζ ^ n = (uζ ^ l) ^ (n / l) * uζ ^ (n % l) := by
    rw [← pow_mul, ← pow_add]
    congr 1
    exact (Nat.div_add_mod n l).symm
  rw [hsplit, hζl, one_pow, one_mul] at hn
  exact hn

/-- ★★★★★★★★★★**`μ_l` の自明でない点では `DX ≠ 0`**（`l` は奇素数）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`c4_velu_tate`・`c6_velu_tate` が受けていた侧条件 `hDX` は、
`μ_l` の設定では**仮説ではなく定理**である。

☆`DX = 0` は `P = −P`（第 893）だから `2P = 0`、
すなわち `ζ^{2i} = 1`、つまり `l ∣ 2i`。
`l` は奇素数で `0 < i < l` なので矛盾する。 -/
theorem tateDXpair_ne_zero_of_mu (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (ζ : R) (uζ : Kˣ)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (i : ℕ) (hi0 : i ≠ 0) (hil : i < l) (hu : IsUnit (1 - ζ ^ i)) :
    tateDXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq ≠ 0 := by
  intro h0
  -- ★準備
  have hζlR : ζ ^ l = 1 := pow_eq_one_of_map S ζ uζ hζu hζl
  have hmap : algebraMap R K (ζ ^ i) = ((uζ ^ i : Kˣ) : K) := by
    rw [map_pow, hζu, Units.val_pow_eq_pow_val]
  have hpowl : (uζ ^ i) ^ l = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have haw : (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1)) = S.q := by
    have h1 : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      have hl' : l - 1 + 1 = l := Nat.succ_pred_eq_of_pos hl.pos
      calc (ζ ^ i) * (ζ ^ i) ^ (l - 1) = (ζ ^ i) ^ (l - 1 + 1) := by ring
        _ = (ζ ^ i) ^ l := by rw [hl']
        _ = (ζ ^ l) ^ i := by rw [← pow_mul, mul_comm, pow_mul]
        _ = 1 := by rw [hζlR, one_pow]
    calc (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1))
        = S.q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = S.q := by rw [h1, mul_one]
  have hwu : IsUnit (1 - S.q * (ζ ^ i) ^ (l - 1)) :=
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ S.hq)
  have hne : algebraMap R K (1 - ζ ^ i) ≠ 0 := (hu.map (algebraMap R K)).ne_zero
  have hmk0 : (QuotientGroup.mk (uζ ^ 0) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
    rw [pow_zero]; rfl
  have hcne : (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 := by
    intro hc
    exact hi0 (mk_pow_injOn S.v S.Q S.hQ hl.pos uζ hζl hord hil hl.pos (by rw [hc, hmk0]))
  -- ★段 1: `P = −P`
  have hP := tatePhi_of_pow_eq_one S hΔ hl.pos (ζ ^ i) (uζ ^ i) hmap hpowl hcne haw hwu hne
  have hneg := eq_neg_of_tateDXpair_eq_zero (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq
    haw hwu hne hΔ hu h0
  rw [← hP] at hneg
  -- ★段 2: 類の側へ戻す
  have hinv : -(tatePhi S hΔ (QuotientGroup.mk (uζ ^ i)))
      = tatePhi S hΔ ((QuotientGroup.mk (uζ ^ i))⁻¹) := by
    rw [← hΦ, ← hΦ, ← map_neg]
    rfl
  rw [hinv, ← hΦ, ← hΦ] at hneg
  have hcl : (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q)
      = (QuotientGroup.mk (uζ ^ i))⁻¹ := Φ.injective hneg
  -- ★段 3: `ζ^{2i} = 1`
  have hsq : (QuotientGroup.mk (uζ ^ (2 * i)) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
    have h2 : (uζ ^ (2 * i)) = (uζ ^ i) * (uζ ^ i) := by
      rw [← pow_add]; congr 1; omega
    rw [h2, QuotientGroup.mk_mul]
    nth_rewrite 1 [hcl]
    rw [inv_mul_cancel]
  have hpow2 : (uζ ^ (2 * i)) ^ l = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have hone := eq_one_of_mk_eq_one_of_pow_eq_one S.v S.Q S.hQ hl.pos (uζ ^ (2 * i)) hpow2 hsq
  -- ★段 4: `l ∣ 2i` は矛盾
  have hdvd := dvd_of_pow_eq_one hl.pos uζ hζl hord (2 * i) hone
  rcases (Nat.Prime.dvd_mul hl).1 hdvd with h2 | hi
  · exact hodd (Nat.le_antisymm (Nat.le_of_dvd (by norm_num) h2) hl.two_le)
  · exact absurd (Nat.le_of_dvd (Nat.pos_of_ne_zero hi0) hi) (by omega)

end Mu

/-! ## ★出典の紐付け(`.src`) -/

def tateDXpair_ne_zero_of_mu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l の自明でない点では DX ≠ 0。★無条件)",
    sectionId := "genell-lemma-3-2" }

def eq_neg_of_tateDXpair_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(DX = 0 は P = −P と同じ。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
