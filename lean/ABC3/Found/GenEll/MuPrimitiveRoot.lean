/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.MuRational

/-!
# 第 947 ブロック —— **★★★★★★★★★★★★★★★★★★★★有理な `l`-捉れ点から
原始 `l` 乗根を取る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——**946 と 927 を繋ぐ**

第 946 は `∃ ζ : Kˣ, ζ^l = 1 ∧ P = tatePhi([ζ])` を与える。
だが `tateParam_quot_velu_dvr`（第 927）が欲しいのは

    `ζ : R`、`IsPrimitiveRoot ζ l`、`uζ : Kˣ`、`algebraMap R K ζ = uζ`、
    `uζ^l = 1`、`∀ n, 0 < n → n < l → uζ^n ≠ 1`

の 6 つである。★本ブロックは**その 6 つを一気に作る**。

☆道は 3 段:

1. `P ≠ 0` なら `ζ ≠ 1`（`tatePhi(1) = 0` だから）
2. `l` は素だから `orderOf ζ = l`、すなわち原始 `l` 乗根
3. 1 の冪根の付値は `0`（第 896）なので `ζ ∈ R`（`TateSetup.hmem0`）

| 定理 | 内容 |
|---|---|
| `isPrimitiveRoot_of_pow_eq_one_of_ne_one` | ★`ζ^l = 1`・`ζ ≠ 1`・`l` 素 ⇒ 原始 |
| `exists_primitiveRoot_of_torsion_point` | ★★★★★★★★★★★★★★★★★★★★**927 が欲しい 6 つを一気に** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GaloisRep

/-! ## ★素数位数の判定 -/

/-- ★**`ζ^l = 1`、`ζ ≠ 1`、`l` が素なら `ζ` は原始 `l` 乗根**。

☆`orderOf ζ ∣ l` で `l` は素だから `orderOf ζ` は `1` か `l`。
`1` なら `ζ = 1` なので除外される。 -/
theorem isPrimitiveRoot_of_pow_eq_one_of_ne_one {G : Type} [CommMonoid G] {l : ℕ} (hl : l.Prime)
    (ζ : G) (hζ : ζ ^ l = 1) (hne : ζ ≠ 1) : IsPrimitiveRoot ζ l := by
  have hord : orderOf ζ = l := by
    rcases hl.eq_one_or_self_of_dvd _ (orderOf_dvd_of_pow_eq_one hζ) with h | h
    · exact absurd (orderOf_eq_one_iff.1 h) hne
    · exact h
  refine ⟨hζ, fun m hm => ?_⟩
  rw [← hord]
  exact orderOf_dvd_of_pow_eq_one hm

/-- ★**原始 `l` 乗根は `l` 未満の正の冪で `1` にならない**。 -/
theorem pow_ne_one_of_isPrimitiveRoot {G : Type} [CommMonoid G] {l : ℕ}
    {ζ : G} (hζ : IsPrimitiveRoot ζ l) (n : ℕ) (hn : 0 < n) (hnl : n < l) : ζ ^ n ≠ 1 := by
  intro h
  exact absurd (Nat.le_of_dvd hn (hζ.dvd_of_pow_eq_one n h)) (not_le.2 hnl)

/-! ## ★★★★★★★★★★★★★★★★★★★★927 が欲しい 6 つを一気に -/

section Dvr

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★★★★★★★★★★★★★★★★★**有理な `l`-捉れ点から
`tateParam_quot_velu_dvr` が欲しい 6 つを一気に作る**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 947）**——これで 946（`ζ` の存在）と
927（`q_{E′} = q_E^l`）が**直接繋がる**。
☆引数は `q`・`Δ ≠ 0`・`l` が素・`l ∤ v(q)`・`l • P = 0`・`P ≠ 0` だけである。 -/
theorem exists_primitiveRoot_of_torsion_point
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := K) q hq hq0).v
      (mkTateSetup (K := K) q hq hq0).Q))
    (P : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0) :
    ∃ (ζ : R) (uζ : Kˣ), IsPrimitiveRoot ζ l
      ∧ algebraMap R K ζ = (uζ : K)
      ∧ uζ ^ l = 1
      ∧ (∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
      ∧ P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ) := by
  -- ★段 1: 第 946
  obtain ⟨uζ, hζl, hPz⟩ := exists_mu_point_dvr q hq hq0 hΔ hl hcop P hP
  -- ★段 2: `P ≠ 0` なら `uζ ≠ 1`
  have hne1 : uζ ≠ 1 := by
    intro h
    apply hP0
    have hone : (QuotientGroup.mk (1 : Kˣ) :
        Kˣ ⧸ Subgroup.zpowers (mkTateSetup (K := K) q hq hq0).Q) = 1 := rfl
    rw [hPz, h, hone, tatePhi_one]
  have hprimU : IsPrimitiveRoot uζ l :=
    isPrimitiveRoot_of_pow_eq_one_of_ne_one hl uζ hζl hne1
  -- ★段 3: 1 の冪根の付値は `0` なので `ζ ∈ R`
  have hv0 : vAdd (mkTateSetup (K := K) q hq hq0).v uζ = 0 :=
    vAdd_eq_zero_of_pow_eq_one _ hl.pos uζ hζl
  obtain ⟨ζR, hζu⟩ := (mkTateSetup (K := K) q hq hq0).hmem0 uζ (le_of_eq hv0.symm)
  -- ★段 4: `R` の側でも原始
  have hcoe : ((uζ : K)) ^ l = 1 := by
    rw [← Units.val_pow_eq_pow_val, hζl, Units.val_one]
  have hζRpow : ζR ^ l = 1 := by
    refine IsFractionRing.injective R K ?_
    rw [map_pow, hζu, hcoe, map_one]
  have hζR : IsPrimitiveRoot ζR l := by
    refine ⟨hζRpow, fun m hm => ?_⟩
    have hmK : ((uζ : K)) ^ m = 1 := by
      rw [← hζu, ← map_pow, hm, map_one]
    have hmU : uζ ^ m = 1 :=
      Units.ext (by rw [Units.val_pow_eq_pow_val, hmK, Units.val_one])
    exact hprimU.dvd_of_pow_eq_one m hmU
  exact ⟨ζR, uζ, hζR, hζu, hζl,
    fun n hn hnl => pow_ne_one_of_isPrimitiveRoot hprimU n hn hnl, hPz⟩

end Dvr

/-! ## ★出典の紐付け(`.src`) -/

def isPrimitiveRoot_of_pow_eq_one_of_ne_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ζ^l = 1 で ζ ≠ 1 なら原始 l 乗根)",
    sectionId := "genell-lemma-3-5" }

def pow_ne_one_of_isPrimitiveRoot.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(原始 l 乗根は l 未満の冪で 1 にならない)",
    sectionId := "genell-lemma-3-5" }

def exists_primitiveRoot_of_torsion_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有理な l-捉れ点から原始 l 乗根を取る)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
