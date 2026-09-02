/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.MuPrimitiveRoot
import ABC3.Found.GenEll.MuOrDeep
import ABC3.Meta.Claim

/-!
# 第 1414 ブロック —— **★★★★★★★★★★★★★★★★原始根か深い代表か**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——第 947 から `hcop` を外す

第 947 の `exists_primitiveRoot_of_torsion_point` は
`l ∤ v(q)`（＝`hcop`）を仮定して「捻れ点から原始 `l` 乗根が取れる」を出していた。
★`hcop` を使うのは**最初の 1 段だけ**（第 946 の `exists_mu_point_dvr`）で、
残りの 4 段（`ζ ≠ 1`・原始性・`ζ ∈ R`・`R` の側の原始性）は `hcop` と無関係である。

☆本ブロックはその 4 段を **`primitiveRoot_of_mu_point`** として切り出し、
第 1410-1413 の二者択一と繋いで

| 場合 | 結論 |
|---|---|
| `μ_l` 型 | ★第 947 と同じ 6 つ（`ζ : R`・`IsPrimitiveRoot ζ l`・…） |
| 深い代表 | ★★`∀ 0 < i < l, v(Q) ∤ i·v(y)`（＝第 1412 の `hdeep`） |

を出す。★★★これが `semistableAt_veluQuotient_bad_ram`（第 1404）の
入口を **`hcop` なし**にするための一本である。

| 定理 | 内容 |
|---|---|
| `primitiveRoot_of_mu_point` | ★★★★`μ_l` 型の点から原始根の 6 つ（`hcop` 不要） |
| `primitiveRoot_or_deep_of_torsion_point` | ★★★★★★★★★★★★★★★★**二者択一** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GaloisRep

section Dvr

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★**`μ_l` 型の点から `tateParam_quot_velu_dvr` が欲しい 6 つを作る**（第 1414）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 947 の段 2-5 をそのまま切り出したものである
——★**`l ∤ v(q)` を使わない**のが要点。

1. `P ≠ 0` なら `ζ ≠ 1`（`tatePhi(1) = 0` だから）
2. `l` は素だから `orderOf ζ = l`、すなわち原始 `l` 乗根
3. 1 の冪根の付値は `0`（第 896）なので `ζ ∈ R`
4. `R → K` は単射なので `R` の側でも原始 -/
theorem primitiveRoot_of_mu_point
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime)
    (P : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Point)
    (hP0 : P ≠ 0) (uζ : Kˣ) (hζl : uζ ^ l = 1)
    (hPz : P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ)) :
    ∃ (ζ : R) (uζ' : Kˣ), IsPrimitiveRoot ζ l
      ∧ algebraMap R K ζ = (uζ' : K)
      ∧ uζ' ^ l = 1
      ∧ (∀ n : ℕ, 0 < n → n < l → uζ' ^ n ≠ 1)
      ∧ P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ') := by
  -- ★段 1: `P ≠ 0` なら `ζ ≠ 1`
  have hne1 : uζ ≠ 1 := by
    intro h
    apply hP0
    have hone : (QuotientGroup.mk (1 : Kˣ) :
        Kˣ ⧸ Subgroup.zpowers (mkTateSetup (K := K) q hq hq0).Q) = 1 := rfl
    rw [hPz, h, hone, tatePhi_one]
  have hprimU : IsPrimitiveRoot uζ l :=
    isPrimitiveRoot_of_pow_eq_one_of_ne_one hl uζ hζl hne1
  -- ★段 2: 1 の冪根の付値は `0` なので `ζ ∈ R`
  have hv0 : vAdd (mkTateSetup (K := K) q hq hq0).v uζ = 0 :=
    vAdd_eq_zero_of_pow_eq_one _ hl.pos uζ hζl
  obtain ⟨ζR, hζu⟩ := (mkTateSetup (K := K) q hq hq0).hmem0 uζ (le_of_eq hv0.symm)
  -- ★段 3: `R` の側でも原始
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

/-- ★★★★★★★★★★★★★★★★**捻れ点は原始根を与えるか、深い代表を持つ**（第 1414）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-02（第 1414）**——第 947 の `exists_primitiveRoot_of_torsion_point`
から **`hcop`（`l ∤ v(q)`）が落ちた**（その代わり結論が二者択一になった）。
☆引数は `q`・`Δ ≠ 0`・`l` が素・`l • P = 0`・`P ≠ 0` の 5 つだけである。

★★★深い側の結論は第 1412 の `veluQuotientFull_tate_deep` の入力 `hdeep` そのもので、
そちらでは `c₄(veluCurve) = c₄(E_q) + 240v` が `v ∈ 𝔪` から単元になる。 -/
theorem primitiveRoot_or_deep_of_torsion_point
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime)
    (P : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0) :
    (∃ (ζ : R) (uζ : Kˣ), IsPrimitiveRoot ζ l
        ∧ algebraMap R K ζ = (uζ : K)
        ∧ uζ ^ l = 1
        ∧ (∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
        ∧ P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))
      ∨ (∃ y : Kˣ, P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk y)
          ∧ ∀ i : ℕ, 0 < i → i < l →
              ¬ (vAdd (mkTateSetup (K := K) q hq hq0).v
                  (mkTateSetup (K := K) q hq hq0).Q
                ∣ (i : ℤ) * vAdd (mkTateSetup (K := K) q hq hq0).v y)) := by
  rcases mu_or_deep_point_dvr q hq hq0 hΔ hl P hP with
    ⟨uζ, hζl, hPz⟩ | ⟨y, hPy, _, _, hdeep⟩
  · exact Or.inl (primitiveRoot_of_mu_point q hq hq0 hΔ hl P hP0 uζ hζl hPz)
  · exact Or.inr ⟨y, hPy, hdeep⟩

end Dvr

/-! ## ★出典の紐付け(`.src`) -/

def primitiveRoot_of_mu_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(μ_l 型の点から原始 l 乗根の 6 つ。★l ∤ v(q) 不要)",
    sectionId := "genell-lemma-3-5" }

def primitiveRoot_or_deep_of_torsion_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(捻れ点は原始根を与えるか深い代表を持つ。★l ∤ v(q) 不要)",
    sectionId := "genell-lemma-3-5" }

def primitiveRoot_or_deep_of_torsion_point.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "mu_or_deep_point_dvr(第 1410・1413、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mu_or_deep_point_dvr") 1,
    .citation "[ABC3]" "exists_primitiveRoot_of_torsion_point(第 947、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1414）**——第 947 が `hcop`（`l ∤ v(q)`）を使うのは" ++
       "**最初の 1 段だけ**（第 946 の `exists_mu_point_dvr`）だと測った。" ++
       "☆残りの 4 段を `primitiveRoot_of_mu_point` として切り出し、" ++
       "第 1410-1413 の二者択一と繋いだ。" ++
       "★★★これで `semistableAt_veluQuotient_bad_ram`（第 1404）の入口が " ++
       "`hcop` なしになる。") 17 ]

end ABC3.Found.GenEll
