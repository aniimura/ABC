/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.SemistableScaled
import ABC3.Meta.Claim

/-!
# 第 1430 ブロック —— **短 Weierstrass 形で `minDeltaExp = 0`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1425 の**良い素点版**

第 1425 は `v(c₄) = 4·v(n)`・`v(c₆) = 6·v(n)` から `v_p(c₄) = 0` の整モデルを作り、
**乗法還元の側**（`SemistableAt` の第 2 の選択肢）を出した。

★良い素点で欲しいのは第 1 の選択肢 **`minDeltaExp p W = 0`** であり、
そちらには**不等式で足りる**:

    `c₄/(48·n⁴) ∈ 𝒪_p`、`c₆/(864·n⁶) ∈ 𝒪_p`、`v_p(Δ) = 12·v_p(n)`
      ⟹ `minDeltaExp p W = 0`

☆機構は第 1425 と同じ——`toShortNF`（`u = 1`）で `a₄ = −c₄/48`・`a₆ = −c₆/864` にし、
`u = n` で割ると係数がちょうど整になり、`Δ` の付値は `12v(n)` 減って `0` になる。
★`c₄ = 0`・`c₆ = 0` の場合も込みで扱える（`0 ∈ 𝒪_p`）。

| 定理 | 内容 |
|---|---|
| `div_mem_primeSubring_of_valAdd_le` | ☆`v(c) ≥ k·v(n)` なら `c/(d·nᵏ) ∈ 𝒪_p` |
| `minDeltaExp_eq_zero_of_short_integral` | ★★★★★★★★**`minDeltaExp = 0`** |
| `semistableAt_of_short_integral` | ★★★★★★★★半安定 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ☆`v(c) ≥ k·v(n)` かつ `v(d) = 0` なら `c/(d·nᵏ)` は付値環の元（`c = 0` も込み）。 -/
theorem div_mem_primeSubring_of_valAdd_le (p : HeightOneSpectrum (𝓞 L))
    {c d : L} (hd : d ≠ 0) (hd0 : valAdd p (Units.mk0 d hd) = 0)
    {n : L} (hn : n ≠ 0) (k : ℕ)
    (hc : ∀ hc0 : c ≠ 0, (k : ℤ) * valAdd p (Units.mk0 n hn)
      ≤ valAdd p (Units.mk0 c hc0)) :
    c / (d * n ^ k) ∈ primeSubring p := by
  rcases eq_or_ne c 0 with rfl | hc0
  · simp
  · have hne : c / (d * n ^ k) ≠ 0 :=
      div_ne_zero hc0 (mul_ne_zero hd (pow_ne_zero k hn))
    refine (mem_primeSubring_iff p _).2 ((valAdd_nonneg_iff p (Units.mk0 _ hne)).1 ?_)
    have hEq : Units.mk0 (c / (d * n ^ k)) hne
        = Units.mk0 c hc0 * (Units.mk0 d hd)⁻¹ * ((Units.mk0 n hn)⁻¹) ^ k := by
      refine Units.ext ?_
      show c / (d * n ^ k) = ((Units.mk0 c hc0 * (Units.mk0 d hd)⁻¹
        * ((Units.mk0 n hn)⁻¹) ^ k : Lˣ) : L)
      push_cast
      rw [Units.val_mk0, Units.val_mk0, Units.val_mk0, inv_pow]
      field_simp
    rw [hEq, valAdd_mul, valAdd_mul, valAdd_pow, valAdd_inv, valAdd_inv, hd0, mul_neg]
    have hh := hc hc0
    omega

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**短 Weierstrass 形の係数が整で `v(Δ) = 12v(n)` なら `minDeltaExp = 0`**
（第 1430）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1425 の**不等式版**である——良い素点で `minDeltaExp p W = 0` を出すのに使う。
★`toShortNF`（`u = 1`）で `a₄ = −c₄/48`・`a₆ = −c₆/864` にし、`u = n` で割る。 -/
theorem minDeltaExp_eq_zero_of_short_integral (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] {n : L} (hn : n ≠ 0)
    (hm4 : W.c₄ / (48 * n ^ 4) ∈ primeSubring p)
    (hm6 : W.c₆ / (864 * n ^ 6) ∈ primeSubring p)
    (hΔ : valAdd p (Units.mk0 W.Δ W.isUnit_Δ.ne_zero)
      = 12 * valAdd p (Units.mk0 n hn)) :
    minDeltaExp p W = 0 := by
  letI : Invertible (2 : L) := invertibleOfNonzero two_ne_zero
  letI : Invertible (3 : L) := invertibleOfNonzero three_ne_zero
  set C0 : WeierstrassCurve.VariableChange L := W.toShortNF with hC0
  haveI hshort : (C0 • W).IsShortNF := W.toShortNF_spec
  have hu0 : C0.u = 1 := by
    rw [hC0]
    simp [WeierstrassCurve.toShortNF, WeierstrassCurve.toCharNeTwoNF,
      WeierstrassCurve.VariableChange.mul_def]
  have hc4W0 : (C0 • W).c₄ = W.c₄ := by
    rw [WeierstrassCurve.variableChange_c₄, hu0]; simp
  have hc6W0 : (C0 • W).c₆ = W.c₆ := by
    rw [WeierstrassCurve.variableChange_c₆, hu0]; simp
  have ha4W0 : (C0 • W).a₄ = -W.c₄ / 48 := by
    have h := WeierstrassCurve.c₄_of_isShortNF (C0 • W)
    rw [hc4W0] at h
    linear_combination h / 48
  have ha6W0 : (C0 • W).a₆ = -W.c₆ / 864 := by
    have h := WeierstrassCurve.c₆_of_isShortNF (C0 • W)
    rw [hc6W0] at h
    linear_combination h / 864
  obtain ⟨hA1, hA2, hA3, hA4, hA6⟩ := short_scale_coeffs (C0 • W) hn
  set Cn : WeierstrassCurve.VariableChange L :=
    { u := Units.mk0 n hn, r := 0, s := 0, t := 0 } with hCn
  have hEqA4 : (Cn • (C0 • W)).a₄ = -(W.c₄ / (48 * n ^ 4)) := by
    rw [hCn] at hA4 ⊢
    rw [hA4, ha4W0, inv_pow]
    field_simp
  have hEqA6 : (Cn • (C0 • W)).a₆ = -(W.c₆ / (864 * n ^ 6)) := by
    rw [hCn] at hA6 ⊢
    rw [hA6, ha6W0, inv_pow]
    field_simp
  haveI hint : WeierstrassCurve.IsIntegral (primeSubring p) ((Cn * C0) • W) := by
    rw [mul_smul]
    refine isIntegral_of_mem _ _ ?_ ?_ ?_ ?_ ?_
    · rw [hCn] at hA1 ⊢; rw [hA1]; exact zero_mem _
    · rw [hCn] at hA2 ⊢; rw [hA2]; exact zero_mem _
    · rw [hCn] at hA3 ⊢; rw [hA3]; exact zero_mem _
    · rw [hEqA4]; exact neg_mem hm4
    · rw [hEqA6]; exact neg_mem hm6
  have hΔW : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  have hle := minDeltaExp_le_of_isIntegral p W hΔW (Cn * C0) hint
  have huu : (Cn * C0).u = Units.mk0 n hn := by
    rw [WeierstrassCurve.VariableChange.mul_def]
    simp [hu0, hCn]
  have hΔ1 : valAdd p (Units.mk0 (((Cn * C0) • W).Δ)
      (variableChange_Delta_ne_zero W hΔW (Cn * C0))) = 0 := by
    have hU : Units.mk0 (((Cn * C0) • W).Δ) (variableChange_Delta_ne_zero W hΔW (Cn * C0))
        = ((Cn * C0).u)⁻¹ ^ 12 * Units.mk0 W.Δ hΔW := by
      refine Units.ext ?_
      show ((Cn * C0) • W).Δ = _
      rw [WeierstrassCurve.variableChange_Δ]
      push_cast
      simp
    rw [hU, valAdd_mul, valAdd_pow, valAdd_inv, huu, hΔ]
    ring
  rw [hΔ1] at hle
  exact le_antisymm hle (minDeltaExp_nonneg p W)

/-- ★★★★★★★★**同じ仮定から半安定**（第 1430）。 -/
theorem semistableAt_of_short_integral (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] {n : L} (hn : n ≠ 0)
    (hm4 : W.c₄ / (48 * n ^ 4) ∈ primeSubring p)
    (hm6 : W.c₆ / (864 * n ^ 6) ∈ primeSubring p)
    (hΔ : valAdd p (Units.mk0 W.Δ W.isUnit_Δ.ne_zero)
      = 12 * valAdd p (Units.mk0 n hn)) :
    SemistableAt p W :=
  Or.inl (minDeltaExp_eq_zero_of_short_integral p W hn hm4 hm6 hΔ)

/-! ## ★出典の紐付け(`.src`) -/

def div_mem_primeSubring_of_valAdd_le.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v(c) ≥ k·v(n) なら c/(d·n^k) は付値環の元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_zero_of_short_integral.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(短 Weierstrass 形の係数が整で v(Δ) = 12v(n) なら minDeltaExp = 0)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_of_short_integral.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同じ仮定から半安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_zero_of_short_integral.needs : List ProofObligation :=
  [ .citation "[mathlib]" "WeierstrassCurve.toShortNF(短 Weierstrass 形への変数変換)"
      (.inMathlib "WeierstrassCurve.toShortNF") 1,
    .citation "[ABC3]" "minDeltaExp_le_of_isIntegral(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_le_of_isIntegral") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1430）**——第 1425 の**不等式版**である。" ++
       "☆良い素点で `p ∣ l` の場合、判別式の恒等式から `v(Δ(E′)) = 12S` は出るので、" ++
       "`v(c₄(E′)) ≥ 4S`・`v(c₆(E′)) ≥ 6S` さえ言えれば本補題で " ++
       "`minDeltaExp p E′ = 0` になる。" ++
       "★その 2 つは **`jExp p E′ ≥ 0`**（`j` の整性）から出る" ++
       "——`3v(c₄) ≥ v(Δ)` と `c₄³ − c₆² = 1728Δ`。") 17 ]

end ABC3.Found.GaloisRep
