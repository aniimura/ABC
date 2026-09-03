/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.SemistableFromC4
import Mathlib.AlgebraicGeometry.EllipticCurve.NormalForms
import ABC3.Meta.Claim

/-!
# 第 1425 ブロック —— **`v(c₄) = 4k`・`v(c₆) = 6k` なら半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——`p ∣ l` の `μ_l` 側の鍵

第 1424 で残ったのは **`p ∣ l` かつ核が `μ_l` 型**の場合である。
そこでは `K` の水準の恒等式（第 1129-1138、`hlu` なし）から

    `valAdd p (c₄ E′) = 4·valAdd p (l)`、`valAdd p (c₆ E′) = 6·valAdd p (l)`

が出るが、`semistableAt_of_c4_valAdd_zero`（第 1322）は `v_p(c₄) = 0` を要求する。
☆`blocked-leaves.json` の `pDivLMuGapNamed2026_09_02` では
「Kraus か弱近似が要る」と見立てていた。

★★★**どちらも要らなかった**——`2, 3` が `p` で単元なら（`l ≥ 5` なので `p ∤ 6`）
**短 Weierstrass 形**（mathlib の `WeierstrassCurve.toShortNF`）が使える:

1. `C₀ := W.toShortNF` は `u = 1` なので `c₄`・`c₆` を変えない。
   `C₀ • W` は `y² = x³ + a₄x + a₆` の形で、`a₄ = −c₄/48`、`a₆ = −c₆/864`。
2. `u = n` で割ると `a₄ ↦ −c₄/(48n⁴)`、`a₆ ↦ −c₆/(864n⁶)`——
   仮定 `v(c₄) = 4v(n)`・`v(c₆) = 6v(n)` と `v(48) = v(864) = 0` から**どちらも付値 `0`**。
3. したがってその模型は `p` 上整で `v_p(c₄) = 0`——第 1322 が効く。

☆`r`・`s`・`t` は mathlib の `toShortNF` が全部作ってくれるので、
**弱近似も Kraus も要らない**。

★注意: `v_p(48) = v_p(864) = 0` は `p ∤ 6` を意味するので、`p ∣ l` の場合は
`l ≥ 5`（`l ≠ 2, 3`）が要る。☆`l = 3` は別扱いである。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ☆道具 2 本 -/

/-- ☆`valAdd = 0` なら付値環の元。 -/
theorem mem_primeSubring_of_valAdd_eq_zero (p : HeightOneSpectrum (𝓞 L))
    {z : L} (hz : z ≠ 0) (h : valAdd p (Units.mk0 z hz) = 0) : z ∈ primeSubring p :=
  (mem_primeSubring_iff p z).2 ((valAdd_nonneg_iff p (Units.mk0 z hz)).1 (le_of_eq h.symm))

/-- ☆`v(c) = k·v(n)`・`v(d) = 0` なら `v(c/(d·nᵏ)) = 0`。 -/
theorem valAdd_div_mul_pow_eq_zero (p : HeightOneSpectrum (𝓞 L))
    {c n d : L} (hc : c ≠ 0) (hn : n ≠ 0) (hd : d ≠ 0) (k : ℕ)
    (hk : valAdd p (Units.mk0 c hc) = (k : ℤ) * valAdd p (Units.mk0 n hn))
    (hd0 : valAdd p (Units.mk0 d hd) = 0) :
    valAdd p (Units.mk0 (c / (d * n ^ k))
      (div_ne_zero hc (mul_ne_zero hd (pow_ne_zero k hn)))) = 0 := by
  have hEq : Units.mk0 (c / (d * n ^ k))
      (div_ne_zero hc (mul_ne_zero hd (pow_ne_zero k hn)))
      = Units.mk0 c hc * (Units.mk0 d hd)⁻¹ * ((Units.mk0 n hn)⁻¹) ^ k := by
    refine Units.ext ?_
    show c / (d * n ^ k) = ((Units.mk0 c hc * (Units.mk0 d hd)⁻¹
      * ((Units.mk0 n hn)⁻¹) ^ k : Lˣ) : L)
    push_cast
    rw [Units.val_mk0, Units.val_mk0, Units.val_mk0, inv_pow]
    field_simp
  rw [hEq, valAdd_mul, valAdd_mul, valAdd_pow, valAdd_inv, valAdd_inv, hd0, hk]
  ring

omit [NumberField L] in
/-- ☆短 Weierstrass 形を `u = n` で割ったときの係数。 -/
theorem short_scale_coeffs (V : WeierstrassCurve L) [V.IsShortNF] {n : L} (hn : n ≠ 0) :
    (({ u := Units.mk0 n hn, r := 0, s := 0, t := 0 }
        : WeierstrassCurve.VariableChange L) • V).a₁ = 0
      ∧ (({ u := Units.mk0 n hn, r := 0, s := 0, t := 0 }
        : WeierstrassCurve.VariableChange L) • V).a₂ = 0
      ∧ (({ u := Units.mk0 n hn, r := 0, s := 0, t := 0 }
        : WeierstrassCurve.VariableChange L) • V).a₃ = 0
      ∧ (({ u := Units.mk0 n hn, r := 0, s := 0, t := 0 }
        : WeierstrassCurve.VariableChange L) • V).a₄ = n⁻¹ ^ 4 * V.a₄
      ∧ (({ u := Units.mk0 n hn, r := 0, s := 0, t := 0 }
        : WeierstrassCurve.VariableChange L) • V).a₆ = n⁻¹ ^ 6 * V.a₆ := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;>
    simp [WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₂,
      WeierstrassCurve.variableChange_a₃, WeierstrassCurve.variableChange_a₄,
      WeierstrassCurve.variableChange_a₆, WeierstrassCurve.a₂_of_isShortNF]

/-! ## ★★★★★★★★★★★★★★★★★★★★本体 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`v_p(c₄) = 4·v_p(n)`・`v_p(c₆) = 6·v_p(n)`・`p ∤ 6` なら半安定**（第 1425）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1322（`c₄` が単元の整モデルは半安定）の**`n` 倍ずれた版**である。

★★★道は短 Weierstrass 形 1 本:
`W.toShortNF` は `u = 1` なので `c₄`・`c₆` を変えず、`a₄ = −c₄/48`・`a₆ = −c₆/864` にする。
そこから `u = n` で割れば `a₄ = −c₄/(48n⁴)`・`a₆ = −c₆/(864n⁶)` で、
仮定よりどちらも付値 `0`——`p` 上整で `v_p(c₄) = 0` になる。

☆`v_p(48) = v_p(864) = 0` は `p ∤ 6` を言っている。 -/
theorem semistableAt_of_valAdd_c4_c6_scaled (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic]
    {n : L} (hn : n ≠ 0)
    (hc4ne : W.c₄ ≠ 0) (hc6ne : W.c₆ ≠ 0)
    (h4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 4 * valAdd p (Units.mk0 n hn))
    (h6 : valAdd p (Units.mk0 W.c₆ hc6ne) = 6 * valAdd p (Units.mk0 n hn))
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0) :
    SemistableAt p W := by
  letI : Invertible (2 : L) := invertibleOfNonzero two_ne_zero
  letI : Invertible (3 : L) := invertibleOfNonzero three_ne_zero
  -- ★段 1: 短 Weierstrass 形（`u = 1` なので `c₄`・`c₆` は変わらない）
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
  -- ★段 2: `u = n` で割る
  obtain ⟨hA1, hA2, hA3, hA4, hA6⟩ := short_scale_coeffs (C0 • W) hn
  set Cn : WeierstrassCurve.VariableChange L :=
    { u := Units.mk0 n hn, r := 0, s := 0, t := 0 } with hCn
  have hz4 : valAdd p (Units.mk0 (W.c₄ / (48 * n ^ 4))
      (div_ne_zero hc4ne (mul_ne_zero (by norm_num) (pow_ne_zero 4 hn)))) = 0 :=
    valAdd_div_mul_pow_eq_zero p hc4ne hn (by norm_num) 4 (by exact_mod_cast h4) h48
  have hz6 : valAdd p (Units.mk0 (W.c₆ / (864 * n ^ 6))
      (div_ne_zero hc6ne (mul_ne_zero (by norm_num) (pow_ne_zero 6 hn)))) = 0 :=
    valAdd_div_mul_pow_eq_zero p hc6ne hn (by norm_num) 6 (by exact_mod_cast h6) h864
  have hm4 : W.c₄ / (48 * n ^ 4) ∈ primeSubring p :=
    mem_primeSubring_of_valAdd_eq_zero p _ hz4
  have hm6 : W.c₆ / (864 * n ^ 6) ∈ primeSubring p :=
    mem_primeSubring_of_valAdd_eq_zero p _ hz6
  have hEqA4 : (Cn • (C0 • W)).a₄ = -(W.c₄ / (48 * n ^ 4)) := by
    rw [hCn] at hA4 ⊢
    rw [hA4, ha4W0, inv_pow]
    field_simp
  have hEqA6 : (Cn • (C0 • W)).a₆ = -(W.c₆ / (864 * n ^ 6)) := by
    rw [hCn] at hA6 ⊢
    rw [hA6, ha6W0, inv_pow]
    field_simp
  -- ★段 3: 整性
  haveI hint : WeierstrassCurve.IsIntegral (primeSubring p) (Cn • (C0 • W)) := by
    refine isIntegral_of_mem _ _ ?_ ?_ ?_ ?_ ?_
    · rw [hCn] at hA1 ⊢; rw [hA1]; exact zero_mem _
    · rw [hCn] at hA2 ⊢; rw [hA2]; exact zero_mem _
    · rw [hCn] at hA3 ⊢; rw [hA3]; exact zero_mem _
    · rw [hEqA4]; exact neg_mem hm4
    · rw [hEqA6]; exact neg_mem hm6
  -- ★段 4: `v_p(c₄) = 0`
  have hΔW : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  have hΔ1 : (Cn • (C0 • W)).Δ ≠ 0 :=
    variableChange_Delta_ne_zero _ (variableChange_Delta_ne_zero W hΔW C0) Cn
  have hc4Eq : (Cn • (C0 • W)).c₄ = W.c₄ / (1 * n ^ 4) := by
    rw [WeierstrassCurve.variableChange_c₄, hc4W0, hCn]
    simp only [Units.val_inv_eq_inv_val, Units.val_mk0, inv_pow, one_mul]
    field_simp
  have hc4ne1 : (Cn • (C0 • W)).c₄ ≠ 0 := by
    rw [hc4Eq]
    exact div_ne_zero hc4ne (mul_ne_zero one_ne_zero (pow_ne_zero 4 hn))
  have hone : valAdd p (Units.mk0 (1 : L) one_ne_zero) = 0 := by
    rw [show Units.mk0 (1 : L) one_ne_zero = 1 from Units.ext rfl, valAdd_one]
  have hval1 : valAdd p (Units.mk0 ((Cn • (C0 • W)).c₄) hc4ne1) = 0 := by
    have hz : valAdd p (Units.mk0 (W.c₄ / (1 * n ^ 4))
        (div_ne_zero hc4ne (mul_ne_zero one_ne_zero (pow_ne_zero 4 hn)))) = 0 :=
      valAdd_div_mul_pow_eq_zero p hc4ne hn one_ne_zero 4 (by exact_mod_cast h4) hone
    have hU : Units.mk0 ((Cn • (C0 • W)).c₄) hc4ne1
        = Units.mk0 (W.c₄ / (1 * n ^ 4))
          (div_ne_zero hc4ne (mul_ne_zero one_ne_zero (pow_ne_zero 4 hn))) := by
      refine Units.ext ?_
      show (Cn • (C0 • W)).c₄ = W.c₄ / (1 * n ^ 4)
      exact hc4Eq
    rw [hU, hz]
  -- ★段 5: 第 1322 ＋ 変数変換で戻す
  have hss1 := semistableAt_of_c4_valAdd_zero p (Cn • (C0 • W)) hΔ1 hc4ne1 hval1
  have hfin := semistableAt_variableChange p _ (C0⁻¹ * Cn⁻¹) hss1
  rwa [show (C0⁻¹ * Cn⁻¹) • (Cn • (C0 • W)) = W by
    rw [mul_smul, inv_smul_smul, inv_smul_smul]] at hfin

/-! ## ★出典の紐付け(`.src`) -/

def mem_primeSubring_of_valAdd_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(valAdd = 0 なら付値環の元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_div_mul_pow_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v(c) = k·v(n)・v(d) = 0 なら v(c/(d·n^k)) = 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def short_scale_coeffs.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(短 Weierstrass 形を u = n で割ったときの係数。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_of_valAdd_c4_c6_scaled.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v(c₄) = 4v(n)・v(c₆) = 6v(n)・p ∤ 6 なら半安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_of_valAdd_c4_c6_scaled.needs : List ProofObligation :=
  [ .citation "[mathlib]" "WeierstrassCurve.toShortNF(短 Weierstrass 形への変数変換)"
      (.inMathlib "WeierstrassCurve.toShortNF") 1,
    .citation "[ABC3]" "semistableAt_of_c4_valAdd_zero(第 1322、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_of_c4_valAdd_zero") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1425）**——`p ∣ l` の `μ_l` 側で要ると見ていた" ++
       "**Kraus も弱近似も要らなかった**。" ++
       "☆`2, 3` が `p` で単元なら mathlib の `toShortNF` が `r, s, t` を全部作ってくれる。" ++
       "★`toShortNF` は `u = 1` なので `c₄`・`c₆` を変えず、" ++
       "`a₄ = −c₄/48`・`a₆ = −c₆/864` にする。そこを `u = n` で割れば" ++
       "仮定 `v(c₄) = 4v(n)`・`v(c₆) = 6v(n)` から係数の付値がちょうど `0` になる。" ++
       "☆制約は `p ∤ 6`——`p ∣ l` の場合は `l ≥ 5` を意味する。") 17 ]

end ABC3.Found.GaloisRep
