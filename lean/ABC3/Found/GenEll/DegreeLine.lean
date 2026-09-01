/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.DegreeBound
import ABC3.Found.GenEll.GalRepContinuity
import ABC3.Meta.Claim

/-!
# 第 1218 ブロック —— **安定直線の定義体は次数 `≤ l−1`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——第 1209 の義務の解消

第 1209 で「残る義務は `[M:L] < l`」となり、第 1214-1217 で道具が揃った:

| 道具 | 内容 | 第 |
|---|---|---|
| 数え上げ | `M →ₐ[L] L̄` が有限集合に単射で入れば `[M:L] ≤ #s` | 1214 |
| 延長 | `φ : M →ₐ[L] L̄` は `Gal(L̄/L)` の元から来る | 1215 |
| 決定 | 不変量が生成元の上で `φ` を決めれば次数が抑えられる | 1216 |
| 生成集合 | `L(T)` の中で `T` は `Algebra.adjoin` としても生成 | 1217 |

★本ブロックはこれを `c(φ)` で繋ぐ——`σ Q = c • Q` の `c ∈ ZMod l` である。

☆`c` が同じなら `σ₁ Q = σ₂ Q`、したがって `Q = (x, y)` の座標の上で `σ₁ = σ₂`、
すなわち `φ₁ = φ₂`。★`c ≠ 0` は `σ Q ≠ 0` から出る。

★★★したがって `[L(x_Q, y_Q) : L] ≤ l − 1 < l` であり、
第 1209 の `l ∤ e(P|p)` が使える。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine
open ABC3.Meta

variable {L : Type} [Field L] [CharZero L]

local notation "Lbar" => AlgebraicClosure L

open scoped Classical in
set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**安定直線の定義体は次数 `≤ l−1`**——★**無条件**（第 1218）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`Q` の座標が生成する `M ≔ L(x_Q, y_Q)` について、
`φ : M →ₐ[L] L̄` を `Gal(L̄/L)` へ延ばし（第 1215）、
`σ Q = c • Q` の `c ∈ ZMod l` を対応させる。

★`c` が同じなら `σ₁ Q = σ₂ Q` なので座標の上で `σ₁ = σ₂`、すなわち `φ₁ = φ₂`。
☆`c ≠ 0` は `σ Q ≠ 0` から出る。したがって `[M:L] ≤ l − 1`。

★★★これで第 1209 が残した義務 `[M:L] < l` が解ける。 -/
theorem finrank_adjoin_ptCoordSet_le (W : WeierstrassCurve L)
    [(W.baseChange Lbar).IsElliptic]
    {l : ℕ} (hl : l.Prime) (Q : (W.baseChange Lbar).toAffine.Point) (hQ : addOrderOf Q = l)
    (hst : ∀ σ : Lbar ≃ₐ[L] Lbar, ∃ k : ℕ, galPoint W σ Q = k • Q) :
    Module.finrank L (IntermediateField.adjoin L (ptCoordSet W Q)) ≤ l - 1 := by
  classical
  haveI : NeZero l := ⟨hl.ne_zero⟩
  have hlQ : l • Q = 0 := by rw [← hQ]; exact addOrderOf_nsmul_eq_zero Q
  have hQne : Q ≠ 0 := by
    intro h0
    rw [h0, addOrderOf_zero] at hQ
    exact absurd hQ.symm hl.ne_one
  have hmodQ : ∀ a : ℕ, a • Q = (a % l) • Q := by
    intro a
    conv_lhs => rw [(Nat.div_add_mod' a l).symm]
    rw [add_smul, mul_smul, hlQ, nsmul_zero, zero_add]
  have hgalne : ∀ σ : Lbar ≃ₐ[L] Lbar, galPoint W σ Q ≠ 0 := by
    intro σ h0
    have h1 : galPoint W σ⁻¹ (galPoint W σ Q) = Q := by
      rw [← galPoint_mul, inv_mul_cancel, galPoint_one]
    rw [h0, map_zero] at h1
    exact hQne h1.symm
  obtain ⟨x, y, hns, rfl⟩ : ∃ x y hns, Q = Point.some x y hns := by
    cases Q with
    | zero => exact absurd rfl hQne
    | some x y hns => exact ⟨x, y, hns, rfl⟩
  have hTalg : ∀ z ∈ ptCoordSet W (Point.some x y hns), IsAlgebraic L z :=
    fun z _ => Algebra.IsAlgebraic.isAlgebraic z
  haveI : Finite (ptCoordSet W (Point.some x y hns) : Set Lbar) :=
    (finite_ptCoordSet W (Point.some x y hns)).to_subtype
  haveI : FiniteDimensional L
      (IntermediateField.adjoin L (ptCoordSet W (Point.some x y hns))) :=
    IntermediateField.finiteDimensional_adjoin (fun z hz => (hTalg z hz).isIntegral)
  set M := IntermediateField.adjoin L (ptCoordSet W (Point.some x y hns)) with hMdef
  have hext : ∀ φ : M →ₐ[L] Lbar,
      ∃ σ : Lbar ≃ₐ[L] Lbar, ∀ w : M, σ (algebraMap M Lbar w) = φ w :=
    fun φ => exists_algEquiv_extend φ
  set σOf : (M →ₐ[L] Lbar) → (Lbar ≃ₐ[L] Lbar) := fun φ => (hext φ).choose with hσOf
  set kOf : (M →ₐ[L] Lbar) → ℕ := fun φ => (hst (σOf φ)).choose with hkOf
  have hσspec : ∀ (φ : M →ₐ[L] Lbar) (w : M), (σOf φ) (algebraMap M Lbar w) = φ w :=
    fun φ => (hext φ).choose_spec
  have hkspec : ∀ φ : M →ₐ[L] Lbar,
      galPoint W (σOf φ) (Point.some x y hns) = (kOf φ) • (Point.some x y hns) :=
    fun φ => (hst (σOf φ)).choose_spec
  refine le_trans (finrank_le_of_determines_on_generators
    (((↑) : M → Lbar) ⁻¹' (ptCoordSet W (Point.some x y hns)))
    (algebra_adjoin_preimage_eq_top _ hTalg)
    (fun φ => ((kOf φ : ℕ) : ZMod l)) ?_ (Finset.univ.erase (0 : ZMod l)) ?_) ?_
  · intro φ₁ φ₂ hc z hz
    have hmod : (kOf φ₁) % l = (kOf φ₂) % l :=
      (ZMod.natCast_eq_natCast_iff' _ _ _).mp hc
    have hkQ : (kOf φ₁) • (Point.some x y hns) = (kOf φ₂) • (Point.some x y hns) := by
      rw [hmodQ (kOf φ₁), hmodQ (kOf φ₂), hmod]
    have hQeq : galPoint W (σOf φ₁) (Point.some x y hns)
        = galPoint W (σOf φ₂) (Point.some x y hns) := by
      rw [hkspec φ₁, hkspec φ₂, hkQ]
    rw [galPoint, galPoint, Point.map_some, Point.map_some] at hQeq
    injection hQeq with hx hy
    have hzmem : (z : Lbar) = x ∨ (z : Lbar) = y := by
      have hz' : (z : Lbar) ∈ ({x, y} : Set Lbar) := hz
      simpa using hz'
    show φ₁ z = φ₂ z
    rw [← hσspec φ₁ z, ← hσspec φ₂ z]
    show (σOf φ₁) (z : Lbar) = (σOf φ₂) (z : Lbar)
    rcases hzmem with h | h
    · rw [h]; exact hx
    · rw [h]; exact hy
  · intro φ
    refine Finset.mem_erase.2 ⟨?_, Finset.mem_univ _⟩
    intro h0
    have h0' : ((kOf φ : ℕ) : ZMod l) = ((0 : ℕ) : ZMod l) := by simpa using h0
    have hmod0 : (kOf φ) % l = 0 % l := (ZMod.natCast_eq_natCast_iff' _ _ _).mp h0'
    have hdvd : l ∣ kOf φ := Nat.dvd_of_mod_eq_zero (by simpa using hmod0)
    obtain ⟨m, hm⟩ := hdvd
    have hzero : (kOf φ) • (Point.some x y hns) = 0 := by
      rw [hm, mul_comm, mul_smul, hlQ, nsmul_zero]
    exact hgalne (σOf φ) (by rw [hkspec φ, hzero])
  · have hcard : (Finset.univ.erase (0 : ZMod l)).card = l - 1 := by
      rw [Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ, ZMod.card]
    exact le_of_eq hcard

/-! ## ★出典の紐付け(`.src`) -/

def finrank_adjoin_ptCoordSet_le.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(安定直線の定義体は次数 ≤ l−1。★無条件)",
    sectionId := "genell-lemma-3-5" }

def finrank_adjoin_ptCoordSet_le.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_algEquiv_extend(第 1215、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_algEquiv_extend") 1,
    .citation "[ABC3]" "finrank_le_of_determines_on_generators(第 1216、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.finrank_le_of_determines_on_generators") 1,
    .citation "[ABC3]" "algebra_adjoin_preimage_eq_top(第 1217、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.algebra_adjoin_preimage_eq_top") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1218）**——第 1209 が残した義務 `[M:L] < l` を解く。" ++
       "☆`φ` に `σ Q = c • Q` の `c ∈ ZMod l` を対応させると、" ++
       "`c` が同じなら座標の上で `σ₁ = σ₂` なので `φ₁ = φ₂`、" ++
       "`c ≠ 0` は `σ Q ≠ 0` から出る。") 3 ]

end ABC3.Found.GenEll
