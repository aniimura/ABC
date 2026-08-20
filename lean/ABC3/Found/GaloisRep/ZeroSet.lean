import ABC3.Found.GaloisRep.Coprime23

/-!
# Galois (G1) 第 60 ブロック —— **★★★★★★零集合は位数の倍数全体**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★Ward の定理が効く

点 `P = (x,y)` に対し `v(j) := normEDS ψ Ψ₃(x) preΨ₄(x) j`(`ψ = y − negY`)とおくと

    v(j)² = ΨSq_j(x)

である。★★Ward の定理(第 58)から、`v(K) = 0` なら

    v(K+r)·v(K−r) = v(K+1)·v(K−1)·v(r)²

がすべての `r` で成り立つ。★★★`r = m − K` を入れると **降下**

    v(m) = 0  ⟹  v(m−K) = 0        (μ := v(K+1)v(K−1) ≠ 0 のとき)

が出て、`K` が最小なら `K ∣ m` になる。

## ★★★★`μ ≠ 0`——連続 2 項は消えない

`v(K) = v(K+1) = 0` を仮定すると `v(K+r)v(K−r) = 0` (∀r) となり、
★`K ≥ 4` では `EllAt (K−1) (K−2)` から `v(K−3) = 0` が出て最小性に反する。
★★残る `K = 2, 3` は:

| `K` | 潰し方 |
|---|---|
| 2 | ★第 59 ブロック(`Ψ₂Sq` と `Ψ₃` の互いに素性) |
| 3 | ★★**乗法公式**——位数 3 の点では `preΨ₄(x) = −ψ⁴ ≠ 0` |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `vSeq` / `vSeq_sq` | ★点での分点列と `v² = ΨSq` |
| `vSeq_ward` | ★★★Ward の帰結 |
| `preP4_of_order_three` | ★★★位数 3 の点で `preΨ₄(x) = −ψ⁴` |
| `not_two_consec` | ★★★★★連続 2 項は消えない |
| `dvd_of_PSq_root` | ★★★★★★**零集合は位数の倍数全体** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-! ## ★点での分点列 -/

theorem preP_eval (j : ℤ) (x : F) :
    (W.preΨ j).eval x
      = preNormEDS (W.Ψ₂Sq.eval x ^ 2) (W.Ψ₃.eval x) (W.preΨ₄.eval x) j := by
  simp [WeierstrassCurve.preΨ, ← Polynomial.coe_evalRingHom]

/-- ★**点での分点列** `v(j) = ψ_j(P)`。 -/
noncomputable def vSeq (W : WeierstrassCurve F) (x y : F) (j : ℤ) : F :=
  normEDS (y - W.toAffine.negY x y) (W.Ψ₃.eval x) (W.preΨ₄.eval x) j

/-- ★★`v(j)² = ΨSq_j(x)`。 -/
theorem vSeq_sq {x y : F} (h : W.toAffine.Equation x y) (j : ℤ) :
    vSeq W x y j ^ 2 = (W.ΨSq j).eval x := by
  have hpsi : (y - W.toAffine.negY x y) ^ 2 = W.Ψ₂Sq.eval x := (psi2Sq_eval W h).symm
  simp only [vSeq, normEDS]
  rw [PSq_eval, preP_eval, ← hpsi]
  by_cases hj : Even j
  · rw [if_pos hj, if_pos hj]
    rw [show ((y - W.toAffine.negY x y) ^ 2) ^ 2 = (y - W.toAffine.negY x y) ^ 4 by ring]
    ring
  · rw [if_neg hj, if_neg hj]
    rw [show ((y - W.toAffine.negY x y) ^ 2) ^ 2 = (y - W.toAffine.negY x y) ^ 4 by ring]
    ring

theorem vSeq_one (x y : F) : vSeq W x y 1 = 1 := by simp [vSeq]
theorem vSeq_two (x y : F) : vSeq W x y 2 = y - W.toAffine.negY x y := by simp [vSeq]
theorem vSeq_three (x y : F) : vSeq W x y 3 = W.Ψ₃.eval x := by simp [vSeq]
theorem vSeq_four (x y : F) :
    vSeq W x y 4 = W.preΨ₄.eval x * (y - W.toAffine.negY x y) := by simp [vSeq]

theorem vSeq_ell (x y : F) (m n : ℤ) :
    vSeq W x y (m + n) * vSeq W x y (m - n)
      = vSeq W x y (m + 1) * vSeq W x y (m - 1) * vSeq W x y n ^ 2
        - vSeq W x y (n + 1) * vSeq W x y (n - 1) * vSeq W x y m ^ 2 :=
  ell _ _ _ m n

/-- ★★★**Ward の帰結**——`v K = 0` なら `v(K+r)v(K−r) = v(K+1)v(K−1)v(r)²`。 -/
theorem vSeq_ward (x y : F) (K : ℤ) (hzero : vSeq W x y K = 0) (r : ℤ) :
    vSeq W x y (K + r) * vSeq W x y (K - r)
      = vSeq W x y (K + 1) * vSeq W x y (K - 1) * vSeq W x y r ^ 2 := by
  have hw := vSeq_ell W x y K r
  rw [hw, hzero]; ring

theorem vSeq_eq_zero_iff {x y : F} (h : W.toAffine.Equation x y) (j : ℤ) :
    vSeq W x y j = 0 ↔ (W.ΨSq j).eval x = 0 := by
  rw [← vSeq_sq W h j]
  exact (pow_eq_zero_iff two_ne_zero).symm

/-! ## ★★★位数 3 の点 -/

/-- ★★★**位数 3 の点では `preΨ₄(x) = −ψ⁴`**(特に `≠ 0`)。

★乗法公式(第 52)の `n = 2` を `2P = −P` に当てるだけである。 -/
theorem preP4_of_order_three {x y : F} (h : W.toAffine.Nonsingular x y)
    (hpsi : y - W.toAffine.negY x y ≠ 0)
    (h3 : (3 : ℕ) • (Point.some x y h) = 0) :
    W.preΨ₄.eval x = -(y - W.toAffine.negY x y) ^ 4 := by
  have hA2 : (W.ΨSq ((2 : ℕ) : ℤ)).eval x ≠ 0 := by
    rw [show (((2 : ℕ)) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two, psi2Sq_eval W h.left]
    exact pow_ne_zero 2 hpsi
  obtain ⟨x2, y2, h2', hP2, hx2, hy2⟩ := mulOK_two W h hA2
  have hneg : (2 : ℕ) • (Point.some x y h) = -(Point.some x y h) := by
    have hh : (3 : ℕ) • (Point.some x y h) = (2 : ℕ) • (Point.some x y h) + Point.some x y h := by
      rw [show (3 : ℕ) = 2 + 1 from rfl, succ_nsmul]
    rw [hh] at h3
    exact eq_neg_of_add_eq_zero_left h3
  rw [hP2, Point.neg_some] at hneg
  obtain ⟨rfl, rfl⟩ := Point.some.inj hneg
  rw [show (((2 : ℕ)) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two, psi2Sq_eval W h.left,
    show (2 : ℤ) * 2 = 4 by norm_num, preP_four W, WeierstrassCurve.Affine.negY_negY] at hy2
  refine mul_right_cancel₀ hpsi ?_
  linear_combination -hy2

/-! ## ★★★★★連続 2 項は消えない -/

/-- ★★★★★**連続 2 項は消えない**。 -/
theorem not_two_consec {x y : F} (hΔ : W.Δ ≠ 0) (h : W.toAffine.Nonsingular x y)
    (K : ℤ) (hK : 1 ≤ K)
    (hzero : vSeq W x y K = 0)
    (hmin : ∀ j : ℤ, 1 ≤ j → j < K → vSeq W x y j ≠ 0)
    (hcons : vSeq W x y (K + 1) = 0) : False := by
  have hK2 : 2 ≤ K := by
    rcases eq_or_lt_of_le hK with h1 | h1
    · exact absurd (h1 ▸ hzero) (by rw [vSeq_one]; exact one_ne_zero)
    · omega
  have hward : ∀ r : ℤ, vSeq W x y (K + r) * vSeq W x y (K - r) = 0 := by
    intro r; rw [vSeq_ward W x y K hzero r, hcons]; ring
  by_cases hK4 : 4 ≤ K
  · have h3ne : vSeq W x y 3 ≠ 0 := hmin 3 (by omega) (by omega)
    have h2K3 : vSeq W x y (2 * K - 3) = 0 := by
      have hh := hward (K - 3)
      rw [show K + (K - 3) = 2 * K - 3 by ring, show K - (K - 3) = 3 by ring] at hh
      exact (mul_eq_zero.1 hh).resolve_right h3ne
    have hKm1 : vSeq W x y (K - 1) ≠ 0 := hmin (K - 1) (by omega) (by omega)
    have hKm3 : vSeq W x y (K - 3) ≠ 0 := hmin (K - 3) (by omega) (by omega)
    have hw2 := vSeq_ell W x y (K - 1) (K - 2)
    rw [show K - 1 + (K - 2) = 2 * K - 3 by ring, show K - 1 - (K - 2) = 1 by ring,
      show K - 1 + 1 = K by ring, show K - 2 + 1 = K - 1 by ring,
      show K - 2 - 1 = K - 3 by ring, h2K3, hzero, vSeq_one] at hw2
    apply hKm3
    have hz : vSeq W x y (K - 1) ^ 3 * vSeq W x y (K - 3) = 0 := by linear_combination hw2
    exact (mul_eq_zero.1 hz).resolve_left (pow_ne_zero 3 hKm1)
  · interval_cases K
    · rw [vSeq_two] at hzero
      rw [show (2 : ℤ) + 1 = 3 by norm_num, vSeq_three] at hcons
      have h2 : W.Ψ₂Sq.eval x = 0 := by rw [psi2Sq_eval W h.left, hzero]; ring
      exact hΔ (pow_eq_zero_iff two_ne_zero |>.1 (no_common_root_23 W x h2 hcons))
    · rw [vSeq_three] at hzero
      rw [show (3 : ℤ) + 1 = 4 by norm_num, vSeq_four] at hcons
      have hpsi : y - W.toAffine.negY x y ≠ 0 := by
        have hh := hmin 2 (by norm_num) (by norm_num); rwa [vSeq_two] at hh
      have hp4 : W.preΨ₄.eval x = 0 := (mul_eq_zero.1 hcons).resolve_right hpsi
      have h3P : (3 : ℕ) • (Point.some x y h) = 0 := by
        refine smul_eq_zero_of_min_root' W h 3 (by norm_num) ?_ ?_
        · rw [show ((3 : ℕ) : ℤ) = 3 from rfl, ← vSeq_eq_zero_iff W h.left, vSeq_three]
          exact hzero
        · intro k hk1 hk2 hk3
          exact hmin (k : ℤ) (by exact_mod_cast hk1) (by exact_mod_cast hk2)
            ((vSeq_eq_zero_iff W h.left (k : ℤ)).mpr hk3)
      have hkey := preP4_of_order_three W h hpsi h3P
      rw [hp4] at hkey
      exact hpsi (pow_eq_zero_iff (n := 4) (by norm_num) |>.1 (by linear_combination hkey))

/-! ## ★★★★★★降下 -/

theorem dvd_of_vSeq_zero_nat {x y : F} (hΔ : W.Δ ≠ 0) (h : W.toAffine.Nonsingular x y)
    (K : ℤ) (hK : 1 ≤ K) (hzero : vSeq W x y K = 0)
    (hmin : ∀ j : ℤ, 1 ≤ j → j < K → vSeq W x y j ≠ 0) :
    ∀ m : ℕ, vSeq W x y (m : ℤ) = 0 → K ∣ (m : ℤ) := by
  have hK2 : 2 ≤ K := by
    rcases eq_or_lt_of_le hK with h1 | h1
    · exact absurd (h1 ▸ hzero) (by rw [vSeq_one]; exact one_ne_zero)
    · omega
  have hmu : vSeq W x y (K + 1) * vSeq W x y (K - 1) ≠ 0 :=
    mul_ne_zero (fun hc => not_two_consec W hΔ h K hK hzero hmin hc)
      (hmin (K - 1) (by omega) (by omega))
  obtain ⟨k, hk⟩ : ∃ k : ℕ, (k : ℤ) = K := ⟨K.toNat, Int.toNat_of_nonneg (by omega)⟩
  intro m
  induction m using Nat.strong_induction_on with
  | _ m ih =>
    intro hm
    by_cases hmk : (m : ℤ) < K
    · rcases Nat.eq_zero_or_pos m with rfl | hpos
      · simp
      · exact absurd hm (hmin m (by exact_mod_cast hpos) hmk)
    · have hmk' : K ≤ (m : ℤ) := by omega
      have hd : vSeq W x y ((m : ℤ) - K) = 0 := by
        have hw := vSeq_ward W x y K hzero ((m : ℤ) - K)
        rw [show K + ((m : ℤ) - K) = (m : ℤ) by ring, hm, zero_mul] at hw
        have hsq : vSeq W x y ((m : ℤ) - K) ^ 2 = 0 :=
          (mul_eq_zero.1 hw.symm).resolve_left hmu
        exact pow_eq_zero_iff two_ne_zero |>.1 hsq
      have hkm : k ≤ m := by omega
      have hlt : m - k < m := by omega
      have hcast : ((m - k : ℕ) : ℤ) = (m : ℤ) - K := by rw [Nat.cast_sub hkm, hk]
      have hdv := ih (m - k) hlt (by rw [hcast]; exact hd)
      rw [hcast] at hdv
      obtain ⟨t, ht⟩ := hdv
      exact ⟨t + 1, by rw [mul_add, mul_one, ← ht]; ring⟩

/-- ★★★★★★**零集合は位数の倍数全体**——`ΨSq_m(x) = 0` なら `K ∣ m`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが `Φ_n` と `ΨSq_n` の互いに素性を与える。 -/
theorem dvd_of_PSq_root {x y : F} (hΔ : W.Δ ≠ 0) (h : W.toAffine.Nonsingular x y)
    (K : ℤ) (hK : 1 ≤ K) (hzero : (W.ΨSq K).eval x = 0)
    (hmin : ∀ j : ℤ, 1 ≤ j → j < K → (W.ΨSq j).eval x ≠ 0)
    (m : ℕ) (hm : (W.ΨSq (m : ℤ)).eval x = 0) : K ∣ (m : ℤ) := by
  refine dvd_of_vSeq_zero_nat W hΔ h K hK ((vSeq_eq_zero_iff W h.left K).mpr hzero)
    (fun j hj1 hj2 hc => hmin j hj1 hj2 ((vSeq_eq_zero_iff W h.left j).mp hc)) m ?_
  exact (vSeq_eq_zero_iff W h.left (m : ℤ)).mpr hm

/-! ## ★出典の紐付け(`.src`) -/

def dvd_of_PSq_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(分点列の零集合は位数の倍数全体)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
