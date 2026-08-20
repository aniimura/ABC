import ABC3.Found.GaloisRep.MulPoint

/-!
# Galois (G1) 第 53 ブロック —— **★★★★★★最小の根は位数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★逆向きも出る

第 52 ブロックの対偶は

    n•P = 0  ⟹  ∃ k ≤ n, ΨSq_k(x) = 0

だが、これだけでは `k` が `n` に依らない多項式に収まらない。
★そこで**逆向き**を作る:

    ΨSq_K(x) = 0 が最小  ⟹  K•P = 0

★★これで「最小の `K`」= **点の位数**となり、`K ∣ n` が出る。

## ★★★逆向きの証明(§9-393)

| `K` | 議論 |
|---|---|
| `K = 1` | ★`ΨSq_1 = 1 ≠ 0` なので起きない |
| `K = 2` | ★`Ψ₂Sq(x) = 0 ⟺ y = negY ⟺ 2P = 0` |
| `K ≥ 3` | ★★`preΨ_K(x) = 0` ⟹ `x((K−1)P) = x` ⟹ `(K−1)P = ±P` |

★`(K−1)P = P` は `(K−2)P = 0` を意味するが、最小性に反する ✅

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_PSq_root` | ★★対偶 |
| `min_root_step` / `smul_eq_zero_of_min_root'` | ★★★★★逆向き |
| `exists_divisor_root` | ★★★★★★**`n` の約数 `k` で `ΨSq_k(x) = 0`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★対偶——`n • P = 0` なら `k ≤ n` で `ΨSq_k(x) = 0` となる `k` がある。 -/
theorem exists_PSq_root {x y : F} (h : W.toAffine.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hP : n • (Point.some x y h) = 0) :
    ∃ k : ℕ, 1 ≤ k ∧ k ≤ n ∧ (W.ΨSq (k : ℤ)).eval x = 0 := by
  rcases Classical.em (∃ k : ℕ, 1 ≤ k ∧ k ≤ n ∧ (W.ΨSq (k : ℤ)).eval x = 0) with hyes | hno
  · exact hyes
  · exfalso
    have hne : ∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval x ≠ 0 :=
      fun k hk1 hk2 hk3 => hno ⟨k, hk1, hk2, hk3⟩
    obtain ⟨x', y', h', hEq, -, -⟩ := mulOK_of_ne W h n hn hne
    rw [hP] at hEq
    exact Point.some_ne_zero h' hEq.symm

/-- ★★★★★逆向きの本体(`K ≥ 3`)。 -/
theorem min_root_step {x y : F} (h : W.toAffine.Nonsingular x y) (j : ℕ)
    (hroot : (W.ΨSq (((j + 3 : ℕ)) : ℤ)).eval x = 0)
    (hmin : ∀ k : ℕ, 1 ≤ k → k < j + 3 → (W.ΨSq (k : ℤ)).eval x ≠ 0) :
    (j + 3) • (Point.some x y h) = 0 := by
  have e3 : ((j + 3 : ℕ) : ℤ) = (j : ℤ) + 3 := by push_cast; ring
  have e2 : ((j + 2 : ℕ) : ℤ) = (j : ℤ) + 2 := by push_cast; ring
  rw [e3] at hroot
  have hpre : (W.preΨ ((j : ℤ) + 3)).eval x = 0 := by
    rw [PSq_eval] at hroot
    rcases mul_eq_zero.1 hroot with h1 | h1
    · exact pow_eq_zero_iff two_ne_zero |>.1 h1
    · exfalso
      by_cases hev : Even ((j : ℤ) + 3)
      · rw [if_pos hev] at h1
        refine hmin 2 (by omega) (by omega) ?_
        rw [show ((2 : ℕ) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two]; exact h1
      · rw [if_neg hev] at h1; exact one_ne_zero h1
  obtain ⟨xk, yk, hk', hPk, hxk, -⟩ := mulOK_of_ne W h (j + 2) (by omega)
    (fun k hk1 hk2 => hmin k hk1 (by omega))
  rw [e2] at hxk
  have hA2 : (W.ΨSq ((j : ℤ) + 2)).eval x ≠ 0 := e2 ▸ hmin (j + 2) (by omega) (by omega)
  have hxkx : xk = x := by
    refine mul_right_cancel₀ hA2 ?_
    rw [hxk]
    have hd := x_mul_PSq_sub_Phi W ((j : ℤ) + 2) x
    rw [show (j : ℤ) + 2 + 1 = (j : ℤ) + 3 by ring, hpre, zero_mul, zero_mul] at hd
    linear_combination -hd
  have hsucc : (j + 3) • (Point.some x y h)
      = Point.some xk yk hk' + Point.some x y h := by
    rw [show j + 3 = (j + 2) + 1 from rfl, succ_nsmul, hPk]
  rcases Point.X_eq_iff.mp hxkx with hcase | hcase
  · exfalso
    obtain ⟨x1, y1, h1', hP1, -, -⟩ := mulOK_of_ne W h (j + 1) (by omega)
      (fun k hk1 hk2 => hmin k hk1 (by omega))
    have hz : (j + 1) • (Point.some x y h) = 0 := by
      have hEq : (j + 2) • (Point.some x y h) = Point.some x y h := by rw [hPk, hcase]
      rw [show j + 2 = (j + 1) + 1 from rfl, succ_nsmul] at hEq
      exact add_right_cancel (b := Point.some x y h) (by rw [hEq, zero_add])
    rw [hz] at hP1
    exact Point.some_ne_zero h1' hP1.symm
  · rw [hsucc, hcase, neg_add_cancel]

/-- ★★★★★**逆向き**——`ΨSq_K(x) = 0` が最小なら `K • P = 0`。 -/
theorem smul_eq_zero_of_min_root' {x y : F} (h : W.toAffine.Nonsingular x y) (K : ℕ)
    (hK : 1 ≤ K) (hroot : (W.ΨSq (K : ℤ)).eval x = 0)
    (hmin : ∀ k : ℕ, 1 ≤ k → k < K → (W.ΨSq (k : ℤ)).eval x ≠ 0) :
    K • (Point.some x y h) = 0 := by
  match K, hK with
  | 1, _ => exact absurd hroot (by norm_num [WeierstrassCurve.ΨSq_one])
  | 2, _ =>
      refine (two_smul_eq_zero_iff W h).2 ?_
      rw [show ((2 : ℕ) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two] at hroot
      exact hroot
  | (j + 3), _ => exact min_root_step W h j hroot hmin

/-- ★★★★★★**`n • P = 0` なら `n` の約数 `k` で `ΨSq_k(x) = 0`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★最小の根 `K` は点の位数に一致し、位数は `n` を割る。
★★これで **`k` を `n` の約数に限れる**ので、`(k : F) ≠ 0` が `(n : F) ≠ 0` から従う。 -/
theorem exists_divisor_root {x y : F} (h : W.toAffine.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hP : n • (Point.some x y h) = 0) :
    ∃ k : ℕ, 1 ≤ k ∧ k ∣ n ∧ (W.ΨSq (k : ℤ)).eval x = 0 := by
  classical
  obtain ⟨k0, hk01, -, hk0root⟩ := exists_PSq_root W h n hn hP
  have hex : ∃ k : ℕ, 1 ≤ k ∧ (W.ΨSq (k : ℤ)).eval x = 0 := ⟨k0, hk01, hk0root⟩
  obtain ⟨hK1, hKroot⟩ := Nat.find_spec hex
  set K := Nat.find hex with hKdef
  have hmin : ∀ k : ℕ, 1 ≤ k → k < K → (W.ΨSq (k : ℤ)).eval x ≠ 0 :=
    fun k hk1 hkK hkroot => Nat.find_min hex hkK ⟨hk1, hkroot⟩
  have hKzero : K • (Point.some x y h) = 0 := smul_eq_zero_of_min_root' W h K hK1 hKroot hmin
  refine ⟨K, hK1, ?_, hKroot⟩
  set d := addOrderOf (Point.some x y h) with hd
  have hdK : d ∣ K := addOrderOf_dvd_of_nsmul_eq_zero hKzero
  have hdn : d ∣ n := addOrderOf_dvd_of_nsmul_eq_zero hP
  have hd1 : 1 ≤ d := Nat.pos_of_ne_zero (fun h0 => by
    rw [h0, Nat.zero_dvd] at hdK; omega)
  have hdle : d ≤ K := Nat.le_of_dvd (by omega) hdK
  have hdeq : d = K := by
    by_contra hne
    have hlt : d < K := lt_of_le_of_ne hdle hne
    obtain ⟨x1, y1, h1', hP1, -, -⟩ := mulOK_of_ne W h d hd1
      (fun k hk1 hk2 => hmin k hk1 (by omega))
    rw [hd, addOrderOf_nsmul_eq_zero] at hP1
    exact Point.some_ne_zero h1' hP1.symm
  exact hdeq ▸ hdn

/-! ## ★出典の紐付け(`.src`) -/

def exists_divisor_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(最小の分点多項式の根は点の位数)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
