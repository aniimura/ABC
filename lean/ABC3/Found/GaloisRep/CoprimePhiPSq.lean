import ABC3.Found.GaloisRep.ZeroSet

/-!
# Galois (G1) 第 61 ブロック —— **★★★★★★`Φ_n` と `ΨSq_n` は共通根を持たない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★零集合が部分群であることから 3 行で出る

共通根 `x₀` があると、第 44 ブロックより

    preΨ_{n+1}(x₀)·preΨ_{n−1}(x₀)·f_n(x₀) = x₀ΨSq_n(x₀) − Φ_n(x₀) = 0

で、第 50 ブロック(`PSq_mul_PSq`)より

    ΨSq_{n+1}(x₀)·ΨSq_{n−1}(x₀) = (preΨ_{n+1}preΨ_{n−1}f_n)²(x₀) = 0

★したがって `ΨSq_{n+1}(x₀) = 0` か `ΨSq_{n−1}(x₀) = 0`。
★★第 60 ブロック(零集合は位数 `K` の倍数全体)より

    K ∣ n  かつ  K ∣ n±1  ⟹  K ∣ 1  ⟹  K = 1  ⟹  ΨSq_1 = 1 = 0  ✗

## ★★これで何が開くか

`Φ_n − c·ΨSq_n` の根はすべて `ΨSq_n(x) ≠ 0` を満たす
——★**乗法公式が全ての根に使える**ようになる。
★★Wronskian(第 55)と合わせて `#E[n] = n²` へ進める。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `Phi_PSq_no_common_root` | ★★★★★★**共通根は無い** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★★**`Φ_n` と `ΨSq_n` は共通根を持たない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`x` が曲線上の点の x 座標であることを仮定する(代数閉体では常に取れる)。 -/
theorem Phi_PSq_no_common_root {x y : F} (hΔ : W.Δ ≠ 0) (h : W.toAffine.Nonsingular x y)
    (n : ℕ) (hn : 1 ≤ n)
    (hPhi : (W.Φ (n : ℤ)).eval x = 0) (hPSq : (W.ΨSq (n : ℤ)).eval x = 0) : False := by
  classical
  rcases Nat.lt_or_ge n 2 with hn1 | hn2
  · have hn' : n = 1 := by omega
    subst hn'
    rw [show ((1 : ℕ) : ℤ) = 1 from rfl, WeierstrassCurve.ΨSq_one] at hPSq
    simp at hPSq
  have hex : ∃ k : ℕ, 1 ≤ k ∧ (W.ΨSq (k : ℤ)).eval x = 0 := ⟨n, hn, hPSq⟩
  obtain ⟨hK1, hKroot⟩ := Nat.find_spec hex
  set K := Nat.find hex with hKdef
  have hmin : ∀ j : ℤ, 1 ≤ j → j < (K : ℤ) → (W.ΨSq j).eval x ≠ 0 := by
    intro j hj1 hj2 hj3
    lift j to ℕ using (by omega : (0 : ℤ) ≤ j) with j'
    exact Nat.find_min hex (by exact_mod_cast hj2) ⟨by exact_mod_cast hj1, hj3⟩
  have hd := x_mul_PSq_sub_Phi W (n : ℤ) x
  rw [hPhi, hPSq, mul_zero, sub_zero] at hd
  have hprod := PSq_mul_PSq W (n : ℤ) x
  rw [← hd] at hprod
  simp only [zero_pow (two_ne_zero)] at hprod
  have hdvdn : (K : ℤ) ∣ (n : ℤ) :=
    dvd_of_PSq_root W hΔ h (K : ℤ) (by exact_mod_cast hK1) (by exact_mod_cast hKroot) hmin n hPSq
  have hone : (K : ℤ) ∣ 1 := by
    rcases mul_eq_zero.1 hprod with hp | hp
    · have hdd : (K : ℤ) ∣ ((n + 1 : ℕ) : ℤ) := by
        refine dvd_of_PSq_root W hΔ h (K : ℤ) (by exact_mod_cast hK1)
          (by exact_mod_cast hKroot) hmin (n + 1) ?_
        rw [show (((n + 1 : ℕ)) : ℤ) = (n : ℤ) + 1 by push_cast; ring]; exact hp
      push_cast at hdd
      simpa using dvd_sub hdd hdvdn
    · have hdd : (K : ℤ) ∣ ((n - 1 : ℕ) : ℤ) := by
        refine dvd_of_PSq_root W hΔ h (K : ℤ) (by exact_mod_cast hK1)
          (by exact_mod_cast hKroot) hmin (n - 1) ?_
        rw [show (((n - 1 : ℕ)) : ℤ) = (n : ℤ) - 1 by
          rw [Nat.cast_sub (by omega)]; norm_num]
        exact hp
      rw [Nat.cast_sub (by omega)] at hdd
      push_cast at hdd
      simpa using dvd_sub hdvdn hdd
  have hKone : K = 1 := by
    have hh := Int.eq_one_of_dvd_one (by exact_mod_cast Nat.zero_le K) hone
    exact_mod_cast hh
  rw [hKone, show ((1 : ℕ) : ℤ) = 1 from rfl, WeierstrassCurve.ΨSq_one] at hKroot
  simp at hKroot

/-! ## ★出典の紐付け(`.src`) -/

def Phi_PSq_no_common_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Φ_n と ΨSq_n の互いに素性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
