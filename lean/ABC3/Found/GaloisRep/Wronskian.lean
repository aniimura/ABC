import ABC3.Found.GaloisRep.TorsionAll

/-!
# Galois (G1) 第 55 ブロック —— **★★★★★Wronskian は 0 でない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`#E[n] = n²` へ向けての第一歩

`x(nP) = Φ_n(x)/ΨSq_n(x)` が**分離的**であることは

    Wr_n := Φ_n' · ΨSq_n − Φ_n · ΨSq_n'  ≠ 0

と同値である。★不変微分の議論では `Wr_n = n·preΨ_{2n}` だが、
★★**それを証明する必要はない**——**次数と主係数だけ**で `≠ 0` が出る:

| 多項式 | 次数 | 主係数 |
|---|---|---|
| `Φ_n` | `n²` | `1` |
| `ΨSq_n` | `n² − 1` | `n²` |
| `Φ_n'` | `n² − 1` | `n²` |
| `ΨSq_n'` | `n² − 2` | `n²(n² − 1)` |

★★★`2n²−2` 次の係数は `n²·n² − 1·n²(n²−1) = n²` ✅
——**`(n : F) ≠ 0` なら `Wr_n ≠ 0`**。

## ★★不変微分を回避できた理由

★`Wr_n = n·preΨ_{2n}` は帰納で示すしかない(微分の漸化式が要る)が、
★★**主係数だけなら mathlib の `coeff_Φ` / `coeff_ΨSq` の在庫で足りる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `coeff_mul_of_le` | ★次数の上界だけで積の係数を取る |
| `wronskian_coeff` | ★★★`Wr_n` の `2n²−2` 次係数は `n²` |
| `wronskian_ne_zero` | ★★★★★**`Wr_n ≠ 0`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

/-- ★**次数の上界だけで積の係数が取れる**。 -/
theorem coeff_mul_of_le {R : Type u} [CommSemiring R] {f g : R[X]} {a b : ℕ}
    (hf : f.natDegree ≤ a) (hg : g.natDegree ≤ b) :
    (f * g).coeff (a + b) = f.coeff a * g.coeff b := by
  rw [Polynomial.coeff_mul]
  refine Finset.sum_eq_single (a, b) ?_ ?_
  · rintro ⟨i, j⟩ hij hne
    simp only [Finset.mem_antidiagonal] at hij
    rcases lt_trichotomy i a with hlt | heq | hgt
    · have hj : b < j := by omega
      rw [Polynomial.coeff_eq_zero_of_natDegree_lt (lt_of_le_of_lt hg hj), mul_zero]
    · exact absurd (Prod.ext heq (by omega)) hne
    · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (lt_of_le_of_lt hf hgt), zero_mul]
  · intro hmem
    simp at hmem

/-- ★★★**`Wr_n` の `2n²−2` 次の係数は `n²`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Φ_n` が monic で次数 `n²`、`ΨSq_n` が次数 `n²−1` で主係数 `n²` であることから出る。 -/
theorem wronskian_coeff {F : Type u} [Field F] (W : WeierstrassCurve F) (n : ℤ)
    (hn : 2 ≤ n.natAbs) :
    (derivative (W.Φ n) * W.ΨSq n - W.Φ n * derivative (W.ΨSq n)).coeff
      (2 * n.natAbs ^ 2 - 2) = (n : F) ^ 2 := by
  have hN2 : 4 ≤ n.natAbs ^ 2 := by nlinarith [hn]
  set N := n.natAbs ^ 2 with hNdef
  have hdPhi : (W.Φ n).natDegree ≤ N := W.natDegree_Φ_le n
  have hdPsi : (W.ΨSq n).natDegree ≤ N - 1 := W.natDegree_ΨSq_le n
  have hdPhi' : (derivative (W.Φ n)).natDegree ≤ N - 1 :=
    le_trans (natDegree_derivative_le _) (Nat.sub_le_sub_right hdPhi 1)
  have hdPsi' : (derivative (W.ΨSq n)).natDegree ≤ N - 2 :=
    le_trans (natDegree_derivative_le _) (by omega)
  have e1 : (N - 1) + (N - 1) = 2 * N - 2 := by omega
  have e2 : N + (N - 2) = 2 * N - 2 := by omega
  have h1 : (derivative (W.Φ n) * W.ΨSq n).coeff (2 * N - 2)
      = (derivative (W.Φ n)).coeff (N - 1) * (W.ΨSq n).coeff (N - 1) := by
    rw [← e1]; exact coeff_mul_of_le hdPhi' hdPsi
  have h2 : (W.Φ n * derivative (W.ΨSq n)).coeff (2 * N - 2)
      = (W.Φ n).coeff N * (derivative (W.ΨSq n)).coeff (N - 2) := by
    rw [← e2]; exact coeff_mul_of_le hdPhi hdPsi'
  have hcast : ((N : ℕ) : F) = (n : F) ^ 2 := by
    have hZ : ((n.natAbs ^ 2 : ℕ) : ℤ) = n ^ 2 := by
      push_cast [Int.natCast_natAbs]; rw [sq_abs]
    rw [hNdef, ← Int.cast_natCast (R := F), hZ]; push_cast; ring
  rw [Polynomial.coeff_sub, h1, h2, Polynomial.coeff_derivative, Polynomial.coeff_derivative,
    show N - 1 + 1 = N by omega, show N - 2 + 1 = N - 1 by omega,
    W.coeff_Φ n, W.coeff_ΨSq n]
  have hc1 : (((N - 1 : ℕ) : F) + 1) = ((N : ℕ) : F) := by
    rw [Nat.cast_sub (by omega)]; push_cast; ring
  have hc2 : (((N - 2 : ℕ) : F) + 1) = ((N : ℕ) : F) - 1 := by
    rw [Nat.cast_sub (by omega)]; push_cast; ring
  rw [hc1, hc2, hcast]
  ring

/-- ★★★★★**Wronskian は 0 でない**——`(n : F) ≠ 0` のとき。

★これが `Φ_n − c·ΨSq_n` の分離性(重根が有限個の `c` にしか現れないこと)を与える。 -/
theorem wronskian_ne_zero {F : Type u} [Field F] (W : WeierstrassCurve F) (n : ℤ)
    (hn : 2 ≤ n.natAbs) (hchar : (n : F) ≠ 0) :
    derivative (W.Φ n) * W.ΨSq n - W.Φ n * derivative (W.ΨSq n) ≠ 0 := by
  intro h
  have hw := wronskian_coeff W n hn
  rw [h, Polynomial.coeff_zero] at hw
  exact hchar (pow_eq_zero_iff two_ne_zero |>.1 hw.symm)

/-! ## ★出典の紐付け(`.src`) -/

def wronskian_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法写像の分離性——Wronskian が 0 でないこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
