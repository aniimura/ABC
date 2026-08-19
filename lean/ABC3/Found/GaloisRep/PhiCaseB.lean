import ABC3.Found.GaloisRep.PhiDegenerate

/-!
# Galois (G1) 第 43 ブロック —— **★★★★★帰納の場合 (b) を一様な形に**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`(n+1)P = O` の場合

`(n+1)P = O` は `nP = −P`、すなわち **`x(nP) = x`** と同値。
乗法公式 `x(nP)·ΨSq_n = Φ_n` を入れると `x·ΨSq_n(x) = Φ_n(x)` になり、
本ブロックの補題で

    ΨSq_{n+1}(x) · ΨSq_{n−1}(x) = 0

が出る。★★どちらが 0 かは分からないが、**積が 0** という形なら偶奇によらず一様である。

## ★★偶奇で機構が違うのに結論は同じ

| `n` | `Φ_n` の余分な因子 | `ΨSq_{n±1}` の因子 |
|---|---|---|
| 偶 | `1` | `1`(`n±1` は奇) |
| 奇 | `Ψ₂Sq` | `Ψ₂Sq`(`n±1` は偶) |

★★★どちらの場合も `preΨ_{n+1}·preΨ_{n−1}·(因子) = 0` から
`(preΨ_{n+1}preΨ_{n−1})²·(因子)² = 0` が出る——**因子が揃う**のが効く。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `PSq_succ_mul_PSq_pred_eq_zero` | ★★★★**場合 (b) の一様形** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- ★★★`x·ΨSq_n(x) = Φ_n(x)` なら `ΨSq_{n+1}(x)·ΨSq_{n−1}(x) = 0`。 -/
theorem PSq_succ_mul_PSq_pred_eq_zero {n : ℤ} {x : F}
    (h : x * (W.ΨSq n).eval x = (W.Φ n).eval x) :
    (W.ΨSq (n + 1)).eval x * (W.ΨSq (n - 1)).eval x = 0 := by
  have hkey : (W.preΨ (n + 1)).eval x * (W.preΨ (n - 1)).eval x
      * (if Even n then (1 : F) else W.Ψ₂Sq.eval x) = 0 := by
    have hd := Phi_def W n
    rw [hd] at h
    by_cases hn : Even n
    · rw [if_pos hn] at h ⊢
      simp only [eval_sub, eval_mul, eval_X, eval_one] at h
      linear_combination h
    · rw [if_neg hn] at h ⊢
      simp only [eval_sub, eval_mul, eval_X] at h
      linear_combination h
  have hp : ∀ m : ℤ, (W.ΨSq m).eval x
      = ((W.preΨ m).eval x) ^ 2 * (if Even m then W.Ψ₂Sq.eval x else 1) := by
    intro m
    rw [WeierstrassCurve.ΨSq]
    by_cases hm : Even m
    · rw [if_pos hm, if_pos hm, eval_mul, eval_pow]
    · rw [if_neg hm, if_neg hm, mul_one, eval_pow, mul_one]
  rw [hp, hp]
  by_cases hn : Even n
  · have h1 : ¬ Even (n + 1) := by simpa [Int.even_add_one] using hn
    have h2 : ¬ Even (n - 1) := by simpa [Int.even_sub] using hn
    rw [if_neg h1, if_neg h2, if_pos hn] at *
    linear_combination ((W.preΨ (n + 1)).eval x * (W.preΨ (n - 1)).eval x) * hkey
  · have h1 : Even (n + 1) := by simpa [Int.even_add_one] using hn
    have h2 : Even (n - 1) := by simpa [Int.even_sub] using hn
    rw [if_pos h1, if_pos h2, if_neg hn] at *
    linear_combination ((W.preΨ (n + 1)).eval x * (W.preΨ (n - 1)).eval x
      * W.Ψ₂Sq.eval x) * hkey


/-! ## ★出典の紐付け(`.src`) -/

def PSq_succ_mul_PSq_pred_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の帰納——(n+1)P = O の場合)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
