import ABC3.Found.GaloisRep.PhiCases

/-!
# Galois (G1) 第 42 ブロック —— **★★★★退化した点での `Φ_{n+1}`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★帰納の場合 (a) を点の言葉で

第 41 は `preΨ_n = 0`(多項式として)の場合だったが、
帰納で実際に手に入るのは **`ΨSq_n(x) = 0`(点での値)** である。

    ΨSq_n(x) = 0  ⟹  Φ_{n+1}(x) = x · ΨSq_{n+1}(x)

★★これは `nP = O` のとき `(n+1)P = P` になることの、多項式側の対応物である。

## ★★偶奇で 2 通り——どちらも通る

| `n` | `ΨSq_n` | `Φ_{n+1}` の因子 |
|---|---|---|
| 偶 | `preΨ_n²·Ψ₂Sq` | ★`Ψ₂Sq`(`Even (n+1)` が偽) |
| 奇 | `preΨ_n²` | ★`1`(`Even (n+1)` が真) |

★★★偶の場合は `preΨ_n(x) = 0` **または** `Ψ₂Sq(x) = 0` のどちらでもよい
——どちらでも `Φ_{n+1}` の余分な項が消える。
★奇の場合は `preΨ_n(x) = 0` が出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `Phi_succ_eval_of_PSq_eval_eq_zero` | ★★★★**点での退化の処理** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

theorem Phi_succ_eval_of_PSq_eval_eq_zero {n : ℤ} {x : F} (h : (W.ΨSq n).eval x = 0) :
    (W.Φ (n + 1)).eval x = x * (W.ΨSq (n + 1)).eval x := by
  by_cases hn : Even n
  · have hodd : ¬ Even (n + 1) := by simpa [Int.even_add_one] using hn
    have hPSq : (W.ΨSq n).eval x = ((W.preΨ n).eval x) ^ 2 * (W.Ψ₂Sq.eval x) := by
      rw [WeierstrassCurve.ΨSq, if_pos hn, eval_mul, eval_pow]
    rw [hPSq] at h
    rw [Phi_def W (n + 1), show n + 1 - 1 = n by ring, if_neg hodd]
    simp only [eval_sub, eval_mul, eval_X]
    rcases mul_eq_zero.1 h with h' | h'
    · rw [pow_eq_zero_iff two_ne_zero] at h'
      rw [h']; ring
    · rw [h']; ring
  · have heven : Even (n + 1) := by simpa [Int.even_add_one] using hn
    have hPSq : (W.ΨSq n).eval x = ((W.preΨ n).eval x) ^ 2 := by
      rw [WeierstrassCurve.ΨSq, if_neg hn, mul_one, eval_pow]
    rw [hPSq, pow_eq_zero_iff two_ne_zero] at h
    rw [Phi_def W (n + 1), show n + 1 - 1 = n by ring, if_pos heven]
    simp only [eval_sub, eval_mul, eval_X, eval_one]
    rw [h]; ring


/-! ## ★出典の紐付け(`.src`) -/

def Phi_succ_eval_of_PSq_eval_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の帰納——点での退化の処理)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
