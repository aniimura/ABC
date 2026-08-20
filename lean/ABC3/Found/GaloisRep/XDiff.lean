import ABC3.Found.GaloisRep.PhiCaseB

/-!
# Galois (G1) 第 44 ブロック —— **`x` と `x(nP)` の差**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★一般の場合に要る入力

加法公式(mathlib `addX_eq_addX_negY_sub`)は

    x((n+1)P) = x((n−1)P) − ψ(nP)·ψ(P) / (x − x(nP))²

★分母 `x − x(nP)` を分点多項式で書くのが本ブロックである:

    x·ΨSq_n(x) − Φ_n(x) = preΨ_{n+1}(x)·preΨ_{n−1}(x)·(if Even n then 1 else Ψ₂Sq(x))

★★乗法公式 `x(nP)·ΨSq_n = Φ_n` を入れると、左辺は `(x − x(nP))·ΨSq_n(x)` である。

## ★★★これが「差が消えるのはいつか」を与える

★`x = x(nP)` ⟺ `preΨ_{n+1}(x)preΨ_{n−1}(x)·(因子) = 0`
——第 43 の場合 (b) はこれの言い換えであった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `x_mul_PSq_sub_Phi` | ★★★**`x` と `x(nP)` の差(分母払い)** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- ★★★`x` と `x(nP)` の差(分母を払った形)。 -/
theorem x_mul_PSq_sub_Phi (n : ℤ) (x : F) :
    x * (W.ΨSq n).eval x - (W.Φ n).eval x
      = (W.preΨ (n + 1)).eval x * (W.preΨ (n - 1)).eval x
        * (if Even n then (1 : F) else W.Ψ₂Sq.eval x) := by
  rw [Phi_def W n]
  by_cases hn : Even n
  · rw [if_pos hn, if_pos hn]
    simp only [eval_sub, eval_mul, eval_X, eval_one]
    ring
  · rw [if_neg hn, if_neg hn]
    simp only [eval_sub, eval_mul, eval_X]
    ring


/-! ## ★出典の紐付け(`.src`) -/

def x_mul_PSq_sub_Phi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の帰納——x 座標の差)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
