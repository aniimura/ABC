import ABC3.Found.GaloisRep.CrossIdentity

/-!
# Galois (G1) 第 41 ブロック —— **乗法公式の帰納の場合分け**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★乗法公式の帰納で「途中で `O` になる」場合を捌く

`x(mP)·ΨSq_m = Φ_m` の帰納では、**`mP = O` の場合**を別に扱う必要がある。
★そのとき `ΨSq_m(x) = 0`、すなわち `preΨ_m(x) = 0` であり、

    Φ_{m±1} = X · ΨSq_{m±1}

になる——`Φ` の定義に `preΨ_m` が因子として現れるからである。

## ★★mathlib の `Φ` の定義がそのまま使える

    Φ n = X · ΨSq n − preΨ(n+1) · preΨ(n−1) · (if Even n then 1 else Ψ₂Sq)

★**`rfl` で取り出せる**(`Phi_def`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `Phi_def` | ★定義の取り出し(`rfl`) |
| `Phi_succ_of_preP_eq_zero` | ★★`preΨ_n = 0 ⟹ Φ_{n+1} = X·ΨSq_{n+1}` |
| `Phi_pred_of_preP_eq_zero` | ★★`preΨ_n = 0 ⟹ Φ_{n−1} = X·ΨSq_{n−1}` |
| `PSq_eq_zero_of_preP_eq_zero` | ★`preΨ_n = 0 ⟹ ΨSq_n = 0` |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- ★`Φ` の定義そのもの。 -/
theorem Phi_def (n : ℤ) :
    W.Φ n = X * W.ΨSq n - W.preΨ (n + 1) * W.preΨ (n - 1) * (if Even n then 1 else W.Ψ₂Sq) :=
  rfl

/-- ★★`preΨ n = 0` なら `Φ (n+1) = X · ΨSq (n+1)`。 -/
theorem Phi_succ_of_preP_eq_zero {n : ℤ} (h : W.preΨ n = 0) :
    W.Φ (n + 1) = X * W.ΨSq (n + 1) := by
  rw [Phi_def, show n + 1 - 1 = n by ring, h]
  ring

/-- ★★`preΨ n = 0` なら `Φ (n−1) = X · ΨSq (n−1)`。 -/
theorem Phi_pred_of_preP_eq_zero {n : ℤ} (h : W.preΨ n = 0) :
    W.Φ (n - 1) = X * W.ΨSq (n - 1) := by
  rw [Phi_def, show n - 1 + 1 = n by ring, h]
  ring

/-- ★`preΨ n = 0` なら `ΨSq n = 0`。 -/
theorem PSq_eq_zero_of_preP_eq_zero {n : ℤ} (h : W.preΨ n = 0) : W.ΨSq n = 0 := by
  rw [WeierstrassCurve.ΨSq, h]
  ring


/-! ## ★出典の紐付け(`.src`) -/

def Phi_succ_of_preP_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の帰納——途中で O になる場合)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
