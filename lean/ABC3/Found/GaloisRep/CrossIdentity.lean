import ABC3.Found.GaloisRep.FourTorsion

/-!
# Galois (G1) 第 40 ブロック —— **★★★★★交叉恒等式 `Φ₃ΨSq₂ − Φ₂ΨSq₃ = −preΨ₅`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`E[5]` への骨

`5P = O` は `3P = −2P`、すなわち **`x(3P) = x(2P)`** と同値。
乗法公式を入れると `Φ₃·ΨSq₂ = Φ₂·ΨSq₃`、そして本ブロックの恒等式で

    preΨ₅(x) = 0    ⟹    ΨSq₅(x) = preΨ₅(x)² = 0

## ★★★b-不変量の関係は**要らなかった**

§9-379 では計算機で `4b₈ = b₂b₆ − b₄²` を法として確認したが、
★Lean で書いてみると **第 36 の `Φ₂ = X·Ψ₂Sq − Ψ₃` だけ**で出た:

    Φ₃ΨSq₂ − Φ₂ΨSq₃ = Ψ₃²·(X·Ψ₂Sq − Φ₂ − Ψ₃) − preΨ₅ + preΨ₅ = −preΨ₅

★★**計算機の測定は「成り立つこと」を教えるが、「最短の道」は教えない**
——Lean で書くときに改めて構造を見ると短くなることがある。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `preP_five` | ★`preΨ₅ = preΨ₄·Ψ₂Sq² − Ψ₃³` |
| `phi_cross_32` | ★★★★**交叉恒等式** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

theorem preP_five : W.preΨ 5 = W.preΨ₄ * W.Ψ₂Sq ^ 2 - W.Ψ₃ ^ 3 := by
  have h := preNormEDS_odd (W.Ψ₂Sq ^ 2) W.Ψ₃ W.preΨ₄ 2
  norm_num at h
  exact h

/-- ★★★`Φ₃·ΨSq₂ − Φ₂·ΨSq₃ = −preΨ₅`。 -/
theorem phi_cross_32 : W.Φ 3 * W.ΨSq 2 - W.Φ 2 * W.ΨSq 3 = -(W.preΨ 5) := by
  rw [WeierstrassCurve.ΨSq_two, WeierstrassCurve.ΨSq_three, WeierstrassCurve.Φ_three,
    preP_five W, Phi_two_eq W]
  ring


/-! ## ★出典の紐付け(`.src`) -/

def phi_cross_32.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(交叉恒等式——E[5] への骨)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
