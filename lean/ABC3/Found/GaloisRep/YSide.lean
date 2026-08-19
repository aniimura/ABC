import ABC3.Found.GaloisRep.PhiCross

/-!
# Galois (G1) 第 49 ブロック —— **★★★★★★y 側の恒等式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★§9-390 で開いた穴が塞がる

mathlib の `addY_sub_negY_addY`

    ψ(P₁+P₂)·(x₂−x₁) = ψ(P₂)(x₁−x₃) − ψ(P₁)(x₂−x₃)

に `P₁ = nP`、`P₂ = P` を入れ、`ΨSq_n²ΨSq_{n+1}²` を掛けると、
★必要な多項式恒等式は

    (Φ_n ΨSq_{n+1} − Φ_{n+1} ΨSq_n)·ΨSq_n ΨSq_{n+1}
      − preΨ_{2n}·(preΨ_{n+2} preΨ_n f_{n+1})·ΨSq_{n+1}
    = preΨ_{2n+2}·(preΨ_{n+1} preΨ_{n−1} f_n)·ΨSq_n          …(Y)

だけになる(`f_k = 1` (k 偶)、`Ψ₂Sq` (k 奇))。

## ★★★★★(Y) は (★) と `preΨ_odd` でちょうど閉じる

§9-394 で Python(sympy)で測った証明書:

| 生成元 | 係数 |
|---|---|
| `preΨ_even n` | ★`−preΨ_n preΨ_{n+1}² preΨ_{n+2} Ψ₂Sq` |
| `preΨ_even (n+1)` | ★`−preΨ_n² preΨ_{n−1} preΨ_{n+1} Ψ₂Sq` |
| **(★)**`preP_star` | ★★`−preΨ_n² preΨ_{n+1}² Ψ₂Sq` |
| `preΨ_odd n` | ★★`−preΨ_n² preΨ_{n+1}² Ψ₂Sq` |

★★★**パリティに依らず同じ係数**であった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `preP_star` | ★★★曲線版の (★)(第 47 ブロックの言い換え) |
| `y_side` | ★★★★★★**(Y)——y 側の帰納段の心臓** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

variable {R : Type} [CommRing R] (W : WeierstrassCurve R)

/-- ★★★**曲線版の (★)**。

`W.preΨ = preNormEDS (Ψ₂Sq²) Ψ₃ preΨ₄` は定義そのものなので、
第 47 ブロックの `eds_star` を代入するだけである。 -/
theorem preP_star (n : ℤ) :
    W.preΨ (n + 3) * W.preΨ (n - 1) * W.preΨ n ^ 2
      - W.preΨ (n + 2) * W.preΨ (n - 2) * W.preΨ (n + 1) ^ 2 = W.preΨ (2 * n + 1) :=
  eds_star (W.Ψ₂Sq ^ 2) W.Ψ₃ W.preΨ₄ n

/-- ★★★★★★**y 側の恒等式 (Y)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これに `addY_sub_negY_addY` を合わせると、乗法公式の y 側が帰納で回る。
★★§9-390 で「多項式恒等式ではない」と判定した穴は、
**(★) を経由すれば恒等式になる**——それが本ブロックである。 -/
theorem y_side (n : ℤ) :
    (W.Φ n * W.ΨSq (n + 1) - W.Φ (n + 1) * W.ΨSq n) * (W.ΨSq n * W.ΨSq (n + 1))
        - W.preΨ (2 * n) * (W.preΨ (n + 2) * W.preΨ n * (if Even (n + 1) then 1 else W.Ψ₂Sq))
          * W.ΨSq (n + 1)
      = W.preΨ (2 * (n + 1))
          * (W.preΨ (n + 1) * W.preΨ (n - 1) * (if Even n then 1 else W.Ψ₂Sq)) * W.ΨSq n := by
  have he0 := W.preΨ_even n
  have he1 := W.preΨ_even (n + 1)
  have hst := preP_star W n
  have hod := W.preΨ_odd n
  rw [show n + 1 - 1 = n by ring, show n + 1 + 2 = n + 3 by ring,
    show n + 1 - 2 = n - 1 by ring, show n + 1 + 1 = n + 2 by ring] at he1
  simp only [WeierstrassCurve.Φ, WeierstrassCurve.ΨSq]
  rw [show n + 1 + 1 = n + 2 by ring, show n + 1 - 1 = n by ring]
  by_cases hn : Even n
  · have h1 : ¬ Even (n + 1) := by rw [Int.even_add_one]; exact not_not.mpr hn
    simp only [if_pos hn, if_neg h1] at hod ⊢
    linear_combination (-(W.preΨ n) * W.preΨ (n + 1) ^ 2 * W.preΨ (n + 2) * W.Ψ₂Sq) * he0
      + (-(W.preΨ n) ^ 2 * W.preΨ (n - 1) * W.preΨ (n + 1) * W.Ψ₂Sq) * he1
      + (-(W.preΨ n) ^ 2 * W.preΨ (n + 1) ^ 2 * W.Ψ₂Sq) * hst
      + (-(W.preΨ n) ^ 2 * W.preΨ (n + 1) ^ 2 * W.Ψ₂Sq) * hod
  · have h1 : Even (n + 1) := Int.even_add_one.mpr hn
    simp only [if_neg hn, if_pos h1] at hod ⊢
    linear_combination (-(W.preΨ n) * W.preΨ (n + 1) ^ 2 * W.preΨ (n + 2) * W.Ψ₂Sq) * he0
      + (-(W.preΨ n) ^ 2 * W.preΨ (n - 1) * W.preΨ (n + 1) * W.Ψ₂Sq) * he1
      + (-(W.preΨ n) ^ 2 * W.preΨ (n + 1) ^ 2 * W.Ψ₂Sq) * hst
      + (-(W.preΨ n) ^ 2 * W.preΨ (n + 1) ^ 2 * W.Ψ₂Sq) * hod

/-! ## ★出典の紐付け(`.src`) -/

def y_side.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の y 側の恒等式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
