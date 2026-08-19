import ABC3.Found.GaloisRep.EdsSomos

/-!
# Galois (G1) 第 48 ブロック —— **★★★★★交叉恒等式の一般形**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★第 40 ブロックの一般化——しかも `preΨ_even` だけで出る

第 40 ブロックで `Φ₃ΨSq₂ − Φ₂ΨSq₃ = −preΨ₅` を出したが、
★これは一般の `n` で成り立つ:

    Φ_{n+1}·ΨSq_{n−1} = Φ_{n−1}·ΨSq_{n+1} − preΨ_{2n}·Ψ₂Sq

## ★★★★★★これは `x(mP) − x(nP)` の公式そのもの

`x((n+1)P) = Φ_{n+1}/ΨSq_{n+1}` を認めれば、この恒等式は

    x((n+1)P) − x((n−1)P) = −preΨ_{2n}·Ψ₂Sq / (ΨSq_{n+1}ΨSq_{n−1})

と読める。★★mathlib の `addX_eq_addX_negY_sub`

    x(P₁+P₂) = x(P₁−P₂) − ψ(P₁)ψ(P₂)/(x(P₂)−x(P₁))²

とちょうど対応する——**帰納段で `x_{n+1}` を決めるのがこれ**である。

## ★★★★驚くべきことに `preΨ_even` 1 本で出る

証明は §9-394 で測ったとおり:

| 段 | 要るもの |
|---|---|
| `ΨSq_{n±1} = preΨ_{n±1}²·f_n` | ★`ΨSq` の定義(パリティが揃う) |
| 残り | ★★**`preΨ_even` ただ 1 本**、係数は `Ψ₂Sq` |

★★★3 項関係も Somos も**要らない**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `Phi_succ_mul_PSq_pred` | ★★★★★**交叉恒等式の一般形** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

variable {R : Type} [CommRing R] (W : WeierstrassCurve R)

/-- ★★★★★**交叉恒等式の一般形**——`Φ_{n+1}ΨSq_{n−1} = Φ_{n−1}ΨSq_{n+1} − preΨ_{2n}Ψ₂Sq`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 40 ブロック `phi_cross_32` は `n = 2` の場合である。
★★証明は `preΨ_even` ただ 1 本、係数は `Ψ₂Sq`(§9-394 で実測)。 -/
theorem Phi_succ_mul_PSq_pred (n : ℤ) :
    W.Φ (n + 1) * W.ΨSq (n - 1)
      = W.Φ (n - 1) * W.ΨSq (n + 1) - W.preΨ (2 * n) * W.Ψ₂Sq := by
  have hrec := W.preΨ_even n
  simp only [WeierstrassCurve.Φ, WeierstrassCurve.ΨSq]
  rw [show n + 1 + 1 = n + 2 by ring, show n + 1 - 1 = n by ring,
    show n - 1 + 1 = n by ring, show n - 1 - 1 = n - 2 by ring]
  by_cases hn : Even n
  · have h1 : ¬ Even (n + 1) := by rw [Int.even_add_one]; exact not_not.mpr hn
    have h2 : ¬ Even (n - 1) := by rw [Int.even_sub_one]; exact not_not.mpr hn
    simp only [if_neg h1, if_neg h2]
    linear_combination W.Ψ₂Sq * hrec
  · have h1 : Even (n + 1) := Int.even_add_one.mpr hn
    have h2 : Even (n - 1) := Int.even_sub_one.mpr hn
    simp only [if_pos h1, if_pos h2]
    linear_combination W.Ψ₂Sq * hrec

/-! ## ★出典の紐付け(`.src`) -/

def Phi_succ_mul_PSq_pred.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(交叉恒等式の一般形——x 側の帰納段)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
