import ABC3.Found.GaloisRep.SquareCompletion

/-!
# Galois (G5) 第 130 ブロック —— **★★★★`Ψ₂Sq` の判別式は `16Δ`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★mathlib に既にあった

第 129 ブロックで `F[W] = F[x][z]`(`z² = Ψ₂Sq(x)`)を得た。
★経路 2 の次は「`Ψ₂Sq` が squarefree」であり、そのために判別式が要る。

★★**2026-08-20 の実測で、mathlib が既に持っていることが分かった**:

| mathlib | 内容 |
|---|---|
| `WeierstrassCurve.twoTorsionPolynomial` | 3 次式 `⟨4, b₂, 2b₄, b₆⟩` |
| `twoTorsionPolynomial_discr` | `discr = 16·Δ` |
| `twoTorsionPolynomial_discr_ne_zero` | `IsUnit 2` と `IsUnit Δ` から `discr ≠ 0` |

★★★**`Ψ₂Sq` は `twoTorsionPolynomial.toPoly` そのものである**——
判別式の計算を自分で書く必要は無かった。

★★★★一度は自分で `discr = 16Δ` を `ring` で証明したが、
**mathlib の重複だった**——記録に残す。

## ★★残っている段

`discr ≠ 0` から `Squarefree Ψ₂Sq` を出す段が残る。
★mathlib の `Cubic.discr_ne_zero_iff_roots_nodup` で根が相異なることは出るが、
★★**`roots.Nodup ⟹ Squarefree`(あるいは `Separable`)の一撃は無い**
(2026-08-20 実測、`exact?` で 0 件)。
`Polynomial.Separable.squarefree` はあるので、`Separable` を経由することになる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `Psi2Sq_eq_twoTorsionPolynomial` | ★★★`Ψ₂Sq = twoTorsionPolynomial.toPoly` |
| `twoTorsionPolynomial_discr_ne_zero_of_isElliptic` | ★★★★`discr ≠ 0` |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

variable {F : Type} [Field F]

/-- ★★★**`Ψ₂Sq` は mathlib の `twoTorsionPolynomial` そのものである**。 -/
theorem Psi2Sq_eq_twoTorsionPolynomial (W : WeierstrassCurve F) :
    W.Ψ₂Sq = W.twoTorsionPolynomial.toPoly := by
  rw [WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.twoTorsionPolynomial, Cubic.toPoly]

/-- ★★★★**楕円曲線(標数 ≠ 2)では `Ψ₂Sq` の判別式は 0 でない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `twoTorsionPolynomial_discr = 16·Δ` と `IsElliptic.isUnit` から出る。 -/
theorem twoTorsionPolynomial_discr_ne_zero_of_isElliptic (W : WeierstrassCurve F) [W.IsElliptic]
    (h2 : IsUnit (2 : F)) : W.twoTorsionPolynomial.discr ≠ 0 :=
  W.twoTorsionPolynomial_discr_ne_zero h2 WeierstrassCurve.IsElliptic.isUnit

/-! ## ★出典の紐付け(`.src`) -/

def Psi2Sq_eq_twoTorsionPolynomial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——Ψ₂Sq の判別式が 16Δ であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
