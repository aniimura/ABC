import ABC3.Found.GaloisRep.CoordinateUnits

/-!
# Galois (G5) 第 129 ブロック —— **★★★★★★平方完成 `z² = Ψ₂Sq(x)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★層 3 の入口(経路 2)の土台

§9-438 で `IsIntegrallyClosed W.CoordinateRing` への経路を 2 本測った。
★経路 2 は**平方完成**である:

    z := 2y + a₁x + a₃   と置くと   z² = Ψ₂Sq(x)

★★これは mathlib の `C_Ψ₂Sq : C Ψ₂Sq = ψ₂² − 4·(Weierstrass 多項式)` と
`AdjoinRoot.mk_self`(Weierstrass 多項式は座標環で 0)から**直ちに出る**。

★★★すなわち

    F[W] = F[x][z]  で  z² = Ψ₂Sq(x)

という**二次拡大の形**が確定した。`Ψ₂Sq` が squarefree なら整閉である、
というのが経路 2 の残りである。

## ★★足場

| mathlib | 役割 |
|---|---|
| `WeierstrassCurve.C_Ψ₂Sq` | `C Ψ₂Sq = ψ₂² − 4·polynomial` |
| `WeierstrassCurve.ψ₂ = polynomialY` | `ψ₂` は `Y` 偏微分 |
| `AdjoinRoot.mk_self` | Weierstrass 多項式は座標環で 0 |
| 自前(第 116) | `eval₂_genX`——`F[X]` の元を生成点で評価 |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mk_C_Psi2Sq` | ★★★★`mk(C Ψ₂Sq) = mk(ψ₂)²` |
| `genZ` | ★平方完成の変数 `z = 2y + a₁x + a₃` |
| `genZ_sq` | ★★★★★★**`z² = Ψ₂Sq(x)`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★★**`mk(C Ψ₂Sq) = mk(ψ₂)²`**——Weierstrass 多項式が座標環で消えることの帰結。 -/
theorem mk_C_Psi2Sq (W : WeierstrassCurve.Affine F) :
    CoordinateRing.mk W (Polynomial.C W.Ψ₂Sq) = (CoordinateRing.mk W W.ψ₂) ^ 2 := by
  rw [WeierstrassCurve.C_Ψ₂Sq, map_sub, map_pow]
  have h : CoordinateRing.mk W W.toAffine.polynomial = 0 := AdjoinRoot.mk_self
  rw [show (4 : Polynomial (Polynomial F)) * W.toAffine.polynomial
      = 4 * W.toAffine.polynomial from rfl, map_mul, h, mul_zero, sub_zero]

/-- ★平方完成の変数 `z = 2y + a₁x + a₃`。 -/
noncomputable def genZ (W : WeierstrassCurve.Affine F) : W.CoordinateRing :=
  2 * genY W + algebraMap F W.CoordinateRing W.a₁ * genX W
    + algebraMap F W.CoordinateRing W.a₃

theorem mk_psi2 (W : WeierstrassCurve.Affine F) :
    CoordinateRing.mk W W.ψ₂ = genZ W := by
  have halg : ∀ a : F, CoordinateRing.mk W (Polynomial.C (Polynomial.C a))
      = algebraMap F W.CoordinateRing a := fun _ => rfl
  rw [WeierstrassCurve.ψ₂, WeierstrassCurve.Affine.polynomialY, genZ, genX, genY]
  simp only [map_add, map_mul, map_ofNat, halg]
  ring

/-- ★★★★★★**平方完成**——`(2y + a₁x + a₃)² = Ψ₂Sq(x)` in `F[W]`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これにより `F[W] = F[x][z]`(`z² = Ψ₂Sq(x)`)という**二次拡大の形**が確定する。
★★`Ψ₂Sq` が squarefree なら整閉である、というのが `IsIntegrallyClosed` への経路 2 である。 -/
theorem genZ_sq (W : WeierstrassCurve.Affine F) :
    (genZ W) ^ 2 = Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W) W.Ψ₂Sq := by
  rw [eval₂_genX, mk_C_Psi2Sq, mk_psi2]

/-! ## ★出典の紐付け(`.src`) -/

def genZ_sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——座標環の平方完成 z² = Ψ₂Sq(x))",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
