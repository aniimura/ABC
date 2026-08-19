import ABC3.Found.GaloisRep.TorsionIdeal

/-!
# Galois (G5) 第 114 ブロック —— **★★★★★★生成点は曲線の点である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★これが (G5) の 2 枚の葉の土台である

`Skeleton/GaloisRep/WeilFunctionField.lean` の葉は
**平行移動 `τ_Q`** と **乗法 `[n]^*`** の関数体への引き戻しであった。
★どちらも「座標関数を別の式に送る環準同型」を作ることであり、
そのためには **`W.polynomial` が消えること**を言わねばならない。

★★ここで効くのが次の事実である:

    座標環 `F[W] = AdjoinRoot(W.polynomial)` の生成元 `(ξ, η)` は、
    **`F[W]` 上に底変換した曲線の点である**

★★★すなわち `AdjoinRoot.mk_self` がそのまま Weierstrass 方程式になる。
★★★★さらに `W` が楕円曲線なら `equation_iff_nonsingular` で**非特異点**になり、
mathlib の `Point.some` が使えて、**群法則がそのまま生成点に効く**。

## ★★★★★これで平行移動が「点の加法」になる

`τ_Q` の像は「生成点 `+` 定数点 `Q`」の座標である。★mathlib の `nonsingular_add` が
その座標が再び方程式を満たすことを保証するので、
**引き戻しが環準同型であることの本体が mathlib から出る**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `equation_gen` | ★★★座標環の生成点が方程式を満たす |
| `equation_coord` | ★★★★関数体でも同じ |
| `nonsingular_coord` | ★★★★★★**生成点は非特異点**(楕円曲線なら) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★座標環の生成元 -/

/-- ★座標環の生成元 `x`。 -/
noncomputable def genX (W : WeierstrassCurve.Affine F) : W.CoordinateRing :=
  CoordinateRing.mk W (Polynomial.C Polynomial.X)

/-- ★座標環の生成元 `y`。 -/
noncomputable def genY (W : WeierstrassCurve.Affine F) : W.CoordinateRing :=
  CoordinateRing.mk W Polynomial.X

/-- ★★★**生成点は Weierstrass 方程式を満たす**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★座標環は `AdjoinRoot W.polynomial` なので、`AdjoinRoot.mk_self` がそのまま方程式になる。 -/
theorem equation_gen (W : WeierstrassCurve.Affine F) :
    (W.map (algebraMap F W.CoordinateRing)).Equation (genX W) (genY W) := by
  rw [WeierstrassCurve.Affine.equation_iff, genX, genY]
  have h : CoordinateRing.mk W (Polynomial.X ^ 2
      + Polynomial.C (Polynomial.C W.a₁ * Polynomial.X + Polynomial.C W.a₃) * Polynomial.X
      - Polynomial.C (Polynomial.X ^ 3 + Polynomial.C W.a₂ * Polynomial.X ^ 2
          + Polynomial.C W.a₄ * Polynomial.X + Polynomial.C W.a₆)) = 0 := by
    rw [← WeierstrassCurve.Affine.polynomial]
    exact AdjoinRoot.mk_self
  simp only [map_add, map_sub, map_mul, map_pow] at h
  have halg : ∀ a : F, algebraMap F W.CoordinateRing a = CoordinateRing.mk W (C (C a)) := by
    intro a
    rfl
  simp only [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₃,
    WeierstrassCurve.map_a₄, WeierstrassCurve.map_a₆, halg]
  linear_combination h

/-! ## ★関数体へ -/

/-- ★関数体の座標関数 `x`。 -/
noncomputable def coordX (W : WeierstrassCurve.Affine F) : W.FunctionField :=
  algebraMap W.CoordinateRing W.FunctionField (genX W)

/-- ★関数体の座標関数 `y`。 -/
noncomputable def coordY (W : WeierstrassCurve.Affine F) : W.FunctionField :=
  algebraMap W.CoordinateRing W.FunctionField (genY W)

/-- ★底変換の合成——`F → F[W] → F(W)`。 -/
theorem map_coord_functionField (W : WeierstrassCurve.Affine F) :
    (W.map (algebraMap F W.CoordinateRing)).map (algebraMap W.CoordinateRing W.FunctionField)
      = W.map (algebraMap F W.FunctionField) :=
  WeierstrassCurve.ext rfl rfl rfl rfl rfl

/-- ★★★★**関数体の生成点も Weierstrass 方程式を満たす**。 -/
theorem equation_coord (W : WeierstrassCurve.Affine F) :
    (W.map (algebraMap F W.FunctionField)).Equation (coordX W) (coordY W) := by
  rw [← map_coord_functionField W]
  exact (equation_gen W).map (algebraMap W.CoordinateRing W.FunctionField)

/-- ★★★★★★**生成点は非特異点である**(楕円曲線なら)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これにより `Point.some _ _ (nonsingular_coord W)` が作れ、
**mathlib の群法則が生成点にそのまま効く**。 -/
theorem nonsingular_coord (W : WeierstrassCurve.Affine F) [W.IsElliptic] :
    (W.map (algebraMap F W.FunctionField)).Nonsingular (coordX W) (coordY W) :=
  (WeierstrassCurve.Affine.equation_iff_nonsingular).1 (equation_coord W)

/-- ★★★★★**生成点**——関数体上の曲線の点として。 -/
noncomputable def genericPoint (W : WeierstrassCurve.Affine F) [W.IsElliptic] :
    (W.map (algebraMap F W.FunctionField)).Point :=
  Point.some (coordX W) (coordY W) (nonsingular_coord W)

/-! ## ★出典の紐付け(`.src`) -/

def nonsingular_coord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——生成点が曲線の非特異点であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
