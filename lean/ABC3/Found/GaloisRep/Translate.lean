import ABC3.Found.GaloisRep.GenericPoint

/-!
# Galois (G5) 第 115 ブロック —— **★★★★★★平行移動の環準同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★葉 1(平行移動)の本体が取れた

第 114 ブロックで**生成点が曲線の非特異点である**ことが分かった。★すると

    平行移動 τ_Q の像 = 生成点 + 定数点 Q

であり、**mathlib の `nonsingular_add` がその座標が再び方程式を満たすことを保証する**。
★★これがそのまま `AdjoinRoot.lift` の仮説になるので、環準同型

    τ_Q : F[W] →+* F(W)

が構成できる(`translateHom`)。★★★生成元の像は加法公式そのものである。

## ★★★★足場の記録

| mathlib | 役割 |
|---|---|
| `nonsingular_add` | 和の座標が非特異点であること |
| `slope_of_X_ne` | `x` 座標が違うときの傾き |
| `CoordinateRing.smul_basis_eq_zero` | 生成点の `x` が定数でないこと |
| `AdjoinRoot.lift` | 座標環からの環準同型 |
| `Polynomial.eval₂_eval₂RingHom_apply` | `eval₂` と `evalEval` の橋 |
| `Affine.map_polynomial` | 底変換した曲線の多項式 |

★どれも 2026-08-20 に実測して見つけたものである。

## ★★残っている段

`τ_Q` を **`F(W)` の自己同型**にするには、`translateHom` が単射であることが要る
(分数体への延長の条件)。★★`F[W]` は `F[X]` 上階数 2 の自由加群なので整拡大であり、
「整拡大で (0) の上にある素イデアルは (0)」から出るはずだが、
**mathlib での経路は未特定**である——`Skeleton/GaloisRep/WeilFunctionField.lean` に記録する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `genX_ne_const` / `coordX_ne_const` | ★★★生成点の `x` は定数でない |
| `slopeFF_eq` | ★★傾きが mathlib の `slope` と一致 |
| `nonsingular_translate` | ★★★★★★**平行移動した座標も非特異点** |
| `translateHom` | ★★★★★★**平行移動の環準同型** |
| `translateHom_genX` / `_genY` | ★★★生成元の像は加法公式 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★★生成点は定数ではない -/

/-- ★★★**生成点の `x` 座標は定数ではない**。

★`F[W]` が `F[X]` 上階数 2 の自由加群であること(`smul_basis_eq_zero`)から出る。 -/
theorem genX_ne_const (W : WeierstrassCurve.Affine F) (c : F) :
    genX W ≠ algebraMap F W.CoordinateRing c := by
  intro hc
  have h0 : (Polynomial.X - Polynomial.C c) • (1 : W.CoordinateRing)
      + (0 : Polynomial F) • CoordinateRing.mk W Polynomial.X = 0 := by
    rw [zero_smul, add_zero, CoordinateRing.smul, mul_one]
    simp only [map_sub]
    exact sub_eq_zero_of_eq hc
  have h1 := CoordinateRing.smul_basis_eq_zero h0
  have h2 : Polynomial.X - Polynomial.C c ≠ (0 : Polynomial F) := by
    intro h
    have hd := congrArg Polynomial.natDegree h
    simp at hd
  exact h2 h1.1

theorem coordX_ne_const (W : WeierstrassCurve.Affine F) (c : F) :
    coordX W ≠ algebraMap F W.FunctionField c := by
  intro hc
  refine genX_ne_const W c ?_
  have hinj : Function.Injective (algebraMap W.CoordinateRing W.FunctionField) :=
    IsFractionRing.injective W.CoordinateRing W.FunctionField
  refine hinj ?_
  rw [← coordX, hc]
  rfl

/-! ## ★★★★★平行移動した座標 -/

/-- ★平行移動の割線の傾き(関数体の元として)。 -/
noncomputable def slopeFF (W : WeierstrassCurve.Affine F) (x₀ y₀ : F) : W.FunctionField :=
  (coordY W - algebraMap F W.FunctionField y₀) / (coordX W - algebraMap F W.FunctionField x₀)

/-- ★★**傾きは mathlib の `slope` と一致する**——`x` 座標が違うので割り算の場合になる。 -/
theorem slopeFF_eq (W : WeierstrassCurve.Affine F) (x₀ y₀ : F) :
    slopeFF W x₀ y₀
      = @WeierstrassCurve.Affine.slope _ _ (W.map (algebraMap F W.FunctionField))
          (Classical.decEq _) (coordX W) (algebraMap F W.FunctionField x₀)
          (coordY W) (algebraMap F W.FunctionField y₀) := by
  rw [@WeierstrassCurve.Affine.slope_of_X_ne _ _ _ (Classical.decEq _) _ _ _ _
    (coordX_ne_const W x₀)]
  rfl

/-- ★★★★★**平行移動した `x` 座標**——生成点 `+` 定数点 `Q`。 -/
noncomputable def translateX (W : WeierstrassCurve.Affine F) (x₀ y₀ : F) : W.FunctionField :=
  (W.map (algebraMap F W.FunctionField)).addX (coordX W)
    (algebraMap F W.FunctionField x₀) (slopeFF W x₀ y₀)

/-- ★★★★★**平行移動した `y` 座標**。 -/
noncomputable def translateY (W : WeierstrassCurve.Affine F) (x₀ y₀ : F) : W.FunctionField :=
  (W.map (algebraMap F W.FunctionField)).addY (coordX W)
    (algebraMap F W.FunctionField x₀) (coordY W) (slopeFF W x₀ y₀)

/-- ★★★★★★**平行移動した座標も曲線の非特異点である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `nonsingular_add` を生成点と定数点に適用しただけである
——**生成点が点であると分かったので群法則がそのまま効く**。 -/
theorem nonsingular_translate (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    (W.map (algebraMap F W.FunctionField)).Nonsingular
      (translateX W x₀ y₀) (translateY W x₀ y₀) := by
  letI : DecidableEq W.FunctionField := Classical.decEq _
  have hQ' : (W.map (algebraMap F W.FunctionField)).Nonsingular
      (algebraMap F W.FunctionField x₀) (algebraMap F W.FunctionField y₀) :=
    WeierstrassCurve.Affine.equation_iff_nonsingular.1
      ((hQ.1).map (algebraMap F W.FunctionField))
  have hne : ¬ (coordX W = algebraMap F W.FunctionField x₀
      ∧ coordY W = (W.map (algebraMap F W.FunctionField)).negY
          (algebraMap F W.FunctionField x₀) (algebraMap F W.FunctionField y₀)) :=
    fun h => coordX_ne_const W x₀ h.1
  have h := WeierstrassCurve.Affine.nonsingular_add (nonsingular_coord W) hQ' hne
  rw [translateX, translateY, slopeFF_eq]
  exact h

/-! ## ★★★★★★環準同型 -/

/-- ★★★★★★**平行移動の環準同型** `τ_Q : F[W] →+* F(W)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`nonsingular_translate` がそのまま `AdjoinRoot.lift` の仮説になる。 -/
noncomputable def translateHom (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) : W.CoordinateRing →+* W.FunctionField :=
  AdjoinRoot.lift
    (Polynomial.eval₂RingHom (algebraMap F W.FunctionField) (translateX W x₀ y₀))
    (translateY W x₀ y₀) (by
      rw [Polynomial.eval₂_eval₂RingHom_apply, ← WeierstrassCurve.Affine.map_polynomial]
      exact (nonsingular_translate W hQ).1)

/-- ★★★生成元の像は加法公式そのものである。 -/
theorem translateHom_genX (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    translateHom W hQ (genX W) = translateX W x₀ y₀ := by
  rw [translateHom, genX, CoordinateRing.mk, AdjoinRoot.lift_mk]
  simp

theorem translateHom_genY (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    translateHom W hQ (genY W) = translateY W x₀ y₀ := by
  rw [translateHom, genY, CoordinateRing.mk, AdjoinRoot.lift_mk]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def translateHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の環準同型)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
