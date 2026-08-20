import ABC3.Found.GaloisRep.PointHom

/-!
# Galois (G5) 第 119 ブロック —— **★★★★★関数体の射は座標関数で決まる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★下流のすべてが要求する道具

Weil 対の残りの葉——**合成則** `τ_{Q+Q'} = τ_Q ∘ τ_{Q'}`、**全単射性**、
`[n]` の性質——はどれも「2 つの射が等しい」という形をしている。
★それを**生成元での一致**に落とす道具が要る:

    f, g : F(W) →ₐ[F] F(W),  f(x) = g(x),  f(y) = g(y)  ⟹  f = g

★★`F[W] = AdjoinRoot(W.polynomial)` は `F[X][Y]` の商なので、
`Polynomial.ringHom_ext` を 2 回使えば座標環の上で決まる。
★★★分数体へは `IsLocalization.ringHom_ext` で上がる。

## ★★★★併せて `F` 代数射として packaging する

第 118 の `pointHom` は `→+*` であった。★単射なら分数体へ延び、
しかも **`F` を固定する**ので `→ₐ[F]` になる(`pointFieldHom`)。
★★スケルトンの statement は `≃ₐ[F]` / `→ₐ[F]` で書かれているので、この形が要る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pointHom_algebraMap` | ★★`pointHom` は `F` を固定する |
| `pointFieldHom` | ★★★★★★**点 ⟹ `F(W) →ₐ[F] F(W)`** |
| `coordinateRing_hom_ext` | ★★★★座標環の射は生成元で決まる |
| `functionField_algHom_ext` | ★★★★★**関数体の `F` 代数射は座標関数で決まる** |
| `translateAlgHom` | ★★★★★★**平行移動の `F` 代数自己準同型** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★`F` を固定すること -/

/-- ★★`pointHom` は `F` 代数射である。 -/
theorem pointHom_algebraMap (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y) (c : F) :
    pointHom W h (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c := by
  have hc : algebraMap F W.CoordinateRing c
      = algebraMap (Polynomial F) W.CoordinateRing (Polynomial.C c) := rfl
  rw [hc, pointHom_polynomial, Polynomial.eval₂_C]

/-- ★★★★★★**点から関数体の `F` 代数自己準同型へ**(単射なとき)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
noncomputable def pointFieldHom (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hinj : Function.Injective (pointHom W h)) :
    W.FunctionField →ₐ[F] W.FunctionField where
  toRingHom := IsFractionRing.lift hinj
  commutes' := by
    intro c
    calc IsFractionRing.lift hinj (algebraMap F W.FunctionField c)
        = IsFractionRing.lift hinj (algebraMap W.CoordinateRing W.FunctionField
            (algebraMap F W.CoordinateRing c)) := by rw [← IsScalarTower.algebraMap_apply]
      _ = pointHom W h (algebraMap F W.CoordinateRing c) := IsFractionRing.lift_algebraMap _ _
      _ = algebraMap F W.FunctionField c := pointHom_algebraMap W h c

theorem pointFieldHom_coordX (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hinj : Function.Injective (pointHom W h)) :
    pointFieldHom W h hinj (coordX W) = x := by
  show IsFractionRing.lift hinj (coordX W) = x
  rw [coordX, IsFractionRing.lift_algebraMap, pointHom_genX]

theorem pointFieldHom_coordY (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hinj : Function.Injective (pointHom W h)) :
    pointFieldHom W h hinj (coordY W) = y := by
  show IsFractionRing.lift hinj (coordY W) = y
  rw [coordY, IsFractionRing.lift_algebraMap, pointHom_genY]

/-! ## ★★★★★射は生成元で決まる -/

/-- ★★★★**座標環からの環準同型は生成元で決まる**。

★`F[W] = AdjoinRoot(W.polynomial)` は `F[X][Y]` の商なので、
`Polynomial.ringHom_ext` を 2 回使う。 -/
theorem coordinateRing_hom_ext {W : WeierstrassCurve.Affine F} {S : Type} [CommRing S]
    (f g : W.CoordinateRing →+* S)
    (hc : ∀ a : F, f (algebraMap F W.CoordinateRing a) = g (algebraMap F W.CoordinateRing a))
    (hx : f (genX W) = g (genX W)) (hy : f (genY W) = g (genY W)) : f = g := by
  have hinner : (f.comp (CoordinateRing.mk W)).comp (Polynomial.C (R := Polynomial F))
      = (g.comp (CoordinateRing.mk W)).comp (Polynomial.C (R := Polynomial F)) := by
    refine Polynomial.ringHom_ext ?_ ?_
    · intro a
      exact hc a
    · exact hx
  have hcomp : f.comp (CoordinateRing.mk W) = g.comp (CoordinateRing.mk W) := by
    refine Polynomial.ringHom_ext ?_ ?_
    · intro q
      exact congrArg (fun h : Polynomial F →+* S => h q) hinner
    · exact hy
  refine RingHom.ext (fun z => ?_)
  obtain ⟨p, rfl⟩ := AdjoinRoot.mk_surjective (g := W.polynomial) z
  exact congrArg (fun h : Polynomial (Polynomial F) →+* S => h p) hcomp

/-- ★★★★★**関数体の `F` 代数自己準同型は座標関数で決まる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★合成則・全単射性はすべてこの補題を経由する。 -/
theorem functionField_algHom_ext {W : WeierstrassCurve.Affine F}
    (f g : W.FunctionField →ₐ[F] W.FunctionField)
    (hx : f (coordX W) = g (coordX W)) (hy : f (coordY W) = g (coordY W)) : f = g := by
  refine AlgHom.coe_ringHom_injective ?_
  refine IsLocalization.ringHom_ext (nonZeroDivisors W.CoordinateRing) ?_
  refine coordinateRing_hom_ext _ _ ?_ ?_ ?_
  · intro a
    show f (algebraMap W.CoordinateRing W.FunctionField (algebraMap F W.CoordinateRing a))
      = g (algebraMap W.CoordinateRing W.FunctionField (algebraMap F W.CoordinateRing a))
    rw [← IsScalarTower.algebraMap_apply, f.commutes, g.commutes]
  · exact hx
  · exact hy

/-! ## ★★★★★★平行移動の `F` 代数射 -/

theorem translateHom_eq_pointHom (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    translateHom W hQ = pointHom W (nonsingular_translate W hQ).1 := rfl

/-- ★★★★★★**平行移動の `F` 代数自己準同型**(`Q` が 2 等分点でないとき)。 -/
noncomputable def translateAlgHom (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    W.FunctionField →ₐ[F] W.FunctionField :=
  pointFieldHom W (nonsingular_translate W hQ).1
    (translateHom_eq_pointHom W hQ ▸ translateHom_injective W hQ h2)

theorem translateAlgHom_coordX (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    translateAlgHom W hQ h2 (coordX W) = translateX W x₀ y₀ :=
  pointFieldHom_coordX W _ _

theorem translateAlgHom_coordY (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    translateAlgHom W hQ h2 (coordY W) = translateY W x₀ y₀ :=
  pointFieldHom_coordY W _ _

/-! ## ★出典の紐付け(`.src`) -/

def functionField_algHom_ext.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——関数体の射が座標関数で決まること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
