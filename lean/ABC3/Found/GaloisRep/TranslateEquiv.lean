import ABC3.Found.GaloisRep.FieldHomExt

/-!
# Galois (G5) 第 120 ブロック —— **★★★★★★平行移動は関数体の自己同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★葉「全単射性」が閉じた

第 119 の `functionField_algHom_ext` により、`τ_{−Q} ∘ τ_Q = id` は
**生成元 `x`・`y` での等式**に落ちた。★それを群法則で片付ける:

    G + Z = (平行移動した点)      (`genericPoint_add_const`)
    (G + Z) + Q = G + (Z + Q) = G  (`Z = −Q`、mathlib の結合則)

★★これを座標で読むと `τ_Z(τ_Q(x)) = x`、`τ_Z(τ_Q(y)) = y` になる。

## ★★★★★鍵は 3 つ

| 補題 | 内容 |
|---|---|
| `map_translateX` / `map_translateY` | ★★★`F` 代数射は加法公式と可換 |
| `genericPoint_add_const` | ★★★★★**生成点 + 定数点 = 平行移動した点** |
| `translate_neg_add'` | ★★★★★★**`(G + Z) + Q = G`** |

★`add_of_X_ne`(mathlib)を使うために `x` 座標が異なることが要るが、
それは第 117 の**超越性**から出る(`translateX_ne_const`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `genericPoint_add_const` | ★★★★★生成点 + 定数点 |
| `translate_neg_add'` | ★★★★★★`(G + Z) + Q = G` |
| `translateAlgHom_comp` | ★★★★★★**合成が恒等** |
| `translateAlgEquiv` | ★★★★★★**平行移動の自己同型** |
| `exists_translateAut_of_not_twoTorsion` | ★★★★★★**葉が閉じた**(2 等分点以外) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★★`F` 代数射は加法公式と可換 -/

/-- ★★★**`F` 代数射は加法公式と可換である**(`x` 側)。 -/
theorem map_translateX (W : WeierstrassCurve.Affine F)
    (τ : W.FunctionField →ₐ[F] W.FunctionField) (x₀ y₀ : F) :
    τ (translateX W x₀ y₀)
      = (W.map (algebraMap F W.FunctionField)).addX (τ (coordX W))
          (algebraMap F W.FunctionField x₀)
          ((τ (coordY W) - algebraMap F W.FunctionField y₀)
            / (τ (coordX W) - algebraMap F W.FunctionField x₀)) := by
  rw [translateX, WeierstrassCurve.Affine.addX, WeierstrassCurve.Affine.addX, slopeFF]
  simp only [map_sub, map_add, map_mul, map_pow, map_div₀, AlgHom.commutes,
    WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂]

/-- ★★★**`F` 代数射は加法公式と可換である**(`y` 側)。 -/
theorem map_translateY (W : WeierstrassCurve.Affine F)
    (τ : W.FunctionField →ₐ[F] W.FunctionField) (x₀ y₀ : F) :
    τ (translateY W x₀ y₀)
      = (W.map (algebraMap F W.FunctionField)).addY (τ (coordX W))
          (algebraMap F W.FunctionField x₀) (τ (coordY W))
          ((τ (coordY W) - algebraMap F W.FunctionField y₀)
            / (τ (coordX W) - algebraMap F W.FunctionField x₀)) := by
  rw [translateY, WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.addY,
    WeierstrassCurve.Affine.negY, WeierstrassCurve.Affine.negY,
    WeierstrassCurve.Affine.negAddY, WeierstrassCurve.Affine.negAddY,
    WeierstrassCurve.Affine.addX, WeierstrassCurve.Affine.addX, slopeFF]
  simp only [map_sub, map_add, map_mul, map_pow, map_div₀, map_neg, AlgHom.commutes,
    WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₃]

/-! ## ★★★★★生成点の加法 -/

theorem slopeFF_eq' (W : WeierstrassCurve.Affine F) [DecidableEq W.FunctionField] (x₀ y₀ : F) :
    slopeFF W x₀ y₀ = (W.map (algebraMap F W.FunctionField)).slope (coordX W)
      (algebraMap F W.FunctionField x₀) (coordY W) (algebraMap F W.FunctionField y₀) := by
  rw [WeierstrassCurve.Affine.slope_of_X_ne (coordX_ne_const W x₀)]
  rfl

theorem mapNonsingular (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x₀ y₀ : F}
    (hQ : W.Nonsingular x₀ y₀) :
    (W.map (algebraMap F W.FunctionField)).Nonsingular
      (algebraMap F W.FunctionField x₀) (algebraMap F W.FunctionField y₀) :=
  WeierstrassCurve.Affine.equation_iff_nonsingular.1 ((hQ.1).map (algebraMap F W.FunctionField))

theorem point_some_congr {K : Type} [Field K] {V : WeierstrassCurve.Affine K} {x y x' y' : K}
    (hx : x = x') (hy : y = y') (h : V.Nonsingular x y) (h' : V.Nonsingular x' y') :
    Point.some x y h = Point.some x' y' h' := by
  subst hx
  subst hy
  rfl

theorem map_negY (W : WeierstrassCurve.Affine F) (x y : F) :
    (W.map (algebraMap F W.FunctionField)).negY
        (algebraMap F W.FunctionField x) (algebraMap F W.FunctionField y)
      = algebraMap F W.FunctionField (W.negY x y) := by
  rw [WeierstrassCurve.Affine.negY, WeierstrassCurve.Affine.negY]
  simp only [map_sub, map_neg, map_mul, WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₃]

theorem translateX_ne_const (W : WeierstrassCurve.Affine F) {x₀ y₀ : F}
    (hQ : W.Equation x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) (c : F) :
    translateX W x₀ y₀ ≠ algebraMap F W.FunctionField c := by
  intro hc
  have hp : (Polynomial.X - Polynomial.C c) ≠ (0 : Polynomial F) := by
    intro h
    have hd := congrArg Polynomial.natDegree h
    simp at hd
  refine translateX_transcendental W hQ h2 hp ?_
  rw [Polynomial.eval₂_sub, Polynomial.eval₂_X, Polynomial.eval₂_C, hc, sub_self]

/-- ★★★★★**生成点 + 定数点 = 平行移動した点**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが加法公式と mathlib の群法則を繋ぐ橋である。 -/
theorem genericPoint_add_const (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x₀ y₀ : F}
    (hQ : W.Nonsingular x₀ y₀) :
    genericPoint W + Point.some (algebraMap F W.FunctionField x₀)
        (algebraMap F W.FunctionField y₀) (mapNonsingular W hQ)
      = Point.some (translateX W x₀ y₀) (translateY W x₀ y₀) (nonsingular_translate W hQ) := by
  rw [genericPoint, WeierstrassCurve.Affine.Point.add_of_X_ne (coordX_ne_const W x₀)]
  refine point_some_congr ?_ ?_ _ _
  · rw [translateX, slopeFF_eq']
  · rw [translateY, slopeFF_eq']

/-- ★★★★★★**`(G + Z) + Q = G`**(`Z = −Q`)。 -/
theorem translate_neg_add' (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x₀ y₀ z₀ : F}
    (hz : W.negY x₀ y₀ = z₀) (hQ : W.Nonsingular x₀ y₀) (hZ : W.Nonsingular x₀ z₀) :
    Point.some (translateX W x₀ z₀) (translateY W x₀ z₀) (nonsingular_translate W hZ)
      + Point.some (algebraMap F W.FunctionField x₀) (algebraMap F W.FunctionField y₀)
        (mapNonsingular W hQ)
      = genericPoint W := by
  have hnegQ : Point.some (algebraMap F W.FunctionField x₀)
      (algebraMap F W.FunctionField z₀) (mapNonsingular W hZ)
      = -Point.some (algebraMap F W.FunctionField x₀) (algebraMap F W.FunctionField y₀)
        (mapNonsingular W hQ) := by
    rw [WeierstrassCurve.Affine.Point.neg_some]
    refine point_some_congr rfl ?_ _ _
    rw [map_negY, hz]
  rw [← genericPoint_add_const W hZ, hnegQ, add_assoc, neg_add_cancel, add_zero]

/-! ## ★★★★★★合成が恒等 -/

theorem translateAlgHom_translateX' (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ z₀ : F} (hz : W.negY x₀ y₀ = z₀) (hQ : W.Nonsingular x₀ y₀)
    (hZ : W.Nonsingular x₀ z₀) (h2Z : W.negY x₀ z₀ ≠ z₀) :
    translateAlgHom W hZ h2Z (translateX W x₀ y₀) = coordX W := by
  rw [map_translateX, translateAlgHom_coordX, translateAlgHom_coordY]
  have hne : translateX W x₀ z₀ ≠ algebraMap F W.FunctionField x₀ :=
    translateX_ne_const W hZ.1 h2Z x₀
  rw [show (translateY W x₀ z₀ - algebraMap F W.FunctionField y₀)
        / (translateX W x₀ z₀ - algebraMap F W.FunctionField x₀)
      = (W.map (algebraMap F W.FunctionField)).slope (translateX W x₀ z₀)
        (algebraMap F W.FunctionField x₀) (translateY W x₀ z₀)
        (algebraMap F W.FunctionField y₀) from
    (WeierstrassCurve.Affine.slope_of_X_ne hne).symm]
  have hsum := translate_neg_add' W hz hQ hZ
  rw [WeierstrassCurve.Affine.Point.add_of_X_ne hne, genericPoint,
    WeierstrassCurve.Affine.Point.some.injEq] at hsum
  exact hsum.1

theorem translateAlgHom_translateY' (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ z₀ : F} (hz : W.negY x₀ y₀ = z₀) (hQ : W.Nonsingular x₀ y₀)
    (hZ : W.Nonsingular x₀ z₀) (h2Z : W.negY x₀ z₀ ≠ z₀) :
    translateAlgHom W hZ h2Z (translateY W x₀ y₀) = coordY W := by
  rw [map_translateY, translateAlgHom_coordX, translateAlgHom_coordY]
  have hne : translateX W x₀ z₀ ≠ algebraMap F W.FunctionField x₀ :=
    translateX_ne_const W hZ.1 h2Z x₀
  rw [show (translateY W x₀ z₀ - algebraMap F W.FunctionField y₀)
        / (translateX W x₀ z₀ - algebraMap F W.FunctionField x₀)
      = (W.map (algebraMap F W.FunctionField)).slope (translateX W x₀ z₀)
        (algebraMap F W.FunctionField x₀) (translateY W x₀ z₀)
        (algebraMap F W.FunctionField y₀) from
    (WeierstrassCurve.Affine.slope_of_X_ne hne).symm]
  have hsum := translate_neg_add' W hz hQ hZ
  rw [WeierstrassCurve.Affine.Point.add_of_X_ne hne, genericPoint,
    WeierstrassCurve.Affine.Point.some.injEq] at hsum
  exact hsum.2

/-- ★★★★★★**`τ_Z ∘ τ_Q = id`**(`Z = −Q`)。 -/
theorem translateAlgHom_comp (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x₀ y₀ z₀ : F}
    (hz : W.negY x₀ y₀ = z₀) (hQ : W.Nonsingular x₀ y₀) (hZ : W.Nonsingular x₀ z₀)
    (h2Q : W.negY x₀ y₀ ≠ y₀) (h2Z : W.negY x₀ z₀ ≠ z₀) :
    (translateAlgHom W hZ h2Z).comp (translateAlgHom W hQ h2Q)
      = AlgHom.id F W.FunctionField := by
  refine functionField_algHom_ext _ _ ?_ ?_
  · simp only [AlgHom.comp_apply, AlgHom.coe_id, id_eq]
    rw [translateAlgHom_coordX]
    exact translateAlgHom_translateX' W hz hQ hZ h2Z
  · simp only [AlgHom.comp_apply, AlgHom.coe_id, id_eq]
    rw [translateAlgHom_coordY]
    exact translateAlgHom_translateY' W hz hQ hZ h2Z

/-! ## ★★★★★★自己同型 -/

/-- ★★★★★★**平行移動は関数体の `F` 自己同型である**(`Q` が 2 等分点でないとき)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★逆は `τ_{−Q}` である。 -/
noncomputable def translateAlgEquiv (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x₀ y₀ : F}
    (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    W.FunctionField ≃ₐ[F] W.FunctionField :=
  AlgEquiv.ofAlgHom (translateAlgHom W hQ h2)
    (translateAlgHom W ((WeierstrassCurve.Affine.nonsingular_neg x₀ y₀).2 hQ)
      (by rw [WeierstrassCurve.Affine.negY_negY]; exact Ne.symm h2))
    (translateAlgHom_comp W (WeierstrassCurve.Affine.negY_negY x₀ y₀) _ _ _ _)
    (translateAlgHom_comp W rfl hQ _ _ _)

theorem translateAlgEquiv_coordX (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x₀ y₀ : F}
    (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    translateAlgEquiv W hQ h2 (coordX W) = translateX W x₀ y₀ :=
  translateAlgHom_coordX W hQ h2

theorem translateAlgEquiv_coordY (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x₀ y₀ : F}
    (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    translateAlgEquiv W hQ h2 (coordY W) = translateY W x₀ y₀ :=
  translateAlgHom_coordY W hQ h2

/-- ★★★★★★**葉が閉じた**——平行移動 `τ_Q` は関数体の `F` 自己同型を誘導する
(`Q` が 2 等分点でないとき)。 -/
theorem exists_translateAut_of_not_twoTorsion (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    ∃ τ : W.FunctionField ≃ₐ[F] W.FunctionField,
      τ (coordX W) = translateX W x₀ y₀ ∧ τ (coordY W) = translateY W x₀ y₀ :=
  ⟨translateAlgEquiv W hQ h2, translateAlgEquiv_coordX W hQ h2, translateAlgEquiv_coordY W hQ h2⟩

/-! ## ★出典の紐付け(`.src`) -/

def translateAlgEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動が関数体の自己同型であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
