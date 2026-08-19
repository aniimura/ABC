import ABC3.Found.GaloisRep.TranslateEquiv

/-!
# Galois (G5) 第 121 ブロック —— **★★★★★★平行移動の合成則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`τ_{Q₁} ∘ τ_{Q₂} = τ_{Q₁+Q₂}`

第 120 で `τ_{−Q} ∘ τ_Q = id` を出した。★同じ道が**一般の合成**にそのまま伸びる:

    τ_{Q₁}(τ_{Q₂}(x)) = x((G + Q₁) + Q₂) = x(G + (Q₁+Q₂)) = τ_{Q₁+Q₂}(x)

★★真ん中の等号は mathlib の**結合則**である。

## ★★★★★これが 2 等分点の壁を回避する

第 117 の超越性は `−Q ≠ Q` を使うので、2 等分点では効かなかった。
★合成則があれば、2 等分点 `Q₃` を `Q₃ = Q₁ + Q₂`(どちらも 2 等分点でない)と
分解して

    τ_{Q₃} = τ_{Q₁} ∘ τ_{Q₂}

と書けるので、**単射性が合成で出る**(`translateHom_injective_of_decomp`)。
★★体の自己準同型は単射だから、`ψ = τ_{Q₁} ∘ τ_{Q₂}` は単射であり、
`translateHom W h₃ = ψ ∘ (F[W] ↪ F(W))` も単射になる。

★★★残るのは**分解の存在**だけである——`Skeleton` に葉として記録する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `translateAlgHom_comp_general` | ★★★★★★**合成則** |
| `translateHom_injective_of_decomp` | ★★★★★★**分解があれば単射**(2 等分点でも) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★★★★**平行移動の合成則**——`τ_{Q₁}(τ_{Q₂}(x)) = τ_{Q₁+Q₂}(x)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`(G + Q₁) + Q₂ = G + (Q₁ + Q₂)` という mathlib の結合則を座標で読んだものである。 -/
theorem translateAlgHom_comp_general (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₁ y₁ x₂ y₂ x₃ y₃ : F}
    (h₁ : W.Nonsingular x₁ y₁) (h₂ : W.Nonsingular x₂ y₂) (h₃ : W.Nonsingular x₃ y₃)
    (k₁ : W.negY x₁ y₁ ≠ y₁)
    (hsum : Point.some (algebraMap F W.FunctionField x₁) (algebraMap F W.FunctionField y₁)
          (mapNonsingular W h₁)
        + Point.some (algebraMap F W.FunctionField x₂) (algebraMap F W.FunctionField y₂)
          (mapNonsingular W h₂)
        = Point.some (algebraMap F W.FunctionField x₃) (algebraMap F W.FunctionField y₃)
          (mapNonsingular W h₃)) :
    translateAlgHom W h₁ k₁ (translateX W x₂ y₂) = translateX W x₃ y₃
    ∧ translateAlgHom W h₁ k₁ (translateY W x₂ y₂) = translateY W x₃ y₃ := by
  have hne : translateX W x₁ y₁ ≠ algebraMap F W.FunctionField x₂ :=
    translateX_ne_const W h₁.1 k₁ x₂
  have hpt : Point.some (translateX W x₁ y₁) (translateY W x₁ y₁) (nonsingular_translate W h₁)
      + Point.some (algebraMap F W.FunctionField x₂) (algebraMap F W.FunctionField y₂)
        (mapNonsingular W h₂)
      = Point.some (translateX W x₃ y₃) (translateY W x₃ y₃) (nonsingular_translate W h₃) := by
    rw [← genericPoint_add_const W h₁, add_assoc, hsum, genericPoint_add_const W h₃]
  rw [WeierstrassCurve.Affine.Point.add_of_X_ne hne,
    WeierstrassCurve.Affine.Point.some.injEq] at hpt
  constructor
  · rw [map_translateX, translateAlgHom_coordX, translateAlgHom_coordY,
      show (translateY W x₁ y₁ - algebraMap F W.FunctionField y₂)
          / (translateX W x₁ y₁ - algebraMap F W.FunctionField x₂)
        = (W.map (algebraMap F W.FunctionField)).slope (translateX W x₁ y₁)
          (algebraMap F W.FunctionField x₂) (translateY W x₁ y₁)
          (algebraMap F W.FunctionField y₂) from
      (WeierstrassCurve.Affine.slope_of_X_ne hne).symm]
    exact hpt.1
  · rw [map_translateY, translateAlgHom_coordX, translateAlgHom_coordY,
      show (translateY W x₁ y₁ - algebraMap F W.FunctionField y₂)
          / (translateX W x₁ y₁ - algebraMap F W.FunctionField x₂)
        = (W.map (algebraMap F W.FunctionField)).slope (translateX W x₁ y₁)
          (algebraMap F W.FunctionField x₂) (translateY W x₁ y₁)
          (algebraMap F W.FunctionField y₂) from
      (WeierstrassCurve.Affine.slope_of_X_ne hne).symm]
    exact hpt.2

/-- ★★★★★★**分解 `Q₃ = Q₁ + Q₂` があれば `τ_{Q₃}` は単射である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Q₃` が 2 等分点でも使える——第 117 の超越性が効かない場合の回避路である。
★★体の自己準同型は単射なので、合成 `ψ = τ_{Q₁} ∘ τ_{Q₂}` は単射であり、
`translateHom W h₃ = ψ ∘ (F[W] ↪ F(W))` も単射になる。 -/
theorem translateHom_injective_of_decomp (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₁ y₁ x₂ y₂ x₃ y₃ : F}
    (h₁ : W.Nonsingular x₁ y₁) (h₂ : W.Nonsingular x₂ y₂) (h₃ : W.Nonsingular x₃ y₃)
    (k₁ : W.negY x₁ y₁ ≠ y₁) (k₂ : W.negY x₂ y₂ ≠ y₂)
    (hsum : Point.some (algebraMap F W.FunctionField x₁) (algebraMap F W.FunctionField y₁)
          (mapNonsingular W h₁)
        + Point.some (algebraMap F W.FunctionField x₂) (algebraMap F W.FunctionField y₂)
          (mapNonsingular W h₂)
        = Point.some (algebraMap F W.FunctionField x₃) (algebraMap F W.FunctionField y₃)
          (mapNonsingular W h₃)) :
    Function.Injective (translateHom W h₃) := by
  set ψ := (translateAlgHom W h₁ k₁).comp (translateAlgHom W h₂ k₂) with hψ
  obtain ⟨hx, hy⟩ := translateAlgHom_comp_general W h₁ h₂ h₃ k₁ hsum
  have hψx : ψ (coordX W) = translateX W x₃ y₃ := by
    rw [hψ, AlgHom.comp_apply, translateAlgHom_coordX]
    exact hx
  have hψy : ψ (coordY W) = translateY W x₃ y₃ := by
    rw [hψ, AlgHom.comp_apply, translateAlgHom_coordY]
    exact hy
  have heq : translateHom W h₃
      = (ψ : W.FunctionField →+* W.FunctionField).comp
        (algebraMap W.CoordinateRing W.FunctionField) := by
    refine coordinateRing_hom_ext _ _ ?_ ?_ ?_
    · intro a
      rw [translateHom_eq_pointHom, pointHom_algebraMap]
      show _ = ψ (algebraMap W.CoordinateRing W.FunctionField (algebraMap F W.CoordinateRing a))
      rw [← IsScalarTower.algebraMap_apply, ψ.commutes]
    · rw [translateHom_genX]
      show _ = ψ (coordX W)
      rw [hψx]
    · rw [translateHom_genY]
      show _ = ψ (coordY W)
      rw [hψy]
  rw [heq, RingHom.coe_comp]
  exact Function.Injective.comp (RingHom.injective (ψ : W.FunctionField →+* W.FunctionField))
    (IsFractionRing.injective W.CoordinateRing W.FunctionField)

/-! ## ★出典の紐付け(`.src`) -/

def translateAlgHom_comp_general.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の合成則)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
