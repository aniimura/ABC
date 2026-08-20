import ABC3.Found.GaloisRep.WeilBilin

/-!
# Galois (G5) 第 182 ブロック —— **★★★★★★和の点の生成元の関係式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★双線型性に要る関係式

第 181 の `aut_div_mul` は `g₁ g₂ = ζ·z·g₃` の形の関係を要求する。★その材料が

    f₁ · f₂ = u · h^n · f₃        (u ∈ F^×、h ∈ F(W)^×)

である。★★これは `toClass` の加法性から出る:

1. `toClass(P₁+P₂) = toClass P₁ + toClass P₂`(mathlib、`toClass` は `→+`)
2. したがって `U := XYIdeal'(P₁)·XYIdeal'(P₂)·XYIdeal'(P₁+P₂)⁻¹` の類は 1
3. 第 140 の `isPrincipal_of_classGroup_eq_one` で `U = (h)`
4. **`n` 乗して**元のレベルに落とす——`XYIdeal(P_i)^n = (f_i)` なので
   `spanSingleton(f₁f₂) = spanSingleton(h^n f₃)`、
   `spanSingleton_eq_spanSingleton` で `f₁f₂ = z·h^n·f₃`(`z` は `F[W]` の単元)
5. 第 128 の**単元は定数**から `z = u ∈ F^×`

★★★これは **Abel–Jacobi の元のレベルでの言い換え**であり、
D2 で使った道具(第 128・第 140)がそのまま効いている。

## ★★残る段

`μ̃` を当てて `μf₁·μf₂ = (定数)·μ̃(h)^n·μf₃` にし、`g_i^n = μf_i` から
`(g₁g₂/(μ̃(h) g₃))^n = 定数` ⟹ 第 177 の `const_of_pow_eq_const` で
`g₁g₂ = ζ·μ̃(h)·g₃`。★そこに第 181 の `aut_div_mul` を当てる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_h_of_add` | ★★★★★★イデアルのレベルの関係式 |
| `elem_relation_of_add` | ★★★★★★**元のレベルの関係式** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)

/-- ★★★★★★**和の点のイデアルは積と主イデアル 1 個だけ違う**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`toClass` の加法性(mathlib)と第 140 の `isPrincipal_of_classGroup_eq_one`。 -/
theorem exists_h_of_add
    {x₁ y₁ x₂ y₂ x₃ y₃ : F} (h₁ : W.Nonsingular x₁ y₁) (h₂ : W.Nonsingular x₂ y₂)
    (h₃ : W.Nonsingular x₃ y₃)
    (hsum : Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂ = Point.some x₃ y₃ h₃) :
    ∃ hh : W.FunctionField,
      (CoordinateRing.XYIdeal' h₁ : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
        * (CoordinateRing.XYIdeal' h₂)
        = FractionalIdeal.spanSingleton W.CoordinateRing⁰ hh
          * (CoordinateRing.XYIdeal' h₃) := by
  have hcls : ClassGroup.mk W.FunctionField (CoordinateRing.XYIdeal' h₁)
      * ClassGroup.mk W.FunctionField (CoordinateRing.XYIdeal' h₂)
      = ClassGroup.mk W.FunctionField (CoordinateRing.XYIdeal' h₃) := by
    have hq := congrArg Point.toClass hsum
    rw [map_add, Point.toClass_some, Point.toClass_some, Point.toClass_some] at hq
    exact hq
  set U := CoordinateRing.XYIdeal' h₁ * CoordinateRing.XYIdeal' h₂
    * (CoordinateRing.XYIdeal' h₃)⁻¹ with hU
  have hUone : ClassGroup.mk W.FunctionField U = 1 := by
    rw [hU, map_mul, map_mul, map_inv, hcls, mul_inv_cancel]
  obtain ⟨hh, hhs⟩ := isPrincipal_of_classGroup_eq_one U hUone
  refine ⟨hh, ?_⟩
  rw [← hhs, hU, Units.val_mul, Units.val_mul, mul_assoc, Units.inv_mul, mul_one]

/-- ★★★★★★**元のレベルの関係式**——`f₁ f₂ = u · h^n · f₃`(`u ∈ F^×`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★イデアルの等式を `n` 乗し、`spanSingleton_eq_spanSingleton` と
第 128 の**単元は定数**で元に落とす。 -/
theorem elem_relation_of_add {n : ℕ}
    {x₁ y₁ x₂ y₂ x₃ y₃ : F} (h₁ : W.Nonsingular x₁ y₁) (h₂ : W.Nonsingular x₂ y₂)
    (h₃ : W.Nonsingular x₃ y₃) {f₁ f₂ f₃ : W.CoordinateRing}
    (hf₁ : (CoordinateRing.XYIdeal W x₁ (Polynomial.C y₁)) ^ n = Ideal.span {f₁})
    (hf₂ : (CoordinateRing.XYIdeal W x₂ (Polynomial.C y₂)) ^ n = Ideal.span {f₂})
    (hf₃ : (CoordinateRing.XYIdeal W x₃ (Polynomial.C y₃)) ^ n = Ideal.span {f₃})
    {hh : W.FunctionField}
    (hid : (CoordinateRing.XYIdeal' h₁ : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
      * (CoordinateRing.XYIdeal' h₂)
      = FractionalIdeal.spanSingleton W.CoordinateRing⁰ hh * (CoordinateRing.XYIdeal' h₃)) :
    ∃ u : F, u ≠ 0 ∧
      (algebraMap W.CoordinateRing W.FunctionField f₁)
        * (algebraMap W.CoordinateRing W.FunctionField f₂)
        = algebraMap F W.FunctionField u * hh ^ n
          * (algebraMap W.CoordinateRing W.FunctionField f₃) := by
  have hcoe : ∀ {x y : F} (h : W.Nonsingular x y) {f : W.CoordinateRing},
      (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {f} →
      ((CoordinateRing.XYIdeal' h : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)) ^ n
        = FractionalIdeal.spanSingleton W.CoordinateRing⁰
          (algebraMap W.CoordinateRing W.FunctionField f) := by
    intro x y h f hf
    rw [CoordinateRing.XYIdeal'_eq, ← FractionalIdeal.coeIdeal_pow, hf,
      FractionalIdeal.coeIdeal_span_singleton]
  have hpow := congrArg (fun I => I ^ n) hid
  simp only [mul_pow, FractionalIdeal.spanSingleton_pow] at hpow
  rw [hcoe h₁ hf₁, hcoe h₂ hf₂, hcoe h₃ hf₃,
    FractionalIdeal.spanSingleton_mul_spanSingleton,
    FractionalIdeal.spanSingleton_mul_spanSingleton] at hpow
  obtain ⟨z, hz⟩ := FractionalIdeal.spanSingleton_eq_spanSingleton.1 hpow
  obtain ⟨u, hu0, hu⟩ := isUnit_coordinateRing z.isUnit
  refine ⟨u⁻¹, inv_ne_zero hu0, ?_⟩
  rw [Units.smul_def, Algebra.smul_def, hu,
    ← IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField u] at hz
  have hum : algebraMap F W.FunctionField u ≠ 0 := by simpa using hu0
  rw [map_inv₀]
  field_simp
  linear_combination hz

/-! ## ★出典の紐付け(`.src`) -/

def elem_relation_of_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——和の点の生成元の関係式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
