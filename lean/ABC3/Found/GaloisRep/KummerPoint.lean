import ABC3.Found.GaloisRep.KummerX

/-!
# Galois (G5) 第 201 ブロック —— **★★★★★★★点の言葉での Kummer 公式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★梯子が消費する形にする

第 200 の 2 公式は `addX`・`slope` の言葉で書いてある。★Montgomery 梯子は
**点の加減**(`A + B` と `A − B`)で回るので、その形に直しておく。

### ★★★`A ± B` はどちらも `some` である

`x₁ ≠ x₂` なら `y₁ = negY x₂ y₂` は起きないので、`Point.add_some` が両方に効く。
★`A − B = A + (−B)` で、`−B = Point.some x₂ (negY x₂ y₂) _` は `rfl`。

### ★★★★★★取れる形

`A + B = (u, ·)`、`A − B = (v, ·)` として

    (u + v)·(x₁−x₂)² = 2x₁x₂(x₁+x₂) + b₂x₁x₂ + b₄(x₁+x₂) + b₆
    (u · v)·(x₁−x₂)² = x₁²x₂² − b₄x₁x₂ − b₆(x₁+x₂) − b₈

★★`x(A−B)` が既知なら **`x(A+B)` が積の式から一意に決まる**——これが梯子の 1 段である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `add_sub_some` | `A ± B` はどちらも `some`(`x₁ ≠ x₂` のとき) |
| `kummer_sum` | ★★★★★★**点の言葉での和の公式** |
| `kummer_prod` | ★★★★★★★**点の言葉での積の公式** |
-/

namespace ABC3.Found.GaloisRep

open Polynomial WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)

/-- `A + B` と `A − B` はどちらも `some` で、`x` 座標は `addX` である。 -/
theorem add_sub_some {x₁ y₁ x₂ y₂ : F} (h₁ : W.Nonsingular x₁ y₁) (h₂ : W.Nonsingular x₂ y₂)
    (hx : x₁ ≠ x₂) :
    (∃ h', Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂
        = Point.some (W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂))
          (W.addY x₁ x₂ y₁ (W.slope x₁ x₂ y₁ y₂)) h')
      ∧ (∃ h'', Point.some x₁ y₁ h₁ - Point.some x₂ y₂ h₂
        = Point.some (W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂)))
          (W.addY x₁ x₂ y₁ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂))) h'') := by
  have hn2 : W.Nonsingular x₂ (W.negY x₂ y₂) := (nonsingular_neg x₂ y₂).mpr h₂
  have hxy : ¬(x₁ = x₂ ∧ y₁ = W.negY x₂ y₂) := fun hc => hx hc.1
  have hxy' : ¬(x₁ = x₂ ∧ y₁ = W.negY x₂ (W.negY x₂ y₂)) := fun hc => hx hc.1
  refine ⟨⟨nonsingular_add h₁ h₂ hxy, Point.add_some hxy⟩,
    ⟨nonsingular_add h₁ hn2 hxy', ?_⟩⟩
  have hneg : -(Point.some x₂ y₂ h₂) = Point.some x₂ (W.negY x₂ y₂) hn2 := rfl
  rw [sub_eq_add_neg, hneg]
  exact Point.add_some hxy'

/-- ★★★★★★**点の言葉での Kummer の和の公式**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem kummer_sum {x₁ y₁ x₂ y₂ u yu v yv : F} (h₁ : W.Nonsingular x₁ y₁)
    (h₂ : W.Nonsingular x₂ y₂) (hx : x₁ ≠ x₂)
    (hu : W.Nonsingular u yu) (hv : W.Nonsingular v yv)
    (heu : Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂ = Point.some u yu hu)
    (hev : Point.some x₁ y₁ h₁ - Point.some x₂ y₂ h₂ = Point.some v yv hv) :
    (u + v) * (x₁ - x₂) ^ 2
      = 2 * x₁ * x₂ * (x₁ + x₂) + W.b₂ * x₁ * x₂ + W.b₄ * (x₁ + x₂) + W.b₆ := by
  obtain ⟨⟨ha, hae⟩, ⟨hb, hbe⟩⟩ := add_sub_some W h₁ h₂ hx
  rw [hae] at heu
  rw [hbe] at hev
  have hu' : W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂) = u := by injection heu
  have hv' : W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂)) = v := by injection hev
  rw [← hu', ← hv']
  exact addX_sum W h₁.1 h₂.1 hx

/-- ★★★★★★★**点の言葉での Kummer の積の公式**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`x(A−B)` が既知なら `x(A+B)` がこれで決まる——梯子の 1 段である。 -/
theorem kummer_prod (h4 : (4 : F) ≠ 0) {x₁ y₁ x₂ y₂ u yu v yv : F} (h₁ : W.Nonsingular x₁ y₁)
    (h₂ : W.Nonsingular x₂ y₂) (hx : x₁ ≠ x₂)
    (hu : W.Nonsingular u yu) (hv : W.Nonsingular v yv)
    (heu : Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂ = Point.some u yu hu)
    (hev : Point.some x₁ y₁ h₁ - Point.some x₂ y₂ h₂ = Point.some v yv hv) :
    (u * v) * (x₁ - x₂) ^ 2
      = x₁ ^ 2 * x₂ ^ 2 - W.b₄ * (x₁ * x₂) - W.b₆ * (x₁ + x₂) - W.b₈ := by
  obtain ⟨⟨ha, hae⟩, ⟨hb, hbe⟩⟩ := add_sub_some W h₁ h₂ hx
  rw [hae] at heu
  rw [hbe] at hev
  have hu' : W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂) = u := by injection heu
  have hv' : W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂)) = v := by injection hev
  rw [← hu', ← hv']
  exact addX_prod W h4 h₁.1 h₂.1 hx

/-! ## ★出典の紐付け(`.src`) -/

def kummer_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——点の言葉での Kummer の和の公式)",
    sectionId := "genell-thm-3-8" }

def kummer_prod.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——点の言葉での Kummer の積の公式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
