/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Ker
import Mathlib.Algebra.Group.Commutator
import Mathlib.GroupTheory.Commutator.Basic
import Mathlib.Tactic.Group

/-!
# [FrdI] `Proposition 4.4, (ii)` —— 中心拡大の交換子が誘導する交代形式

原文 (FrdI p.83):
> is a commutative monoid, and that the natural homomorphism

★★`Proposition 4.4, (ii)` は `𝒪^×(A^birat)` が**可換**であることを言う。
在庫で次の形まで来ている(`Prop44Ker.lean`):

```
1 → (𝒪^×(A) の像) → 𝒪^×(A^birat) → Φ^birat(A) → 0
```

★`Φ^birat` は**加法群の部分群**なので**可換**である。
したがって `𝒪^×(A^birat)` の**交換子は核に落ちる**。

## ★★本ファイルの担当 —— 一般論の側を先に固める

★★★交換子が交代双線形形式を誘導するのは**純粋な群論**であり、
Frobenioid とは無関係である。先にそちらを一般の形で取っておけば、
`pairing-vanishes`(原典が「routine exercise」と書いた所)に
**Frobenioid 固有の議論だけ**が残る。

## ★測ったこと(2026-08-19)

★台帳は `image-central`(核が中心に入る)を「済」としているが、
★★**Lean のポインタが記録されていない**。本ファイルの
`commutatorElement_mul_left_of_central` などは
**「交換子が中心に入る」を仮定として明示的に受ける**形にしてあるので、
`image-central` の実体が確定した時点でそのまま繋がる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped commutatorElement

/-! ## ★1. 一般の群論 —— 交換子が中心に落ちるときの双線形性 -/

section CentralPairing

variable {G : Type*} [Group G]

/-- ★**可換群への準同型の交換子は核に落ちる**。 -/
theorem commutatorElement_mem_ker {Q : Type*} [CommGroup Q] (f : G →* Q) (a b : G) :
    ⁅a, b⁆ ∈ f.ker := by
  rw [MonoidHom.mem_ker, commutatorElement_def]
  rw [map_mul, map_mul, map_mul, map_inv, map_inv]
  rw [mul_comm (f a) (f b)]
  group

/-- ★**交代性** —— `⁅a, a⁆ = 1`。 -/
theorem commutatorElement_self' (a : G) : ⁅a, a⁆ = 1 := by
  rw [commutatorElement_def]; group

/-- ★★**中心の元は交換子から消える**(左)。 -/
theorem commutatorElement_center_left {z b : G} (hz : z ∈ Subgroup.center G) :
    ⁅z, b⁆ = 1 := by
  have h : b * z = z * b := Subgroup.mem_center_iff.mp hz b
  rw [commutatorElement_def, ← h]
  group

/-- ★★**中心の元は交換子から消える**(右)。 -/
theorem commutatorElement_center_right {a z : G} (hz : z ∈ Subgroup.center G) :
    ⁅a, z⁆ = 1 := by
  have h : a * z = z * a := Subgroup.mem_center_iff.mp hz a
  rw [commutatorElement_def, h]
  group

/-- ★★★**第 1 変数の乗法性** —— `⁅b, c⁆` が中心に入るときに成り立つ。

★★これが「交換子が**双線形**である」の第 1 の柱である。 -/
theorem commutatorElement_mul_left_of_central {a b c : G}
    (hbc : ⁅b, c⁆ ∈ Subgroup.center G) :
    ⁅a * b, c⁆ = ⁅a, c⁆ * ⁅b, c⁆ := by
  have hw : ∀ g : G, g * ⁅b, c⁆ = ⁅b, c⁆ * g := Subgroup.mem_center_iff.mp hbc
  have hbcb : b * c * b⁻¹ = ⁅b, c⁆ * c := by rw [commutatorElement_def]; group
  have e1 : ⁅a * b, c⁆ = a * (b * c * b⁻¹) * a⁻¹ * c⁻¹ := by
    rw [commutatorElement_def]; group
  have e2 : a * (⁅b, c⁆ * c) * a⁻¹ * c⁻¹ = ⁅b, c⁆ * (a * c * a⁻¹ * c⁻¹) := by
    rw [← mul_assoc a ⁅b, c⁆ c, hw a]
    simp only [mul_assoc]
  calc ⁅a * b, c⁆ = a * (b * c * b⁻¹) * a⁻¹ * c⁻¹ := e1
    _ = a * (⁅b, c⁆ * c) * a⁻¹ * c⁻¹ := by rw [hbcb]
    _ = ⁅b, c⁆ * (a * c * a⁻¹ * c⁻¹) := e2
    _ = (a * c * a⁻¹ * c⁻¹) * ⁅b, c⁆ := (hw _).symm
    _ = ⁅a, c⁆ * ⁅b, c⁆ := by simp only [commutatorElement_def]

/-- ★`⁅x, y⁆ = ⁅y, x⁆⁻¹`。 -/
theorem commutatorElement_symm_inv (x y : G) : ⁅x, y⁆ = (⁅y, x⁆)⁻¹ := by
  rw [commutatorElement_def, commutatorElement_def]; group

/-- ★★★**第 2 変数の乗法性**。★`⁅a, b⁆ = ⁅b, a⁆⁻¹` で第 1 変数へ帰着する。 -/
theorem commutatorElement_mul_right_of_central {a b c : G}
    (hab : ⁅a, b⁆ ∈ Subgroup.center G) (hac : ⁅a, c⁆ ∈ Subgroup.center G) :
    ⁅a, b * c⁆ = ⁅a, b⁆ * ⁅a, c⁆ := by
  have hca : ⁅c, a⁆ ∈ Subgroup.center G := by
    rw [commutatorElement_symm_inv c a]
    exact (Subgroup.center G).inv_mem hac
  rw [commutatorElement_symm_inv a (b * c),
    commutatorElement_mul_left_of_central hca, mul_inv_rev,
    ← commutatorElement_symm_inv a c, ← commutatorElement_symm_inv a b]
  exact (Subgroup.mem_center_iff.mp hac _).symm

/-- ★★**中心の元をずらしても交換子は変わらない**(第 1 変数)。

★★★これが「交換子が**商の上で well-defined**」の中身である。 -/
theorem commutatorElement_mul_center_left {a z b : G} (hz : z ∈ Subgroup.center G) :
    ⁅a * z, b⁆ = ⁅a, b⁆ := by
  have h1 : ⁅z, b⁆ = 1 := commutatorElement_center_left hz
  have h2 : ⁅z, b⁆ ∈ Subgroup.center G := by rw [h1]; exact Subgroup.one_mem _
  rw [commutatorElement_mul_left_of_central h2, h1, mul_one]

/-- ★★**中心の元をずらしても交換子は変わらない**(第 2 変数)。 -/
theorem commutatorElement_mul_center_right {a b z : G} (hz : z ∈ Subgroup.center G)
    (hab : ⁅a, b⁆ ∈ Subgroup.center G) :
    ⁅a, b * z⁆ = ⁅a, b⁆ := by
  have h1 : ⁅a, z⁆ = 1 := commutatorElement_center_right hz
  have h2 : ⁅a, z⁆ ∈ Subgroup.center G := by rw [h1]; exact Subgroup.one_mem _
  rw [commutatorElement_mul_right_of_central hab h2, h1, mul_one]

/-! ### ★★★まとめ —— 核が中心に入る中心拡大では交換子は交代双線形

★`hZ : f.ker ≤ Subgroup.center G` を 1 本置けば、
`commutatorElement_mem_ker` と合わせて**すべての交換子が中心に入る**ので、
上の 4 本の仮定が自動的に満たされる。 -/

variable {Q : Type*} [CommGroup Q] (f : G →* Q)

/-- ★★★**核が中心なら、すべての交換子が中心に入る**。 -/
theorem commutatorElement_mem_center_of_ker_le_center
    (hZ : f.ker ≤ Subgroup.center G) (a b : G) : ⁅a, b⁆ ∈ Subgroup.center G :=
  hZ (commutatorElement_mem_ker f a b)

variable (hZ : f.ker ≤ Subgroup.center G)
include hZ

/-- ★★★★**交代双線形性 (1)** —— 第 1 変数で乗法的。 -/
theorem commutatorElement_mul_left (a b c : G) : ⁅a * b, c⁆ = ⁅a, c⁆ * ⁅b, c⁆ :=
  commutatorElement_mul_left_of_central
    (commutatorElement_mem_center_of_ker_le_center f hZ b c)

/-- ★★★★**交代双線形性 (2)** —— 第 2 変数で乗法的。 -/
theorem commutatorElement_mul_right (a b c : G) : ⁅a, b * c⁆ = ⁅a, b⁆ * ⁅a, c⁆ :=
  commutatorElement_mul_right_of_central
    (commutatorElement_mem_center_of_ker_le_center f hZ a b)
    (commutatorElement_mem_center_of_ker_le_center f hZ a c)

omit [CommGroup Q] hZ in
/-- ★★★★**交代双線形性 (3)** —— 交代的。 -/
theorem commutatorElement_alternating (a : G) : ⁅a, a⁆ = 1 := commutatorElement_self' a

/-- ★★★★**商の上で well-defined** —— `f a₁ = f a₂`、`f b₁ = f b₂` なら交換子は等しい。

★★★これが「交換子が `Q` 上の交代双線形形式を誘導する」の中身である。 -/
theorem commutatorElement_congr {a₁ a₂ b₁ b₂ : G}
    (ha : f a₁ = f a₂) (hb : f b₁ = f b₂) : ⁅a₁, b₁⁆ = ⁅a₂, b₂⁆ := by
  have hza : a₁⁻¹ * a₂ ∈ f.ker := by
    rw [MonoidHom.mem_ker, map_mul, map_inv, ha, inv_mul_cancel]
  have hzb : b₁⁻¹ * b₂ ∈ f.ker := by
    rw [MonoidHom.mem_ker, map_mul, map_inv, hb, inv_mul_cancel]
  have hA : a₂ = a₁ * (a₁⁻¹ * a₂) := by group
  have hB : b₂ = b₁ * (b₁⁻¹ * b₂) := by group
  rw [hA, hB, commutatorElement_mul_center_left (hZ hza),
    commutatorElement_mul_center_right (hZ hzb)
      (commutatorElement_mem_center_of_ker_le_center f hZ _ _)]

/-- ★★★**交換子は核に値を取る** —— 誘導される形式の行き先。 -/
theorem commutatorElement_mem_ker' (a b : G) : ⁅a, b⁆ ∈ f.ker :=
  commutatorElement_mem_ker f a b

end CentralPairing

def commutatorElement_congr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (ii) — 中心拡大の交換子が誘導する交代双線形形式",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
