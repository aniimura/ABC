/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.JScale
import Mathlib.AlgebraicGeometry.EllipticCurve.NormalForms

/-!
# 第 920 ブロック —— **★★★★★★★★★★★★二次の捧り（quadratic twist）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★これは何か

`Lemma 3.5` に残っているのは「`j` が同じ**分裂**乗法還元の対を 1 つ作る」
ことである（第 919 `minDeltaExp_eq_mul_of_twist`）。
★それを与えるのが**二次の捧り**である。

`a₁ = a₃ = 0` の形（mathlib の `IsCharNeTwoNF`、`2` が可逆ならいつでも取れる）で

    `W^d : a₂ ↦ d·a₂,  a₄ ↦ d²·a₄,  a₆ ↦ d³·a₆`

とすると、`c₄ ↦ d²c₄`・`c₆ ↦ d³c₆`・`Δ ↦ d⁶Δ` となり、
★★**`j` は変わらない**。

☆`j = Δ⁻¹c₄³`（第 913）に直しておけば、`d⁻⁶·d⁶ = 1` の計算だけである。

| 定義・定理 | 内容 |
|---|---|
| `quadTwist` | ★二次の捧り |
| `quadTwist_c₄` / `_c₆` / `_Δ` | 重み `d²`・`d³`・`d⁶` |
| `quadTwist_j` | ★★★★**`j` は変わらない** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve

variable {F : Type} [Field F]

/-- ★**二次の捧り**（`a₁ = a₃ = 0` の形で）。 -/
def quadTwist (W : WeierstrassCurve F) (d : F) : WeierstrassCurve F where
  a₁ := 0
  a₂ := d * W.a₂
  a₃ := 0
  a₄ := d ^ 2 * W.a₄
  a₆ := d ^ 3 * W.a₆

instance quadTwist_isCharNeTwoNF (W : WeierstrassCurve F) (d : F) :
    (quadTwist W d).IsCharNeTwoNF := ⟨rfl, rfl⟩

@[simp] theorem quadTwist_a₂ (W : WeierstrassCurve F) (d : F) :
    (quadTwist W d).a₂ = d * W.a₂ := rfl

@[simp] theorem quadTwist_a₄ (W : WeierstrassCurve F) (d : F) :
    (quadTwist W d).a₄ = d ^ 2 * W.a₄ := rfl

@[simp] theorem quadTwist_a₆ (W : WeierstrassCurve F) (d : F) :
    (quadTwist W d).a₆ = d ^ 3 * W.a₆ := rfl

/-- ★★`c₄` は `d²` 倍。 -/
theorem quadTwist_c₄ (W : WeierstrassCurve F) [W.IsCharNeTwoNF] (d : F) :
    (quadTwist W d).c₄ = d ^ 2 * W.c₄ := by
  rw [WeierstrassCurve.c₄_of_isCharNeTwoNF, WeierstrassCurve.c₄_of_isCharNeTwoNF,
    quadTwist_a₂, quadTwist_a₄]
  ring

/-- ★★`c₆` は `d³` 倍。 -/
theorem quadTwist_c₆ (W : WeierstrassCurve F) [W.IsCharNeTwoNF] (d : F) :
    (quadTwist W d).c₆ = d ^ 3 * W.c₆ := by
  rw [WeierstrassCurve.c₆_of_isCharNeTwoNF, WeierstrassCurve.c₆_of_isCharNeTwoNF,
    quadTwist_a₂, quadTwist_a₄, quadTwist_a₆]
  ring

/-- ★★`Δ` は `d⁶` 倍。 -/
theorem quadTwist_Δ (W : WeierstrassCurve F) [W.IsCharNeTwoNF] (d : F) :
    (quadTwist W d).Δ = d ^ 6 * W.Δ := by
  rw [WeierstrassCurve.Δ_of_isCharNeTwoNF, WeierstrassCurve.Δ_of_isCharNeTwoNF,
    quadTwist_a₂, quadTwist_a₄, quadTwist_a₆]
  ring

/-- ★★★捧りも楼円曲線である（`d ≠ 0` なら）。 -/
theorem quadTwist_isElliptic (W : WeierstrassCurve F) [W.IsCharNeTwoNF] [W.IsElliptic]
    {d : F} (hd : d ≠ 0) : (quadTwist W d).IsElliptic := by
  refine ⟨?_⟩
  rw [quadTwist_Δ]
  refine (isUnit_iff_ne_zero).2 ?_
  exact mul_ne_zero (pow_ne_zero 6 hd) W.isUnit_Δ.ne_zero

/-- ★★★★**捧りで `j` は変わらない**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`j = Δ⁻¹c₄³`（第 913）なので `(d⁶Δ)⁻¹(d²c₄)³ = Δ⁻¹c₄³` である。 -/
theorem quadTwist_j (W : WeierstrassCurve F) [W.IsCharNeTwoNF] [W.IsElliptic]
    {d : F} (hd : d ≠ 0) [(quadTwist W d).IsElliptic] :
    (quadTwist W d).j = W.j := by
  rw [j_eq_inv_Delta_mul, j_eq_inv_Delta_mul, quadTwist_Δ, quadTwist_c₄]
  have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  field_simp

/-! ## ★出典の紐付け(`.src`) -/

def quadTwist.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(二次の捧りの定義。★無条件)",
    sectionId := "genell-lemma-3-5" }

def quadTwist_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(捧りで j は変わらない。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
