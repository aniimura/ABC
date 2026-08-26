/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63Model

/-!
# `Example 6.3` の射の記述

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★原文が「one verifies immediately」で畳んだ所

原文は `Example 6.3` の Frobenioid の射を
「(a) `Spec L → Spec M`、(b) `d ∈ ℕ≥1`、(c) 非零な元」の組だと述べる。
★これは `Theorem 5.2, (i)` の `ModelData.Hom` の 4 成分そのものである:

| 原文 | `ModelData.Hom` の成分 | `Example 6.3` での中身 |
|---|---|---|
| (a) `Spec L → Spec M` | `base` | `𝒟 = B(G)⁰` の射(= `F`-代数射 `M → L`) |
| (b) `d ∈ ℕ≥1` | `deg` | `ℕ+` |
| (c) 非零な元 | `u` | `L^×`(加法的に書いたもの) |
| — | `div` | 有効な算術因子 `Φ(L)` の元 |
| — | `cond` | `d·α_A + Div(φ) = Base(φ)^*(α_B) + div(u)` |

★★**`ex63HomEquiv` がその 1 対 1 対応**である。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `ex63HomEquiv` | ★★★**射 ⟺ (底の射, 因子, 次数, 単元) の組で関係式を満たすもの** |
| `ex63_hom_unit_ne_zero` | ★(c) の元は非零 |
-/

namespace ABC3.Found.Divisor

open CategoryTheory NumberField ABC3.Found.FrdI

open scoped NumberField

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

/-- ★`Example 6.3` の Frobenioid の対象 —— `(L, α)`(`α ∈ Φ(L)^gp`)。 -/
abbrev Ex63Obj (F Kbar : Type) [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
    [IsGalois F Kbar] : Type :=
  ModelData.Obj (ex63ModelData (F := F) (Kbar := Kbar))

/-- ★★★**`Example 6.3` の射の記述** —— 原文の (a)(b)(c) と 1 対 1。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model] -/
def ex63HomEquiv (A B : Ex63Obj F Kbar) :
    ModelData.Hom (ex63ModelData (F := F) (Kbar := Kbar)) A B
      ≃ {t : (A.base ⟶ B.base) × (ex63ModelData (F := F) (Kbar := Kbar)).phi.val A.base
              × ℕ+ × (ex63ModelData (F := F) (Kbar := Kbar)).bmon.val A.base //
          (t.2.2.1 : ℕ) • A.cls + toGpHom _ t.2.1
            = (ex63ModelData (F := F) (Kbar := Kbar)).phi.gpMapOn t.1 B.cls
              + (ex63ModelData (F := F) (Kbar := Kbar)).divB _ t.2.2.2} where
  toFun φ := ⟨(φ.base, φ.div, φ.deg, φ.u), φ.cond⟩
  invFun t := ⟨t.1.1, t.1.2.1, t.1.2.2.1, t.1.2.2.2, t.2⟩
  left_inv _ := rfl
  right_inv _ := rfl

/-- ★**(c) の元は `L` の非零元**。 -/
theorem ex63_hom_unit_ne_zero (A B : Ex63Obj F Kbar)
    (φ : ModelData.Hom (ex63ModelData (F := F) (Kbar := Kbar)) A B) :
    ((Additive.toMul (φ.u : Additive ((A.base.unop.toIF)ˣ)) : (A.base.unop.toIF)ˣ)
      : A.base.unop.toIF) ≠ 0 :=
  Units.ne_zero _

/-- ★**(b) の次数は `1` 以上**。 -/
theorem ex63_hom_deg_pos (A B : Ex63Obj F Kbar)
    (φ : ModelData.Hom (ex63ModelData (F := F) (Kbar := Kbar)) A B) :
    1 ≤ (φ.deg : ℕ) := φ.deg.property

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.3` の射の記述((a)(b)(c) の組)。 -/
def ex63HomEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 射は (Spec L → Spec M, d, 非零元) の組である",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
