/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex61Model

/-!
# `Example 6.1` の射の記述

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109–110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★原文が畳んだ所

`Example 6.1` の Frobenioid の射は
「(a) `V[M] → V[L]`、(b) `d ∈ ℕ≥1`、(c) 有理函数」の組である。
★これは `Theorem 5.2, (i)` の `ModelData.Hom` の 4 成分そのものであり、
`Ex63Morph.lean`(算術)と同じ形をしている。

| 原文 | `ModelData.Hom` の成分 | `Example 6.1` での中身 |
|---|---|---|
| (a) 底の射 | `base` | `𝒟 = B(G)⁰` の射 |
| (b) `d ∈ ℕ≥1` | `deg` | `ℕ+` |
| (c) 有理函数 | `u` | `B(L)`(零点・極が `D_L` にある有理函数) |
| — | `div` | 有効な Cartier 因子 `Φ(L)` の元 |
| — | `cond` | `d·α_A + Div(φ) = Base(φ)^*(α_B) + div(u)` |

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `ex61HomEquiv` | ★★★**射 ⟺ (底の射, 因子, 次数, 有理函数) の組で関係式を満たすもの** |
| `ex61_hom_unit_ne_zero` | ★(c) の有理函数は非零 |
| `ex61_hom_deg_pos` | ★(b) の次数は `1` 以上 |
-/

namespace ABC3.Found.Divisor

open CategoryTheory AlgebraicGeometry ABC3.Found.FrdI

universe u

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar] [IsGalois V.functionField Kbar]
  (DK : Set (PrimeDivisorPt V))
  [∀ L : FinSub V.functionField Kbar, IsLocallyNoetherian (normObj V L)]
  [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)]
  (hkq : IsKQCartier V DK
    (fun (L : FinSub V.functionField Kbar) _ => normObj_isNormalScheme V L))

/-- ★`Example 6.1` の Frobenioid の対象 —— `(L, α)`(`α ∈ Φ(L)^gp`)。 -/
abbrev Ex61Obj : Type u := ModelData.Obj (ex61ModelData V DK hkq)

/-- ★★★**`Example 6.1` の射の記述** —— 原文の (a)(b)(c) と 1 対 1。 -/
def ex61HomEquiv (A B : Ex61Obj V DK hkq) :
    ModelData.Hom (ex61ModelData V DK hkq) A B
      ≃ {t : (A.base ⟶ B.base) × (ex61ModelData V DK hkq).phi.val A.base
              × ℕ+ × (ex61ModelData V DK hkq).bmon.val A.base //
          (t.2.2.1 : ℕ) • A.cls + toGpHom _ t.2.1
            = (ex61ModelData V DK hkq).phi.gpMapOn t.1 B.cls
              + (ex61ModelData V DK hkq).divB _ t.2.2.2} where
  toFun φ := ⟨(φ.base, φ.div, φ.deg, φ.u), φ.cond⟩
  invFun t := ⟨t.1.1, t.1.2.1, t.1.2.2.1, t.1.2.2.2, t.2⟩
  left_inv _ := rfl
  right_inv _ := rfl

/-- ★**(c) の有理函数は非零**。 -/
theorem ex61_hom_unit_ne_zero (A B : Ex61Obj V DK hkq)
    (φ : ModelData.Hom (ex61ModelData V DK hkq) A B) :
    (((Additive.toMul (φ.u) : BSubgroup V DK A.base.unop (normObj_isNormalScheme V _))
      : ((normObj V A.base.unop).functionField)ˣ)
        : (normObj V A.base.unop).functionField) ≠ 0 :=
  Units.ne_zero _

/-- ★**(b) の次数は `1` 以上**。 -/
theorem ex61_hom_deg_pos (A B : Ex61Obj V DK hkq)
    (φ : ModelData.Hom (ex61ModelData V DK hkq) A B) :
    1 ≤ (φ.deg : ℕ) := φ.deg.property

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の射の記述((a)(b)(c) の組)。 -/
def ex61HomEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — 射は (V[M] → V[L], d, 有理函数) の組である",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
