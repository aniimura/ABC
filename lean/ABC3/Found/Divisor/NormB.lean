/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormKQC

/-!
# `B(L)` —— 零点・極が `D_L` にある有理函数(鎖 `normalize` の `B-functor` の第 1 段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★`ord` の 3 本だけで部分群になる

`B(L)` は「`V[L]` 上の有理函数で、零点と極が `D_L` に台を持つもの」である。
★スキームの言葉では **`D_L` の外のすべての素因子で `ord` が消えること**:

  `B(L) = {u ∈ K(V[L])^× | ∀ x ∉ D_L, ord_x(u) = 0}`

★★部分群であることの中身は `SchemeWeilOrd.lean` / `SchemeCartier.lean` の
**`ordPt_one` / `ordPt_mul` / `ordPt_inv`** の 3 本だけである。
★★★**`div(f)` の台の有限性は要らない** —— 台の条件を「外で消える」と書けば
有限性を経由しないで済む。

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `BSubgroup` | ★★**`B(L)`**(零点・極が `D_L` にある有理函数の群) |
| `mem_BSubgroup` | 判定条件 |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

variable (V : Scheme.{u}) [IsIntegral V] [IsLocallyNoetherian V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-- ★★★**`B(L)`** —— `V[L]` 上の有理函数で、零点と極が `D_L` に台を持つもの。

★部分群であることの中身は `ordPt_one` / `ordPt_mul` / `ordPt_inv` の 3 本だけ。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported -/
noncomputable def BSubgroup (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar)
    [IsLocallyNoetherian (normObj V L)] (hnorm : IsNormalScheme (normObj V L)) :
    Subgroup ((normObj V L).functionField)ˣ where
  carrier := {u | ∀ x : PrimeDivisorPt (normObj V L), x ∉ DLSet V DK L →
    ordPt (normObj V L) hnorm x (u : (normObj V L).functionField) = 0}
  one_mem' := by
    intro x _
    show ordPt (normObj V L) hnorm x (1 : (normObj V L).functionField) = 0
    exact ordPt_one hnorm x
  mul_mem' := by
    intro u v hu hv x hx
    show ordPt (normObj V L) hnorm x
      ((u : (normObj V L).functionField) * (v : (normObj V L).functionField)) = 0
    rw [ordPt_mul hnorm x u.ne_zero v.ne_zero, hu x hx, hv x hx, add_zero]
  inv_mem' := by
    intro u hu x hx
    show ordPt (normObj V L) hnorm x ((u⁻¹ : ((normObj V L).functionField)ˣ) :
      (normObj V L).functionField) = 0
    rw [Units.val_inv_eq_inv_val, ordPt_inv hnorm x u.ne_zero, hu x hx, neg_zero]

omit [IsLocallyNoetherian V] in
@[simp] theorem mem_BSubgroup (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar)
    [IsLocallyNoetherian (normObj V L)] (hnorm : IsNormalScheme (normObj V L))
    (u : ((normObj V L).functionField)ˣ) :
    u ∈ BSubgroup V DK L hnorm ↔ ∀ x : PrimeDivisorPt (normObj V L), x ∉ DLSet V DK L →
      ordPt (normObj V L) hnorm x (u : (normObj V L).functionField) = 0 := Iff.rfl

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `B(L)`。 -/
def BSubgroup.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — B(L)(零点・極が D_L にある有理函数の群)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
