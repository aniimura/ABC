/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottHyperplane
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★高さは射に沿って関手的である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★これは何か —— 一般の `X` へ渡る道

`§9-878` で `ℙᴺ_ℤ` の超平面因子について Northcott が出た。
★一般の `X`（`L_ℚ` が豊富）へ渡すには、射影埋め込み `ψ : X ⟶ ℙᴺ` に沿って

    `ht_{ψ^*E}(x) = ht_E(ψ ∘ x)`

が要る。★★本ファイルはそれを取る——**両辺とも定義を展開すれば同じ**である。

## ★★★機構 —— 有限側とアルキメデス側で 1 行ずつ

| 側 | 中身 |
|---|---|
| 有限 | `(E.divisor.comap ψ).comap xF = E.divisor.comap (xF ≫ ψ)`（`IdealSheafData.comap_comp`） |
| アルキメデス | `archPoint (xF ≫ ψ) v = archPoint xF v ≫ ψ`（結合則） |

★これで「`ℙᴺ` で測った高さ」と「`X` で測った高さ」が繋がる。

## ★残っている段（明示）

★★一般の `X` への還元には、あと

1. `L^n` が**非常に豊富**（`ψ^*O(1) ≅ L^n`）——段 E3
2. `ht_{L^n} = n·ht_L + O(1)`——`htArith_tensor_unconditional` の反復

が要る。★★★本ファイルはその 2 つを繋ぐ**土台**である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★複素点は射に沿って動く -/

/-- ★**`archPoint` は後合成と可換である**（結合則そのもの）。 -/
theorem archPoint_comp_right (F : Type) [Field F] [NumberField F] {X Y : Scheme.{0}}
    (xF : specRingOfIntegers F ⟶ X) (ψ : X ⟶ Y) (v : InfinitePlace F) :
    archPoint (xF ≫ ψ) v = archPoint xF v ≫ ψ := by
  rw [archPoint, archPoint, Category.assoc]

/-! ## ★★算術因子の引き戻し -/

/-- ★★**算術因子の引き戻し** —— 因子は `comap`、Green 関数は前合成。 -/
noncomputable def ArithCartier.comap {X Y : Scheme.{0}} (E : ArithCartier Y) (ψ : X ⟶ Y) :
    ArithCartier X :=
  { divisor := E.divisor.comap ψ
    green := fun p => E.green (p ≫ ψ) }

@[simp] theorem ArithCartier.comap_divisor {X Y : Scheme.{0}} (E : ArithCartier Y)
    (ψ : X ⟶ Y) : (E.comap ψ).divisor = E.divisor.comap ψ := rfl

@[simp] theorem ArithCartier.comap_green {X Y : Scheme.{0}} (E : ArithCartier Y)
    (ψ : X ⟶ Y) (p : complexPoints X) : (E.comap ψ).green p = E.green (p ≫ ψ) := rfl

/-- ★**有限側の関手性** —— `IdealSheafData.comap_comp` だけである。 -/
theorem pullbackIdeal_comap (F : Type) [Field F] [NumberField F] {X Y : Scheme.{0}}
    (E : ArithCartier Y) (ψ : X ⟶ Y) (xF : specRingOfIntegers F ⟶ X) :
    pullbackIdeal F (E.comap ψ).divisor xF = pullbackIdeal F E.divisor (xF ≫ ψ) := by
  show pullbackIdeal F (E.divisor.comap ψ) xF = _
  rw [← pullbackIdealOf_eq_pullbackIdeal, ← pullbackIdealOf_eq_pullbackIdeal,
    pullbackIdealOf, pullbackIdealOf, ← Scheme.IdealSheafData.comap_comp]

/-! ## ★★★★★★★★★★★★高さの関手性 -/

/-- ★★★★★★★★★★★★**高さは射に沿って関手的である**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `ht_{ψ^*E}(x) = ht_E(ψ ∘ x)`

★これで「`ℙᴺ` で測った高さ」と「`X` で測った高さ」が繋がる
——`§9-878`（`ℙᴺ` の Northcott）を一般の `X` へ渡す土台である。 -/
theorem htArith_comap (F : Type) [Field F] [NumberField F] {X Y : Scheme.{0}}
    (E : ArithCartier Y) (ψ : X ⟶ Y) (xF : specRingOfIntegers F ⟶ X) :
    htArith F (E.comap ψ) xF = htArith F E (xF ≫ ψ) := by
  rw [htArith_eq_add, htArith_eq_add, pullbackIdeal_comap, archADiv_sum_eq, archADiv_sum_eq]
  congr 2

/-- ★★**高さで有界な点の集合は引き戻しでも同じである**。 -/
theorem setOf_htArith_comap (F : Type) [Field F] [NumberField F] {X Y : Scheme.{0}}
    (E : ArithCartier Y) (ψ : X ⟶ Y) (C : ℝ) :
    {xF : specRingOfIntegers F ⟶ X | htArith F (E.comap ψ) xF ≤ C}
      = {xF | htArith F E (xF ≫ ψ) ≤ C} := by
  ext xF
  rw [Set.mem_setOf_eq, Set.mem_setOf_eq, htArith_comap]

/-! ## ★出典の紐付け(`.src`) -/

def archPoint_comp_right.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(archPoint は後合成と可換)",
    sectionId := "genell-prop-1-4" }

def ArithCartier.comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(算術因子の引き戻し)",
    sectionId := "genell-prop-1-4" }

def htArith_comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(高さは射に沿って関手的である)",
    sectionId := "genell-prop-1-4" }

def htArith_comap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.IdealSheafData.comap_comp"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.comap_comp") 2,
    .citation "[ABC3]" "htArith_eq_add(高さ = 有限側 + アルキメデス側)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_eq_add") 2,
    .implicitStep
      ("★機構は 2 行である: 有限側は IdealSheafData.comap_comp、" ++
       "アルキメデス側は archPoint (xF ≫ ψ) v = archPoint xF v ≫ ψ(結合則)") 2,
    .implicitStep
      ("★★一般の X への還元には、あと (1) L^n が**非常に豊富**(ψ^*O(1) ≅ L^n)——段 E3、" ++
       "(2) ht_{L^n} = n·ht_L + O(1)——htArith_tensor_unconditional の反復、が要る。" ++
       "★本ファイルはその 2 つを繋ぐ**土台**である") 4 ]

end ABC3.Found.GenEll
