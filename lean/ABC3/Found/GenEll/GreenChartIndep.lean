/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GreenFubiniStudy
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★Fubini–Study はチャートに依らない（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★これは何か —— 大域的な `green` を作るための本体

`§9-868` の `greenChartOf` は**チャートを固定した**形である。
★`ArithCartier` の欄（`green : complexPoints X → ℝ`）に入れるには、
**別のチャートで測っても同じ値になる**ことが要る。

★★本ファイルはその**代数の側**を取る:

| 補題 | 内容 |
|---|---|
| `toOverlap` / `toOverlap'` | `A⁰_{x_i}`・`A⁰_{x_{i'}}` から重なり `A⁰_{x_i x_{i'}}` へ |
| `toOverlap_eq_mul` | ★★**`x_k/x_i = (x_k/x_{i'})·(x_{i'}/x_i)`**（重なりの上で） |
| `log_iSup_norm_ratio_congr` | ★共通因子は `sup` と比の両方から消える |
| `greenChartOf_overlap` | ★★★**Fubini–Study はチャートに依らない** |

## ★★★機構 —— 遷移単元が約分される

★重なり `A⁰_{x_i x_{i'}}` の中では

    `x_k x_{i'}/(x_i x_{i'}) = (x_k x_i/(x_i x_{i'})) · (x_{i'}²/(x_i x_{i'}))`

が成り立つ（両辺とも `x_k/x_i`）。★★右の第 2 因子 `t = x_{i'}/x_i` は **`k` に依らない**ので、
`sup_k` からも比の分母からも同じだけ出て**約分される**。
★★★これは `§9-868` の `log_iSup_norm_eq`（`v(x_{i₀})` が約分される）と**同じ形**である。

## ★残っている段（明示）

★★★★大域的な `green` を書き下すには、あと 2 つ要る:

1. 「同じ複素点が 2 つのチャートに入るなら、**重なりのチャートを通る**」という
   スキームの側の段（`Proj.SpecMap_awayMap_awayι` と開埋め込みの持ち上げ）。
2. 得られた `green` の**連続性**（`htArith_sub_abs_le` が要求する `ArcModel` の位相）。
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★重なりのチャート -/

/-- ★**`x_i x_{i'}` は 2 次の斉次式である**。 -/
theorem X_mul_X_mem (N : ℕ) (i i' : Fin (N+1)) :
    (MvPolynomial.X i * MvPolynomial.X i' : MvPolynomial (Fin (N+1)) ℤ)
      ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ 2 := by
  have h := (MvPolynomial.isHomogeneous_X ℤ i).mul (MvPolynomial.isHomogeneous_X ℤ i')
  simpa using h

/-- ★**`A⁰_{x_i} → A⁰_{x_i x_{i'}}`**。 -/
noncomputable abbrev toOverlap (N : ℕ) (i i' : Fin (N+1)) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i)
      →+* HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i * MvPolynomial.X i') :=
  awayMap _ (MvPolynomial.isHomogeneous_X ℤ i') rfl

/-- ★**`A⁰_{x_{i'}} → A⁰_{x_i x_{i'}}`**。 -/
noncomputable abbrev toOverlap' (N : ℕ) (i i' : Fin (N+1)) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i')
      →+* HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i * MvPolynomial.X i') :=
  awayMap _ (MvPolynomial.isHomogeneous_X ℤ i) (mul_comm _ _)

/-- ★★**遷移単元** `x_{i'}²/(x_i x_{i'}) = x_{i'}/x_i`。 -/
noncomputable abbrev overlapUnit (N : ℕ) (i i' : Fin (N+1)) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
      (MvPolynomial.X i * MvPolynomial.X i') :=
  Away.mk _ (X_mul_X_mem N i i') 1 (MvPolynomial.X i' * MvPolynomial.X i')
    (X_mul_X_mem N i' i')

/-- ★★★★★**重なりの上で `x_k/x_i = (x_k/x_{i'})·(x_{i'}/x_i)`**。

★これがチャートの遷移そのものである。★★右の第 2 因子は **`k` に依らない**。 -/
theorem toOverlap_eq_mul (N : ℕ) (i i' k : Fin (N+1)) :
    toOverlap N i i' (projCoord N ℤ i k)
      = toOverlap' N i i' (projCoord N ℤ i' k) * overlapUnit N i i' := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [projCoord, projCoord, toOverlap, toOverlap', awayMap_mk, awayMap_mk,
    HomogeneousLocalization.val_mul, overlapUnit, Away.val_mk, Away.val_mk, Away.val_mk,
    Localization.mk_mul, Localization.mk_eq_mk_iff, Localization.r_iff_exists]
  exact ⟨1, by push_cast; ring⟩

/-! ## ★★共通因子は約分される -/

/-- ★★**共通因子は `sup` と比の両方から消える**。

★`§9-868` の `log_iSup_norm_eq`（`v(x_{i₀})` が約分される）と**同じ形**である。 -/
theorem log_iSup_norm_ratio_congr {ι : Type} (c d : ι → ℂ) (t : ℂ) (j : ι)
    (hc : ∀ k, c k = d k * t) (ht : t ≠ 0) :
    Real.log ((⨆ k, ‖c k‖) / ‖c j‖) = Real.log ((⨆ k, ‖d k‖) / ‖d j‖) := by
  have hnt : (0:ℝ) < ‖t‖ := norm_pos_iff.2 ht
  have hck : ∀ k, ‖c k‖ = ‖d k‖ * ‖t‖ := fun k => by rw [hc k, norm_mul]
  simp only [hck]
  rw [← Real.iSup_mul_of_nonneg hnt.le (fun k => ‖d k‖),
    mul_div_mul_right _ _ hnt.ne']

/-! ## ★★★★★★★★★★★Fubini–Study はチャートに依らない -/

/-- ★★★★★★★★★★★**Fubini–Study はチャートに依らない**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★重なりのチャートを通る複素点について、`i` で測っても `i'` で測っても同じ値になる。
★★機構は `toOverlap_eq_mul`（遷移単元が `k` に依らない）と
`log_iSup_norm_ratio_congr`（共通因子は約分される）だけである。 -/
theorem greenChartOf_overlap (N : ℕ) (i i' : Fin (N+1))
    (χ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i * MvPolynomial.X i')) ⟶ CommRingCat.of ℂ)
    (hne : χ.hom (overlapUnit N i i') ≠ 0) :
    greenChartOf N i (CommRingCat.ofHom (toOverlap N i i') ≫ χ)
      = greenChartOf N i' (CommRingCat.ofHom (toOverlap' N i i') ≫ χ) := by
  have hc : ∀ a, (CommRingCat.ofHom (toOverlap N i i') ≫ χ).hom a
      = χ.hom (toOverlap N i i' a) := by
    intro a; rw [CommRingCat.hom_comp, CommRingCat.hom_ofHom]; rfl
  have hc' : ∀ a, (CommRingCat.ofHom (toOverlap' N i i') ≫ χ).hom a
      = χ.hom (toOverlap' N i i' a) := by
    intro a; rw [CommRingCat.hom_comp, CommRingCat.hom_ofHom]; rfl
  rw [greenChartOf, greenChartOf]
  simp only [hc, hc']
  refine log_iSup_norm_ratio_congr
    (fun k => χ.hom (toOverlap N i i' (projCoord N ℤ i k)))
    (fun k => χ.hom (toOverlap' N i i' (projCoord N ℤ i' k)))
    (χ.hom (overlapUnit N i i')) 0 (fun k => ?_) hne
  rw [toOverlap_eq_mul, map_mul]

/-! ## ★出典の紐付け(`.src`) -/

def toOverlap_eq_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(重なりの上で x_k/x_i = (x_k/x_{i'})·(x_{i'}/x_i))",
    sectionId := "genell-prop-1-4" }

def log_iSup_norm_ratio_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(共通因子は sup と比の両方から消える)",
    sectionId := "genell-prop-1-4" }

def greenChartOf_overlap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(Fubini–Study はチャートに依らない)",
    sectionId := "genell-prop-1-4" }

def greenChartOf_overlap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "greenChartOf(チャートを固定した Fubini–Study、§9-868)"
      (.inProject "ABC3" "ABC3.Found.GenEll.greenChartOf") 2,
    .citation "[mathlib]" "HomogeneousLocalization.awayMap_mk(重なりへの写像)"
      (.inMathlib "HomogeneousLocalization.awayMap_mk") 2,
    .implicitStep
      ("★重なり A⁰_{x_i x_{i'}} の中では " ++
       "x_k x_{i'}/(x_i x_{i'}) = (x_k x_i/(x_i x_{i'}))·(x_{i'}²/(x_i x_{i'})) が成り立つ" ++
       "(両辺とも x_k/x_i)。★★右の第 2 因子 t = x_{i'}/x_i は **k に依らない**ので、" ++
       "sup_k からも比の分母からも同じだけ出て**約分される**") 3,
    .implicitStep
      ("★★★大域的な green を書き下すには、あと 2 つ要る: " ++
       "(1)「同じ複素点が 2 つのチャートに入るなら**重なりのチャートを通る**」という" ++
       "スキームの側の段(Proj.SpecMap_awayMap_awayι と開埋め込みの持ち上げ); " ++
       "(2) 得られた green の**連続性**(htArith_sub_abs_le が要求する ArcModel の位相)") 5 ]

end ABC3.Found.GenEll
