/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GreenGlobal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★座標の組から `ℙᴺ` の複素点を作る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★これは何か —— `ArcModel` の逆向き

`§9-871` までで、複素点から**座標**を取る向き（`projPointCoord`）と、
そこから Fubini–Study を測る向き（`greenFS`）ができた。
★`ArcModel (ℙᴺ_ℤ) (ℂ^{N+1})` を作るには**逆向き**——
座標の組から複素点を作る向き——が要る（`emb_range` の証明に効く）。

★★本ファイルはそれを取る:

| 定義・補題 | 内容 |
|---|---|
| `evalCoords` | `ℤ[x] →+* ℂ`、`x_k ↦ a_k/a_i` |
| `awayHomOfCoords` | ★**`A⁰_{x_i} →+* ℂ`**（`x_i − 1` を殺すので降りる） |
| `awayHomOfCoords_projCoord` | `x_k/x_i ↦ a_k/a_i` |
| `projPointOfCoords` | ★★座標の組が定める ℂ-点 |
| `greenFS_projPointOfCoords` | ★★★**`greenFS` の同次座標での表示** |

## ★★★機構 —— `A⁰_{x_i} ≅ ℤ[x]/(x_i − 1)` を使う

★`§9-859`（段 C2a-2）の `awayQuotHom : A⁰_{x_i} →+* ℤ[x]/(x_i − 1)` に
`Ideal.Quotient.lift` を繋ぐだけである。★★`x_i − 1` が核に入るのは
`a_i/a_i = 1` だからで、ここで **`a_i ≠ 0`** を使う。

## ★★★★★到達点 —— 同次座標での Fubini–Study

    `greenFS N (projPointOfCoords N a i hi) = log( sup_k ‖a_k‖ / ‖a_0‖ )`

★これが `ArcModel` の位相（`ℙ(ℂ^{N+1})` からの引き戻し）で連続性を論じるときに読む形である
——右辺は `a` の**同次座標だけ**で書けている。
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★座標の組が定める環準同型 -/

/-- ★**座標の組から `ℤ[x] →+* ℂ`**（`x_k ↦ a_k/a_i`）。 -/
noncomputable def evalCoords (N : ℕ) (a : Fin (N+1) → ℂ) (i : Fin (N+1)) :
    MvPolynomial (Fin (N+1)) ℤ →+* ℂ :=
  MvPolynomial.eval₂Hom (Int.castRingHom ℂ) (fun k => a k / a i)

theorem evalCoords_X (N : ℕ) (a : Fin (N+1) → ℂ) (i k : Fin (N+1)) :
    evalCoords N a i (MvPolynomial.X k) = a k / a i := by
  rw [evalCoords, MvPolynomial.eval₂Hom_X']

/-- ★**`x_i − 1` は核に入る** —— ここで `a_i ≠ 0` を使う。 -/
theorem span_le_ker_evalCoords (N : ℕ) (a : Fin (N+1) → ℂ) (i : Fin (N+1)) (hi : a i ≠ 0) :
    Ideal.span {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) ℤ)}
      ≤ RingHom.ker (evalCoords N a i) := by
  rw [Ideal.span_le]
  rintro z rfl
  rw [SetLike.mem_coe, RingHom.mem_ker, map_sub, map_one, evalCoords_X, div_self hi, sub_self]

/-- ★★★★★**座標の組が定める `A⁰_{x_i} →+* ℂ`**。

★`§9-859`（段 C2a-2）の `A⁰_{x_i} ≅ ℤ[x]/(x_i − 1)` に `Ideal.Quotient.lift` を繋ぐだけ。 -/
noncomputable def awayHomOfCoords (N : ℕ) (a : Fin (N+1) → ℂ) (i : Fin (N+1)) (hi : a i ≠ 0) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
      (MvPolynomial.X i) →+* ℂ :=
  (Ideal.Quotient.lift _ (evalCoords N a i) (fun _ hz =>
    RingHom.mem_ker.1 (span_le_ker_evalCoords N a i hi hz))).comp (awayQuotHom N ℤ i)

/-- ★★**正規化座標は `a_k/a_i` へ行く**。 -/
theorem awayHomOfCoords_projCoord (N : ℕ) (a : Fin (N+1) → ℂ) (i : Fin (N+1)) (hi : a i ≠ 0)
    (k : Fin (N+1)) :
    awayHomOfCoords N a i hi (projCoord N ℤ i k) = a k / a i := by
  have hq : awayQuotHom N ℤ i (projCoord N ℤ i k)
      = Ideal.Quotient.mk _ (MvPolynomial.X k) := by
    have h := congrArg (fun f : MvPolynomial (Fin (N+1)) ℤ →+*
        (MvPolynomial (Fin (N+1)) ℤ ⧸ Ideal.span
          {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) ℤ)}) => f (MvPolynomial.X k))
      (awayQuotHom_comp_awayEval N ℤ i)
    simp only [RingHom.comp_apply] at h
    rw [show awayEval N ℤ i (MvPolynomial.X k) = projCoord N ℤ i k from by
      rw [awayEval, MvPolynomial.eval₂Hom_X']] at h
    exact h
  rw [awayHomOfCoords, RingHom.comp_apply, hq, Ideal.Quotient.lift_mk, evalCoords_X]

/-! ## ★★★★★★★★★★★座標の組が定める複素点 -/

/-- ★★**座標の組が定める ℂ-点**。 -/
noncomputable def projPointOfCoords (N : ℕ) (a : Fin (N+1) → ℂ) (i : Fin (N+1)) (hi : a i ≠ 0) :
    complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :=
  Spec.map (CommRingCat.ofHom (awayHomOfCoords N a i hi)) ≫ chartA N ℤ i

/-- ★★★★★★★★★★★**`greenFS` の同次座標での表示**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `greenFS N (projPointOfCoords N a i hi) = log( sup_k ‖a_k‖ / ‖a_0‖ )`

★右辺は `a` の**同次座標だけ**で書けている——`i` にも `hi` にも依らない。
★★これが `ArcModel` の位相（`ℙ(ℂ^{N+1})` からの引き戻し）で
連続性を論じるときに読む形である。 -/
theorem greenFS_projPointOfCoords (N : ℕ) (a : Fin (N+1) → ℂ) (i : Fin (N+1)) (hi : a i ≠ 0) :
    greenFS N (projPointOfCoords N a i hi)
      = Real.log ((⨆ k, ‖a k‖) / ‖a 0‖) := by
  rw [projPointOfCoords, greenFS_eq_greenChartOf, greenChartOf]
  have hval : ∀ k, (CommRingCat.ofHom (awayHomOfCoords N a i hi)).hom (projCoord N ℤ i k)
      = a k * (a i)⁻¹ := by
    intro k
    rw [CommRingCat.hom_ofHom, awayHomOfCoords_projCoord, div_eq_mul_inv]
  simp only [hval]
  exact log_iSup_norm_ratio_congr (fun k => a k * (a i)⁻¹) a (a i)⁻¹ 0
    (fun _ => rfl) (inv_ne_zero hi)

/-! ## ★出典の紐付け(`.src`) -/

def awayHomOfCoords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(座標の組が定める A⁰_{x_i} →+* ℂ)",
    sectionId := "genell-prop-1-4" }

def projPointOfCoords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(座標の組が定める ℂ-点)",
    sectionId := "genell-prop-1-4" }

def greenFS_projPointOfCoords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(greenFS の同次座標での表示)",
    sectionId := "genell-prop-1-4" }

def greenFS_projPointOfCoords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "awayQuotHom / ker_awayEval(A⁰_{x_i} ≅ ℤ[x]/(x_i − 1)、段 C2a-2、§9-859)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayQuotHom") 2,
    .citation "[ABC3]" "greenFS_eq_greenChartOf(チャートに依らない、§9-871)"
      (.inProject "ABC3" "ABC3.Found.GenEll.greenFS_eq_greenChartOf") 2,
    .citation "[ABC3]" "log_iSup_norm_ratio_congr(共通因子は約分される、§9-869)"
      (.inProject "ABC3" "ABC3.Found.GenEll.log_iSup_norm_ratio_congr") 2,
    .implicitStep
      ("★到達点の右辺 log( sup_k ‖a_k‖ / ‖a_0‖ ) は a の**同次座標だけ**で書けている" ++
       "——i にも hi にも依らない。これが ArcModel の位相" ++
       "(ℙ(ℂ^{N+1}) からの引き戻し)で連続性を論じるときに読む形である") 3,
    .implicitStep
      ("★★段 C2c-3 に残るのは ArcModel (ℙᴺ_ℤ) (ℂ^{N+1}) の構成である: " ++
       "emb は projPointCoord、emb_range の全射性は本ファイルの projPointOfCoords、" ++
       "cone は univ(閉錐)。★emb_injective は §9-850 の ext_of_projCoord から出る") 4 ]

end ABC3.Found.GenEll
