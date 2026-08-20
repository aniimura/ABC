import ABC3.Found.GaloisRep.PsiDouble

/-!
# Galois (G1) 第 46 ブロック —— **★★★★★★同時帰納の述語と基底**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★乗法公式を `x` 側と `y` 側で**同時に**述べる

    MulFormula n :=
      (n•P = O → ΨSq_n(x) = 0)
      ∧ (n•P = (x', y') → x'·ΨSq_n(x) = Φ_n(x)
                        ∧ ψ(n•P)·ΨSq_n(x)² = preΨ_{2n}(x)·ψ(P))

★`ψ(Q) := y_Q − negY x_Q y_Q`(= `2y_Q + a₁x_Q + a₃`)。

★★**除算も二変数多項式も出てこない**——§9-386 で形を決めたとおりである。

## ★★★基底 3 つは在庫だけで出た

| `n` | `x` 側 | `y` 側 |
|---|---|---|
| 0 | ★`ΨSq_0 = 0`(mathlib) | 空虚 |
| 1 | ★`ΨSq_1 = 1`, `Φ_1 = X` | ★`preΨ_2 = 1` |
| 2 | ★★**第 30**(倍化公式) | ★★**第 45**(`ψ(2P)d³ = preΨ₄`) |

★★★★`n = 2` は `2P = O` かどうかで場合分けし、
どちらも第 29・31 の 2-捩れの判定で捌ける。

## ★★残るのは帰納段だけ

| 場合 | 在庫 |
|---|---|
| `nP = O` | ★第 42 |
| `(n+1)P = O` | ★第 43 |
| 一般 | ★第 44(分母)+ 加法公式(mathlib) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★同時帰納の述語。 -/
def MulFormula {x y : F} (h : W.toAffine.Nonsingular x y) (n : ℕ) : Prop :=
  ((n : ℕ) • (Point.some x y h) = 0 → (W.ΨSq (n : ℤ)).eval x = 0)
  ∧ (∀ x' y' h', (n : ℕ) • (Point.some x y h) = Point.some x' y' h' →
      x' * (W.ΨSq (n : ℤ)).eval x = (W.Φ (n : ℤ)).eval x
      ∧ (y' - W.toAffine.negY x' y') * ((W.ΨSq (n : ℤ)).eval x) ^ 2
        = (W.preΨ (2 * (n : ℤ))).eval x * (y - W.toAffine.negY x y))

theorem mulFormula_zero {x y : F} (h : W.toAffine.Nonsingular x y) :
    MulFormula W h 0 := by
  refine ⟨fun _ => ?_, fun x' y' h' hP => ?_⟩
  · norm_num [WeierstrassCurve.ΨSq_zero]
  · rw [zero_smul] at hP
    exact absurd hP.symm (Point.some_ne_zero h')

theorem mulFormula_one {x y : F} (h : W.toAffine.Nonsingular x y) :
    MulFormula W h 1 := by
  refine ⟨fun hP => ?_, fun x' y' h' hP => ?_⟩
  · rw [one_smul] at hP
    exact absurd hP (Point.some_ne_zero h)
  · rw [one_smul] at hP
    obtain ⟨rfl, rfl⟩ := Point.some.inj hP
    norm_num [WeierstrassCurve.ΨSq_one, WeierstrassCurve.Φ_one]

theorem preP_four : W.preΨ 4 = W.preΨ₄ := by
  have h := preNormEDS_four (W.Ψ₂Sq ^ 2) W.Ψ₃ W.preΨ₄
  exact h

theorem mulFormula_two {x y : F} (h : W.toAffine.Nonsingular x y) :
    MulFormula W h 2 := by
  by_cases hy : y = W.toAffine.negY x y
  · refine ⟨fun _ => ?_, fun x' y' h' hP => ?_⟩
    · rw [show ((2:ℕ):ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two]
      exact (psi2Sq_eval_eq_zero_iff W h.left).2 hy
    · exfalso
      have hz := (two_smul_eq_zero_iff W h).2 ((psi2Sq_eval_eq_zero_iff W h.left).2 hy)
      rw [hz] at hP
      exact Point.some_ne_zero h' hP.symm
  · have hQ : (2 : ℕ) • (Point.some x y h)
        = Point.some _ _ (nonsingular_add h h fun hxy => hy hxy.right) := by
      rw [two_nsmul]; exact Point.add_self_of_Y_ne hy
    refine ⟨fun hP => ?_, fun x' y' h' hP => ?_⟩
    · exfalso; rw [hQ] at hP; exact Point.some_ne_zero _ hP
    · rw [hQ] at hP
      obtain ⟨rfl, rfl⟩ := Point.some.inj hP
      rw [show ((2:ℕ):ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two,
        show (2:ℤ) * 2 = 4 by norm_num, preP_four W]
      refine ⟨doubling_x W h.left hy, ?_⟩
      rw [psi2Sq_eval W h.left]
      linear_combination (y - W.toAffine.negY x y) * psi_double W h.left hy


/-! ## ★出典の紐付け(`.src`) -/

def MulFormula.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の同時帰納——述語と基底)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
