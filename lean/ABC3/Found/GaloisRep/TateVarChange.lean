import ABC3.Found.GaloisRep.TateDvrSetup
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point

/-!
# Galois (G6) 第 303 ブロック —— **★★★★★★★★変数変換は点の加法同型を与える**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> `W.Point ≃+ (C • W).Point`(`vcPointAddEquiv`)

★★★**mathlib に 0 件**である(2026-08-26 実測)。係数側の
`variableChange_c₄`・`variableChange_Δ`・`variableChange_j` は在るが、
**点への作用が無い**。★G6 の残り 4 段のうちの第 1 段である。

## ★★★★★★変数変換は `(x, y)` を動かす

mathlib の `C • W` は `x = u²x' + r`、`y = u³y' + u²sx' + t` の置き換えである。
したがって点は**逆向き**に動く:

    x' = u⁻²(x − r)                        (`vcX`)
    y' = u⁻³(y − s(x − r) − t)              (`vcY`)

## ★★★★★★★重さで揃う——4 本の恒等式

| 対象 | 重み | 恒等式 |
|---|---|---|
| 方程式 `F` | `u⁻⁶` | `F'(x', y') = u⁻⁶ F(x, y)`(`equation_vc`) |
| `∂F/∂x` | `u⁻⁴` | `∂F'/∂x' = u⁻⁴(∂F/∂x + s ∂F/∂y)`(`nonsingular_vc`) |
| `∂F/∂y` | `u⁻³` | `∂F'/∂y' = u⁻³ ∂F/∂y` |
| 傾き `ℓ` | `u⁻¹` | `ℓ' = u⁻¹(ℓ − s)`(`slope_vc_of_*`) |

★★★★**`x` 微分だけが `s` で混ざる**——置き換えの Jacobi 行列が三角だからである。
そのため**単独の非退化条件は保たれない**が、
**2 つ一組(勾配が消えないこと)は保たれる**。ここが唯一の注意点だった。

## ★★★★★傾きの退化だけは別扱い

mathlib は `y₁ = negY x₂ y₂` のとき `slope := 0` と**約束**している。
★このとき `ℓ' = u⁻¹(ℓ − s)` は**成り立たない**(`0 ≠ −u⁻¹s`)。
★★しかしその場合は和が `0` になる場合(`add_of_Y_eq`)なので、傾きを使わない。
★★★したがって傾きの恒等式は **`x₁ ≠ x₂` と「`x₁ = x₂` かつ非退化」の 2 つに分ける**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `equation_vc`・`nonsingular_vc` | ★★★★★★方程式と非特異性の移り方 |
| `negY_vc`・`addX_vc`・`addY_vc` | ★★★★★加法公式の移り方 |
| `slope_vc_of_X_ne`・`slope_vc_of_X_eq` | ★★★★★★傾きの移り方 |
| `vcPoint_add` | ★★★★★★★**加法を保つ** |
| `vcPointAddEquiv` | ★★★★★★★★**`W.Point ≃+ (C • W).Point`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

section VarChange

variable {F : Type} [Field F]

/-! ## ★★★★★動いた座標 -/

/-- ★★変数変換後の `x` 座標 `u⁻²(x − r)`。 -/
def vcX (C : VariableChange F) (x : F) : F := ((C.u⁻¹ : Fˣ) : F) ^ 2 * (x - C.r)

/-- ★★変数変換後の `y` 座標 `u⁻³(y − s(x − r) − t)`。 -/
def vcY (C : VariableChange F) (x y : F) : F :=
  ((C.u⁻¹ : Fˣ) : F) ^ 3 * (y - C.s * (x - C.r) - C.t)

theorem vcX_inj (C : VariableChange F) {x x' : F} (h : vcX C x = vcX C x') : x = x' := by
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  simp only [vcX] at h
  have h2 := mul_left_cancel₀ (pow_ne_zero 2 hu0) h
  linear_combination h2

theorem vcY_inj (C : VariableChange F) {x y y' : F} (h : vcY C x y = vcY C x y') : y = y' := by
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  simp only [vcY] at h
  have h2 := mul_left_cancel₀ (pow_ne_zero 3 hu0) h
  linear_combination h2

/-- ★逆向きの `x` 座標を入れると戻る。 -/
theorem vcX_apply_inv (C : VariableChange F) (x' : F) :
    vcX C (((C.u : Fˣ) : F) ^ 2 * x' + C.r) = x' := by
  have h : ((C.u⁻¹ : Fˣ) : F) * ((C.u : Fˣ) : F) = 1 := by
    rw [← Units.val_mul]
    simp
  have h2 : ((C.u⁻¹ : Fˣ) : F) ^ 2 * ((C.u : Fˣ) : F) ^ 2 = 1 := by
    rw [← mul_pow, h, one_pow]
  simp only [vcX]
  linear_combination x' * h2

/-- ★逆向きの `y` 座標を入れると戻る。 -/
theorem vcY_apply_inv (C : VariableChange F) (x' y' : F) :
    vcY C (((C.u : Fˣ) : F) ^ 2 * x' + C.r)
        (((C.u : Fˣ) : F) ^ 3 * y' + ((C.u : Fˣ) : F) ^ 2 * C.s * x' + C.t) = y' := by
  have h : ((C.u⁻¹ : Fˣ) : F) * ((C.u : Fˣ) : F) = 1 := by
    rw [← Units.val_mul]
    simp
  have h3 : ((C.u⁻¹ : Fˣ) : F) ^ 3 * ((C.u : Fˣ) : F) ^ 3 = 1 := by
    rw [← mul_pow, h, one_pow]
  simp only [vcY]
  linear_combination y' * h3

/-! ## ★★★★★★方程式と非特異性 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**方程式は `u⁻⁶` 倍で移る**。 -/
theorem equation_vc (W : WeierstrassCurve F) (C : VariableChange F) (x y : F) :
    (C • W).toAffine.Equation (vcX C x) (vcY C x y) ↔ W.toAffine.Equation x y := by
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hE : ((vcY C x y) ^ 2 + (C • W).a₁ * (vcX C x) * (vcY C x y) + (C • W).a₃ * (vcY C x y))
        - ((vcX C x) ^ 3 + (C • W).a₂ * (vcX C x) ^ 2 + (C • W).a₄ * (vcX C x) + (C • W).a₆)
      = ((C.u⁻¹ : Fˣ) : F) ^ 6
        * ((y ^ 2 + W.a₁ * x * y + W.a₃ * y)
          - (x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆)) := by
    simp only [vcX, vcY, variableChange_a₁, variableChange_a₂, variableChange_a₃,
      variableChange_a₄, variableChange_a₆]
    ring
  rw [WeierstrassCurve.Affine.equation_iff, WeierstrassCurve.Affine.equation_iff,
    ← sub_eq_zero (a := (vcY C x y) ^ 2 + _ * _ * _ + _ * _),
    ← sub_eq_zero (a := y ^ 2 + W.a₁ * x * y + W.a₃ * y), hE, mul_eq_zero]
  constructor
  · rintro (h | h)
    · exact absurd h (pow_ne_zero _ hu0)
    · exact h
  · intro h
    exact Or.inr h

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**非特異性は保たれる**——`x` 微分は `s` で混ざるが、勾配が消えないことは同値。 -/
theorem nonsingular_vc (W : WeierstrassCurve F) (C : VariableChange F) (x y : F) :
    (C • W).toAffine.Nonsingular (vcX C x) (vcY C x y) ↔ W.toAffine.Nonsingular x y := by
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have h1 : (C • W).a₁ * vcY C x y
        - (3 * (vcX C x) ^ 2 + 2 * (C • W).a₂ * (vcX C x) + (C • W).a₄)
      = ((C.u⁻¹ : Fˣ) : F) ^ 4
        * ((W.a₁ * y - (3 * x ^ 2 + 2 * W.a₂ * x + W.a₄))
          + C.s * (y - (-y - W.a₁ * x - W.a₃))) := by
    simp only [vcX, vcY, variableChange_a₁, variableChange_a₂, variableChange_a₄]
    ring
  have h2 : vcY C x y - (-(vcY C x y) - (C • W).a₁ * (vcX C x) - (C • W).a₃)
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * (y - (-y - W.a₁ * x - W.a₃)) := by
    simp only [vcX, vcY, variableChange_a₁, variableChange_a₃]
    ring
  have hAB : (((C • W).a₁ * vcY C x y
          - (3 * (vcX C x) ^ 2 + 2 * (C • W).a₂ * (vcX C x) + (C • W).a₄) = 0)
        ∧ (vcY C x y - (-(vcY C x y) - (C • W).a₁ * (vcX C x) - (C • W).a₃) = 0))
      ↔ ((W.a₁ * y - (3 * x ^ 2 + 2 * W.a₂ * x + W.a₄) = 0)
        ∧ (y - (-y - W.a₁ * x - W.a₃) = 0)) := by
    constructor
    · rintro ⟨hA, hB⟩
      rw [h2, mul_eq_zero] at hB
      have hB' : y - (-y - W.a₁ * x - W.a₃) = 0 := hB.resolve_left (pow_ne_zero _ hu0)
      rw [h1, hB', mul_eq_zero] at hA
      have hA' := hA.resolve_left (pow_ne_zero _ hu0)
      exact ⟨by linear_combination hA', hB'⟩
    · rintro ⟨hA, hB⟩
      refine ⟨?_, ?_⟩
      · rw [h1, hA, hB]
        ring
      · rw [h2, hB]
        ring
  rw [WeierstrassCurve.Affine.nonsingular_iff, WeierstrassCurve.Affine.nonsingular_iff,
    equation_vc]
  refine and_congr_right fun _ => ?_
  rw [ne_eq, ne_eq, ne_eq, ne_eq,
    ← sub_eq_zero (a := (C • W).a₁ * vcY C x y),
    ← sub_eq_zero (a := vcY C x y),
    ← sub_eq_zero (a := W.a₁ * y),
    ← sub_eq_zero (a := y),
    ← not_and_or, ← not_and_or]
  exact not_congr hAB

/-! ## ★★★★★加法公式の移り方 -/

/-- ★★★★★`negY` は座標の移動と可換。 -/
theorem negY_vc (W : WeierstrassCurve F) (C : VariableChange F) (x y : F) :
    (C • W).toAffine.negY (vcX C x) (vcY C x y) = vcY C x (W.toAffine.negY x y) := by
  show -(vcY C x y) - (C • W).a₁ * (vcX C x) - (C • W).a₃ = _
  simp only [vcX, vcY, variableChange_a₁, variableChange_a₃, WeierstrassCurve.Affine.negY]
  ring

/-- ★★★★★`addX` は座標の移動と可換(傾きは `u⁻¹(ℓ − s)`)。 -/
theorem addX_vc (W : WeierstrassCurve F) (C : VariableChange F) (x₁ x₂ L : F) :
    (C • W).toAffine.addX (vcX C x₁) (vcX C x₂) (((C.u⁻¹ : Fˣ) : F) * (L - C.s))
      = vcX C (W.toAffine.addX x₁ x₂ L) := by
  simp only [vcX, WeierstrassCurve.Affine.addX, variableChange_a₁, variableChange_a₂]
  ring

/-- ★★★★★`addY` は座標の移動と可換。 -/
theorem addY_vc (W : WeierstrassCurve F) (C : VariableChange F) (x₁ x₂ y₁ L : F) :
    (C • W).toAffine.addY (vcX C x₁) (vcX C x₂) (vcY C x₁ y₁) (((C.u⁻¹ : Fˣ) : F) * (L - C.s))
      = vcY C (W.toAffine.addX x₁ x₂ L) (W.toAffine.addY x₁ x₂ y₁ L) := by
  simp only [vcX, vcY, WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negAddY,
    WeierstrassCurve.Affine.addX, WeierstrassCurve.Affine.negY,
    variableChange_a₁, variableChange_a₂, variableChange_a₃]
  ring

/-! ## ★★★★★★傾きの移り方 -/

variable [DecidableEq F]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**傾きは `u⁻¹(ℓ − s)` で移る**(`x` が異なる場合)。 -/
theorem slope_vc_of_X_ne (W : WeierstrassCurve F) (C : VariableChange F) {x₁ x₂ : F} (y₁ y₂ : F)
    (hx : x₁ ≠ x₂) :
    (C • W).toAffine.slope (vcX C x₁) (vcX C x₂) (vcY C x₁ y₁) (vcY C x₂ y₂)
      = ((C.u⁻¹ : Fˣ) : F) * (W.toAffine.slope x₁ x₂ y₁ y₂ - C.s) := by
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hx' : vcX C x₁ ≠ vcX C x₂ := fun h => hx (vcX_inj C h)
  rw [WeierstrassCurve.Affine.slope_of_X_ne hx', WeierstrassCurve.Affine.slope_of_X_ne hx]
  have hne : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hx
  have hnum : vcY C x₁ y₁ - vcY C x₂ y₂
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * ((y₁ - y₂) - C.s * (x₁ - x₂)) := by
    simp only [vcY]
    ring
  have hden : vcX C x₁ - vcX C x₂ = ((C.u⁻¹ : Fˣ) : F) ^ 2 * (x₁ - x₂) := by
    simp only [vcX]
    ring
  rw [hnum, hden]
  field_simp

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**傾きは `u⁻¹(ℓ − s)` で移る**(`x` が等しく、非退化な場合)。 -/
theorem slope_vc_of_X_eq (W : WeierstrassCurve F) (C : VariableChange F) {x₁ x₂ y₁ y₂ : F}
    (h₁ : W.toAffine.Equation x₁ y₁) (h₂ : W.toAffine.Equation x₂ y₂) (hx : x₁ = x₂)
    (hy : y₁ ≠ W.toAffine.negY x₂ y₂) :
    (C • W).toAffine.slope (vcX C x₁) (vcX C x₂) (vcY C x₁ y₁) (vcY C x₂ y₂)
      = ((C.u⁻¹ : Fˣ) : F) * (W.toAffine.slope x₁ x₂ y₁ y₂ - C.s) := by
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hyy : y₁ = y₂ := WeierstrassCurve.Affine.Y_eq_of_Y_ne h₁ h₂ hx hy
  subst hx
  subst hyy
  have hD : y₁ - W.toAffine.negY x₁ y₁ ≠ 0 := sub_ne_zero.2 hy
  have hy' : vcY C x₁ y₁ ≠ (C • W).toAffine.negY (vcX C x₁) (vcY C x₁ y₁) := by
    rw [negY_vc]
    intro hc
    exact hy (vcY_inj C hc)
  rw [WeierstrassCurve.Affine.slope_of_Y_ne rfl hy', WeierstrassCurve.Affine.slope_of_Y_ne rfl hy]
  have hN : 3 * (vcX C x₁) ^ 2 + 2 * (C • W).a₂ * (vcX C x₁) + (C • W).a₄
        - (C • W).a₁ * (vcY C x₁ y₁)
      = ((C.u⁻¹ : Fˣ) : F) ^ 4
        * ((3 * x₁ ^ 2 + 2 * W.a₂ * x₁ + W.a₄ - W.a₁ * y₁)
          - C.s * (y₁ - W.toAffine.negY x₁ y₁)) := by
    simp only [vcX, vcY, variableChange_a₁, variableChange_a₂, variableChange_a₄,
      WeierstrassCurve.Affine.negY]
    ring
  have hDD : vcY C x₁ y₁ - (C • W).toAffine.negY (vcX C x₁) (vcY C x₁ y₁)
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * (y₁ - W.toAffine.negY x₁ y₁) := by
    rw [negY_vc]
    simp only [vcY, WeierstrassCurve.Affine.negY]
    ring
  rw [hN, hDD]
  field_simp

/-! ## ★★★★★★★★点の加法同型 -/

/-- ★★★★★★**変数変換で点を動かす**。 -/
noncomputable def vcPoint (W : WeierstrassCurve F) (C : VariableChange F) :
    W.toAffine.Point → (C • W).toAffine.Point
  | .zero => 0
  | .some x y h => .some (vcX C x) (vcY C x y) ((nonsingular_vc W C x y).2 h)

theorem vcPoint_zero (W : WeierstrassCurve F) (C : VariableChange F) :
    vcPoint W C 0 = 0 := rfl

theorem vcPoint_some (W : WeierstrassCurve F) (C : VariableChange F) (x y : F)
    (h : W.toAffine.Nonsingular x y) :
    vcPoint W C (.some x y h) = .some (vcX C x) (vcY C x y) ((nonsingular_vc W C x y).2 h) := rfl

/-- ★★★★★符号を保つ。 -/
theorem vcPoint_neg (W : WeierstrassCurve F) (C : VariableChange F) (P : W.toAffine.Point) :
    vcPoint W C (-P) = -(vcPoint W C P) := by
  cases P with
  | zero => rfl
  | some x y h =>
    rw [WeierstrassCurve.Affine.Point.neg_some, vcPoint_some, vcPoint_some,
      WeierstrassCurve.Affine.Point.neg_some]
    simp only [negY_vc]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**加法を保つ**。

★3 つの場合(`x` が異なる・`x` が等しく退化・`x` が等しく非退化)に分かれ、
そのそれぞれで mathlib の `add_of_X_ne`・`add_of_Y_eq`・`add_of_Y_ne` に対応する。 -/
theorem vcPoint_add (W : WeierstrassCurve F) (C : VariableChange F) (P Q : W.toAffine.Point) :
    vcPoint W C (P + Q) = vcPoint W C P + vcPoint W C Q := by
  cases P with
  | zero =>
    show vcPoint W C (0 + Q) = vcPoint W C 0 + vcPoint W C Q
    rw [zero_add, vcPoint_zero, zero_add]
  | some x₁ y₁ h₁ =>
    cases Q with
    | zero =>
      show vcPoint W C (_ + 0) = vcPoint W C _ + vcPoint W C 0
      rw [add_zero, vcPoint_zero, add_zero]
    | some x₂ y₂ h₂ =>
      by_cases hx : x₁ = x₂
      · subst hx
        by_cases hy : y₁ = W.toAffine.negY x₁ y₂
        · have hy' : vcY C x₁ y₁ = (C • W).toAffine.negY (vcX C x₁) (vcY C x₁ y₂) := by
            rw [negY_vc, hy]
          rw [WeierstrassCurve.Affine.Point.add_of_Y_eq rfl hy, vcPoint_zero,
            vcPoint_some, vcPoint_some,
            WeierstrassCurve.Affine.Point.add_of_Y_eq rfl hy']
        · have hy' : vcY C x₁ y₁ ≠ (C • W).toAffine.negY (vcX C x₁) (vcY C x₁ y₂) := by
            rw [negY_vc]
            intro hc
            exact hy (vcY_inj C hc)
          rw [WeierstrassCurve.Affine.Point.add_of_Y_ne hy, vcPoint_some, vcPoint_some,
            vcPoint_some, WeierstrassCurve.Affine.Point.add_of_Y_ne hy']
          simp only [WeierstrassCurve.Affine.Point.some.injEq]
          refine ⟨?_, ?_⟩
          · rw [slope_vc_of_X_eq W C h₁.1 h₂.1 rfl hy, addX_vc]
          · rw [slope_vc_of_X_eq W C h₁.1 h₂.1 rfl hy, addY_vc]
      · have hx' : vcX C x₁ ≠ vcX C x₂ := fun hc => hx (vcX_inj C hc)
        rw [WeierstrassCurve.Affine.Point.add_of_X_ne hx, vcPoint_some, vcPoint_some,
          vcPoint_some, WeierstrassCurve.Affine.Point.add_of_X_ne hx']
        simp only [WeierstrassCurve.Affine.Point.some.injEq]
        refine ⟨?_, ?_⟩
        · rw [slope_vc_of_X_ne W C y₁ y₂ hx, addX_vc]
        · rw [slope_vc_of_X_ne W C y₁ y₂ hx, addY_vc]

theorem vcPoint_injective (W : WeierstrassCurve F) (C : VariableChange F) :
    Function.Injective (vcPoint W C) := by
  intro P Q hPQ
  cases P with
  | zero =>
    cases Q with
    | zero => rfl
    | some x y h => simp [vcPoint] at hPQ
  | some x₁ y₁ h₁ =>
    cases Q with
    | zero => simp [vcPoint] at hPQ
    | some x₂ y₂ h₂ =>
      rw [vcPoint_some, vcPoint_some] at hPQ
      simp only [WeierstrassCurve.Affine.Point.some.injEq] at hPQ
      have hx : x₁ = x₂ := vcX_inj C hPQ.1
      subst hx
      have hy : y₁ = y₂ := vcY_inj C hPQ.2
      subst hy
      rfl

theorem vcPoint_surjective (W : WeierstrassCurve F) (C : VariableChange F) :
    Function.Surjective (vcPoint W C) := by
  intro P'
  cases P' with
  | zero => exact ⟨0, rfl⟩
  | some x' y' h' =>
    refine ⟨.some (((C.u : Fˣ) : F) ^ 2 * x' + C.r)
      (((C.u : Fˣ) : F) ^ 3 * y' + ((C.u : Fˣ) : F) ^ 2 * C.s * x' + C.t) ?_, ?_⟩
    · refine (nonsingular_vc W C _ _).1 ?_
      rw [vcX_apply_inv, vcY_apply_inv]
      exact h'
    · rw [vcPoint_some]
      simp only [WeierstrassCurve.Affine.Point.some.injEq]
      exact ⟨vcX_apply_inv C x', vcY_apply_inv C x' y'⟩

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**変数変換は点の加法同型を与える**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def vcPointAddEquiv (W : WeierstrassCurve F) (C : VariableChange F) :
    W.toAffine.Point ≃+ (C • W).toAffine.Point :=
  AddEquiv.mk'
    (Equiv.ofBijective (vcPoint W C) ⟨vcPoint_injective W C, vcPoint_surjective W C⟩)
    (vcPoint_add W C)

end VarChange

/-! ## ★出典の紐付け(`.src`) -/

def equation_vc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(変数変換——方程式の移り方)",
    sectionId := "genell-def-3-3" }

def vcPoint_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(変数変換——加法を保つこと)",
    sectionId := "genell-def-3-3" }

def vcPointAddEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(変数変換が与える点の加法同型)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
