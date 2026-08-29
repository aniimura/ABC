/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.EllipticCurve.VariableChange
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Formula
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point
import ABC3.Found.GenEll.Velu
import ABC3.Meta.Claim

/-!
# 変数変換の点への作用 —— `Equation` の側

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.17。

## ★★このファイルで測ったこと

★**mathlib には `VariableChange` の点への作用が無い**（2026-08-29 実測）。
`Mathlib/AlgebraicGeometry/EllipticCurve/Affine/Point.lean` の `Point.map` は
**環準同型**（正確には `F →ₐ[S] K`）に対するもので、変数変換に対するものは無い。
`equation_iff_variableChange` は `(1, x, 0, y)` という特別な変数変換で
`(0,0)` へ移す形しかない。

## ★このファイルで建てるもの

変数変換 `C = (u, r, s, t)` は座標を

    x = u²·x′ + r,   y = u³·y′ + u²·s·x′ + t

で移す。逆に解くと

    x′ = u⁻²·(x − r),   y′ = u⁻³·(y − s·(x − r) − t)

これを `vcX`・`vcY` と書く。★★**`Equation` はこの対応で保たれる**——
これが `Found/GenEll/Uniformization.lean` の解析側の結論を
`E.map σ` の側へ戻すために要る。

☆群法則が保たれること（`AddEquiv` になること）は本ファイルでは示していない。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve

variable {R : Type*} [CommRing R]

/-! ## ★★★★★座標の変換 -/

/-- ★★★★★変数変換 `C = (u, r, s, t)` による `x` 座標の移動

    `x′ = u⁻²·(x − r)` -/
def vcX (C : VariableChange R) (x : R) : R := ((C.u⁻¹ : Rˣ) : R) ^ 2 * (x - C.r)

/-- ★★★★★変数変換 `C = (u, r, s, t)` による `y` 座標の移動

    `y′ = u⁻³·(y − s·(x − r) − t)` -/
def vcY (C : VariableChange R) (x y : R) : R :=
  ((C.u⁻¹ : Rˣ) : R) ^ 3 * (y - C.s * (x - C.r) - C.t)

/-! ## ★★★★★★★★`Equation` は保たれる -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**変数変換の方程式の側の恒等式**

    `F_{C•W}(x′, y′) = u⁻⁶ · F_W(x, y)`

（`F_W(x,y) = y² + a₁xy + a₃y − x³ − a₂x² − a₄x − a₆`）

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★`variableChange_a*` の各式は `u⁻¹` **だけ**で書かれているので、
`u⁻¹` を不定元と見た多項式恒等式であり、`ring` で閉じる。 -/
theorem variableChange_equation_poly (C : VariableChange R) (W : WeierstrassCurve R)
    (x y : R) :
    (vcY C x y) ^ 2 + (C • W).a₁ * (vcX C x) * (vcY C x y) + (C • W).a₃ * (vcY C x y)
        - ((vcX C x) ^ 3 + (C • W).a₂ * (vcX C x) ^ 2 + (C • W).a₄ * (vcX C x)
          + (C • W).a₆)
      = ((C.u⁻¹ : Rˣ) : R) ^ 6
        * (y ^ 2 + W.a₁ * x * y + W.a₃ * y
          - (x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆)) := by
  simp only [vcX, vcY, WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₂,
    WeierstrassCurve.variableChange_a₃, WeierstrassCurve.variableChange_a₄,
    WeierstrassCurve.variableChange_a₆]
  ring

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Equation` は変数変換で保たれる**

    `(C • W).Equation (vcX C x) (vcY C x y) ↔ W.Equation x y`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★☆**mathlib に無い**（`Point.map` は環準同型に対するものだけ）。
☆これで `latticeCurve (P σ) = C σ • (E.map σ)` の両側の点が対応する。 -/
theorem equation_variableChange (C : VariableChange R) (W : WeierstrassCurve R)
    (x y : R) :
    (C • W).toAffine.Equation (vcX C x) (vcY C x y) ↔ W.toAffine.Equation x y := by
  have hu : ((C.u : Rˣ) : R) * ((C.u⁻¹ : Rˣ) : R) = 1 := C.u.mul_inv
  have hcancel : ∀ z : R, ((C.u⁻¹ : Rˣ) : R) ^ 6 * z = 0 → z = 0 := by
    intro z hz
    calc z = (((C.u : Rˣ) : R) * ((C.u⁻¹ : Rˣ) : R)) ^ 6 * z := by rw [hu]; ring
      _ = ((C.u : Rˣ) : R) ^ 6 * (((C.u⁻¹ : Rˣ) : R) ^ 6 * z) := by ring
      _ = 0 := by rw [hz, mul_zero]
  rw [WeierstrassCurve.Affine.equation_iff', WeierstrassCurve.Affine.equation_iff']
  constructor
  · intro h
    have hkey := variableChange_equation_poly C W x y
    rw [h] at hkey
    exact hcancel _ hkey.symm
  · intro h
    have hkey := variableChange_equation_poly C W x y
    rw [h, mul_zero] at hkey
    exact hkey

/-! ## ★★★★★★★★群法則の部品の変換則 -/

/-- ★★★★★変数変換による傾きの移動 `L′ = u⁻¹·(L − s)`。 -/
def vcSlope (C : VariableChange R) (L : R) : R := ((C.u⁻¹ : Rˣ) : R) * (L - C.s)

/-- ★★★★★★★★★★★★**`negY` は変数変換と両立する**

    `negY_{C•W}(x′, y′) = vcY C x (negY_W(x, y))` -/
theorem negY_variableChange (C : VariableChange R) (W : WeierstrassCurve R) (x y : R) :
    (C • W).toAffine.negY (vcX C x) (vcY C x y)
      = vcY C x (W.toAffine.negY x y) := by
  simp only [WeierstrassCurve.Affine.negY, vcX, vcY,
    WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₃]
  ring

/-- ★★★★★★★★★★★★**`addX` は変数変換と両立する**

    `addX_{C•W}(x₁′, x₂′, L′) = vcX C (addX_W(x₁, x₂, L))` -/
theorem addX_variableChange (C : VariableChange R) (W : WeierstrassCurve R)
    (x₁ x₂ L : R) :
    (C • W).toAffine.addX (vcX C x₁) (vcX C x₂) (vcSlope C L)
      = vcX C (W.toAffine.addX x₁ x₂ L) := by
  simp only [WeierstrassCurve.Affine.addX, vcX, vcSlope,
    WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₂]
  ring

/-- ★★★★★★★★★★★★**`negAddY` は変数変換と両立する**

    `negAddY_{C•W}(x₁′, x₂′, y₁′, L′) = vcY C (addX_W) (negAddY_W)` -/
theorem negAddY_variableChange (C : VariableChange R) (W : WeierstrassCurve R)
    (x₁ x₂ y₁ L : R) :
    (C • W).toAffine.negAddY (vcX C x₁) (vcX C x₂) (vcY C x₁ y₁) (vcSlope C L)
      = vcY C (W.toAffine.addX x₁ x₂ L) (W.toAffine.negAddY x₁ x₂ y₁ L) := by
  simp only [WeierstrassCurve.Affine.negAddY, addX_variableChange, vcX, vcY, vcSlope,
    WeierstrassCurve.Affine.addX, WeierstrassCurve.variableChange_a₁,
    WeierstrassCurve.variableChange_a₂]
  ring

/-- ★★★★★★★★★★★★★★★★**`addY` は変数変換と両立する**

    `addY_{C•W}(x₁′, x₂′, y₁′, L′) = vcY C (addX_W) (addY_W)`

★`addY = negY(addX, negAddY)` なので `negY` と `negAddY` の変換則から出る。 -/
theorem addY_variableChange (C : VariableChange R) (W : WeierstrassCurve R)
    (x₁ x₂ y₁ L : R) :
    (C • W).toAffine.addY (vcX C x₁) (vcX C x₂) (vcY C x₁ y₁) (vcSlope C L)
      = vcY C (W.toAffine.addX x₁ x₂ L) (W.toAffine.addY x₁ x₂ y₁ L) := by
  rw [WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.addY,
    ← negY_variableChange C W (W.toAffine.addX x₁ x₂ L),
    ← negAddY_variableChange, addX_variableChange]

def negY_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(negY は変数変換と両立する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def addX_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(addX は変数変換と両立する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def addY_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(addY は変数変換と両立する。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★`slope` の変換則 -/

section Field

variable {F : Type*} [Field F]

/-- ★★★★`vcX C` は単射。 -/
theorem vcX_injective (C : VariableChange F) : Function.Injective (vcX C) := by
  intro a b hab
  have hu : ((C.u⁻¹ : Fˣ) : F) ^ 2 ≠ 0 := pow_ne_zero 2 (Units.ne_zero _)
  have h := mul_left_cancel₀ hu hab
  linear_combination h

/-- ★★★★`vcY C x` は単射。 -/
theorem vcY_injective (C : VariableChange F) (x : F) : Function.Injective (vcY C x) := by
  intro a b hab
  have hu : ((C.u⁻¹ : Fˣ) : F) ^ 3 ≠ 0 := pow_ne_zero 3 (Units.ne_zero _)
  have h := mul_left_cancel₀ hu hab
  linear_combination h

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**`slope` の変換則——`x₁ ≠ x₂` の場合**

    `slope_{C•W}(x₁′, x₂′, y₁′, y₂′) = u⁻¹·(slope_W(x₁,x₂,y₁,y₂) − s)` -/
theorem slope_variableChange_of_ne (C : VariableChange F) (W : WeierstrassCurve F)
    {x₁ x₂ : F} (y₁ y₂ : F) (hx : x₁ ≠ x₂) :
    (C • W).toAffine.slope (vcX C x₁) (vcX C x₂) (vcY C x₁ y₁) (vcY C x₂ y₂)
      = vcSlope C (W.toAffine.slope x₁ x₂ y₁ y₂) := by
  have hx' : vcX C x₁ ≠ vcX C x₂ := fun hc => hx (vcX_injective C hc)
  rw [WeierstrassCurve.Affine.slope, if_neg hx', WeierstrassCurve.Affine.slope, if_neg hx]
  have hu : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hd : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hx
  have hnum : vcY C x₁ y₁ - vcY C x₂ y₂
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * ((y₁ - y₂) - C.s * (x₁ - x₂)) := by
    simp only [vcY]; ring
  have hden : vcX C x₁ - vcX C x₂ = ((C.u⁻¹ : Fˣ) : F) ^ 2 * (x₁ - x₂) := by
    simp only [vcX]; ring
  rw [hnum, hden, vcSlope, div_eq_iff (mul_ne_zero (pow_ne_zero 2 hu) hd)]
  field_simp

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★**`slope` の変換則——接線の場合**

★`x₁ = x₂` かつ `y₁ ≠ negY(x, y₂)` かつ `y₁ ≠ negY(x, y₁)`（＝接線が定義される）。 -/
theorem slope_variableChange_of_eq (C : VariableChange F) (W : WeierstrassCurve F)
    (x y₁ y₂ : F) (hy : y₁ ≠ W.toAffine.negY x y₂)
    (hy' : y₁ ≠ W.toAffine.negY x y₁) :
    (C • W).toAffine.slope (vcX C x) (vcX C x) (vcY C x y₁) (vcY C x y₂)
      = vcSlope C (W.toAffine.slope x x y₁ y₂) := by
  have hy2 : vcY C x y₁ ≠ (C • W).toAffine.negY (vcX C x) (vcY C x y₂) := by
    rw [negY_variableChange]
    exact fun hc => hy (vcY_injective C x hc)
  rw [WeierstrassCurve.Affine.slope, if_pos rfl, if_neg hy2,
    WeierstrassCurve.Affine.slope, if_pos rfl, if_neg hy]
  have hu : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hd : y₁ - W.toAffine.negY x y₁ ≠ 0 := sub_ne_zero.2 hy'
  have hd2 : vcY C x y₁ - (C • W).toAffine.negY (vcX C x) (vcY C x y₁)
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * (y₁ - W.toAffine.negY x y₁) := by
    rw [negY_variableChange]
    simp only [vcY, WeierstrassCurve.Affine.negY]
    ring
  have hnum : 3 * (vcX C x) ^ 2 + 2 * (C • W).a₂ * (vcX C x) + (C • W).a₄
        - (C • W).a₁ * (vcY C x y₁)
      = ((C.u⁻¹ : Fˣ) : F) ^ 4
        * ((3 * x ^ 2 + 2 * W.a₂ * x + W.a₄ - W.a₁ * y₁)
          - C.s * (y₁ - W.toAffine.negY x y₁)) := by
    simp only [vcX, vcY, WeierstrassCurve.Affine.negY,
      WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₂,
      WeierstrassCurve.variableChange_a₄]
    ring
  rw [hd2, hnum, vcSlope, div_eq_iff (mul_ne_zero (pow_ne_zero 3 hu) hd)]
  field_simp

end Field

def slope_variableChange_of_ne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(slope の変換則——x₁ ≠ x₂ の場合。★無条件)",
    sectionId := "genell-lemma-3-5" }

def slope_variableChange_of_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(slope の変換則——接線の場合。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★`Nonsingular` も保たれる -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Nonsingular` は変数変換で保たれる**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★2 つの偏微分は

    `A′ = u⁻⁴·(A + s·B)`,   `B′ = u⁻³·B`
    （`A = a₁y − (3x² + 2a₂x + a₄)`、`B = 2y + a₁x + a₃`）

と変換される。`B = 0` なら `A′ = u⁻⁴A` なので、
`(A′ ≠ 0 ∨ B′ ≠ 0) ⟺ (A ≠ 0 ∨ B ≠ 0)`。

★★★☆**`Equation`（第 682）と合わせて、点が点へ移ることが取れた。** -/
theorem nonsingular_variableChange {F : Type*} [Field F] (C : VariableChange F)
    (W : WeierstrassCurve F) (x y : F) :
    (C • W).toAffine.Nonsingular (vcX C x) (vcY C x y) ↔ W.toAffine.Nonsingular x y := by
  have hu : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hA : (C • W).a₁ * (vcY C x y)
        - (3 * (vcX C x) ^ 2 + 2 * (C • W).a₂ * (vcX C x) + (C • W).a₄)
      = ((C.u⁻¹ : Fˣ) : F) ^ 4
        * ((W.a₁ * y - (3 * x ^ 2 + 2 * W.a₂ * x + W.a₄))
          + C.s * (2 * y + W.a₁ * x + W.a₃)) := by
    simp only [vcX, vcY, WeierstrassCurve.variableChange_a₁,
      WeierstrassCurve.variableChange_a₂, WeierstrassCurve.variableChange_a₄]
    ring
  have hB : 2 * (vcY C x y) + (C • W).a₁ * (vcX C x) + (C • W).a₃
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * (2 * y + W.a₁ * x + W.a₃) := by
    simp only [vcX, vcY, WeierstrassCurve.variableChange_a₁,
      WeierstrassCurve.variableChange_a₃]
    ring
  have key : ∀ a b : F,
      (((C.u⁻¹ : Fˣ) : F) ^ 4 * (a + C.s * b) ≠ 0 ∨ ((C.u⁻¹ : Fˣ) : F) ^ 3 * b ≠ 0)
        ↔ (a ≠ 0 ∨ b ≠ 0) := by
    intro a b
    have h4 : ((C.u⁻¹ : Fˣ) : F) ^ 4 ≠ 0 := pow_ne_zero 4 hu
    have h3 : ((C.u⁻¹ : Fˣ) : F) ^ 3 ≠ 0 := pow_ne_zero 3 hu
    rw [mul_ne_zero_iff, mul_ne_zero_iff]
    simp only [h4, h3, true_and, ne_eq, not_false_eq_true]
    constructor
    · rintro (h | h)
      · rcases eq_or_ne b 0 with hb | hb
        · left; simpa [hb] using h
        · right; exact hb
      · right; exact h
    · rintro (h | h)
      · rcases eq_or_ne b 0 with hb | hb
        · left; simpa [hb] using h
        · right; exact hb
      · right; exact h
  rw [WeierstrassCurve.Affine.nonsingular_iff', WeierstrassCurve.Affine.nonsingular_iff',
    hA, hB, equation_variableChange, key]

def nonsingular_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Nonsingular も変数変換で保たれる——mathlib に無い。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★点の写像 -/

section Point

variable {F : Type*} [Field F]

/-- ★★★★★`Point.some` の合同。 -/
theorem point_some_congr' {W : WeierstrassCurve F} {x₁ y₁ x₂ y₂ : F}
    {h₁ : W.toAffine.Nonsingular x₁ y₁} {h₂ : W.toAffine.Nonsingular x₂ y₂}
    (hx : x₁ = x₂) (hy : y₁ = y₂) :
    (WeierstrassCurve.Affine.Point.some x₁ y₁ h₁ : W.toAffine.Point)
      = WeierstrassCurve.Affine.Point.some x₂ y₂ h₂ := by
  subst hx; subst hy; rfl

/-- ★★★★★★★★★★★★★★★★**変数変換による点の写像**

    `O ↦ O`,  `(x, y) ↦ (u⁻²(x − r), u⁻³(y − s(x − r) − t))`

★`Nonsingular` が保たれる（第 685）ので well-defined。 -/
noncomputable def vcPoint (C : VariableChange F) (W : WeierstrassCurve F) :
    W.toAffine.Point → (C • W).toAffine.Point
  | .zero => 0
  | .some x y h => .some (vcX C x) (vcY C x y) ((nonsingular_variableChange C W x y).2 h)

@[simp] theorem vcPoint_zero (C : VariableChange F) (W : WeierstrassCurve F) :
    vcPoint C W 0 = 0 := rfl

@[simp] theorem vcPoint_some (C : VariableChange F) (W : WeierstrassCurve F)
    {x y : F} (h : W.toAffine.Nonsingular x y) :
    vcPoint C W (.some x y h)
      = .some (vcX C x) (vcY C x y) ((nonsingular_variableChange C W x y).2 h) := rfl

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**変数変換による点の写像は加法を保つ**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★☆**これで `W.Point ≃+ (C • W).Point` が取れる**——mathlib には
`VariableChange` の点への作用が無いので、これは新しく建てたものである。
☆位数が保たれるので、`latticeCurve (P σ) = C σ • (E.map σ)` の両側で
「位数 `l` の点」が対応する。

★場合分けは 2 つだけ:

* 垂直（`x₁ = x₂` かつ `y₁ = negY(x₂,y₂)`）—— 両側とも `O`（`add_of_Y_eq`）
* それ以外 —— `add_some` で座標が出るので、
  `slope`（第 684）・`addX`・`addY`（第 683）の変換則で一致する -/
theorem vcPoint_add (C : VariableChange F) (W : WeierstrassCurve F)
    (Pt Qt : W.toAffine.Point) :
    vcPoint C W (Pt + Qt) = vcPoint C W Pt + vcPoint C W Qt := by
  rcases Pt with _ | ⟨x₁, y₁, h₁⟩
  · change vcPoint C W (0 + Qt) = vcPoint C W 0 + vcPoint C W Qt
    rw [zero_add, vcPoint_zero, zero_add]
  rcases Qt with _ | ⟨x₂, y₂, h₂⟩
  · change vcPoint C W (_ + 0) = vcPoint C W _ + vcPoint C W 0
    rw [add_zero, vcPoint_zero, add_zero]
  by_cases hxy : x₁ = x₂ ∧ y₁ = W.toAffine.negY x₂ y₂
  · have hx' : vcX C x₁ = vcX C x₂ := by rw [hxy.1]
    have hy' : vcY C x₁ y₁ = (C • W).toAffine.negY (vcX C x₂) (vcY C x₂ y₂) := by
      rw [negY_variableChange, hxy.1, hxy.2]
    rw [WeierstrassCurve.Affine.Point.add_of_Y_eq hxy.1 hxy.2, vcPoint_zero,
      vcPoint_some, vcPoint_some,
      WeierstrassCurve.Affine.Point.add_of_Y_eq hx' hy']
  · have hxy' : ¬(vcX C x₁ = vcX C x₂
        ∧ vcY C x₁ y₁ = (C • W).toAffine.negY (vcX C x₂) (vcY C x₂ y₂)) := by
      rintro ⟨ha, hb⟩
      have hx12 : x₁ = x₂ := vcX_injective C ha
      subst hx12
      rw [negY_variableChange] at hb
      exact hxy ⟨rfl, vcY_injective C x₁ hb⟩
    have hslope : (C • W).toAffine.slope (vcX C x₁) (vcX C x₂) (vcY C x₁ y₁)
        (vcY C x₂ y₂) = vcSlope C (W.toAffine.slope x₁ x₂ y₁ y₂) := by
      by_cases hx : x₁ = x₂
      · subst hx
        have hy : y₁ ≠ W.toAffine.negY x₁ y₂ := fun hc => hxy ⟨rfl, hc⟩
        have hy2 : y₁ = y₂ :=
          WeierstrassCurve.Affine.Y_eq_of_Y_ne h₁.left h₂.left rfl hy
        subst hy2
        exact slope_variableChange_of_eq C W x₁ y₁ y₁ hy hy
      · exact slope_variableChange_of_ne C W y₁ y₂ hx
    rw [WeierstrassCurve.Affine.Point.add_some hxy, vcPoint_some, vcPoint_some,
      vcPoint_some, WeierstrassCurve.Affine.Point.add_some hxy']
    refine point_some_congr' ?_ ?_
    · rw [hslope, addX_variableChange]
    · rw [hslope, addY_variableChange]

/-- ★★★★★★★★★★**点の写像は単射**——`vcX`・`vcY` が単射だから。 -/
theorem vcPoint_injective (C : VariableChange F) (W : WeierstrassCurve F) :
    Function.Injective (vcPoint C W) := by
  rintro (_ | ⟨x₁, y₁, h₁⟩) (_ | ⟨x₂, y₂, h₂⟩) hab
  · rfl
  · exact absurd hab (by simp [vcPoint])
  · exact absurd hab (by simp [vcPoint])
  · simp only [vcPoint_some] at hab
    have hx : vcX C x₁ = vcX C x₂ := by
      injection hab with hx hy
    have hy : vcY C x₁ y₁ = vcY C x₂ y₂ := by
      injection hab with hx hy
    have hx12 : x₁ = x₂ := vcX_injective C hx
    subst hx12
    exact point_some_congr' rfl (vcY_injective C x₁ hy)

open scoped Classical in
/-- ★★★★★★★★`n • ` と可換。 -/
theorem vcPoint_nsmul (C : VariableChange F) (W : WeierstrassCurve F)
    [W.IsElliptic] [(C • W).IsElliptic] (Pt : W.toAffine.Point) (n : ℕ) :
    vcPoint C W (n • Pt) = n • vcPoint C W Pt := by
  induction n with
  | zero => simp only [zero_nsmul, vcPoint_zero]
  | succ k ih => rw [succ_nsmul, succ_nsmul, vcPoint_add, ih]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**変数変換は点の位数を保つ**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★☆**これが Galois 降下に要るものである**——
`latticeCurve (P σ) = C σ • (E.map σ)` の両側で「位数 `l` の点」が対応する。 -/
theorem addOrderOf_vcPoint (C : VariableChange F) (W : WeierstrassCurve F)
    [W.IsElliptic] [(C • W).IsElliptic] (Pt : W.toAffine.Point) :
    addOrderOf (vcPoint C W Pt) = addOrderOf Pt := by
  refine Nat.dvd_antisymm ?_ ?_
  · refine addOrderOf_dvd_of_nsmul_eq_zero ?_
    rw [← vcPoint_nsmul, addOrderOf_nsmul_eq_zero, vcPoint_zero]
  · refine addOrderOf_dvd_of_nsmul_eq_zero ?_
    refine vcPoint_injective C W ?_
    rw [vcPoint_nsmul, vcPoint_zero, addOrderOf_nsmul_eq_zero]

/-- ★★★★★★逆向きの座標 `x = u²x′ + r`。 -/
theorem vcX_inv (C : VariableChange F) (x' : F) :
    vcX C (((C.u : Fˣ) : F) ^ 2 * x' + C.r) = x' := by
  have hu : ((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) = 1 := C.u.mul_inv
  simp only [vcX]
  linear_combination (x' * (((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) + 1)) * hu

/-- ★★★★★★逆向きの座標 `y = u³y′ + u²s·x′ + t`。 -/
theorem vcY_inv (C : VariableChange F) (x' y' : F) :
    vcY C (((C.u : Fˣ) : F) ^ 2 * x' + C.r)
        (((C.u : Fˣ) : F) ^ 3 * y' + ((C.u : Fˣ) : F) ^ 2 * C.s * x' + C.t) = y' := by
  have hu : ((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) = 1 := C.u.mul_inv
  simp only [vcY]
  linear_combination (y' * ((((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F)) ^ 2
      + ((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) + 1)
    + ((C.u : Fˣ) : F) ^ 2 * C.s * x' * ((C.u⁻¹ : Fˣ) : F) ^ 3
      * (((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) + 1) * 0) * hu

/-- ★★★★★★★★★★**点の写像は全射**——逆向きの変数変換で戻せる。 -/
theorem vcPoint_surjective (C : VariableChange F) (W : WeierstrassCurve F) :
    Function.Surjective (vcPoint C W) := by
  rintro (_ | ⟨x', y', h'⟩)
  · exact ⟨0, rfl⟩
  · refine ⟨.some (((C.u : Fˣ) : F) ^ 2 * x' + C.r)
      (((C.u : Fˣ) : F) ^ 3 * y' + ((C.u : Fˣ) : F) ^ 2 * C.s * x' + C.t) ?_, ?_⟩
    · refine (nonsingular_variableChange C W _ _).1 ?_
      rw [vcX_inv, vcY_inv]
      exact h'
    · rw [vcPoint_some]
      exact point_some_congr' (vcX_inv C x') (vcY_inv C x' y')

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**変数変換は点の加法群の同型を与える**

    `W.Point ≃+ (C • W).Point`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★☆**mathlib にはこれが無い**（`Point.map` は環準同型に対するものだけ）。
☆位数が保たれる（第 686）ので、`latticeCurve (P σ) = C σ • (E.map σ)` の
両側で「位数 `l` の点」と「その生成する部分群」が 1 対 1 に対応する。 -/
noncomputable def vcEquiv (C : VariableChange F) (W : WeierstrassCurve F)
    [W.IsElliptic] [(C • W).IsElliptic] :
    W.toAffine.Point ≃+ (C • W).toAffine.Point :=
  AddEquiv.mk'
    (Equiv.ofBijective (vcPoint C W) ⟨vcPoint_injective C W, vcPoint_surjective C W⟩)
    (vcPoint_add C W)

open scoped Classical in
@[simp] theorem vcEquiv_apply (C : VariableChange F) (W : WeierstrassCurve F)
    [W.IsElliptic] [(C • W).IsElliptic] (Pt : W.toAffine.Point) :
    vcEquiv C W Pt = vcPoint C W Pt := rfl

end Point

def vcPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換による点の写像——mathlib に無い)",
    sectionId := "genell-lemma-3-5" }

def vcEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換は点の加法群の同型を与える——mathlib に無い)",
    sectionId := "genell-lemma-3-5" }

def vcPoint_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(点の写像は全射——逆向きの変数変換で戻せる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def addOrderOf_vcPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換は点の位数を保つ——Galois 降下の要。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vcPoint_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換による点の写像は加法を保つ——mathlib に無い。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★Vélu の量の変換則 -/

section VeluVC

variable {F : Type*} [Field F]

/-- ★★★★★★`B ≔ 2y + a₁x + a₃`（`= y − negY(x,y)`）は `±` で打ち消し合う。

★★☆これが `Σ_{w ∈ T∖{0}} ℘′(w) = 0`（第 670）の代数版である。 -/
theorem veluB_add_negY (W : WeierstrassCurve F) (x y : F) :
    (2 * y + W.a₁ * x + W.a₃) + (2 * (W.toAffine.negY x y) + W.a₁ * x + W.a₃) = 0 := by
  simp only [WeierstrassCurve.Affine.negY]
  ring

/-- ★★★★★★★★★★★★★★**`g^x_Q` の一般の変数変換による変換則**

    `g^x_{C•W}(x′, y′) = u⁻⁴·(g^x_W(x, y) − s·B)`,  `B = 2y + a₁x + a₃`

★純粋なスケーリング（`veluGx_scale`）と違い、`s ≠ 0` のときは `−s·B` の補正が付く。
☆だが `S` が `±` で閉じていれば `Σ_S B = 0`（`veluB_add_negY`）なので、
**和を取れば補正は消える**。 -/
theorem veluV2_variableChange (C : VariableChange F) (W : WeierstrassCurve F)
    (x y : F) :
    veluV2 (C • W) (vcX C x) (vcY C x y)
      = ((C.u⁻¹ : Fˣ) : F) ^ 4
        * (veluV2 W x y - C.s * (2 * y + W.a₁ * x + W.a₃)) := by
  simp only [veluV2, veluGx, vcX, vcY, WeierstrassCurve.variableChange_a₁,
    WeierstrassCurve.variableChange_a₂, WeierstrassCurve.variableChange_a₄]
  ring

/-- ★★★★★★★★★★★★**`g^y_Q` の変換則**——`g^y = −B` は重さ 3。 -/
theorem veluGy_variableChange (C : VariableChange F) (W : WeierstrassCurve F)
    (x y : F) :
    veluGy (C • W) (vcX C x) (vcY C x y)
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * veluGy W x y := by
  simp only [veluGy, vcX, vcY, WeierstrassCurve.variableChange_a₁,
    WeierstrassCurve.variableChange_a₃]
  ring

/-- ★★★★★★★★★★★★**`u_Q` の変換則**——重さ 6（補正なし）。 -/
theorem veluU_variableChange (C : VariableChange F) (W : WeierstrassCurve F)
    (x y : F) :
    veluU (C • W) (vcX C x) (vcY C x y)
      = ((C.u⁻¹ : Fˣ) : F) ^ 6 * veluU W x y := by
  rw [veluU, veluGy_variableChange, veluU]
  ring

end VeluVC

def veluV2_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(g^x_Q の一般の変数変換による変換則。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluB_add_negY.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(B = 2y + a₁x + a₃ は ± で打ち消し合う——Σ℘′(w) = 0 の代数版。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★商曲線の一般変数変換 -/

section VeluFullVC

variable {F : Type*} [Field F]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Vélu の商曲線の一般の変数変換による変換則**

    `veluCurve (C•W) (u⁻⁴v) (u⁻⁶(w − r·v)) = C • veluCurve W v w`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★純粋なスケーリング（`veluCurve_scale`）では `w ↦ u⁻⁶w` だったが、
一般には `r` の分だけ `w ↦ u⁻⁶(w − r·v)` とずれる。
☆`b₂(C•W) = u⁻²(b₂ + 12r)` が効く。 -/
theorem veluCurve_variableChange (C : VariableChange F) (W : WeierstrassCurve F)
    (v w : F) :
    veluCurve (C • W) (((C.u⁻¹ : Fˣ) : F) ^ 4 * v)
        (((C.u⁻¹ : Fˣ) : F) ^ 6 * (w - C.r * v))
      = C • veluCurve W v w := by
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · simp only [veluCurve, WeierstrassCurve.variableChange_a₁]
  · simp only [veluCurve, WeierstrassCurve.variableChange_a₂]
  · simp only [veluCurve, WeierstrassCurve.variableChange_a₃]
  · simp only [veluCurve, WeierstrassCurve.variableChange_a₁,
      WeierstrassCurve.variableChange_a₂, WeierstrassCurve.variableChange_a₃,
      WeierstrassCurve.variableChange_a₄]
    ring
  · simp only [veluCurve, WeierstrassCurve.variableChange_a₁,
      WeierstrassCurve.variableChange_a₂, WeierstrassCurve.variableChange_a₃,
      WeierstrassCurve.variableChange_a₄, WeierstrassCurve.variableChange_a₆,
      WeierstrassCurve.b₂]
    ring

/-- ★★★★点の対の変数変換は単射。 -/
theorem vcPair_injective (C : VariableChange F) :
    Function.Injective (fun Q : F × F => (vcX C Q.1, vcY C Q.1 Q.2)) := by
  intro a b hab
  have h1 : vcX C a.1 = vcX C b.1 := congrArg Prod.fst hab
  have ha : a.1 = b.1 := vcX_injective C h1
  have h2 : vcY C a.1 a.2 = vcY C b.1 b.2 := congrArg Prod.snd hab
  rw [ha] at h2
  exact Prod.ext ha (vcY_injective C b.1 h2)

open scoped Classical in
/-- ★★★★★★★★★★★★★★`veluVFull` の一般変数変換——`Σ_S B = 0` なら補正は消える。 -/
theorem veluVFull_variableChange (C : VariableChange F) (W : WeierstrassCurve F)
    (S : Finset (F × F))
    (hB : ∑ Q ∈ S, (2 * Q.2 + W.a₁ * Q.1 + W.a₃) = 0) :
    veluVFull (C • W) (S.image (fun Q => (vcX C Q.1, vcY C Q.1 Q.2)))
      = ((C.u⁻¹ : Fˣ) : F) ^ 4 * veluVFull W S := by
  rw [veluVFull, veluVFull, Finset.sum_image (fun a _ b _ hab => vcPair_injective C hab)]
  have hstep : ∀ Q ∈ S, veluV2 (C • W) (vcX C Q.1) (vcY C Q.1 Q.2)
      = ((C.u⁻¹ : Fˣ) : F) ^ 4
        * (veluV2 W Q.1 Q.2 - C.s * (2 * Q.2 + W.a₁ * Q.1 + W.a₃)) :=
    fun Q _ => veluV2_variableChange C W Q.1 Q.2
  rw [Finset.sum_congr rfl hstep, ← Finset.mul_sum, Finset.sum_sub_distrib,
    ← Finset.mul_sum, hB, mul_zero, sub_zero]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★`veluWFull` の一般変数変換

    `veluWFull (C•W) (vc S) = u⁻⁶·(veluWFull W S − r·veluVFull W S)`

★`Σ_S B = 0` と `Σ_S B·x = 0` が要る——どちらも `S` が `±` で閉じていれば従う。 -/
theorem veluWFull_variableChange (C : VariableChange F) (W : WeierstrassCurve F)
    (S : Finset (F × F))
    (hB : ∑ Q ∈ S, (2 * Q.2 + W.a₁ * Q.1 + W.a₃) = 0)
    (hBx : ∑ Q ∈ S, (2 * Q.2 + W.a₁ * Q.1 + W.a₃) * Q.1 = 0) :
    veluWFull (C • W) (S.image (fun Q => (vcX C Q.1, vcY C Q.1 Q.2)))
      = ((C.u⁻¹ : Fˣ) : F) ^ 6
        * (veluWFull W S - C.r * veluVFull W S) := by
  rw [veluWFull, veluWFull, veluVFull,
    Finset.sum_image (fun a _ b _ hab => vcPair_injective C hab)]
  have hstep : ∀ Q ∈ S,
      veluU (C • W) (vcX C Q.1) (vcY C Q.1 Q.2) / 2
        + veluV2 (C • W) (vcX C Q.1) (vcY C Q.1 Q.2) * vcX C Q.1
      = ((C.u⁻¹ : Fˣ) : F) ^ 6
        * ((veluU W Q.1 Q.2 / 2 + veluV2 W Q.1 Q.2 * Q.1)
          - C.r * veluV2 W Q.1 Q.2
          - C.s * ((2 * Q.2 + W.a₁ * Q.1 + W.a₃) * Q.1)
          + C.s * C.r * (2 * Q.2 + W.a₁ * Q.1 + W.a₃)) := by
    intro Q _
    rw [veluU_variableChange, veluV2_variableChange, vcX]
    ring
  rw [Finset.sum_congr rfl hstep, ← Finset.mul_sum]
  congr 1
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_sub_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, hB, hBx,
    mul_zero, mul_zero, sub_zero, add_zero]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`H∖{O}` 全体で書いた Vélu の商は一般の変数変換と両立する**

    `veluQuotientFull (C•W) (vc S) = C • veluQuotientFull W S`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★☆**これで `L` 上の `E/H` と `ℂ` 上の `latticeCurve Λ′` が
一意化の変数変換を挟んで結びつく**——`Lemma 3.5` の最後の橋である。 -/
theorem veluQuotientFull_variableChange (C : VariableChange F) (W : WeierstrassCurve F)
    (S : Finset (F × F))
    (hB : ∑ Q ∈ S, (2 * Q.2 + W.a₁ * Q.1 + W.a₃) = 0)
    (hBx : ∑ Q ∈ S, (2 * Q.2 + W.a₁ * Q.1 + W.a₃) * Q.1 = 0) :
    veluQuotientFull (C • W) (S.image (fun Q => (vcX C Q.1, vcY C Q.1 Q.2)))
      = C • veluQuotientFull W S := by
  rw [veluQuotientFull, veluVFull_variableChange C W S hB,
    veluWFull_variableChange C W S hB hBx, veluQuotientFull,
    veluCurve_variableChange]

end VeluFullVC

def veluCurve_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商曲線の一般の変数変換による変換則。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(H∖{O} 全体で書いた商は一般の変数変換と両立する——最後の橋)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★点集合の引き戻し -/

section Pullback

variable {F : Type*} [Field F]

/-- ★★★★★★逆向きの座標変換 `(x′, y′) ↦ (u²x′ + r, u³y′ + u²s·x′ + t)`。 -/
def vcInvPair (C : VariableChange F) (Q : F × F) : F × F :=
  (((C.u : Fˣ) : F) ^ 2 * Q.1 + C.r,
    ((C.u : Fˣ) : F) ^ 3 * Q.2 + ((C.u : Fˣ) : F) ^ 2 * C.s * Q.1 + C.t)

/-- ★★★★★★`vcPair ∘ vcInvPair = id`。 -/
@[simp] theorem vcPair_vcInvPair (C : VariableChange F) (Q : F × F) :
    (vcX C (vcInvPair C Q).1, vcY C (vcInvPair C Q).1 (vcInvPair C Q).2) = Q :=
  Prod.ext (vcX_inv C Q.1) (vcY_inv C Q.1 Q.2)

open scoped Classical in
/-- ★★★★★★★★★★★★★★**点集合は引き戻せる**

    `(S.image (vcInvPair C)).image (vcPair C) = S`

★★☆これで「`ℂ` 上で得た代表点集合 `S` は、`C` を外した曲線 `W` の上の
点集合の像である」と言える——`L` 上の `E/H` を作るための引き戻しである。 -/
theorem image_vcPair_vcInvPair (C : VariableChange F) (S : Finset (F × F)) :
    (S.image (vcInvPair C)).image (fun Q => (vcX C Q.1, vcY C Q.1 Q.2)) = S := by
  rw [Finset.image_image,
    show ((fun Q : F × F => (vcX C Q.1, vcY C Q.1 Q.2)) ∘ vcInvPair C) = id from
      funext (vcPair_vcInvPair C)]
  exact Finset.image_id

/-- ★★★★★★★★`B` は逆向きの変換で `u³` 倍になる。 -/
theorem veluB_vcInvPair (C : VariableChange F) (W : WeierstrassCurve F) (Q : F × F) :
    2 * (vcInvPair C Q).2 + W.a₁ * (vcInvPair C Q).1 + W.a₃
      = ((C.u : Fˣ) : F) ^ 3
        * (2 * Q.2 + (C • W).a₁ * Q.1 + (C • W).a₃) := by
  have hu : ((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) = 1 := C.u.mul_inv
  simp only [vcInvPair, WeierstrassCurve.variableChange_a₁,
    WeierstrassCurve.variableChange_a₃]
  linear_combination (-((W.a₁ + 2 * C.s) * Q.1 * ((C.u : Fˣ) : F) ^ 2)
    - (W.a₃ + C.r * W.a₁ + 2 * C.t)
      * ((((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F)) ^ 2
        + ((C.u : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) + 1)) * hu

end Pullback

def image_vcPair_vcInvPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(点集合は引き戻せる——L 上の E/H を作るため。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluB_vcInvPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(B は逆向きの変換で u³ 倍になる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def equation_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Equation は変数変換で保たれる——mathlib に無い。★無条件)",
    sectionId := "genell-lemma-3-5" }

def variableChange_equation_poly.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換の方程式の側の恒等式 F_{C•W} = u⁻⁶·F_W。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
