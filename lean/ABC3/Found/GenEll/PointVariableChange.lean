/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.EllipticCurve.VariableChange
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Formula
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

def equation_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Equation は変数変換で保たれる——mathlib に無い。★無条件)",
    sectionId := "genell-lemma-3-5" }

def variableChange_equation_poly.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換の方程式の側の恒等式 F_{C•W} = u⁻⁶·F_W。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
