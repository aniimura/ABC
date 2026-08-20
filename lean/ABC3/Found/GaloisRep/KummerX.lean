import ABC3.Found.GaloisRep.PhiTwo

/-!
# Galois (G5) 第 200 ブロック —— **★★★★★★★Kummer の 2 公式(`y` が消える)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★`ω_n` を定義せずに済ませる道

残る 1 本 `Φ_n(x) = x_n·ΨSq_n(x)` を一般の `n` へ進めるには、古典的には `y` 側の多項式
`ω_n` が要る。★**mathlib はまだ `ω_n` を持たない**(docstring に `TODO`)。
★★そこで **`x` だけの(Kummer)公式**で帰納法を回す道を測った。

### ★★★★★★2 つの対称式は `y` が消える

`x₁ ≠ x₂` の 2 点 `A = (x₁,y₁)`・`B = (x₂,y₂)` について、
`x(A+B)` と `x(A−B)` の**和と積**は `y` を含まない:

    (x(A+B) + x(A−B))·(x₁−x₂)² = 2x₁x₂(x₁+x₂) + b₂x₁x₂ + b₄(x₁+x₂) + b₆
    (x(A+B) · x(A−B))·(x₁−x₂)² = x₁²x₂² − b₄x₁x₂ − b₆(x₁+x₂) − b₈

★★これで Montgomery 梯子が回る——`x(A−B)` が既知なら `x(A+B)` が積の式から出る。

### ★★★★★積の式は「差の式」経由で出した

和の式は `linear_combination 2h₁ + 2h₂` で 1 行。★積の式を直接やると
`y` の 4 次になって係数が見つけにくい。★★そこで**差**を先に出した:

    (x(A+B) − x(A−B))·(x₁−x₂)² = −(2y₁+a₁x₁+a₃)(2y₂+a₁x₂+a₃)

これは `ring` で出る。★★★二乗すると第 199 の `ΨSq₂(x) = (2y+a₁x+a₃)²` で
**`y` が消える**。★★★★あとは `4uv = (u+v)² − (u−v)²` と、`b` だけの多項式恒等式

    M² − ΨSq₂(x₁)·ΨSq₂(x₂) = 4N(x₁−x₂)²

(`ring` で出る)を合わせるだけ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `ΨSq_two_eval` | `ΨSq₂(x) = (2y+a₁x+a₃)²`(曲線上で) |
| `eval_ΨSq_two` | `ΨSq₂(x) = 4x³+b₂x²+2b₄x+b₆` |
| `addX_sum` | ★★★★★★**Kummer の和の公式** |
| `addX_sub` | ★★★★★差の公式(`y` は残るが二乗で消える) |
| `addX_prod` | ★★★★★★★**Kummer の積の公式** |
-/

namespace ABC3.Found.GaloisRep

open Polynomial WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)

/-- 曲線上では `ΨSq₂(x) = (2y + a₁x + a₃)²`。 -/
theorem ΨSq_two_eval {x y : F} (h : W.Equation x y) :
    (W.ΨSq 2).eval x = (2 * y + W.a₁ * x + W.a₃) ^ 2 := by
  rw [equation_iff] at h
  rw [WeierstrassCurve.ΨSq_two]
  simp only [WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆, eval_add, eval_mul, eval_pow, eval_C, eval_X]
  linear_combination -4 * h

theorem eval_ΨSq_two (x : F) :
    (W.ΨSq 2).eval x = 4 * x ^ 3 + W.b₂ * x ^ 2 + 2 * W.b₄ * x + W.b₆ := by
  rw [WeierstrassCurve.ΨSq_two]
  simp only [WeierstrassCurve.Ψ₂Sq, eval_add, eval_mul, eval_pow, eval_C, eval_X]

/-- ★★★★★★**Kummer の和の公式**——`y` が消える。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem addX_sum {x₁ y₁ x₂ y₂ : F} (h₁ : W.Equation x₁ y₁) (h₂ : W.Equation x₂ y₂)
    (hx : x₁ ≠ x₂) :
    (W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂) + W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂)))
        * (x₁ - x₂) ^ 2
      = 2 * x₁ * x₂ * (x₁ + x₂) + W.b₂ * x₁ * x₂ + W.b₄ * (x₁ + x₂) + W.b₆ := by
  have hne : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hx
  rw [equation_iff] at h₁ h₂
  rw [slope_of_X_ne hx, slope_of_X_ne hx]
  simp only [negY, addX, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆]
  obtain ⟨d, hdd⟩ : ∃ d : F, d = x₁ - x₂ := ⟨_, rfl⟩
  rw [← hdd]
  rw [← hdd] at hne
  field_simp
  subst hdd
  linear_combination 2 * h₁ + 2 * h₂

/-- ★★★★★差の公式——`y` は残るが、二乗すると第 199 で消える。 -/
theorem addX_sub {x₁ y₁ x₂ y₂ : F} (hx : x₁ ≠ x₂) :
    (W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂) - W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂)))
        * (x₁ - x₂) ^ 2
      = -((2 * y₁ + W.a₁ * x₁ + W.a₃) * (2 * y₂ + W.a₁ * x₂ + W.a₃)) := by
  have hne : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hx
  rw [slope_of_X_ne hx, slope_of_X_ne hx]
  simp only [negY, addX]
  obtain ⟨d, hdd⟩ : ∃ d : F, d = x₁ - x₂ := ⟨_, rfl⟩
  rw [← hdd]
  rw [← hdd] at hne
  field_simp
  subst hdd
  ring

/-- ★★★★★★★**Kummer の積の公式**——`y` が消える。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`4uv = (u+v)² − (u−v)²` と差の公式の二乗、そして `b` だけの多項式恒等式から。 -/
theorem addX_prod (h4 : (4 : F) ≠ 0) {x₁ y₁ x₂ y₂ : F}
    (h₁ : W.Equation x₁ y₁) (h₂ : W.Equation x₂ y₂) (hx : x₁ ≠ x₂) :
    (W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂) * W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂)))
        * (x₁ - x₂) ^ 2
      = x₁ ^ 2 * x₂ ^ 2 - W.b₄ * (x₁ * x₂) - W.b₆ * (x₁ + x₂) - W.b₈ := by
  have hne : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hx
  have hsum := addX_sum W h₁ h₂ hx
  have hdiff := addX_sub W (y₁ := y₁) (y₂ := y₂) hx
  have hΨ1 := ΨSq_two_eval W h₁
  have hΨ2 := ΨSq_two_eval W h₂
  rw [eval_ΨSq_two] at hΨ1 hΨ2
  have hpoly : (2 * x₁ * x₂ * (x₁ + x₂) + W.b₂ * x₁ * x₂ + W.b₄ * (x₁ + x₂) + W.b₆) ^ 2
      - (4 * x₁ ^ 3 + W.b₂ * x₁ ^ 2 + 2 * W.b₄ * x₁ + W.b₆)
        * (4 * x₂ ^ 3 + W.b₂ * x₂ ^ 2 + 2 * W.b₄ * x₂ + W.b₆)
      = 4 * (x₁ ^ 2 * x₂ ^ 2 - W.b₄ * (x₁ * x₂) - W.b₆ * (x₁ + x₂) - W.b₈) * (x₁ - x₂) ^ 2 := by
    simp only [WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
      WeierstrassCurve.b₈]
    ring
  have hkey : ((W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂)
        * W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂))) * (x₁ - x₂) ^ 2)
        * (4 * (x₁ - x₂) ^ 2)
      = (x₁ ^ 2 * x₂ ^ 2 - W.b₄ * (x₁ * x₂) - W.b₆ * (x₁ + x₂) - W.b₈)
        * (4 * (x₁ - x₂) ^ 2) := by
    have hh : (4 : F) * ((W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂)
        * W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂))) * (x₁ - x₂) ^ 2) * (x₁ - x₂) ^ 2
      = 4 * (x₁ ^ 2 * x₂ ^ 2 - W.b₄ * (x₁ * x₂) - W.b₆ * (x₁ + x₂) - W.b₈)
        * (x₁ - x₂) ^ 2 := by
      rw [← hpoly]
      linear_combination
        ((W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂)
          + W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂))) * (x₁ - x₂) ^ 2
          + (2 * x₁ * x₂ * (x₁ + x₂) + W.b₂ * x₁ * x₂ + W.b₄ * (x₁ + x₂) + W.b₆)) * hsum
        - ((W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂)
          - W.addX x₁ x₂ (W.slope x₁ x₂ y₁ (W.negY x₂ y₂))) * (x₁ - x₂) ^ 2
          - (2 * y₁ + W.a₁ * x₁ + W.a₃) * (2 * y₂ + W.a₁ * x₂ + W.a₃)) * hdiff
        + (4 * x₂ ^ 3 + W.b₂ * x₂ ^ 2 + 2 * W.b₄ * x₂ + W.b₆) * hΨ1
        + (2 * y₁ + W.a₁ * x₁ + W.a₃) ^ 2 * hΨ2
    linear_combination hh
  exact mul_right_cancel₀ (mul_ne_zero h4 (pow_ne_zero 2 hne)) hkey

/-! ## ★出典の紐付け(`.src`) -/

def addX_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——Kummer の和の公式)",
    sectionId := "genell-thm-3-8" }

def addX_prod.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——Kummer の積の公式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
