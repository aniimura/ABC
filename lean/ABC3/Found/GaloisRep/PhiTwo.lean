import ABC3.Found.GaloisRep.PhiDegree

/-!
# Galois (G5) 第 199 ブロック —— **★★★★★★★`x(2P) = Φ₂(x)/ΨSq₂(x)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★分点多項式と群法則が実際に一致することを確かめた

残る 1 本は **`Φ_n(x) = x_n · ΨSq_n(x)`** である(第 198)。★その **`n = 2` の場合**を証明した。
★★これは帰納法の底であると同時に、**mathlib の `Φ`・`ΨSq` の定義が我々の群法則
(`addX`・`slope`)と本当に噛み合うことの確認**でもある。

### ★★★★★機構

`ℓ = (3x² + 2a₂x + a₄ − a₁y)/d`(`d = 2y + a₁x + a₃`)として

    x(2P) = ℓ² + a₁ℓ − a₂ − 2x

★まず **`ΨSq₂(x) = d²`**(Weierstrass 方程式を `linear_combination -4·h` で使う)。
★★次に分母を払って **`linear_combination −(a₁² + 4a₂ + 8x)·h`**。
★★★係数 `−(a₁² + 4a₂ + 8x)` は手で出した:`N·(N + a₁d) = A² + AB − a₁²(y²+a₁xy+a₃y)` と
`d² = ΨSq₂ + 4E` から `E` の係数を集めたものである。

### ★Lean 上の注意

`field_simp` は `d = 2y + a₁x + a₃` を**展開してしまう**と `d ≠ 0` と噛み合わなくなる。
★`obtain ⟨d, hdd⟩ : ∃ d, d = …` で **`d` を不透明な変数にしてから** `field_simp` し、
そのあと `subst` して `linear_combination` する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `addX_self_eq_Φ_div` | ★★★★★★**2 倍公式 = `Φ₂/ΨSq₂`** |
| `two_nsmul_some` | ★★★★★★★**`2•P` の座標** |
-/

namespace ABC3.Found.GaloisRep

open Polynomial WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)

/-- ★★★★★★**2 倍公式が mathlib の分点多項式と一致する**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ΨSq₂(x) = (2y + a₁x + a₃)²` を出してから分母を払う。 -/
theorem addX_self_eq_Φ_div {x y : F} (h : W.Equation x y) (hd : W.negY x y ≠ y) :
    W.addX x x (W.slope x x y y) = (W.Φ 2).eval x / (W.ΨSq 2).eval x := by
  have hdne : 2 * y + W.a₁ * x + W.a₃ ≠ 0 := by
    simp only [negY] at hd
    intro h0
    exact hd (by linear_combination -h0)
  rw [equation_iff] at h
  have hΨ : (W.ΨSq 2).eval x = (2 * y + W.a₁ * x + W.a₃) ^ 2 := by
    rw [WeierstrassCurve.ΨSq_two]
    simp only [WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
      WeierstrassCurve.b₆, eval_add, eval_mul, eval_pow, eval_C, eval_X]
    linear_combination -4 * h
  rw [hΨ, WeierstrassCurve.Φ_two, slope, if_pos rfl, if_neg (fun hc => hd hc.symm)]
  simp only [negY, addX, WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈,
    eval_sub, eval_mul, eval_pow, eval_C, eval_X]
  have hd2 : y - (-y - W.a₁ * x - W.a₃) = 2 * y + W.a₁ * x + W.a₃ := by ring
  rw [hd2]
  obtain ⟨d, hdd⟩ : ∃ d : F, d = 2 * y + W.a₁ * x + W.a₃ := ⟨_, rfl⟩
  rw [← hdd]
  rw [← hdd] at hdne
  field_simp
  subst hdd
  linear_combination (-(W.a₁ ^ 2 + 4 * W.a₂ + 8 * x)) * h

/-- ★★★★★★★**`x(2P) = Φ₂(x)/ΨSq₂(x)`**——点の言葉で。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem two_nsmul_some {x y : F} (h : W.Nonsingular x y) (hd : W.negY x y ≠ y) :
    ∃ hns : W.Nonsingular ((W.Φ 2).eval x / (W.ΨSq 2).eval x)
      (W.addY x x y (W.slope x x y y)),
      (2 : ℕ) • Point.some x y h
        = Point.some ((W.Φ 2).eval x / (W.ΨSq 2).eval x)
          (W.addY x x y (W.slope x x y y)) hns := by
  have hxy : ¬(x = x ∧ y = W.negY x y) := fun hc => hd hc.2.symm
  have hadd : Point.some x y h + Point.some x y h
      = Point.some (W.addX x x (W.slope x x y y)) (W.addY x x y (W.slope x x y y))
        (nonsingular_add h h hxy) := Point.add_some hxy
  refine ⟨?_, ?_⟩
  · rw [← addX_self_eq_Φ_div W h.1 hd]
    exact nonsingular_add h h hxy
  · rw [two_nsmul, hadd]
    exact point_some_congr (addX_self_eq_Φ_div W h.1 hd) rfl _ _

/-! ## ★出典の紐付け(`.src`) -/

def addX_self_eq_Φ_div.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——2 倍公式と分点多項式)",
    sectionId := "genell-thm-3-8" }

def two_nsmul_some.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——x(2P) = Φ₂/ΨSq₂)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
