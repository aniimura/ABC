import ABC3.Found.GaloisRep.TatePair

/-!
# Galois (G6) 第 212 ブロック —— **★★★★★★★Weierstrass 方程式が `I` を法として成り立つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★葉 (b) に入った

§9-413 が挙げた (G6) の残りの葉は
**(b) Weierstrass 方程式・(c) 準同型性・(d) 核・(e) 全射性**の 4 つだった。
★本ブロックは (b) の**最下段**——`I` を法とした場合——を取る。

### ★★★★★★特異 Tate 曲線の等式は有理式の恒等式だった

`q ≡ 0` では `w = q/u ≡ 0` なので `X ≡ f(a)`・`Y ≡ g(a)` になり、`a₄ ≡ a₆ ≡ 0` である。
★残るのは **`g(t)² + f(t)g(t) = f(t)³`**——すなわち節点付き三次曲線 `y² + xy = x³` の
有理パラメータ表示である。★★`f(t) = t/(1−t)²`・`g(t) = t²/(1−t)³` を入れると

    (1−t)⁶ g² = t⁴,  (1−t)⁶ f g = t³(1−t),  (1−t)⁶ f³ = t³

なので `t⁴ + t³(1−t) = t³` で閉じる。★★★`1 − t` の逆元 `v` を取って
`linear_combination (−t³v⁵) * ((1−t)v = 1)` の 1 行になった。

### ★★★合同の段

| 段 | 内容 |
|---|---|
| `evalAdic_mem` | 定数項 0 の冪級数の `I` 進値は `I` に入る(`evalAdic_spec` を `n = 1` で使う) |
| `tateCurveAt_a4_mem` / `_a6_mem` | ★`a₄(q), a₆(q) ∈ I` |
| `tateXpair_sub_mem` / `tateYpair_sub_mem` | ★★`X − f(a) ∈ I`・`Y − g(a) ∈ I` |
| `tate_equation_mem` | ★★★★★★★**方程式が `I` を法として成り立つ** |

★差の分解 `Y² − g² = (Y−g)(Y+g)`・`XY − fg = (X−f)Y + f(Y−g)`・
`X³ − f³ = (X−f)(X²+Xf+f²)` を `ring` で作り、各項が `I` に入ることを言うだけである。

## ★★残っているもの

`I` を法とした先——`I^n` を法とした帰納、すなわち**厳密な等式**——には
第 112 ブロックの `I` 進 Cauchy 積(`adicSum_mul`)が要る。
★(c) 準同型性・(d) 核・(e) 全射性は手つかずである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_singular_equation` | ★★★★★★**`g(t)² + f(t)g(t) = f(t)³`** |
| `evalAdic_mem` | ★★★定数項 0 なら値は `I` に入る |
| `tateCurveAt_a4_mem` / `_a6_mem` | ★`a₄(q), a₆(q) ∈ I` |
| `tateXpair_sub_mem` / `tateYpair_sub_mem` | ★★`X ≡ f(a)`・`Y ≡ g(a)` |
| `tate_equation_mem` | ★★★★★★★**Weierstrass 方程式(`I` 法)** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★★特異 Tate 曲線の等式 -/

/-- ★★★★★★**特異 Tate 曲線の方程式**——`g(t)² + f(t)g(t) = f(t)³`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★節点付き三次曲線 `y² + xy = x³` の有理パラメータ表示である。 -/
theorem tate_singular_equation {t : R} (hu : IsUnit (1 - t)) :
    tateYterm t ^ 2 + tateXterm t * tateYterm t = tateXterm t ^ 3 := by
  have hv : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  set v := Ring.inverse (1 - t) with hvdef
  show (t ^ 2 * v ^ 3) ^ 2 + (t * v ^ 2) * (t ^ 2 * v ^ 3) = (t * v ^ 2) ^ 3
  linear_combination (-(t ^ 3 * v ^ 5)) * hv

/-! ## ★★★合同の段 -/

/-- ★★★定数項が 0 の冪級数の `I` 進値は `I` に入る。 -/
theorem evalAdic_mem [IsAdicComplete I R] (f : PowerSeries ℤ) (hf : PowerSeries.coeff 0 f = 0)
    (q : R) (hq : q ∈ I) : evalAdic f q hq ∈ I := by
  have h := evalAdic_spec (I := I) f q hq 1
  have h0 : partialEval f q 1 = 0 := by
    simp [partialEval, hf]
  rw [h0] at h
  have hmem : evalAdic f q hq ∈ (I ^ 1 • (⊤ : Submodule R R)) := SModEq.zero.1 h.symm
  simpa [Ideal.smul_top_eq_map] using hmem

theorem tateCurveAt_a4_mem [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₄ ∈ I := by
  have h : (tateCurveAt (I := I) q hq).a₄ = evalAdic tateA4 q hq := by
    simp [tateCurveAt, tateCurve, WeierstrassCurve.map, evalAdicHom]
  rw [h]
  exact evalAdic_mem _ coeff_zero_tateA4 q hq

theorem tateCurveAt_a6_mem [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₆ ∈ I := by
  have h : (tateCurveAt (I := I) q hq).a₆ = evalAdic tateA6 q hq := by
    simp [tateCurveAt, tateCurve, WeierstrassCurve.map, evalAdicHom]
  rw [h]
  exact evalAdic_mem _ coeff_zero_tateA6 q hq

/-- ★★`X(u,q) ≡ f(a)`(`I` 法)。 -/
theorem tateXpair_sub_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I) :
    tateXpair a w q hq - tateXterm a ∈ I := by
  have h1 : tateXtail (I := I) a q hq ∈ I := tateXtail_mem a q hq
  have h2 : tateXterm w ∈ I := by
    simpa using tateXterm_mem_pow (I := I) (k := 1) (by simpa using hw)
  have h3 : tateXtail (I := I) w q hq ∈ I := tateXtail_mem w q hq
  have h4 : evalAdic (sigmaSeries 1) q hq ∈ I := evalAdic_mem _ (coeff_zero_sigmaSeries 1) q hq
  have hkey : tateXpair a w q hq - tateXterm a
      = tateXtail a q hq + tateXterm w + tateXtail w q hq
        - 2 * evalAdic (sigmaSeries 1) q hq := by
    simp only [tateXpair]; ring
  rw [hkey]
  exact I.sub_mem (I.add_mem (I.add_mem h1 h2) h3) (I.mul_mem_left _ h4)

/-- ★★`Y(u,q) ≡ g(a)`(`I` 法)。 -/
theorem tateYpair_sub_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I) :
    tateYpair a w q hq - tateYterm a ∈ I := by
  have h1 : tateYtail (I := I) a q hq ∈ I := tateYtail_mem a q hq
  have h2 : tateXterm w ∈ I := by
    simpa using tateXterm_mem_pow (I := I) (k := 1) (by simpa using hw)
  have h3 : tateXtail (I := I) w q hq ∈ I := tateXtail_mem w q hq
  have h4 : evalAdic (sigmaSeries 1) q hq ∈ I := evalAdic_mem _ (coeff_zero_sigmaSeries 1) q hq
  have h5 : tateYterm w ∈ I := by
    simpa using tateYterm_mem_pow (I := I) (k := 1) (by simpa using hw)
  have h6 : tateYtail (I := I) w q hq ∈ I := tateYtail_mem w q hq
  have hkey : tateYpair a w q hq - tateYterm a
      = tateYtail a q hq - (tateXterm w + tateXtail w q hq)
        - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq := by
    simp only [tateYpair]; ring
  rw [hkey]
  exact I.add_mem (I.sub_mem (I.sub_mem h1 (I.add_mem h2 h3)) (I.add_mem h5 h6)) h4

/-! ## ★★★★★★★方程式(`I` 法) -/

/-- ★★★★★★★**Tate 級数は Weierstrass 方程式を `I` を法として満たす**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★Tate 曲線は `a₁ = 1`, `a₂ = a₃ = 0` なので方程式は `y² + xy = x³ + a₄x + a₆` である。
★★`I` を法とすると `a₄ ≡ a₆ ≡ 0`・`X ≡ f(a)`・`Y ≡ g(a)` になり、
`tate_singular_equation` で閉じる。 -/
theorem tate_equation_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I)
    (hu : IsUnit (1 - a)) :
    (tateYpair a w q hq) ^ 2 + tateXpair a w q hq * tateYpair a w q hq
      - ((tateXpair a w q hq) ^ 3 + (tateCurveAt q hq).a₄ * tateXpair a w q hq
        + (tateCurveAt q hq).a₆) ∈ I := by
  have hX := tateXpair_sub_mem (I := I) a w q hq hw
  have hY := tateYpair_sub_mem (I := I) a w q hq hw
  have ha4 := tateCurveAt_a4_mem (I := I) q hq
  have ha6 := tateCurveAt_a6_mem (I := I) q hq
  have hsing : tateYterm a ^ 2 + tateXterm a * tateYterm a - tateXterm a ^ 3 = 0 := by
    rw [tate_singular_equation hu]; ring
  set X := tateXpair a w q hq with hXd
  set Y := tateYpair a w q hq with hYd
  set f := tateXterm a with hfd
  set g := tateYterm a with hgd
  set A4 := (tateCurveAt (I := I) q hq).a₄ with hA4
  set A6 := (tateCurveAt (I := I) q hq).a₆ with hA6
  have hkey : Y ^ 2 + X * Y - (X ^ 3 + A4 * X + A6)
      = (Y - g) * (Y + g) + (X - f) * Y + f * (Y - g)
        - (X - f) * (X ^ 2 + X * f + f ^ 2) - A4 * X - A6
        + (g ^ 2 + f * g - f ^ 3) := by ring
  rw [hkey, hsing, add_zero]
  exact I.sub_mem (I.sub_mem (I.sub_mem
    (I.add_mem (I.add_mem (I.mul_mem_right _ hY) (I.mul_mem_right _ hX)) (I.mul_mem_left _ hY))
    (I.mul_mem_right _ hX)) (I.mul_mem_right _ ha4)) ha6

/-! ## ★出典の紐付け(`.src`) -/

def tate_singular_equation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——特異 Tate 曲線の有理パラメータ表示)",
    sectionId := "genell-def-3-3" }

def tate_equation_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Weierstrass 方程式が I を法として成り立つこと)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
