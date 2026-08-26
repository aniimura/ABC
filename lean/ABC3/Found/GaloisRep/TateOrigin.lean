import ABC3.Found.GaloisRep.TateOriginUniv

/-!
# Galois (G6) 第 275 ブロック —— **★★★★★★★★★★原点近傍でも方程式が成り立つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点——`1 − a` の単元性が要らない方程式

    YE² + XE·YE·(1−a) = XE³ + a₄·XE·(1−a)⁴ + a₆·(1−a)⁶      (`tate_equationE`)

ここで `XE = (1−a)²·X`、`YE = (1−a)³·Y` である(`Ring.inverse (1−a)` を含まない)。
★★★**`u ≡ 1 mod 𝔪` でもそのまま成り立つ**。

## ★★★★★★★★`K` を経由しなくてよかった

原点近傍は「`X` が `K` に極を持つから `K` の水準で定式化が要る」と見ていたが、
**分母を先に払えば `R` の中の等式で済む**:

| | 極を持つ形 | 分母を払った形 |
|---|---|---|
| 座標 | `X = XE/(1−a)²` | `XE ∈ R` |
| 方程式 | `Y² + XY = X³ + a₄X + a₆` | `YE² + XE·YE·(1−a) = …` |
| 住む場所 | `K` | **`R`** |

★★★★`(1−a)⁶ × (方程式の差)` が `I^n` に入る(第 274)——これを `n → ∞` すれば
`IsHausdorff` で 0。**付値も極限も要らない**。`K` は最後に商を取るだけの場所になった。

## ★★★★★★収束は尾だけ

`XE − XE_n = (1−a)²·(尾の残り)` なので、`(1−a)²` は**外に出たまま**である。
★したがって残りは `I^n` の元の定数倍で、`n` に依らない有界性が自動的に効く。

## ★★★★★★`K` の水準の座標

    X_K := XE / (1−a)²、  Y_K := YE / (1−a)³        (`algebraMap R K` の像で割る)

★`1 − a` が単元なら `X_K = algebraMap (X)` に一致する(`tateXK_eq`)ので、
これまでの領域の議論と**接ぎ木できる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXpairE`・`tateYpairE` | ★★★★★★分母を払った座標 |
| `tateDefectE_eq_zero` | ★★★★★★★★★★**原点近傍でも差は 0** |
| `tate_equationE` | ★★★★★★★★★★**原点近傍でも Weierstrass 方程式** |
| `tateXK`・`tateYK` | ★★★★★★`K` の水準の座標 |
| `tateXK_eq`・`tateYK_eq` | ★★★★既存の領域との接ぎ木 |
| `tate_equation_mapK`・`nonsingular_tateK` | ★★★★★★★★★`K` の点は曲線の上にある |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★★★★★分母を払った座標 -/

section Pair

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★**`(1−a)²` を掛けた `X`**——`Ring.inverse (1−a)` を含まない。 -/
noncomputable def tateXpairE [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  a + (1 - a) ^ 2 * (tateXtail a q hq + (tateXterm w + tateXtail w q hq)
    - 2 * evalAdic (sigmaSeries 1) q hq)

/-- ★★★★★★**`(1−a)³` を掛けた `Y`**。 -/
noncomputable def tateYpairE [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  a ^ 2 + (1 - a) ^ 3 * (tateYtail a q hq - (tateXterm w + tateXtail w q hq)
    - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq)

theorem tateXpairE_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) :
    tateXpairE a w q hq = (1 - a) ^ 2 * tateXpair a w q hq := by
  rw [tateXpairE, tateXpair]
  linear_combination -mul_tateXterm' ha

theorem tateYpairE_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) :
    tateYpairE a w q hq = (1 - a) ^ 3 * tateYpair a w q hq := by
  rw [tateYpairE, tateYpair]
  linear_combination -mul_tateYterm' ha

/-- ★★★★★★収束は尾だけ——`(1−a)²` は外に出たまま。 -/
theorem tateXpairE_sub_trunc [IsAdicComplete I R] (n : ℕ) (a w q : R) (hq : q ∈ I) :
    tateXpairE a w q hq - tateXtruncE n a w q ∈ I ^ n := by
  have hkey : tateXpairE a w q hq - tateXtruncE n a w q
      = (1 - a) ^ 2 * ((tateXtail a q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * a)) n)
        + (tateXtail w q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
        - 2 * (evalAdic (sigmaSeries 1) q hq - partialEval (sigmaSeries 1) q n)) := by
    rw [tateXpairE, tateXtruncE]; ring
  rw [hkey]
  exact Ideal.mul_mem_left _ _ (Ideal.sub_mem _ (Ideal.add_mem _
    (tateXtail_sub_partialSum n a q hq) (tateXtail_sub_partialSum n w q hq))
    (Ideal.mul_mem_left _ _ (evalAdic_sub_partialEval n (sigmaSeries 1) q hq)))

theorem tateYpairE_sub_trunc [IsAdicComplete I R] (n : ℕ) (a w q : R) (hq : q ∈ I) :
    tateYpairE a w q hq - tateYtruncE n a w q ∈ I ^ n := by
  have hkey : tateYpairE a w q hq - tateYtruncE n a w q
      = (1 - a) ^ 3 * ((tateYtail a q hq - partialSum (fun m => tateYterm (q ^ (m + 1) * a)) n)
        - (tateXtail w q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
        - (tateYtail w q hq - partialSum (fun m => tateYterm (q ^ (m + 1) * w)) n)
        + (evalAdic (sigmaSeries 1) q hq - partialEval (sigmaSeries 1) q n)) := by
    rw [tateYpairE, tateYtruncE]; ring
  rw [hkey]
  exact Ideal.mul_mem_left _ _ (Ideal.add_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
    (tateYtail_sub_partialSum n a q hq) (tateXtail_sub_partialSum n w q hq))
    (tateYtail_sub_partialSum n w q hq)) (evalAdic_sub_partialEval n (sigmaSeries 1) q hq))

/-! ## ★★★★★★★★★★原点近傍でも方程式が成り立つ -/

/-- ★★★★★★**分母を払った方程式の差**。 -/
noncomputable def tateDefectE [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  tateYpairE a w q hq ^ 2 + tateXpairE a w q hq * tateYpairE a w q hq * (1 - a)
    - (tateXpairE a w q hq ^ 3 + (tateCurveAt q hq).a₄ * tateXpairE a w q hq * (1 - a) ^ 4
      + (tateCurveAt q hq).a₆ * (1 - a) ^ 6)

theorem tateDefectE_sub_trunc [IsAdicComplete I R] (n : ℕ) (a w q : R) (hq : q ∈ I) :
    tateDefectE a w q hq - tateDefectTruncE n a w q ∈ I ^ n := by
  refine Ideal.Quotient.eq.1 ?_
  have hX : (Ideal.Quotient.mk (I ^ n)) (tateXpairE a w q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtruncE n a w q) :=
    Ideal.Quotient.eq.2 (tateXpairE_sub_trunc n a w q hq)
  have hY : (Ideal.Quotient.mk (I ^ n)) (tateYpairE a w q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtruncE n a w q) :=
    Ideal.Quotient.eq.2 (tateYpairE_sub_trunc n a w q hq)
  have h4 : (Ideal.Quotient.mk (I ^ n)) ((tateCurveAt q hq).a₄)
      = (Ideal.Quotient.mk (I ^ n)) (partialEval tateA4 q n) :=
    Ideal.Quotient.eq.2 (evalAdic_sub_partialEval n tateA4 q hq)
  have h6 : (Ideal.Quotient.mk (I ^ n)) ((tateCurveAt q hq).a₆)
      = (Ideal.Quotient.mk (I ^ n)) (partialEval tateA6 q n) :=
    Ideal.Quotient.eq.2 (evalAdic_sub_partialEval n tateA6 q hq)
  simp only [tateDefectE, tateDefectTruncE, map_sub, map_add, map_mul, map_pow, map_one,
    hX, hY, h4, h6]

/-- ★★★★★★★★★★**原点近傍でも方程式の差は 0**——`1 − a` の単元性が要らない。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateDefectE_eq_zero [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hw : IsUnit (1 - w)) : tateDefectE a w q hq = 0 := by
  refine eq_zero_of_mem_pow (I := I) fun n => ?_
  have h1 := tateDefectE_sub_trunc (I := I) n a w q hq
  have h2 := tateDefectTruncE_mem (I := I) n a w q hq haw hw
  have h3 := Ideal.add_mem (I ^ n) h1 h2
  simpa using h3

/-- ★★★★★★★★★★**原点近傍でも Weierstrass 方程式が成り立つ**(分母を払った形)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_equationE [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hw : IsUnit (1 - w)) :
    tateYpairE a w q hq ^ 2 + tateXpairE a w q hq * tateYpairE a w q hq * (1 - a)
      = tateXpairE a w q hq ^ 3 + (tateCurveAt q hq).a₄ * tateXpairE a w q hq * (1 - a) ^ 4
        + (tateCurveAt q hq).a₆ * (1 - a) ^ 6 := by
  have h := tateDefectE_eq_zero a w q hq haw hw
  rw [tateDefectE] at h
  linear_combination h

end Pair

/-! ## ★★★★★★`K` の水準の座標 -/

section Field

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K] [Algebra R K]

/-- ★★★★★★**`K` の水準の `X`**——原点近傍では極を持つ。 -/
noncomputable def tateXK (a w q : R) (hq : q ∈ I) : K :=
  algebraMap R K (tateXpairE a w q hq) * (algebraMap R K (1 - a))⁻¹ ^ 2

/-- ★★★★★★**`K` の水準の `Y`**。 -/
noncomputable def tateYK (a w q : R) (hq : q ∈ I) : K :=
  algebraMap R K (tateYpairE a w q hq) * (algebraMap R K (1 - a))⁻¹ ^ 3

theorem tateXK_eq (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) :
    (tateXK a w q hq : K) = algebraMap R K (tateXpair a w q hq) := by
  have hne : algebraMap R K (1 - a) ≠ 0 := (ha.map (algebraMap R K)).ne_zero
  rw [tateXK, tateXpairE_eq a w q hq ha, map_mul, map_pow]
  field_simp

theorem tateYK_eq (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) :
    (tateYK a w q hq : K) = algebraMap R K (tateYpair a w q hq) := by
  have hne : algebraMap R K (1 - a) ≠ 0 := (ha.map (algebraMap R K)).ne_zero
  rw [tateYK, tateYpairE_eq a w q hq ha, map_mul, map_pow]
  field_simp

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★★★**原点近傍でも `K` の点は曲線の上にある**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_equation_mapK (a w q : R) (hq : q ∈ I) (haw : a * w = q) (hw : IsUnit (1 - w))
    (ha : algebraMap R K (1 - a) ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Equation
      (tateXK a w q hq) (tateYK (K := K) a w q hq) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  have h := congrArg (algebraMap R K) (tate_equationE a w q hq haw hw)
  simp only [map_add, map_mul, map_pow] at h
  have e1 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt q hq).a₁) = 1
    rw [show (tateCurveAt q hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e2 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₂ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₂) = 0
    rw [show (tateCurveAt q hq).a₂ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have e3 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₃) = 0
    rw [show (tateCurveAt q hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have e4 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₄
      = algebraMap R K ((tateCurveAt q hq).a₄) := rfl
  have e6 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₆
      = algebraMap R K ((tateCurveAt q hq).a₆) := rfl
  rw [e1, e2, e3, e4, e6, tateXK, tateYK]
  field_simp
  linear_combination (algebraMap R K (1 - a)) ^ 0 * h

theorem nonsingular_tateK (a w q : R) (hq : q ∈ I) (haw : a * w = q) (hw : IsUnit (1 - w))
    (ha : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (tateXK a w q hq) (tateYK (K := K) a w q hq) :=
  (WeierstrassCurve.Affine.equation_iff_nonsingular_of_Δ_ne_zero hΔ).mp
    (tate_equation_mapK a w q hq haw hw ha)

end Field

/-! ## ★出典の紐付け(`.src`) -/

def tate_equationE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍でも Weierstrass 方程式)",
    sectionId := "genell-def-3-3" }

def tate_equation_mapK.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍でも K の点は曲線の上)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
