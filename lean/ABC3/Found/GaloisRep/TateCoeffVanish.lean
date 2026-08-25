import ABC3.Found.GaloisRep.TateUWFree

/-!
# Galois (G6) 第 238 ブロック —— **★★★★★★★★分子の低次係数は消える**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★三つの部品を噛み合わせる

| 部品 | 出どころ |
|---|---|
| `tateEval u w P = tateDefectTrunc (n+1) u w (uw) · tateEval u w d` | 第 225 |
| `‖tateDefectTrunc (n+1) u w (uw)‖ ≤ C_n·50·(64‖uw‖)^{n+1}` | 第 237 |
| `‖f(w)‖ ≤ K‖w‖ᵐ` なら `f` の低次係数は 0 | 第 236 |

★あと一つ、**分母 `d` の値が有界**であることが要る。

## ★★★★分母の有界性は帰納法で出る

`‖u‖ ≤ 1`、`‖w‖ ≤ 1` の上では `tateEval u w Q` は有界である。
★`MvPolynomial.induction_on` で:

| 段 | 上界 |
|---|---|
| `C r` | `‖(r : ℂ)‖` |
| `p + q` | `B_p + B_q` |
| `p · X i` | `B_p`(`‖u‖, ‖w‖ ≤ 1` なので変数を掛けても増えない) |

★★**上界を存在量化しておく**のがよい——係数の和として明示すると
`MvPolynomial.support` の扱いが要るが、帰納法なら三段で済む。

## ★★★仕上げ

`‖u‖ ≤ 1` なので `‖uw‖ ≤ ‖w‖`、したがって `(64‖uw‖)^{n+1} ≤ 64^{n+1}‖w‖^{n+1}`。
`‖w‖ < 1/128` なら第 237 の仮定(`‖w‖ ≤ 1/2` と `‖uw‖ ≤ 1/128`)がどちらも満たされる。
よって

    ‖(evalW u P).eval w‖ ≤ (C_n · 50 · 64^{n+1} · B) · ‖w‖^{n+1}

となり、第 236 から **`j ≤ n` の係数が 0** が出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateEval_C` | ★定数の評価 |
| `exists_bound_tateEval` | ★★★★分母は有界 |
| `coeff_evalW_eq_zero` | ★★★★★★★★**分子の低次係数は消える** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★★分母の有界性 -/

theorem tateEval_C (u w : ℂ) (r : ℤ) : tateEval u w (MvPolynomial.C r) = (r : ℂ) := by
  rw [tateEval, MvPolynomial.eval₂Hom_C]
  rfl

/-- ★★★★**`‖u‖, ‖w‖ ≤ 1` の上で `tateEval u w Q` は有界**。 -/
theorem exists_bound_tateEval (Q : TateBase) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ u w : ℂ, ‖u‖ ≤ 1 → ‖w‖ ≤ 1 → ‖tateEval u w Q‖ ≤ B := by
  refine MvPolynomial.induction_on (motive := fun p =>
    ∃ B : ℝ, 0 ≤ B ∧ ∀ u w : ℂ, ‖u‖ ≤ 1 → ‖w‖ ≤ 1 → ‖tateEval u w p‖ ≤ B) Q ?_ ?_ ?_
  · intro r
    refine ⟨‖((r : ℤ) : ℂ)‖, norm_nonneg _, fun u w _ _ => ?_⟩
    rw [tateEval_C]
  · rintro p q ⟨Bp, hBp0, hBp⟩ ⟨Bq, hBq0, hBq⟩
    refine ⟨Bp + Bq, by linarith, fun u w hu hw => ?_⟩
    rw [map_add]
    refine (norm_add_le _ _).trans ?_
    have h1 := hBp u w hu hw
    have h2 := hBq u w hu hw
    linarith
  · rintro p i ⟨Bp, hBp0, hBp⟩
    refine ⟨Bp, hBp0, fun u w hu hw => ?_⟩
    rw [map_mul, norm_mul]
    have hxi : ‖tateEval u w (MvPolynomial.X i)‖ ≤ 1 := by
      fin_cases i
      · show ‖tateEval u w univA‖ ≤ 1
        rw [tateEval_A]
        exact hu
      · show ‖tateEval u w univW‖ ≤ 1
        rw [tateEval_W]
        exact hw
    have hp := hBp u w hu hw
    nlinarith [norm_nonneg (tateEval u w p), norm_nonneg (tateEval u w (MvPolynomial.X i))]

/-! ## ★★★★★★★★分子の低次係数は消える -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**分子の低次係数は消える**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem coeff_evalW_eq_zero (z : UpperHalfPlane) (n : ℕ) (P : TateBase) (d : tateDenoms)
    (hPd : tateDefectTrunc (n + 1) uA uW (uA * uW)
      * algebraMap TateBase TateUniv (d : TateBase) = algebraMap TateBase TateUniv P)
    (j : ℕ) (hj : j < n + 1) :
    (evalW (Complex.exp (2 * ↑π * I * (z : ℂ))) P).coeff j = 0 := by
  set u := Complex.exp (2 * ↑π * I * (z : ℂ)) with hudef
  have hu : ‖u‖ < 1 := UpperHalfPlane.norm_exp_two_pi_I_lt_one z
  have hu1 : ‖u‖ ≤ 1 := hu.le
  obtain ⟨B, hB0, hB⟩ := exists_bound_tateEval (d : TateBase)
  set Cn : ℝ := 3 * (‖tateXterm u‖ + ‖tateYterm u‖ + 30 * ((n : ℝ) + 1) + 80) ^ 2
      + 6 * (‖tateXterm u‖ + ‖tateYterm u‖ + 30 * ((n : ℝ) + 1) + 80) + 1 with hCndef
  have hCn0 : (0 : ℝ) ≤ Cn := by
    simp only [hCndef]
    have h1 := norm_nonneg (tateXterm u)
    have h2 := norm_nonneg (tateYterm u)
    have h3 : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
    positivity
  refine coeff_eq_zero_of_norm_le _ (Cn * 50 * 64 ^ (n + 1) * B) (1 / 128) (by norm_num)
    (n + 1) ?_ j hj
  intro w hw0 hwr
  have hw2 : ‖w‖ ≤ 1 / 2 := by linarith
  have hw1 : ‖w‖ < 1 := by linarith
  have huw : ‖u * w‖ ≤ 1 / 128 := by
    rw [norm_mul]
    nlinarith [norm_nonneg w, norm_nonneg u]
  have hnum := tateEval_numerator (n + 1) P d hPd u w hu hw1
  have hbound := norm_tateDefectTrunc_uw_le z w hw0 hw2 huw n
  have hd := hB u w hu1 hw1.le
  rw [eval_evalW, hnum, norm_mul]
  have hmono : ((64 : ℝ) * ‖u * w‖) ^ (n + 1) ≤ 64 ^ (n + 1) * ‖w‖ ^ (n + 1) := by
    have hle : (64 : ℝ) * ‖u * w‖ ≤ 64 * ‖w‖ := by
      rw [norm_mul]
      nlinarith [norm_nonneg w, norm_nonneg u]
    calc ((64 : ℝ) * ‖u * w‖) ^ (n + 1) ≤ ((64 : ℝ) * ‖w‖) ^ (n + 1) :=
          pow_le_pow_left₀ (by positivity) hle _
      _ = 64 ^ (n + 1) * ‖w‖ ^ (n + 1) := by rw [mul_pow]
  have hdefnn : (0 : ℝ) ≤ ‖tateDefectTrunc (n + 1) u w (u * w)‖ := norm_nonneg _
  have hwpow : (0 : ℝ) ≤ ‖w‖ ^ (n + 1) := by positivity
  calc ‖tateDefectTrunc (n + 1) u w (u * w)‖ * ‖tateEval u w (d : TateBase)‖
      ≤ (Cn * (50 * ((64 : ℝ) * ‖u * w‖) ^ (n + 1))) * B :=
        mul_le_mul hbound hd (norm_nonneg _) (by positivity)
    _ ≤ (Cn * (50 * (64 ^ (n + 1) * ‖w‖ ^ (n + 1)))) * B := by
        have hstep : Cn * (50 * ((64 : ℝ) * ‖u * w‖) ^ (n + 1))
            ≤ Cn * (50 * (64 ^ (n + 1) * ‖w‖ ^ (n + 1))) := by
          nlinarith [hCn0, hmono]
        nlinarith [hB0, hstep]
    _ = Cn * 50 * 64 ^ (n + 1) * B * ‖w‖ ^ (n + 1) := by ring

/-! ## ★出典の紐付け(`.src`) -/

def exists_bound_tateEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分母の有界性)",
    sectionId := "genell-def-3-3" }

def coeff_evalW_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分子の低次係数は消える)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
