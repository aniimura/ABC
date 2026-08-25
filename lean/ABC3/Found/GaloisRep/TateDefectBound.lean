import ABC3.Found.GaloisRep.TateBounds

/-!
# Galois (G6) 第 235 ブロック —— **★★★★★★★★切り詰めた差は `O(q^{n+1})`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★解析の段の本体

第 220–234 で積み上げた部品を一つにまとめる。`u = exp(2πiz)`、`w = exp(2πi(τ−z))`、
`q = uw = exp(2πiτ)` として、`‖q‖ ≤ 1/128` かつ `‖w‖ ≤ 1/2` なら:

    ‖tateDefectTrunc (n+1) u w (uw)‖ ≤ (3M² + 6M + 1) · 50 · (64‖q‖)^{n+1}

ここで `M = ‖f(u)‖ + ‖g(u)‖ + 30(n+1) + 80`。★`(64‖q‖)^{n+1} = 64^{n+1}‖q‖^{n+1}` なので
これは **`C_n · ‖q‖^{n+1}`** の形である(定数は `n` と `u` に依る)。

## ★★★組み立て

| 段 | 使うもの |
|---|---|
| 解析側の差は厳密に 0 | 第 232 `tate_equation_uw` |
| `X`・`Y` の差 | 第 231 |
| `a₄`・`a₆` の差 | 第 233 |
| 値の大きさ | 第 234 |
| 差の差を足す | 第 234 `norm_defect_diff_le` |

★★**四つの差をすべて `ε := 50(64‖q‖)^{n+1}` に揃える**のが要点である。
`(4‖q‖)^{n+1} ≤ (16‖q‖)^{n+1} ≤ (64‖q‖)^{n+1}` なので、底を一番大きいものに合わせれば
`norm_defect_diff_le` が一度で当たる。

★★★`‖w‖ ≤ 1/2` は仮定として持つ。`w = q/u` は `u` を固定して `q → 0` とすれば
いずれ満たされる——最終段で `u` を固定して `q → 0` を見るときの前提である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_defect_analytic_zero` | ★★★解析側の差は厳密に 0 |
| `norm_tateDefectTrunc_le` | ★★★★★★★★**切り詰めた差は `O(q^{n+1})`** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-- ★★★**解析側の差は厳密に 0**——第 232 の方程式の言い換え。 -/
theorem tate_defect_analytic_zero (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    tateYanalytic (Complex.exp (2 * ↑π * I * z))
        (Complex.exp (2 * ↑π * I * (subUH z τ him))) ^ 2
      + tateXanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him)))
        * tateYanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him)))
      - (tateXanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him))) ^ 3
        + (-5 * sigmaSum 3 τ) * tateXanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him)))
        + (-(5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12)) = 0 := by
  have h := tate_equation_uw z τ him
  linear_combination h

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★**切り詰めた差は `O(q^{n+1})` である**。

★これが段 6 の解析の本体である。`u` を固定して `q → 0` を見るとき、各 `n` について
分子が `q^{n+1}` の位数で消えることがこれで出る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_tateDefectTrunc_le (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im)
    (hq : ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ ≤ 1 / 128)
    (hw : ‖Complex.exp (2 * ↑π * I * ((subUH z τ him : UpperHalfPlane) : ℂ))‖ ≤ 1 / 2)
    (n : ℕ) :
    ‖tateDefectTrunc (n + 1) (Complex.exp (2 * ↑π * I * z))
        (Complex.exp (2 * ↑π * I * (subUH z τ him)))
        (Complex.exp (2 * ↑π * I * z) * Complex.exp (2 * ↑π * I * (subUH z τ him)))‖
      ≤ (3 * (‖tateXterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖
          + ‖tateYterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖ + 30 * ((n : ℝ) + 1) + 80) ^ 2
        + 6 * (‖tateXterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖
          + ‖tateYterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖ + 30 * ((n : ℝ) + 1) + 80) + 1)
        * (50 * ((64 : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖) ^ (n + 1)) := by
  set u := Complex.exp (2 * ↑π * I * (z : ℂ)) with hudef
  set w := Complex.exp (2 * ↑π * I * ((subUH z τ him : UpperHalfPlane) : ℂ)) with hwdef
  have hmul : u * w = Complex.exp (2 * ↑π * I * (τ : ℂ)) := exp_mul_exp_sub z τ him
  have hu : ‖u‖ ≤ 1 := (UpperHalfPlane.norm_exp_two_pi_I_lt_one z).le
  have hw1 : ‖w‖ ≤ 1 := by linarith
  have hquw : ‖u * w‖ ≤ 1 / 128 := by rw [hmul]; exact hq
  have hq8 : ‖u * w‖ ≤ 1 / 8 := by linarith
  have hq0 : (0 : ℝ) ≤ ‖u * w‖ := norm_nonneg _
  have hnn : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
  set M : ℝ := ‖tateXterm u‖ + ‖tateYterm u‖ + 30 * ((n : ℝ) + 1) + 80 with hMdef
  have hfu : (0 : ℝ) ≤ ‖tateXterm u‖ := norm_nonneg _
  have hgu : (0 : ℝ) ≤ ‖tateYterm u‖ := norm_nonneg _
  have hMnn : (0 : ℝ) ≤ M := by simp only [hMdef]; linarith
  have h4le64 : ((4 : ℝ) * ‖u * w‖) ^ (n + 1) ≤ ((64 : ℝ) * ‖u * w‖) ^ (n + 1) :=
    pow_le_pow_left₀ (by positivity) (by linarith) _
  have h16le64 : ((16 : ℝ) * ‖u * w‖) ^ (n + 1) ≤ ((64 : ℝ) * ‖u * w‖) ^ (n + 1) :=
    pow_le_pow_left₀ (by positivity) (by linarith) _
  have h64le1 : ((64 : ℝ) * ‖u * w‖) ^ (n + 1) ≤ 1 :=
    pow_le_one₀ (by positivity) (by linarith)
  have h64nn : (0 : ℝ) ≤ ((64 : ℝ) * ‖u * w‖) ^ (n + 1) := by positivity
  have dX0 := norm_tateXanalytic_sub_trunc u w hu hw1 hq8 n
  have dY0 := norm_tateYanalytic_sub_trunc u w hu hw1 hq8 n
  have dA40 : ‖(-5 * sigmaSum 3 τ) - partialEval tateA4 (u * w) (n + 1)‖
      ≤ 10 * ((16 : ℝ) * ‖u * w‖) ^ (n + 1) := by
    rw [hmul]
    exact norm_a4_sub_partialEval_le τ (by linarith) n
  have dA60 : ‖(-(5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12)
      - partialEval tateA6 (u * w) (n + 1)‖ ≤ 2 * ((64 : ℝ) * ‖u * w‖) ^ (n + 1) := by
    rw [hmul]
    exact norm_a6_sub_partialEval_le τ (by linarith) n
  have bXtr := norm_tateXtrunc_le u w hu hw hq8 n
  have bYtr := norm_tateYtrunc_le u w hu hw hq8 n
  have bA4tr : ‖partialEval tateA4 (u * w) (n + 1)‖ ≤ 5 * (n : ℝ) :=
    norm_partialEval_tateA4_le (u * w) hquw n
  have bXan : ‖tateXanalytic u w‖ ≤ M := by
    have hb := norm_add_le (tateXanalytic u w - tateXtrunc (n + 1) u w (u * w))
      (tateXtrunc (n + 1) u w (u * w))
    simp only [sub_add_cancel] at hb
    simp only [hMdef]
    nlinarith [dX0, bXtr, h4le64, h64le1, hb]
  have bYan : ‖tateYanalytic u w‖ ≤ M := by
    have hb := norm_add_le (tateYanalytic u w - tateYtrunc (n + 1) u w (u * w))
      (tateYtrunc (n + 1) u w (u * w))
    simp only [sub_add_cancel] at hb
    simp only [hMdef]
    nlinarith [dY0, bYtr, h4le64, h64le1, hb]
  have bA4an : ‖(-5 * sigmaSum 3 τ : ℂ)‖ ≤ M := by
    have hb := norm_add_le ((-5 * sigmaSum 3 τ) - partialEval tateA4 (u * w) (n + 1))
      (partialEval tateA4 (u * w) (n + 1))
    simp only [sub_add_cancel] at hb
    simp only [hMdef]
    nlinarith [dA40, bA4tr, h16le64, h64le1, hb]
  have bXtr' : ‖tateXtrunc (n + 1) u w (u * w)‖ ≤ M := by simp only [hMdef]; linarith
  have bYtr' : ‖tateYtrunc (n + 1) u w (u * w)‖ ≤ M := by simp only [hMdef]; linarith
  have bA4tr' : ‖partialEval tateA4 (u * w) (n + 1)‖ ≤ M := by simp only [hMdef]; linarith
  set ε : ℝ := 50 * ((64 : ℝ) * ‖u * w‖) ^ (n + 1) with hεdef
  have hεnn : (0 : ℝ) ≤ ε := by simp only [hεdef]; positivity
  have dX : ‖tateXtrunc (n + 1) u w (u * w) - tateXanalytic u w‖ ≤ ε := by
    rw [norm_sub_rev]
    simp only [hεdef]
    linarith
  have dY : ‖tateYtrunc (n + 1) u w (u * w) - tateYanalytic u w‖ ≤ ε := by
    rw [norm_sub_rev]
    simp only [hεdef]
    linarith
  have dA4 : ‖partialEval tateA4 (u * w) (n + 1) - (-5 * sigmaSum 3 τ)‖ ≤ ε := by
    rw [norm_sub_rev]
    simp only [hεdef]
    linarith
  have dA6 : ‖partialEval tateA6 (u * w) (n + 1)
      - (-(5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12)‖ ≤ ε := by
    rw [norm_sub_rev]
    simp only [hεdef]
    linarith
  have hzero := tate_defect_analytic_zero z τ him
  have hkey : tateDefectTrunc (n + 1) u w (u * w)
      = (tateYtrunc (n + 1) u w (u * w) ^ 2
          + tateXtrunc (n + 1) u w (u * w) * tateYtrunc (n + 1) u w (u * w)
          - (tateXtrunc (n + 1) u w (u * w) ^ 3
            + partialEval tateA4 (u * w) (n + 1) * tateXtrunc (n + 1) u w (u * w)
            + partialEval tateA6 (u * w) (n + 1)))
        - (tateYanalytic u w ^ 2 + tateXanalytic u w * tateYanalytic u w
          - (tateXanalytic u w ^ 3 + (-5 * sigmaSum 3 τ) * tateXanalytic u w
            + (-(5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12))) := by
    rw [tateDefectTrunc, hzero]
    ring
  rw [hkey, ← hmul]
  exact norm_defect_diff_le _ _ _ _ _ _ _ _ M ε hMnn hεnn bXtr' bXan bYtr' bYan
    bA4tr' bA4an dX dY dA4 dA6

/-! ## ★出典の紐付け(`.src`) -/

def norm_tateDefectTrunc_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——切り詰めた差は O(q^{n+1}))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
