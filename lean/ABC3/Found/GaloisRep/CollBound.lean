import ABC3.Found.GaloisRep.CollDescent

/-!
# Galois (G6) 第 254 ブロック —— **★★★★★★★切り詰めた共線性の差の評価**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★在庫がそのまま 3 点に効いた

第 231 の `norm_tateXanalytic_sub_trunc`・`norm_tateYanalytic_sub_trunc` と
第 234 の `norm_tateXtrunc_le`・`norm_tateYtrunc_le` は `(u, w)` の対に対する評価である。
★★3 点はどれも `(u, vw)`・`(v, uw)`・`(w, uv)` という **`(点, 相方)` の対**なので、
**そのまま 3 回当てられる**。しかも相方との積はどれも `q = uvw` なので、
**評価はどの点でも同じ `ε = 50(4‖q‖)^{n+1}`** になる。

★`‖u‖, ‖v‖, ‖w‖ ≤ 1/8` と取れば `‖tateXterm t‖ ≤ 4‖t‖ ≤ 1/2`、
`‖tateYterm t‖ ≤ 8‖t‖ ≤ 1` なので、**切り詰めの大きさ `M` が `n` だけで決まる**。
第 237 では `u` が固定されていたので `‖tateXterm u‖` を残したが、
ここは 3 変数とも動くので数値で潰す。

## ★★★★★★行列式の差

    D − D' = Σᵢ [ (Xᵢ − Xᵢ')(Yⱼ − Yₖ) + Xᵢ'((Yⱼ − Yⱼ') − (Yₖ − Yₖ')) ]

★各項は `ε·2M` と `M·2ε` で押さえられ、合計 `12Mε`(`norm_coll_diff_le`)。
第 234 の `norm_defect_diff_le` が `(3M²+6M+1)ε` だったのに比べ、
**行列式は双線形なので `M` の 1 次で済む**。

## ★★★★★★★到達点

解析側は 0(第 249)なので、

    ‖collDefectTrunc (n+1) u v w (uvw)‖ ≤ 12(25(n+1)+8) · 50(4‖uvw‖)^{n+1}

★右辺は `‖uvw‖^{n+1} = ‖uv‖^{n+1}‖w‖^{n+1}` を含むので、`u, v` を固定して `w` を
動かせば `C‖w‖^{n+1}` の形になり、第 236 の `X_pow_dvd_of_norm_le` が効く。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `norm_sum6` | ★6 項の三角不等式 |
| `norm_coll_diff_le` | ★★★★★★**行列式の差は `12Mε`** |
| `norm_collDefectTrunc_le` | ★★★★★★★**切り詰めた共線性の差の評価** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★6 項の三角不等式 -/

theorem norm_sum6 (a b c d e f : ℂ) :
    ‖a + b + c + d + e + f‖ ≤ ‖a‖ + ‖b‖ + ‖c‖ + ‖d‖ + ‖e‖ + ‖f‖ := by
  have h2 := norm_add_le (a + b + c + d + e) f
  have h3 := norm_add_le (a + b + c + d) e
  have h4 := norm_add_le (a + b + c) d
  have h5 := norm_add_le (a + b) c
  have h6 := norm_add_le a b
  linarith

/-! ## ★★★★★★行列式の差 -/

/-- ★★★★★★**共線性の行列式の差**——値の大きさ `M` と差 `ε` で押さえる。

★行列式は双線形なので `M` の 1 次で済む(第 234 の `norm_defect_diff_le` は
3 次式が要った)。 -/
theorem norm_coll_diff_le (X1 Y1 X2 Y2 X3 Y3 X1' Y1' X2' Y2' X3' Y3' : ℂ) (M ε : ℝ)
    (hM : 0 ≤ M) (hε : 0 ≤ ε)
    (hY1 : ‖Y1‖ ≤ M) (hY2 : ‖Y2‖ ≤ M) (hY3 : ‖Y3‖ ≤ M)
    (hX1' : ‖X1'‖ ≤ M) (hX2' : ‖X2'‖ ≤ M) (hX3' : ‖X3'‖ ≤ M)
    (dX1 : ‖X1 - X1'‖ ≤ ε) (dX2 : ‖X2 - X2'‖ ≤ ε) (dX3 : ‖X3 - X3'‖ ≤ ε)
    (dY1 : ‖Y1 - Y1'‖ ≤ ε) (dY2 : ‖Y2 - Y2'‖ ≤ ε) (dY3 : ‖Y3 - Y3'‖ ≤ ε) :
    ‖(X1 * (Y2 - Y3) + X2 * (Y3 - Y1) + X3 * (Y1 - Y2))
        - (X1' * (Y2' - Y3') + X2' * (Y3' - Y1') + X3' * (Y1' - Y2'))‖ ≤ 12 * M * ε := by
  have hkey : (X1 * (Y2 - Y3) + X2 * (Y3 - Y1) + X3 * (Y1 - Y2))
      - (X1' * (Y2' - Y3') + X2' * (Y3' - Y1') + X3' * (Y1' - Y2'))
      = (X1 - X1') * (Y2 - Y3) + X1' * ((Y2 - Y2') - (Y3 - Y3'))
        + (X2 - X2') * (Y3 - Y1) + X2' * ((Y3 - Y3') - (Y1 - Y1'))
        + (X3 - X3') * (Y1 - Y2) + X3' * ((Y1 - Y1') - (Y2 - Y2')) := by ring
  rw [hkey]
  have h := norm_sum6 ((X1 - X1') * (Y2 - Y3)) (X1' * ((Y2 - Y2') - (Y3 - Y3')))
    ((X2 - X2') * (Y3 - Y1)) (X2' * ((Y3 - Y3') - (Y1 - Y1')))
    ((X3 - X3') * (Y1 - Y2)) (X3' * ((Y1 - Y1') - (Y2 - Y2')))
  have b1 : ‖(X1 - X1') * (Y2 - Y3)‖ ≤ ε * (2 * M) := by
    rw [norm_mul]
    have hs : ‖Y2 - Y3‖ ≤ 2 * M := by have := norm_sub_le Y2 Y3; linarith
    exact mul_le_mul dX1 hs (norm_nonneg _) hε
  have b2 : ‖X1' * ((Y2 - Y2') - (Y3 - Y3'))‖ ≤ M * (2 * ε) := by
    rw [norm_mul]
    have hs : ‖(Y2 - Y2') - (Y3 - Y3')‖ ≤ 2 * ε := by
      have := norm_sub_le (Y2 - Y2') (Y3 - Y3'); linarith
    exact mul_le_mul hX1' hs (norm_nonneg _) hM
  have b3 : ‖(X2 - X2') * (Y3 - Y1)‖ ≤ ε * (2 * M) := by
    rw [norm_mul]
    have hs : ‖Y3 - Y1‖ ≤ 2 * M := by have := norm_sub_le Y3 Y1; linarith
    exact mul_le_mul dX2 hs (norm_nonneg _) hε
  have b4 : ‖X2' * ((Y3 - Y3') - (Y1 - Y1'))‖ ≤ M * (2 * ε) := by
    rw [norm_mul]
    have hs : ‖(Y3 - Y3') - (Y1 - Y1')‖ ≤ 2 * ε := by
      have := norm_sub_le (Y3 - Y3') (Y1 - Y1'); linarith
    exact mul_le_mul hX2' hs (norm_nonneg _) hM
  have b5 : ‖(X3 - X3') * (Y1 - Y2)‖ ≤ ε * (2 * M) := by
    rw [norm_mul]
    have hs : ‖Y1 - Y2‖ ≤ 2 * M := by have := norm_sub_le Y1 Y2; linarith
    exact mul_le_mul dX3 hs (norm_nonneg _) hε
  have b6 : ‖X3' * ((Y1 - Y1') - (Y2 - Y2'))‖ ≤ M * (2 * ε) := by
    rw [norm_mul]
    have hs : ‖(Y1 - Y1') - (Y2 - Y2')‖ ≤ 2 * ε := by
      have := norm_sub_le (Y1 - Y1') (Y2 - Y2'); linarith
    exact mul_le_mul hX3' hs (norm_nonneg _) hM
  nlinarith [h, b1, b2, b3, b4, b5, b6]

/-! ## ★★★★★★★切り詰めた共線性の差の評価 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★**切り詰めた共線性の差の評価**——解析側が 0(第 249)であることから。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_collDefectTrunc_le (u v w : ℂ) (hu0 : u ≠ 0) (hv0 : v ≠ 0) (hw0 : w ≠ 0)
    (hu : ‖u‖ ≤ 1 / 8) (hv : ‖v‖ ≤ 1 / 8) (hw : ‖w‖ ≤ 1 / 8) (n : ℕ) :
    ‖collDefectTrunc (n + 1) u v w (u * v * w)‖
      ≤ 12 * (25 * ((n : ℝ) + 1) + 8) * (50 * (4 * ‖u * v * w‖) ^ (n + 1)) := by
  have nu := norm_nonneg u
  have nv := norm_nonneg v
  have nw := norm_nonneg w
  have hvw : ‖v * w‖ ≤ 1 / 64 := by rw [norm_mul]; nlinarith
  have huw : ‖u * w‖ ≤ 1 / 64 := by rw [norm_mul]; nlinarith
  have huv : ‖u * v‖ ≤ 1 / 64 := by rw [norm_mul]; nlinarith
  have e1 : u * (v * w) = u * v * w := by ring
  have e2 : v * (u * w) = u * v * w := by ring
  have e3 : w * (u * v) = u * v * w := by ring
  have hq : ‖u * v * w‖ ≤ 1 / 512 := by
    rw [← e1, norm_mul]; nlinarith [norm_nonneg (v * w)]
  have hq8 : ‖u * v * w‖ ≤ 1 / 8 := by linarith
  set M : ℝ := 25 * ((n : ℝ) + 1) + 7 with hMdef
  have hn0 : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
  have hM0 : (0 : ℝ) ≤ M := by simp only [hMdef]; linarith
  have hxt : ∀ t : ℂ, ‖t‖ ≤ 1 / 8 → ‖tateXterm t‖ ≤ 1 / 2 := by
    intro t ht
    have := norm_tateXterm_le_of_small (t := t) (by linarith)
    have h2 := norm_nonneg t
    linarith
  have hyt : ∀ t : ℂ, ‖t‖ ≤ 1 / 8 → ‖tateYterm t‖ ≤ 1 := by
    intro t ht
    have := norm_tateYterm_le_of_small (t := t) (by linarith)
    have h2 := norm_nonneg t
    linarith
  have bX1 : ‖tateXtrunc (n + 1) u (v * w) (u * v * w)‖ ≤ M := by
    have h := norm_tateXtrunc_le u (v * w) (by linarith) (by linarith) (by rw [e1]; linarith) n
    rw [e1] at h; have := hxt u hu; simp only [hMdef]; linarith
  have bX2 : ‖tateXtrunc (n + 1) v (u * w) (u * v * w)‖ ≤ M := by
    have h := norm_tateXtrunc_le v (u * w) (by linarith) (by linarith) (by rw [e2]; linarith) n
    rw [e2] at h; have := hxt v hv; simp only [hMdef]; linarith
  have bX3 : ‖tateXtrunc (n + 1) w (u * v) (u * v * w)‖ ≤ M := by
    have h := norm_tateXtrunc_le w (u * v) (by linarith) (by linarith) (by rw [e3]; linarith) n
    rw [e3] at h; have := hxt w hw; simp only [hMdef]; linarith
  have bY1 : ‖tateYtrunc (n + 1) u (v * w) (u * v * w)‖ ≤ M := by
    have h := norm_tateYtrunc_le u (v * w) (by linarith) (by linarith) (by rw [e1]; linarith) n
    rw [e1] at h; have := hyt u hu; simp only [hMdef]; linarith
  have bY2 : ‖tateYtrunc (n + 1) v (u * w) (u * v * w)‖ ≤ M := by
    have h := norm_tateYtrunc_le v (u * w) (by linarith) (by linarith) (by rw [e2]; linarith) n
    rw [e2] at h; have := hyt v hv; simp only [hMdef]; linarith
  have bY3 : ‖tateYtrunc (n + 1) w (u * v) (u * v * w)‖ ≤ M := by
    have h := norm_tateYtrunc_le w (u * v) (by linarith) (by linarith) (by rw [e3]; linarith) n
    rw [e3] at h; have := hyt w hw; simp only [hMdef]; linarith
  set ε : ℝ := 50 * (4 * ‖u * v * w‖) ^ (n + 1) with hεdef
  have hpowpos : (0 : ℝ) ≤ (4 * ‖u * v * w‖) ^ (n + 1) := by positivity
  have hε0 : (0 : ℝ) ≤ ε := by simp only [hεdef]; linarith
  have hε1 : ε ≤ 1 := by
    have h4 : 4 * ‖u * v * w‖ ≤ 1 / 128 := by linarith
    have hpow : (4 * ‖u * v * w‖) ^ (n + 1) ≤ ((1 : ℝ) / 128) ^ (n + 1) :=
      pow_le_pow_left₀ (by positivity) h4 _
    have hpow2 : ((1 : ℝ) / 128) ^ (n + 1) ≤ ((1 : ℝ) / 128) ^ 1 :=
      pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega)
    simp only [hεdef]
    norm_num at hpow2
    linarith
  have dX1 : ‖tateXanalytic u (v * w) - tateXtrunc (n + 1) u (v * w) (u * v * w)‖ ≤ ε := by
    have h := norm_tateXanalytic_sub_trunc u (v * w) (by linarith) (by linarith)
      (by rw [e1]; linarith) n
    rw [e1] at h; simp only [hεdef]; linarith
  have dX2 : ‖tateXanalytic v (u * w) - tateXtrunc (n + 1) v (u * w) (u * v * w)‖ ≤ ε := by
    have h := norm_tateXanalytic_sub_trunc v (u * w) (by linarith) (by linarith)
      (by rw [e2]; linarith) n
    rw [e2] at h; simp only [hεdef]; linarith
  have dX3 : ‖tateXanalytic w (u * v) - tateXtrunc (n + 1) w (u * v) (u * v * w)‖ ≤ ε := by
    have h := norm_tateXanalytic_sub_trunc w (u * v) (by linarith) (by linarith)
      (by rw [e3]; linarith) n
    rw [e3] at h; simp only [hεdef]; linarith
  have dY1 : ‖tateYanalytic u (v * w) - tateYtrunc (n + 1) u (v * w) (u * v * w)‖ ≤ ε := by
    have h := norm_tateYanalytic_sub_trunc u (v * w) (by linarith) (by linarith)
      (by rw [e1]; linarith) n
    rw [e1] at h; simp only [hεdef]; linarith
  have dY2 : ‖tateYanalytic v (u * w) - tateYtrunc (n + 1) v (u * w) (u * v * w)‖ ≤ ε := by
    have h := norm_tateYanalytic_sub_trunc v (u * w) (by linarith) (by linarith)
      (by rw [e2]; linarith) n
    rw [e2] at h; simp only [hεdef]; linarith
  have dY3 : ‖tateYanalytic w (u * v) - tateYtrunc (n + 1) w (u * v) (u * v * w)‖ ≤ ε := by
    have h := norm_tateYanalytic_sub_trunc w (u * v) (by linarith) (by linarith)
      (by rw [e3]; linarith) n
    rw [e3] at h; simp only [hεdef]; linarith
  have habs : ∀ (A T : ℂ), ‖A - T‖ ≤ ε → ‖T‖ ≤ M → ‖A‖ ≤ M + ε := by
    intro A T h1 h2
    have h := norm_add_le (A - T) T
    simp only [sub_add_cancel] at h
    linarith
  have hzero := collinear_analytic_uvw u v w hu0 hv0 hw0 (by linarith) (by linarith) (by linarith)
  have hdiff := norm_coll_diff_le
    (tateXtrunc (n + 1) u (v * w) (u * v * w)) (tateYtrunc (n + 1) u (v * w) (u * v * w))
    (tateXtrunc (n + 1) v (u * w) (u * v * w)) (tateYtrunc (n + 1) v (u * w) (u * v * w))
    (tateXtrunc (n + 1) w (u * v) (u * v * w)) (tateYtrunc (n + 1) w (u * v) (u * v * w))
    (tateXanalytic u (v * w)) (tateYanalytic u (v * w))
    (tateXanalytic v (u * w)) (tateYanalytic v (u * w))
    (tateXanalytic w (u * v)) (tateYanalytic w (u * v))
    (M + ε) ε (by linarith) hε0
    (by linarith) (by linarith) (by linarith)
    (habs _ _ dX1 bX1) (habs _ _ dX2 bX2) (habs _ _ dX3 bX3)
    (by rw [norm_sub_rev]; exact dX1) (by rw [norm_sub_rev]; exact dX2)
    (by rw [norm_sub_rev]; exact dX3)
    (by rw [norm_sub_rev]; exact dY1) (by rw [norm_sub_rev]; exact dY2)
    (by rw [norm_sub_rev]; exact dY3)
  rw [hzero, sub_zero] at hdiff
  have hform : collDefectTrunc (n + 1) u v w (u * v * w)
      = tateXtrunc (n + 1) u (v * w) (u * v * w)
          * (tateYtrunc (n + 1) v (u * w) (u * v * w) - tateYtrunc (n + 1) w (u * v) (u * v * w))
        + tateXtrunc (n + 1) v (u * w) (u * v * w)
          * (tateYtrunc (n + 1) w (u * v) (u * v * w) - tateYtrunc (n + 1) u (v * w) (u * v * w))
        + tateXtrunc (n + 1) w (u * v) (u * v * w)
          * (tateYtrunc (n + 1) u (v * w) (u * v * w)
            - tateYtrunc (n + 1) v (u * w) (u * v * w)) := rfl
  rw [hform]
  rw [hMdef] at hdiff
  nlinarith [hdiff, hε0, hε1]

/-! ## ★出典の紐付け(`.src`) -/

def norm_coll_diff_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性の行列式の差)",
    sectionId := "genell-def-3-3" }

def norm_collDefectTrunc_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——切り詰めた共線性の差の評価)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
