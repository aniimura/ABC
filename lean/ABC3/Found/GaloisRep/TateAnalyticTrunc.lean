import ABC3.Found.GaloisRep.TateSigmaTail

/-!
# Galois (G6) 第 231 ブロック —— **★★★★★★解析値と切り詰めの差をまとめる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★部品を合わせる

第 227(`X`・`Y` の尾)と第 230(`s₁` の尾)を合わせて、
**解析側の `X`・`Y` と形式側の切り詰めの差**を一つの評価にまとめる。

    ‖X_an(u,w) − tateXtrunc (n+1) u w (uw)‖ ≤ 20·(4‖uw‖)^{n+1}
    ‖Y_an(u,w) − tateYtrunc (n+1) u w (uw)‖ ≤ 50·(4‖uw‖)^{n+1}

★解析側の `X_an`・`Y_an` は、形式側の切り詰めの**部分和を `tsum` に、
`partialEval (sigmaSeries 1)` を `∑_{m≥1}f(qᵐ)` に置き換えた**ものである
(第 226 でこの形が `∑_{n∈ℤ}` の和と一致することは確かめた)。

## ★★切り詰めの添字は `n+1` で取る

`s₁` の側の評価(第 230)が `partialEval … (n+1)` の形で出るので、
切り詰めも `n+1` で揃える。
★これは損にならない——`tateDefect − tateDefectTrunc (n+1) ∈ I^{n+1} ⊆ I^n` なので、
`n+1` 次の切り詰めで `I^{n+1}` に入ることを示せば `∀ n, tateDefect ∈ I^n` が出る。

## ★★★指数の帳尻

`‖q‖^{n+2} ≤ ‖q‖^{n+1} ≤ (4‖q‖)^{n+1}`(`‖q‖ ≤ 1/8` より)。
これで尾の評価(`8‖q‖^{n+2}` 型)も `s₁` の評価(`2(4‖q‖)^{n+1}` 型)も
**同じ `(4‖q‖)^{n+1}` の定数倍**に揃う(`pow_le_four_pow_mul`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXanalytic` / `tateYanalytic` | ★★★★★★解析側の `X`・`Y`(`(u,w)` の形) |
| `pow_le_four_pow_mul` | ★指数の帳尻 |
| `norm_three_comb` / `norm_four_comb` | ★三角不等式の組み立て |
| `norm_tateXanalytic_sub_trunc` | ★★★★★★**`X` の差の評価** |
| `norm_tateYanalytic_sub_trunc` | ★★★★★★**`Y` の差の評価** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★★★★★解析側の `X`・`Y` -/

/-- ★★★★★★**解析側の `X(u,w)`**(`q = uw`)——形式側 `tateXtrunc` と同じ形。 -/
noncomputable def tateXanalytic (u w : ℂ) : ℂ :=
  (tateXterm u + ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * u))
    + (tateXterm w + ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w))
    - 2 * ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1))

/-- ★★★★★★**解析側の `Y(u,w)`**。 -/
noncomputable def tateYanalytic (u w : ℂ) : ℂ :=
  (tateYterm u + ∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * u))
    - (tateXterm w + ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w))
    - (tateYterm w + ∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * w))
    + ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1))

/-! ## ★指数と三角不等式の下ごしらえ -/

theorem pow_le_four_pow_mul (r : ℝ) (hr0 : 0 ≤ r) (hr : r ≤ 1 / 8) (n : ℕ) :
    r ^ (n + 2) ≤ (4 * r) ^ (n + 1) := by
  have h1 : r ^ (n + 2) ≤ r ^ (n + 1) :=
    pow_le_pow_of_le_one hr0 (by linarith) (by omega)
  have h2 : r ^ (n + 1) ≤ (4 * r) ^ (n + 1) :=
    pow_le_pow_left₀ hr0 (by linarith) (n + 1)
  linarith

theorem norm_three_comb (Au Aw As : ℂ) :
    ‖Au + Aw + (-2 : ℂ) * As‖ ≤ ‖Au‖ + ‖Aw‖ + 2 * ‖As‖ := by
  have h1 := norm_add_le (Au + Aw) ((-2 : ℂ) * As)
  have h2 := norm_add_le Au Aw
  have h3 : ‖(-2 : ℂ) * As‖ = 2 * ‖As‖ := by
    rw [norm_mul]
    simp
  linarith

theorem norm_four_comb (Au Aw Bw As : ℂ) :
    ‖Au + (-1 : ℂ) * Aw + (-1 : ℂ) * Bw + As‖ ≤ ‖Au‖ + ‖Aw‖ + ‖Bw‖ + ‖As‖ := by
  have h1 := norm_add_le (Au + (-1 : ℂ) * Aw + (-1 : ℂ) * Bw) As
  have h2 := norm_add_le (Au + (-1 : ℂ) * Aw) ((-1 : ℂ) * Bw)
  have h3 := norm_add_le Au ((-1 : ℂ) * Aw)
  have h4 : ‖(-1 : ℂ) * Aw‖ = ‖Aw‖ := by rw [norm_mul]; simp
  have h5 : ‖(-1 : ℂ) * Bw‖ = ‖Bw‖ := by rw [norm_mul]; simp
  linarith

/-! ## ★★★★★★差の評価 -/

/-- ★★★★★★**`X` の切り詰めと解析値の差**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_tateXanalytic_sub_trunc (u w : ℂ) (hu : ‖u‖ ≤ 1) (hw : ‖w‖ ≤ 1)
    (hq : ‖u * w‖ ≤ 1 / 8) (n : ℕ) :
    ‖tateXanalytic u w - tateXtrunc (n + 1) u w (u * w)‖
      ≤ 20 * (4 * ‖u * w‖) ^ (n + 1) := by
  have hq2 : ‖u * w‖ ≤ 1 / 2 := by linarith
  have hq0 : (0 : ℝ) ≤ ‖u * w‖ := norm_nonneg _
  have hA := norm_tateXtail_sub_partialSum_le (u * w) u hq2 hu (n + 1)
  have hB := norm_tateXtail_sub_partialSum_le (u * w) w hq2 hw (n + 1)
  have hC := norm_tateXtail_one_sub_partialEval_le (u * w) hq n
  have hkey : tateXanalytic u w - tateXtrunc (n + 1) u w (u * w)
      = ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * u))
          - partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * u)) (n + 1))
        + ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w))
          - partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
        + (-2) * ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1)))
          - partialEval (sigmaSeries 1) (u * w) (n + 1)) := by
    rw [tateXanalytic, tateXtrunc]
    ring
  rw [hkey]
  have hcomb := norm_three_comb
    ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * u))
      - partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * u)) (n + 1))
    ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w))
      - partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
    ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1)))
      - partialEval (sigmaSeries 1) (u * w) (n + 1))
  have hpow : ‖u * w‖ ^ (n + 1 + 1) ≤ (4 * ‖u * w‖) ^ (n + 1) :=
    pow_le_four_pow_mul _ hq0 hq n
  have hpos : (0 : ℝ) ≤ (4 * ‖u * w‖) ^ (n + 1) := by positivity
  linarith

/-- ★★★★★★**`Y` の切り詰めと解析値の差**。 -/
theorem norm_tateYanalytic_sub_trunc (u w : ℂ) (hu : ‖u‖ ≤ 1) (hw : ‖w‖ ≤ 1)
    (hq : ‖u * w‖ ≤ 1 / 8) (n : ℕ) :
    ‖tateYanalytic u w - tateYtrunc (n + 1) u w (u * w)‖
      ≤ 50 * (4 * ‖u * w‖) ^ (n + 1) := by
  have hq2 : ‖u * w‖ ≤ 1 / 2 := by linarith
  have hq0 : (0 : ℝ) ≤ ‖u * w‖ := norm_nonneg _
  have hA := norm_tateYtail_sub_partialSum_le (u * w) u hq2 hu (n + 1)
  have hB := norm_tateXtail_sub_partialSum_le (u * w) w hq2 hw (n + 1)
  have hB' := norm_tateYtail_sub_partialSum_le (u * w) w hq2 hw (n + 1)
  have hC := norm_tateXtail_one_sub_partialEval_le (u * w) hq n
  have hkey : tateYanalytic u w - tateYtrunc (n + 1) u w (u * w)
      = ((∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * u))
          - partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * u)) (n + 1))
        + (-1) * ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w))
          - partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
        + (-1) * ((∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * w))
          - partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * w)) (n + 1))
        + ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1)))
          - partialEval (sigmaSeries 1) (u * w) (n + 1)) := by
    rw [tateYanalytic, tateYtrunc]
    ring
  rw [hkey]
  have hcomb := norm_four_comb
    ((∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * u))
      - partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * u)) (n + 1))
    ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w))
      - partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
    ((∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * w))
      - partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * w)) (n + 1))
    ((∑' m : ℕ, tateXterm ((u * w) ^ (m + 1)))
      - partialEval (sigmaSeries 1) (u * w) (n + 1))
  have hpow : ‖u * w‖ ^ (n + 1 + 1) ≤ (4 * ‖u * w‖) ^ (n + 1) :=
    pow_le_four_pow_mul _ hq0 hq n
  have hpos : (0 : ℝ) ≤ (4 * ‖u * w‖) ^ (n + 1) := by positivity
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def tateXanalytic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——解析側の X(u,w))",
    sectionId := "genell-def-3-3" }

def norm_tateXanalytic_sub_trunc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X の切り詰めと解析値の差)",
    sectionId := "genell-def-3-3" }

def norm_tateYanalytic_sub_trunc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Y の切り詰めと解析値の差)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
