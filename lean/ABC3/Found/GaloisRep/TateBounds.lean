import ABC3.Found.GaloisRep.TateCoeffTail

/-!
# Galois (G6) 第 234 ブロック —— **★★★★★★方程式の差の差と大きさの評価**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★差を足すには「大きさ」も要る

第 231・233 で四つの部品(`X`・`Y`・`a₄`・`a₆`)の**差**が揃った。
方程式の差は二次・三次の項を含むので、差を足すには**値の大きさ**も要る:

    Y² − Y'² = (Y−Y')(Y+Y'),  X³ − X'³ = (X−X')(X² + XX' + X'²),  …

★これをまとめたのが `norm_defect_diff_le` である:

    ‖D(X,Y,A₄,A₆) − D(X',Y',A₄',A₆')‖ ≤ (3M² + 6M + 1)·ε

(`M` は値の大きさの上界、`ε` は差の上界。`D` は Weierstrass の左辺 − 右辺。)

## ★★★★切り詰めの大きさ

`‖q‖ ≤ 1/8` かつ `‖u‖ ≤ 1`、`‖w‖ ≤ 1/2` なら、切り詰めの各項は絶対定数で押さえられる:

| 項 | 上界 |
|---|---|
| `f(q^{m+1}c)`(`‖c‖ ≤ 1`) | `4` |
| `g(q^{m+1}c)` | `8` |
| `partialEval (sigmaSeries k) q (n+1)` | `n`(各項が `1` 以下) |
| `f(w)`(`‖w‖ ≤ 1/2`) | `2` |
| `g(w)` | `4` |

★★**`u` に依存するのは `f(u)`・`g(u)` だけ**である。`w = q/u` は `q → 0` で 0 に行くので
`f(w)`・`g(w)` は絶対定数で押さえられる。`u` を固定して `q → 0` を見る、という
最終段の構図がここに現れている。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `norm_sum7` | ★7 項の三角不等式 |
| `norm_defect_diff_le` | ★★★★★★**方程式の差の差** |
| `norm_partialSum_le` | ★部分和の大きさ |
| `norm_tateXtrunc_le` / `norm_tateYtrunc_le` | ★★★★★切り詰めの大きさ |
| `norm_partialEval_tateA4_le` / `_tateA6_le` | ★★★★係数の切り詰めの大きさ |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★三角不等式の下ごしらえ -/

theorem norm_sum7 (a b c d e f g : ℂ) :
    ‖a + b + c + d + e + f + g‖ ≤ ‖a‖ + ‖b‖ + ‖c‖ + ‖d‖ + ‖e‖ + ‖f‖ + ‖g‖ := by
  have h1 := norm_add_le (a + b + c + d + e + f) g
  have h2 := norm_add_le (a + b + c + d + e) f
  have h3 := norm_add_le (a + b + c + d) e
  have h4 := norm_add_le (a + b + c) d
  have h5 := norm_add_le (a + b) c
  have h6 := norm_add_le a b
  linarith

/-! ## ★★★★★★方程式の差の差 -/

/-- ★★★★★★**方程式の差の差**——二組の値の差を、値の大きさ `M` と差 `ε` で押さえる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_defect_diff_le (X Y A4 A6 X' Y' A4' A6' : ℂ) (M ε : ℝ)
    (hM : 0 ≤ M) (hε : 0 ≤ ε)
    (hX : ‖X‖ ≤ M) (hX' : ‖X'‖ ≤ M) (hY : ‖Y‖ ≤ M) (hY' : ‖Y'‖ ≤ M)
    (_hA4 : ‖A4‖ ≤ M) (hA4' : ‖A4'‖ ≤ M)
    (dX : ‖X - X'‖ ≤ ε) (dY : ‖Y - Y'‖ ≤ ε) (dA4 : ‖A4 - A4'‖ ≤ ε)
    (dA6 : ‖A6 - A6'‖ ≤ ε) :
    ‖(Y ^ 2 + X * Y - (X ^ 3 + A4 * X + A6))
        - (Y' ^ 2 + X' * Y' - (X' ^ 3 + A4' * X' + A6'))‖
      ≤ (3 * M ^ 2 + 6 * M + 1) * ε := by
  have hkey : (Y ^ 2 + X * Y - (X ^ 3 + A4 * X + A6))
      - (Y' ^ 2 + X' * Y' - (X' ^ 3 + A4' * X' + A6'))
      = (Y - Y') * (Y + Y') + (X - X') * Y + X' * (Y - Y')
        + (-((X - X') * (X ^ 2 + X * X' + X' ^ 2))) + (-((A4 - A4') * X))
        + (-(A4' * (X - X'))) + (-(A6 - A6')) := by ring
  rw [hkey]
  have h := norm_sum7 ((Y - Y') * (Y + Y')) ((X - X') * Y) (X' * (Y - Y'))
    (-((X - X') * (X ^ 2 + X * X' + X' ^ 2))) (-((A4 - A4') * X))
    (-(A4' * (X - X'))) (-(A6 - A6'))
  have b1 : ‖(Y - Y') * (Y + Y')‖ ≤ ε * (2 * M) := by
    rw [norm_mul]
    have hsum : ‖Y + Y'‖ ≤ 2 * M := by
      have hy := norm_add_le Y Y'
      linarith
    exact mul_le_mul dY hsum (norm_nonneg _) hε
  have b2 : ‖(X - X') * Y‖ ≤ ε * M := by
    rw [norm_mul]
    exact mul_le_mul dX hY (norm_nonneg _) hε
  have b3 : ‖X' * (Y - Y')‖ ≤ M * ε := by
    rw [norm_mul]
    exact mul_le_mul hX' dY (norm_nonneg _) hM
  have b4 : ‖-((X - X') * (X ^ 2 + X * X' + X' ^ 2))‖ ≤ ε * (3 * M ^ 2) := by
    rw [norm_neg, norm_mul]
    have hsum : ‖X ^ 2 + X * X' + X' ^ 2‖ ≤ 3 * M ^ 2 := by
      have h1 := norm_add_le (X ^ 2 + X * X') (X' ^ 2)
      have h2 := norm_add_le (X ^ 2) (X * X')
      have e1 : ‖X ^ 2‖ = ‖X‖ ^ 2 := by rw [norm_pow]
      have e2 : ‖X' ^ 2‖ = ‖X'‖ ^ 2 := by rw [norm_pow]
      have e3 : ‖X * X'‖ = ‖X‖ * ‖X'‖ := by rw [norm_mul]
      nlinarith [norm_nonneg X, norm_nonneg X']
    exact mul_le_mul dX hsum (norm_nonneg _) hε
  have b5 : ‖-((A4 - A4') * X)‖ ≤ ε * M := by
    rw [norm_neg, norm_mul]
    exact mul_le_mul dA4 hX (norm_nonneg _) hε
  have b6 : ‖-(A4' * (X - X'))‖ ≤ M * ε := by
    rw [norm_neg, norm_mul]
    exact mul_le_mul hA4' dX (norm_nonneg _) hM
  have b7 : ‖-(A6 - A6')‖ ≤ ε := by
    rw [norm_neg]
    exact dA6
  nlinarith [h, b1, b2, b3, b4, b5, b6, b7]

/-! ## ★部分和と項の大きさ -/

theorem norm_partialSum_le (f : ℕ → ℂ) (B : ℝ) (hB : ∀ i, ‖f i‖ ≤ B) (n : ℕ) :
    ‖partialSum f n‖ ≤ n * B := by
  rw [partialSum]
  refine (norm_sum_le _ _).trans ?_
  calc ∑ i ∈ Finset.range n, ‖f i‖ ≤ ∑ _i ∈ Finset.range n, B :=
        Finset.sum_le_sum (fun i _ => hB i)
    _ = n * B := by rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]

theorem norm_tateXterm_pow_le (q c : ℂ) (hq : ‖q‖ ≤ 1 / 2) (hc : ‖c‖ ≤ 1) (m : ℕ) :
    ‖tateXterm (q ^ (m + 1) * c)‖ ≤ 4 := by
  have h := norm_tateXterm_le_of_small (norm_pow_succ_mul_le_half hq hc m)
  have h2 := norm_pow_succ_mul_le_half hq hc m
  linarith

theorem norm_tateYterm_pow_le (q c : ℂ) (hq : ‖q‖ ≤ 1 / 2) (hc : ‖c‖ ≤ 1) (m : ℕ) :
    ‖tateYterm (q ^ (m + 1) * c)‖ ≤ 8 := by
  have h := norm_tateYterm_le_of_small (norm_pow_succ_mul_le_half hq hc m)
  have h2 := norm_pow_succ_mul_le_half hq hc m
  linarith

theorem norm_partialEval_sigma_le (k : ℕ) (q : ℂ)
    (hq : (2 ^ (k + 1) : ℝ) * ‖q‖ ≤ 1 / 2) (n : ℕ) :
    ‖partialEval (sigmaSeries k) q (n + 1)‖ ≤ (n : ℝ) * 1 := by
  rw [partialEval_sigmaSeries_k_succ]
  refine norm_partialSum_le _ 1 (fun N => ?_) n
  have h := norm_sigma_k_term_le k q N
  have hbase : (0 : ℝ) ≤ (2 ^ (k + 1) : ℝ) * ‖q‖ := by positivity
  have hpow : ((2 ^ (k + 1) : ℝ) * ‖q‖) ^ (N + 1) ≤ 1 :=
    pow_le_one₀ hbase (by linarith)
  linarith

theorem norm_partialEval_sigma1_le (q : ℂ) (hq : ‖q‖ ≤ 1 / 8) (n : ℕ) :
    ‖partialEval (sigmaSeries 1) q (n + 1)‖ ≤ (n : ℝ) * 1 := by
  refine norm_partialEval_sigma_le 1 q ?_ n
  have h4 : ((2 : ℝ) ^ (1 + 1)) = 4 := by norm_num
  rw [h4]
  linarith

/-! ## ★★★★★切り詰めの大きさ -/

/-- ★★★★★**`X` の切り詰めの大きさ**——`u` に依存するのは `f(u)` だけ。 -/
theorem norm_tateXtrunc_le (u w : ℂ) (hu : ‖u‖ ≤ 1) (hw : ‖w‖ ≤ 1 / 2)
    (hq : ‖u * w‖ ≤ 1 / 8) (n : ℕ) :
    ‖tateXtrunc (n + 1) u w (u * w)‖ ≤ ‖tateXterm u‖ + 10 * ((n : ℝ) + 1) + 2 := by
  have hq2 : ‖u * w‖ ≤ 1 / 2 := by linarith
  have hw1 : ‖w‖ ≤ 1 := by linarith
  have hs1 : ‖partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * u)) (n + 1)‖
      ≤ ((n : ℝ) + 1) * 4 := by
    have h := norm_partialSum_le (fun m => tateXterm ((u * w) ^ (m + 1) * u)) 4
      (fun m => norm_tateXterm_pow_le (u * w) u hq2 hu m) (n + 1)
    push_cast at h
    linarith
  have hs2 : ‖partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1)‖
      ≤ ((n : ℝ) + 1) * 4 := by
    have h := norm_partialSum_le (fun m => tateXterm ((u * w) ^ (m + 1) * w)) 4
      (fun m => norm_tateXterm_pow_le (u * w) w hq2 hw1 m) (n + 1)
    push_cast at h
    linarith
  have hs3 := norm_partialEval_sigma1_le (u * w) hq n
  have hfw : ‖tateXterm w‖ ≤ 2 := by
    have h := norm_tateXterm_le_of_small hw
    linarith
  rw [tateXtrunc]
  have hb := norm_sum7 (tateXterm u)
    (partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * u)) (n + 1))
    (tateXterm w) (partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
    ((-2) * partialEval (sigmaSeries 1) (u * w) (n + 1)) 0 0
  have heq : (tateXterm u + partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * u)) (n + 1))
      + (tateXterm w + partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
      - 2 * partialEval (sigmaSeries 1) (u * w) (n + 1)
      = tateXterm u + partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * u)) (n + 1)
        + tateXterm w + partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1)
        + ((-2) * partialEval (sigmaSeries 1) (u * w) (n + 1)) + 0 + 0 := by ring
  rw [heq]
  have hlast : ‖(-2 : ℂ) * partialEval (sigmaSeries 1) (u * w) (n + 1)‖
      ≤ 2 * ((n : ℝ) * 1) := by
    rw [norm_mul]
    have h2 : ‖(-2 : ℂ)‖ = 2 := by simp
    rw [h2]
    linarith
  simp only [norm_zero] at hb
  linarith

/-- ★★★★★**`Y` の切り詰めの大きさ**。 -/
theorem norm_tateYtrunc_le (u w : ℂ) (hu : ‖u‖ ≤ 1) (hw : ‖w‖ ≤ 1 / 2)
    (hq : ‖u * w‖ ≤ 1 / 8) (n : ℕ) :
    ‖tateYtrunc (n + 1) u w (u * w)‖ ≤ ‖tateYterm u‖ + 25 * ((n : ℝ) + 1) + 6 := by
  have hq2 : ‖u * w‖ ≤ 1 / 2 := by linarith
  have hw1 : ‖w‖ ≤ 1 := by linarith
  have hs1 : ‖partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * u)) (n + 1)‖
      ≤ ((n : ℝ) + 1) * 8 := by
    have h := norm_partialSum_le (fun m => tateYterm ((u * w) ^ (m + 1) * u)) 8
      (fun m => norm_tateYterm_pow_le (u * w) u hq2 hu m) (n + 1)
    push_cast at h
    linarith
  have hs2 : ‖partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1)‖
      ≤ ((n : ℝ) + 1) * 4 := by
    have h := norm_partialSum_le (fun m => tateXterm ((u * w) ^ (m + 1) * w)) 4
      (fun m => norm_tateXterm_pow_le (u * w) w hq2 hw1 m) (n + 1)
    push_cast at h
    linarith
  have hs3 : ‖partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * w)) (n + 1)‖
      ≤ ((n : ℝ) + 1) * 8 := by
    have h := norm_partialSum_le (fun m => tateYterm ((u * w) ^ (m + 1) * w)) 8
      (fun m => norm_tateYterm_pow_le (u * w) w hq2 hw1 m) (n + 1)
    push_cast at h
    linarith
  have hs4 := norm_partialEval_sigma1_le (u * w) hq n
  have hfw : ‖tateXterm w‖ ≤ 2 := by
    have h := norm_tateXterm_le_of_small hw
    linarith
  have hgw : ‖tateYterm w‖ ≤ 4 := by
    have h := norm_tateYterm_le_of_small hw
    linarith
  rw [tateYtrunc]
  have hb := norm_sum7 (tateYterm u)
    (partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * u)) (n + 1))
    (-tateXterm w) (-partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
    (-tateYterm w) (-partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * w)) (n + 1))
    (partialEval (sigmaSeries 1) (u * w) (n + 1))
  have heq : (tateYterm u + partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * u)) (n + 1))
      - (tateXterm w + partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
      - (tateYterm w + partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * w)) (n + 1))
      + partialEval (sigmaSeries 1) (u * w) (n + 1)
      = tateYterm u + partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * u)) (n + 1)
        + (-tateXterm w) + (-partialSum (fun m => tateXterm ((u * w) ^ (m + 1) * w)) (n + 1))
        + (-tateYterm w) + (-partialSum (fun m => tateYterm ((u * w) ^ (m + 1) * w)) (n + 1))
        + partialEval (sigmaSeries 1) (u * w) (n + 1) := by ring
  rw [heq]
  simp only [norm_neg] at hb
  linarith

/-! ## ★★★★係数の切り詰めの大きさ -/

theorem norm_partialEval_tateA4_le (q : ℂ) (hq : ‖q‖ ≤ 1 / 128) (n : ℕ) :
    ‖partialEval tateA4 q (n + 1)‖ ≤ 5 * (n : ℝ) := by
  have h3 : ‖partialEval (sigmaSeries 3) q (n + 1)‖ ≤ (n : ℝ) * 1 := by
    refine norm_partialEval_sigma_le 3 q ?_ n
    have h16 : ((2 : ℝ) ^ (3 + 1)) = 16 := by norm_num
    rw [h16]
    linarith
  rw [partialEval_tateA4, norm_mul]
  have h5 : ‖(-5 : ℂ)‖ = 5 := by simp
  rw [h5]
  linarith

theorem norm_partialEval_tateA6_le (q : ℂ) (hq : ‖q‖ ≤ 1 / 128) (n : ℕ) :
    ‖partialEval tateA6 q (n + 1)‖ ≤ (n : ℝ) := by
  have h3 : ‖partialEval (sigmaSeries 3) q (n + 1)‖ ≤ (n : ℝ) * 1 := by
    refine norm_partialEval_sigma_le 3 q ?_ n
    have h16 : ((2 : ℝ) ^ (3 + 1)) = 16 := by norm_num
    rw [h16]
    linarith
  have h5 : ‖partialEval (sigmaSeries 5) q (n + 1)‖ ≤ (n : ℝ) * 1 := by
    refine norm_partialEval_sigma_le 5 q ?_ n
    have h64 : ((2 : ℝ) ^ (5 + 1)) = 64 := by norm_num
    rw [h64]
    linarith
  have hpe := partialEval_tateA6 q (n + 1)
  have hkey : ‖(12 : ℂ) * partialEval tateA6 q (n + 1)‖ ≤ 12 * (n : ℝ) := by
    rw [hpe, norm_neg]
    have hb := norm_add_le ((5 : ℂ) * partialEval (sigmaSeries 3) q (n + 1))
      ((7 : ℂ) * partialEval (sigmaSeries 5) q (n + 1))
    rw [norm_mul, norm_mul] at hb
    have e5 : ‖(5 : ℂ)‖ = 5 := by simp
    have e7 : ‖(7 : ℂ)‖ = 7 := by simp
    rw [e5, e7] at hb
    linarith
  rw [norm_mul] at hkey
  have e12 : ‖(12 : ℂ)‖ = 12 := by simp
  rw [e12] at hkey
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def norm_defect_diff_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——方程式の差の差)",
    sectionId := "genell-def-3-3" }

def norm_tateXtrunc_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——切り詰めの大きさ)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
