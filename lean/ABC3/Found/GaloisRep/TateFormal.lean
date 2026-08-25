import ABC3.Found.GaloisRep.TateOrigin

/-!
# Galois (G6) 第 276 ブロック —— **★★★★★★★★★原点近傍の単射性と全射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点——葉 (e) の 3 つ目の領域

> `1 − u ∈ I` の領域で **`z = (u−1)·XE/YE` が `u` と 1 対 1**
> (`tateZ_inj`・`exists_tateZ_eq`)

★★★`z` は **`−X_K/Y_K`**、すなわち形式群の座標そのものである(`tateZ_map`)。

## ★★★★★★★★母数の選び方

| 領域 | 母数 | 主要項 | 出典 |
|---|---|---|---|
| 単元 | `Λ = Y/(X+Y)` | `Λ ≡ a` | 第 263 |
| 環帯 | `(a,w)` の 2 変数 | `X ≡ a+w`、`Y ≡ −w` | 第 269 |
| 原点近傍 | `z = (u−1)·XE/YE` | `z ≡ u−1` | **本ブロック** |

★★★どれも「主要項の微分が 1」になるように選んである。それが縮小写像定理の
入り口である。原点近傍では **`XE ≡ u`、`YE ≡ u²`** なので
`XE/YE ≡ 1/u ≡ 1` となり、`z ≡ u−1` が出る。

## ★★★★★`XE − YE` に `1 − u` が括り出せる

    XE − YE = (u − u²) + (1−u)²A − (1−u)³B = (1−u)·[u + (1−u)A − (1−u)²B]

★これで `κ = z − (u−1) = (u−1)(XE−YE)·YE⁻¹` が **`I²` の元**になり、
差の評価が 1 つ位を上げる(`tateZ_diff_mem`)。

★★`YE` が単元であること(`YE = u² + I` の元)も同じ観察から出る。

## ★★葉 (e) の領域がそろった

| 領域 | 全射性 | 単射性 |
|---|---|---|
| 単元 | 第 266 | 第 267 |
| 環帯 | 第 271 | 第 259 |
| 原点近傍 | **本ブロック** | **本ブロック** |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXpairE_sub_tateYpairE_mem` | ★★★★★`XE − YE ∈ I` |
| `isUnit_tateYpairE` | ★★★★`YE` は単元 |
| `tateXpairE_diff_mem`・`tateYpairE_diff_mem` | ★★★★★★`u` についての差 |
| `tateZ`・`tateZ_sub_mem` | ★★★★★★★原点近傍の母数 |
| `tateZ_diff_mem` | ★★★★★★★★**原点近傍の縮小性** |
| `tateZ_inj` | ★★★★★★★★**原点近傍での単射性** |
| `exists_tateZ_eq` | ★★★★★★★★★**原点近傍での全射性** |
| `tateZ_map` | ★★★★★★`z = −X_K/Y_K` |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★`XE − YE` に `1 − u` が括り出せる -/

theorem one_sub_form_mem {u A B : R} (hu : 1 - u ∈ I) :
    u + (1 - u) ^ 2 * A - (u ^ 2 + (1 - u) ^ 3 * B) ∈ I := by
  have hexp : u + (1 - u) ^ 2 * A - (u ^ 2 + (1 - u) ^ 3 * B)
      = (1 - u) * (u + (1 - u) * A - (1 - u) ^ 2 * B) := by ring
  rw [hexp]
  exact Ideal.mul_mem_right _ _ hu

theorem isUnit_of_one_sub_mem [IsAdicComplete I R] {u : R} (hu : 1 - u ∈ I) : IsUnit u := by
  simpa using isUnit_one_sub hu

/-- ★★★★★**`XE − YE ∈ I`**——`1 − u` が括り出せるから。 -/
theorem tateXpairE_sub_tateYpairE_mem [IsAdicComplete I R] (u w q : R) (hq : q ∈ I)
    (hu : 1 - u ∈ I) : tateXpairE u w q hq - tateYpairE u w q hq ∈ I := by
  rw [tateXpairE, tateYpairE]
  exact one_sub_form_mem hu

/-- ★★★★**`YE` は単元**——`u²` に `I` の元を足したもの。 -/
theorem isUnit_tateYpairE [IsAdicComplete I R] (u w q : R) (hq : q ∈ I) (hu : 1 - u ∈ I) :
    IsUnit (tateYpairE u w q hq) := by
  have huU : IsUnit u := isUnit_of_one_sub_mem hu
  have hmem : (1 - u) ^ 3 * (tateYtail u q hq - (tateXterm w + tateXtail w q hq)
      - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq) ∈ I := by
    refine Ideal.mul_mem_right _ _ ?_
    rw [pow_succ]
    exact Ideal.mul_mem_left _ _ hu
  rw [tateYpairE]
  exact isUnit_add_mem (huU.pow 2) hmem

/-! ## ★★★★★★`u` についての差 -/

theorem tateXpairE_diff_mem [IsAdicComplete I R] {k : ℕ} (u u' w w' q : R) (hq : q ∈ I)
    (hw : w ∈ I) (hw' : w' ∈ I) (huu : u - u' ∈ I ^ k) (hww : w - w' ∈ I ^ k) :
    tateXpairE u w q hq - tateXpairE u' w' q hq ∈ I ^ k := by
  have hq1 : q ∈ I ^ 1 := by simpa using hq
  have huu' : u' - u ∈ I ^ k := by simpa using neg_mem huu
  have hT : (tateXtail u q hq + (tateXterm w + tateXtail w q hq)
        - 2 * evalAdic (sigmaSeries 1) q hq)
      - (tateXtail u' q hq + (tateXterm w' + tateXtail w' q hq)
        - 2 * evalAdic (sigmaSeries 1) q hq) ∈ I ^ k := by
    have h1 : tateXtail u q hq - tateXtail u' q hq ∈ I ^ k :=
      Ideal.pow_le_pow_right (by omega) (tateXtail_diff_mem (m := 1) u u' q hq hq1 huu)
    have h2 : tateXterm w - tateXterm w' ∈ I ^ k := tateXterm_diff_mem' hw hw' hww
    have h3 : tateXtail w q hq - tateXtail w' q hq ∈ I ^ k :=
      Ideal.pow_le_pow_right (by omega) (tateXtail_diff_mem (m := 1) w w' q hq hq1 hww)
    have hexp : (tateXtail u q hq + (tateXterm w + tateXtail w q hq)
          - 2 * evalAdic (sigmaSeries 1) q hq)
        - (tateXtail u' q hq + (tateXterm w' + tateXtail w' q hq)
          - 2 * evalAdic (sigmaSeries 1) q hq)
        = (tateXtail u q hq - tateXtail u' q hq) + (tateXterm w - tateXterm w')
          + (tateXtail w q hq - tateXtail w' q hq) := by ring
    rw [hexp]
    exact Ideal.add_mem _ (Ideal.add_mem _ h1 h2) h3
  rw [tateXpairE, tateXpairE]
  have hexp : (u + (1 - u) ^ 2 * (tateXtail u q hq + (tateXterm w + tateXtail w q hq)
        - 2 * evalAdic (sigmaSeries 1) q hq))
      - (u' + (1 - u') ^ 2 * (tateXtail u' q hq + (tateXterm w' + tateXtail w' q hq)
        - 2 * evalAdic (sigmaSeries 1) q hq))
      = (u - u') + ((u' - u) * ((1 - u) + (1 - u')))
          * (tateXtail u q hq + (tateXterm w + tateXtail w q hq)
            - 2 * evalAdic (sigmaSeries 1) q hq)
        + (1 - u') ^ 2 * ((tateXtail u q hq + (tateXterm w + tateXtail w q hq)
            - 2 * evalAdic (sigmaSeries 1) q hq)
          - (tateXtail u' q hq + (tateXterm w' + tateXtail w' q hq)
            - 2 * evalAdic (sigmaSeries 1) q hq)) := by ring
  rw [hexp]
  exact Ideal.add_mem _ (Ideal.add_mem _ huu
    (Ideal.mul_mem_right _ _ (Ideal.mul_mem_right _ _ huu')))
    (Ideal.mul_mem_left _ _ hT)

theorem tateYpairE_diff_mem [IsAdicComplete I R] {k : ℕ} (u u' w w' q : R) (hq : q ∈ I)
    (hw : w ∈ I) (hw' : w' ∈ I) (huu : u - u' ∈ I ^ k) (hww : w - w' ∈ I ^ k) :
    tateYpairE u w q hq - tateYpairE u' w' q hq ∈ I ^ k := by
  have hq1 : q ∈ I ^ 1 := by simpa using hq
  have huu' : u' - u ∈ I ^ k := by simpa using neg_mem huu
  have hsq : u ^ 2 - u' ^ 2 ∈ I ^ k := by
    have h : u ^ 2 - u' ^ 2 = (u - u') * (u + u') := by ring
    rw [h]
    exact Ideal.mul_mem_right _ _ huu
  have hT : (tateYtail u q hq - (tateXterm w + tateXtail w q hq)
        - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq)
      - (tateYtail u' q hq - (tateXterm w' + tateXtail w' q hq)
        - (tateYterm w' + tateYtail w' q hq) + evalAdic (sigmaSeries 1) q hq) ∈ I ^ k := by
    have h1 : tateYtail u q hq - tateYtail u' q hq ∈ I ^ k :=
      Ideal.pow_le_pow_right (by omega) (tateYtail_diff_mem (m := 1) u u' q hq hq1 huu)
    have h2 : tateXterm w - tateXterm w' ∈ I ^ k := tateXterm_diff_mem' hw hw' hww
    have h3 : tateXtail w q hq - tateXtail w' q hq ∈ I ^ k :=
      Ideal.pow_le_pow_right (by omega) (tateXtail_diff_mem (m := 1) w w' q hq hq1 hww)
    have h4 : tateYterm w - tateYterm w' ∈ I ^ k := tateYterm_diff_mem' hw hw' hww
    have h5 : tateYtail w q hq - tateYtail w' q hq ∈ I ^ k :=
      Ideal.pow_le_pow_right (by omega) (tateYtail_diff_mem (m := 1) w w' q hq hq1 hww)
    have hexp : (tateYtail u q hq - (tateXterm w + tateXtail w q hq)
          - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq)
        - (tateYtail u' q hq - (tateXterm w' + tateXtail w' q hq)
          - (tateYterm w' + tateYtail w' q hq) + evalAdic (sigmaSeries 1) q hq)
        = (tateYtail u q hq - tateYtail u' q hq) - (tateXterm w - tateXterm w')
          - (tateXtail w q hq - tateXtail w' q hq) - (tateYterm w - tateYterm w')
          - (tateYtail w q hq - tateYtail w' q hq) := by ring
    rw [hexp]
    exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ h1 h2) h3) h4) h5
  rw [tateYpairE, tateYpairE]
  have hexp : (u ^ 2 + (1 - u) ^ 3 * (tateYtail u q hq - (tateXterm w + tateXtail w q hq)
        - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq))
      - (u' ^ 2 + (1 - u') ^ 3 * (tateYtail u' q hq - (tateXterm w' + tateXtail w' q hq)
        - (tateYterm w' + tateYtail w' q hq) + evalAdic (sigmaSeries 1) q hq))
      = (u ^ 2 - u' ^ 2) + ((u' - u) * ((1 - u) ^ 2 + (1 - u) * (1 - u') + (1 - u') ^ 2))
          * (tateYtail u q hq - (tateXterm w + tateXtail w q hq)
            - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq)
        + (1 - u') ^ 3 * ((tateYtail u q hq - (tateXterm w + tateXtail w q hq)
            - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq)
          - (tateYtail u' q hq - (tateXterm w' + tateXtail w' q hq)
            - (tateYterm w' + tateYtail w' q hq) + evalAdic (sigmaSeries 1) q hq)) := by ring
  rw [hexp]
  exact Ideal.add_mem _ (Ideal.add_mem _ hsq
    (Ideal.mul_mem_right _ _ (Ideal.mul_mem_right _ _ huu')))
    (Ideal.mul_mem_left _ _ hT)

/-! ## ★★★★★★★原点近傍の母数 -/

/-- ★★★★★★★**原点近傍の母数** `z = (u−1)·XE/YE`——形式群の座標 `−x/y` にあたる。 -/
noncomputable def tateZ [IsAdicComplete I R] (u w q : R) (hq : q ∈ I) : R :=
  (u - 1) * tateXpairE u w q hq * Ring.inverse (tateYpairE u w q hq)

theorem tateZ_sub_eq [IsAdicComplete I R] (u w q : R) (hq : q ∈ I) (hu : 1 - u ∈ I) :
    tateZ u w q hq - (u - 1)
      = (u - 1) * (tateXpairE u w q hq - tateYpairE u w q hq)
        * Ring.inverse (tateYpairE u w q hq) := by
  have hY := isUnit_tateYpairE u w q hq hu
  have h : tateYpairE u w q hq * Ring.inverse (tateYpairE u w q hq) = 1 :=
    Ring.mul_inverse_cancel _ hY
  rw [tateZ]
  linear_combination (u - 1) * h

/-- ★★★★★**`z ≡ u − 1` は 2 次の精度**。 -/
theorem tateZ_sub_mem [IsAdicComplete I R] (u w q : R) (hq : q ∈ I) (hu : 1 - u ∈ I) :
    tateZ u w q hq - (u - 1) ∈ I ^ 2 := by
  have hu1 : u - 1 ∈ I := by simpa using neg_mem hu
  have hQ := tateXpairE_sub_tateYpairE_mem u w q hq hu
  rw [tateZ_sub_eq u w q hq hu, pow_two]
  exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_mul hu1 hQ)

/-- ★★★★★★★★**原点近傍の縮小性**——`κ = z − (u−1)` は位を 1 つ上げる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateZ_diff_mem [IsAdicComplete I R] {k : ℕ} (u u' w w' q : R) (hq : q ∈ I)
    (hu : 1 - u ∈ I) (hu' : 1 - u' ∈ I) (hw : w ∈ I) (hw' : w' ∈ I)
    (huu : u - u' ∈ I ^ k) (hww : w - w' ∈ I ^ k) :
    (tateZ u w q hq - (u - 1)) - (tateZ u' w' q hq - (u' - 1)) ∈ I ^ (k + 1) := by
  have hu1' : u' - 1 ∈ I := by simpa using neg_mem hu'
  have hQ := tateXpairE_sub_tateYpairE_mem u w q hq hu
  have hdX := tateXpairE_diff_mem u u' w w' q hq hw hw' huu hww
  have hdY := tateYpairE_diff_mem u u' w w' q hq hw hw' huu hww
  have hdQ : (tateXpairE u w q hq - tateYpairE u w q hq)
      - (tateXpairE u' w' q hq - tateYpairE u' w' q hq) ∈ I ^ k := by
    have hexp : (tateXpairE u w q hq - tateYpairE u w q hq)
        - (tateXpairE u' w' q hq - tateYpairE u' w' q hq)
        = (tateXpairE u w q hq - tateXpairE u' w' q hq)
          - (tateYpairE u w q hq - tateYpairE u' w' q hq) := by ring
    rw [hexp]
    exact Ideal.sub_mem _ hdX hdY
  have hdS : Ring.inverse (tateYpairE u w q hq) - Ring.inverse (tateYpairE u' w' q hq) ∈ I ^ k :=
    inverse_diff_mem (isUnit_tateYpairE u w q hq hu) (isUnit_tateYpairE u' w' q hq hu') hdY
  rw [tateZ_sub_eq u w q hq hu, tateZ_sub_eq u' w' q hq hu']
  set P := u - 1
  set P' := u' - 1
  set Q := tateXpairE u w q hq - tateYpairE u w q hq
  set Q' := tateXpairE u' w' q hq - tateYpairE u' w' q hq
  set S := Ring.inverse (tateYpairE u w q hq)
  set S' := Ring.inverse (tateYpairE u' w' q hq)
  have hPP : P - P' ∈ I ^ k := by simpa [P, P'] using huu
  have hexp : P * Q * S - P' * Q' * S'
      = ((P - P') * Q) * S + ((Q - Q') * P') * S + ((S - S') * P') * Q' := by ring
  rw [hexp, pow_succ]
  exact Ideal.add_mem _ (Ideal.add_mem _
    (Ideal.mul_mem_right _ _ (Ideal.mul_mem_mul hPP hQ))
    (Ideal.mul_mem_right _ _ (Ideal.mul_mem_mul hdQ hu1')))
    (Ideal.mul_mem_right _ _ (Ideal.mul_mem_mul hdS hu1'))

/-! ## ★★★★★★★★★原点近傍の単射性と全射性 -/

/-- ★★★★★★★★**原点近傍での単射性**——母数 `z` が `u` を決める。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateZ_inj [IsAdicComplete I R] (u u' q : R) (hq : q ∈ I)
    (hu : 1 - u ∈ I) (hu' : 1 - u' ∈ I)
    (hZ : tateZ u (wOf q u) q hq = tateZ u' (wOf q u') q hq) : u = u' := by
  have hUu : IsUnit u := isUnit_of_one_sub_mem hu
  have hUu' : IsUnit u' := isUnit_of_one_sub_mem hu'
  have hstep : ∀ k : ℕ, u - u' ∈ I ^ k := by
    intro k
    induction k with
    | zero => simp
    | succ n ih =>
      have hww : wOf q u - wOf q u' ∈ I ^ n :=
        Ideal.pow_le_pow_right (by omega) (wOf_diff_mem hq hUu hUu' ih)
      have h := tateZ_diff_mem u u' (wOf q u) (wOf q u') q hq hu hu'
        (wOf_mem u hq) (wOf_mem u' hq) ih hww
      have heq : (tateZ u (wOf q u) q hq - (u - 1)) - (tateZ u' (wOf q u') q hq - (u' - 1))
          = -(u - u') := by rw [hZ]; ring
      rw [heq] at h
      simpa using neg_mem h
  exact sub_eq_zero.1 (eq_zero_of_mem_pow (I := I) hstep)

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★★**原点近傍での全射性**——どの `z ∈ I` にも `u` が在る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_tateZ_eq [IsAdicComplete I R] (z q : R) (hz : z ∈ I) (hq : q ∈ I) :
    ∃ u : R, 1 - u ∈ I ∧ tateZ u (wOf q u) q hq = z := by
  set F : R → R := fun s => z - (tateZ (1 + s) (wOf q (1 + s)) q hq - ((1 + s) - 1)) with hF
  have hone : ∀ s : R, s ∈ I → 1 - (1 + s) ∈ I := by
    intro s hs
    have h : 1 - (1 + s) = -s := by ring
    rw [h]
    exact neg_mem hs
  have hFI : ∀ s ∈ I, F s ∈ I := by
    intro s hs
    refine Ideal.sub_mem _ hz ?_
    have h2 := tateZ_sub_mem (1 + s) (wOf q (1 + s)) q hq (hone s hs)
    simpa using Ideal.pow_le_pow_right (I := I) (show (1:ℕ) ≤ 2 by omega) h2
  have hcon : ∀ x ∈ I, ∀ y ∈ I, ∀ k : ℕ, x - y ∈ I ^ k → F x - F y ∈ I ^ (k + 1) := by
    intro x hx y hy k hxy
    have hUx : IsUnit (1 + x) := isUnit_of_one_sub_mem (hone x hx)
    have hUy : IsUnit (1 + y) := isUnit_of_one_sub_mem (hone y hy)
    have hxy' : (1 + x) - (1 + y) ∈ I ^ k := by
      have h : (1 + x) - (1 + y) = x - y := by ring
      rw [h]; exact hxy
    have hww : wOf q (1 + x) - wOf q (1 + y) ∈ I ^ k :=
      Ideal.pow_le_pow_right (by omega) (wOf_diff_mem hq hUx hUy hxy')
    have h := tateZ_diff_mem (1 + x) (1 + y) (wOf q (1 + x)) (wOf q (1 + y)) q hq
      (hone x hx) (hone y hy) (wOf_mem (1 + x) hq) (wOf_mem (1 + y) hq) hxy' hww
    have heq : F x - F y
        = -((tateZ (1 + x) (wOf q (1 + x)) q hq - ((1 + x) - 1))
          - (tateZ (1 + y) (wOf q (1 + y)) q hq - ((1 + y) - 1))) := by
      rw [hF]; ring
    rw [heq]
    exact neg_mem h
  obtain ⟨s, hs, hfix⟩ :=
    exists_fixedPoint_of_contraction F hFI (hFI 0 (Submodule.zero_mem I)) hcon
  refine ⟨1 + s, hone s hs, ?_⟩
  rw [hF] at hfix
  linear_combination -hfix

/-! ## ★★★★★★`z = −X_K/Y_K` -/

section Field

variable {K : Type} [IsAdicComplete I R] [Field K] [Algebra R K]

/-- ★★★★★★**`z = −X_K/Y_K`**——形式群の座標そのもの。 -/
theorem tateZ_map (u w q : R) (hq : q ∈ I) (hu : 1 - u ∈ I)
    (hne : algebraMap R K (1 - u) ≠ 0) :
    algebraMap R K (tateZ u w q hq)
      = -(tateXK u w q hq / tateYK (K := K) u w q hq) := by
  have hYE := isUnit_tateYpairE u w q hq hu
  have hYne : algebraMap R K (tateYpairE u w q hq) ≠ 0 := (hYE.map (algebraMap R K)).ne_zero
  have hinv : algebraMap R K (Ring.inverse (tateYpairE u w q hq))
      = (algebraMap R K (tateYpairE u w q hq))⁻¹ := by
    have h2 := congrArg (algebraMap R K) (Ring.mul_inverse_cancel _ hYE)
    rw [map_mul, map_one] at h2
    exact (inv_eq_of_mul_eq_one_right h2).symm
  have hu1 : algebraMap R K (u - 1) = -algebraMap R K (1 - u) := by
    rw [← map_neg]
    congr 1
    ring
  rw [tateZ, map_mul, map_mul, hinv, hu1, tateXK, tateYK]
  field_simp

end Field

/-! ## ★出典の紐付け(`.src`) -/

def tateZ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍の母数)",
    sectionId := "genell-def-3-3" }

def tateZ_inj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍での単射性)",
    sectionId := "genell-def-3-3" }

def exists_tateZ_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍での全射性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
