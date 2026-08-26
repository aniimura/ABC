import ABC3.Found.GaloisRep.TateComplexSpec

/-!
# Galois (G6) 第 226 ブロック —— **★★★★★★ℤ 和を頭と尾に分ける**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★解析側と形式側の形を揃える

解析側(第 218・219)の `X`・`Y` は **`n ∈ ℤ` にわたる和**で書かれている。
形式側(第 105-107)の `tateXpair`・`tateYpair` は **頭 + `ℕ` の尾**で書かれている。
本ブロックで**同じものであることを型で確かめる**。

`q = uw` とすると `q^{−(m+1)}·u = (q^m·w)⁻¹` である(`zpow_neg_succ_mul`)。よって:

| 段 `n` | `X` 側 | `Y` 側 |
|---|---|---|
| `n = 0` | `f(u)` | `g(u)` |
| `n = −(m+1)` | `f(q^{m+1}u)` | `g(q^{m+1}u)` |
| `n = m+1` | `f((q^m w)⁻¹) = f(q^m w)` | `g((q^m w)⁻¹) = −(f(q^m w) + g(q^m w))` |

★★`f` は反転で不変、`g` は `g(1/t) = −(f(t)+g(t))`(`tateYterm_inv'`)——
**偶数べき・奇数べきの違いがそのまま `X` と `Y` の形の違いになる**。

    ∑_{n∈ℤ} f(q^{−n}u) = [f(u) + ∑_{m≥1}f(qᵐu)] + [f(w) + ∑_{m≥1}f(qᵐw)]
    ∑_{n∈ℤ} g(q^{−n}u) = [g(u) + ∑_{m≥1}g(qᵐu)] − [f(w) + ∑_{m≥1}f(qᵐw)]
                                                  − [g(w) + ∑_{m≥1}g(qᵐw)]

★★★これは `tateXpair`・`tateYpair` の定義そのものの形である。

## ★★収束の条件は `‖u‖ < 1`、`‖w‖ < 1` だけ

`q = uw` なので `‖q‖ < 1`。`‖q^m w‖ ≤ ‖w‖ < 1` なので、反転則に要る `t ≠ 1` も自動で出る。
★`t ≠ 0` には `u ≠ 0`、`w ≠ 0` が要る(`‖·‖ < 1` からは出ない)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `zpow_neg_succ_mul` | ★`q^{−(m+1)}u = (q^m w)⁻¹` |
| `norm_mul_lt_one` / `norm_pow_mul_lt_one` | ★ノルムの評価 |
| `summable_nat_tateXterm` 他 | ★★`ℕ` 側の可和性 |
| `tateYterm_inv'` | ★★★★`g(1/t) = −(f(t)+g(t))` |
| `tsum_int_tateXterm_split` | ★★★★★★**`X` の ℤ 和は頭と二本の尾** |
| `tsum_int_tateYterm_split` | ★★★★★★**`Y` の ℤ 和は頭と三本の尾** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★指数と収束の下ごしらえ -/

/-- ★`q^{−(m+1)}·u = (q^m·w)⁻¹`(`q = uw`)。 -/
theorem zpow_neg_succ_mul (u w : ℂ) (hu0 : u ≠ 0) (hw0 : w ≠ 0) (m : ℕ) :
    (u * w) ^ (-((m : ℤ) + 1)) * u = ((u * w) ^ m * w)⁻¹ := by
  have huw : u * w ≠ 0 := mul_ne_zero hu0 hw0
  rw [zpow_neg, ← zpow_natCast (u * w) m]
  field_simp
  rw [zpow_add_one₀ huw]
  ring

theorem norm_mul_lt_one {u w : ℂ} (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) : ‖u * w‖ < 1 := by
  rw [norm_mul]
  nlinarith [norm_nonneg u, norm_nonneg w]

theorem norm_pow_mul_lt_one {q c : ℂ} (hq : ‖q‖ < 1) (hc : ‖c‖ < 1) (m : ℕ) :
    ‖q ^ m * c‖ < 1 := by
  rw [norm_mul, norm_pow]
  have h1 : ‖q‖ ^ m ≤ 1 := pow_le_one₀ (norm_nonneg q) hq.le
  nlinarith [norm_nonneg c, pow_nonneg (norm_nonneg q) m]

theorem ne_one_of_norm_lt_one {t : ℂ} (h : ‖t‖ < 1) : t ≠ 1 := by
  intro hc
  rw [hc] at h
  simp at h

/-! ## ★★`ℕ` 側の可和性 -/

theorem summable_nat_tateXterm (q c : ℂ) (hq : ‖q‖ < 1) :
    Summable fun m : ℕ => tateXterm (q ^ m * c) := by
  refine (summable_tateXterm_small q c hq).congr fun k => ?_
  rw [zpow_natCast]

theorem summable_nat_tateXterm_succ (q c : ℂ) (hq : ‖q‖ < 1) :
    Summable fun m : ℕ => tateXterm (q ^ (m + 1) * c) :=
  (summable_nat_add_iff 1).2 (summable_nat_tateXterm q c hq)

theorem summable_nat_tateDterm (q c : ℂ) (hq : ‖q‖ < 1) :
    Summable fun m : ℕ => tateDterm (q ^ m * c) := by
  refine (summable_tateDterm_small q c hq).congr fun k => ?_
  rw [zpow_natCast]

theorem summable_nat_tateYterm (q c : ℂ) (hq : ‖q‖ < 1) :
    Summable fun m : ℕ => tateYterm (q ^ m * c) := by
  have h := ((summable_nat_tateDterm q c hq).sub (summable_nat_tateXterm q c hq)).mul_right
    (1 / 2 : ℂ)
  refine h.congr fun m => ?_
  rw [tateYterm_eq_half]
  ring

theorem summable_nat_tateYterm_succ (q c : ℂ) (hq : ‖q‖ < 1) :
    Summable fun m : ℕ => tateYterm (q ^ (m + 1) * c) :=
  (summable_nat_add_iff 1).2 (summable_nat_tateYterm q c hq)

/-! ## ★★★★反転則(`Y` の側) -/

/-- ★★★★**`g(1/t) = −(f(t) + g(t))`**——`Y` の段が `X` と `Y` の両方を生む理由。 -/
theorem tateYterm_inv' {t : ℂ} (ht0 : t ≠ 0) (ht1 : t ≠ 1) :
    tateYterm t⁻¹ = -(tateXterm t + tateYterm t) := by
  have h1 : (1 : ℂ) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht1)
  rw [tateYterm_inv ht0 ht1, tateXterm_field, tateYterm_field]
  field_simp
  ring

/-! ## ★★★★★★ℤ 和を頭と尾に分ける -/

/-- ★★★★★★**`X` の ℤ 和を頭と二本の尾に分ける**——形式側 `tateXpair` と同じ形。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tsum_int_tateXterm_split (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1)
    (hu0 : u ≠ 0) (hw0 : w ≠ 0) :
    (∑' n : ℤ, tateXterm ((u * w) ^ (-n) * u))
      = (tateXterm u + ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * u))
        + (tateXterm w + ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w)) := by
  have hq : ‖u * w‖ < 1 := norm_mul_lt_one hu hw
  have hq0 : u * w ≠ 0 := mul_ne_zero hu0 hw0
  have hpos : ∀ m : ℕ,
      tateXterm ((u * w) ^ (-((m : ℤ) + 1)) * u) = tateXterm ((u * w) ^ m * w) := by
    intro m
    have hne0 : (u * w) ^ m * w ≠ 0 := mul_ne_zero (pow_ne_zero _ hq0) hw0
    have hne1 : (u * w) ^ m * w ≠ 1 := ne_one_of_norm_lt_one (norm_pow_mul_lt_one hq hw m)
    rw [zpow_neg_succ_mul u w hu0 hw0 m, tateXterm_inv hne0 hne1]
  have hneg : ∀ m : ℕ,
      tateXterm ((u * w) ^ (-(-((m : ℤ) + 1))) * u) = tateXterm ((u * w) ^ (m + 1) * u) := by
    intro m
    rw [neg_neg, show ((m : ℤ) + 1) = ((m + 1 : ℕ) : ℤ) by push_cast; ring, zpow_natCast]
  have h1 : Summable fun m : ℕ => tateXterm ((u * w) ^ (-((m : ℤ) + 1)) * u) := by
    refine (summable_nat_tateXterm (u * w) w hq).congr fun m => ?_
    rw [hpos m]
  have h2 : Summable fun m : ℕ => tateXterm ((u * w) ^ (-(-((m : ℤ) + 1))) * u) := by
    refine (summable_nat_tateXterm_succ (u * w) u hq).congr fun m => ?_
    rw [hneg m]
  have key := tsum_of_add_one_of_neg_add_one
    (f := fun n : ℤ => tateXterm ((u * w) ^ (-n) * u)) h1 h2
  rw [key]
  simp only [hpos, hneg]
  rw [(summable_nat_tateXterm (u * w) w hq).tsum_eq_zero_add]
  simp only [pow_zero, one_mul, neg_zero, zpow_zero]
  ring

/-- ★★★★★★**`Y` の ℤ 和を頭と三本の尾に分ける**——形式側 `tateYpair` と同じ形。

★`g(1/t) = −(f(t)+g(t))` なので、負の段が `X` の尾と `Y` の尾の**両方**を生む。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tsum_int_tateYterm_split (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1)
    (hu0 : u ≠ 0) (hw0 : w ≠ 0) :
    (∑' n : ℤ, tateYterm ((u * w) ^ (-n) * u))
      = (tateYterm u + ∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * u))
        - (tateXterm w + ∑' m : ℕ, tateXterm ((u * w) ^ (m + 1) * w))
        - (tateYterm w + ∑' m : ℕ, tateYterm ((u * w) ^ (m + 1) * w)) := by
  have hq : ‖u * w‖ < 1 := norm_mul_lt_one hu hw
  have hq0 : u * w ≠ 0 := mul_ne_zero hu0 hw0
  have hpos : ∀ m : ℕ, tateYterm ((u * w) ^ (-((m : ℤ) + 1)) * u)
      = -(tateXterm ((u * w) ^ m * w) + tateYterm ((u * w) ^ m * w)) := by
    intro m
    have hne0 : (u * w) ^ m * w ≠ 0 := mul_ne_zero (pow_ne_zero _ hq0) hw0
    have hne1 : (u * w) ^ m * w ≠ 1 := ne_one_of_norm_lt_one (norm_pow_mul_lt_one hq hw m)
    rw [zpow_neg_succ_mul u w hu0 hw0 m, tateYterm_inv' hne0 hne1]
  have hneg : ∀ m : ℕ,
      tateYterm ((u * w) ^ (-(-((m : ℤ) + 1))) * u) = tateYterm ((u * w) ^ (m + 1) * u) := by
    intro m
    rw [neg_neg, show ((m : ℤ) + 1) = ((m + 1 : ℕ) : ℤ) by push_cast; ring, zpow_natCast]
  have hsX := summable_nat_tateXterm (u * w) w hq
  have hsY := summable_nat_tateYterm (u * w) w hq
  have h1 : Summable fun m : ℕ => tateYterm ((u * w) ^ (-((m : ℤ) + 1)) * u) := by
    refine ((hsX.add hsY).neg).congr fun m => ?_
    rw [hpos m]
  have h2 : Summable fun m : ℕ => tateYterm ((u * w) ^ (-(-((m : ℤ) + 1))) * u) := by
    refine (summable_nat_tateYterm_succ (u * w) u hq).congr fun m => ?_
    rw [hneg m]
  have key := tsum_of_add_one_of_neg_add_one
    (f := fun n : ℤ => tateYterm ((u * w) ^ (-n) * u)) h1 h2
  rw [key]
  simp only [hpos, hneg]
  rw [tsum_neg, hsX.tsum_add hsY, hsX.tsum_eq_zero_add, hsY.tsum_eq_zero_add]
  simp only [pow_zero, one_mul, neg_zero, zpow_zero]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def tsum_int_tateXterm_split.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X の ℤ 和を頭と尾に分ける)",
    sectionId := "genell-def-3-3" }

def tsum_int_tateYterm_split.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Y の ℤ 和を頭と尾に分ける)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
