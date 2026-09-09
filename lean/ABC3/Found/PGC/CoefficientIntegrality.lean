import ABC3.Found.PGC.ComponentBasisExpansion

/-!
# [pGC] 係数の整数性 —— `𝒪_M = 𝒪_K[π]` を経由しないで出る

## 持ち場（前波で私が挙げた最後の 1 点）

前波の `ComponentBasisExpansion.exists_expansion_of_valK` は展開 `x = Σ_{j<n} f_j π^j` を
出すが、★**係数が整である（`‖f_j‖ ≤ 1`）とは言っていなかった**。下流の
`ComponentExpansion` / `LossExponentMatch` はそれを仮説で受けている。本ファイルが供給する。

## ★開いて分かったこと —— **構造論は要らなかった**

私は前波で「一般に `𝒪_M ≠ 𝒪_K[π]`。全分岐なら等しいが、その段は未証明」と書き、
`PowerBasis` / `integralClosure` / `IsIntegrallyClosed` あたりの重い段を予期していた。
★**要らなかった。** 全分岐のときは**ノルムだけで**出る:

1. `‖f_j π^j‖ = ‖π‖^{n·m_j + j}` で、指数は `j` を法として ★**全部違う剰余類**にある。
2. ⇒ `‖Σ_j f_j π^j‖ = max_j ‖f_j π^j‖`（超距離、`nnnorm_sum_eq_sup_of_pairwise_ne`）。
3. ⇒ `‖x‖ ≤ 1` から各項が `≤ 1`、すなわち `n·m_j + j ≥ 0`。
4. ⇒ ★**`ℤ` の算術 1 本**（`0 ≤ j < n`）で `m_j ≥ 0`、すなわち `‖f_j‖ ≤ 1`。

★段 (4) が「係数が整である」ことの正体である（§2 `exponent_nonneg_of_lt`）。

## ★在庫の測定（索引がまた「無い」と嘘をついた —— 8 例目）

`IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne` は ★**mathlib に在る**が
`grep -n "nnnorm_sum_eq_sup_of_pairwise_ne" .cache/mathlib-index.txt` は **0 件**を返す
（`to_additive` 生成名は索引に載らない）。★見つけた経路は索引ではなく、
**木の `TotallyRamifiedValueGroup.lean:113` が実際に使っているのを読んだ**こと。

★同ファイルには `exists_zpow_norm`（「底が全分岐なら `M` の値群は `‖π‖^ℤ` に入る」）も在る。
本ファイルでは使わなかった（`hvalK` から直接出るため）。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `norm_term_le_sum_of_pairwise_ne` | ★抽象核。ノルムが相異なる項は和より大きくならない |
| §2 | `exponent_nonneg_of_lt` | ★抽象核。`ℤ` の算術 1 本。★整数性の正体 |
| §3 | `norm_coeff_le_one_of_valK` | ★係数の整数性 |
| §3 | `exists_integral_expansion_of_valK` | ★★**到達点**。`f n = 0` ∧ `∀ j, ‖f j‖ ≤ 1` ∧ `x = Σ f_j π^j` が**全部定理から出る** |
| §4 | `norm_sum_sub_deriv_le_scaled_tail` ほか | ★`j = 0` を仮定から外した版（下記の訂正） |

## ★★訂正 —— 前波の `LossExponentMatch.lean` §6 は `j = 0` に過剰要求していた

あちらの `norm_sum_sub_deriv_le_scaled` は `hf : ∀ j, ‖f j‖ ≤ C` と**全部の `j`** に要求する。
★これでは `C = ‖f_{j₀}‖` を代入できない —— `f_0`（`E₁` 成分）は `f_{j₀}` より
**大きいのが普通**だから。★`j = 0` の項は誤差に寄与しない（`(π+ρ)^0 − π^0 = 0`）ので、
仮定は `1 ≤ j` に制限してよい。§4 がその形で、
`MaxMinIndex.exists_no_cancel_index` が出す「`∀ j ∈ [1,p−1], ‖f j‖ ≤ ‖f j₀‖」がそのまま入る。

## ★測定（`tools/numerology-check.py` (n11)、5 設定）

`x` が整な標本で `v(f_j) ≥ 0` を数えた: ★**1590/1590 件で成立**、最小の `v(f_j)` は
ちょうど **0**（＝ 評価は sharp、単数係数が実際に現れる）。

| p | n | 標本 | `v(f_j) ≥ 0` |
|---|---|---|---|
| 3 | 3 | 200 | 600/600 |
| 3 | 4 | 60 | 180/180 |
| 5 | 3 | 40 | 200/200 |
| 7 | 2 | 30 | 210/210 |
| 2 | 4 | 200 | 400/400 |

## 逸脱の記録

- ★環としての `𝒪_M = 𝒪_K[π]` は**証明していない**。証明したのは**ノルムの主張**
  （整な元の係数は整）で、下流が要るのはこちらである。
- `hvalK`（`Γ_K ⊆ ‖π‖^{nℤ}`）は仮説のまま。★`TotallyRamifiedValueGroup.lean` の
  docstring が「`hvalK` 自身は落ちない（落としてはいけない）」と言っているのと整合する
  ——`hvalK` は「全分岐かつ `π` が素元」そのものだからである。
-/

namespace ABC3.Found.PGC

namespace CoefficientIntegrality

open Finset

/-! ## §1 抽象核 —— ノルムが相異なる項は和より大きくならない -/

section Kernel

/-- ★抽象核（分岐も付値も出ない）。超距離で**ノルムが相異なる**項の族なら、
各項のノルムは和のノルム以下。★`IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne`
（★mathlib に在るが `.cache/mathlib-index.txt` には**出てこない** —— `to_additive`
生成名は索引に載らない。木の `TotallyRamifiedValueGroup.lean:113` が使っているのを見て確認した）。 -/
theorem norm_term_le_sum_of_pairwise_ne {M ι : Type*} [NormedAddCommGroup M]
    [IsUltrametricDist M] {t : Finset ι} {y : ι → M}
    (hpair : (t : Set ι).Pairwise fun i j => ‖y i‖₊ ≠ ‖y j‖₊) {i : ι} (hi : i ∈ t) :
    ‖y i‖ ≤ ‖∑ j ∈ t, y j‖ := by
  have hsup : ‖∑ j ∈ t, y j‖₊ = t.sup fun j => ‖y j‖₊ :=
    IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne hpair
  have hle : ‖y i‖₊ ≤ ‖∑ j ∈ t, y j‖₊ := by
    rw [hsup]
    exact Finset.le_sup (f := fun j => ‖y j‖₊) hi
  exact_mod_cast hle

end Kernel

/-! ## §2 抽象核 2 —— 指数の非負性 -/

section Exponent

/-- ★抽象核（`ℤ` の算術 1 本）。`n·m + j ≥ 0` と `0 ≤ j < n` なら `m ≥ 0`。
★ここが「係数が整である」ことの正体である。 -/
theorem exponent_nonneg_of_lt {n j : ℕ} {m : ℤ} (hj : j < n)
    (h : 0 ≤ (n : ℤ) * m + (j : ℤ)) : 0 ≤ m := by
  by_contra hcon
  have hm : m ≤ -1 := by omega
  have hjn : (j : ℤ) < (n : ℤ) := by exact_mod_cast hj
  have hn0 : (0 : ℤ) < (n : ℤ) := by omega
  have : (n : ℤ) * m ≤ (n : ℤ) * (-1) := by
    exact mul_le_mul_of_nonneg_left hm (le_of_lt hn0)
  omega

end Exponent

/-! ## §3 ★係数の整数性 —— `‖x‖ ≤ 1` なら `‖f_j‖ ≤ 1` -/

section Integral

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★**係数の整数性**。底が全分岐（`Γ_K ⊆ ‖π‖^{nℤ}`）なら、整な元の展開の係数は整。

★構造論（`𝒪_M = 𝒪_K[π]`）を経由しない。★スロットの剰余が全部違うので
`‖Σ f_j π^j‖ = max_j ‖f_j π^j‖` となり、そこから **`ℤ` の算術 1 本**で出る。 -/
theorem norm_coeff_le_one_of_valK {π : M} {n : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (c : ℕ → K) (hx : ‖∑ i ∈ range n, algebraMap K M (c i) * π ^ i‖ ≤ 1)
    {j : ℕ} (hj : j < n) :
    ‖algebraMap K M (c j)‖ ≤ 1 := by
  classical
  by_cases hc : c j = 0
  · simp [hc]
  set y : ℕ → M := fun i => algebraMap K M (c i) * π ^ i with hy
  set t : Finset ℕ := (range n).filter (fun i => c i ≠ 0) with ht
  have hmem : ∀ i, i ∈ t ↔ (i ∈ range n ∧ c i ≠ 0) := fun i => Finset.mem_filter
  have hexp : ∀ i ∈ t, ∃ m : ℤ, ‖algebraMap K M (c i)‖ = ‖π‖ ^ ((n : ℤ) * m) ∧
      ‖y i‖ = ‖π‖ ^ ((n : ℤ) * m + (i : ℤ)) := by
    intro i hi
    obtain ⟨m, hm⟩ := hvalK (c i) ((hmem i).mp hi).2
    refine ⟨m, hm, ?_⟩
    rw [hy]
    simp only [norm_mul, norm_pow, hm]
    rw [zpow_add₀ (ne_of_gt hπ0), zpow_natCast]
  have hpair : (t : Set ℕ).Pairwise fun i i' => ‖y i‖₊ ≠ ‖y i'‖₊ := by
    intro i hi i' hi' hii hcon
    have hreal : ‖y i‖ = ‖y i'‖ := by simpa using congrArg NNReal.toReal hcon
    obtain ⟨m, -, hm⟩ := hexp i hi
    obtain ⟨m', -, hm'⟩ := hexp i' hi'
    rw [hm, hm'] at hreal
    have hzeq := (zpow_right_inj₀ hπ0 (ne_of_lt hπ1)).mp hreal
    have hin : i < n := Finset.mem_range.mp ((hmem i).mp hi).1
    have hin' : i' < n := Finset.mem_range.mp ((hmem i').mp hi').1
    refine TotallyRamifiedLayer.not_dvd_sub_of_lt hin hin' hii ⟨m' - m, ?_⟩
    linear_combination hzeq
  have hjt : j ∈ t := (hmem j).mpr ⟨Finset.mem_range.mpr hj, hc⟩
  have hsum : ∑ i ∈ range n, y i = ∑ i ∈ t, y i := by
    refine (Finset.sum_subset (Finset.filter_subset _ _) ?_).symm
    intro i hi hit
    have hci : c i = 0 := by
      by_contra hne
      exact hit ((hmem i).mpr ⟨hi, hne⟩)
    simp [hy, hci]
  have hle : ‖y j‖ ≤ 1 := by
    refine le_trans (norm_term_le_sum_of_pairwise_ne hpair hjt) ?_
    rw [← hsum]
    exact hx
  obtain ⟨m, hmc, hm⟩ := hexp j hjt
  rw [hm] at hle
  have hnn : 0 ≤ (n : ℤ) * m + (j : ℤ) := (zpow_le_one_iff_right_of_lt_one₀ hπ0 hπ1).mp hle
  have hm0 : 0 ≤ m := exponent_nonneg_of_lt hj hnn
  rw [hmc]
  refine (zpow_le_one_iff_right_of_lt_one₀ hπ0 hπ1).mpr ?_
  positivity

/-- ★★★**本ファイルの到達点** —— 整な元は**整な係数で**展開できる。

★これで `ComponentExpansion` / `LossExponentMatch` が要求する 3 つ
（`f n = 0`、`∀ j, ‖f j‖ ≤ 1`、`x = Σ f_j π^j`）が**すべて定理から出る**。 -/
theorem exists_integral_expansion_of_valK [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {x : M} (hx : ‖x‖ ≤ 1) :
    ∃ f : ℕ → M, (∀ j, n ≤ j → f j = 0) ∧ (∀ j, ‖f j‖ ≤ 1) ∧
      x = ∑ j ∈ range n, f j * π ^ j := by
  obtain ⟨f, hf0, hfa, hfsum⟩ :=
    ComponentBasisExpansion.exists_expansion_of_valK hπ0 (ne_of_lt hπ1) hn hvalK x
  choose c hc using hfa
  have hsum' : ∑ i ∈ range n, algebraMap K M (c i) * π ^ i = x := by
    rw [hfsum]
    exact Finset.sum_congr rfl fun i _ => by rw [hc i]
  refine ⟨f, hf0, ?_, hfsum⟩
  intro j
  by_cases hj : j < n
  · rw [hc j]
    refine norm_coeff_le_one_of_valK hπ0 hπ1 hvalK c ?_ hj
    rw [hsum']
    exact hx
  · rw [hf0 j (by omega)]
    simp

end Integral

/-! ## §4 ★`j = 0` を仮定から外す（★これが無いと `C = ‖f_{j₀}‖` を代入できない） -/

section DropZero

/-- ★★訂正（前波の私のファイル `LossExponentMatch.lean` §6 について）。

あちらの `norm_sum_sub_deriv_le_scaled` は `hf : ∀ j, ‖f j‖ ≤ C` と**全部の `j`** に
要求している。★これでは `C = ‖f_{j₀}‖` を代入できない ——
`f_0`（`E₁` 成分）は `f_{j₀}` より**大きいのが普通**だからである。

★実は `j = 0` の項は誤差に**寄与しない**（`(π+ρ)^0 − π^0 = 0`）。証明の
`match j with | 0 => …` の枝がそれを言っている。⇒ 仮定を `1 ≤ j` に制限してよい。
本定理はその形。★`MaxMinIndex.exists_no_cancel_index` が出す
「`∀ j ∈ [1,p−1], ‖f j‖ ≤ ‖f j₀‖`」がそのまま入る。 -/
theorem norm_sum_sub_deriv_le_scaled_tail {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (f : ℕ → M) {C : ℝ} (hC : 0 ≤ C)
    (hf : ∀ j, 1 ≤ j → ‖f j‖ ≤ C) (n : ℕ) :
    ‖(∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        - ρ * ∑ j ∈ range n, (j : M) * f j * π ^ (j - 1)‖ ≤ C * ‖ρ‖ ^ 2 := by
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
  intro j _
  match j with
  | 0 => simpa using by positivity
  | (m + 1) =>
      have hrw : f (m + 1) * ((π + ρ) ^ (m + 1) - π ^ (m + 1))
            - ρ * ((((m + 1 : ℕ)) : M) * f (m + 1) * π ^ (m + 1 - 1))
          = f (m + 1) * ((π + ρ) ^ (m + 1) - π ^ (m + 1) - ((m : M) + 1) * π ^ m * ρ) := by
        simp only [Nat.add_sub_cancel]
        push_cast
        ring
      rw [hrw, norm_mul]
      exact mul_le_mul (hf _ (by omega))
        (BinomialFirstOrder.norm_add_pow_sub_linear_le hπ hρ m) (norm_nonneg _) hC

/-- ★上と `ComponentExpansion.sum_deriv_shift` を合わせた形（`j = 0` を外した版）。 -/
theorem norm_expansion_main_le_scaled_tail {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ w : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (hw : ρ = w * (1 + π))
    (f : ℕ → M) {C : ℝ} (hC : 0 ≤ C) (hf : ∀ j, 1 ≤ j → ‖f j‖ ≤ C) {n : ℕ} (hfn : f n = 0) :
    ‖(∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j‖
      ≤ C * ‖ρ‖ ^ 2 := by
  have key : w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
      = ρ * ∑ j ∈ range n, (j : M) * f j * π ^ (j - 1) := by
    rw [ComponentExpansion.sum_deriv_shift, hfn, hw]
    ring
  rw [key]
  exact norm_sum_sub_deriv_le_scaled_tail hπ hρ f hC hf n

end DropZero

/-! ## §5 使っている公理の一覧 -/

#print axioms norm_term_le_sum_of_pairwise_ne
#print axioms exponent_nonneg_of_lt
#print axioms norm_coeff_le_one_of_valK
#print axioms exists_integral_expansion_of_valK
#print axioms norm_sum_sub_deriv_le_scaled_tail
#print axioms norm_expansion_main_le_scaled_tail

end CoefficientIntegrality

end ABC3.Found.PGC
