import ABC3.Meta.Claim
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.Order.BigOperators.Ring.Finset

/-!
# [GenEll] Lemma 3.6 / Lemma 4.2 —— 2 理論を要しない初等的な評価(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.17・p.21。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.17):
> Lemma 3.6. (An Elementary Estimate) Let

原文 (GenEll p.21):
> Lemma 4.2. (Some Elementary Estimates) Let n be a positive integer;

## ★★なぜこの 2 件だけ先に実装できるのか

[GenEll] の必要 24 件のうち、**Arakelov 理論も l 進 Galois 表現も Tate 曲線も
要求しない**のは 4 件だけである(2026-08-16 実測):

| 項目 | 要るもの |
|---|---|
| `Lemma 3.1` | 群論。(i)(ii)(iii) は `Lemma31.lean` に実装済み、(iv) が [Serre] 待ち |
| **`Lemma 3.6`** | ★**実解析だけ**(本ファイル) |
| `Lemma 4.1` | 背理法と 6 行の不等式計算。条件 (i)(ii) は**仮説**である |
| **`Lemma 4.2`** | ★**初等的な不等式だけ**(本ファイル) |

★**`Lemma 4.1` の条件 (i)(ii) が仮説であること**は
`Remark 4.1.1` が「条件 (i) は素数定理の帰結」と書いて**明示的に圏外へ置いている**
ことから読める。素数定理は mathlib に無い(2026-08-16 実測)が、**それは効かない**。

## ★★原文の証明より短い道を取った箇所(`Lemma 4.2`)

原文 (GenEll p.21) は第 3 不等式を
「任意の正整数 `H` について `log(H+1) ≤ (3H/2)·log(2)`」から出す。

★**我々は `H + 1 ≤ 2^H`(`Nat.lt_two_pow_self`)を使った。**
これは `log(H+1) ≤ H·log 2` を与え、`(3H/2)·log 2` より**強い**。
★したがって結論 `≤ 3h/2` は**そのまま従い、しかも `≤ h` まで言える**。

★**主張は弱めていない**(G5)——原文の `3h/2` をそのまま証明している。
強い中間評価を経由しただけである。
-/

namespace ABC3.Found.GenEll

open Real

/-! ## Lemma 3.6 —— An Elementary Estimate -/

/-- **[GenEll] Lemma 3.6**(An Elementary Estimate)。

原文 (GenEll p.17):
> Lemma 3.6. (An Elementary Estimate) Let

`ϵ ∈ ℝ_{>0}` に対し**定数 `C₀ ∈ ℝ_{>0}` が存在して**、`y ≥ 1` かつ `x ≥ C₀·y^{1+ϵ}` なる
全ての `x, y ∈ ℝ` について **`x ≥ y·log(x)`**。

## ★定数を明示的に構成した

原文 (GenEll p.18):
> This follows immediately from the well-known elementary fact that x

——続きは `x^{1/(1+ϵ)}·log(x)/x = log(x)·x^{−ϵ/(1+ϵ)} → 0`(`x → ∞`)。
★原文は**極限**で片づけているが、**それでは定数が出ない**。

★我々は `δ ≝ ϵ/(1+ϵ)` と置き、**`C₀ ≝ (1/δ)^{1+ϵ}` を明示的に取った**。
道筋は次の 3 つの不等式である:

1. `log(x) ≤ x^δ/δ` —— `log(x^δ) ≤ x^δ − 1`(`Real.log_le_sub_one_of_pos`)から。
2. `y ≤ (x/C₀)^{1/(1+ϵ)}` —— 仮定を `1/(1+ϵ)` 乗して。
3. `(x/C₀)^{1/(1+ϵ)}·x^δ/δ = x` —— `C₀^{1/(1+ϵ)} = 1/δ` と
   `1/(1+ϵ) + δ = 1` から**きっかり `x` になる**。

★★**3 の等号がきっかり成り立つように `C₀` を選んである**。
`δ < 1` ゆえ `C₀ > 1`、したがって `x > 1` で `log x > 0` も従う——
この符号が無いと 2 の両辺に `log x` を掛けられない。 -/
theorem lemma_3_6 (eps : ℝ) (heps : 0 < eps) :
    ∃ C₀ : ℝ, 0 < C₀ ∧ ∀ x y : ℝ, 1 ≤ y → C₀ * y ^ (1 + eps) ≤ x → y * Real.log x ≤ x := by
  have h1e : (0 : ℝ) < 1 + eps := by linarith
  set δ : ℝ := eps / (1 + eps) with hδdef
  have hδ : 0 < δ := div_pos heps h1e
  have hδ1 : δ < 1 := by
    rw [hδdef, div_lt_one h1e]; linarith
  have hinvδ : 1 < 1 / δ := by
    rw [lt_div_iff₀ hδ]; linarith
  -- ★`C₀ ≝ (1/δ)^{1+ϵ}`
  set C₀ : ℝ := (1 / δ) ^ (1 + eps) with hC₀def
  have hC₀1 : 1 < C₀ := by
    have := Real.rpow_lt_rpow (by norm_num) hinvδ h1e
    rwa [Real.one_rpow] at this
  have hC₀0 : 0 < C₀ := lt_trans one_pos hC₀1
  have hmul : (1 + eps) * (1 / (1 + eps)) = 1 := by field_simp
  -- ★`C₀^{1/(1+ϵ)} = 1/δ`
  have hC₀root : C₀ ^ (1 / (1 + eps)) = 1 / δ := by
    rw [hC₀def, ← Real.rpow_mul (by positivity), hmul, Real.rpow_one]
  refine ⟨C₀, hC₀0, ?_⟩
  intro x y hy hx
  have hy0 : (0 : ℝ) < y := lt_of_lt_of_le one_pos hy
  -- `1 ≤ y^{1+ϵ}`
  have hyp1 : (1 : ℝ) ≤ y ^ (1 + eps) := by
    have := Real.rpow_le_rpow (by norm_num) hy h1e.le
    rwa [Real.one_rpow] at this
  -- `1 < C₀ ≤ C₀·y^{1+ϵ} ≤ x`
  have hx1 : 1 < x := by
    have : C₀ ≤ C₀ * y ^ (1 + eps) := le_mul_of_one_le_right hC₀0.le hyp1
    linarith
  have hx0 : (0 : ℝ) < x := lt_trans one_pos hx1
  have hlogpos : 0 ≤ Real.log x := Real.log_nonneg hx1.le
  -- ★① `log(x) ≤ x^δ/δ`
  have hlog : Real.log x ≤ x ^ δ / δ := by
    have hxδ : (0 : ℝ) < x ^ δ := Real.rpow_pos_of_pos hx0 δ
    have h := Real.log_le_sub_one_of_pos hxδ
    rw [Real.log_rpow hx0] at h
    rw [le_div_iff₀ hδ]
    nlinarith
  -- ★② `y ≤ (x/C₀)^{1/(1+ϵ)}`
  have hyle : y ≤ (x / C₀) ^ (1 / (1 + eps)) := by
    have hdiv : y ^ (1 + eps) ≤ x / C₀ := by
      rw [le_div_iff₀ hC₀0]; linarith [hx]
    have := Real.rpow_le_rpow (by positivity) hdiv (by positivity : (0:ℝ) ≤ 1 / (1 + eps))
    rwa [← Real.rpow_mul hy0.le, hmul, Real.rpow_one] at this
  -- ★③ `(x/C₀)^{1/(1+ϵ)} · x^δ/δ = x`
  have hkey : (x / C₀) ^ (1 / (1 + eps)) * (x ^ δ / δ) = x := by
    rw [Real.div_rpow hx0.le hC₀0.le, hC₀root]
    have hsum : 1 / (1 + eps) + δ = 1 := by
      rw [hδdef]; field_simp
    have : x ^ (1 / (1 + eps)) * x ^ δ = x := by
      rw [← Real.rpow_add hx0, hsum, Real.rpow_one]
    field_simp
    linarith [this]
  calc y * Real.log x
      ≤ (x / C₀) ^ (1 / (1 + eps)) * (x ^ δ / δ) := by
        exact mul_le_mul hyle hlog hlogpos (by positivity)
    _ = x := hkey

/-! ## Lemma 4.2 —— Some Elementary Estimates -/

/-- 原文が使う初等的事実の**強い版**。

原文 (GenEll p.21):
> follows from the easily verified fact that log(H +1) ≤ (3H/2)·log(2) for any positive integer H.

原文は「任意の正整数 `H` について `log(H+1) ≤ (3H/2)·log(2)`」を使う。
★**我々は `H + 1 ≤ 2^H`(`Nat.lt_two_pow_self`)から `log(H+1) ≤ H·log 2` を出す。**
係数 `3/2` が要らない分だけ強い。 -/
theorem log_succ_le_mul_log_two (H : ℕ) (hH : 0 < H) :
    Real.log ((H : ℝ) + 1) ≤ (H : ℝ) * Real.log 2 := by
  have hlt : H < 2 ^ H := Nat.lt_two_pow_self
  have hle : (H : ℝ) + 1 ≤ (2 : ℝ) ^ H := by
    have : (H : ℝ) + 1 ≤ ((2 ^ H : ℕ) : ℝ) := by exact_mod_cast hlt
    simpa using this
  have h0 : (0 : ℝ) < (H : ℝ) + 1 := by positivity
  calc Real.log ((H : ℝ) + 1)
      ≤ Real.log ((2 : ℝ) ^ H) := Real.log_le_log h0 hle
    _ = (H : ℝ) * Real.log 2 := by rw [Real.log_pow]

/-- **[GenEll] Lemma 4.2**(Some Elementary Estimates)。

原文 (GenEll p.21):
> Lemma 4.2. (Some Elementary Estimates) Let n be a positive integer;

`n` 正整数、`p_1,…,p_n` 素数、`h_1,…,h_n` 正整数、`h ≝ Σ h_j·log(p_j)` とすると
**`Σ log(p_j) ≤ h`** かつ **`Σ log(h_j) ≤ Σ log(h_j+1) ≤ 3h/2`**。

★3 つの不等式の道筋:
1. `h_j ≥ 1` と `log(p_j) ≥ 0` から各項ごとに `log(p_j) ≤ h_j·log(p_j)`。
2. `log` の単調性(`h_j ≤ h_j + 1`)。
3. ★`log(h_j+1) ≤ h_j·log 2 ≤ h_j·log(p_j)`(`p_j ≥ 2`)。総和して `≤ h ≤ 3h/2`。

★**3 で得ているのは原文より強い `≤ h` である**が、
主張は原文どおり `≤ 3h/2` のまま書いてある(G5——主張を強めも弱めもしない)。 -/
theorem lemma_4_2 (n : ℕ) (_hn : 0 < n) (p : Fin n → ℕ) (hp : ∀ j, (p j).Prime)
    (hh : Fin n → ℕ) (hhpos : ∀ j, 0 < hh j) :
    (∑ j, Real.log (p j)) ≤ (∑ j, (hh j : ℝ) * Real.log (p j))
  ∧ (∑ j, Real.log (hh j)) ≤ (∑ j, Real.log ((hh j : ℝ) + 1))
  ∧ (∑ j, Real.log ((hh j : ℝ) + 1))
      ≤ 3 / 2 * (∑ j, (hh j : ℝ) * Real.log (p j)) := by
  -- 各素数について `log 2 ≤ log(p j)`(とくに `0 ≤ log(p j)`)
  have hlog2 : ∀ j, Real.log 2 ≤ Real.log (p j) := fun j =>
    Real.log_le_log (by norm_num) (by exact_mod_cast (hp j).two_le)
  have hlogpos : ∀ j, 0 ≤ Real.log (p j) := fun j =>
    le_trans (Real.log_nonneg (by norm_num)) (hlog2 j)
  have hh1 : ∀ j, (1 : ℝ) ≤ (hh j : ℝ) := fun j => by exact_mod_cast hhpos j
  -- ★第 1 不等式
  have first : (∑ j, Real.log (p j)) ≤ (∑ j, (hh j : ℝ) * Real.log (p j)) :=
    Finset.sum_le_sum fun j _ => le_mul_of_one_le_left (hlogpos j) (hh1 j)
  -- ★第 3 不等式の主部: `Σ log(h_j+1) ≤ Σ h_j·log(p_j)`
  have third' : (∑ j, Real.log ((hh j : ℝ) + 1)) ≤ (∑ j, (hh j : ℝ) * Real.log (p j)) :=
    Finset.sum_le_sum fun j _ =>
      le_trans (log_succ_le_mul_log_two (hh j) (hhpos j))
        (mul_le_mul_of_nonneg_left (hlog2 j) (by positivity))
  have hsum0 : 0 ≤ (∑ j, (hh j : ℝ) * Real.log (p j)) :=
    Finset.sum_nonneg fun j _ => mul_nonneg (by positivity) (hlogpos j)
  refine ⟨first, ?_, ?_⟩
  · -- ★第 2 不等式
    exact Finset.sum_le_sum fun j _ =>
      Real.log_le_log (by exact_mod_cast hhpos j) (by linarith)
  · -- ★第 3 不等式(原文どおり `3h/2` で書く)
    linarith

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_6.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17, item := "Lemma 3.6",
    sectionId := "genell-lemma-3-6" }

def lemma_4_2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Lemma 4.2",
    sectionId := "genell-lemma-4-2" }

def log_succ_le_mul_log_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Lemma 4.2",
    sectionId := "genell-lemma-4-2" }

end ABC3.Found.GenEll
