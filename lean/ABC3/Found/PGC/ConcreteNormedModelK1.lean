import ABC3.Found.PGC.ConcreteNormedModelP3

/-!
# [pGC] ★★★★★ `k = 1`（`p²` 次）の具体模型 —— `ℚ₂(i)((1+i)^{1/4}) / ℚ₂(i)`

`JumpFromValueGroup.exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree` を
★**`k = 1`（4 次巡回・全分岐）の具体的な体に、仮説を 1 つも残さずに代入した**もの
（`model_k1_exists`）。

## ★★これまでとの違い —— **Hasse–Arf の合同が初めて具体例で試された**

`Found/PGC/ConcreteNormedModel.lean`（`p = 2`, `k = 0`）と
`Found/PGC/ConcreteNormedModelP3.lean`（`p = 3`, `k = 0`）は
どちらも `k = 0` なので、`harith` の中段 2 条件

* `∀ m, m < k → u m < u (m+1)`（増加）
* `∀ m, m < k → p^{m+1} ∣ u (m+1) − u m`（★**Hasse–Arf の合同**）

が `m < 0` で**空虚**であった。本ファイルは `k = 1` なので `m = 0` が実際に現れ、
`2 ∣ u₁ − u₀ = 8 − 4 = 4` を**証明して**通している（`harith_k1`）。

## 数値の勘定（★2 通りの独立な方法で検算した）

`K = ℚ₂(i)`, `ϖ = 1 + i`（`K` の素元、`ϖ² = 2i`）, `M = K(π)`, `π⁴ = ϖ`,
`p = 2`, `k = 1`, `e = v_K(2) = 2`, `v_M(π) = 1`, `v_M(2) = 8`。

| 仮説 | 値 | 本ファイル |
|---|---|---|
| `hnK` | `[M:K] = 4 = p^{k+1}` | `finrank_M4` |
| `hnormp` | `‖2‖ = 1/2` | `norm_two_M4` |
| `heM` | `‖2‖ = ‖π‖^{p^{k+1}·e} = ‖π‖⁸` | `heM4` |
| `hvalK` | `‖a‖ ∈ ‖π‖^{4ℤ}`（`a ∈ K`） | `valK4` |
| `htop` | `K[π] = M` | `adjoin_root_top` |
| `hg` | `orderOf g = 4`（`g : π ↦ i·π`） | `orderOf_g4` |
| 第 0 跳び | `gπ − π = (i−1)π`, `‖i−1‖ = ‖ϖ‖ = ‖π‖⁴` ⇒ `‖π‖⁵` ⇒ `u₀ = 4` | `g4_pi4_sub` |
| 第 1 跳び | `g²π − π = −2π` ⇒ `‖π‖⁹` ⇒ `u₁ = 8` | `g4_sq_pi4_sub` |
| ★`harith` | `1 ≤ 4`, `4 < 8`, `2 ∣ 4`, `1·8 ≤ 4·2 = 8`（等号） | `harith_k1` |

★**独立な検算**（形式化はしていないが手計算で 1 度）: 相対差積は
`Σ_{s≠1} i_G(s) = 5 + 9 + 5 = 19` で、`f = X⁴ − ϖ` の `v_M(f'(π)) = v_M(4π³) = 16 + 3 = 19`
と一致する。上付き番号は `φ(4) = 4`, `φ(8) = 4 + (8−4)/2 = 6` でどちらも整数
（＝ Hasse–Arf の定理そのもの）。

## ★★配られた見立ての検算（3 点のうち 2 点が当たり、1 点は理由が違った）

1. 「`p = 2`（4 次）が一番安いはず」→ ★**当たり**。ただし理由は違う。
   `p = 2` は Kummer 塔が **2 段**（`ℚ₂ → ℚ₂(i) → M`）で済むのに対し、
   `p = 3` の `k = 1` は `μ₉ ⊂ K` が要るので `K = ℚ₃(ζ₉)`（`ℚ₃` 上 6 次）となり
   **3 段**必要になる。
2. 「★`ζ_{p²}` が要る可能性がある」→ ★★**当たり**。`k = 1` の Kummer には
   `μ_{p²} ⊂ K` が要り、`p = 2` では `ζ₄ = i`。だから底が `ℚ₂` ではなく `ℚ₂(i)` になる。
3. 「★Hasse–Arf の合同で偽が出る可能性がある」→ ★**外れ（真だった）**。
   `2 ∣ 8 − 4` は成り立つ。★偽は出ていない。

## ★★★配られた字面からの逸脱（記録）

* 持ち場は「`k = 1`（`p²` 次）の模型を 1 つ」であり、`p` は指定されていない。
  ★`p = 2`, `e = 2`, `K = ℚ₂(i)` を選んだ。出口定理は `[Fact p.Prime]` だけを要求し
  `p ≠ 2` を要求していないので、原典の設定の正当な一例である。
* ★`k = 1` の Kummer 拡大に `μ₄ ⊂ K` が要るため、**底が `ℚ_p` そのものではない**。
  これは `k = 0` の 2 本（`ℚ₂` / `ℚ₃(ζ₃)`）と違う点で、`e = 2` の由来でもある。

## ★★在庫の測定（コマンドを残す）

```
grep -n "X_pow_sub_C" .cache/mathlib-index.txt
  → X_pow_sub_C_irreducible_of_prime_pow は ★**(hp' : p ≠ 2) を要求する**。
    mathlib 側にも "-- TODO: generalize to `p = 2`" が書かれている
    (FieldTheory/KummerExtension.lean:145)。
    ⇒ ★`p = 2`, `n = 4` は mathlib の判定法では**扱えない**。
      本ファイルの `X_pow_four_sub_C_irreducible`（抽象核 1）がその穴を
      十分条件の形で埋める（`a ∉ K²` かつ `−a ∉ K²` ⇒ `X⁴ − a` 既約）。
      ★これは古典的判定「`a ∉ K²` かつ `a ∉ −4K⁴`」より強い条件だが、
      `a = −4b⁴ ⇒ −a = (2b²)²` なので**整合している**。
grep -n "toZModPow" .cache/mathlib-index.txt
  → PadicInt.toZModPow (n) : ℤ_[p] →+* ZMod (p ^ n)。
    ★`−1` が `ℚ₂` の平方でないことを `ZMod 4` に落として `decide` で閉じた
    （`neg_one_not_sq`）。`exact?` では出ない形。
grep -n "orderOf_eq_prime_pow" .cache/mathlib-index.txt
  → GroupTheory/OrderOfElement.lean:529 に在る。★位数 4 を「4 の約数」から
    絞り込む手作業が要らなくなった。
```

★**在るが使えなかったもの**: `ConcreteNormedModel.exists_zpow_norm_of_quadratic`
（`k = 0` の 2 本が使った値群の抽象核）は `0 < ‖π‖`, `‖π‖ ≠ 1` を要求するので、
★**生成元が単数の場合（`ℚ₂(i)` の `i`、`‖i‖ = 1`）には使えない**。
そこで抽象核 2（`exists_zpow_norm_of_sq_rel`、2 次の関係式 `z² = bz − c` だけを使う）
を新しく切り出した。

## ★抽象核（分岐・付値・Galois の語彙が 1 語も出ない）

1. `X_pow_four_sub_C_irreducible` —— 一般の体 `K`。
2. `exists_zpow_norm_of_sq_rel` —— 一般のノルム体 `L`（超距離）。
3. `not_pow_of_norm_zpow` —— 一般のノルム体 `L`。

具体層（`K1` / `M4` / `g4` / `harith_k1` / `model_k1_exists`）はこの 3 本に代入するだけ。
-/

open Polynomial

namespace ABC3.Found.PGC
namespace ConcreteNormedModelK1

section AbstractIrred

variable {K : Type*} [Field K]

/-- ★★**抽象核 1** —— 分岐・付値・Galois の語彙が 1 語も出ない。
`a` も `−a` も `K` の平方でなければ `X⁴ − a` は既約。

★mathlib は `X_pow_sub_C_irreducible_of_prime_pow` に
`-- TODO: generalize to `p = 2`` を残しており、`p = 2` の冪は**扱えない**。
本補題はその `n = 4` の場合を（十分条件の形で）埋める。 -/
theorem X_pow_four_sub_C_irreducible {a : K}
    (h1 : ∀ b : K, b ^ 2 ≠ a) (h2 : ∀ b : K, b ^ 2 ≠ -a) :
    Irreducible (X ^ 4 - C a) := by
  have h4 : (4 : ℕ) = 2 * 2 := by norm_num
  rw [h4]
  refine X_pow_mul_sub_C_irreducible
    (X_pow_sub_C_irreducible_of_prime Nat.prime_two h1) ?_
  intro E _ _ x hx
  refine X_pow_sub_C_irreducible_of_prime Nat.prime_two ?_
  intro b hb
  have hint : IsIntegral K x := by
    refine not_not.mp fun h => ?_
    simpa only [degree_zero, degree_X_pow_sub_C (by norm_num : 0 < 2),
      WithBot.natCast_ne_bot] using congr_arg degree (hx.symm.trans (dif_neg h))
  refine h2 (Algebra.norm K b) ?_
  rw [← map_pow, hb, ← IntermediateField.adjoin.powerBasis_gen hint,
    Algebra.PowerBasis.norm_gen_eq_coeff_zero_minpoly]
  simp [IntermediateField.adjoin.powerBasis_gen, IntermediateField.minpoly_gen, hx]

end AbstractIrred

section AbstractValue

/-- 実数の平方根の一意性（非負の範囲）。 -/
theorem eq_of_sq_eq_sq_of_nonneg {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) (h : x ^ 2 = y ^ 2) :
    x = y := by
  rcases lt_trichotomy x y with hlt | heq | hgt
  · nlinarith
  · exact heq
  · nlinarith

/-- ★★**抽象核 2** —— 分岐・付値・Galois の語彙が 1 語も出ない。

ノルム体 `L` の元 `z` が **2 次の関係式** `z² = b·z − c` をみたし、
`b`, `c` のノルムがどちらも `‖ϖ‖^{2ℤ}` に入るなら、`‖z‖` は `‖ϖ‖^ℤ` に入る。

★これは `ConcreteNormedModel.exists_zpow_norm_of_quadratic`（分解 `x + yπ` を使う）が
**使えない場合**（生成元 `π` が単数、たとえば `ℚ₂(i)` の `i`）を覆う。 -/
theorem exists_zpow_norm_of_sq_rel {L : Type*} [NormedField L] [IsUltrametricDist L]
    {ϖ : L} (h0 : 0 < ‖ϖ‖) {z b c : L} (hz : z ≠ 0)
    (hrel : z ^ 2 = b * z - c)
    (hb : b ≠ 0 → ∃ m : ℤ, ‖b‖ = ‖ϖ‖ ^ (2 * m))
    (hc : c ≠ 0 → ∃ m : ℤ, ‖c‖ = ‖ϖ‖ ^ (2 * m)) :
    ∃ m : ℤ, ‖z‖ = ‖ϖ‖ ^ m := by
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz
  by_cases hc0 : c = 0
  · subst hc0
    have : z * (z - b) = 0 := by rw [sub_zero] at hrel; ring_nf; linear_combination hrel
    have hzb : z = b := by
      rcases mul_eq_zero.mp this with h | h
      · exact absurd h hz
      · exact sub_eq_zero.mp h
    obtain ⟨m, hm⟩ := hb (hzb ▸ hz)
    exact ⟨2 * m, by rw [hzb, hm]⟩
  obtain ⟨mc, hmc⟩ := hc hc0
  by_cases hb0 : b = 0
  · subst hb0
    have h2 : ‖z‖ ^ 2 = ‖ϖ‖ ^ mc * ‖ϖ‖ ^ mc := by
      rw [← norm_pow, hrel, zero_mul, zero_sub, norm_neg, hmc, two_mul, zpow_add₀ (ne_of_gt h0)]
    refine ⟨mc, eq_of_sq_eq_sq_of_nonneg (norm_nonneg _) (le_of_lt (zpow_pos h0 mc)) ?_⟩
    rw [h2]; ring
  obtain ⟨mb, hmb⟩ := hb hb0
  by_cases heq : ‖b * z‖ = ‖c‖
  · rw [norm_mul, hmb, hmc] at heq
    refine ⟨2 * mc - 2 * mb, ?_⟩
    rw [zpow_sub₀ (ne_of_gt h0), eq_div_iff (ne_of_gt (zpow_pos h0 (2 * mb))), mul_comm]
    exact heq
  · have hne : ‖b * z‖ ≠ ‖(-c : L)‖ := by rwa [norm_neg]
    have hmax : ‖z‖ ^ 2 = max ‖b * z‖ ‖(-c : L)‖ := by
      rw [← norm_pow, hrel, sub_eq_add_neg]
      exact IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm hne
    rw [norm_neg, norm_mul] at hmax
    rcases le_total (‖b‖ * ‖z‖) ‖c‖ with hle | hle
    · rw [max_eq_right hle] at hmax
      refine ⟨mc, eq_of_sq_eq_sq_of_nonneg (norm_nonneg _) (le_of_lt (zpow_pos h0 mc)) ?_⟩
      rw [hmax, hmc, two_mul, zpow_add₀ (ne_of_gt h0)]; ring
    · rw [max_eq_left hle] at hmax
      refine ⟨2 * mb, ?_⟩
      rw [← hmb]
      have : ‖z‖ * ‖z‖ = ‖b‖ * ‖z‖ := by rw [← sq]; exact hmax
      exact (mul_right_cancel₀ (ne_of_gt hzn) this)

end AbstractValue

section QTwoBase

/-- `−1` は `ℚ₂` の平方でない（`ZMod 4` に落として `decide`）。
★`ℚ₂(i)` を作るのに要る唯一の非自明な入力。 -/
theorem neg_one_not_sq (b : ℚ_[2]) : b ^ 2 ≠ (-1 : ℚ_[2]) := by
  intro hb2
  have hnb : ‖b‖ ≤ 1 := by
    have : ‖b‖ ^ 2 = 1 := by rw [← norm_pow, hb2, norm_neg, norm_one]
    nlinarith [norm_nonneg b]
  set B : ℤ_[2] := ⟨b, hnb⟩ with hB
  have hcoe : ((B ^ 2 : ℤ_[2]) : ℚ_[2]) = ((-1 : ℤ_[2]) : ℚ_[2]) := by
    simp only [PadicInt.coe_pow, PadicInt.coe_neg, PadicInt.coe_one]
    exact hb2
  have h2 : B ^ 2 = (-1 : ℤ_[2]) := Subtype.ext hcoe
  have h3 : (PadicInt.toZModPow 2 B) ^ 2 = (-1 : ZMod (2 ^ 2)) := by
    rw [← map_pow, h2, map_neg, map_one]
  revert h3
  generalize (PadicInt.toZModPow 2 B) = t
  revert t
  decide

end QTwoBase

section BaseField

open ConcreteNormedModel

/-- `X² + 1` は `ℚ₂` 上既約。 -/
instance factIrrK1 : Fact (Irreducible (X ^ 2 - C (-1 : ℚ_[2]))) :=
  ⟨X_pow_sub_C_irreducible_of_prime Nat.prime_two neg_one_not_sq⟩

/-- ★底の体 `K1 = ℚ₂(i)`（`μ₄ ⊂ K1` にするために必要）。 -/
abbrev K1 : Type := AdjoinRoot (X ^ 2 - C (-1 : ℚ_[2]))

instance : FiniteDimensional ℚ_[2] K1 :=
  (AdjoinRoot.powerBasis (poly_ne_zero 2 (-1 : ℚ_[2]))).finite

attribute [local instance] Algebra.IsAlgebraic.of_finite

@[implicit_reducible] noncomputable instance normedFieldK1 : NormedField K1 :=
  spectralNorm.normedField ℚ_[2] K1

noncomputable instance ultraK1 : IsUltrametricDist K1 :=
  IsUltrametricDist.isUltrametricDist_of_forall_norm_add_le_max_norm
    (isNonarchimedean_spectralNorm (K := ℚ_[2]) (L := K1))

theorem norm_algebraMap_K1 (x : ℚ_[2]) : ‖algebraMap ℚ_[2] K1 x‖ = ‖x‖ := spectralNorm_extends x

/-- `i = √−1`。 -/
noncomputable def iota : K1 := AdjoinRoot.root (X ^ 2 - C (-1 : ℚ_[2]))

theorem iota_sq' : iota ^ 2 = algebraMap ℚ_[2] K1 (-1) := root_pow_eq (n := 2) (-1 : ℚ_[2])

theorem iota_sq : iota ^ 2 = -1 := by
  rw [iota_sq']
  simp

theorem norm_two_K1 : ‖(2 : K1)‖ = (2 : ℝ)⁻¹ := by
  rw [show (2 : K1) = ((2 : ℕ) : K1) by norm_num, ← map_natCast (algebraMap ℚ_[2] K1) 2,
    norm_algebraMap_K1]
  simpa using Padic.norm_p (p := 2)

theorem two_ne_zero_K1 : (2 : K1) ≠ 0 := by
  intro h
  have := norm_two_K1
  rw [h, norm_zero] at this
  norm_num at this

theorem norm_iota : ‖iota‖ = 1 := by
  have h : ‖iota‖ ^ 2 = 1 := by rw [← norm_pow, iota_sq, norm_neg, norm_one]
  nlinarith [norm_nonneg iota]

/-- ★素元 `ϖ = 1 + i`（`ϖ² = 2i`）。 -/
noncomputable def varpi : K1 := 1 + iota

theorem varpi_sq : varpi ^ 2 = 2 * iota := by
  unfold varpi
  linear_combination iota_sq

theorem norm_varpi_sq : ‖varpi‖ ^ 2 = 1 / 2 := by
  rw [← norm_pow, varpi_sq, norm_mul, norm_iota, norm_two_K1]
  norm_num

theorem norm_varpi_pos : 0 < ‖varpi‖ := by nlinarith [norm_varpi_sq, norm_nonneg varpi]

theorem norm_varpi_lt_one : ‖varpi‖ < 1 := by nlinarith [norm_varpi_sq, norm_varpi_pos]

theorem norm_varpi_ne_one : ‖varpi‖ ≠ 1 := ne_of_lt norm_varpi_lt_one

theorem varpi_ne_zero : varpi ≠ 0 := norm_pos_iff.mp norm_varpi_pos

end BaseField

section ValueGroupK1

open ConcreteNormedModel

/-- ★★**抽象核 3** —— 値群の指数が `n` で割れなければ `n` 乗ではない。
分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem not_pow_of_norm_zpow {L : Type*} [NormedField L] {ϖ : L} (h0 : 0 < ‖ϖ‖) (h1 : ‖ϖ‖ ≠ 1)
    (hval : ∀ z : L, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖ϖ‖ ^ m)
    {a : L} {r : ℤ} (ha : ‖a‖ = ‖ϖ‖ ^ r) {n : ℕ} (hn0 : n ≠ 0) (hn : ¬ (n : ℤ) ∣ r) (b : L) :
    b ^ n ≠ a := by
  intro hb
  have ha0 : a ≠ 0 := by
    intro h
    rw [h, norm_zero] at ha
    exact absurd ha.symm (ne_of_gt (zpow_pos h0 r))
  have hb0 : b ≠ 0 := by
    intro h
    rw [h, zero_pow hn0] at hb
    exact ha0 hb.symm
  obtain ⟨m, hm⟩ := hval b hb0
  have hmn : ‖ϖ‖ ^ (m * (n : ℤ)) = ‖ϖ‖ ^ r := by
    rw [zpow_mul, ← hm, zpow_natCast, ← norm_pow, hb, ha]
  exact hn ⟨m, by rw [← zpow_right_injective₀ h0 h1 hmn]; ring⟩

theorem base_val_K1 (x : ℚ_[2]) (hx : x ≠ 0) :
    ∃ m : ℤ, ‖algebraMap ℚ_[2] K1 x‖ = ‖varpi‖ ^ ((2 : ℤ) * m) := by
  refine ⟨Padic.valuation x, ?_⟩
  rw [norm_algebraMap_K1, Padic.norm_eq_zpow_neg_valuation hx, zpow_mul]
  have h2 : ‖varpi‖ ^ (2 : ℤ) = ((2 : ℕ) : ℝ) ^ (-1 : ℤ) := by
    rw [show (2 : ℤ) = ((2 : ℕ) : ℤ) by norm_num, zpow_natCast, norm_varpi_sq]
    norm_num
  rw [h2, ← zpow_mul]
  ring_nf

theorem decomp_K1 (z : K1) :
    ∃ x y : ℚ_[2], z = algebraMap ℚ_[2] K1 x + algebraMap ℚ_[2] K1 y * iota :=
  quadratic_decomp (-1 : ℚ_[2]) z

/-- ★`K1 = ℚ₂(i)` の値群は `‖ϖ‖^ℤ`（＝ `K1/ℚ₂` は全分岐）。
★★生成元 `i` は**単数**なので `exists_zpow_norm_of_quadratic`（分解を使う核）は
使えず、抽象核 2（2 次の関係式を使う核）で通す。 -/
theorem val_K1 {z : K1} (hz : z ≠ 0) : ∃ m : ℤ, ‖z‖ = ‖varpi‖ ^ m := by
  obtain ⟨x, y, rfl⟩ := decomp_K1 z
  refine exists_zpow_norm_of_sq_rel (ϖ := varpi) (b := algebraMap ℚ_[2] K1 (2 * x))
    (c := algebraMap ℚ_[2] K1 (x ^ 2 + y ^ 2)) norm_varpi_pos hz ?_ ?_ ?_
  · have h1 : algebraMap ℚ_[2] K1 (2 * x) = 2 * algebraMap ℚ_[2] K1 x := by
      rw [map_mul, map_ofNat]
    have h2 : algebraMap ℚ_[2] K1 (x ^ 2 + y ^ 2)
        = algebraMap ℚ_[2] K1 x ^ 2 + algebraMap ℚ_[2] K1 y ^ 2 := by
      rw [map_add, map_pow, map_pow]
    rw [h1, h2]
    linear_combination (algebraMap ℚ_[2] K1 y ^ 2) * iota_sq
  · intro hb0
    refine base_val_K1 (2 * x) (fun h => hb0 ?_)
    rw [h, map_zero]
  · intro hc0
    refine base_val_K1 (x ^ 2 + y ^ 2) (fun h => hc0 ?_)
    rw [h, map_zero]

/-- ★`ϖ` は `K1` の平方でない（値群の指数 1 が 2 で割れない）。 -/
theorem varpi_not_sq (b : K1) : b ^ 2 ≠ varpi :=
  not_pow_of_norm_zpow norm_varpi_pos norm_varpi_ne_one (fun _ hz => val_K1 hz)
    (a := varpi) (r := 1) (by rw [zpow_one]) two_ne_zero (by decide) b

/-- ★`−ϖ` も `K1` の平方でない。 -/
theorem neg_varpi_not_sq (b : K1) : b ^ 2 ≠ -varpi :=
  not_pow_of_norm_zpow norm_varpi_pos norm_varpi_ne_one (fun _ hz => val_K1 hz)
    (a := -varpi) (r := 1) (by rw [zpow_one, norm_neg]) two_ne_zero (by decide) b

end ValueGroupK1

section Zeta4

/-- ★`i` は `K1` の原始 4 乗根（`k = 1` の Kummer 拡大に必要）。 -/
theorem iota_prim : IsPrimitiveRoot iota 4 := by
  refine IsPrimitiveRoot.mk_of_lt iota (by norm_num) ?_ ?_
  · have h := iota_sq
    linear_combination (iota ^ 2 - 1) * h
  · intro l hl0 hl4
    interval_cases l
    · intro h
      exact two_ne_zero_K1 (by linear_combination iota_sq - (iota + 1) * h)
    · intro h
      exact two_ne_zero_K1 (by linear_combination iota_sq - h)
    · intro h
      exact two_ne_zero_K1 (by
        linear_combination (1 - iota ^ 2 + iota) * iota_sq + (iota - 1) * h)

end Zeta4

section TopField

open ConcreteNormedModel

attribute [local instance] Algebra.IsAlgebraic.of_finite

@[implicit_reducible] noncomputable instance normedAlgebraK1 : NormedAlgebra ℚ_[2] K1 :=
  spectralNorm.normedAlgebra ℚ_[2] K1

@[implicit_reducible] noncomputable instance nontriviallyK1 : NontriviallyNormedField K1 where
  toNormedField := normedFieldK1
  non_trivial := by
    refine ⟨algebraMap ℚ_[2] K1 ((2 : ℚ_[2])⁻¹), ?_⟩
    rw [norm_algebraMap_K1, norm_inv, show (2 : ℚ_[2]) = ((2 : ℕ) : ℚ_[2]) by norm_num,
      Padic.norm_p, inv_inv]
    norm_num

noncomputable instance properK1 : ProperSpace K1 := FiniteDimensional.proper ℚ_[2] K1

example : CompleteSpace K1 := by infer_instance

/-- ★`X⁴ − ϖ` は `K1` 上既約（抽象核 1）。 -/
instance factIrrM4 : Fact (Irreducible (X ^ 4 - C varpi)) :=
  ⟨X_pow_four_sub_C_irreducible varpi_not_sq neg_varpi_not_sq⟩

/-- ★★模型の台 `M4 = ℚ₂(i)((1+i)^{1/4})`（`K1` 上 4 次巡回・全分岐）。 -/
abbrev M4 : Type := AdjoinRoot (X ^ 4 - C varpi)

instance : FiniteDimensional K1 M4 := (AdjoinRoot.powerBasis (poly_ne_zero 4 varpi)).finite

@[implicit_reducible] noncomputable instance normedFieldM4 : NormedField M4 :=
  spectralNorm.normedField K1 M4

noncomputable instance ultraM4 : IsUltrametricDist M4 :=
  IsUltrametricDist.isUltrametricDist_of_forall_norm_add_le_max_norm
    (isNonarchimedean_spectralNorm (K := K1) (L := M4))

theorem norm_algebraMap_M4 (x : K1) : ‖algebraMap K1 M4 x‖ = ‖x‖ := spectralNorm_extends x

end TopField

section Uniformizer

open ConcreteNormedModel

/-- ★素元 `π = ϖ^{1/4} = (1+i)^{1/4}`。 -/
noncomputable def pi4 : M4 := AdjoinRoot.root (X ^ 4 - C varpi)

theorem pi4_pow : pi4 ^ 4 = algebraMap K1 M4 varpi := root_pow_eq (n := 4) varpi

theorem pi4_ne_zero : pi4 ≠ 0 := root_X_pow_sub_C_ne_zero (by norm_num) _

theorem norm_pi4_pow : ‖pi4‖ ^ (4 : ℕ) = ‖varpi‖ := by
  rw [← norm_pow, pi4_pow, norm_algebraMap_M4]

theorem norm_pi4_pos : 0 < ‖pi4‖ := norm_pos_iff.mpr pi4_ne_zero

theorem norm_pi4_lt_one : ‖pi4‖ < 1 := by
  by_contra hcon
  rw [not_lt] at hcon
  have h4 : (1 : ℝ) ≤ ‖pi4‖ ^ (4 : ℕ) := one_le_pow₀ hcon
  rw [norm_pi4_pow] at h4
  linarith [norm_varpi_lt_one]

theorem norm_pi4_ne_one : ‖pi4‖ ≠ 1 := ne_of_lt norm_pi4_lt_one

theorem norm_two_M4 : ‖((2 : ℕ) : M4)‖ = ((2 : ℕ) : ℝ)⁻¹ := by
  rw [← map_natCast (algebraMap K1 M4) 2, norm_algebraMap_M4]
  simpa using norm_two_K1

/-- ★`heM` —— `‖2‖ = ‖π‖^{p^{k+1}·e} = ‖π‖^{4·2} = ‖π‖^8`（`e = 2`）。 -/
theorem heM4 : ‖((2 : ℕ) : M4)‖ = ‖pi4‖ ^ (2 ^ (1 + 1) * 2) := by
  rw [norm_two_M4]
  have h8 : ‖pi4‖ ^ (2 ^ (1 + 1) * 2) = (‖pi4‖ ^ (4 : ℕ)) ^ (2 : ℕ) := by
    rw [← pow_mul]
    norm_num
  rw [h8, norm_pi4_pow, ← norm_pow, varpi_sq, norm_mul, norm_iota, norm_two_K1]
  norm_num

theorem finrank_M4 : Module.finrank K1 M4 = 2 ^ (1 + 1) := by
  rw [finrank_kummer (n := 4) varpi]
  norm_num

/-- ★`hvalK` —— `K1` の値群は `‖π‖^{4ℤ}`。 -/
theorem valK4 (z : K1) (hz : z ≠ 0) :
    ∃ m : ℤ, ‖algebraMap K1 M4 z‖ = ‖pi4‖ ^ (((2 ^ (1 + 1) : ℕ) : ℤ) * m) := by
  obtain ⟨m, hm⟩ := val_K1 hz
  refine ⟨m, ?_⟩
  rw [norm_algebraMap_M4, hm, ← norm_pi4_pow, ← zpow_natCast (‖pi4‖) 4, ← zpow_mul]
  norm_num

end Uniformizer

section Aut4

open ConcreteNormedModel GainedTowerModel.GaloisTower JumpFromValueGroup

/-- ★★位数 4 の自己同型 `g : π ↦ i·π`（`k = 1` の巡回性）。 -/
noncomputable def g4 : M4 ≃ₐ[K1] M4 := kummerAut iota_prim varpi

theorem g4_pi4 : g4 pi4 = algebraMap K1 M4 iota * pi4 := kummerAut_root iota_prim varpi

theorem two_ne_zero_M4 : ((2 : ℕ) : M4) ≠ 0 := by
  intro h
  have := norm_two_M4
  rw [h, norm_zero] at this
  norm_num at this

/-- ★`g² : π ↦ −π`（`i² = −1`）。 -/
theorem g4_sq_pi4 : (g4 ^ 2) pi4 = -pi4 := by
  have h : (g4 ^ 2) pi4 = g4 (g4 pi4) := by
    rw [show (g4 ^ 2) = g4 * g4 from sq g4, AlgEquiv.mul_apply]
  rw [h, g4_pi4, map_mul, AlgEquiv.commutes, g4_pi4, ← mul_assoc, ← map_mul, ← sq, iota_sq]
  simp

theorem g4_sq_ne_one : ¬ (g4 ^ 2 ^ 1 = 1) := by
  intro h
  simp only [pow_one] at h
  have h2 : (g4 ^ 2) pi4 = pi4 := by rw [h]; rfl
  rw [g4_sq_pi4] at h2
  refine pi4_ne_zero ?_
  have h3 : ((2 : ℕ) : M4) * pi4 = 0 := by push_cast; linear_combination -h2
  rcases mul_eq_zero.mp h3 with h4 | h4
  · exact absurd h4 two_ne_zero_M4
  · exact h4

/-- ★★`orderOf g = 4 = p^{k+1}`（`k = 1`）。 -/
theorem orderOf_g4 : orderOf g4 = 2 ^ (1 + 1) :=
  orderOf_eq_prime_pow g4_sq_ne_one (by
    rw [show (2 : ℕ) ^ (1 + 1) = 4 from rfl]
    exact kummerAut_pow iota_prim varpi)

theorem norm_iota_sub_one : ‖iota - 1‖ = ‖varpi‖ := by
  have h : iota - 1 = iota * varpi := by unfold varpi; linear_combination -iota_sq
  rw [h, norm_mul, norm_iota, one_mul]

/-- ★第 0 跳び: `‖gπ − π‖ = ‖π‖⁵` ⇒ `u₀ = 4`。 -/
theorem g4_pi4_sub : ‖g4 pi4 - pi4‖ = ‖pi4‖ ^ (5 : ℕ) := by
  have h : g4 pi4 - pi4 = algebraMap K1 M4 (iota - 1) * pi4 := by
    rw [g4_pi4, map_sub, map_one]; ring
  rw [h, norm_mul, norm_algebraMap_M4, norm_iota_sub_one, ← norm_pi4_pow]
  ring

/-- ★第 1 跳び: `‖g²π − π‖ = ‖π‖⁹` ⇒ `u₁ = 8`。 -/
theorem g4_sq_pi4_sub : ‖(g4 ^ 2) pi4 - pi4‖ = ‖pi4‖ ^ (9 : ℕ) := by
  have h : (g4 ^ 2) pi4 - pi4 = -(((2 : ℕ) : M4) * pi4) := by
    rw [g4_sq_pi4]; push_cast; ring
  rw [h, norm_neg, norm_mul, heM4]
  norm_num
  ring

theorem u0_eq_four (u : ℕ → ℤ) (hu : ‖g4 pi4 - pi4‖ = ‖pi4‖ ^ (u 0 + 1)) : u 0 = 4 := by
  have h5 : ‖pi4‖ ^ (u 0 + 1) = ‖pi4‖ ^ ((5 : ℕ) : ℤ) := by
    rw [← hu, g4_pi4_sub, zpow_natCast]
  have h6 := zpow_right_injective₀ norm_pi4_pos norm_pi4_ne_one h5
  omega

theorem u1_eq_eight (u : ℕ → ℤ) (hu : ‖(g4 ^ 2) pi4 - pi4‖ = ‖pi4‖ ^ (u 1 + 1)) : u 1 = 8 := by
  have h9 : ‖pi4‖ ^ (u 1 + 1) = ‖pi4‖ ^ ((9 : ℕ) : ℤ) := by
    rw [← hu, g4_sq_pi4_sub, zpow_natCast]
  have h10 := zpow_right_injective₀ norm_pi4_pos norm_pi4_ne_one h9
  omega

noncomputable def tau4 : M4 →+ M4 := AddMonoidHom.mk' (fun z => g4 z) (fun x y => map_add g4 x y)

noncomputable def s4 (j : ℕ) : M4 →+* M4 := ((g4 ^ 2 ^ j).toAlgHom : M4 →ₐ[K1] M4).toRingHom

end Aut4

section Instantiate

open ConcreteNormedModel GainedTowerModel.GaloisTower JumpFromValueGroup

/-- ★★★★**`harith` の 4 条件が `k = 1` で成り立つ** —— 中段 2 条件が
`k = 0` と違って**空虚でない**。

* `u₀ = 4`, `u₁ = 8`
* 増加: `4 < 8`
* ★★**Hasse–Arf の合同**: `p^{m+1} = 2 ∣ u₁ − u₀ = 4`（★具体例で初めて試された）
* 上界: `(p−1)·u₁ = 8 ≤ p^{k+1}·e = 4·2 = 8`（★等号でぎりぎり通る） -/
theorem harith_k1 : ∀ u : ℕ → ℤ, (∀ j, j ≤ 1 → ‖(s4 j) pi4 - pi4‖ = ‖pi4‖ ^ (u j + 1)) →
    1 ≤ u 0 ∧ (∀ m, m < 1 → u m < u (m + 1)) ∧
      (∀ m, m < 1 → (2 : ℤ) ^ (m + 1) ∣ u (m + 1) - u m) ∧
      ((2 : ℤ) - 1) * u 1 ≤ (2 : ℤ) ^ (1 + 1) * (2 : ℤ) := by
  intro u hu
  have h0 : ‖g4 pi4 - pi4‖ = ‖pi4‖ ^ (u 0 + 1) := by
    simpa [s4] using hu 0 (by norm_num)
  have h1 : ‖(g4 ^ 2) pi4 - pi4‖ = ‖pi4‖ ^ (u 1 + 1) := by
    simpa [s4] using hu 1 le_rfl
  have hv0 := u0_eq_four u h0
  have hv1 := u1_eq_eight u h1
  refine ⟨by omega, ?_, ?_, ?_⟩
  · intro m hm
    interval_cases m
    simp only [zero_add]
    omega
  · intro m hm
    interval_cases m
    simp only [zero_add, pow_one]
    rw [hv0, hv1]
    norm_num
  · rw [hv1]
    norm_num

/-- ★★★★★**`k = 1`（`p²` 次）の具体模型の出口** ——
`ℚ₂(i)((1+i)^{1/4}) / ℚ₂(i)`（4 次巡回・全分岐、`p = 2`, `k = 1`, `e = 2`）。 -/
theorem model_k1_exists (x : M4) :
    ∃ y : twr g4 2 1, ‖x - algebraMap (twr g4 2 1) M4 y‖
      ≤ (∏ j ∈ Finset.Icc 1 (1 + 1), axDecay 2 j) * ‖tau4 x - x‖ :=
  exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree
    (p := 2) (e := 2) (k := 1) (π := pi4) g4 orderOf_g4
    tau4 (fun _ => rfl) s4 (fun _ _ => rfl) (by norm_num)
    harith_k1
    finrank_M4 (by simpa using valK4) norm_two_M4 heM4 (adjoin_root_top (n := 4) varpi) x

end Instantiate

section Teeth

/-! ## ★★中段 2 条件が「歯を持つ」ことの確認

`k = 0` では `m < 0` が空なので中段 2 条件は**どんな `u` でも通る**。
`k = 1` では `m = 0` が現れ、実際の跳び `(u₀, u₁) = (4, 8)` に対して
`2 ∣ u₁ − u₀` を証明する必要がある。 -/

/-- ★`k = 0` では中段 2 条件は空虚（`k = 0` の 2 本が試していなかったもの）。 -/
theorem middle_vacuous_at_k_zero (u : ℕ → ℤ) :
    (∀ m, m < 0 → u m < u (m + 1)) ∧ (∀ m, m < 0 → (2 : ℤ) ^ (m + 1) ∣ u (m + 1) - u m) :=
  ⟨fun m hm => absurd hm (Nat.not_lt_zero m), fun m hm => absurd hm (Nat.not_lt_zero m)⟩

/-- ★★合同は**空でない制約**である —— 偶奇の違う値（例 `u₁ = 7`）なら落ちる。 -/
theorem congruence_has_teeth : ¬ ((2 : ℤ) ∣ (7 : ℤ) - 4) := by decide

/-- ★実際の跳び `(u₀, u₁) = (4, 8)` が満たすもの。★上界は**等号**で乗る。 -/
theorem jumps_k1_numbers :
    (1 : ℤ) ≤ 4 ∧ (4 : ℤ) < 8 ∧ (2 : ℤ) ^ (0 + 1) ∣ (8 : ℤ) - 4 ∧
      ((2 : ℤ) - 1) * 8 = (2 : ℤ) ^ (1 + 1) * 2 := by
  refine ⟨by norm_num, by norm_num, ?_, by norm_num⟩
  decide

end Teeth

/-! ## 使っている公理の一覧 -/

#print axioms X_pow_four_sub_C_irreducible
#print axioms exists_zpow_norm_of_sq_rel
#print axioms not_pow_of_norm_zpow
#print axioms neg_one_not_sq
#print axioms val_K1
#print axioms iota_prim
#print axioms orderOf_g4
#print axioms g4_pi4_sub
#print axioms g4_sq_pi4_sub
#print axioms harith_k1
#print axioms model_k1_exists
#print axioms jumps_k1_numbers

end ConcreteNormedModelK1
end ABC3.Found.PGC
