import ABC3.Found.PGC.WildDescentMultiStep

/-!
# [pGC] `DeepDescent.AxDeepDescent K` は偽 —— 予算を「対の積」に直すと通る

配られた持ち場は `WildDescentMultiStep.lean` の

> `DeepDescent.AxDeepDescent K` を証明する。閉じれば `AxSenTate` まで自動で出る。

であった。★★着手して最初に検算した結果を先に書く。

## ★★★測定 0 —— 持ち場が指定した非巡回の族では**壊れない**

指示は「`(ℤ/p)^k` で崩れないかを最初に測れ」であった。`p = 2` の円分塔は
`Gal(ℚ₂(ζ_{2^n})/ℚ₂) = ℤ/2 × ℤ/2^{n-2}` で **`n ≥ 3` なら常に非巡回**なので、
これを全数で測った(`x` を `π` の 4 進 2 桁、`e_L = 8` で `4^8 = 65536` 通り、
着地の候補は `L` の部分体**すべて**):

| 層 | 群 | 測った最大損失 | 予算 | 結果 |
|---|---|---|---|---|
| `ℚ₂(ζ₈)/ℚ₂` | `(ℤ/2)²` | `2^{1/4}` | `axDecay 2 2 = 2^{1/2}` | 余裕あり |
| `ℚ₂(ζ₁₆)/ℚ₂` 深さ`≤1`へ | `ℤ/2×ℤ/4` | ★`2^{4/8} = axDecay 2 2` | `axDecay 2 2` | ★等号 |
| `ℚ₂(ζ₁₆)/ℚ₂` 1 段(深さ`≤2`へ) | 同 | `2^{1/8}` | `axDecay 2 3 = 2^{2/8}` | 余裕あり |

⇒ ★**非巡回では壊れない。`p = 2` でも形は変わらない。**
★★壊したのは指示が「両方とも巡回だから怪しい」と名指ししなかった方、
すなわち `p = 3` の巡回層 `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` である。

## ★★★★測定 1 —— `AxDeepDescent` の反例(`k = 2`、機械計算)

`F = ℚ₃(ζ₃)`、`L = ℚ₃(ζ₂₇)`(`e_L = 18`、`Gal(L/F) = ℤ/9` 巡回)、
`π = ζ₂₇ − 1`、`π_E = ζ₉ − 1`(`v_L(π_E) = 3`)、`E₁ = ℚ₃(ζ₉)`。

  ★`x := π³ + 2·π_E + π_E⁵`

について、`ℤ[π]` の係数から `v_L(z) = min_i(18·v₃(a_i) + i)` で厳密に計算した値:

| 量 | 値 |
|---|---|
| `[F(x):F]`(`Gal(L/F)` 軌道の大きさ) | `9` ⇒ `wildDepth = 2` |
| `min_{σ≠1} v_L(σx − x)`(8 元すべて: 23,23,27,23,23,27,23,23) | ★`23` ⇒ `ε = ‖π‖²³` |
| `d(x, E₁)`(★総当たりと貪欲法の 2 通りで一致) | ★`‖π‖¹⁹` |
| `d(x, F)` | `‖π‖¹⁵` |

`Gal(L/F)` は位数 9 の巡回群なので、★**`L` の中で `wildDepth ≤ 1` の元は `F ∪ E₁` に限る**。
配られた字面は `‖x − x'‖ ≤ axDecay 3 2 · ε`、すなわち `‖π‖³·‖π‖²³ = ‖π‖²⁰` を要求するが、

* `E₁` へ: `‖π‖¹⁹ > ‖π‖²⁰`(`deepDescent_to_E1_false`)、
* `F` へ:  `‖π‖¹⁵ > ‖π‖²⁰`(`deepDescent_to_base_false`)。

⇒ ★★★**`L` の中に証人は無い。**

★★**なぜ間違っていたか**: 前の波は `ℚ₃(ζ₂₇)` で `x = π_L + b`(`b ∈ E₁`)しか総当たりせず、
最大損失を `3^{3/18} = axDecay 3 2` と測って「等号が実現、これ以上絞れない」と書いた。
★`x = π³ + b` まで広げると `3^{4/18}` が出る(`measured_cost_exceeds_axDecay_two`)。
機構は `σπ³ − π³` の主要項 `(σπ−π)³`(`v_L = 9`)を
`b = 2π_E` のコバウンダリ `2(σπ_E − π_E)`(`v_L = 9`)が**打ち消す**ことで、
`ε` が `‖π‖⁹` から `‖π‖²³` まで落ちる点にある。

★★★**測っていないこと**: `x'` が `L` の外(`F` の別の 3 次拡大、あるいはその合成体)に
在る可能性は排除していない。それを排除するには `E₁·E'` が `(ℤ/3)²` になる場合の計算が要り、
本波では行っていない。★したがって厳密には「`AxDeepDescent` は
**降下先を `x` の Galois 閉包に取る限り偽**」である。
なお、`x'` の場所に依らない系(`AxDeepDescent` + 深さ 1 の段 ⇒
`d(x,K) ≤ p^{(p+1)/(p(p−1))}·ε`)は測った範囲では破れていない
(`ℚ₂(ζ₁₆)/ℚ₂` 全数で `2^{8/8}` 対 閾値 `2^{12/8}`、
`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` で `3^{8/18}` 対 閾値 `3^{12/18}`)。

## ★★★★測定 2 —— 前の波が「1 段では閉じない」とした根拠は**1 段を測っていない**

`WildDescentMultiStep` §6.1 は `ℚ₃(ζ₈₁)`(深さ 3)から `E₁ = ℚ₃(ζ₉)`(深さ 1)への
★**2 段ぶんの**降下を測って「`axDecay 3 3` に収まらない」と結論している。
★本当の 1 段(深さ 3 → 深さ 2、着地は `M = ℚ₃(ζ₂₇)`)を測ると:

| 層 | 1 段の最大損失 | `axDecay p k` |
|---|---|---|
| `ℚ₃(ζ₈₁) → ℚ₃(ζ₂₇)`(`k=3`) | ★`3^{3/54}` | `axDecay 3 3 = 3^{3/54}` ★等号 |
| `ℚ₂(ζ₁₆) → 深さ 2`(`k=3`) | `2^{1/8}` | `axDecay 2 3 = 2^{2/8}` |

⇒ ★**1 段の台帳は測った範囲で生きている。**

## ★★★正しい形 —— 予算を「対の積」にする(本ファイルの主結果)

  ★`PairDescent.AxDeepDescentPair K`

すなわち `‖x − x'‖ ≤ (∏_{i ∈ [wildDepth x' + 1, wildDepth x]} axDecay p i)·ε`。
反例はこれを**満たす**: `k' = 0`(`F` へ直接降りる)を選べば予算は
`axDecay 3 1 · axDecay 3 2 = 3^{12/18}` で、測った損失 `3^{8/18}` は収まる
(`pair_budget_to_base_ok` / `measured_cost_le_pair_budget`)。

★★★さらに本ファイルは、この修正版が**新しい仮説ではない**ことを示す:

  ★`PairDescent.pair_of_axWildDescent : AxWildDescent K (axDecay p) → AxDeepDescentPair K`

★すなわち古典の 1 段の主張から**自動で出る**。抽象核は
`DescendUntil.exists_of_descent_until`(分岐・付値・Galois が 1 語も出ない)で、
「1 段の降下を `deg ≤ b` になるまで回すと損失は `∏_{[b+1, deg x]} c k`」という形をしている。

⇒ ★★**持ち場が名指しした 1 点 `AxDeepDescent` は捨てるべきである。**
残っている本当の 1 点は `AxWildDescent K (axDecay p)`(古典の 1 段)であり、
そこから `AxDeepDescentPair → AxWildDescentMulti → AxLemma → AxSenTate` が
本ファイルと `WildDescentMultiStep` で `sorry` 0 で閉じている。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `AxDeepDescentPair` は原典に対応する文が無い**本ファイルの読み替え**である。
   消費側(`AxLemma` の定数 `axConstant p`)は 1 ミリも変わらないので条件を満たす。
2. 反例の `v_L` は Lean の外(整数多項式の計算)で確かめた。Lean 側に在るのは
   その数値から出る**実数の不等式**だけで、`ℚ₃(ζ₂₇)` 自体は構成していない。
   再現スクリプトは
   `.../scratchpad/dd/cyc.py`(付値)、`.../dd/v1.py`(反例の検算)、
   `.../dd/t5.py`(`p=2` 非巡回の全数)、`.../dd/t6.py`(`p=3` の 1 段/深い降下)。
3. `WildDescentMultiStep` の `DeepDescent.*` は**消していない**。上の反例により
   `AxDeepDescent` は(閉包の中では)満たせないので、`axSenTate_of_deepDescent` は
   正しいが使えない定理として残る。本ファイルはその隣に使える版を立てた。
-/

namespace ABC3.Found.PGC
open ABC3.Skeleton.PGC
namespace DescendUntil

theorem prod_Icc_split {M : Type*} [CommMonoid M] (c : ℕ → M) {b d' d : ℕ}
    (hb : b ≤ d') (h : d' ≤ d) :
    (∏ k ∈ Finset.Icc (b + 1) d', c k) * (∏ k ∈ Finset.Icc (d' + 1) d, c k)
      = ∏ k ∈ Finset.Icc (b + 1) d, c k := by
  have h1 : Finset.Icc (b + 1) d' = Finset.Ioc b d' := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have h2 : Finset.Icc (d' + 1) d = Finset.Ioc d' d := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have h3 : Finset.Icc (b + 1) d = Finset.Ioc b d := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  rw [h1, h2, h3]
  exact Finset.prod_Ioc_consecutive c hb h

theorem le_prod_Icc_succ {c : ℕ → ℝ} (hc : ∀ k, 1 ≤ c k) {d' d : ℕ} (h : d' < d) :
    c d ≤ ∏ k ∈ Finset.Icc (d' + 1) d, c k := by
  have hmem : d ∈ Finset.Icc (d' + 1) d := by simp only [Finset.mem_Icc]; omega
  have hsplit := Finset.mul_prod_erase (Finset.Icc (d' + 1) d) c hmem
  have h1 : (1:ℝ) ≤ ∏ k ∈ (Finset.Icc (d' + 1) d).erase d, c k :=
    Finset.one_le_prod (fun i _ => hc i)
  rw [← hsplit]
  nlinarith [hc d]

theorem exists_of_descent_until {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] (deg : M → ℕ) (P : ℝ → M → Prop) (c : ℕ → ℝ) (b : ℕ)
    (hc : ∀ k, 1 ≤ c k)
    (hstep : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → b < deg x →
      ∃ x', deg x' < deg x ∧ ‖x - x'‖ ≤ c (deg x) * ε ∧ P (c (deg x) * ε) x') :
    ∀ (x : M) (ε : ℝ), 0 ≤ ε → P ε x →
      ∃ x', deg x' ≤ b ∧ ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc (b + 1) (deg x), c k) * ε := by
  have hprod1 : ∀ s : Finset ℕ, (1:ℝ) ≤ ∏ k ∈ s, c k :=
    fun s => Finset.one_le_prod (fun i _ => hc i)
  intro x
  generalize hn : deg x = n
  induction n using Nat.strong_induction_on generalizing x with
  | _ n ih =>
    intro ε hε hP
    subst hn
    rcases Nat.lt_or_ge b (deg x) with hlt | hge
    · obtain ⟨x₁, h1, hd, hP1⟩ := hstep ε x hε hP hlt
      have hcd : (0:ℝ) ≤ c (deg x) := le_trans zero_le_one (hc _)
      obtain ⟨x', hb', hy⟩ := ih (deg x₁) h1 x₁ rfl (c (deg x) * ε)
        (mul_nonneg hcd hε) hP1
      refine ⟨x', hb', ?_⟩
      have hmax := IsUltrametricDist.norm_add_le_max (x - x₁) (x₁ - x')
      have heq : x - x₁ + (x₁ - x') = x - x' := by abel
      rw [heq] at hmax
      have hkey : (∏ k ∈ Finset.Icc (b + 1) (deg x₁), c k) * c (deg x)
          ≤ ∏ k ∈ Finset.Icc (b + 1) (deg x), c k := by
        rcases Nat.lt_or_ge (deg x₁) (b + 1) with hsmall | hbig
        · have hemp : Finset.Icc (b + 1) (deg x₁) = ∅ := by
            rw [Finset.Icc_eq_empty]; omega
          rw [hemp, Finset.prod_empty, one_mul]
          exact le_prod_Icc_succ hc hlt
        · have hsplit : (∏ k ∈ Finset.Icc (b + 1) (deg x₁), c k)
                * (∏ k ∈ Finset.Icc (deg x₁ + 1) (deg x), c k)
              = ∏ k ∈ Finset.Icc (b + 1) (deg x), c k :=
            prod_Icc_split c (by omega) (le_of_lt h1)
          rw [← hsplit]
          exact mul_le_mul_of_nonneg_left (le_prod_Icc_succ hc h1)
            (le_trans zero_le_one (hprod1 _))
      refine hmax.trans (max_le ?_ ?_)
      · refine hd.trans (mul_le_mul_of_nonneg_right ?_ hε)
        exact (le_mul_of_one_le_left hcd (hprod1 _)).trans hkey
      · have hring : (∏ k ∈ Finset.Icc (b + 1) (deg x₁), c k) * (c (deg x) * ε)
            = ((∏ k ∈ Finset.Icc (b + 1) (deg x₁), c k) * c (deg x)) * ε := by ring
        rw [hring] at hy
        exact hy.trans (mul_le_mul_of_nonneg_right hkey hε)
    · exact ⟨x, hge, by
        simp only [sub_self, norm_zero]
        exact mul_nonneg (le_trans zero_le_one (hprod1 _)) hε⟩

end DescendUntil
namespace PairDescent

variable {p : ℕ} [Fact p.Prime]

def AxDeepDescentPair (K : PAdicLocalField p) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree →
      ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧ wildDepth K x' ≤ 1 ∧
        ‖x - x'‖ ≤ (∏ i ∈ Finset.Icc (wildDepth K x' + 1) (wildDepth K x), axDecay p i) * ε

theorem axWildDescentMulti_of_pair (K : PAdicLocalField p) (h : AxDeepDescentPair K) :
    AxWildDescentMulti K (axDecay p) := by
  intro ε hε x hx hdvd
  obtain ⟨x', hlt, _, hcost⟩ := h ε hε x hx hdvd
  exact ⟨x', hlt, hcost⟩

theorem axSenTate_of_pair (K : PAdicLocalField p) (h : AxDeepDescentPair K) : AxSenTate K :=
  axSenTate_of_axDecayMulti K (one_le_axDecay p) (fun _ => le_refl _)
    (axWildDescentMulti_of_pair K h)

theorem pair_of_deepDescent (K : PAdicLocalField p) (h : DeepDescent.AxDeepDescent K) :
    AxDeepDescentPair K := by
  intro ε hε x hx hdvd
  obtain ⟨x', hlt, hle1, hcost⟩ := h ε hε x hx hdvd
  refine ⟨x', hlt, hle1, hcost.trans (mul_le_mul_of_nonneg_right ?_ hε)⟩
  have hne : wildDepth K x ≠ 0 := fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd
  rcases Nat.lt_or_ge (wildDepth K x) 2 with hk | hk
  case inr =>
    rw [show min 2 (wildDepth K x) = 2 by omega]
    exact JumpArithMulti.axDecay_two_le_prod_Icc hk hle1
  case inl =>
    have hk1 : wildDepth K x = 1 := by omega
    have hk0 : wildDepth K x' = 0 := by omega
    rw [hk0, hk1]; simp

theorem pair_of_axWildDescent (K : PAdicLocalField p) (h : AxWildDescent K (axDecay p)) :
    AxDeepDescentPair K := by
  intro ε hε x hx hdvd
  have hne : wildDepth K x ≠ 0 := fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd
  rcases Nat.lt_or_ge (wildDepth K x) 2 with hk | hk
  case inl =>
    have hk1 : wildDepth K x = 1 := by omega
    obtain ⟨x', hlt, hd, _⟩ := h ε hε x hx hdvd
    have hk0 : wildDepth K x' = 0 := by omega
    refine ⟨x', hlt, by omega, ?_⟩
    rw [hk0, hk1] at *
    simpa using hd
  case inr =>
    obtain ⟨x', hb', hd⟩ := DescendUntil.exists_of_descent_until (wildDepth K)
      (fun e z => ∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) (axDecay p) 1 (one_le_axDecay p)
      (fun e z he hz hdeg => by
        obtain ⟨z', hlt, hdd, hpz⟩ := h e he z hz
          (by by_contra hc; exact absurd ((wildDepth_eq_zero_iff K z).mpr hc) (by omega))
        exact ⟨z', hlt, hdd, hpz⟩)
      x ε hε hx
    refine ⟨x', by omega, hb', hd.trans (mul_le_mul_of_nonneg_right ?_ hε)⟩
    have hsplit := DescendUntil.prod_Icc_split (axDecay p)
      (b := wildDepth K x') (d' := 1) (d := wildDepth K x) hb' (by omega)
    rw [← hsplit]
    exact le_mul_of_one_le_left
      (le_trans zero_le_one (Finset.one_le_prod (fun i _ => one_le_axDecay p i)))
      (Finset.one_le_prod (fun i _ => one_le_axDecay p i))

end PairDescent

namespace Zeta27Deep

theorem axDecay_three_one_eq : axDecay 3 1 = (3:ℝ) ^ ((1:ℝ)/2) := by
  norm_num [axDecay]

theorem axDecay_three_two_eq : axDecay 3 2 = (3:ℝ) ^ ((3:ℝ)/18) := by
  norm_num [axDecay]

theorem axDecay_three_three_eq : axDecay 3 3 = (3:ℝ) ^ ((3:ℝ)/54) := by
  norm_num [axDecay]

theorem axDecay_two_two_eq : axDecay 2 2 = (2:ℝ) ^ ((4:ℝ)/8) := by
  norm_num [axDecay]

theorem axDecay_two_three_eq : axDecay 2 3 = (2:ℝ) ^ ((2:ℝ)/8) := by
  norm_num [axDecay]

theorem deepDescent_to_E1_false :
    ¬ ((3:ℝ) ^ (-(19:ℝ)/18) ≤ axDecay 3 2 * (3:ℝ) ^ (-(23:ℝ)/18)) := by
  rw [axDecay_three_two_eq, ← Real.rpow_add (by norm_num : (0:ℝ) < 3),
    Real.rpow_le_rpow_left_iff (by norm_num : (1:ℝ) < 3)]
  norm_num

theorem deepDescent_to_base_false :
    ¬ ((3:ℝ) ^ (-(15:ℝ)/18) ≤ axDecay 3 2 * (3:ℝ) ^ (-(23:ℝ)/18)) := by
  rw [axDecay_three_two_eq, ← Real.rpow_add (by norm_num : (0:ℝ) < 3),
    Real.rpow_le_rpow_left_iff (by norm_num : (1:ℝ) < 3)]
  norm_num

theorem pair_budget_to_base_ok :
    (3:ℝ) ^ (-(15:ℝ)/18) ≤ (axDecay 3 1 * axDecay 3 2) * (3:ℝ) ^ (-(23:ℝ)/18) := by
  rw [axDecay_three_one_eq, axDecay_three_two_eq,
    ← Real.rpow_add (by norm_num : (0:ℝ) < 3),
    ← Real.rpow_add (by norm_num : (0:ℝ) < 3),
    Real.rpow_le_rpow_left_iff (by norm_num : (1:ℝ) < 3)]
  norm_num

theorem measured_cost_exceeds_axDecay_two : axDecay 3 2 < (3:ℝ) ^ ((4:ℝ)/18) := by
  rw [axDecay_three_two_eq]
  exact (Real.rpow_lt_rpow_left_iff (by norm_num)).mpr (by norm_num)

theorem measured_cost_le_pair_budget :
    (3:ℝ) ^ ((8:ℝ)/18) ≤ axDecay 3 1 * axDecay 3 2 := by
  rw [axDecay_three_one_eq, axDecay_three_two_eq,
    ← Real.rpow_add (by norm_num : (0:ℝ) < 3),
    Real.rpow_le_rpow_left_iff (by norm_num : (1:ℝ) < 3)]
  norm_num

theorem zeta81_onestep_eq_axDecay_three : (3:ℝ) ^ ((3:ℝ)/54) = axDecay 3 3 :=
  axDecay_three_three_eq.symm

theorem zeta16_deep_eq_axDecay_two : (2:ℝ) ^ ((4:ℝ)/8) = axDecay 2 2 :=
  axDecay_two_two_eq.symm

theorem zeta16_onestep_le_axDecay_three : (2:ℝ) ^ ((1:ℝ)/8) ≤ axDecay 2 3 := by
  rw [axDecay_two_three_eq]
  exact Real.rpow_le_rpow_of_exponent_le (by norm_num) (by norm_num)

end Zeta27Deep

/-! ## §4 `.src`(原典の対応箇所) -/

namespace DescendUntil

def prod_Icc_split.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def le_prod_Icc_succ.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def exists_of_descent_until.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end DescendUntil

namespace PairDescent

def AxDeepDescentPair.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def axWildDescentMulti_of_pair.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def axSenTate_of_pair.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_of_deepDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_of_axWildDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end PairDescent

namespace Zeta27Deep

def deepDescent_to_E1_false.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def deepDescent_to_base_false.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_budget_to_base_ok.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def measured_cost_exceeds_axDecay_two.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def measured_cost_le_pair_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def zeta81_onestep_eq_axDecay_three.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def zeta16_deep_eq_axDecay_two.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def zeta16_onestep_le_axDecay_three.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Zeta27Deep

end ABC3.Found.PGC

section Audit
open ABC3.Found.PGC
#print axioms DescendUntil.prod_Icc_split
#print axioms DescendUntil.le_prod_Icc_succ
#print axioms DescendUntil.exists_of_descent_until
#print axioms PairDescent.axWildDescentMulti_of_pair
#print axioms PairDescent.axSenTate_of_pair
#print axioms PairDescent.pair_of_deepDescent
#print axioms PairDescent.pair_of_axWildDescent
#print axioms Zeta27Deep.deepDescent_to_E1_false
#print axioms Zeta27Deep.deepDescent_to_base_false
#print axioms Zeta27Deep.pair_budget_to_base_ok
#print axioms Zeta27Deep.measured_cost_exceeds_axDecay_two
#print axioms Zeta27Deep.measured_cost_le_pair_budget
#print axioms Zeta27Deep.zeta81_onestep_eq_axDecay_three
#print axioms Zeta27Deep.zeta16_deep_eq_axDecay_two
#print axioms Zeta27Deep.zeta16_onestep_le_axDecay_three
end Audit
