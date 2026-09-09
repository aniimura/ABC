import ABC3.Found.PGC.DeepLossExit
import ABC3.Found.PGC.GainedDescentBridge

/-!
# [pGC] 正規化トレースの「得」の抽象核 ＋ ★★前波の「残りは 1 本」の訂正

## ★★訂正（自己訂正 16 度目）—— 私の前波の言葉が「空白」を誇張していた

前波の報告で私はこう書いた（逐語）:

> **残りはちょうど 1 点**: … `AxWildDescent K c` を作ること ——
> ★**降下そのもの 1 本だけ**です。

★★**「1 本だけ」は正しいが、「そこが空白だ」という含みは誤り。**
本波で `grep` したところ、木は既に**体の層で `axDecay` の積まで到達した定理**を持っている:

* `GainedBridge.exists_norm_sub_algebraMap_le_prod_axDecay`
  （`GainedDescentBridge.lean:527`、名前空間 `ABC3.Found.PGC.GainedBridge`）
  結論（逐語）:
  `∃ y : F, ‖x - algebraMap F M y‖ ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖`
* `GainedBridge.exists_norm_sub_algebraMap_le_gainedLoss`（`:495`）
* `GainedBridgeSupply.exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps`
  （`GainedBridgeSupply.lean:474`）

★これは `AxLemma` の形（`∃ y ∈ F, ‖x − y‖ ≤ C·ε`）を**1 つの層 `M/F`（`[M:F] = p`）で**
達成している。★私の前波の 3 つのファイル（`BaseLayerConstant` / `DeepLossExit` /
`FiniteExceptionsIcc`）は `AxWildDescent`（**段**の形）の側から組んでいた。
★★**木には 2 本目の道が既にある**（`DeepDescentPairDirect.lean:306`
`axSenTate_of_axLemmaGraded (K) (h : AxLemmaGraded K) : AxSenTate K`）。

★測ったコマンド:
`grep -n '^namespace\|^theorem ' lean/ABC3/Found/PGC/GainedDescentBridge.lean`
`grep -n '^namespace' lean/ABC3/Found/PGC/GainedDescentBridge.lean` → `:142` / `:144`（#359）。

★★**私は前波まで、この 2 本のファイル（603 行 ＋ 561 行）を一度も測っていなかった。**
「残りは降下 1 本」と書いたとき、私は `AxWildDescent` 側の在庫しか見ていなかった。

## ★★★訂正 2 —— 本ファイルの下書きの断定も偽だった（commit 前に測り直した）

下書きの私はこう書いていた:

> ★★**正規化トレース `Σ_{i<n} σ^i x`（`n` 個の和）については木に無い**。
> 測定: `grep -rn "^theorem norm_sum_iterate" lean/ABC3/Found/PGC/*.lean` → **0 件**。

★★**そのコマンドを実際に叩いたら 0 件ではなく 4 件だった。** 私は grep を書いただけで
**実行していなかった**（★「測っていない」を「測って無かった」と書く典型の誤り）。

在るもの（`CyclicLayerDescent.lean`、名前空間 `ABC3.Found.PGC`）:

* `:412 norm_sum_iterate_sub_self_le_max`
  `‖∑_{j<p}(((D+1)^j) x − x)‖ ≤ max ‖(D^{p−1}) x‖ (c * ‖D x‖)`
* `:461 norm_sum_iterate_sub_self_le_of_contract`
  `‖∑_{j<p}(((D+1)^j) x − x)‖ ≤ max (θ^{p−2}) c * ‖D x‖`

★`σ = D + 1` と置けば `∑_{j<p}(σ^j x − x) = (∑_{j<p} σ^j x) − p•x` で、
★★**本ファイル §1 と同じ量**である。

★★しかも**木の方が良い**: 木は Newton の二項展開（`p ∣ C(p,l)`, `0 < l < p`）で
主要項を `θ^{p−2}` まで落としており、本ファイルのガウス和 `p(p−1)/2` の道は
`max (1/p) b` までしか落ちない。`θ ≤ 1` かつ `p ≥ 3` なら `θ^{p−2} ≤ θ` なので
★**木の評価は本ファイルの評価以上に良い**（§4 `tree_contract_le_gauss` で証明した）。

⇒ ★★**本ファイル §1〜§3 は「2 本目の証明」であって新しい数学ではない。**
★残す理由は 2 つだけ: ① ガウス和の道が二項展開の道と**独立に**同じ結論を出すことの記録、
② §4 の比較（どちらが良いかを**測った**もの）。★使うべきは木の方である。

## 本ファイルの内容（新規性の申告つき）

  `‖(Σ_{i<n} f^[i] x) − n•x‖ ≤ max a b · ‖f x − x‖`   （§1、★新規性なし）

* `a` = `‖(Σ_{i<n} i)•y‖/‖y‖` の上界（★`n = p` 奇素数なら `Σ_{i<p} i = p(p−1)/2` が
  `p` で割れるので **`a = ‖(p:M)‖ = 1/p`**、§2。★この「ガウス和の得」は木に無い形）
* `b` = `‖f^[j]y − y‖/‖y‖` の上界

正規化トレース `x' = (1/p)·Σ_{i<p} σ^i x` の 1 段の損失（§3）:

  `‖x' − x‖ ≤ max 1 (p·b) · ‖σ x − x‖`
  ⇒ ★`b ≤ 1/p` なら損失は `1`（ただ）。

## 逸脱の記録

- §1 は `M →+ M` の反復についての純粋な超距離評価で、★体も Galois も素数も出ない。
  ★`f` が全単射である必要すら無い。
- §2 は `p ≠ 2` を使う（`p = 2` では `Σ_{i<2} i = 1` で `2 ∤ 1`、★得が無い）。
  ★これは `p = 2` が特別だという既知の事実と整合する。§2 の `two_no_gain` に記録した。
- §3 は `NormedField` の層で `x' = (p:M)⁻¹ • (Σ …)` を作る。
  ★`GainedBridge` の道（桁展開）とは**別の道**であり、どちらが安いかは測っていない。
- ★★**新規性の申告**: §1〜§3 は `CyclicLayerDescent.lean:412/:461` と**同じ量**の評価で、
  ★そちらの方が良い（§4）。★本ファイルを新規結果として数えないこと。
  ★`prime_dvd_sum_range`（ガウス和が `p` で割れる）だけは木に無い形だが、
  それ単独では使い道が無い。
-/

namespace ABC3.Found.PGC

namespace TraceGainKernel

open Finset

/-! ## §1 ★抽象核 —— 正規化トレースの分解（体も Galois も素数も出ない） -/

section Kernel

variable {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]

/-- ★★**正規化トレースの抽象核**。

`y := f x − x` として

  `‖(Σ_{i<n} f^[i] x) − n•x‖ ≤ max a b · ‖y‖`

* `ha` : `‖(Σ_{i<n} i)•y‖ ≤ a·‖y‖` ——★**ガウス和 `n(n−1)/2` の「得」**
* `hb` : `∀ j < n, ‖f^[j] y − y‖ ≤ b·‖y‖` ——★分岐が要るのはここだけ

★木の `AxEpsilonDecay.norm_iterate_sub_self_le`（`:267`）と**同じ `max a b` の形**だが、
あちらは `f^[n] x − x`（1 個の像）、こちらは `Σ_{i<n} f^[i] x`（**トレース**）。
★★**トレース版は木にも在る**（`CyclicLayerDescent.lean:412` / `:461`、二項展開の道）。
★そちらの方が良い評価を出す（§4 `tree_contract_le_gauss`）。★本補題は独立な 2 本目の証明。 -/
theorem norm_sum_iterate_sub_nsmul_le (f : M →+ M) (x : M) (n : ℕ) {a b : ℝ}
    (ha : ‖(∑ i ∈ Finset.range n, i) • (f x - x)‖ ≤ a * ‖f x - x‖)
    (hb : ∀ j < n, ‖f^[j] (f x - x) - (f x - x)‖ ≤ b * ‖f x - x‖) :
    ‖(∑ i ∈ Finset.range n, f^[i] x) - n • x‖ ≤ max a b * ‖f x - x‖ := by
  set y := f x - x with hy
  have hy0 : (0 : ℝ) ≤ ‖y‖ := norm_nonneg _
  have hay : 0 ≤ a * ‖y‖ := le_trans (norm_nonneg _) ha
  have hmax : a * ‖y‖ ≤ max a b * ‖y‖ :=
    mul_le_mul_of_nonneg_right (le_max_left a b) hy0
  have hC : 0 ≤ max a b * ‖y‖ := le_trans hay hmax
  have h1 : (∑ i ∈ Finset.range n, f^[i] x) - n • x
      = ∑ i ∈ Finset.range n, (f^[i] x - x) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range]
  have h2 : ∀ i ∈ Finset.range n,
      f^[i] x - x = i • y + ∑ j ∈ Finset.range i, (f^[j] y - y) :=
    fun i _ => iterate_sub_self_eq_nsmul_add f x i
  have h3 : (∑ i ∈ Finset.range n, f^[i] x) - n • x
      = (∑ i ∈ Finset.range n, i) • y
        + ∑ i ∈ Finset.range n, ∑ j ∈ Finset.range i, (f^[j] y - y) := by
    rw [h1, Finset.sum_congr rfl h2, Finset.sum_add_distrib, ← Finset.sum_smul]
  rw [h3]
  refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le (ha.trans hmax) ?_)
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hC ?_
  intro i hi
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hC ?_
  intro j hj
  exact (hb j (lt_trans (Finset.mem_range.mp hj) (Finset.mem_range.mp hi))).trans
    (mul_le_mul_of_nonneg_right (le_max_right a b) hy0)

end Kernel

/-! ## §2 ガウス和の「得」—— `p` 奇素数なら `p ∣ Σ_{i<p} i` -/

section Gauss

/-- ★奇素数なら `p ∣ Σ_{i<p} i = p(p−1)/2`。★これが `a = ‖(p:L)‖` の出どころ。 -/
theorem prime_dvd_sum_range {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) :
    p ∣ ∑ i ∈ Finset.range p, i := by
  have h2 : (∑ i ∈ Finset.range p, i) * 2 = p * (p - 1) := Finset.sum_range_id_mul_two p
  have hdvd : p ∣ (∑ i ∈ Finset.range p, i) * 2 := ⟨p - 1, h2⟩
  rcases (Nat.Prime.dvd_mul hp).mp hdvd with h | h
  · exact h
  · exact absurd ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp h) hp2

/-- ★★`p = 2` では得が無い（`Σ_{i<2} i = 1` で `2 ∤ 1`）。

★「奇素数」の条件は本質的で、`p = 2` は別扱いが要る。
★`DeepLossExit` で測った「`p = 2` の最小 `K₀` は 3」と同じ向きの現象である。 -/
theorem two_no_gain : ¬ (2 ∣ ∑ i ∈ Finset.range 2, i) := by
  norm_num [Finset.sum_range_succ]

end Gauss

/-! ## §3 体の層 —— 正規化トレースの 1 段の損失 -/

section FieldLayer

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- §1 の `a` を ★**`‖(p:L)‖`** に確定させた形（`p` 奇素数）。 -/
theorem norm_sum_iterate_sub_nsmul_le_prime {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (f : L →+ L) (x : L) {b : ℝ}
    (hb : ∀ j < p, ‖f^[j] (f x - x) - (f x - x)‖ ≤ b * ‖f x - x‖) :
    ‖(∑ i ∈ Finset.range p, f^[i] x) - p • x‖ ≤ max ‖(p : L)‖ b * ‖f x - x‖ := by
  refine norm_sum_iterate_sub_nsmul_le f x p ?_ hb
  rw [nsmul_eq_mul, norm_mul]
  exact mul_le_mul_of_nonneg_right
    (GainedBridge.norm_natCast_le_of_dvd (prime_dvd_sum_range hp hp2)) (norm_nonneg _)

/-- ★★**正規化トレース `x' = (p:L)⁻¹·Σ_{i<p} f^[i] x` の 1 段の損失**。 -/
theorem norm_normalizedTrace_sub_le {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hp0 : (p : L) ≠ 0) (f : L →+ L) (x : L) {b : ℝ}
    (hb : ∀ j < p, ‖f^[j] (f x - x) - (f x - x)‖ ≤ b * ‖f x - x‖) :
    ‖(p : L)⁻¹ * (∑ i ∈ Finset.range p, f^[i] x) - x‖
      ≤ ‖(p : L)‖⁻¹ * max ‖(p : L)‖ b * ‖f x - x‖ := by
  have hrw : (p : L)⁻¹ * (∑ i ∈ Finset.range p, f^[i] x) - x
      = (p : L)⁻¹ * ((∑ i ∈ Finset.range p, f^[i] x) - p • x) := by
    rw [mul_sub, nsmul_eq_mul, ← mul_assoc, inv_mul_cancel₀ hp0, one_mul]
  rw [hrw, norm_mul, norm_inv, mul_assoc]
  exact mul_le_mul_of_nonneg_left (norm_sum_iterate_sub_nsmul_le_prime hp hp2 f x hb)
    (by positivity)

/-- ★★★**閉じた形** —— `‖(p:L)‖ = 1/p` なら 1 段の損失は `max 1 (p·b)`。

⇒ ★★**`b ≤ 1/p` なら損失は `1`（ただ）**。
⇒ 目標 `p^{(2p−2)/(m·e)}` は `b ≤ p^{(2p−2)/(m·e) − 1}` で達成される。
★これが前波で「残っている 1 点」と書いたものの**閉じた形**である。 -/
theorem norm_normalizedTrace_sub_le_max_one {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hnp : ‖(p : L)‖ = ((p : ℝ))⁻¹) (f : L →+ L) (x : L) {b : ℝ}
    (hb : ∀ j < p, ‖f^[j] (f x - x) - (f x - x)‖ ≤ b * ‖f x - x‖) :
    ‖(p : L)⁻¹ * (∑ i ∈ Finset.range p, f^[i] x) - x‖
      ≤ max 1 ((p : ℝ) * b) * ‖f x - x‖ := by
  have hpR : (0 : ℝ) < p := by exact_mod_cast hp.pos
  have hp0 : (p : L) ≠ 0 := by
    intro h
    rw [h, norm_zero] at hnp
    exact (inv_ne_zero (ne_of_gt hpR)) hnp.symm
  refine (norm_normalizedTrace_sub_le hp hp2 hp0 f x hb).trans ?_
  rw [hnp, inv_inv]
  refine mul_le_mul_of_nonneg_right (le_of_eq ?_) (norm_nonneg _)
  rw [mul_max_of_nonneg _ _ (le_of_lt hpR), mul_inv_cancel₀ (ne_of_gt hpR)]

/-- ★★**`b ≤ 1/p` なら 1 段はただ**（損失 `1`）。 -/
theorem norm_normalizedTrace_sub_le_of_small {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hnp : ‖(p : L)‖ = ((p : ℝ))⁻¹) (f : L →+ L) (x : L) {b : ℝ}
    (hbs : b ≤ ((p : ℝ))⁻¹)
    (hb : ∀ j < p, ‖f^[j] (f x - x) - (f x - x)‖ ≤ b * ‖f x - x‖) :
    ‖(p : L)⁻¹ * (∑ i ∈ Finset.range p, f^[i] x) - x‖ ≤ ‖f x - x‖ := by
  have hpR : (0 : ℝ) < p := by exact_mod_cast hp.pos
  refine (norm_normalizedTrace_sub_le_max_one hp hp2 hnp f x hb).trans ?_
  have hle : (p : ℝ) * b ≤ 1 := by
    calc (p : ℝ) * b ≤ (p : ℝ) * ((p : ℝ))⁻¹ :=
          mul_le_mul_of_nonneg_left hbs (le_of_lt hpR)
      _ = 1 := mul_inv_cancel₀ (ne_of_gt hpR)
  rw [max_eq_left hle, one_mul]

end FieldLayer

/-! ## §4 ★どちらの評価が良いかを測る -/

section Compare

/-- ★★**木の評価は本ファイルの評価以上に良い**。

`CyclicLayerDescent.norm_sum_iterate_sub_self_le_of_contract`（`:461`）が出すのは
`max (θ^{p−2}) c`、本ファイル §3 が出すのは（`b = θ` と読むと）`max θ c` の形。
`0 ≤ θ ≤ 1` かつ `3 ≤ p` なら `θ^{p−2} ≤ θ` なので木の方が小さい。

⇒ ★**使うべきは木の `CyclicLayerDescent` の方である。**
★本ファイル §1〜§3 は独立な 2 本目の証明として残す。 -/
theorem tree_contract_le_gauss {θ c : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) {p : ℕ} (hp : 3 ≤ p) :
    max (θ ^ (p - 2)) c ≤ max θ c := by
  refine max_le_max ?_ (le_refl c)
  calc θ ^ (p - 2) ≤ θ ^ 1 := pow_le_pow_of_le_one hθ0 hθ1 (by omega)
    _ = θ := pow_one θ

end Compare

/-! ## §5 使っている公理の一覧 -/

#print axioms norm_sum_iterate_sub_nsmul_le
#print axioms prime_dvd_sum_range
#print axioms two_no_gain
#print axioms norm_sum_iterate_sub_nsmul_le_prime
#print axioms norm_normalizedTrace_sub_le
#print axioms norm_normalizedTrace_sub_le_max_one
#print axioms norm_normalizedTrace_sub_le_of_small
#print axioms tree_contract_le_gauss

end TraceGainKernel

end ABC3.Found.PGC
