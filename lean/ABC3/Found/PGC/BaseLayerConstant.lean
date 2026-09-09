import ABC3.Found.PGC.FiniteExceptionsIcc
import ABC3.Found.PGC.WildDepthFieldDescent

/-!
# [pGC] 基底段の定数 `B` を**測った** —— 前波の見込み `p^{K₀}` は粗すぎた

## ★★訂正（自己訂正 14 度目）—— 私自身の前波の言葉

前波の報告で私はこう書いた（逐語）:

> 出口 `axLemma_of_eventually_Icc` の仮説のうち、**未測定は `B` ただ 1 つ**。
> ★`axWildDescent_pow` から `B = p^{K₀}` が取れそうだが、これは**まだ測っていない**。

★測った結果、**`B = p^{K₀}` は在庫より `p^{K₀−1}` 倍粗い**。
★正しい値は ★★**`B = p`**（`k` に依らない定数）である。

根拠（本波で読んだ）: `WildDepthFieldDescent.lean:141`

```
theorem axWildDescent_prime (K : PAdicLocalField p) :
    AxWildDescent K (fun _ => (p : ℝ))
```

そして同ファイル `:137` の docstring が自分で言っている:
「★`AxTowerDecay.axWildDescent_pow`(`c k = p^k`)より真に良い」。
★私は `axWildDescent_pow`（本体が名指しした方）だけを見て `p^{K₀}` と見込んでいた。
★★**本体の名指しは 20 波連続で外れ、実装者は毎回「その 1 つ手前の部品」で通す** ——
本波もその形である（名指し `axWildDescent_pow` → 実際に効いた部品 `axWildDescent_prime`）。

⇒ 出口の定数は `(p^{K₀})^{K₀} = p^{K₀²}` ではなく ★**`p^{K₀} · axConstant p`**。

## ★在庫の測定（コマンドを添える）

- `grep -rn "^theorem axWildDescent_mono" lean/ABC3/Found/PGC/*.lean`
  → ★**2 件**（`FirstJumpLedger.lean:458`（名前空間 `ABC3.Found.PGC`）と
  `PGroupDescentToAxWild.lean:202`（名前空間 `ABC3.Found.PGC.PGroupToAxWild`））。
  ★**同じ主張が木に 2 度書かれている。** 本ファイルは 3 度目を書かず、
  **真に強い形**（各 `k` ごとに別の `c'` を選んでよい＝貼り合わせ）だけを足す。
- `grep -n '^namespace' lean/ABC3/Found/PGC/WildDepthFieldDescent.lean` → `:65` の 1 行だけ
  （`ABC3.Found.PGC`）。★#359 を書く前に確認したので修飾なしで引けた。
- `hD`（`∏_{Icc 1 n} axDecay p k ≤ axConstant p`）は木に**独立した宣言としては無い** ——
  `AxEpsilonDecay.axLemma_of_axDecay`(`:532`)の**証明の中**に埋まっている。
  ★§4 で切り出した（`prod_axDecay_le`）。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `select_pointwise` | ★**抽象核**（`ι`・`m : ι → ℕ`・`P : ι → ℝ → Prop` だけ。体も群も出ない） |
| §2 | `wildStep` / `wildStep_mono` / `wildStep_iff` / `axWildDescent_glue` | 抽象核への代入 |
| §3 | `baseGlue` / `axWildDescent_baseGlue` | ★浅い段は `p`、深い段は `axDecay p` |
| §4 | `prod_axDecay_le` | `∏_{Icc 1 n} axDecay p k ≤ axConstant p`（木の証明から切り出し） |
| §5 | `axLemma_of_deep_bound` / `axSenTate_of_deep_bound` | ★★**出口** |
| §6 | `prev_guess_worse` | ★訂正の証拠 |

## ★残っているのはちょうど 1 点

`axSenTate_of_deep_bound` の仮説は **2 つだけ**:

* `hdeep : AxWildDescent K cdeep`
* `hbig : ∀ k, K₀ ≤ k → cdeep k ≤ axDecay p k`（★**深い段だけ**）

★`hc`（`1 ≤ c k`）も `hc1` も要らなくなった —— 貼り合わせ後の `c` は
浅い段で `p`、深い段で `axDecay p k` なので下界は自動である。
★前波の出口が要求していた 6 仮説のうち 4 つが消えた。

## 逸脱の記録

- `axWildDescent_mono` は木に既に 2 件あるので**再宣言しない**（#348）。
  本ファイルの `axWildDescent_glue` は `c' := c` と置けば `mono` を含む。
- §4 は `AxEpsilonDecay.axLemma_of_axDecay` の証明本体と**同じ 4 行**である。
  木がそれを宣言として持っていなかったので切り出した（内容の新規性は無い）。
-/

namespace ABC3.Found.PGC

namespace BaseLayerConstant

open Finset

/-! ## §1 ★抽象核 —— 各点ごとに別の証人を選んでよい

★★**分岐・付値・Galois の語が 1 語も出ない。** 出てくるのは
添字型 `ι`、深さ `m : ι → ℕ`、そして「第 2 引数について単調な述語」`P : ι → ℝ → Prop` だけ。 -/

section Kernel

/-- ★★**抽象核**。

`P i a` が `a` について単調なとき、「結論が `c` に **`c (m i)` の形でしか依存しない**」
主張は**各点ごとに別の `c'` を選んで貼り合わせられる**。

★これが「浅い段は定数 `p`、深い段は `axDecay p k`」を 1 本の `c` にまとめる仕組みである。
★木の `axWildDescent_mono`（2 箇所にある）は `c' := c` と置いた特別な場合。 -/
theorem select_pointwise {ι : Type*} (m : ι → ℕ) (P : ι → ℝ → Prop)
    (hmono : ∀ i : ι, ∀ a b : ℝ, a ≤ b → P i a → P i b)
    {c : ℕ → ℝ} (h : ∀ k : ℕ, ∃ c' : ℕ → ℝ, (∀ i : ι, P i (c' (m i))) ∧ c' k ≤ c k)
    (i : ι) : P i (c (m i)) := by
  obtain ⟨c', hc', hle⟩ := h (m i)
  exact hmono i _ _ hle (hc' i)

end Kernel

/-! ## §2 抽象核への代入 —— `AxWildDescent` の貼り合わせ -/

section Glue

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- `AxWildDescent K c` の「`ε` と `x` を 1 つ固定した 1 段分」。

★`AxWildDescent K c` は `∀ i, wildStep K i (c (wildDepth K i.2))` と**同値**
（`wildStep_iff`）で、右辺は §1 の抽象核がそのまま当たる形をしている。 -/
def wildStep (K : PAdicLocalField p) (i : ℝ × K.closure) (a : ℝ) : Prop :=
  0 ≤ i.1 → (∀ σ : K.absGal, ‖σ • i.2 - i.2‖ ≤ i.1) →
    p ∣ (minpoly K.carrier i.2).natDegree →
      ∃ x' : K.closure, wildDepth K x' < wildDepth K i.2 ∧
        ‖i.2 - x'‖ ≤ a * i.1 ∧ ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ a * i.1

/-- `wildStep` は定数について単調（★`0 ≤ ε` は `wildStep` の中から取れる）。 -/
theorem wildStep_mono (K : PAdicLocalField p) (i : ℝ × K.closure) (a b : ℝ) (hab : a ≤ b)
    (h : wildStep K i a) : wildStep K i b := by
  intro hε hx hdvd
  obtain ⟨x', hlt, h1, h2⟩ := h hε hx hdvd
  exact ⟨x', hlt, h1.trans (mul_le_mul_of_nonneg_right hab hε),
    fun σ => (h2 σ).trans (mul_le_mul_of_nonneg_right hab hε)⟩

/-- `AxWildDescent` を「1 段分の述語」の形に書き直す（★カリー化を外すだけ）。 -/
theorem wildStep_iff (K : PAdicLocalField p) (c : ℕ → ℝ) :
    AxWildDescent K c ↔ ∀ i : ℝ × K.closure, wildStep K i (c (wildDepth K i.2)) := by
  constructor
  · intro h i hε hx hdvd
    exact h i.1 hε i.2 hx hdvd
  · intro h ε hε x hx hdvd
    exact h (ε, x) hε hx hdvd

/-- ★★**`AxWildDescent` の貼り合わせ** —— 深さ `k` ごとに別の証人 `c'` を選んでよい。

★木の `axWildDescent_mono`（`FirstJumpLedger.lean:458` と
`PGroupDescentToAxWild.lean:202` の **2 箇所**にある）は `c' := c` と置いた場合。
★本補題は「浅い段は `axWildDescent_prime`、深い段は分岐の議論」を
**1 本の `c` にまとめる**ために要る。 -/
theorem axWildDescent_glue (K : PAdicLocalField p) {c : ℕ → ℝ}
    (h : ∀ k : ℕ, ∃ c' : ℕ → ℝ, AxWildDescent K c' ∧ c' k ≤ c k) :
    AxWildDescent K c :=
  (wildStep_iff K c).mpr fun i =>
    select_pointwise (fun j : ℝ × K.closure => wildDepth K j.2) (wildStep K)
      (wildStep_mono K)
      (fun k => by
        obtain ⟨c', hc', hle⟩ := h k
        exact ⟨c', fun j => (wildStep_iff K c').mp hc' j, hle⟩) i

end Glue

/-! ## §3 ★測った `B` —— 浅い段は `p`、深い段は `axDecay p` -/

section BaseGlue

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★**貼り合わせた定数列**: 深さ `< K₀` では `p`、それ以上では `axDecay p k`。

★浅い側の値が **`p`（`k` に依らない）** であることが本波の測定結果。
前波の見込み `p^k`（`axWildDescent_pow`）より `p^{k−1}` 倍良い。 -/
noncomputable def baseGlue (p : ℕ) (K₀ : ℕ) : ℕ → ℝ :=
  fun k => if k < K₀ then (p : ℝ) else axDecay p k

/-- `AxWildDescent K (baseGlue p K₀)` —— ★浅い側は `axWildDescent_prime`
（`WildDepthFieldDescent.lean:141`）、深い側は仮説 `hdeep` から。 -/
theorem axWildDescent_baseGlue (K : PAdicLocalField p) {cdeep : ℕ → ℝ} (K₀ : ℕ)
    (hbig : ∀ k, K₀ ≤ k → cdeep k ≤ axDecay p k)
    (hdeep : AxWildDescent K cdeep) :
    AxWildDescent K (baseGlue p K₀) := by
  refine axWildDescent_glue K fun k => ?_
  by_cases hk : k < K₀
  · exact ⟨fun _ => (p : ℝ), axWildDescent_prime K, by simp [baseGlue, hk]⟩
  · exact ⟨cdeep, hdeep, by simpa [baseGlue, hk] using hbig k (not_lt.mp hk)⟩

/-- `1 ≤ baseGlue p K₀ k`（★浅い側は `1 ≤ p`、深い側は `one_le_axDecay`）。 -/
theorem one_le_baseGlue (K₀ k : ℕ) : 1 ≤ baseGlue p K₀ k := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  by_cases hk : k < K₀
  · simpa [baseGlue, hk] using h1.le
  · simpa [baseGlue, hk] using one_le_axDecay p k

end BaseGlue

/-! ## §4 `axDecay` の有限積 —— 木の証明から切り出す

★`grep -rn "prod.*axDecay" lean/ABC3/Found/PGC/*.lean` では**宣言としては出ない**。
中身は `AxEpsilonDecay.axLemma_of_axDecay`(`:532`) の証明の 4 行に埋まっている。
★内容の新規性は無い（同じ 4 行）。★宣言になっていなかったので次のノードから引けなかった。 -/

section ProdBound

theorem prod_axDecay_le (p : ℕ) [Fact p.Prime] (n : ℕ) :
    ∏ k ∈ Finset.Icc 1 n, axDecay p k ≤ axConstant p := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (0 : ℝ) < (p : ℝ) := by linarith
  have hbound := prod_le_rpow_of_geometric_shift (D := (p : ℝ)) (A := 1 / ((p : ℝ) - 1))
    (r := 1 / (p : ℝ)) (le_of_lt h1) (by positivity) (by positivity)
    (by rw [div_lt_one hp0]; linarith)
    (fun k => le_trans zero_le_one (one_le_axDecay p k)) (fun k => le_of_eq rfl) n
  rwa [axExponent_eq h1] at hbound

end ProdBound

/-! ## §5 ★★出口 —— 仮説は 2 つだけ -/

section Exit

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**基底段を別扱いした出口**。

仮説は **2 つだけ**:

* `hdeep : AxWildDescent K cdeep`
* `hbig : ∀ k, K₀ ≤ k → cdeep k ≤ axDecay p k`（★**深さ `K₀` 以上でだけ**）

⇒ `AxLemma K (p^{K₀} · axConstant p)`。

★★定数の測定: 前波の私は `B = p^{K₀}` を見込んでいたので `(p^{K₀})^{K₀} = p^{K₀²}` に
なるはずだった。★在庫の `axWildDescent_prime` が `B = p` を与えるので **`p^{K₀}`** で済む。 -/
theorem axLemma_of_deep_bound (K : PAdicLocalField p) {cdeep : ℕ → ℝ} (K₀ : ℕ)
    (hbig : ∀ k, K₀ ≤ k → cdeep k ≤ axDecay p k)
    (hdeep : AxWildDescent K cdeep) :
    AxLemma K ((p : ℝ) ^ K₀ * axConstant p) := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  refine FiniteExceptionsIcc.axLemma_of_eventually_Icc K (d := axDecay p) (K₀ := K₀)
    (B := (p : ℝ)) (C := axConstant p)
    (fun k => one_le_baseGlue K₀ k) (one_le_axDecay p) h1.le
    (fun k hk => le_of_eq (by simp [baseGlue, hk]))
    (fun k hk => le_of_eq (by simp [baseGlue, Nat.not_lt.mpr hk]))
    (prod_axDecay_le p) (axWildDescent_baseGlue K K₀ hbig hdeep)

/-- ★★★**そこから Ax–Sen–Tate**。★これで `cor_3_1` の鎖が
「**深さ `K₀` 以上の 1 段が `axDecay p k` で済む**」ただ 1 点に落ちる
（★浅い `K₀` 段は在庫の `axWildDescent_prime` で無条件に埋まる）。 -/
theorem axSenTate_of_deep_bound (K : PAdicLocalField p) {cdeep : ℕ → ℝ} (K₀ : ℕ)
    (hbig : ∀ k, K₀ ≤ k → cdeep k ≤ axDecay p k)
    (hdeep : AxWildDescent K cdeep) : AxSenTate K := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  refine axSenTate_of_axLemma ?_ (axLemma_of_deep_bound K K₀ hbig hdeep)
  have := one_le_axConstant p
  positivity

end Exit

/-! ## §6 ★訂正の証拠 -/

section Evidence

variable {p : ℕ} [Fact p.Prime]

/-- ★前波の私の見込み `B = p^{K₀}` は、在庫の `B = p` より粗い。

★差は `p^{K₀−1}` 倍、出口の定数では `p^{K₀²}` 対 `p^{K₀}`。 -/
theorem prev_guess_worse (K₀ : ℕ) (hK : 1 ≤ K₀) : (p : ℝ) ≤ (p : ℝ) ^ K₀ := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  simpa using pow_le_pow_right₀ (le_of_lt h1) hK

end Evidence

/-! ## §7 使っている公理の一覧 -/

#print axioms select_pointwise
#print axioms wildStep_mono
#print axioms wildStep_iff
#print axioms axWildDescent_glue
#print axioms axWildDescent_baseGlue
#print axioms one_le_baseGlue
#print axioms prod_axDecay_le
#print axioms axLemma_of_deep_bound
#print axioms axSenTate_of_deep_bound
#print axioms prev_guess_worse

end BaseLayerConstant

end ABC3.Found.PGC
