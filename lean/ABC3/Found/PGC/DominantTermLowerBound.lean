import ABC3.Found.PGC.DigitLossFactor

/-!
# [pGC] ★★`p ∣ j₀` の逃げ道 —— 二項展開の**主役の項**による**下界**

## 持ち場（前波で「次の 1 点」とした点）

前波の私の言葉（逐語）:

> `p ∣ j₀` の場合の扱い。★測っていない具体的な問いは 1 つ: **`x` の主役の桁 `j₀` が `p ∣ j₀` のとき、
> `x` は「`π^p` で生成される部分体の元」にどれだけ近いか**。

## ★測った結果 —— 問いの立て方を変えた（自己訂正 20 度目）

★**「部分体に近いか」は測らなかった。** 代わりに測ったのは
★**「`p ∣ j₀` でも `‖u^{j₀} − π^{j₀}‖` の下界が別の項から出るか」**である。理由:

`MaxMinIndex.lean:143`（★本体が `grep` で先に測ってくれた）:

```lean
theorem not_dvd_max_min_index {p j₀ : ℕ} (hp : 2 ≤ p) (h1 : 1 ≤ j₀) (h2 : j₀ ≤ p - 1) :
    ¬ p ∣ j₀
```

⇒ ★旧設定（`j₀ ∈ [1, p−1]`）では `p ∤ j₀` は**自動**だった。
★新設定（次数 `p^{k+1}`）では `j₀ ∈ [1, p^{k+1}−1]` なので `p ∣ j₀` が起こりうる。

★★**前波の `DigitLossFactor` は「上界がどれだけ落ちるか」を測った。本波はその裏、
「**下界がどこから出るか**」である。** 前波は `l = 0`（形式微分）の項だけを主役と見ていたが、
★二項展開には `j` 個の項があり、★**`p ∣ j` のときは別の項が主役になりうる**。

## ★★★本ファイルの中身

`u = π + δ`（`δ = u − π`）として `add_pow` から

  `(π + δ)^j − π^j = Σ_{k ∈ range j} π^k · δ^{j−k} · C(j,k)`   （§1、★恒等式）

* `k = j−1` の項 = `π^{j−1}·δ·j` ——★**形式微分**（前波が見ていた項）
* `k = 0` の項 = `δ^j` ——★★**係数が `C(j,0) = 1` で必ず単数**

★★超距離では「**1 つの項が他をすべて真に上回れば和のノルムはその項に等しい**」（§2、
★これは**下界**を与える）。⇒ §3 で `k = 0` の項を主役にすると

  `‖u^j − π^j‖ = ‖δ‖^j`

★★★**`j = p^a` のとき `hdom` は機械的に確かめられる**（§4）:
`Nat.Prime.dvd_choose_pow` により `0 < k < p^a` では `p ∣ C(p^a, k)` なので
`‖C(p^a,k)‖ ≤ ‖(p:L)‖ < 1`。

⇒ ★**`x = π^{p^a}`（私が第 1135 で「下界が出ない」と書いた正準例）では、
下界が `‖δ‖^{p^a}` として出る。**

## ★在庫の測定（予防 4 本目を書く前に実行）

* `grep -rn "norm_sum_eq_of\|map_sum_eq_of_lt\|主役" lean/ABC3/Found/PGC/*.lean`
  → 木にあるのは **付値版** `map_sum_eq_of_lt`（`UniformizerExpansion.lean:74`）で、
  ★docstring が「mathlib に **`AddValuation` 版が無い**（乗法版 `Valuation.map_sum_eq_of_lt` はある）」
  と言っている。★**ノルム版（`‖·‖`）は木にも無い**。§2 がそれである。
* `IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne` は mathlib に在る
  （`CyclicJumpNorm.lean:274` が使っている）が、★**「相異なる」より「1 つが勝つ」の方が弱くて使いやすい**。
* `Nat.Prime.dvd_choose_pow`（`Data/Nat/Multiplicity.lean:260`）を `.cache/mathlib-index.txt` で確認。

## ★残っているちょうど 1 点

`j = p^a·m`（`p ∤ m`, `m ≥ 2`）の場合。★Kummer の定理により `p ∤ C(j, p^a)` なので
`k = j − p^a` の項が主役の候補だが、★**mathlib の Kummer は `C(p^n, k)` の形しか無い**
（`emultiplicity_choose_prime_pow`）。★一般の `C(j, p^a)` は**測っていない**。

## 逸脱の記録

- §1–§3 は超距離ノルム体だけで、★分岐・付値・Galois・桁展開の語彙が 1 語も出ない。
- ★§4 は `hdom` を**仮定として残していない**（`Nat.Prime.dvd_choose_pow` で消した）が、
  「`‖δ‖^j` が他を上回る」ことは `‖δ‖ < ‖π‖` からは出ない ——
  ★**`‖π‖^k‖δ‖^{j−k}‖C‖ < ‖δ‖^j` は `k ≥ 1` で `‖π‖^k` が効くので逆向き**である。
  ⇒ §4 は `‖π‖ ≤ ‖δ‖`（★`δ` が大きい側）を仮定する。★これは `σ` を**深く取る**ことに当たる。
  ★どの `σ` でそれが実現するかは本ファイルでは**測っていない**。
-/

namespace ABC3.Found.PGC

namespace DominantTermLowerBound

/-! ## §1 ★恒等式 —— 二項展開 -/

section Identity

variable {L : Type*} [CommRing L]

/-- ★`(π + δ)^j − π^j = Σ_{k ∈ range j} π^k · δ^{j−k} · C(j,k)`。

★`k = j−1` の項が**形式微分** `j π^{j−1} δ`、★`k = 0` の項が `δ^j`（係数は `C(j,0) = 1`）。 -/
theorem pow_sub_pow_eq_sum_choose (π δ : L) (j : ℕ) :
    (π + δ) ^ j - π ^ j
      = ∑ k ∈ Finset.range j, π ^ k * δ ^ (j - k) * ((j.choose k : ℕ) : L) := by
  have h := add_pow π δ j
  rw [Finset.sum_range_succ] at h
  simp only [Nat.sub_self, pow_zero, Nat.choose_self, Nat.cast_one, mul_one] at h
  rw [h]
  ring

end Identity

/-! ## §2 ★★抽象核 —— 超距離で「主役の項」が和を決める（★下界） -/

section Dominant

variable {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]

/-- ★★★**主役の項が和のノルムを決める**（超距離）。

`i₀ ∈ s` の項が他のすべてを**真に**上回れば `‖Σ_{i ∈ s} f i‖ = ‖f i₀‖`。

★★これは**下界**を与える補題である（上界だけなら `norm_sum_le_of_forall_le` で足りる）。
★木にあるのは**付値版** `UniformizerExpansion.map_sum_eq_of_lt` だけで、
★**ノルム版は木にも mathlib にも見当たらない**（モジュール docstring に測定を書いた）。
★mathlib の `IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne` は
「相異なる」を要求するので、こちらの方が仮説が弱い。 -/
theorem norm_sum_eq_of_dominant {ι : Type*} (s : Finset ι) (f : ι → M) {i₀ : ι}
    (hi₀ : i₀ ∈ s) (hdom : ∀ i ∈ s, i ≠ i₀ → ‖f i‖ < ‖f i₀‖) :
    ‖∑ i ∈ s, f i‖ = ‖f i₀‖ := by
  classical
  rcases Finset.eq_empty_or_nonempty (s.erase i₀) with he | hne
  · have hsum : ∑ i ∈ s, f i = f i₀ := by
      rw [← Finset.add_sum_erase s f hi₀, he, Finset.sum_empty, add_zero]
    rw [hsum]
  · obtain ⟨i₁, hi₁, hmax⟩ := Finset.exists_max_image (s.erase i₀) (fun i => ‖f i‖) hne
    have hlt : ‖∑ i ∈ s.erase i₀, f i‖ < ‖f i₀‖ := by
      refine lt_of_le_of_lt
        (IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (norm_nonneg (f i₁)) ?_) ?_
      · intro i hi
        exact hmax i hi
      · exact hdom i₁ (Finset.mem_of_mem_erase hi₁) (Finset.ne_of_mem_erase hi₁)
    rw [← Finset.add_sum_erase s f hi₀,
      IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (ne_of_lt hlt).symm]
    exact max_eq_left hlt.le

end Dominant

/-! ## §3 ★二項展開に代入 -/

section Apply

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- ★★`k₀` の項が主役なら `‖(π+δ)^j − π^j‖` はその項のノルムに**等しい**。 -/
theorem norm_pow_sub_pow_eq_of_dominant (π δ : L) {j k₀ : ℕ} (hk₀ : k₀ ∈ Finset.range j)
    (hdom : ∀ k ∈ Finset.range j, k ≠ k₀ →
      ‖π ^ k * δ ^ (j - k) * ((j.choose k : ℕ) : L)‖
        < ‖π ^ k₀ * δ ^ (j - k₀) * ((j.choose k₀ : ℕ) : L)‖) :
    ‖(π + δ) ^ j - π ^ j‖ = ‖π ^ k₀ * δ ^ (j - k₀) * ((j.choose k₀ : ℕ) : L)‖ := by
  rw [pow_sub_pow_eq_sum_choose π δ j]
  exact norm_sum_eq_of_dominant _ _ hk₀ hdom

/-- ★★★**`k₀ = 0` の項（`δ^j`、係数 `C(j,0) = 1`）が主役のときの下界**。

  `‖(π+δ)^j − π^j‖ = ‖δ‖^j`

★★これが `p ∣ j` でも消えない項である（`C(j,0) = 1` は常に単数）。 -/
theorem norm_pow_sub_pow_eq_pow_of_dominant (π δ : L) {j : ℕ} (hj : 0 < j)
    (hdom : ∀ k ∈ Finset.range j, k ≠ 0 →
      ‖π ^ k * δ ^ (j - k) * ((j.choose k : ℕ) : L)‖ < ‖δ‖ ^ j) :
    ‖(π + δ) ^ j - π ^ j‖ = ‖δ‖ ^ j := by
  have hmem : (0 : ℕ) ∈ Finset.range j := Finset.mem_range.mpr hj
  have hzero : ‖π ^ (0 : ℕ) * δ ^ (j - 0) * ((j.choose 0 : ℕ) : L)‖ = ‖δ‖ ^ j := by
    simp [norm_pow]
  rw [norm_pow_sub_pow_eq_of_dominant π δ hmem (by simpa [hzero] using hdom), hzero]

end Apply

/-! ## §4 ★★`j = p^a` —— `hdom` が機械的に確かめられる場合 -/

section PrimePow

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- ★`0 < k < p^a` なら `‖C(p^a, k)‖ < 1`（`Nat.Prime.dvd_choose_pow` ＋
`GainedBridge.norm_natCast_le_of_dvd`）。 -/
theorem norm_choose_prime_pow_lt_one {p a k : ℕ} (hp : p.Prime) (hk : k ≠ 0)
    (hkp : k ≠ p ^ a) (hple : ‖(p : L)‖ < 1) :
    ‖(((p ^ a).choose k : ℕ) : L)‖ < 1 :=
  lt_of_le_of_lt (GainedBridge.norm_natCast_le_of_dvd (hp.dvd_choose_pow hk hkp)) hple

/-- ★★★**正準例 `j = p^a` での下界**。

`‖π‖ ≤ ‖δ‖`（★`δ` が大きい側 ＝ `σ` を深く取る）かつ `0 < ‖δ‖` なら

  `‖(π+δ)^{p^a} − π^{p^a}‖ = ‖δ‖^{p^a}`

★★これが第 1135 で「下界が出ない」と書いた `x = π^p` の場合の**答え**である。
★`‖π‖ ≤ ‖δ‖` がどの `σ` で実現するかは本ファイルでは**測っていない**。 -/
theorem norm_pow_sub_pow_prime_pow (π δ : L) {p a : ℕ} (hp : p.Prime)
    (hple : ‖(p : L)‖ < 1) (hδ0 : 0 < ‖δ‖) (hπδ : ‖π‖ ≤ ‖δ‖) :
    ‖(π + δ) ^ (p ^ a) - π ^ (p ^ a)‖ = ‖δ‖ ^ (p ^ a) := by
  have hjpos : 0 < p ^ a := pow_pos hp.pos a
  refine norm_pow_sub_pow_eq_pow_of_dominant π δ hjpos ?_
  intro k hk hk0
  have hklt : k < p ^ a := Finset.mem_range.mp hk
  have hcoef : ‖(((p ^ a).choose k : ℕ) : L)‖ < 1 :=
    norm_choose_prime_pow_lt_one hp hk0 (by omega) hple
  have hπk : ‖π‖ ^ k ≤ ‖δ‖ ^ k := pow_le_pow_left₀ (norm_nonneg _) hπδ k
  have hsplit : ‖δ‖ ^ k * ‖δ‖ ^ (p ^ a - k) = ‖δ‖ ^ (p ^ a) := by
    rw [← pow_add]
    congr 1
    omega
  calc ‖π ^ k * δ ^ (p ^ a - k) * (((p ^ a).choose k : ℕ) : L)‖
      = ‖π‖ ^ k * ‖δ‖ ^ (p ^ a - k) * ‖(((p ^ a).choose k : ℕ) : L)‖ := by
        rw [norm_mul, norm_mul, norm_pow, norm_pow]
    _ ≤ ‖δ‖ ^ k * ‖δ‖ ^ (p ^ a - k) * ‖(((p ^ a).choose k : ℕ) : L)‖ := by
        have h1 : (0 : ℝ) ≤ ‖δ‖ ^ (p ^ a - k) * ‖(((p ^ a).choose k : ℕ) : L)‖ := by positivity
        nlinarith [hπk, h1]
    _ = ‖δ‖ ^ (p ^ a) * ‖(((p ^ a).choose k : ℕ) : L)‖ := by rw [hsplit]
    _ < ‖δ‖ ^ (p ^ a) * 1 := by
        have hpos : (0 : ℝ) < ‖δ‖ ^ (p ^ a) := by positivity
        exact mul_lt_mul_of_pos_left hcoef hpos
    _ = ‖δ‖ ^ (p ^ a) := mul_one _

end PrimePow

/-! ## §5 使っている公理の一覧 -/

#print axioms pow_sub_pow_eq_sum_choose
#print axioms norm_sum_eq_of_dominant
#print axioms norm_pow_sub_pow_eq_of_dominant
#print axioms norm_pow_sub_pow_eq_pow_of_dominant
#print axioms norm_choose_prime_pow_lt_one
#print axioms norm_pow_sub_pow_prime_pow

end DominantTermLowerBound

end ABC3.Found.PGC
