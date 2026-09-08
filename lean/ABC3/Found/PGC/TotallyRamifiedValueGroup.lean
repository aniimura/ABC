import ABC3.Found.PGC.GainedTowerModel

/-!
# [pGC] `hvalj`(層ごとの全分岐 `k+1` 本)を **底の 1 本** `hvalK` に落とす

`Found/PGC/GainedTowerModel.lean` の出口
`exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic` に残っていた唯一の仮説は

    hvalj : ∀ j ≤ k, ∀ d ∈ (E j)^×, ∃ m, ‖d‖ = ‖π‖^(p^{k+1−j}·m)

すなわち「**構成された固定体 `E j = M^{⟨g^{p^j}⟩}` のすべてで全分岐**」であった。
本ファイルはこれを

    hvalK : ∀ a ∈ K^×, ∃ m, ‖a‖ = ‖π‖^(p^{k+1}·m)      (★底 `K` だけ)
    hnK   : [M : K] = p^{k+1}

の 2 本に落とす(`valj_of_valK` / 出口は
`exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_totallyRamified`)。

## ★★前波の見立てのどこが違ったか

前波は「落とすのに 2 本要る」と書いた:

1. `d_j ≤ [M : E_j]`(1 層で済む)
2. `d_j / d_{j+1} ≤ [E_{j+1} : E_j]`(★**2 層をまたぐ ＝ #59 の危険帯**)

★**2 は要らなかった。** 層どうし(`E_j` と `E_{j+1}`)を比べる代わりに、
**`E_j` を底 `K` と直接比べる**と同じ結論が出る:

    d_j · [Γ_{E_j} : Γ_K] = [Γ_M : Γ_K] = p^{k+1},   [Γ_{E_j} : Γ_K] ≤ [E_j : K] = p^j
    ⟹ d_j ≥ p^{k+1−j}

★これは #314 の構え(**層どうしを比較せず全部底 `K` から測る**)そのもので、
本ファイルは `IntermediateField` を使いながら #59/#69 に**一度も当たっていない**。

## 構成(★抽象核 → 具体層の順)

* §1 `linearIndependent_of_norm_pairwise_notMem_coset` —— ★**抽象核**。
  超距離ノルム体で「ノルムが `‖K^×‖` の**相異なる剰余類**にある」族は 1 次独立。
  ★分岐・付値・Galois の語彙が 1 語も出ない。
* §2 `linearIndependent_of_ne_mod` / `card_le_finrank_of_ne_mod` —— 指数版。
  `exists_zpow_norm` —— ★`Γ_M ⊆ ‖π‖^ℤ`(`z, π^0, …, π^{n−1}` の `n+1` 個で矛盾)。
* §3 `int_subgroup_dvd` —— ★**抽象核(ℤ だけ)**。`ℤ` の部分群は `cℤ`。
* §4 `exists_zpow_norm_intermediate` —— ★★**主定理**。挟み撃ち
  `c ≤ [M:E]`(`π^l` を `E` 上で)と `gcd(n,c) ≥ [M:E]`(`d₀^i` を `K` 上で、★`E` の中で数える)。
* §5 塔への代入(`twr g p j` は既存の `GainedTowerModel.GaloisTower`)。
* §6 ★検算: `hvalK` が**不分岐を排除する**こと(`valK_forces_ramified`)と、数値側の非空虚性。

## ★在庫の測定(コマンドを残す)

```
grep -n "Int.subgroup_cyclic|mem_closure_singleton" .cache/mathlib-index.txt
  → Int.subgroup_cyclic は在る。★AddSubgroup.mem_closure_singleton は
    索引に**出ない**(to_additive 生成名)が**在る**。索引の「無い」は嘘(既知の 7 例に 1 つ追加)。
grep -n "fintype_card_le_finrank" .cache/mathlib-index.txt
  → LinearIndependent.fintype_card_le_finrank が在る。
grep -n "zpow_right_inj" .cache/mathlib-index.txt
  → zpow_right_inj₀ (0 < a) (a ≠ 1) : a ^ m = a ^ n ↔ m = n。★これが指数の一意性の要。
grep -rn "IsTotallyRamified" lean/ABC3/ --include=*.lean ; grep -in "totallyramified" .cache/mathlib-index.txt
  → mathlib に 0 件(前波の測定を追認)。★**自前の述語は作らなかった** ——
    `∀ a ≠ 0, ∃ m, ‖a‖ = ‖π‖^(n·m)` という**素の命題**のままで 4 本の定理が通る。
```

## ★★閉じていないもの

`hvalK` 自身は**落ちない**(落としてはいけない)。§6 `valK_forces_ramified` が示す通り、
`hvalK` は「`M/K` が全分岐で `π` が素元」そのものであり、不分岐拡大では偽になる。
⇒ 残りは**点 1(体・ノルム込みの模型 `ℚ₃` の 18 次全分岐拡大の構成)**だけである。
-/

namespace ABC3.Found.PGC

namespace TotallyRamified

open Finset

section Core

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★**抽象核 1** —— ノルムが `‖K^×‖` の**相異なる剰余類**にある族は `K` 上 1 次独立。 -/
theorem linearIndependent_of_norm_pairwise_notMem_coset
    {ι : Type*} [Fintype ι] (x : ι → M) (hx : ∀ i, x i ≠ 0)
    (h : ∀ i j, i ≠ j → ∀ a : K, a ≠ 0 → ‖x i‖ ≠ ‖algebraMap K M a‖ * ‖x j‖) :
    LinearIndependent K x := by
  classical
  rw [Fintype.linearIndependent_iff]
  intro g hg i₀
  by_contra hgi
  set t : Finset ι := Finset.univ.filter (fun i => g i ≠ 0) with ht
  have hi₀t : i₀ ∈ t := by simp [ht, hgi]
  have hsum : ∑ i ∈ t, g i • x i = 0 := by
    rw [Finset.sum_subset (f := fun i => g i • x i) (Finset.subset_univ t)
      (fun i _ hit => by
        have hz : g i = 0 := by simpa [ht] using hit
        simp [hz])]
    exact hg
  have hpair : (t : Set ι).Pairwise fun i j => ‖g i • x i‖₊ ≠ ‖g j • x j‖₊ := by
    intro i hi j hj hij hcon
    have hgi' : g i ≠ 0 := by simpa [ht] using hi
    have hgj' : g j ≠ 0 := by simpa [ht] using hj
    have hreal : ‖g i • x i‖ = ‖g j • x j‖ := by
      simpa using congrArg NNReal.toReal hcon
    rw [Algebra.smul_def, Algebra.smul_def, norm_mul, norm_mul] at hreal
    have hne0 : ‖algebraMap K M (g i)‖ ≠ 0 := by
      simpa using (map_ne_zero (algebraMap K M)).mpr hgi'
    refine h i j hij (g j / g i) (div_ne_zero hgj' hgi') ?_
    rw [map_div₀, norm_div]
    field_simp
    linarith [hreal]
  have hzero : ‖∑ i ∈ t, g i • x i‖₊ = t.sup fun i => ‖g i • x i‖₊ :=
    IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne hpair
  rw [hsum] at hzero
  have hle : ‖g i₀ • x i₀‖₊ ≤ 0 := by
    have hs := Finset.le_sup (f := fun i => ‖g i • x i‖₊) hi₀t
    rw [← hzero] at hs
    simpa using hs
  have hz : g i₀ • x i₀ = 0 := by
    have : ‖g i₀ • x i₀‖₊ = 0 := le_antisymm hle (by simp)
    simpa using this
  rw [Algebra.smul_def] at hz
  rcases mul_eq_zero.mp hz with h1 | h2
  · exact hgi ((map_eq_zero (algebraMap K M)).mp h1)
  · exact hx i₀ h2

end Core
section Exponent

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★**抽象核 2** —— ノルムが `β^{w i}` で、指数 `w i` が `n` を法として相異なる族は
`K` 上 1 次独立。★`β` は `‖π‖` のつもりだが、ここでは体の元でなくてよい。 -/
theorem linearIndependent_of_ne_mod
    {β : ℝ} (hβ0 : 0 < β) (hβ1 : β ≠ 1) {n : ℤ}
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = β ^ (n * m))
    {ι : Type*} [Fintype ι] (x : ι → M) (w : ι → ℤ)
    (hxw : ∀ i, ‖x i‖ = β ^ w i)
    (hne : ∀ i j, i ≠ j → ¬ (n ∣ (w i - w j))) :
    LinearIndependent K x := by
  refine linearIndependent_of_norm_pairwise_notMem_coset x (fun i => ?_) (fun i j hij a ha => ?_)
  · have hpos : 0 < ‖x i‖ := by rw [hxw i]; exact zpow_pos hβ0 _
    exact norm_pos_iff.mp hpos
  · intro hcon
    obtain ⟨m, hm⟩ := hvalK a ha
    rw [hxw i, hxw j, hm, ← zpow_add₀ (ne_of_gt hβ0)] at hcon
    have hw := (zpow_right_inj₀ hβ0 hβ1).mp hcon
    exact hne i j hij ⟨m, by rw [hw]; ring⟩

/-- 抽象核 2 の系 —— 族の個数は `[M : K]` 以下。 -/
theorem card_le_finrank_of_ne_mod [FiniteDimensional K M]
    {β : ℝ} (hβ0 : 0 < β) (hβ1 : β ≠ 1) {n : ℤ}
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = β ^ (n * m))
    {ι : Type*} [Fintype ι] (x : ι → M) (w : ι → ℤ)
    (hxw : ∀ i, ‖x i‖ = β ^ w i)
    (hne : ∀ i j, i ≠ j → ¬ (n ∣ (w i - w j))) :
    Fintype.card ι ≤ Module.finrank K M :=
  LinearIndependent.fintype_card_le_finrank
    (linearIndependent_of_ne_mod hβ0 hβ1 hvalK x w hxw hne)
/-- ★★★**底が全分岐なら `M` の値群は `‖π‖^ℤ` に入る**(整数性)。

`n = [M : K]` 本の `π^l`(`l < n`)と `z` の合計 `n+1` 個は、`‖z‖` が `‖π‖` の
整数冪でなければ `K` 上 1 次独立になり、`n+1 ≤ n` で矛盾。 -/
theorem exists_zpow_norm [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {z : M} (hz : z ≠ 0) :
    ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m := by
  classical
  by_contra hcon
  simp only [not_exists] at hcon
  set x : Option (Fin n) → M := fun i => i.elim z (fun l => π ^ (l : ℕ)) with hxdef
  have hxnorm : ∀ l : Fin n, ‖x (some l)‖ = ‖π‖ ^ ((l : ℕ) : ℤ) := by
    intro l; simp [hxdef, zpow_natCast]
  have hx0 : ∀ i, x i ≠ 0 := by
    intro i
    cases i with
    | none => exact hz
    | some l =>
        have : (0:ℝ) < ‖x (some l)‖ := by rw [hxnorm l]; exact zpow_pos hπ0 _
        exact norm_pos_iff.mp this
  have hkey : ∀ i j : Option (Fin n), i ≠ j → ∀ a : K, a ≠ 0 →
      ‖x i‖ ≠ ‖algebraMap K M a‖ * ‖x j‖ := by
    intro i j hij a ha hEq
    obtain ⟨m, hm⟩ := hvalK a ha
    rw [hm] at hEq
    cases i with
    | none =>
        cases j with
        | none => exact hij rfl
        | some l =>
            rw [hxnorm l, ← zpow_add₀ (ne_of_gt hπ0)] at hEq
            exact hcon _ hEq
    | some l =>
        cases j with
        | none =>
            rw [hxnorm l] at hEq
            refine hcon (((l : ℕ) : ℤ) - (n : ℤ) * m) ?_
            have hne0 : ‖π‖ ^ ((n : ℤ) * m) ≠ 0 := (zpow_pos hπ0 _).ne'
            rw [zpow_sub₀ (ne_of_gt hπ0), eq_div_iff hne0, mul_comm]
            exact hEq.symm
        | some l' =>
            rw [hxnorm l, hxnorm l', ← zpow_add₀ (ne_of_gt hπ0)] at hEq
            have hw := (zpow_right_inj₀ hπ0 hπ1).mp hEq
            have hdvd : (n : ℤ) ∣ ((l : ℕ) : ℤ) - ((l' : ℕ) : ℤ) := ⟨m, by rw [hw]; ring⟩
            have hlt : ((l : ℕ) : ℤ) - ((l' : ℕ) : ℤ) = 0 := by
              have h1 : ((l : ℕ) : ℤ) < (n : ℤ) := by exact_mod_cast l.isLt
              have h2 : ((l' : ℕ) : ℤ) < (n : ℤ) := by exact_mod_cast l'.isLt
              have h3 : (0:ℤ) ≤ ((l : ℕ) : ℤ) := Int.natCast_nonneg _
              have h4 : (0:ℤ) ≤ ((l' : ℕ) : ℤ) := Int.natCast_nonneg _
              exact Int.eq_zero_of_abs_lt_dvd hdvd (abs_lt.mpr ⟨by omega, by omega⟩)
            exact hij (by
              have : (l : ℕ) = (l' : ℕ) := by omega
              simp [Fin.ext_iff, this])
  have hcard := LinearIndependent.fintype_card_le_finrank
    (linearIndependent_of_norm_pairwise_notMem_coset (K := K) x hx0 hkey)
  rw [hn] at hcard
  simp at hcard

end Exponent

/-! ## §3 値群は `ℤ` の部分群 -/

section ValueSub

/-- ★**抽象核(ℤ だけ)** —— `ℤ` の部分群は「`c` の倍数全体」。 -/
theorem int_subgroup_dvd (H : AddSubgroup ℤ) : ∃ c : ℕ, ∀ m : ℤ, m ∈ H ↔ (c : ℤ) ∣ m := by
  obtain ⟨a, ha⟩ := Int.subgroup_cyclic H
  refine ⟨a.natAbs, fun m => ?_⟩
  have hmem : m ∈ H ↔ ∃ k : ℤ, k • a = m := by
    rw [ha]; exact AddSubgroup.mem_closure_singleton
  rw [hmem, Int.natAbs_dvd]
  constructor
  · rintro ⟨k, rfl⟩; exact ⟨k, by rw [smul_eq_mul]; ring⟩
  · rintro ⟨k, rfl⟩; exact ⟨k, by rw [smul_eq_mul]; ring⟩

variable {E M : Type*} [Field E] [NormedField M] [Algebra E M]

/-- `E ⊆ M` の値群の指数集合。★`‖π‖ > 0` だけを使う(超距離も次数も要らない)。 -/
def valSub (E : Type*) [Field E] {M : Type*} [NormedField M] [Algebra E M] {π : M}
    (hπ0 : 0 < ‖π‖) : AddSubgroup ℤ where
  carrier := {m : ℤ | ∃ d : E, d ≠ 0 ∧ ‖algebraMap E M d‖ = ‖π‖ ^ m}
  zero_mem' := ⟨1, one_ne_zero, by simp⟩
  add_mem' := by
    rintro a b ⟨d, hd, hda⟩ ⟨d', hd', hdb⟩
    exact ⟨d * d', mul_ne_zero hd hd', by
      rw [map_mul, norm_mul, hda, hdb, zpow_add₀ (ne_of_gt hπ0)]⟩
  neg_mem' := by
    rintro a ⟨d, hd, hda⟩
    exact ⟨d⁻¹, inv_ne_zero hd, by rw [map_inv₀, norm_inv, hda, ← zpow_neg]⟩

/-- ★値群の生成元 `c`。★`c = 0` は「`E` のノルムが自明」の場合。 -/
theorem exists_valSub_gen {π : M} (hπ0 : 0 < ‖π‖) :
    ∃ c : ℕ, ∀ m : ℤ,
      (∃ d : E, d ≠ 0 ∧ ‖algebraMap E M d‖ = ‖π‖ ^ m) ↔ (c : ℤ) ∣ m := by
  obtain ⟨c, hc⟩ := int_subgroup_dvd (valSub E hπ0)
  exact ⟨c, fun m => hc m⟩

end ValueSub

/-! ## §4 主定理 —— 中間層の値群 -/

section Main

variable {K E M : Type*} [Field K] [Field E] [NormedField M] [IsUltrametricDist M]
  [Algebra K M] [Algebra K E] [Algebra E M] [IsScalarTower K E M]

/-- ★★★★**主定理** —— 底 `K` の値群が `‖π‖^{n·ℤ}`(`n = [M:K]`)に入るなら、
**どの中間層 `E` の値群も** `‖π‖^{[M:E]·ℤ}` に入る。

★★これが `hvalj`(層ごとの全分岐)を**底の 1 本**に落とす中身である。
証明は 3 本の 1 次独立(すべて「1 層ぶん」)の挟み撃ち:

* `M` の値群は `‖π‖^ℤ`(`exists_zpow_norm`、`n+1` 個の族)
* `c := ` `E` の値群の生成元 について `c ≤ [M:E]`(`π^l`, `l < c` を `E` 上で)
* `gcd(n,c) ≥ [M:E]`(`d₀^i`, `i < n/gcd` を `K` 上で。★`E` の中で数える) -/
theorem exists_zpow_norm_intermediate
    [FiniteDimensional K M] [FiniteDimensional E M] [FiniteDimensional K E]
    {π : M} {n q : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hn : Module.finrank K M = n) (hq : Module.finrank E M = q)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (d : E) (hd : d ≠ 0) :
    ∃ m : ℤ, ‖algebraMap E M d‖ = ‖π‖ ^ ((q : ℤ) * m) := by
  classical
  set r := Module.finrank K E with hrdef
  have hrq : r * q = n := by rw [hrdef, ← hq, ← hn]; exact Module.finrank_mul_finrank K E M
  have hn0 : 0 < n := by rw [← hn]; exact Module.finrank_pos
  have hq0 : 0 < q := by rw [← hq]; exact Module.finrank_pos
  obtain ⟨c, hc⟩ := exists_valSub_gen (E := E) (M := M) hπ0
  have hint : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m :=
    fun z hz => exists_zpow_norm (K := K) hπ0 hπ1 hn hvalK hz
  have hvalE : ∀ e : E, e ≠ 0 → ∃ m : ℤ, ‖algebraMap E M e‖ = ‖π‖ ^ ((c : ℤ) * m) := by
    intro e he
    obtain ⟨m, hm⟩ := hint (algebraMap E M e) ((map_ne_zero (algebraMap E M)).mpr he)
    obtain ⟨m', rfl⟩ := (hc m).mp ⟨e, he, hm⟩
    exact ⟨m', hm⟩
  rcases Nat.eq_zero_or_pos c with hc0 | hcpos
  · obtain ⟨m, hm⟩ := hvalE d hd
    exact ⟨0, by rw [hm, hc0]; simp⟩
  -- (b) `c ≤ q`
  have hb : c ≤ q := by
    have hcard := card_le_finrank_of_ne_mod (K := E) (M := M) hπ0 hπ1 (n := (c : ℤ)) hvalE
      (ι := Fin c) (fun l => π ^ (l : ℕ)) (fun l => ((l : ℕ) : ℤ))
      (fun l => by simp [zpow_natCast]) (fun i j hij hdvd => ?_)
    · simpa [hq] using hcard
    · refine hij (Fin.ext ?_)
      have h1 : ((i : ℕ) : ℤ) - ((j : ℕ) : ℤ) = 0 := by
        refine Int.eq_zero_of_abs_lt_dvd hdvd ?_
        have hi := i.isLt; have hj := j.isLt
        rw [abs_lt]; omega
      omega
  -- (c) `q ≤ gcd n c`
  obtain ⟨d₀, hd₀, hd₀n⟩ := (hc (c : ℤ)).mpr dvd_rfl
  set G := Nat.gcd n c with hGdef
  have hG0 : 0 < G := Nat.gcd_pos_of_pos_left c hn0
  set t := n / G with htdef
  set c' := c / G with hc'def
  have hnt : G * t = n := Nat.mul_div_cancel' (Nat.gcd_dvd_left n c)
  have hct : G * c' = c := Nat.mul_div_cancel' (Nat.gcd_dvd_right n c)
  have hcop : Nat.Coprime t c' := Nat.coprime_div_gcd_div_gcd hG0
  have ht0 : 0 < t := by
    rcases Nat.eq_zero_or_pos t with h | h
    · rw [h, Nat.mul_zero] at hnt; omega
    · exact h
  have hLI : LinearIndependent K (((IsScalarTower.toAlgHom K E M).toLinearMap : E →ₗ[K] M) ∘
      (fun i : Fin t => d₀ ^ (i : ℕ))) := by
    refine linearIndependent_of_ne_mod (K := K) hπ0 hπ1 hvalK _
      (fun i => (c : ℤ) * ((i : ℕ) : ℤ)) (fun i => ?_) (fun i j hij hdvd => ?_)
    · simp only [Function.comp_apply, AlgHom.toLinearMap_apply, IsScalarTower.coe_toAlgHom']
      rw [map_pow, norm_pow, hd₀n, ← zpow_natCast (‖π‖ ^ (c : ℤ)) (i : ℕ), ← zpow_mul]
    · -- `n ∣ c * (i − j)` から `i = j`
      have hfac : ((c : ℤ) * ((i : ℕ) : ℤ) - (c : ℤ) * ((j : ℕ) : ℤ))
          = (G : ℤ) * ((c' : ℤ) * (((i : ℕ) : ℤ) - ((j : ℕ) : ℤ))) := by
        have : ((c : ℤ)) = (G : ℤ) * (c' : ℤ) := by exact_mod_cast hct.symm
        rw [this]; ring
      have hnfac : ((n : ℤ)) = (G : ℤ) * (t : ℤ) := by exact_mod_cast hnt.symm
      rw [hfac, hnfac] at hdvd
      have hGne : ((G : ℤ)) ≠ 0 := by exact_mod_cast hG0.ne'
      have hdvd' : ((t : ℤ)) ∣ (c' : ℤ) * (((i : ℕ) : ℤ) - ((j : ℕ) : ℤ)) :=
        (mul_dvd_mul_iff_left hGne).mp hdvd
      have hcopZ : IsCoprime ((t : ℤ)) ((c' : ℤ)) := Nat.isCoprime_iff_coprime.mpr hcop
      have hdvd'' : ((t : ℤ)) ∣ (((i : ℕ) : ℤ) - ((j : ℕ) : ℤ)) :=
        hcopZ.dvd_of_dvd_mul_left hdvd'
      refine hij (Fin.ext ?_)
      have h1 : ((i : ℕ) : ℤ) - ((j : ℕ) : ℤ) = 0 := by
        refine Int.eq_zero_of_abs_lt_dvd hdvd'' ?_
        have hi := i.isLt; have hj := j.isLt
        rw [abs_lt]; omega
      omega
  have hter : t ≤ r := by
    have hLI' : LinearIndependent K (fun i : Fin t => d₀ ^ (i : ℕ)) :=
      LinearIndependent.of_comp ((IsScalarTower.toAlgHom K E M).toLinearMap) hLI
    simpa [hrdef] using hLI'.fintype_card_le_finrank
  have hGq : q ≤ G := by
    have h1 : q * t ≤ q * r := Nat.mul_le_mul_left q hter
    have h2 : q * r = G * t := by rw [mul_comm q r, hrq, ← hnt]
    have h3 : q * t ≤ G * t := by omega
    exact Nat.le_of_mul_le_mul_right h3 ht0
  have hGc : G ≤ c := Nat.le_of_dvd hcpos (Nat.gcd_dvd_right n c)
  have hcq : c = q := le_antisymm hb (le_trans hGq hGc)
  obtain ⟨m, hm⟩ := hvalE d hd
  exact ⟨m, by rw [hm, hcq]⟩

end Main

/-! ## §5 塔への代入 —— `hvalj`(層ごと `k+1` 本)を `hvalK`(底の 1 本)に落とす -/

section Tower

open GainedTowerModel GainedTowerModel.GaloisTower

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★★`hvalj` の**供給** —— 底 `K` が全分岐(`Γ_K ⊆ ‖π‖^{p^{k+1}ℤ}`)なら、
`⟨g^{p^j}⟩` の固定体 `E j` はすべて `Γ_{E j} ⊆ ‖π‖^{p^{k+1−j}ℤ}` を満たす。 -/
theorem valj_of_valK {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p k : ℕ} [Fact p.Prime] {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m)) :
    ∀ j, j ≤ k → ∀ d : twr g p j, d ≠ 0 →
      ∃ m : ℤ, ‖algebraMap (twr g p j) M d‖ = ‖π‖ ^ (((p ^ (k + 1 - j) : ℕ) : ℤ) * m) := by
  intro j hj d hd
  have hp' : p.Prime := Fact.out
  exact exists_zpow_norm_intermediate (K := K) (E := twr g p j) (M := M) hπ0 hπ1 hnK
    (finrank_twr hp'.pos g hg (by omega)) hvalK d hd

omit [IsUltrametricDist M] in
/-- `heM` と `hnormp` から `0 < ‖π‖`。 -/
theorem norm_pos_of_normp {p k e : ℕ} [Fact p.Prime] {π : M} (he : 0 < e)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e)) : 0 < ‖π‖ := by
  have hp' : p.Prime := Fact.out
  have hexp : 0 < p ^ (k + 1) * e := Nat.mul_pos (pow_pos hp'.pos _) he
  have hpow_pos : 0 < ‖π‖ ^ (p ^ (k + 1) * e) := by
    rw [← heM, hnormp]
    exact inv_pos.mpr (by exact_mod_cast hp'.pos)
  by_contra hcon
  have hz : ‖π‖ = 0 := le_antisymm (not_lt.mp hcon) (norm_nonneg π)
  rw [hz, zero_pow hexp.ne'] at hpow_pos
  exact lt_irrefl 0 hpow_pos

omit [IsUltrametricDist M] in
/-- `heM` と `hnormp` から `‖π‖ ≠ 1`。 -/
theorem norm_ne_one_of_normp {p k e : ℕ} [Fact p.Prime] {π : M}
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e)) : ‖π‖ ≠ 1 := by
  have hp' : p.Prime := Fact.out
  intro h
  have h1 : ((p : ℝ))⁻¹ = 1 := by rw [← hnormp, heM, h, one_pow]
  have hp0 : ((p : ℝ)) ≠ 0 := by
    have : (0:ℝ) < p := by exact_mod_cast hp'.pos
    exact this.ne'
  have h2 : (p : ℝ) = 1 := by field_simp at h1; exact h1.symm
  exact hp'.one_lt.ne' (by exact_mod_cast h2)

/-- ★★★★★★★**点 2 の最終形** —— 入力は
**巡回群 1 個**(`orderOf g = p^{k+1}`)、`π` の原始性、`[M:K] = p^{k+1}`、
そして★**底 `K` が全分岐であること 1 本**(`hvalK`)だけ。

★★前の形(`exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic`)の `hvalj` は
**構成された体の族 `E j` すべて(`k+1` 本)についての全分岐**を要求していた。
本定理はそれを**底 `K` の 1 本**に落とす(§4 の `exists_zpow_norm_intermediate`)。

★逸脱の記録: `hnK`(`[M:K] = p^{k+1}`)を**足した**。原典の設定では `M/K` は
`⟨τ⟩ ≅ ℤ/p^{k+1}` を Galois 群に持つので自動だが、前の形では `hvalj 0` の中に
隠れていた。★これがあると `E 0 = K` になり、**層どうしを比較せず全部底 `K` から
測れる**(#314 の構え)ので `IntermediateField` の 2 層をまたぐ `rfl`(#59)を踏まない。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_totallyRamified
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (τ : M →+ M) (hτ : ∀ z : M, τ z = g z)
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (σ : M →ₐ[twr g p k] M) (hσg : ∀ z : M, σ z = (g ^ p ^ k) z)
    (he : 0 < e)
    (hu0 : 1 ≤ u 0)
    (hult : ∀ m, m < k → u m < u (m + 1))
    (hudvd : ∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m)
    (hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hjump : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) (x : M) :
    ∃ y : twr g p k, ‖x - algebraMap (twr g p k) M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ :=
  exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic g hg τ hτ s hsg σ hσg he hu0 hult
    hudvd hbnd hjump
    (valj_of_valK g hg (norm_pos_of_normp he hnormp heM) (norm_ne_one_of_normp hnormp heM)
      hnK hvalK)
    hti hnormp heM hbreak htop x

end Tower

/-! ## §6 ★配られた字面の検算 —— 非空虚性と、`hvalK` が不分岐を排除すること -/

section Sanity

/-- ★★**`hvalK` は不分岐を排除する**(前波が挙げた反例の側の検算)。

`M/K` が不分岐なら `Γ_K = Γ_M = ‖π‖^ℤ` なので、`‖a‖ = ‖π‖` となる `a ∈ K^×` が在る。
そのとき `hvalK` は `1 = n·m` を強い、`n = 1` になる。
⇒ `n = p^{k+1} > 1` の塔では `hvalK` は**全分岐であることそのもの**であり、
`hdegj` / `htopj` / `hfixj` からは出ない(前波の診断は正しい)。 -/
theorem valK_forces_ramified {K M : Type*} [Field K] [NormedField M] [Algebra K M]
    {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1) {n : ℕ}
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (a : K) (ha : a ≠ 0) (haπ : ‖algebraMap K M a‖ = ‖π‖) : n = 1 := by
  obtain ⟨m, hm⟩ := hvalK a ha
  rw [haπ] at hm
  have h1 : (1 : ℤ) = (n : ℤ) * m :=
    (zpow_right_inj₀ hπ0 hπ1).mp (by rw [zpow_one]; exact hm)
  have hn : ((n : ℤ)) = 1 := Int.eq_one_of_dvd_one (Int.natCast_nonneg n) ⟨m, h1⟩
  exact_mod_cast hn

/-- 数値側の跳びの列 `(u 0, u 1) = (2, 8)`。 -/
def uz : ℕ → ℤ := fun m => if m = 0 then 2 else 8

/-- ★★**非空虚性(数値側)** —— `(p,k,e,i,u) = (3,1,2,8,(2,8))` は
`hu0` / `hult` / `hudvd` / `hbnd` / `hti` を**同時に**満たす。
★`k = 1`(すなわち `k ≥ 1`)で空虚でない。 -/
theorem numeric_nonvacuous :
    (1 ≤ uz 0) ∧ (∀ m, m < 1 → uz m < uz (m + 1)) ∧
      (∀ m, m < 1 → (3 : ℤ) ^ (m + 1) ∣ uz (m + 1) - uz m) ∧
      (((3 : ℤ) - 1) * uz 1 ≤ (3 : ℤ) ^ (1 + 1) * (2 : ℤ)) ∧ (((8 : ℕ) : ℤ) = uz 1) := by
  refine ⟨by norm_num [uz], fun m hm => ?_, fun m hm => ?_, by norm_num [uz], by norm_num [uz]⟩
  · interval_cases m; norm_num [uz]
  · interval_cases m; norm_num [uz]

/-- ★`hvalK`(`Γ_K ⊆ ‖π‖^{p^{k+1}ℤ}`)と `heM`(`‖p‖ = ‖π‖^{p^{k+1}·e}`)は**両立する**:
`heM` の指数はちょうど `p^{k+1}` の `e` 倍だから、`hvalK` を `a = p` で使った形と一致する。 -/
theorem valK_heM_compatible (p k e : ℕ) :
    ((p : ℤ) ^ (k + 1)) ∣ ((p : ℤ) ^ (k + 1) * (e : ℤ)) := Dvd.intro _ rfl

end Sanity

/-! ## §7 `.src`(原典の対応箇所) -/

def linearIndependent_of_norm_pairwise_notMem_coset.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_zpow_norm.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_zpow_norm_intermediate.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_totallyRamified.src :
    ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §8 使っている公理の一覧 -/

#print axioms linearIndependent_of_norm_pairwise_notMem_coset
#print axioms linearIndependent_of_ne_mod
#print axioms card_le_finrank_of_ne_mod
#print axioms exists_zpow_norm
#print axioms int_subgroup_dvd
#print axioms exists_valSub_gen
#print axioms exists_zpow_norm_intermediate
#print axioms valj_of_valK
#print axioms norm_pos_of_normp
#print axioms norm_ne_one_of_normp
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_totallyRamified
#print axioms valK_forces_ramified
#print axioms numeric_nonvacuous

end TotallyRamified

end ABC3.Found.PGC
