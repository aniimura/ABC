import ABC3.Found.PGC.JumpFromValueGroup
import ABC3.Found.PGC.ConcreteNormedModel
import ABC3.Found.PGC.WildDepthDescent
import Mathlib.GroupTheory.Nilpotent

/-!
# [pGC] 一般の `K` から出口に渡す層 —— **`htop` / `e` / `heM` は仮説ではなく定理**

`Found/PGC/JumpFromValueGroup.lean` の出口

    exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree   (明示引数 13)

に渡すべき組 `(M, g, π)` を、一般の底から作る仕事の**測定と、閉じた部分**である。

## ★★測ったこと(先に結論)

1. ★★**Sylow 降下が供給する層は「全分岐」ではない。**
   `WildDepthDescent.exists_pgroup_descent` が返すのは `P ≤ Q`・`[Q:P] = p`、すなわち
   Galois 対応で**次数 `p` の層** `M^P / M^Q` である(`CyclicLayerDescent.lean` の測定 2)。
   局所体の `p` 次拡大は `e·f = p` より**不分岐か全分岐のどちらか**であって、
   ★どちらかは `x` に依る。不分岐の側では出口の `hvalK` は**偽**である
   (`TotallyRamified.valK_forces_ramified`、本ファイル §6 で対偶を明示した)。
   ⇒ ★**場合分けは避けられない。本ファイルは全分岐の側だけを閉じる。**
2. ★層は 1 段(`[Q:P] = p`)しか出ないので、出口に渡すのは `k = 0` である。
   ★一般の `k` の塔を一度に取る道は Sylow 降下からは出ない(§5)。
3. ★★残るノルム側の仮説 3 本(`hnormp` / `heM` / `hvalK`)のうち、
   ★**`heM`(と `e`・`he`)は `hvalK` + `hnormp` + `‖π‖ < 1` から出る**(§3)。
   ★**`htop`(`M = K(π)`)は `hvalK` から出る**(§2)。
   ⇒ 出口の明示引数は ★**13 → 7** に落ちる(§4 `..._of_uniformizer`)。
4. ★★**群論の段は閉じた**(§8)—— 前波(`CyclicLayerDescent.lean` 冒頭)が
   「★測ったが書いていない」と残した「`p` 群の指数 `p` の部分群は正規」を埋め、
   `exists_pgroup_descent` の出力を★「`P ⊴ Q` かつ `Q/P` は位数 `p` の**巡回**群」まで
   強めた(`exists_pgroup_descent_cyclic`)。★これで `orderOf g = p` の**群の側**は揃う。

## ★閉じていないもの(正確に)

* ★**不分岐側**(`f = p`)の 1 段は本ファイルでは閉じていない。
  `SenLemma.descentStep_of_natDegree_tame` は「`p ∤ deg minpoly`」を扱うのであって、
  「層が不分岐」とは**別の条件**である(前者は `x` の次数、後者は層の分岐)。★同一視しないこと。
* ★`harith`(実際の跳びが `1 ≤ u` と `(p−1)u ≤ p^{k+1}e` を満たすこと)は**残る**。
  `k = 0` ではこれは古典的な「暴分岐の跳びの上下界 `1 ≤ u ≤ pe/(p−1)`」そのもので、
  ★下界 `1 ≤ u` の証明には**剰余体**(`ū^p = 1 ⇒ ū = 1`)と `g` が等長であることが要る。
  ★本ファイルの設定(`NormedField M` と `g : M ≃ₐ[K] M`)は `‖g z‖ = ‖z‖` を**含まない**ので、
  そこで止まる。★「測ったが書いていない」と読むこと。
* ★`hnormp`(`‖p‖ = 1/p`)は落ちない。ノルムの正規化そのものである。
* ★★**体の側の翻訳**(Galois 対応で `M^P / M^Q` を `IntermediateField` として作り、
  §8 の `Q/P` の生成元を `g : M^P ≃ₐ[M^Q] M^P` に移すこと)は**書いていない**。
  ★`lean-idioms.md` #59(中間体 2 層をまたぐ `rfl` が kernel を止める)に当たる区間なので、
  ★`MulAction.stabilizer` の指数で測る #296 の構えが要る。★測っただけである。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. §4 `..._of_uniformizer` は `e` / `he` / `heM` / `htop` を落とす代わりに
   ★`hπlt : ‖π‖ < 1`(π が素元であること)を**足している**。
   これは追加の仮定ではなく**入れ替え**である: 元の形では `he` + `heM` + `hnormp` から
   `0 < ‖π‖`・`‖π‖ ≠ 1` が出ていた(`TotallyRamified.norm_pos_of_normp` /
   `norm_ne_one_of_normp`)が、`‖π‖ < 1` 自体は `heM` の指数が自然数であることから従う。
   逆向き(本ファイル)は `‖π‖ < 1` から `e` を作る。★両者は同値である。
2. §4 の `harith` は `e` を**全称量化**した形にした(`e` を落とした代償)。
   `‖π‖ < 1` のとき `e` は一意なので強さは変わらない。
3. `τ`(加法群準同型)と `s`(環準同型の族)は `g` から**構成した**ので仮説から消えた
   (`g.toRingEquiv.toRingHom` を渡すだけ。★`fun _ _ => rfl` で埋まる)。
4. `JumpFromValueGroup` / `TotallyRamifiedValueGroup` / `WildDepthDescent` は
   **読むだけ**で 1 行も書き換えていない。

## ★在庫の測定(コマンドを残す)

```
grep -n "span_eq_top_of_card_eq_finrank" .cache/mathlib-index.txt
  → LinearIndependent.span_eq_top_of_card_eq_finrank' [FiniteDimensional K V] が在る。
grep -n "toSubmodule_eq_top" .cache/mathlib-index.txt
  → Algebra.toSubmodule_eq_top : Subalgebra.toSubmodule S = ⊤ ↔ S = ⊤。
grep -n "self_mem_adjoin_singleton" .cache/mathlib-index.txt
  → Algebra.self_mem_adjoin_singleton R s が在る(★索引には用例の行でしか出ない)。
grep -n "eq_zero_of_abs_lt_dvd" .cache/mathlib-index.txt
  → Int.eq_zero_of_abs_lt_dvd (h1 : m ∣ x) (h2 : |x| < m) : x = 0。
grep -n "zpow_lt_one_iff_right_of_lt_one₀" .cache/mathlib-index.txt
  → (ha₀ : 0 < a) (ha₁ : a < 1) : a ^ n < 1 ↔ 0 < n。★§3 の要。
```

## ★配管(実測した罠、`lean-idioms.md` に節を足した)

`M` が `NormedField` のとき `Algebra.adjoin K {π}` への**所属**を `exact` で埋めると

    error: (deterministic) timeout at `isDefEq`, maximum number of heartbeats (200000)

になる(21 秒)。★同じ証明を `[Field M]` だけの補題(`adjoin_eq_top_of_span_powers`)に
切り出して**代入**すると 10 秒で通る。★これが本ファイルで抽象核を切った実利上の理由である。
-/

namespace ABC3.Found.PGC

namespace TotallyRamifiedLayer

open TotallyRamified

/-! ## §1 抽象核 —— 分岐・付値・Galois が 1 語も出ない 2 本 -/

/-- ★**抽象核 1**(純 `ℕ`/`ℤ`) —— `i, j < n`, `i ≠ j` なら `n ∤ i − j`。 -/
theorem not_dvd_sub_of_lt {n i j : ℕ} (hi : i < n) (hj : j < n) (hij : i ≠ j) :
    ¬ ((n : ℤ) ∣ ((i : ℤ) - (j : ℤ))) := by
  intro hdvd
  have habs : |((i : ℤ) - (j : ℤ))| < (n : ℤ) := by rw [abs_lt]; omega
  have := Int.eq_zero_of_abs_lt_dvd hdvd habs
  omega

/-- ★★**抽象核 2**(ノルムが 1 語も出ない) —— `π` の冪 `n` 本が `M` を張れば `M = K(π)`。

★`M` を `NormedField` のまま書くと `isDefEq` が 200000 heartbeat で落ちる(冒頭の配管)。
`[Field M]` に落として書き、§2 で代入する。 -/
theorem adjoin_eq_top_of_span_powers {K M : Type*} [Field K] [Field M] [Algebra K M]
    [FiniteDimensional K M] {π : M} {n : ℕ}
    (hspan : Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ))) = ⊤) :
    Algebra.adjoin K ({π} : Set M) = ⊤ := by
  have hle : Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ)))
      ≤ Subalgebra.toSubmodule (Algebra.adjoin K ({π} : Set M)) := by
    refine Submodule.span_le.2 ?_
    rintro z ⟨l, rfl⟩
    exact pow_mem (Algebra.self_mem_adjoin_singleton K π) (l : ℕ)
  rw [hspan] at hle
  exact Algebra.toSubmodule_eq_top.mp (top_le_iff.mp hle)

/-! ## §2 ★★`htop` は仮説ではない —— 全分岐なら素元は生成元 -/

section Top

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★**`htop` の供給** —— 底が全分岐(`Γ_K ⊆ ‖π‖^{nℤ}`, `n = [M:K]`)なら `M = K(π)`。

`π^0, …, π^{n−1}` はノルムの指数が `n` を法として相異なるので
`TotallyRamified.linearIndependent_of_ne_mod` で 1 次独立、本数が `[M:K]` に等しいので
`M` を張る。★出口の仮説 `htop` はこれで**落ちる**。 -/
theorem adjoin_eq_top_of_valK [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    Algebra.adjoin K ({π} : Set M) = ⊤ := by
  have hli : LinearIndependent K (fun l : Fin n => π ^ (l : ℕ)) := by
    refine linearIndependent_of_ne_mod hπ0 hπ1 hvalK _ (fun l => ((l : ℕ) : ℤ)) ?_ ?_
    · intro l; simp [zpow_natCast]
    · intro i j hij
      exact not_dvd_sub_of_lt i.isLt j.isLt (fun h => hij (Fin.ext h))
  have hcard : Fintype.card (Fin n) = Module.finrank K M := by simp [hn]
  exact adjoin_eq_top_of_span_powers (hli.span_eq_top_of_card_eq_finrank' hcard)

end Top

/-! ## §3 ★★`e` / `he` / `heM` は仮説ではない —— 絶対分岐指数は値群が決める -/

section AbsRam

variable {K M : Type*} [Field K] [NormedField M] [Algebra K M]

omit [Algebra K M] in
/-- `hnormp` から `(p : M) ≠ 0`。 -/
theorem nat_ne_zero_of_normp {p : ℕ} [Fact p.Prime] (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹) :
    (p : M) ≠ 0 := by
  have hp' : p.Prime := Fact.out
  have : (0:ℝ) < (p:ℝ) := by exact_mod_cast hp'.pos
  intro h
  rw [h, norm_zero] at hnormp
  exact (inv_pos.mpr this).ne hnormp

/-- ★`0 < ‖π‖` は `hvalK` + `hnormp` から出る(`π = 0` なら `‖p‖ ∈ {0, 1}` になる)。 -/
theorem norm_pos_of_valK {p n : ℕ} [Fact p.Prime] {π : M}
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    0 < ‖π‖ := by
  have hp' : p.Prime := Fact.out
  have hpM : (p : M) ≠ 0 := nat_ne_zero_of_normp (p := p) hnormp
  have hpK : (p : K) ≠ 0 := by
    intro h
    exact hpM (by rw [← map_natCast (algebraMap K M) p, h, map_zero])
  obtain ⟨m, hm⟩ := hvalK (p : K) hpK
  rw [map_natCast, hnormp] at hm
  rcases (norm_nonneg π).lt_or_eq with h | h
  · exact h
  · exfalso
    rw [← h] at hm
    have hpR : (1:ℝ) < (p:ℝ) := by exact_mod_cast hp'.one_lt
    rcases eq_or_ne ((n : ℤ) * m) 0 with h0 | h0
    · rw [h0, zpow_zero] at hm
      have : (p:ℝ) = 1 := by field_simp at hm; linarith [hm]
      linarith
    · rw [zero_zpow _ h0] at hm
      have hpos : (0:ℝ) < (p:ℝ)⁻¹ := by positivity
      rw [hm] at hpos
      linarith

/-- ★★★**`e` と `heM` の供給** —— 底が全分岐(`Γ_K ⊆ ‖π‖^{nℤ}`)でノルムが正規化されていて
`π` が素元(`‖π‖ < 1`)なら、★**絶対分岐指数 `e > 0` が存在して `‖p‖ = ‖π‖^{n·e}`**。

★出口の `e` / `he` / `heM` はこれで**落ちる**(仮説ではなく結論になる)。 -/
theorem exists_absRamIndex {p n : ℕ} [Fact p.Prime] {π : M} (hn0 : 0 < n)
    (hπlt : ‖π‖ < 1)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    ∃ e : ℕ, 0 < e ∧ ‖(p : M)‖ = ‖π‖ ^ (n * e) := by
  have hp' : p.Prime := Fact.out
  have hπ0 : 0 < ‖π‖ := norm_pos_of_valK (p := p) (n := n) hnormp hvalK
  have hpM : (p : M) ≠ 0 := nat_ne_zero_of_normp (p := p) hnormp
  have hpK : (p : K) ≠ 0 := by
    intro h
    exact hpM (by rw [← map_natCast (algebraMap K M) p, h, map_zero])
  obtain ⟨m, hm⟩ := hvalK (p : K) hpK
  rw [map_natCast, hnormp] at hm
  have hlt1 : ‖π‖ ^ ((n : ℤ) * m) < 1 := by
    rw [← hm]
    have hpR : (1:ℝ) < (p:ℝ) := by exact_mod_cast hp'.one_lt
    rw [inv_lt_one_iff₀]; right; exact hpR
  have hpos : 0 < (n : ℤ) * m := (zpow_lt_one_iff_right_of_lt_one₀ hπ0 hπlt).mp hlt1
  have hm0 : 0 < m := by
    by_contra hcon
    rw [not_lt] at hcon
    have hnZ : (0:ℤ) < (n : ℤ) := by exact_mod_cast hn0
    nlinarith [hpos, hnZ]
  refine ⟨m.toNat, by omega, ?_⟩
  rw [hnormp, hm, ← zpow_natCast ‖π‖ (n * m.toNat)]
  congr 1
  push_cast [Int.toNat_of_nonneg hm0.le]
  ring

end AbsRam

/-! ## §4 ★★★出口 —— 明示引数 13 → 7 -/

section Exit

open GainedTowerModel GainedTowerModel.GaloisTower TotallyRamified JumpFromValueGroup

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★**`htop` / `τ` / `s` を落とした形**(明示引数 8、`x` を除いて 7)。

`τ`(加法群準同型)と `s`(環準同型の族)は `g` から作れるので仮説から消え、
`htop` は §2 が供給する。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_totallyRamified
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {e k : ℕ} {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1)) (he : 0 < e)
    (harith : ∀ u : ℕ → ℤ, (∀ j, j ≤ k → ‖(g ^ p ^ j) π - π‖ = ‖π‖ ^ (u j + 1)) →
      1 ≤ u 0 ∧ (∀ m, m < k → u m < u (m + 1)) ∧
        (∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m) ∧
        ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e)) (x : M) :
    ∃ y : twr g p k, ‖x - algebraMap (twr g p k) M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖g x - x‖ := by
  have htop : Algebra.adjoin K ({π} : Set M) = ⊤ :=
    adjoin_eq_top_of_valK (norm_pos_of_normp he hnormp heM)
      (norm_ne_one_of_normp hnormp heM) hnK hvalK
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree
    (p := p) (e := e) (k := k) (π := π) g hg
    (g.toRingEquiv.toRingHom.toAddMonoidHom) (fun _ => rfl)
    (fun j => (g ^ p ^ j).toRingEquiv.toRingHom) (fun _ _ => rfl)
    he harith hnK hvalK hnormp heM htop x

/-- ★★★★**素元だけを渡す形** —— `e` / `he` / `heM` / `htop` / `τ` / `s` が全部落ちた。

★明示引数は `g, hg, harith, hnK, hvalK, hnormp, hπlt`(と `x`)の **7 本**。
残るのは「`M/K` は `g` で巡回・次数 `p^{k+1}`」「全分岐で `π` が素元」「`‖p‖ = 1/p`」
「跳びが Hasse–Arf の合同と上界を満たす」だけである。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {k : ℕ} {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (harith : ∀ (e : ℕ) (u : ℕ → ℤ), 0 < e → ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e) →
      (∀ j, j ≤ k → ‖(g ^ p ^ j) π - π‖ = ‖π‖ ^ (u j + 1)) →
      1 ≤ u 0 ∧ (∀ m, m < k → u m < u (m + 1)) ∧
        (∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m) ∧
        ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹) (hπlt : ‖π‖ < 1) (x : M) :
    ∃ y : twr g p k, ‖x - algebraMap (twr g p k) M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖g x - x‖ := by
  have hp' : p.Prime := Fact.out
  obtain ⟨e, he, heM⟩ :=
    exists_absRamIndex (p := p) (n := p ^ (k + 1)) (pow_pos hp'.pos _) hπlt hnormp hvalK
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_totallyRamified g hg he
    (fun u hu => harith e u he heM hu) hnK hvalK hnormp heM x

end Exit

/-! ## §5 ★★★1 段(`k = 0`)—— Sylow 降下が実際に供給する形

`WildDepthDescent.exists_pgroup_descent` は `[Q:P] = p` の**1 段**しか返さないので、
出口に渡すのは `k = 0` である。★このとき `harith` の中間 2 条件(狭義単調・Hasse–Arf の合同)は
空虚になり、残るのは★**跳びの上下界 `1 ≤ u ≤ p·e/(p−1)` の 1 本だけ**である。 -/

section DegP

open GainedTowerModel GainedTowerModel.GaloisTower TotallyRamified JumpFromValueGroup

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★★**次数 `p` の全分岐巡回層 1 枚**に対する出口。

★仮説は `g` が位数 `p`・`[M:K] = p`・全分岐(`hvalK`)・正規化(`hnormp`)・素元(`hπlt`)、
そして★**跳びの上下界 `hbreak` 1 本**だけである。 -/
theorem exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p)
    (hbreak : ∀ (e : ℕ) (t : ℤ), 0 < e → ‖(p : M)‖ = ‖π‖ ^ (p * e) →
      ‖g π - π‖ = ‖π‖ ^ (t + 1) → 1 ≤ t ∧ ((p : ℤ) - 1) * t ≤ (p : ℤ) * (e : ℤ))
    (hnK : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹) (hπlt : ‖π‖ < 1) (x : M) :
    ∃ y : twr g p 0, ‖x - algebraMap (twr g p 0) M y‖ ≤ axDecay p 1 * ‖g x - x‖ := by
  have hprod : (∏ j ∈ Finset.Icc 1 (0 + 1), axDecay p j) = axDecay p 1 := by simp
  rw [← hprod]
  refine exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer (k := 0) (π := π)
    g (by simpa using hg) ?_ (by simpa using hnK) (by simpa using hvalK) hnormp hπlt x
  intro e u he heM hu
  have hu0 : ‖g π - π‖ = ‖π‖ ^ (u 0 + 1) := by simpa using hu 0 le_rfl
  have heM' : ‖(p : M)‖ = ‖π‖ ^ (p * e) := by simpa using heM
  obtain ⟨h1, h2⟩ := hbreak e (u 0) he heM' hu0
  exact ⟨h1, fun m hm => absurd hm (Nat.not_lt_zero m),
    fun m hm => absurd hm (Nat.not_lt_zero m), by simpa using h2⟩

end DegP

/-! ## §6 ★★測定 —— 「全分岐」は自動ではない(不分岐側では出口の仮説が**偽**)

`WildDepthDescent.exists_pgroup_descent` が返す層 `M^P / M^Q` は次数 `p` であり、
局所体の `p` 次拡大は `e·f = p` から**不分岐か全分岐のどちらか**である。
★不分岐の側では `π` が `‖K^×‖` に入ってしまい、出口の `hvalK` は偽になる。
以下はその対偶(`TotallyRamified.valK_forces_ramified` の言い換え)。 -/

section Obstruction

/-- ★★**不分岐なら出口は使えない** —— `‖a‖ = ‖π‖` なる `a ∈ K^×` が在れば、
`n > 1` のとき `hvalK` は**偽**である。

⇒ ★「Sylow 降下が層を供給するから一般化は済む」は**外れ**。
供給される `p` 次の層が全分岐かどうかは `x` に依り、不分岐の側は**別の議論**が要る。 -/
theorem not_valK_of_norm_eq {K M : Type*} [Field K] [NormedField M] [Algebra K M]
    {π : M} {n : ℕ} (hn : 1 < n) (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    {a : K} (ha : a ≠ 0) (haπ : ‖algebraMap K M a‖ = ‖π‖) :
    ¬ (∀ b : K, b ≠ 0 → ∃ m : ℤ, ‖algebraMap K M b‖ = ‖π‖ ^ ((n : ℤ) * m)) := by
  intro h
  have hn1 : n = 1 := TotallyRamified.valK_forces_ramified hπ0 hπ1 h a ha haπ
  omega

end Obstruction

/-! ## §7 ★★非空虚性 —— 7 引数の形が具体的な体で満たされる

`Found/PGC/ConcreteNormedModel.lean` の `ℚ₂(√2)/ℚ₂`(`p = 2`, `k = 0`, `e = 1`)を
★**§5 の 1 段の出口に代入する**。★仮説を 1 つも残さない。 -/

section Nonvacuous

open ConcreteNormedModel GainedTowerModel.GaloisTower

/-- ★★**§5 の形は空虚でない** —— `ℚ₂(√2)/ℚ₂` が仮説を全部満たす。 -/
theorem model_deg_p_exists (x : M2) :
    ∃ y : twr g2 2 0, ‖x - algebraMap (twr g2 2 0) M2 y‖ ≤ axDecay 2 1 * ‖g2 x - x‖ := by
  refine exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p (p := 2) (π := pi2)
    g2 orderOf_g2 ?_ (by simpa using finrank_M2) valK norm_two_M2 norm_pi2_lt_one x
  intro e t he _ ht
  have htt : t = 2 := u0_eq_two (fun _ => t) ht
  subst htt
  refine ⟨by norm_num, ?_⟩
  have : (1 : ℤ) ≤ (e : ℤ) := by exact_mod_cast he
  norm_num
  linarith

end Nonvacuous

/-! ## §8 ★★★群論の段 —— 「`p` 群の指数 `p` の部分群は正規」

`Found/PGC/CyclicLayerDescent.lean` の冒頭が
「mathlib に `Sylow.exists_subgroup_card_pow_succ` は在るが★正規性を返さない★ので、
『`p` 群の指数 `p` の部分群は正規』を自前で書く必要がある。★測ったが書いていない。」
と残していた穴を★**埋める**。★分岐・付値・体の語彙が 1 語も出ない。

★これで `WildDepthDescent.exists_pgroup_descent` の出力は
「`P ⊴ Q` かつ `Q/P` は位数 `p` の**巡回**群」まで強くなる ——
Galois 対応で `M^P/M^Q` が★**巡回 `p` 次**の拡大であること(群の側)が閉じた。
★残るのは体の側(`IntermediateField` への翻訳)と、★**その層が全分岐かどうか**(§6)である。 -/

section PGroupLayer

/-- ★**抽象核**(純群論) —— 素数指数の部分群は極大(coatom)。 -/
theorem isCoatom_of_index_prime {G : Type*} [Group G] {H : Subgroup G} {p : ℕ}
    (hp : p.Prime) (h : H.index = p) : IsCoatom H := by
  constructor
  · intro hcon
    rw [hcon, Subgroup.index_top] at h
    exact hp.one_lt.ne h
  · intro K hK
    have hle : H ≤ K := le_of_lt hK
    have hmul : H.index = H.relIndex K * K.index := (Subgroup.relIndex_mul_index hle).symm
    rw [h] at hmul
    rcases (Nat.Prime.eq_one_or_self_of_dvd hp K.index ⟨_, by rw [hmul]; ring⟩) with h1 | h1
    · exact Subgroup.index_eq_one.mp h1
    · exfalso
      rw [h1] at hmul
      have hrel : H.relIndex K = 1 := by
        have hp0 : 0 < p := hp.pos
        nlinarith [hmul, hp0]
      have : K ≤ H := Subgroup.relIndex_eq_one.mp hrel
      exact absurd (le_antisymm hle this) (ne_of_lt hK)

/-- ★`[Q:P] = p` かつ `P` が `p` 群なら `Q` も `p` 群。 -/
theorem isPGroup_of_relIndex_prime {G : Type*} [Group G] [Finite G] {p : ℕ} [Fact p.Prime]
    {P Q : Subgroup G} (hPQ : P ≤ Q) (hPp : IsPGroup p P) (h : P.relIndex Q = p) :
    IsPGroup p Q := by
  obtain ⟨a, ha⟩ := IsPGroup.iff_card.mp hPp
  have hcard : Nat.card (P.subgroupOf Q) * (P.subgroupOf Q).index = Nat.card Q :=
    Subgroup.card_mul_index _
  have h1 : Nat.card (P.subgroupOf Q) = Nat.card P :=
    Nat.card_congr (Subgroup.subgroupOfEquivOfLe hPQ).toEquiv
  have h2 : (P.subgroupOf Q).index = p := h
  refine IsPGroup.of_card (n := a + 1) ?_
  rw [← hcard, h1, h2, ha, pow_succ]

/-- ★★★**`p` 群の指数 `p` の部分群は正規**(前波が「測ったが書いていない」と残した穴)。

`p` 群は冪零(`IsPGroup.isNilpotent`)、冪零群は normalizer 条件を満たす
(`Group.normalizerCondition_of_isNilpotent`)、normalizer 条件下では極大部分群は正規
(`Subgroup.NormalizerCondition.normal_of_coatom`)。

★在庫の測定: `Subgroup.NormalizerCondition.normal_of_coatom` は索引の行に `H` が出ないが
**明示引数**である(section の `variable (H : Subgroup G)`。`lean-idioms.md` #297 の形)。
`_` を 1 つ足して `normal_of_coatom _ hnc hmax` と書く。 -/
theorem normal_of_relIndex_prime {G : Type*} [Group G] {p : ℕ} [Fact p.Prime]
    {P Q : Subgroup G} [Finite Q] (hQ : IsPGroup p Q) (h : P.relIndex Q = p) :
    (P.subgroupOf Q).Normal := by
  haveI := IsPGroup.isNilpotent (p := p) hQ
  exact Subgroup.NormalizerCondition.normal_of_coatom _ Group.normalizerCondition_of_isNilpotent
    (isCoatom_of_index_prime Fact.out h)

/-- `[Q:P] = p` は `|Q/P| = p` そのもの(`Subgroup.relIndex` の定義)。 -/
theorem card_quotient_of_relIndex {G : Type*} [Group G] {p : ℕ} {P Q : Subgroup G}
    (h : P.relIndex Q = p) : Nat.card (↥Q ⧸ P.subgroupOf Q) = p := h

/-- ★★**商は位数 `p` の巡回群** —— 生成元 `g` の位数はちょうど `p`。

★出口 `..._of_uniformizer_deg_p` が要求する `orderOf g = p` の**群の側**がこれで揃う。 -/
theorem exists_orderOf_eq_of_relIndex_prime {G : Type*} [Group G] {p : ℕ} [Fact p.Prime]
    {P Q : Subgroup G} [Finite Q] (hQ : IsPGroup p Q) (h : P.relIndex Q = p) :
    letI := normal_of_relIndex_prime hQ h
    ∃ g : (↥Q ⧸ P.subgroupOf Q), orderOf g = p ∧ ∀ z, z ∈ Subgroup.zpowers g := by
  haveI := normal_of_relIndex_prime hQ h
  haveI : IsCyclic (↥Q ⧸ P.subgroupOf Q) := isCyclic_of_prime_card (card_quotient_of_relIndex h)
  obtain ⟨g, hg⟩ := IsCyclic.exists_generator (α := (↥Q ⧸ P.subgroupOf Q))
  exact ⟨g, (orderOf_eq_card_of_forall_mem_zpowers hg).trans (card_quotient_of_relIndex h), hg⟩

/-- ★★★★**Sylow 降下の出力を強めた形** —— `exists_pgroup_descent` の `(P, Q)` について
★`P ⊴ Q` かつ `|Q/P| = p`。⇒ Galois 対応で `M^P/M^Q` は**巡回 `p` 次**である。

★★ただし「全分岐」は**言えない**(§6)。★そこが一般化の止まる場所である。 -/
theorem exists_pgroup_descent_cyclic {G : Type*} [Group G] [Finite G] {p : ℕ} [Fact p.Prime]
    (H : Subgroup G) (k : ℕ) (hk : padicValNat p H.index = k + 1) :
    ∃ P Q : Subgroup G, P ≤ H ∧ P ≤ Q ∧ padicValNat p P.index = k + 1 ∧
      padicValNat p Q.index = k ∧ (P.subgroupOf Q).Normal ∧
      Nat.card (↥Q ⧸ P.subgroupOf Q) = p := by
  obtain ⟨P, Q, hPH, hPp, hPQ, hrel, hPv, hQv⟩ :=
    ABC3.Found.PGC.exists_pgroup_descent (p := p) H k hk
  have hQp : IsPGroup p Q := isPGroup_of_relIndex_prime hPQ hPp hrel
  exact ⟨P, Q, hPH, hPQ, hPv, hQv, normal_of_relIndex_prime hQp hrel, hrel⟩

end PGroupLayer

/-! ## §9 `.src`(原典の対応箇所) -/

def adjoin_eq_top_of_valK.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_absRamIndex.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §10 使っている公理の一覧 -/

#print axioms not_dvd_sub_of_lt
#print axioms adjoin_eq_top_of_span_powers
#print axioms adjoin_eq_top_of_valK
#print axioms norm_pos_of_valK
#print axioms exists_absRamIndex
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_totallyRamified
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer
#print axioms exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p
#print axioms not_valK_of_norm_eq
#print axioms model_deg_p_exists
#print axioms isCoatom_of_index_prime
#print axioms isPGroup_of_relIndex_prime
#print axioms normal_of_relIndex_prime
#print axioms exists_orderOf_eq_of_relIndex_prime
#print axioms exists_pgroup_descent_cyclic

end TotallyRamifiedLayer

end ABC3.Found.PGC
