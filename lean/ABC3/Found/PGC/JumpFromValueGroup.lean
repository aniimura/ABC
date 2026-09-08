import ABC3.Found.PGC.TotallyRamifiedValueGroup

/-!
# [pGC] `hjump`(層ごとの跳び)を**落とす** —— 跳びの列は仮説ではなく**値群から作れる**

`Found/PGC/TotallyRamifiedValueGroup.lean` の出口
`exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_totallyRamified` に残っていた
ノルム側の仮説のうち、★**跳びに関わるもの 3 本**

    hjump  : ∀ j < k, ‖s j π − π‖ ≤ ‖π‖^{u j + 1}     (層ごと k 本)
    hbreak : ‖σ π − π‖ = ‖π‖^{i + 1}                  (頂点)
    hti    : (i : ℤ) = u k

と、整数 `i`・列 `u`・自己同型 `σ`(と `hσg`)を**すべて落とす**。

## ★★配られた見立ては**外れ**であった(§6 で機械検算)

持ち場は「`GainedTowerModel.jump_succ_of_jump_of_step`(跳びの伝播 1 段)で
`hjump` を落とせ」であった。★**これは実データで否定される。**

`M = ℚ₃(ζ₈₁)`, `K = ℚ₃(ζ₃)`(`k = 2`, `e = 2`, `E = v_M(3) = p^{k+1}·e = 54`)、
実際の跳びは `u = (2, 8, 26)`(`GainedTowerModel.Zeta81.hOrd27_breaks`)であり

| | 指数 |
|---|---|
| 伝播が出す上界 `min(E + u₀ + 1, p·u₀ + 1)` | ★`min(57, 7) = 7` |
| `hjump 1` が要求する `u₁ + 1` | ★`9` |
| `GainedJumpSeq.seq` まで緩めた要求 `min(u₁, u₂/p) + 1` | ★やはり `9` |

`‖π‖ < 1` なので `‖π‖^9 < ‖π‖^7`(`propagated_bound_strictly_weaker`)。
⇒ ★★伝播の結論は要求より**真に弱い**。さらに
「`u (j+1) ≤ min(E + u j, p·u j)`」を仮説に足す道は、実データで
`8 ≤ min(56, 6)` が**偽**(`chain_upper_false`)なので★**空虚**になる。

★伝播が緩い理由は §8 の二項展開にある: `D^p π` の評価に使えるのは
`v(D^{r+1}π) ≥ v(D^rπ) + u_j` の繰り返しだけで、`v(D^pπ) ≥ p·u_j + 1` までしか出ない。
実際の `ℚ₃(ζ₂₇)` では `D³π = (τ³π − π) − 3(τ(τπ−π))` の**打ち消し**で
`v(D³π) = 9`(`> 3·2+1 = 7`)になっており、この打ち消しは超距離だけでは見えない。

## ★★★通った道 —— 「跳びの列は**定理**である」

跳び `u j` は「`‖s j π − π‖` が `‖π‖` の何乗か」でしかない。
★底 `K` が全分岐(`hvalK`)なら `Γ_M ⊆ ‖π‖^ℤ`(`TotallyRamified.exists_zpow_norm`)なので、
`s j π − π ≠ 0` でありさえすれば**指数は自動的に存在する**。
そして `s j π = π`(`j ≤ k`)は `orderOf g = p^{k+1}` と `htop` から**排除できる**(§2)。

* §1 `exists_exponent_seq` —— ★**抽象核**(選択だけ)。各点の存在から列を作る。
* §2 `eq_one_of_apply_eq_self` / `pow_pow_ne_one` —— ★**抽象核**(体の生成元 / 純群論)。
  ★分岐・付値の語彙が 1 語も出ない。
* §3 `exists_jump_seq` —— ★★跳びの列の**構成**(等号つき)。
* §4 `..._of_cyclic_jumpEq` —— `hjump`/`hbreak`/`hti`/`i`/`σ` が落ちた出口。
* §5 `..._of_cyclic_jumpFree` —— ★★★列 `u` すら仮説から消えた最終形。
* §6 上の**否定的な検算**と、`harith` の結論が実データで成り立つこと。

## 仮説の数(実測)

| 出口 | 明示引数(`x` を除く) |
|---|---|
| `..._of_cyclic_totallyRamified`(前波) | 21 |
| ★`..._of_cyclic_jumpEq`(§4) | 17 |
| ★★`..._of_cyclic_jumpFree`(§5) | ★**13** |

★残るノルム側の仮説は `hnormp`(`p` の正規化)・`heM`(`π` の位置)・`hvalK`(全分岐)の
**3 本だけ**で、跳びに関するものは 1 本も無い。

## ★在庫の測定(コマンドを残す)

```
grep -n "ext_of_adjoin_eq_top" .cache/mathlib-index.txt
  → AlgHom.ext_of_adjoin_eq_top (h : adjoin R s = ⊤) (hs : s.EqOn φ₁ φ₂) : φ₁ = φ₂ が在る。
    ★これで「生成元を固定する自己同型は恒等」が 5 行で書ける(自作しなかった)。
grep -n "zpow_lt_zpow_right_of_lt_one₀" .cache/mathlib-index.txt
  → (ha₀ : 0 < a) (ha₁ : a < 1) (hmn : m < n) : a ^ n < a ^ m。★§6 の「真に弱い」の要。
grep -n "Nat.pow_dvd_pow_iff_le_right" .cache/mathlib-index.txt   → 在る(1 < p で ↔)。
grep -n "exists_zpow_norm" lean/ABC3/Found/PGC/TotallyRamifiedValueGroup.lean
  → 前波が作った Γ_M ⊆ ‖π‖^ℤ。★本ファイルはこれを**跳びの定義**に使う(前波は値群にしか使っていない)。
```

## ★★閉じていないもの(正確に)

1. `harith`(実際の跳びが Hasse–Arf の合同と `hbnd` を満たすこと)は**残る**。
   これは `ℤ` だけの条件で、`GainedTowerModel` §2 が `hu` に変換する。
   ★Hasse–Arf そのものを分岐理論から証明する仕事は本ファイルの外である。
2. `hvalK`(底が全分岐)は落とせない(`TotallyRamified.valK_forces_ramified`)。
3. ★**体・ノルム込みの模型**(`ℚ₃` の 18 次全分岐拡大の構成)は依然として未構成。
   §6 で検算したのは数値側だけである。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. §4 の `hjumpEq` は `j ≤ k` の**等号**である。前の形は `j < k` が不等号(`hjump`)、
   `j = k` が等号(`hbreak`)に分かれていた。★`j < k` で等号に**強めた**が、
   §5 でその等号を値群から**無料で**作るので、最終形では損をしていない。
2. `σ`(頂点の `E k`-代数準同型)を仮説で受け取る代わりに `topHom` として**構成した**。
   `fix_twr`(既存)で `commutes'` が埋まる。★原典は `σ` を Galois 群の元として与える。
3. `GainedTowerModel` / `TotallyRamifiedValueGroup` / `GainedTowerStep` は
   **読むだけ**で 1 行も書き換えていない。
-/

namespace ABC3.Found.PGC

namespace JumpFromValueGroup

open Finset

/-! ## §1 抽象核(選択だけ) —— 指数の族は「各点で存在する」から作れる -/

section ExponentSeq

variable {M : Type*} [NormedField M]

/-- ★★**抽象核 1** —— `j ≤ k` の各点でノルムが `‖π‖` の整数冪なら、
指数の**列** `w : ℕ → ℤ` が取れる。★分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem exists_exponent_seq {π : M} {k : ℕ} (f : ℕ → M)
    (h : ∀ j, j ≤ k → ∃ m : ℤ, ‖f j‖ = ‖π‖ ^ m) :
    ∃ w : ℕ → ℤ, ∀ j, j ≤ k → ‖f j‖ = ‖π‖ ^ w j := by
  classical
  refine ⟨fun j => if hj : j ≤ k then (h j hj).choose else 0, fun j hj => ?_⟩
  simp only [dif_pos hj]
  exact (h j hj).choose_spec

end ExponentSeq

/-! ## §2 抽象核(群作用だけ) —— 生成元を固定する自己同型は恒等 -/

section GroupCore

variable {K M : Type*} [Field K] [Field M] [Algebra K M]

/-- ★★**抽象核 2** —— `π` が `K` 上の生成元なら、`π` を固定する `K`-自己同型は恒等。 -/
theorem eq_one_of_apply_eq_self {π : M} (htop : Algebra.adjoin K ({π} : Set M) = ⊤)
    (h : M ≃ₐ[K] M) (hπ : h π = π) : h = 1 := by
  have hA : (h : M →ₐ[K] M) = (AlgHom.id K M) := by
    refine AlgHom.ext_of_adjoin_eq_top htop ?_
    intro x hx
    rw [Set.mem_singleton_iff] at hx
    subst hx
    simpa using hπ
  ext x
  have := AlgHom.congr_fun hA x
  simpa using this

/-- ★★**抽象核 3**(純群論) —— `orderOf g = p^{k+1}` かつ `j ≤ k` なら `g^{p^j} ≠ 1`。 -/
theorem pow_pow_ne_one {G : Type*} [Group G] {p k j : ℕ} (hp : 1 < p) (g : G)
    (hg : orderOf g = p ^ (k + 1)) (hj : j ≤ k) : g ^ p ^ j ≠ 1 := by
  intro hone
  have hdvd : p ^ (k + 1) ∣ p ^ j := by
    rw [← hg]; exact orderOf_dvd_of_pow_eq_one hone
  have := (Nat.pow_dvd_pow_iff_le_right hp).1 hdvd
  omega

/-- ★★ `π` の上で層 `j` が動く: `j ≤ k` なら `(g^{p^j}) π ≠ π`。 -/
theorem apply_pow_ne_self {p k j : ℕ} (hp : 1 < p) {g : M ≃ₐ[K] M} {π : M}
    (hg : orderOf g = p ^ (k + 1)) (hj : j ≤ k)
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) :
    (g ^ p ^ j) π ≠ π := fun hfix =>
  pow_pow_ne_one hp g hg hj (eq_one_of_apply_eq_self htop _ hfix)

end GroupCore

/-! ## §3 具体層 —— 全分岐の底から跳びの列そのものを**作る** -/

section Jumps

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★**跳びの列は仮説ではなく定理である**。

底 `K` が全分岐(`hvalK`)で `π` が生成元なら、各層 `j ≤ k` の
`s j π − π` は `0` でなく(§2)、そのノルムは `‖π‖` の整数冪(`exists_zpow_norm`)。
⇒ `‖s j π − π‖ = ‖π‖^{u j + 1}` を満たす列 `u` が**取れる**。 -/
theorem exists_jump_seq {p k : ℕ} (hp : 1 < p) {π : M} [FiniteDimensional K M]
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) :
    ∃ u : ℕ → ℤ, ∀ j, j ≤ k → ‖s j π - π‖ = ‖π‖ ^ (u j + 1) := by
  have hex : ∀ j, j ≤ k → ∃ m : ℤ, ‖s j π - π‖ = ‖π‖ ^ m := by
    intro j hj
    refine TotallyRamified.exists_zpow_norm hπ0 hπ1 hnK hvalK ?_
    refine sub_ne_zero.mpr ?_
    rw [hsg j π]
    exact apply_pow_ne_self hp hg hj htop
  obtain ⟨w, hw⟩ := exists_exponent_seq (π := π) (k := k) (fun j => s j π - π) hex
  refine ⟨fun j => w j - 1, fun j hj => ?_⟩
  simpa using hw j hj

end Jumps

/-! ## §4 出口 —— `hjump` / `hbreak` / `hti` / `i` が**全部落ちた**形 -/

section Exit

open GainedTowerModel GainedTowerModel.GaloisTower TotallyRamified

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- 頂点の自己同型 `g^{p^k}` を、その固定体 `E k = twr g p k` 上の代数準同型として見た形。 -/
def topHom {K : Type*} [Field K] [Algebra K M] (g : M ≃ₐ[K] M) (p k : ℕ) :
    M →ₐ[twr g p k] M :=
  { ((g ^ p ^ k : M ≃ₐ[K] M) : M →ₐ[K] M).toRingHom with
    commutes' := fun d => fix_twr g p k d }

omit [IsUltrametricDist M] in
@[simp] theorem topHom_apply {K : Type*} [Field K] [Algebra K M] (g : M ≃ₐ[K] M) (p k : ℕ)
    (z : M) : topHom g p k z = (g ^ p ^ k) z := rfl

/-- ★★★★★★★**点 2 の跳び無し形** —— `hjump`(層ごと `k` 本)と
`hbreak` / `hti` と整数 `i` が**すべて落ちた**。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpEq
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {e k : ℕ} {π : M} {u : ℕ → ℤ}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (τ : M →+ M) (hτ : ∀ z : M, τ z = g z)
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (he : 0 < e)
    (hu0 : 1 ≤ u 0)
    (hult : ∀ m, m < k → u m < u (m + 1))
    (hudvd : ∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m)
    (hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hjumpEq : ∀ j, j ≤ k → ‖s j π - π‖ = ‖π‖ ^ (u j + 1))
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) (x : M) :
    ∃ y : twr g p k, ‖x - algebraMap (twr g p k) M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have huk : (p : ℤ) ^ k ≤ u k :=
    pow_le_of_strictMono_of_dvd (by positivity) hu0 hult hudvd k le_rfl
  have hk0 : (0 : ℤ) ≤ u k := by
    have : (1 : ℤ) ≤ (p : ℤ) ^ k := one_le_pow₀ (by exact_mod_cast hp'.one_le)
    omega
  have hti : (((u k).toNat : ℕ) : ℤ) = u k := Int.toNat_of_nonneg hk0
  have hbreak : ‖topHom g p k π - π‖ = ‖π‖ ^ ((((u k).toNat : ℕ) : ℤ) + 1) := by
    rw [hti, topHom_apply, ← hsg k π]
    exact hjumpEq k le_rfl
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_totallyRamified
    g hg τ hτ s hsg (topHom g p k) (fun z => rfl) he hu0 hult hudvd hbnd
    (fun j hj => le_of_eq (hjumpEq j (le_of_lt hj))) hnK hvalK hti hnormp heM hbreak htop x

end Exit

/-! ## §5 ★★跳びの列すら仮説から消した形 —— 残るのは `ℤ` の条件だけ -/

section JumpFree

open GainedTowerModel GainedTowerModel.GaloisTower TotallyRamified

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★★★★★**点 2 の最終形** —— 跳びの列 `u` は仮説から**消えた**。

残る `harith` は「**実際の**跳びの列(値群から一意に決まる)が
Hasse–Arf の合同と頂点の上界を満たす」という★`ℤ` だけの条件である。
★ノルムの仮説は `hnormp` / `heM` / `hvalK` の 3 本(＝ `p` の正規化・`π` の位置・全分岐)に減った。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {e k : ℕ} {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (τ : M →+ M) (hτ : ∀ z : M, τ z = g z)
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (he : 0 < e)
    (harith : ∀ u : ℕ → ℤ, (∀ j, j ≤ k → ‖s j π - π‖ = ‖π‖ ^ (u j + 1)) →
      1 ≤ u 0 ∧ (∀ m, m < k → u m < u (m + 1)) ∧
        (∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m) ∧
        ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) (x : M) :
    ∃ y : twr g p k, ‖x - algebraMap (twr g p k) M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  obtain ⟨u, hu⟩ := exists_jump_seq hp'.one_lt g hg s hsg
    (norm_pos_of_normp he hnormp heM) (norm_ne_one_of_normp hnormp heM) hnK hvalK htop
  obtain ⟨hu0, hult, hudvd, hbnd⟩ := harith u hu
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpEq g hg τ hτ s hsg he
    hu0 hult hudvd hbnd hu hnK hvalK hnormp heM htop x

end JumpFree

/-! ## §6 ★★配られた見立ての検算 —— **伝播だけでは `hjump` は落ちない**

持ち場は「`jump_succ_of_jump_of_step`(跳びの伝播 1 段)で `hjump` を落とせ」であった。
★**これは実データで機械的に否定される。**

`M = ℚ₃(ζ₈₁)`, `K = ℚ₃(ζ₃)`(`k = 2`, `e = 2`, `E = v_M(3) = p^{k+1}·e = 54`)、
実際の跳びは `u = (2, 8, 26)`(`GainedTowerModel.Zeta81.hOrd27_breaks`)。

* 伝播が出すのは `‖s 1 π − π‖ ≤ ‖π‖^{min(E + u₀ + 1, p·u₀ + 1)} = ‖π‖^{min(57,7)} = ‖π‖^7`
* `hjump 1` が要求するのは `‖s 1 π − π‖ ≤ ‖π‖^{u₁ + 1} = ‖π‖^9`
* `‖π‖ < 1` なので `‖π‖^9 < ‖π‖^7` ⇒ ★**伝播の結論は要求より真に弱い**

さらに `GainedJumpSeq.seq` が実際に使う値まで緩めても足りない:
`t 2 = min(u₁, u₂ / p) = min(8, 26/3) = 8` なので、要求はやはり指数 `9` である。

⇒ ★★`hjump` を「底 1 本 + 伝播」に落とす道は**閉じない**。
本ファイルが取ったのは**別の道**(値群から跳びの列を作る、§3)で、
そちらは `hjump` を丸ごと落とす(§5)。 -/

section Weak

open GainedTowerModel

/-- ★実データ(`ℚ₃(ζ₈₁)/ℚ₃(ζ₃)`, `k = 2`)で、伝播の出す指数は `7`。 -/
theorem propagated_exponent : min ((54 : ℤ) + 2 + 1) (3 * 2 + 1) = 7 := by norm_num

/-- ★`GainedJumpSeq.seq` が層 1 で要求する指数は `u₁ = 8`(`min(8, 26/3) = 8`)。 -/
theorem required_exponent : min (8 : ℤ) ((26 : ℤ) / 3) = 8 := by decide

/-- ★★**配られた字面は偽** —— 「`u (j+1) ≤ min (E + u j) (p·u j)`」を仮説に足すと
実データ(`E = 54`, `u = (2,8,26)`)で**偽**になる ⇒ その形の定理は空虚である。 -/
theorem chain_upper_false : ¬ ((8 : ℤ) ≤ min ((54 : ℤ) + 2) (3 * 2)) := by norm_num

/-- ★★伝播の結論が要求より**真に弱い**ことのノルム側の確認。 -/
theorem propagated_bound_strictly_weaker {M : Type*} [NormedField M] {π : M}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) : ‖π‖ ^ (9 : ℤ) < ‖π‖ ^ (7 : ℤ) :=
  zpow_lt_zpow_right_of_lt_one₀ hπ0 hπ1 (by norm_num)

/-- ★★`harith` の結論は実データで**成り立つ**(`k = 2`, `u = (2,8,26)`) ⇒ §5 は空虚でない。 -/
theorem harith_zeta81 :
    1 ≤ Zeta81.uz 0 ∧ (∀ m, m < 2 → Zeta81.uz m < Zeta81.uz (m + 1)) ∧
      (∀ m, m < 2 → (3 : ℤ) ^ (m + 1) ∣ Zeta81.uz (m + 1) - Zeta81.uz m) ∧
      ((3 : ℤ) - 1) * Zeta81.uz 2 ≤ (3 : ℤ) ^ (2 + 1) * (2 : ℤ) :=
  ⟨Zeta81.uz_hu0, Zeta81.uz_hult, Zeta81.uz_hudvd, Zeta81.uz_hbnd⟩

end Weak

/-! ## §7 `.src`(原典の対応箇所) -/

def exists_jump_seq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §8 使っている公理の一覧 -/

#print axioms exists_exponent_seq
#print axioms eq_one_of_apply_eq_self
#print axioms pow_pow_ne_one
#print axioms apply_pow_ne_self
#print axioms exists_jump_seq
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpEq
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree
#print axioms chain_upper_false
#print axioms propagated_bound_strictly_weaker
#print axioms harith_zeta81

end JumpFromValueGroup

end ABC3.Found.PGC
