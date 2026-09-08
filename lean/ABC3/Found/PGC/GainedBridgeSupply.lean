import ABC3.Found.PGC.GainedDescentBridge

/-!
# [pGC] ★★★★★`GainedDescentBridge` の仮説を落とす —— 6 本が定義と正規化から出る

`Found/PGC/GainedDescentBridge.lean` は `AxLemmaGraded` の形

```
exists_norm_sub_algebraMap_le_prod_axDecay :
  ∃ y ∈ F,  ‖x − y‖ ≤ (∏_{j∈[1,k+1]} axDecay p j) · ‖τ x − x‖
```

まで閉じたが、docstring の「閉じていないもの」に

> 1. ★**`hstep` / `hlayer` / `hbreak` / `hval` / `hchar` は仮説のままである。**
> 2. ★**`hlayer`(`(p−1)t_j ≤ e_M`)と `hlayerZ`(`(p−1)t_j ≤ p^j e`)を別々に取っている。**

の 2 点が残っていた。★★**本ファイルはこのうち 6 本を落とす**(点 2 の重複も消える)。

## ★★★測定 —— 木に在ったもの / 無かったもの(コマンドつき)

| 部品 | 在庫 | 測定 |
|---|---|---|
| `hchar`(`0<j<p` で `‖j‖=1`) | ★**木に在った** | `grep -n "norm_natCast_eq_one_of_lt_prime" lean/ABC3/Found/PGC/MinpolyOrbitSplit.lean` |
| `‖σx−x‖·‖π‖ = ‖σπ−π‖·‖x−a₀‖` | ★★**木に在った** | `grep -n "norm_algHom_sub_mul_norm_eq" lean/ABC3/Found/PGC/CyclicJumpNorm.lean`(504 行) |
| `‖x−a₀‖ ≤ ‖x−c‖`(最良近似) | ★木に在った | 同上 354 行 `norm_digitSum_sub_digit_zero_le` |
| `σ ∈ G_1 ⟹ σ^{p^k} ∈ G_{1+k}` | 木に在るが**使えない** | `RamificationJumpDivisibility.lean:344`。★下の「使わなかった理由」 |
| 上付き/下付き分岐群・Herbrand | mathlib に**無い** | 木が自前で持つ(`LowerRamificationGroup.lean:265`) |
| `AlgHom.coe_pow` | mathlib に在った | `grep -n "AlgHom.coe_pow" .cache/mathlib-index.txt` |

## ★★★落ちた 6 本と、その代わりに入った 2 本

| 落ちた仮説 | どこから出るか |
|---|---|
| `hchar : ∀ j, 0<j → j<p → ‖(j:M)‖ = 1` | `hnormp` ⟹ `‖(p:M)‖<1` ⟹ Bézout(`MinpolyOrbitSplit`) |
| ★`hlayer : ∀ j, ‖(p:M)‖ ≤ (‖π‖^{t_{j+1}})^{p−1}` | ★`heM` ＋ `hlayerZ`(★`j ≤ k` に限る) |
| `hπE : ‖π‖ = p^{−1/(p^{k+1}e)}` | `hnormp` ＋ `heM`(`E` 乗根を取るだけ) |
| `hπ0 : 0 < ‖π‖` / `hπ1 : ‖π‖ < 1` | `heM` ＋ `hnormp`(`‖π‖^E = p⁻¹ < 1`) |
| ★★`hstep : ∀ j z, ‖τ^{p^j}z − z‖ ≤ ‖π‖^{t_{j+1}}‖z‖` | ★★`hjump`(★**生成元 `π` の上だけ**の条件) |

入ったのは★**正規化 2 本だけ**である:

* `hnormp : ‖(p:M)‖ = (p:ℝ)⁻¹` —— `M` のノルムが `ℚ_p` のノルムを延長する(約束事)、
* `heM : ‖(p:M)‖ = ‖π‖^{p^{k+1}e}` —— ★`π` が `M` の素元で `e_M = p^{k+1}e`(塔が全分岐)。

## ★★★★`hstep` が落ちた仕組み(★これが本ファイルの核)

下付き分岐群の定義は `𝒪_M = 𝒪_F[π]` のとき**生成元 1 個で書ける**:

```
τ^{p^j} ∈ Γ_{t}  ⟺  v_M(τ^{p^j}π − π) ≥ t + 1.
```

★木の `LowerRamificationGroup.mem_lowerRamificationGroup_iff`(277 行)＋
`mem_pow_maximalIdeal_iff_lt_addVal`(285 行)がちょうどこの同値を出している。

「全元 `z ∈ M` について `v_M(τ z − z) − v_M(z) ≥ t`」はその**帰結**であって仮定ではない。
★★しかも木の定義側 `mem_lowerRamificationGroup_iff_forall` は
`∀ x ∈ 𝒪_M, σx − x ∈ 𝔪^{t+1}`(**加法的**、整数環の上だけ)であり、
`hstep` が要求する `‖σz − z‖ ≤ ‖π‖^t‖z‖`(**斉次**、`M` 全体)はそれより強い。
★本ファイルはその強い方を生成元の条件から出す。
★木の `CyclicJumpNorm.norm_algHom_sub_mul_norm_eq` が既に**等式**

```
‖σ x − x‖ · ‖π‖ = ‖σ π − π‖ · ‖x − a₀‖      (x = Σ_{j<p} a_j π^j)
```

を持っていたので、`‖x − a₀‖ ≤ ‖x‖`(`norm_digitSum_sub_digit_zero_le` に `c = 0`)を掛けて
`‖π‖` で割るだけで出る(`norm_algHom_sub_le_mul_of_break`)。
★★**新しい数学は 1 行もない。木の等式の右辺を 1 度緩めただけである。**

## ★★`RamificationJumpDivisibility` を使わなかった理由(★持ち場の見立ては外れた)

持ち場は「`hstep` の部品は `RamificationJumpDivisibility.lean` にあるはず」と見ていた。
★測ると `pow_char_pow_mem_lowerRamificationGroup`(344 行)は

```
σ ∈ G_1  ⟹  σ^{p^k} ∈ G_{1+k}
```

しか出さない。これは `t_{j+1} = 1 + j` に当たり、★**`hstepZ : p·t_{j+1} ≤ t_{j+2}` を満たさない**
(`p·(1+j) ≤ 2+j` は `j ≥ 0`, `p ≥ 2` で偽)。
★★**塔の跳びが `p` 倍で伸びること(`hstepZ`)は `G_n/G_{n+1}` が指数 `p` であることからは出ない。**
それは Hasse–Arf 側の内容で、`GainedTowerDescent` が `ℤ` の側の仮定として持っている。
⇒ 本ファイルは `hstep` を**跳びの列 `t` を与えられた上で**生成元の条件に落とすところまでを担当する。

## ★★★落とせなかった仮説と、その理由(★正確に)

1. ★**`hval : ∀ c ∈ F^×, ‖c‖ ∈ ‖π‖^{pℤ}`** は落ちない。これは「`F` が離散付値体で
   `M/F` が `e = p` の全分岐」という**入力**であり、`M` の側のデータ(`heM`, `hdeg`, `htop`)からは
   出ない(`‖(p:M)‖ = ‖π‖^{p^{k+1}e}` は `c = p^n` の場合しか決めない)。
   ★木の `MinpolyOrbitSplit.norm_algebraMap_mem_zpow_of_uniformizer`(481 行)は
   `hval` を「`‖ϖ‖ = ‖π‖^p`」＋「`F` の値群が `‖ϖ‖^ℤ`」の 2 本に**分ける**だけで、減らさない。
2. ★**`hbreak : ‖τ^{p^k}π − π‖ = ‖π‖^{i+1}` は `i` の定義**なので落ちない(`hti` と対)。
   ★`hjump k` は不等号でこれを含むが、降下が sharp であるためには**等号**が要る。
3. ★**`hstepZ` / `hlayerZ` は `ℤ` の側の跳びの条件**で、`GainedTowerDescent` の仮定である。
4. ★`hjump` / `hstepZ` / `hlayerZ` を `∀ j`(`j > k` も)で取っているのは
   `GainedDescentBridge` からの継承である。★`j > k` では `τ^{p^j} = 1` で `hjump` は自明に成り立ち、
   `t` を `t_{j+1} = p^{j−k}t_{k+1}` と延長すれば `hstepZ`/`hlayerZ` も満たせる(数値検算で確認)。

5. ★★**`ht0 : 0 ≤ t_{j+1}` を `ht1 : 1 ≤ t_{j+1}` に強めた**(★逸脱、正直に記録する)。
   `norm_algHom_sub_le_mul_of_break` が縮小率 `θ = ‖π‖^{t_{j+1}}` に `θ < 1` を要求するためで、
   `t_{j+1} = 0` だと「`σ` は `π` を `‖π‖` しか動かさない」までしか言えず
   `‖σπ − π‖ < ‖π‖`(桁の分離に要る)が出ない。
   ★暴分岐の跳びは `t ≥ 1`(`σ ∈ G_1`)なので、意図した設定では自由に置ける。
6. ★`hσ : ∀ z, σ z = τ^{p^k} z` は不要になった(`σ := τ^{p^k}` を中で作るため)。

## ★★指数の勘定(★`GainedDescentBridge` と同じ例で検算した)

`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`: `p = 3`, `e = 2`, `k = 1`, `t = (t₁,t₂) = (2,8)`, `e_M = p²e = 18`。

* `hnormp`: `‖3‖ = 1/3` ✓、`heM`: `‖3‖ = ‖π‖^{18}` ✓(`e(M/ℚ₃) = φ(27) = 18`、全分岐)
* `ht1`: `2 ≥ 1`, `8 ≥ 1` ✓  `hjump`: `v(τπ−π) = 3 ≥ t₁+1 = 3` ✓、`v(τ³π−π) = 9 ≥ t₂+1 = 9` ✓
* `hlayerZ`: `2·2 = 4 ≤ 3·2 = 6` ✓、`2·8 = 16 ≤ 9·2 = 18` ✓   `hstepZ`: `3·2 = 6 ≤ 8` ✓
* 結論の左 `3^{4/9}` ≤ 右 `∏_{j∈[1,2]} axDecay 3 j = 3^{2/3}` ✓

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `GainedDescentBridge` / `CyclicJumpNorm` / `MinpolyOrbitSplit` / `GainedTowerDescent` は
   **読むだけ**で 1 行も書き換えていない。弱めた版は本ファイルに別名で置いた
   (`..._of_lt` は `hlayer` を `∀ j < k` にしたもの)。
2. ★`hlayer` を `∀ j < k` に弱めた証明は、元の `∀ j` 版に
   `θ' j = if j < k then θ j else 1` を代入するだけで済ませた(元の帰納を書き直していない)。
3. 原典(Serre *Corps Locaux* IV)は付値言語。本ファイルはノルム言語
   (`v_M(z) = m` は `‖z‖ = ‖π‖^m`)。`GainedDescentBridge` と同じ約束である。
-/
open Finset

namespace ABC3.Found.PGC

namespace GainedBridgeSupply

/-! ## §1 抽象核(実数のみ) -/

section RealCore

/-- 抽象核(実数のみ): `0 < b`, `E ≠ 0`, `b^E = c⁻¹` なら `b = c^{−1/E}`。 -/
theorem eq_rpow_neg_one_div_of_pow_eq_inv {b c : ℝ} {E : ℕ} (hb0 : 0 < b) (hc0 : 0 ≤ c)
    (hE : E ≠ 0) (h : b ^ E = c⁻¹) : b = c ^ (-(1 : ℝ) / (E : ℝ)) := by
  have h1 : (b ^ E) ^ ((E : ℝ))⁻¹ = b := Real.pow_rpow_inv_natCast hb0.le hE
  rw [h] at h1
  rw [← h1, ← Real.rpow_neg_one c, ← Real.rpow_mul hc0]
  congr 1

/-- 抽象核(実数のみ): `b^E < 1` なら `b < 1`。 -/
theorem lt_one_of_pow_lt_one {b : ℝ} {E : ℕ} (h : b ^ E < 1) : b < 1 := by
  by_contra hc
  exact absurd h (not_lt.mpr (one_le_pow₀ (not_lt.mp hc)))

end RealCore

/-! ## §2 抽象核(ノルムのみ) —— 塔の得は `j < k` の層しか使わない -/

section NormCore

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- ★★`GainedBridge.norm_iterate_pow_sub_self_le` の `hlayer` を `∀ j < k` に弱めた形。

★証明は `θ' j = if j < k then θ j else 1` を代入するだけ
(`j ≥ k` では `θ' j ^ (p−1) = 1` で `‖(p:L)‖ ≤ 1` が自明に成り立ち、
`hstep` も `θ j ≤ 1` から出る)。★積は `range k` 上なので結論は変わらない。 -/
theorem norm_iterate_pow_sub_self_le_of_lt {p : ℕ} (hp : p.Prime) {τ : L →+ L} {θ : ℕ → ℝ}
    (hθ0 : ∀ j, 0 ≤ θ j) (hθ1 : ∀ j, θ j ≤ 1)
    (hstep : ∀ (j : ℕ) (z : L), ‖(⇑τ)^[p ^ j] z - z‖ ≤ θ j * ‖z‖)
    (k : ℕ) (hlayer : ∀ j, j < k → ‖(p : L)‖ ≤ θ j ^ (p - 1)) (x : L) :
    ‖(⇑τ)^[p ^ k] x - x‖ ≤ (∏ j ∈ Finset.range k, θ j ^ (p - 1)) * ‖τ x - x‖ := by
  classical
  set θ' : ℕ → ℝ := fun j => if j < k then θ j else 1 with hθ'def
  have hlt : ∀ j, j < k → θ' j = θ j := by intro j hj; simp [hθ'def, hj]
  have hge : ∀ j, ¬ j < k → θ' j = 1 := by intro j hj; simp [hθ'def, hj]
  have h0 : ∀ j, 0 ≤ θ' j := by
    intro j; by_cases hj : j < k
    · rw [hlt j hj]; exact hθ0 j
    · rw [hge j hj]; exact zero_le_one
  have h1 : ∀ j, θ' j ≤ 1 := by
    intro j; by_cases hj : j < k
    · rw [hlt j hj]; exact hθ1 j
    · rw [hge j hj]
  have hs : ∀ (j : ℕ) (z : L), ‖(⇑τ)^[p ^ j] z - z‖ ≤ θ' j * ‖z‖ := by
    intro j z; by_cases hj : j < k
    · rw [hlt j hj]; exact hstep j z
    · rw [hge j hj, one_mul]
      refine le_trans (hstep j z) ?_
      nlinarith [norm_nonneg z, hθ1 j]
  have hl : ∀ j, ‖(p : L)‖ ≤ θ' j ^ (p - 1) := by
    intro j; by_cases hj : j < k
    · rw [hlt j hj]; exact hlayer j hj
    · rw [hge j hj, one_pow]; exact IsUltrametricDist.norm_natCast_le_one L p
  have hmain := GainedBridge.norm_iterate_pow_sub_self_le hp h0 h1 hs hl k x
  refine le_trans hmain (le_of_eq ?_)
  congr 1
  exact Finset.prod_congr rfl fun j hj => by rw [hlt j (Finset.mem_range.mp hj)]

end NormCore

/-! ## §3 供給 —— `hlayer` / `hchar` / `hπE` を素元の正規化から出す -/

section Supply

variable {M : Type*} [NormedField M]

/-- 抽象核(ノルムのみ): `b = ‖π‖`, `‖p‖ = b^E`, `c ≤ E` なら `‖p‖ ≤ (b^c)^{p−1}`。

★★これが `hlayer` の供給である。`c = ((p:ℤ)−1)·t_{j+1}`、`E = e_M` と読む。 -/
theorem norm_natCast_le_zpow_pow {p : ℕ} (hp1 : 1 ≤ p) {π : M} {E : ℕ} {c : ℤ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1) (hE : ‖(p : M)‖ = ‖π‖ ^ E)
    (hc : ((p : ℤ) - 1) * c ≤ (E : ℤ)) :
    ‖(p : M)‖ ≤ (‖π‖ ^ c) ^ (p - 1) := by
  have hcast : ((p - 1 : ℕ) : ℤ) = (p : ℤ) - 1 := by
    push_cast [Nat.cast_sub hp1]; ring
  have hterm : (‖π‖ ^ c) ^ (p - 1) = ‖π‖ ^ (((p : ℤ) - 1) * c) := by
    rw [← zpow_natCast (‖π‖ ^ c) (p - 1), ← zpow_mul, hcast, mul_comm]
  rw [hterm, hE, ← zpow_natCast ‖π‖ E]
  exact zpow_le_zpow_right_of_le_one₀ hπ0 hπ1 hc

/-- ★★`hlayer`(層 `j` の sharp な上界)を `hlayerZ`(`ℤ` の側)と
`‖(p:M)‖ = ‖π‖^{e_M}`(★`π` が `M` の素元で `e_M = p^{k+1}e`)から供給する。

★`j ≤ k` に制限しているのが要点で、`GainedDescentBridge` が `hlayer` と `hlayerZ` を
**別々に**取っていた重複はここで消える(`j > k` では `p^{j+1}e > e_M` なので成り立たない)。 -/
theorem norm_natCast_le_zpow_pow_of_layerZ {p e k : ℕ} (hp1 : 1 ≤ p) {π : M} {t : ℕ → ℤ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hlayerZ : ∀ j, ((p : ℤ) - 1) * t (j + 1) ≤ (p : ℤ) ^ (j + 1) * (e : ℤ)) :
    ∀ j, j ≤ k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1) := by
  intro j hj
  refine norm_natCast_le_zpow_pow hp1 hπ0 hπ1 heM (le_trans (hlayerZ j) ?_)
  have hpow : (p : ℤ) ^ (j + 1) ≤ (p : ℤ) ^ (k + 1) :=
    pow_le_pow_right₀ (by exact_mod_cast hp1) (by omega)
  have he0 : (0 : ℤ) ≤ (e : ℤ) := Int.natCast_nonneg e
  calc (p : ℤ) ^ (j + 1) * (e : ℤ) ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) :=
        mul_le_mul_of_nonneg_right hpow he0
    _ = ((p ^ (k + 1) * e : ℕ) : ℤ) := by push_cast; ring

end Supply

/-! ## §4 弱めた `hlayer` で橋を張り直す(`GainedDescentBridge` §3〜§5 の写し) -/

section Bridge

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- `GainedBridge.norm_sub_digit_zero_le_gain` の `hlayer` を `∀ j < k` に弱めた形。 -/
theorem norm_sub_digit_zero_le_gain_of_lt {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) {θ : ℕ → ℝ}
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z)
    (hθ0 : ∀ j, 0 ≤ θ j) (hθ1 : ∀ j, θ j ≤ 1)
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ θ j * ‖z‖)
    (hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ θ j ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(i : ℤ)) * (∏ j ∈ Finset.range k, θ j ^ (p - 1))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  have heq := norm_sub_digit_zero_eq_zpow_mul (a := a) σ hp.pos hπ0 hπ1 hi hbreak hchar hval
  have hchain := norm_iterate_pow_sub_self_le_of_lt hp hθ0 hθ1 hstep k hlayer (digitSum p π a)
  rw [heq, hσ (digitSum p π a)]
  have hz0 : (0 : ℝ) ≤ ‖π‖ ^ (-(i : ℤ)) := le_of_lt (zpow_pos hπ0 _)
  calc ‖π‖ ^ (-(i : ℤ)) * ‖(⇑τ)^[p ^ k] (digitSum p π a) - digitSum p π a‖
      ≤ ‖π‖ ^ (-(i : ℤ)) * ((∏ j ∈ Finset.range k, θ j ^ (p - 1))
          * ‖τ (digitSum p π a) - digitSum p π a‖) := mul_le_mul_of_nonneg_left hchain hz0
    _ = _ := by ring

/-- `GainedBridge.norm_sub_digit_zero_le_topDefect` の `hlayer` を `∀ j < k` に弱めた形。 -/
theorem norm_sub_digit_zero_le_topDefect_of_lt {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(t (k + 1) - ((p : ℤ) - 1) * GainedDescent.jumpSum t k))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  have h := norm_sub_digit_zero_le_gain_of_lt (a := a) (θ := fun j => ‖π‖ ^ (t (j + 1)))
    hp τ σ hσ (fun j => le_of_lt (zpow_pos hπ0 _))
    (fun j => zpow_le_one₀ hπ0 (le_of_lt hπ1) (ht0 j)) hstep hlayer hπ0 hπ1 hi hbreak hchar hval
  refine le_trans h (le_of_eq ?_)
  rw [GainedBridge.prod_theta_eq hπ0 hp.one_lt.le t k, ← zpow_add₀ (ne_of_gt hπ0), hti]
  congr 1
  ring

/-- `GainedBridge.norm_sub_digit_zero_le_gainedLoss` の `hlayer` を `∀ j < k` に弱めた形。 -/
theorem norm_sub_digit_zero_le_gainedLoss_of_lt {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(GainedDescent.gainedLoss (p : ℤ) t (k + 1)))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  refine le_trans (norm_sub_digit_zero_le_topDefect_of_lt hp τ σ t hσ ht0 hstep hlayer hπ0 hπ1 hi
    hti hbreak hchar hval) (mul_le_mul_of_nonneg_right ?_ (norm_nonneg _))
  refine zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) ?_
  have := GainedBridge.topDefect_le_gainedLoss (p : ℤ) t k
  omega

/-- `GainedBridge.exists_norm_sub_algebraMap_le_gainedLoss` の `hlayer` を `∀ j < k` に弱めた形。 -/
theorem exists_norm_sub_algebraMap_le_gainedLoss_of_lt {p i k : ℕ} (hp : p.Prime) {π : M}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ ‖π‖ ^ (-(GainedDescent.gainedLoss (p : ℤ) t (k + 1))) * ‖τ x - x‖ := by
  have hint : IsIntegral F π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp.ne_zero hdeg.symm
  obtain ⟨a, rfl⟩ := exists_digitSum_of_adjoin_eq_top hint hdeg htop x
  exact ⟨a 0, norm_sub_digit_zero_le_gainedLoss_of_lt hp τ σ t hσ ht0 hstep hlayer hπ0 hπ1 hi hti
    hbreak hchar hval⟩

end Bridge

/-! ## §5 出口 —— 落とした仮説は 5 本 -/

section Exit

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★★**`GainedDescentBridge` の出口から仮説を 5 本落とした形**。

落ちたのは `hchar` / `hlayer` / `hπE` / `hπ0` / `hπ1` の 5 本で、
代わりに入るのは★**正規化 2 本**だけである:

* `hnormp : ‖(p:M)‖ = p⁻¹` —— `M` のノルムが `ℚ_p` のノルムを延長する(正規化の約束)、
* `heM : ‖(p:M)‖ = ‖π‖^{p^{k+1}e}` —— ★`π` が `M` の素元で `e_M = p^{k+1}e`(塔が全分岐)。

★★これで `hlayer`(ノルム版)と `hlayerZ`(`ℤ` 版)の**重複が消えた**
(`GainedDescentBridge` の docstring「閉じていないもの 2」)。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ) (he : 0 < e)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hstepZ : ∀ j, (p : ℤ) * t (j + 1) ≤ t (j + 2))
    (hlayerZ : ∀ j, ((p : ℤ) - 1) * t (j + 1) ≤ (p : ℤ) ^ (j + 1) * (e : ℤ))
    (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hplt : ‖(p : M)‖ < 1 := by
    rw [hnormp, inv_lt_one_iff₀]; exact Or.inr hpR
  have hπ0 : 0 < ‖π‖ := by
    rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hEne] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  have hπ1 : ‖π‖ < 1 := by
    refine lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  have hπE : ‖π‖ = (p : ℝ) ^ (-(1 : ℝ) / ((p ^ (k + 1) * e : ℕ) : ℝ)) :=
    eq_rpow_neg_one_div_of_pow_eq_inv hπ0 (le_of_lt hpR0) hEne hπpow
  have hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1 :=
    fun j hj0 hjp => norm_natCast_eq_one_of_lt_prime hp' hplt hj0 hjp
  have hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1) :=
    fun j hj => norm_natCast_le_zpow_pow_of_layerZ hp'.one_lt.le hπ0 (le_of_lt hπ1) heM hlayerZ j
      (le_of_lt hj)
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_gainedLoss_of_lt hp' τ σ t hσ ht0
    hstep hlayer hπ0 hπ1 hi hti hbreak hchar hval hdeg htop x
  refine ⟨y, le_trans hy (mul_le_mul_of_nonneg_right ?_ (norm_nonneg _))⟩
  have hΛ0 : 0 ≤ GainedDescent.gainedLoss (p : ℤ) t (k + 1) :=
    GainedBridge.gainedLoss_nonneg (by simpa using ht0 0) hstepZ k
  have hJ : ((GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat : ℤ)
      = GainedDescent.gainedLoss (p : ℤ) t (k + 1) := Int.toNat_of_nonneg hΛ0
  have hmain := GainedDescent.rpow_le_prod_axDecay (p := p) (e := e) (k := k + 1)
    (J := (GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat) (t := t) he ht0 hstepZ hlayerZ
    (le_of_eq hJ)
  have hcast : ((GainedDescent.gainedLoss (p : ℤ) t (k + 1) : ℤ) : ℝ)
      = (((GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat : ℕ) : ℝ) := by
    exact_mod_cast hJ.symm
  rw [GainedBridge.zpow_eq_rpow_div hπE, hcast]
  exact hmain

end Exit

/-! ## §6 `hstep` の供給 —— 生成元 `π` の上の条件(＝分岐群の定義)から全元の縮小率へ -/

section Step

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★**`σ` が `π` を `θ·‖π‖` しか動かさないなら、`M` の全元を `θ` 倍しか動かさない**。

★★これが「下付き分岐群 `Γ_t` の定義は生成元 `π` の上の条件で書けば十分」という主張であり、
`GainedDescentBridge` の `hstep`(`∀ z ∈ M` の条件)を
★**`π` だけの条件**に落とす部品である。

★中身は木に既に在った `CyclicJumpNorm.norm_algHom_sub_mul_norm_eq` の**等式**

  `‖σx − x‖ · ‖π‖ = ‖σπ − π‖ · ‖x − a₀‖`

と `norm_digitSum_sub_digit_zero_le`(`‖x − a₀‖ ≤ ‖x − c‖`、`c = 0` を取る)だけである。 -/
theorem norm_algHom_sub_le_mul_of_break {p : ℕ} (hp : 0 < p) {π : M} (σ : M →ₐ[F] M) {θ : ℝ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hθ1 : θ < 1)
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hbrk : ‖σ π - π‖ ≤ θ * ‖π‖)
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (z : M) :
    ‖σ z - z‖ ≤ θ * ‖z‖ := by
  have hint : IsIntegral F π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    omega
  have hu : ‖σ π - π‖ < ‖π‖ := by
    refine lt_of_le_of_lt hbrk ?_
    calc θ * ‖π‖ < 1 * ‖π‖ := mul_lt_mul_of_pos_right hθ1 hπ0
      _ = ‖π‖ := one_mul _
  obtain ⟨a, rfl⟩ := exists_digitSum_of_adjoin_eq_top hint hdeg htop z
  have heq := norm_algHom_sub_mul_norm_eq (a := a) σ hp hπ1 hu hchar hval
  have hle : ‖digitSum p π a - algebraMap F M (a 0)‖ ≤ ‖digitSum p π a‖ := by
    have h0 := norm_digitSum_sub_digit_zero_le hp hπ0 hπ1 hval a 0
    simpa using h0
  have hθ0 : (0 : ℝ) ≤ θ * ‖π‖ := le_trans (norm_nonneg _) hbrk
  refine le_of_mul_le_mul_right ?_ hπ0
  calc ‖σ (digitSum p π a) - digitSum p π a‖ * ‖π‖
      = ‖σ π - π‖ * ‖digitSum p π a - algebraMap F M (a 0)‖ := heq
    _ ≤ (θ * ‖π‖) * ‖digitSum p π a‖ :=
        mul_le_mul hbrk hle (norm_nonneg _) hθ0
    _ = θ * ‖digitSum p π a‖ * ‖π‖ := by ring

end Step

/-! ## §7 出口 2 —— `hstep` も生成元の条件に落とした形 -/

section Exit2

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★★★**`GainedDescentBridge` の出口から仮説を 6 本落とした形**。

落ちたのは `hchar` / `hlayer` / `hπE` / `hπ0` / `hπ1` / `hstep` の 6 本。
★★とくに `hstep`(`∀ j`, `∀ z ∈ M` の縮小率)が

```
hjump : ∀ j, ‖τ^{p^j} π − π‖ ≤ ‖π‖^{t_{j+1}+1}
```

すなわち★**生成元 `π` の上の条件だけ**になった。これは下付き分岐群の定義
`τ^{p^j} ∈ Γ_{t_{j+1}} ⟺ v_M(τ^{p^j}π − π) ≥ t_{j+1}+1` そのものである
(`𝒪_M = 𝒪_F[π]` なので生成元で測れば十分、というのが `norm_algHom_sub_le_mul_of_break`)。

★残る仮説は `hstepZ` / `hlayerZ`(★どちらも `ℤ` の側の跳びの条件)、
`hbreak` / `hti`(★`i` の定義)、`hval`(値群 = 全分岐)、`hnormp` / `heM`(正規化)、
`hdeg` / `htop`(`M = F(π)` が `p` 次)である。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M}
    (τ : M →ₐ[F] M) (t : ℕ → ℤ) (he : 0 < e) (ht1 : ∀ j, 1 ≤ t (j + 1))
    (hjump : ∀ j, ‖(τ ^ p ^ j) π - π‖ ≤ ‖π‖ ^ (t (j + 1) + 1))
    (hstepZ : ∀ j, (p : ℤ) * t (j + 1) ≤ t (j + 2))
    (hlayerZ : ∀ j, ((p : ℤ) - 1) * t (j + 1) ≤ (p : ℤ) ^ (j + 1) * (e : ℤ))
    (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖(τ ^ p ^ k) π - π‖ = ‖π‖ ^ (i + 1))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hplt : ‖(p : M)‖ < 1 := by rw [hnormp, inv_lt_one_iff₀]; exact Or.inr hpR
  have hπ0 : 0 < ‖π‖ := by
    rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hEne] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  have hπ1 : ‖π‖ < 1 := by
    refine lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  have hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1 :=
    fun j hj0 hjp => norm_natCast_eq_one_of_lt_prime hp' hplt hj0 hjp
  set τ' : M →+ M := AddMonoidHom.mk' (⇑τ) (fun a b => map_add τ a b) with hτ'
  have hcoe : ⇑τ' = ⇑τ := rfl
  have hpowcoe : ∀ j : ℕ, (⇑τ')^[p ^ j] = ⇑(τ ^ p ^ j) := by
    intro j; rw [hcoe, AlgHom.coe_pow]
  have hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ')^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖ := by
    intro j z
    rw [hpowcoe j]
    refine norm_algHom_sub_le_mul_of_break hp'.pos (τ ^ p ^ j) hπ0 hπ1
      (zpow_lt_one₀ hπ0 hπ1 (lt_of_lt_of_le zero_lt_one (ht1 j))) hchar hval ?_ hdeg htop z
    refine le_trans (hjump j) (le_of_eq ?_)
    rw [zpow_add_one₀ (ne_of_gt hπ0)]
  have hσ : ∀ z : M, (τ ^ p ^ k) z = (⇑τ')^[p ^ k] z := by
    intro z; rw [hpowcoe k]
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer
    (F := F) (M := M) (p := p) (e := e) (i := i) (k := k) (π := π) τ' (τ ^ p ^ k) t he hσ
    (fun j => le_trans zero_le_one (ht1 j)) hstep hstepZ hlayerZ hi hti hnormp heM hbreak hval
    hdeg htop x
  exact ⟨y, hy⟩

end Exit2

/-! ## §8 `.src`(原典の対応箇所) -/

def norm_iterate_pow_sub_self_le_of_lt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_natCast_le_zpow_pow_of_layerZ.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_algHom_sub_le_mul_of_break.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §9 使っている公理の一覧 -/

#print axioms eq_rpow_neg_one_div_of_pow_eq_inv
#print axioms lt_one_of_pow_lt_one
#print axioms norm_iterate_pow_sub_self_le_of_lt
#print axioms norm_natCast_le_zpow_pow
#print axioms norm_natCast_le_zpow_pow_of_layerZ
#print axioms norm_sub_digit_zero_le_gain_of_lt
#print axioms norm_sub_digit_zero_le_topDefect_of_lt
#print axioms norm_sub_digit_zero_le_gainedLoss_of_lt
#print axioms exists_norm_sub_algebraMap_le_gainedLoss_of_lt
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer
#print axioms norm_algHom_sub_le_mul_of_break
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps

end GainedBridgeSupply

end ABC3.Found.PGC
