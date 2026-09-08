import ABC3.Found.PGC.PureStepSetup

/-!
# [pGC] ★★★★★★`k ≥ 1` の塔 —— `τ : M →+ M` の正しい設定で `hstep` を供給する

## ★★★出発点 —— 直前の波が見つけた空虚と、その正確な境界

`Found/PGC/PureStepSetup.lean` の `conj_pure_absurd_of_one_le` / `break_absurd_of_one_le` は
**`τ : M →ₐ[F] M`**(＝`F`-線形)を取る出口 **6 本**が `1 ≤ k` で空虚だと示した:

> `hdeg`(`(minpoly F π).natDegree = p`)＋ `htop`(`Algebra.adjoin F {π} = ⊤`)
> ⟹ `[M : F] = p` ⟹ `orderOf τ ∣ p` ⟹ `τ ^ p = 1` ⟹ `k ≥ 1` で `τ ^ (p^k) = 1`
> ⟹ `hbreak : ‖(τ^{p^k})π − π‖ = ‖π‖^{i+1} > 0` と矛盾。

★★**本ファイルはその境界の**「向こう側」**を埋める。**
`GainedJumpSeq.exists_norm_sub_algebraMap_le_prod_axDecay_lt` は
★**`τ : M →+ M`(加法だけ)**を取り、`F`-線形なのは `σ = τ^{p^k}` の方だけである。
この設定では上の推論は**第 1 段で止まる**(加法的な `τ` の「位数」を `[M:F]` は縛らない)。

★★★**空虚性の分かれ目を形式化した**のが §5 である:

| | `F`-線形版(空虚) | ★本ファイル(塔) |
|---|---|---|
| `τ` に効く体 | `F`(`[M:F] = p`) | ★`E 0`(`[M : E 0] = p^{k+1}`) |
| 出る関係式 | `τ ^ p = 1` | ★`τ ^ (p^{k+1}) = 1`(`iterate_pow_succ_eq_self`) |
| `τ^{p^k} ≠ 1` と両立するか | ★**しない**(`k ≥ 1`) | ★★**する**(`pow_ne_one_of_orderOf_eq`) |

★`pow_ne_one_of_orderOf_eq` は**純群論**(`orderOf τ = p^{k+1}` なら `τ^{p^k} ≠ 1`)で、
`orderOf_ofAdd_one_zmod`(`Multiplicative (ZMod n)` の生成元の位数は `n`)により
★**その位数の型が実在する**ことまで機械検算した。★これで「`k ≥ 1` は空虚ではない」の
群論的な芯は閉じている(体・ノルムまで込めた模型の構成は §6 の「閉じていないもの」に書く)。

## ★★★抽象核 —— 「桁展開さえあれば体は要らない」(§1)

`hstep`(層 `j` の縮小率)は結局

  `‖σ z − z‖ ≤ ‖π‖^c · ‖z‖`   (`σ = τ^{p^j}`、`c = t_{j+1}`)

である。既存の供給 `GainedBridgeSupply.norm_algHom_sub_le_mul_of_break` は
`σ : M →ₐ[F] M` と `hchar`(`0 < j < p` で `‖(j:M)‖ = 1`)と `[M:F] = p` を要求した。
★**層 `j` で要る体は `F_j` であって `F_k = F` ではなく、`[M : F_j] = p^{k+1−j}` は素数でない**ので
`hchar` は**偽**になる(`j = p` で `‖p‖ < 1`)。

★★そこで**等式を捨てて不等式にする**と `hchar` も素数性も要らなくなる:

```
‖σ z − z‖ = ‖Σ_m a_m ((σπ)^m − π^m)‖ ≤ max_m ‖a_m‖·‖σπ − π‖·‖π‖^{m−1}
          ≤ ‖π‖^c · max_m ‖a_m π^m‖ = ‖π‖^c · ‖z‖
```

使うのは `(σπ)^m − π^m = (Σ_r (σπ)^r π^{m−1−r})(σπ − π)`(`geom_sum₂_mul`)と
超距離だけである。★★結果として抽象核 `norm_ringHom_sub_le_zpow_mul` は
**体も Galois も分岐も付値も現れない**:

* `σ : M →+* M`(環準同型。★`F`-線形性は不要)、
* `z = Σ_{m<n} a_m π^m` で `σ a_m = a_m`、`‖a_m π^m‖ ≤ ‖z‖`(桁の分離)、
* `‖σπ − π‖ ≤ ‖π‖^{c+1}`、`0 ≤ c`、`0 < ‖π‖ < 1`

だけから `‖σ z − z‖ ≤ ‖π‖^c ‖z‖` が出る。★`n` は任意で、素数である必要がない。

★★★**#59/#69 の回避**: 抽象核が要求するのは「`σ` が固定する係数での桁展開」だけで、
★**中間体を `IntermediateField` として作る必要が無い**。
`WildDepthDescent` が `MulAction.stabilizer` の指数で次数を測って #59 を避けた
(`lean-idioms` #296)のと同じ趣旨で、本ファイルは
★**塔を「体の族 `E : ℕ → Type*`」として仮説に置く**ことで層をまたぐ `rfl` を 1 度も起こさない。
実測: 本ファイルは 1 回目の `leanfile` で §1–§3 が、2 回目で §4 が通り、#59/#69 に**一度も当たっていない**。

## ★★出口(§4)

| 定理 | 落ちた仮説 |
|---|---|
| `exists_norm_sub_algebraMap_le_prod_axDecay_of_layers` | `hstep`(∀ z ∈ M の縮小率)→ 層ごとの桁展開 `hexp` と素元 1 点の `hbrk` |
| `exists_norm_sub_algebraMap_le_prod_axDecay_of_tower` | `hexp` → 体の塔 `E j`(`[M : E j] = p^{k+1−j}`、全分岐、`s j` が `E j` を固定) |
| ★`exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_jumps` | ★さらに `ℤ` 側 3 本(`ht0`/`hstepZ`/`hlayerZ`)が `hu`＋`hbnd` に落ちる |

★★どれも `τ : M →+ M` のままで、`k ≥ 1` を許す。★これが Ax–Sen–Tate の実質である。

## ★★★数値の検算(§5 `Numeric`) —— `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` で全部の仮説が成り立つ

`p = 3`, `k = 1`, `M = ℚ₃(ζ₂₇)`, `E 0 = ℚ₃(ζ₃)`, `F = E 1 = ℚ₃(ζ₉)`, `π = ζ₂₇ − 1`。

* `e(M/ℚ₃) = 18`, `heM` は `v_M(3) = p^{k+1}·e = 9·e` なので ★`e = 2`
  ——これは**塔の底** `E 0 = ℚ₃(ζ₃)` の絶対分岐指数 `e = 2` に一致する(`[M:E 0] = 18/2 = 9 = p^{k+1}` ✓)。
* `τ : ζ₂₇ ↦ ζ₂₇^4` は `Gal(M/E 0)`(位数 9)の生成元。`τπ − π = ζ(ζ^3−1)` で
  `v_M(τπ−π) = v_M(ζ₉−1) = 3` ⟹ `u₀ + 1 = 3`、★`u₀ = 2`。
* `σ = τ^3 : ζ ↦ ζ^{10}`、`σπ − π = ζ(ζ^9−1) = ζ(ζ₃−1)` で `v_M(ζ₃−1) = 9`
  ⟹ `i + 1 = 9`、★`i = u₁ = 8`。

★★この `(u₀,u₁) = (2,8)` は `GainedJumpSeq.Numeric.seq_zeta27` が別ルートで置いた値と**一致**する
(本ファイルは分岐理論から手計算し、向こうは列の構成から取った)。
`tz` はこのデータから作った列で、`tz_nonneg` / `tz_step` / `tz_layer` / `tz_top` が
★`ht0` / `hstepZ` / `hlayerZ` / `hti` を**すべての `j` について**満たすことを機械検算する
(`t₁ = 2`, `t₂ = 8`, `t_{j+2} = 8·3^j`)。
★とくに `hlayerZ` の一番きつい所は `j = 1` の `2·8 = 16 ≤ 3²·2 = 18` で、**余裕は 2 しかない**。

## ★★★在庫の測定(コマンドつき)

| 部品 | 在庫 | 測定 |
|---|---|---|
| `geom_sum₂_mul` | mathlib に在った | `#check @geom_sum₂_mul` ⟹ `(∑ i ∈ range n, x^i*y^(n-1-i))*(x-y) = x^n - y^n` |
| `IsUltrametricDist.nnnorm_sum_le_of_forall_le` | mathlib に在った | `#check`。★**`norm_` 版は無い**(`Unknown constant IsUltrametricDist.norm_sum_le_of_forall_le`)ので §1 で `lift C to NNReal` で作った |
| 桁の分離 `nnnorm_sum_digit_eq_sup` | ★木に在り、★**`n` が素数である必要が無い**(仮定は `hval` だけ) | `CyclicJumpNorm.lean:312`。★これが本ファイルの鍵 |
| `exists_digitSum_of_adjoin_eq_top` | ★木に在り、★`p` は `natDegree` の名前でしかない | `MinpolyOrbitSplit.lean:528` |
| `ZMod.addOrderOf_one` / `orderOf_ofAdd_eq_addOrderOf` | mathlib に在った | `#check`(位数ちょうど `n` の群元の実在) |
| `Nat.pos_pow_of_pos` | ★mathlib に**無い** | `Unknown constant `Nat.pos_pow_of_pos``。`pow_pos` を使う |
| `PureStepSetup.algHom_pow_apply_eq_self_of_natDegree` | 木に在った | `PureStepSetup.lean:188`。§5 で `M →+* M` から `M →ₐ[E] M` を作って流用 |

## ★★★閉じていないもの(★正確に)

1. ★**体・ノルムまで込めた `k ≥ 1` の模型は構成していない。**
   §5 が閉じたのは「位数の型が実在する」という**群論の部分**だけである
   (`orderOf_ofAdd_one_zmod`)。`hnormp` / `heM` / `hval` / `hdegj` / `htopj` を
   同時に満たす具体的な `M`(例: `ℚ₃(ζ₂₇)`)を Lean で作るには
   `Found/PGC` に円分体の完備化と分岐の計算が要る。★本ファイルはそれを**仮説のまま**置く。
   ★ただし `ℤ` 側の数値は `ℚ₃(ζ₂₇)` の実データで機械検算した(§5 `Numeric`)。
2. ★`hjump`(層 `j` の跳び `‖s j π − π‖ ≤ ‖π‖^{u_j+1}`)は仮説のままである。
   これは下付き分岐群の定義そのもので、`GainedJumpSeq` と同じ状況。
3. ★`hu : ∀ m ≤ k, p^m ≤ u m`(Hasse–Arf の中身)も `GainedJumpSeq` から引き継いだ仮説のままである。
4. `E : ℕ → Type*` の族は「塔がある」ことを**仮定**する形で、
   `M/F₀` が巡回 `p^{k+1}` 次全分岐であることから族を**作る**部分は書いていない。
   ★これは中間体の構成なので #59/#69 の危険がある所であり、本ファイルは意図的に避けた。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `GainedJumpSeq` / `GainedJumpFree` / `PureStepSetup` は**読むだけ**で 1 行も書き換えていない。
2. ★§1 の抽象核は `GainedBridgeSupply.norm_algHom_sub_le_mul_of_break` の**等式**を
   捨てて不等式にした。結論(`hstep` の形)は同じで、仮説から `hchar` と `p` の素数性が落ちる。
   ★原典(Serre *Corps Locaux* IV / Tate 1967 §3.2)は付値言語。本ファイルはノルム言語
   (`v_M(z) = m` は `‖z‖ = ‖π‖^m`)で、`GainedDescentBridge` と同じ約束である。
3. ★塔の層 `j` の体を `IntermediateField F₀ M` ではなく**独立した型 `E j`** として置いた。
   原典は中間体だが、`Algebra (E j) M` と `hdegj`/`htopj`/`hvalj` が同じ内容を表す。
   ★意図は #59(中間体 2 層をまたぐ `rfl` が kernel を止める)の回避である。
-/
open Finset

namespace ABC3.Found.PGC

namespace GainedTowerStep

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-! ## §1 抽象核 -/

section Core

/-- 超距離の有限和の上界(実数版)。
★mathlib にあるのは `IsUltrametricDist.nnnorm_sum_le_of_forall_le`(NNReal 版)だけで、
`IsUltrametricDist.norm_sum_le_of_forall_le` は無い(`Unknown constant`)。`lift C to NNReal` で作る。 -/
theorem norm_sum_le_of_forall_le {ι : Type*} {s : Finset ι} {f : ι → M} {C : ℝ} (hC : 0 ≤ C)
    (h : ∀ i ∈ s, ‖f i‖ ≤ C) : ‖∑ i ∈ s, f i‖ ≤ C := by
  lift C to NNReal using hC with C'
  exact_mod_cast IsUltrametricDist.nnnorm_sum_le_of_forall_le
    (fun i hi => by exact_mod_cast h i hi)

/-- `‖w‖ ≤ ‖π‖` なら `‖w^m − π^m‖ ≤ ‖π‖^{m−1}·‖w − π‖`。

`geom_sum₂_mul`(`(∑_r x^r y^{m−1−r})(x−y) = x^m − y^m`)と超距離だけ。
★ここに `m < p` や `‖(m:M)‖ = 1` は一切入らないのが本ファイルの鍵である。 -/
theorem norm_pow_sub_pow_le {π w : M} (hπ0 : 0 < ‖π‖) (hw : ‖w‖ ≤ ‖π‖) (m : ℕ) :
    ‖w ^ m - π ^ m‖ ≤ ‖π‖ ^ (m - 1) * ‖w - π‖ := by
  rw [← geom_sum₂_mul w π m, norm_mul]
  refine mul_le_mul_of_nonneg_right ?_ (norm_nonneg _)
  refine norm_sum_le_of_forall_le (pow_nonneg (le_of_lt hπ0) _) ?_
  intro i hi
  have him : i < m := Finset.mem_range.mp hi
  rw [norm_mul, norm_pow, norm_pow]
  calc ‖w‖ ^ i * ‖π‖ ^ (m - 1 - i) ≤ ‖π‖ ^ i * ‖π‖ ^ (m - 1 - i) := by gcongr
    _ = ‖π‖ ^ (i + (m - 1 - i)) := (pow_add _ _ _).symm
    _ = ‖π‖ ^ (m - 1) := by congr 1; omega

/-- ★★★抽象核 —— 桁展開さえあれば、`σ` が素元を動かす量だけで全元の縮小率が出る。

★分岐・付値・Galois・体の語彙が 1 語も出ない。入力は

* `σ : M →+* M`(★環準同型だけ。`F`-線形性は不要)、
* `z = ∑_{m<n} a_m π^m` で `σ a_m = a_m` かつ `‖a_m π^m‖ ≤ ‖z‖`(★桁の分離)、
* `‖σπ − π‖ ≤ ‖π‖^{c+1}`、`0 ≤ c`、`0 < ‖π‖ < 1`

だけで、`n` は任意(★素数でなくてよい)。

★`GainedBridgeSupply.norm_algHom_sub_le_mul_of_break` は同じ結論を
等式(`‖σx−x‖·‖π‖ = ‖σπ−π‖·‖x−a₀‖`)経由で出すため
`hchar`(`0<j<p` で `‖(j:M)‖ = 1`)と `[M:F] = p` を要求した。
★★層 `j` では `[M : F_j] = p^{k+1−j}` で素数でないので `hchar` は偽になる
(`j = p` で `‖p‖ < 1`)。本補題は等式を捨てて不等式にすることでその 2 本を落とす。 -/
theorem norm_ringHom_sub_le_zpow_mul {n : ℕ} {π : M} {c : ℤ} (σ : M →+* M)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hc0 : 0 ≤ c)
    (hbrk : ‖σ π - π‖ ≤ ‖π‖ ^ (c + 1))
    {z : M} {a : ℕ → M} (hfix : ∀ m, m < n → σ (a m) = a m)
    (hz : z = ∑ m ∈ Finset.range n, a m * π ^ m)
    (hle : ∀ m, m < n → ‖a m * π ^ m‖ ≤ ‖z‖) :
    ‖σ z - z‖ ≤ ‖π‖ ^ c * ‖z‖ := by
  have hπne : ‖π‖ ≠ 0 := ne_of_gt hπ0
  have hbrk1 : ‖σ π - π‖ ≤ ‖π‖ := by
    refine le_trans hbrk ?_
    calc ‖π‖ ^ (c + 1) ≤ ‖π‖ ^ (1 : ℤ) :=
          zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (by omega)
      _ = ‖π‖ := zpow_one _
  have hσπ : ‖σ π‖ ≤ ‖π‖ := by
    have h2 := IsUltrametricDist.norm_add_le_max (σ π - π) π
    rw [sub_add_cancel] at h2
    exact le_trans h2 (max_le hbrk1 le_rfl)
  have hshift : ∀ m : ℕ, 1 ≤ m →
      ‖π‖ ^ (m - 1) * ‖π‖ ^ (c + 1) = ‖π‖ ^ c * ‖π‖ ^ m := by
    intro m hm
    rw [← zpow_natCast ‖π‖ (m - 1), ← zpow_natCast ‖π‖ m, ← zpow_add₀ hπne, ← zpow_add₀ hπne]
    congr 1
    omega
  have hmap : σ z = ∑ m ∈ Finset.range n, a m * (σ π) ^ m := by
    rw [hz, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_mul, map_pow, hfix m (Finset.mem_range.mp hm)]
  have hdiff : σ z - z = ∑ m ∈ Finset.range n, a m * ((σ π) ^ m - π ^ m) := by
    rw [hmap]
    conv_lhs => rw [hz]
    rw [← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun m _ => (mul_sub _ _ _).symm
  rw [hdiff]
  refine norm_sum_le_of_forall_le
    (mul_nonneg (le_of_lt (zpow_pos hπ0 _)) (norm_nonneg _)) ?_
  intro m hm
  have hmn : m < n := Finset.mem_range.mp hm
  rcases Nat.eq_zero_or_pos m with hm0 | hm1
  · subst hm0
    simp only [pow_zero, sub_self, mul_zero, norm_zero]
    exact mul_nonneg (le_of_lt (zpow_pos hπ0 _)) (norm_nonneg _)
  · calc ‖a m * ((σ π) ^ m - π ^ m)‖
        = ‖a m‖ * ‖(σ π) ^ m - π ^ m‖ := norm_mul _ _
      _ ≤ ‖a m‖ * (‖π‖ ^ (m - 1) * ‖σ π - π‖) :=
          mul_le_mul_of_nonneg_left (norm_pow_sub_pow_le hπ0 hσπ m) (norm_nonneg _)
      _ ≤ ‖a m‖ * (‖π‖ ^ (m - 1) * ‖π‖ ^ (c + 1)) := by gcongr
      _ = ‖π‖ ^ c * (‖a m‖ * ‖π‖ ^ m) := by rw [hshift m hm1]; ring
      _ = ‖π‖ ^ c * ‖a m * π ^ m‖ := by rw [norm_mul, norm_pow]
      _ ≤ ‖π‖ ^ c * ‖z‖ :=
          mul_le_mul_of_nonneg_left (hle m hmn) (le_of_lt (zpow_pos hπ0 _))

end Core

/-! ## §2 桁展開の供給 -/

section FieldBridge

variable {E : Type*} [Field E] [Algebra E M]

/-- ★抽象核の `hexp` を体から供給する。

`M = E(π)`、`[M:E] = n`、`hval`(全分岐 —— `E^×` の値群が `‖π‖^{nℤ}`)、
`σ` が `E` を点ごとに固定 ⟹ 桁展開が取れて桁が分離する。
★中身は木にあった `exists_digitSum_of_adjoin_eq_top`(`MinpolyOrbitSplit.lean:528`)と
`nnnorm_sum_digit_eq_sup`(`CyclicJumpNorm.lean:312`)だけで、どちらも `n` の素数性を使わない。 -/
theorem exists_expansion_of_adjoin {n : ℕ} {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : 0 < n)
    (hval : ∀ d : E, d ≠ 0 → ∃ m : ℤ, ‖algebraMap E M d‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hdeg : (minpoly E π).natDegree = n)
    (htop : Algebra.adjoin E ({π} : Set M) = ⊤)
    (σ : M →+* M) (hσ : ∀ d : E, σ (algebraMap E M d) = algebraMap E M d) (z : M) :
    ∃ a : ℕ → M, (∀ m, m < n → σ (a m) = a m) ∧
      z = ∑ m ∈ Finset.range n, a m * π ^ m ∧
      ∀ m, m < n → ‖a m * π ^ m‖ ≤ ‖z‖ := by
  have hint : IsIntegral E π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    omega
  obtain ⟨b, hb⟩ := exists_digitSum_of_adjoin_eq_top hint hdeg htop z
  refine ⟨fun m => algebraMap E M (b m), fun m _ => hσ _, by rw [hb, digitSum], ?_⟩
  intro m hm
  have hsup : ‖z‖₊ = (Finset.range n).sup fun j => ‖algebraMap E M (b j) * π ^ j‖₊ := by
    rw [hb, digitSum]
    exact nnnorm_sum_digit_eq_sup hπ0 hπ1 hval b fun j hj => Finset.mem_range.mp hj
  have hx := Finset.le_sup (f := fun j => ‖algebraMap E M (b j) * π ^ j‖₊)
    (Finset.mem_range.mpr hm)
  rw [← hsup] at hx
  exact_mod_cast hx

end FieldBridge

omit [IsUltrametricDist M] in
/-- `hnormp`(`‖p‖ = 1/p`)と `heM`(`‖p‖ = ‖π‖^K`, `K ≠ 0`)から `0 < ‖π‖ < 1`。 -/
theorem norm_pos_and_lt_one {p K : ℕ} {π : M} (hp1 : 1 < p) (hK : K ≠ 0)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹) (heM : ‖(p : M)‖ = ‖π‖ ^ K) :
    0 < ‖π‖ ∧ ‖π‖ < 1 := by
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp1
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hπpow : ‖π‖ ^ K = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  refine ⟨?_, ?_⟩
  · rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hK] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  · refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := K) ?_
    rw [hπpow, inv_lt_one_iff₀]
    exact Or.inr hpR

/-! ## §3 層の供給 -/

section Layers

/-- ★★これが持ち場の本体 —— `GainedJumpSeq...lt` の `hstep`(層ごとの縮小率)の供給。

`τ : M →+ M` は加法的なだけでよい。★層 `j` で要るのは
「`s j`(★`τ^{p^j}` と同じ写像である環準同型)の固定する係数での桁展開」だけで、
★`F_k`-線形性は要らない。これが直前の波が見つけた空虚を避ける点である。 -/
theorem norm_iterate_sub_le_of_layers {p k : ℕ} {π : M} {t : ℕ → ℤ} {n : ℕ → ℕ}
    (τ : M →+ M) (s : ℕ → (M →+* M))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hs : ∀ j, j < k → ∀ z : M, s j z = (⇑τ)^[p ^ j] z)
    (hbrk : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (t (j + 1) + 1))
    (hexp : ∀ j, j < k → ∀ z : M, ∃ a : ℕ → M, (∀ m, m < n j → s j (a m) = a m) ∧
        z = ∑ m ∈ Finset.range (n j), a m * π ^ m ∧
        ∀ m, m < n j → ‖a m * π ^ m‖ ≤ ‖z‖) :
    ∀ j, j < k → ∀ z : M, ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖ := by
  intro j hj z
  obtain ⟨a, hfix, hz, hle⟩ := hexp j hj z
  rw [← hs j hj z]
  exact norm_ringHom_sub_le_zpow_mul (s j) hπ0 hπ1 (ht0 j) (hbrk j hj) hfix hz hle

end Layers


/-! ## §4 出口 -/

section Exit

variable {F : Type*} [Field F] [Algebra F M]

/-- ★★★★`GainedJumpSeq...lt` の `hstep` を層ごとの桁展開に落とした形。

`hstep`(「すべての `z ∈ M`」の縮小率)が
`hbrk`(★素元 1 点の条件)と `hexp`(層ごとの桁展開)に分解される。
★`τ : M →+ M` のままで、`1 ≤ k` を許す。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_layers
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {n : ℕ → ℕ}
    (τ : M →+ M) (σ : M →ₐ[F] M) (s : ℕ → (M →+* M)) (t : ℕ → ℤ) (he : 0 < e)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hs : ∀ j, j < k → ∀ z : M, s j z = (⇑τ)^[p ^ j] z)
    (hbrk : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (t (j + 1) + 1))
    (hexp : ∀ j, j < k → ∀ z : M, ∃ a : ℕ → M, (∀ m, m < n j → s j (a m) = a m) ∧
        z = ∑ m ∈ Finset.range (n j), a m * π ^ m ∧
        ∀ m, m < n j → ‖a m * π ^ m‖ ≤ ‖z‖)
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
  obtain ⟨hπ0, hπ1⟩ := norm_pos_and_lt_one (π := π) hp'.one_lt
    (Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne') hnormp heM
  exact GainedJumpSeq.exists_norm_sub_algebraMap_le_prod_axDecay_lt τ σ t he hσ ht0
    (norm_iterate_sub_le_of_layers τ s hπ0 hπ1 ht0 hs hbrk hexp)
    hstepZ hlayerZ hi hti hnormp heM hbreak hval hdeg htop x

/-- ★★★★★塔の形。層 `j` の体 `E j` を与えると `hexp` が自動に出る。

`[M : E j] = p^{k+1−j}`、`M = (E j)(π)`、全分岐(`hvalj`)、
`s j`(`= τ^{p^j}`)が `E j` を点ごとに固定すること(`hfixj`)。

★★`E : ℕ → Type*` は `IntermediateField` ではない。
中間体を作らないのは #59(2 層をまたぐ `rfl` が kernel を止める)の回避であり、
実測で本ファイルは #59/#69 に一度も当たっていない。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_tower
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M}
    {E : ℕ → Type*} [∀ j, Field (E j)] [∀ j, Algebra (E j) M]
    (τ : M →+ M) (σ : M →ₐ[F] M) (s : ℕ → (M →+* M)) (t : ℕ → ℤ) (he : 0 < e)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hs : ∀ j, j < k → ∀ z : M, s j z = (⇑τ)^[p ^ j] z)
    (hdegj : ∀ j, j < k → (minpoly (E j) π).natDegree = p ^ (k + 1 - j))
    (htopj : ∀ j, j < k → Algebra.adjoin (E j) ({π} : Set M) = ⊤)
    (hvalj : ∀ j, j < k → ∀ d : E j, d ≠ 0 →
        ∃ m : ℤ, ‖algebraMap (E j) M d‖ = ‖π‖ ^ (((p ^ (k + 1 - j) : ℕ) : ℤ) * m))
    (hfixj : ∀ j, j < k → ∀ d : E j, s j (algebraMap (E j) M d) = algebraMap (E j) M d)
    (hbrk : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (t (j + 1) + 1))
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
  obtain ⟨hπ0, hπ1⟩ := norm_pos_and_lt_one (π := π) hp'.one_lt
    (Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne') hnormp heM
  refine exists_norm_sub_algebraMap_le_prod_axDecay_of_layers (n := fun j => p ^ (k + 1 - j))
    τ σ s t he hσ ht0 hs hbrk ?_ hstepZ hlayerZ hi hti hnormp heM hbreak hval hdeg htop x
  intro j hj z
  exact exists_expansion_of_adjoin hπ0 hπ1 (pow_pos hp'.pos _) (hvalj j hj) (hdegj j hj)
    (htopj j hj) (s j) (hfixj j hj) z


/-- ★★★★★★本ファイルの最終出口。塔の形から `ℤ` 側の仮説 3 本を落とした。

`ht0` / `hstepZ` / `hlayerZ` は `GainedJumpSeq.exists_jump_seq` が
`hu`(`p^m ≤ u m`)と `hbnd`(頂点 1 本)から作る列 `t` に吸収される。
★入力は実際の跳びの列 `u` と塔のデータだけで、`τ : M →+ M`、`1 ≤ k` 可。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_jumps
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    {E : ℕ → Type*} [∀ j, Field (E j)] [∀ j, Algebra (E j) M]
    (τ : M →+ M) (σ : M →ₐ[F] M) (s : ℕ → (M →+* M)) (he : 0 < e)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z)
    (hs : ∀ j, j < k → ∀ z : M, s j z = (⇑τ)^[p ^ j] z)
    (hu : ∀ m, m ≤ k → (p : ℤ) ^ m ≤ u m)
    (hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hjump : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hdegj : ∀ j, j < k → (minpoly (E j) π).natDegree = p ^ (k + 1 - j))
    (htopj : ∀ j, j < k → Algebra.adjoin (E j) ({π} : Set M) = ⊤)
    (hvalj : ∀ j, j < k → ∀ d : E j, d ≠ 0 →
        ∃ m : ℤ, ‖algebraMap (E j) M d‖ = ‖π‖ ^ (((p ^ (k + 1 - j) : ℕ) : ℤ) * m))
    (hfixj : ∀ j, j < k → ∀ d : E j, s j (algebraMap (E j) M d) = algebraMap (E j) M d)
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hpZ : (2 : ℤ) ≤ (p : ℤ) := by exact_mod_cast hp'.two_le
  obtain ⟨hπ0, hπ1⟩ := norm_pos_and_lt_one (π := π) hp'.one_lt
    (Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne') hnormp heM
  obtain ⟨t, htk, ht1, htu, hstepZ, hlayerZ⟩ :=
    GainedJumpSeq.exists_jump_seq hpZ (e := (e : ℤ)) k hu hbnd
  have ht0 : ∀ j, 0 ≤ t (j + 1) := fun j => le_trans zero_le_one (ht1 j)
  have hti' : (i : ℤ) = t (k + 1) := by rw [htk]; exact hti
  have hi : 0 < i := by
    have h1 : (1 : ℤ) ≤ u k := le_trans (one_le_pow₀ (by linarith)) (hu k (le_refl k))
    omega
  have hbrk : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (t (j + 1) + 1) := by
    intro j hj
    refine le_trans (hjump j hj) ?_
    refine zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) ?_
    have := htu j (le_of_lt hj)
    omega
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_tower τ σ s t he hσ ht0 hs
    hdegj htopj hvalj hfixj hbrk hstepZ hlayerZ hi hti' hnormp heM hbreak hval hdeg htop x

end Exit


/-! ## §5 非空虚性 -/

section Vacuity

/-- ★★★空虚性の分かれ目を純群論で書いたもの。

直前の波が見つけた 6 本は `τ` が `F`-線形で `[M:F] = p` だったため
`orderOf τ ∣ p` となり、`k ≥ 1` で `τ^{p^k} = 1` が強制されて空虚になった。
★塔では `τ` に効く体は底 `E 0` で `[M : E 0] = p^{k+1}` なので
`orderOf τ = p^{k+1}` が取れ、本補題の通り `τ^{p^k} ≠ 1` と両立する。 -/
theorem pow_ne_one_of_orderOf_eq {G : Type*} [Group G] {p k : ℕ} (hp : 1 < p) (τ : G)
    (hord : orderOf τ = p ^ (k + 1)) : τ ^ p ^ k ≠ 1 := by
  intro h
  have hdvd : orderOf τ ∣ p ^ k := orderOf_dvd_of_pow_eq_one h
  rw [hord] at hdvd
  have hlt : p ^ k < p ^ (k + 1) := Nat.pow_lt_pow_right hp (Nat.lt_succ_self k)
  have hpos : 0 < p ^ k := pow_pos (lt_trans Nat.zero_lt_one hp) k
  exact absurd (Nat.le_of_dvd hpos hdvd) (not_le.mpr hlt)

/-- ★位数がちょうど `n` の群元は実在する(`Multiplicative (ZMod n)` の生成元)。
★これで `orderOf τ = p^{k+1}` という要求が空虚でないことが機械検算される。 -/
theorem orderOf_ofAdd_one_zmod (n : ℕ) :
    orderOf (Multiplicative.ofAdd (1 : ZMod n)) = n := by
  rw [orderOf_ofAdd_eq_addOrderOf, ZMod.addOrderOf_one]

/-- ★★`k ≥ 1` でも矛盾しないことの群論的な証拠。
`Multiplicative (ZMod (p^{k+1}))` の生成元 `τ` は `τ^{p^{k+1}} = 1` かつ `τ^{p^k} ≠ 1`。 -/
theorem exists_group_pow_ne_one (p k : ℕ) (hp : 1 < p) :
    (Multiplicative.ofAdd (1 : ZMod (p ^ (k + 1)))) ^ p ^ k ≠ 1 :=
  pow_ne_one_of_orderOf_eq hp _ (orderOf_ofAdd_one_zmod _)

end Vacuity

section BottomOrder

variable {E : Type*} [Field E] [Algebra E M]

/-- 環準同型で `E` を点ごとに固定するものは `E`-代数準同型。 -/
def toAlgHom (σ : M →+* M) (h : ∀ d : E, σ (algebraMap E M d) = algebraMap E M d) :
    M →ₐ[E] M :=
  { σ with commutes' := h }

omit [IsUltrametricDist M] in
/-- `[M:E] = n` なら `E` を固定する環準同型の `n` 回反復は恒等。
★`PureStepSetup.algHom_pow_apply_eq_self_of_natDegree` を `M →+* M` 向けに言い換えただけ。 -/
theorem iterate_natDegree_eq_self {n : ℕ} {π : M} (hn : 0 < n)
    (hdeg : (minpoly E π).natDegree = n) (htop : Algebra.adjoin E ({π} : Set M) = ⊤)
    (σ : M →+* M) (h : ∀ d : E, σ (algebraMap E M d) = algebraMap E M d) (z : M) :
    (⇑σ)^[n] z = z := by
  have hh := PureStepSetup.algHom_pow_apply_eq_self_of_natDegree hn hdeg htop
    (toAlgHom (E := E) σ h) z
  rw [AlgHom.coe_pow] at hh
  exact hh

omit [IsUltrametricDist M] in
/-- ★★★直前の波の「`τ^p = 1`」の正しい置き換え。

塔の底 `E 0` で `[M : E 0] = p^{k+1}` なので、出るのは `τ^{p^{k+1}} = 1` だけであり、
★`τ^{p^k} ≠ 1`(⇐ `hbreak`)と矛盾しない。これが `k ≥ 1` を生かす点である。 -/
theorem iterate_pow_succ_eq_self {p k : ℕ} {π : M} (hp : 0 < p) (τ : M →+ M) (s0 : M →+* M)
    (hs0 : ∀ z : M, s0 z = (⇑τ)^[p ^ 0] z)
    (hdeg0 : (minpoly E π).natDegree = p ^ (k + 1))
    (htop0 : Algebra.adjoin E ({π} : Set M) = ⊤)
    (hfix0 : ∀ d : E, s0 (algebraMap E M d) = algebraMap E M d) (z : M) :
    (⇑τ)^[p ^ (k + 1)] z = z := by
  have hfun : (⇑s0 : M → M) = ⇑τ := by
    funext w
    rw [hs0 w]
    simp
  have hh := iterate_natDegree_eq_self (E := E) (pow_pos hp (k + 1)) hdeg0 htop0 s0 hfix0 z
  rwa [hfun] at hh

end BottomOrder

namespace Numeric

/-- `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` の実データから作った列。
`p = 3`, `k = 1`, `e = 2`, `u = (u₀,u₁) = (2,8)`, `i = 8`。
`t₁ = 2`, `t₂ = 8`, `t_{j+2} = 8·3^j`。 -/
def tz : ℕ → ℤ
  | 0 => 0
  | 1 => 2
  | (n + 2) => 8 * 3 ^ n

theorem tz_nonneg : ∀ j, 0 ≤ tz (j + 1)
  | 0 => by norm_num [tz]
  | (n + 1) => by simp only [tz]; positivity

theorem tz_step : ∀ j, (3 : ℤ) * tz (j + 1) ≤ tz (j + 2)
  | 0 => by norm_num [tz]
  | 1 => by norm_num [tz]
  | (n + 2) => by simp only [tz]; ring_nf; exact le_rfl

theorem tz_layer : ∀ j, ((3 : ℤ) - 1) * tz (j + 1) ≤ (3 : ℤ) ^ (j + 1) * (2 : ℤ)
  | 0 => by norm_num [tz]
  | (n + 1) => by
      have h : (0 : ℤ) < 3 ^ n := pow_pos (by norm_num) n
      have h9 : (3 : ℤ) ^ (n + 1 + 1) = 9 * 3 ^ n := by ring
      simp only [tz, h9]
      linarith

theorem tz_top : (8 : ℤ) = tz (1 + 1) := by norm_num [tz]

end Numeric


/-! ## §6 `.src`(原典の対応箇所) -/

def norm_ringHom_sub_le_zpow_mul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_expansion_of_adjoin.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_iterate_sub_le_of_layers.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_tower.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_jumps.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def pow_ne_one_of_orderOf_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §7 使っている公理の一覧 -/

#print axioms norm_sum_le_of_forall_le
#print axioms norm_pow_sub_pow_le
#print axioms norm_ringHom_sub_le_zpow_mul
#print axioms exists_expansion_of_adjoin
#print axioms norm_pos_and_lt_one
#print axioms norm_iterate_sub_le_of_layers
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_layers
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_tower
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_jumps
#print axioms pow_ne_one_of_orderOf_eq
#print axioms orderOf_ofAdd_one_zmod
#print axioms exists_group_pow_ne_one
#print axioms iterate_natDegree_eq_self
#print axioms iterate_pow_succ_eq_self
#print axioms Numeric.tz_nonneg
#print axioms Numeric.tz_step
#print axioms Numeric.tz_layer
#print axioms Numeric.tz_top

end GainedTowerStep
end ABC3.Found.PGC
