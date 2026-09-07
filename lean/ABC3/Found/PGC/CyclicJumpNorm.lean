import ABC3.Meta.Claim
import Mathlib.Analysis.Normed.Field.Ultra
import Mathlib.Analysis.Normed.Group.Ultra
import Mathlib.Algebra.Ring.GeomSum

/-!
# [pGC] 巡回 `p` 次拡大の跳び `i` と `‖σx − x‖` の関係

`Found/PGC/SenLemma.lean` は Ax–Sen–Tate を **`AxDescentStep K C`(1 段の降下)ただ 1 点**に
落とし、tame(`p ∤ [K(x):K]`)な側は `descentStep_of_natDegree_tame` で埋めた。
残っていたのは wild(`p ∣ [K(x):K]`)な 1 段で、その古典的な定数は
**`|π_L|^{−i}`(`i` は下付き分岐群の跳び)**である。本ファイルはその **`|π_L|^{−i}` を
厳密な等式として出す**(不等式ではない)。

## 在庫調査(自分で測った。★MCP は 0 回。コマンドを残す)

```
grep -n "norm_add_eq_max\|nnnorm_add_eq_max" .cache/mathlib-index.txt
  → PadicInt.norm_add_eq_max_of_ne など。一般の超距離版は名前が違う
grep -n "Analysis/Normed/Group/Ultra.lean" .cache/mathlib-index.txt
  → ★★IsUltrametricDist.nnnorm_prod_eq_sup_of_pairwise_ne (Ultra.lean:358)
      「ノルムが相異なる有限積のノルム = sup」。★**乗法版しか索引に出ない**
  → IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm / norm_add_le_max /
     norm_sum_le_of_forall_le_of_nonneg
★★ここが本波の最大の当たり: 索引は `to_additive` の生成名を持たないので
   `grep` では **加法版が 0 件に見える**。`leanfile.mjs` に `#check` を投げると
   ★**`IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne` は在る**。
   ⇒「桁の分離(相異なるノルム)ならば和のノルムは最大値」は**自前で書く必要が無かった**。
   ★教訓:「索引に無い」は「mathlib に無い」ではない(`to_additive` の生成名)。
grep -n "mul_finset_sup\|Finset.sup_mul" .cache/mathlib-index.txt
  → ★NNReal.mul_finset_sup (Data/NNReal/Basic.lean:115)。sup とスカラー倍の交換
grep -n "geom_sum₂_mul" .cache/mathlib-index.txt
  → ★geom_sum₂_mul (Algebra/Ring/GeomSum.lean:306) `(Σ x^i y^{n-1-i})(x−y) = x^n − y^n`
grep -n "zpow_right_injective\|eq_zero_of_abs_lt_dvd" .cache/mathlib-index.txt
  → ★zpow_right_injective₀ / Int.eq_zero_of_abs_lt_dvd
grep -c "PGC\.<名前>\b" .cache/decl-index.txt (本ファイルの全宣言名について)
  → すべて 0(衝突なし)。★ただし `ABC3.Found.PGC.norm_pow_sub_pow_le`
     (Found/PGC/PadicLogSurjective.lean:39)が**既に在る**(`K.carrier` 版の不等式)ので、
     本ファイルの補題は `norm_pow_succ_sub_pow_succ_le` に改名した。
```

★木の在庫で**使わなかったもの**とその理由:

* `Found/PGC/LowerRamificationGroup.lean` / `UniformizerExpansion.lean` は
  跳びを **`addVal`(`ℕ∞`)言語**で持っている(`ramIndex α σ = addVal (σ•α − α)`)。
  本ファイルは**ノルム言語**で書いた。`AxDescentStep` がノルム言語だからである。
  ★橋は `LubinTateRamificationBreak.lean` の
  `norm_eq_pow_of_addVal_eq` / `addVal_eq_of_norm_eq_pow` が既に持っている。
* ★`UniformizerExpansion.exists_digits` の桁は `∏_{i<n} σ^i π`(Yoshida の基底)であって
  `π^j` ではない。本ファイルは `π^j` を使う。**貼り合わせるときはここが食い違う。**

## 何が言えたか

### §1 抽象核(★分岐・付値・Galois・`p` 進の語彙が **1 語も出てこない**)

| 宣言 | 内容 |
|---|---|
| `norm_eq_of_norm_sub_lt` | `‖u−π‖ < ‖π‖ ⇒ ‖u‖ = ‖π‖` |
| `norm_pow_succ_sub_pow_succ_le` | `‖u^{i+1} − π^{i+1}‖ ≤ ‖u−π‖·‖π‖^i` |
| ★`norm_pow_sub_pow_eq` | **`‖(k:L)‖ = 1` なら `‖u^k − π^k‖ = ‖u−π‖·‖π‖^{k−1}`(等号)** |
| `zpow_mul_pow_ne` | 実数だけ。`c^{nm+j} = c^{nm′+l}`・`j,l<n` ⇒ `j = l` |
| ★`nnnorm_sum_mul_eq_mul_nnnorm_sum` | ノルムが**相似**な 2 族の和のノルムも相似 |

★`norm_pow_sub_pow_eq` が心臓である。超距離では一般に `‖u^k − π^k‖ ≤ ‖u−π‖‖π‖^{k−1}`
しか言えないが、**`k` が `L` で単数なら等号**になる。証明は
`u^k − π^k = (Σ_{i<k} u^i π^{k−1−i})(u−π)` の第 1 因子が
`k·π^{k−1}` と `‖π‖^{k−1}` より真に小さいところしか違わない、という 1 点。

### §2 桁展開の層

`digitSum n π a := Σ_{j<n} a_j π^j`(`a_j ∈ K`)。仮定は 2 つだけ:

* `hval`: `c ∈ K^×` なら `‖c‖ = ‖π‖^{n·m}`(★**「`L/K` は全分岐で `e = n`」の言い換え**)、
* `hchar`: `0 < j < n` なら `‖(j:L)‖ = 1`(★**「`n = p`、剰余標数 `p`」の言い換え**)。

| 宣言 | 主張 |
|---|---|
| `nnnorm_sum_digit_eq_sup` | **桁展開のノルム = 各桁のノルムの最大値**(桁の分離) |
| ★`norm_digitSum_sub_digit_zero_le` | **定数項 `a_0` が最良近似**(`∀c∈K`, `‖x−a_0‖ ≤ ‖x−c‖`) |
| ★★`norm_digitSum_sub_mul_norm_eq` | **`‖x^u − x‖·‖π‖ = ‖u−π‖·‖x − a_0‖`** |

★最後の 1 本が本ファイルの主結果である。`x^u` は `π` を `u` に置き換えた「共役」。
★★**等号**であることに注意(`≤` ではない)。

### §3 `σ` と跳び `i`

`σ : L →ₐ[K] L`(★**Galois も同型も等長性も要らない。`K`-代数準同型だけ**)を入れると
`σ x = x^{σπ}` なので §2 がそのまま使えて

| 宣言 | 主張 |
|---|---|
| `norm_algHom_sub_mul_norm_eq` | `‖σx − x‖·‖π‖ = ‖σπ − π‖·‖x − a_0‖` |
| ★`eq_algebraMap_digit_zero_of_algHom_eq` | **退化検査**: `σx = x` かつ `σπ ≠ π` なら `x ∈ K` |
| ★★`norm_sub_digit_zero_eq_zpow_mul` | 跳び `i`(`‖σπ−π‖ = ‖π‖^{i+1}`)のとき **`‖x − a_0‖ = ‖π‖^{−i}·‖σx − x‖`** |
| ★★`exists_norm_sub_algebraMap_eq_zpow_mul` | **`d(x,K) = ‖π‖^{−i}·‖σx − x‖`**(最良近似込み) |

★★これが持ち場の目標「**1 段の定数は `|π_L|^{−i}`**」そのものである。
★しかも**上からも下からも**押さえた等式なので、次の波が
「もっと良い定数が取れないか」を掘る必要は無い(`|π_L|^{−i}` が**最良**である)。

## ★★★測って分かったこと —— 「1 段が一様なら Ax が出る」は **偽**

☆★**本波の当初の見立ては外れた。**「1 段の定数が `x` に依らず一様(`p^{1/(p−1)}` 程度)なら
`AxDescentStep` が埋まる」と見ていたが、★**塔では定数が掛かり算で積み上がる**。

`x_0 = x`、`x_{j+1}` を 1 段降ろすと
`‖σ x_{j+1} − x_{j+1}‖ ≤ max(‖σx_j − x_j‖, ‖x_j − x_{j+1}‖) = C_j·ε_j`
なので `ε_{j+1} = C_j ε_j`。★超距離は**和**を `max` に潰すが、**積**は潰さない。
よって `a = v_p([K(x):K])` 段で `Π_j C_j` になり、`C_j ≥ 1` が一定なら**発散する**。

★Ax の定数の指数 `p/(p−1)^2` は `Σ_{k≥1} (1/(p−1))·p^{−(k−1)}` に等しい。
⇒ ★★**各段の損失 `i_j/e_j` が幾何級数的に減る**ことを使っている。
これは**上付き番号(Herbrand `φ`/`ψ`)**の話であって、1 段の話ではない。
★★**次の波はここ(塔に沿った `i_j/e_j` の減衰)を掘ること。**
★木には `HerbrandFunction.lean` / `HerbrandComposition.lean` /
`UpperRamificationGroup.lean` / `HasseArf*.lean` が既に在る。

## ★もう 1 つ測ったこと —— `i ≤ e_L` では既存より良くならない

本結果と木の `AxLemma.exists_norm_sub_algebraMap_le_div_norm_natDegree`
(`d(x,K) ≤ ε/‖n‖`、`n = p` なら定数 `|p|^{−1}`)を突き合わせると `‖p‖ ≤ ‖π‖^i`
すなわち **`i ≤ e_L`** が出る。★しかしそれは `‖π‖^{−i} ≤ |p|^{−1}` を意味するだけで、
**既存の定数と同じ**である。★価値があるのは sharp な **`(p−1)·i ≤ e_L`** の方で、
それには different の評価(`d = (p−1)(i+1) ≤ e_L + p − 1`)が要る。★**本波は未着手**。
★この段落は「次の波が `i ≤ e_L` を証明して満足するのを止めるため」に書いてある。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★原典(Serre, *Corps Locaux* IV / Ax 1970)は **付値言語**で `i_G(σ) = v_L(σπ − π)` を
   使う。本ファイルは **ノルム言語**で書いた(`AxDescentStep` がノルム言語だから)。
   跳び `i` は `‖σπ − π‖ = ‖π‖^{i+1}` として現れる。★同値である。
2. ★「`L/K` は全分岐な巡回 `p` 次拡大」という**分岐理論の仮定を一切置いていない**。
   代わりに `hval`(値群が `n` 乗部分群)・`hchar`(`0<j<n` で `‖j‖=1`)・
   `x = Σ_{j<n} a_j π^j` という**ノルムと代数だけの仮定**に分解した。
   ★★その結果 `σ` に **Galois性・全単射性・等長性のどれも要らない**ことが分かった
   (`K`-代数準同型で十分)。★これは原典より弱い仮定である。
3. `.src` は `SenLemma.lean` と**同じ項目**(pGC 物理 p.6 Corollary 3.1)を指す。
   本ファイルの中身は `AxSenTate` の入力であって、原典が独立に立てた項目ではない。
4. `n` は素数だと仮定していない(`hchar` が実質それを担う)。`p` という文字も使わない。
-/

namespace ABC3.Found.PGC

open Finset

/-! ## §1 抽象核

★以下の 5 本には**分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない**。 -/

section AbstractCore

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- **抽象核** —— `‖u − π‖ < ‖π‖` なら `‖u‖ = ‖π‖`。超距離の等号条件そのもの。 -/
theorem norm_eq_of_norm_sub_lt {u π : L} (h : ‖u - π‖ < ‖π‖) : ‖u‖ = ‖π‖ := by
  have hu : u = (u - π) + π := by ring
  rw [hu, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (ne_of_lt h), max_eq_right h.le]

/-- **抽象核** —— `‖u^{i+1} − π^{i+1}‖ ≤ ‖u − π‖·‖π‖^i`。★不等号の側は `k` の仮定が要らない。

★木には `ABC3.Found.PGC.norm_pow_sub_pow_le`(`Found/PGC/PadicLogSurjective.lean`)が
既に在るが、あちらは `K.carrier` 専用なので、ここでは一般の超距離ノルム体で書き直した。 -/
theorem norm_pow_succ_sub_pow_succ_le {u π : L} (h : ‖u - π‖ < ‖π‖) (i : ℕ) :
    ‖u ^ (i + 1) - π ^ (i + 1)‖ ≤ ‖u - π‖ * ‖π‖ ^ i := by
  have hu : ‖u‖ = ‖π‖ := norm_eq_of_norm_sub_lt h
  induction i with
  | zero => simp
  | succ j ih =>
    have hid : u ^ (j + 1 + 1) - π ^ (j + 1 + 1)
        = u * (u ^ (j + 1) - π ^ (j + 1)) + (u - π) * π ^ (j + 1) := by ring
    rw [hid]
    refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ ?_)
    · rw [norm_mul, hu]
      calc ‖π‖ * ‖u ^ (j + 1) - π ^ (j + 1)‖ ≤ ‖π‖ * (‖u - π‖ * ‖π‖ ^ j) :=
            mul_le_mul_of_nonneg_left ih (norm_nonneg _)
        _ = ‖u - π‖ * ‖π‖ ^ (j + 1) := by ring
    · rw [norm_mul, norm_pow]

/-- ★★★**抽象核(本ファイルの心臓)** —— `‖(k : L)‖ = 1` ならば

`‖u^k − π^k‖ = ‖u − π‖ · ‖π‖^{k−1}` (★**等号**)。

超距離では一般に `≤` しか言えない。等号になるのは
`u^k − π^k = (Σ_{i<k} u^i π^{k−1−i})·(u − π)` の第 1 因子が
`(k:L)·π^{k−1}` と `‖π‖^{k−1}` より**真に小さい**ところしか違わないからで、
その「真に小さい」に `‖u − π‖ < ‖π‖` を使う。
★`‖(k:L)‖ = 1` は具体層では「`0 < k < p` は `p` 進単数」の言い換えである。
★`k = 0` は `‖(0:L)‖ = 0 ≠ 1` により仮定から排除される。 -/
theorem norm_pow_sub_pow_eq {u π : L} (h : ‖u - π‖ < ‖π‖) {k : ℕ} (hk : ‖(k : L)‖ = 1) :
    ‖u ^ k - π ^ k‖ = ‖u - π‖ * ‖π‖ ^ (k - 1) := by
  have hπ0 : (0 : ℝ) < ‖π‖ := lt_of_le_of_lt (norm_nonneg _) h
  have hk0 : k ≠ 0 := by rintro rfl; simp at hk
  rcases eq_or_ne k 1 with rfl | hk1
  · simp
  have hk2 : 2 ≤ k := by omega
  set S : L := ∑ i ∈ range k, u ^ i * π ^ (k - 1 - i) with hSdef
  have hgeom : S * (u - π) = u ^ k - π ^ k := geom_sum₂_mul u π k
  have hconst : ∑ i ∈ range k, π ^ i * π ^ (k - 1 - i) = (k : L) * π ^ (k - 1) := by
    rw [Finset.sum_congr rfl (g := fun _ => π ^ (k - 1)) (fun i hi => by
      rw [← pow_add]
      congr 1
      have := Finset.mem_range.mp hi
      omega)]
    simp
  have hsplit : S - (k : L) * π ^ (k - 1)
      = ∑ i ∈ range k, (u ^ i - π ^ i) * π ^ (k - 1 - i) := by
    rw [hSdef, ← hconst, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  have hC : (0 : ℝ) ≤ ‖u - π‖ * ‖π‖ ^ (k - 2) :=
    mul_nonneg (norm_nonneg _) (pow_nonneg hπ0.le _)
  have hdiff : ‖S - (k : L) * π ^ (k - 1)‖ ≤ ‖u - π‖ * ‖π‖ ^ (k - 2) := by
    rw [hsplit]
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hC ?_
    intro i hi
    have hik : i < k := Finset.mem_range.mp hi
    rcases Nat.eq_zero_or_pos i with rfl | hipos
    · simp [hC]
    obtain ⟨j, rfl⟩ : ∃ j, i = j + 1 := ⟨i - 1, by omega⟩
    rw [norm_mul, norm_pow]
    calc ‖u ^ (j + 1) - π ^ (j + 1)‖ * ‖π‖ ^ (k - 1 - (j + 1))
        ≤ (‖u - π‖ * ‖π‖ ^ j) * ‖π‖ ^ (k - 1 - (j + 1)) :=
          mul_le_mul_of_nonneg_right (norm_pow_succ_sub_pow_succ_le h j) (pow_nonneg hπ0.le _)
      _ = ‖u - π‖ * ‖π‖ ^ (k - 2) := by
          rw [mul_assoc, ← pow_add]; congr 2; omega
  have hlt : ‖u - π‖ * ‖π‖ ^ (k - 2) < ‖π‖ ^ (k - 1) := by
    have hpow : ‖π‖ ^ (k - 1) = ‖π‖ * ‖π‖ ^ (k - 2) := by
      rw [← pow_succ']; congr 1; omega
    rw [hpow]
    exact mul_lt_mul_of_pos_right h (pow_pos hπ0 _)
  have hknorm : ‖(k : L) * π ^ (k - 1)‖ = ‖π‖ ^ (k - 1) := by
    rw [norm_mul, hk, one_mul, norm_pow]
  have hSnorm : ‖S‖ = ‖π‖ ^ (k - 1) := by
    have hid : S = (k : L) * π ^ (k - 1) + (S - (k : L) * π ^ (k - 1)) := by ring
    have hne : ‖(k : L) * π ^ (k - 1)‖ ≠ ‖S - (k : L) * π ^ (k - 1)‖ := by
      rw [hknorm]; exact ne_of_gt (lt_of_le_of_lt hdiff hlt)
    rw [hid, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm hne, hknorm,
      max_eq_left (le_of_lt (lt_of_le_of_lt hdiff hlt))]
  rw [← hgeom, norm_mul, hSnorm]
  ring

end AbstractCore

/-- **抽象核(実数だけ)** —— `0 < c ≠ 1` のとき、`j, l < n` で
`c^{nm}·c^j = c^{nm′}·c^l` なら `j = l`。

★これが「桁の分離」の全部である。`c = ‖π‖`、`c^{nm}` が `K^×` のノルムのとき、
指数が `n` を法として `j` を決めるので、桁どうしのノルムは相異なる。 -/
theorem zpow_mul_pow_ne {c : ℝ} (hc0 : 0 < c) (hc1 : c ≠ 1) {n j l : ℕ} {m m' : ℤ}
    (hj : j < n) (hl : l < n) (hjl : j ≠ l) :
    c ^ ((n : ℤ) * m) * c ^ j ≠ c ^ ((n : ℤ) * m') * c ^ l := by
  intro heq
  rw [← zpow_natCast c j, ← zpow_natCast c l, ← zpow_add₀ (ne_of_gt hc0),
    ← zpow_add₀ (ne_of_gt hc0)] at heq
  have hexp : (n : ℤ) * m + j = (n : ℤ) * m' + l := zpow_right_injective₀ hc0 hc1 heq
  have hdvd : (n : ℤ) ∣ ((l : ℤ) - j) := ⟨m - m', by linarith⟩
  have habs : |(l : ℤ) - j| < (n : ℤ) := abs_lt.mpr ⟨by omega, by omega⟩
  have hzero : (l : ℤ) - j = 0 := Int.eq_zero_of_abs_lt_dvd hdvd habs
  exact hjl (by omega)

/-- ★★**抽象核** —— 有限和のノルムの「相似」の保存。

各項のノルムが `‖v i‖·r = t·‖w i‖` という比例関係にあり、どちらの族も
**ノルムが互いに相異なる**なら、和についても `‖Σv‖·r = t·‖Σw‖`。

★超距離だから両辺とも `sup` に潰れ、`sup` はスカラー倍と交換する
(`NNReal.mul_finset_sup`)。★`ℝ≥0` で書くのは `Finset.sup` に `OrderBot` が要るから。
★★mathlib の `IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne` を使う。 -/
theorem nnnorm_sum_mul_eq_mul_nnnorm_sum {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {ι : Type*} (s : Finset ι) (v w : ι → M) (r t : NNReal)
    (hv : (s : Set ι).Pairwise fun i j => ‖v i‖₊ ≠ ‖v j‖₊)
    (hw : (s : Set ι).Pairwise fun i j => ‖w i‖₊ ≠ ‖w j‖₊)
    (h : ∀ i ∈ s, ‖v i‖₊ * r = t * ‖w i‖₊) :
    ‖∑ i ∈ s, v i‖₊ * r = t * ‖∑ i ∈ s, w i‖₊ := by
  rw [IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne hv,
    IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne hw, mul_comm, NNReal.mul_finset_sup,
    NNReal.mul_finset_sup]
  exact Finset.sup_congr rfl fun i hi => by rw [mul_comm]; exact h i hi

/-! ## §2 桁展開の層 -/

/-- **桁展開** `Σ_{j<n} a_j π^j`。★係数は `K`、`π` は `L` の元。 -/
def digitSum {K : Type*} [CommRing K] {L : Type*} [CommRing L] [Algebra K L]
    (n : ℕ) (π : L) (a : ℕ → K) : L :=
  ∑ j ∈ Finset.range n, algebraMap K L (a j) * π ^ j

section Digits

variable {K L : Type*} [Field K] [NormedField L] [IsUltrametricDist L] [Algebra K L]

omit [IsUltrametricDist L] in
/-- 桁展開から定数項を引くと 1 桁目以降の和になる。 -/
theorem digitSum_sub_digit_zero {n : ℕ} (hn : 0 < n) (π : L) (a : ℕ → K) :
    digitSum n π a - algebraMap K L (a 0)
      = ∑ j ∈ Finset.Ico 1 n, algebraMap K L (a j) * π ^ j := by
  rw [digitSum, Finset.range_eq_Ico,
    Finset.sum_eq_sum_Ico_succ_bot hn (fun j => algebraMap K L (a j) * π ^ j)]
  simp

omit [IsUltrametricDist L] in
/-- `K`-代数準同型は桁展開を `π ↦ σπ` に写す。★`σ` に全単射性も等長性も要らない。 -/
theorem map_digitSum (n : ℕ) (π : L) (a : ℕ → K) (σ : L →ₐ[K] L) :
    σ (digitSum n π a) = digitSum n (σ π) a := by
  simp [digitSum, map_sum, AlgHom.commutes]

/-- ★★**桁の分離** —— 桁展開のノルムは各桁のノルムの**最大値**である。

仮定 `hval`(`K^×` のノルムが `‖π‖` の `n` 乗部分群に入る)は
「`L/K` が全分岐で分岐指数が `n`」の言い換えで、これにより
`‖a_j π^j‖ = ‖π‖^{nm_j + j}` の指数が `n` を法として `j` を決め、桁どうしのノルムが
相異なる。★あとは mathlib の
`IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne` に流し込むだけ。 -/
theorem nnnorm_sum_digit_eq_sup {n : ℕ} {π : L} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m))
    (b : ℕ → K) {s : Finset ℕ} (hs : ∀ j ∈ s, j < n) :
    ‖∑ j ∈ s, algebraMap K L (b j) * π ^ j‖₊
      = s.sup fun j => ‖algebraMap K L (b j) * π ^ j‖₊ := by
  classical
  set t : Finset ℕ := s.filter (fun j => b j ≠ 0) with htdef
  have hsub : t ⊆ s := Finset.filter_subset _ _
  have hzero : ∀ j ∈ s, j ∉ t → algebraMap K L (b j) * π ^ j = 0 := by
    intro j hj hjt
    have hbj : b j = 0 := by
      by_contra hc
      exact hjt (Finset.mem_filter.mpr ⟨hj, hc⟩)
    simp [hbj]
  have h1 : ∑ j ∈ t, algebraMap K L (b j) * π ^ j = ∑ j ∈ s, algebraMap K L (b j) * π ^ j :=
    Finset.sum_subset hsub hzero
  have h2 : (t.sup fun j => ‖algebraMap K L (b j) * π ^ j‖₊)
      = s.sup fun j => ‖algebraMap K L (b j) * π ^ j‖₊ := by
    refine le_antisymm (Finset.sup_mono hsub) (Finset.sup_le fun j hj => ?_)
    by_cases hjt : j ∈ t
    · apply Finset.le_sup hjt
    · rw [hzero j hj hjt]; simp
  rw [← h1, ← h2]
  refine IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne ?_
  intro j hj l hl hjl
  have hjm : j ∈ t := by simpa using hj
  have hlm : l ∈ t := by simpa using hl
  obtain ⟨hjs, hjb⟩ := Finset.mem_filter.mp hjm
  obtain ⟨hls, hlb⟩ := Finset.mem_filter.mp hlm
  obtain ⟨m, hm⟩ := hval (b j) hjb
  obtain ⟨m', hm'⟩ := hval (b l) hlb
  intro hcontra
  have hreal : ‖algebraMap K L (b j) * π ^ j‖ = ‖algebraMap K L (b l) * π ^ l‖ := by
    simpa using congrArg NNReal.toReal hcontra
  rw [norm_mul, norm_mul, norm_pow, norm_pow, hm, hm'] at hreal
  exact zpow_mul_pow_ne hπ0 (ne_of_lt hπ1) (hs j hjs) (hs l hls) hjl hreal

/-- ★★**定数項が最良近似** —— `‖x − a_0‖ ≤ ‖x − c‖`(`∀ c ∈ K`)。

★これにより `‖x − a_0‖` は `d(x, K)`(`K` への距離)**そのもの**である。
★証明は「`c` を引くのは 0 桁目を `a_0 − c` に取り替えるだけ」で、
桁の分離(`nnnorm_sum_digit_eq_sup`)から `sup` の単調性で出る。 -/
theorem norm_digitSum_sub_digit_zero_le {n : ℕ} (hn : 0 < n) {π : L} (hπ0 : 0 < ‖π‖)
    (hπ1 : ‖π‖ < 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m))
    (a : ℕ → K) (c : K) :
    ‖digitSum n π a - algebraMap K L (a 0)‖ ≤ ‖digitSum n π a - algebraMap K L c‖ := by
  classical
  set a' : ℕ → K := Function.update a 0 (a 0 - c) with ha'
  have hupd : ∀ j ∈ Finset.Ico 1 n,
      algebraMap K L (a' j) * π ^ j = algebraMap K L (a j) * π ^ j := by
    intro j hj
    have hj0 : j ≠ 0 := by have := (Finset.mem_Ico.mp hj).1; omega
    rw [ha', Function.update_of_ne hj0]
  have hrepr : digitSum n π a - algebraMap K L c
      = ∑ j ∈ Finset.range n, algebraMap K L (a' j) * π ^ j := by
    rw [Finset.range_eq_Ico,
      Finset.sum_eq_sum_Ico_succ_bot hn (fun j => algebraMap K L (a' j) * π ^ j),
      Finset.sum_congr rfl hupd, ha', Function.update_self, digitSum, Finset.range_eq_Ico,
      Finset.sum_eq_sum_Ico_succ_bot hn (fun j => algebraMap K L (a j) * π ^ j)]
    simp [map_sub]
    ring
  have hnn1 : ‖digitSum n π a - algebraMap K L (a 0)‖₊
      = (Finset.Ico 1 n).sup fun j => ‖algebraMap K L (a j) * π ^ j‖₊ := by
    rw [digitSum_sub_digit_zero hn]
    exact nnnorm_sum_digit_eq_sup hπ0 hπ1 hval a fun j hj => (Finset.mem_Ico.mp hj).2
  have hnn2 : ‖digitSum n π a - algebraMap K L c‖₊
      = (Finset.range n).sup fun j => ‖algebraMap K L (a' j) * π ^ j‖₊ := by
    rw [hrepr]
    exact nnnorm_sum_digit_eq_sup hπ0 hπ1 hval a' fun j hj => Finset.mem_range.mp hj
  have hle : ‖digitSum n π a - algebraMap K L (a 0)‖₊
      ≤ ‖digitSum n π a - algebraMap K L c‖₊ := by
    rw [hnn1, hnn2]
    refine Finset.sup_le fun j hj => ?_
    rw [← hupd j hj]
    apply Finset.le_sup (Finset.mem_range.mpr (Finset.mem_Ico.mp hj).2)
  simpa using NNReal.coe_le_coe.mpr hle

/-- ★★★**本ファイルの主結果(§2)** —— 桁展開の「共役」との差の厳密な関係:

`‖x^u − x‖ · ‖π‖ = ‖u − π‖ · ‖x − a_0‖`   (`x = Σ_{j<n} a_j π^j`, `x^u = Σ_{j<n} a_j u^j`)

★★**等号**である。左辺の `‖π‖` と右辺の `‖u − π‖` を移項すれば
`‖x − a_0‖ = (‖π‖/‖u−π‖)·‖x^u − x‖` で、`u = σπ`・`‖σπ−π‖ = ‖π‖^{i+1}` のとき
係数はちょうど **`‖π‖^{−i}`** になる(§3)。

証明は 3 行:
1. 各桁で `‖a_j(u^j − π^j)‖·‖π‖ = ‖u−π‖·‖a_j π^j‖`(★`norm_pow_sub_pow_eq` の等号)、
2. 桁どうしのノルムは相異なる(`zpow_mul_pow_ne`)、
3. よって和でも比例が保たれる(`nnnorm_sum_mul_eq_mul_nnnorm_sum`)。

★`u = π` のときは両辺 `0` で真(退化検査)。 -/
theorem norm_digitSum_sub_mul_norm_eq {n : ℕ} {π u : L} {a : ℕ → K}
    (hn : 0 < n) (hπ1 : ‖π‖ < 1) (hu : ‖u - π‖ < ‖π‖)
    (hchar : ∀ j, 0 < j → j < n → ‖(j : L)‖ = 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    ‖digitSum n u a - digitSum n π a‖ * ‖π‖
      = ‖u - π‖ * ‖digitSum n π a - algebraMap K L (a 0)‖ := by
  classical
  have hπ0 : (0 : ℝ) < ‖π‖ := lt_of_le_of_lt (norm_nonneg _) hu
  rcases eq_or_ne u π with rfl | hune
  · simp
  have hupos : (0 : ℝ) < ‖u - π‖ := norm_pos_iff.mpr (sub_ne_zero.mpr hune)
  set S : Finset ℕ := (Finset.Ico 1 n).filter (fun j => a j ≠ 0) with hSdef
  set v : ℕ → L := fun j => algebraMap K L (a j) * (u ^ j - π ^ j) with hvdef
  set w : ℕ → L := fun j => algebraMap K L (a j) * π ^ j with hwdef
  have hSsub : S ⊆ Finset.Ico 1 n := Finset.filter_subset _ _
  have hSmem : ∀ j ∈ S, 1 ≤ j ∧ j < n ∧ a j ≠ 0 := by
    intro j hj
    have h1 := Finset.mem_filter.mp hj
    have h2 := Finset.mem_Ico.mp h1.1
    exact ⟨h2.1, h2.2, h1.2⟩
  have hsum1 : digitSum n u a - digitSum n π a = ∑ j ∈ S, v j := by
    rw [digitSum, digitSum, ← Finset.sum_sub_distrib,
      Finset.sum_congr rfl (g := v) (fun j _ => by rw [hvdef]; ring)]
    refine (Finset.sum_subset (fun j hj => ?_) ?_).symm
    · exact Finset.mem_range.mpr (hSmem j hj).2.1
    · intro j _ hj
      rcases Nat.eq_zero_or_pos j with rfl | hjpos
      · simp
      by_cases haj : a j = 0
      · simp [haj]
      · exact absurd (Finset.mem_filter.mpr ⟨Finset.mem_Ico.mpr ⟨hjpos,
          Finset.mem_range.mp (by assumption)⟩, haj⟩) hj
  have hsum2 : digitSum n π a - algebraMap K L (a 0) = ∑ j ∈ S, w j := by
    rw [digitSum_sub_digit_zero hn]
    refine (Finset.sum_subset hSsub ?_).symm
    intro j hj hjS
    have hja : a j = 0 := by
      by_contra hc
      exact hjS (Finset.mem_filter.mpr ⟨hj, hc⟩)
    simp [hja]
  have hvw : ∀ j ∈ S, ‖v j‖ * ‖π‖ = ‖u - π‖ * ‖w j‖ := by
    intro j hj
    obtain ⟨hj1, hjn, -⟩ := hSmem j hj
    have hpow : ‖π‖ ^ (j - 1) * ‖π‖ = ‖π‖ ^ j := by rw [← pow_succ]; congr 1; omega
    show ‖algebraMap K L (a j) * (u ^ j - π ^ j)‖ * ‖π‖
        = ‖u - π‖ * ‖algebraMap K L (a j) * π ^ j‖
    rw [norm_mul, norm_mul, norm_pow, norm_pow_sub_pow_eq hu (hchar j (by omega) hjn)]
    calc ‖algebraMap K L (a j)‖ * (‖u - π‖ * ‖π‖ ^ (j - 1)) * ‖π‖
        = ‖u - π‖ * (‖algebraMap K L (a j)‖ * (‖π‖ ^ (j - 1) * ‖π‖)) := by ring
      _ = ‖u - π‖ * (‖algebraMap K L (a j)‖ * ‖π‖ ^ j) := by rw [hpow]
  have hwdist : ∀ j ∈ S, ∀ l ∈ S, j ≠ l → ‖w j‖ ≠ ‖w l‖ := by
    intro j hj l hl hjl
    obtain ⟨-, hjn, haj⟩ := hSmem j hj
    obtain ⟨-, hln, hal⟩ := hSmem l hl
    obtain ⟨m, hm⟩ := hval (a j) haj
    obtain ⟨m', hm'⟩ := hval (a l) hal
    show ‖algebraMap K L (a j) * π ^ j‖ ≠ ‖algebraMap K L (a l) * π ^ l‖
    rw [norm_mul, norm_mul, norm_pow, norm_pow, hm, hm']
    exact zpow_mul_pow_ne hπ0 (ne_of_lt hπ1) hjn hln hjl
  have hvdist : ∀ j ∈ S, ∀ l ∈ S, j ≠ l → ‖v j‖ ≠ ‖v l‖ := by
    intro j hj l hl hjl hEq
    refine hwdist j hj l hl hjl ?_
    have h1 := hvw j hj
    have h2 := hvw l hl
    rw [hEq] at h1
    exact mul_left_cancel₀ (ne_of_gt hupos) (by rw [← h1, h2])
  have hnn : ∀ x y : L, ‖x‖ ≠ ‖y‖ → ‖x‖₊ ≠ ‖y‖₊ := by
    intro x y h he
    exact h (by simpa using congrArg NNReal.toReal he)
  have hpv : (S : Set ℕ).Pairwise fun i j => ‖v i‖₊ ≠ ‖v j‖₊ := fun i hi j hj hij =>
    hnn _ _ (hvdist i (by simpa using hi) j (by simpa using hj) hij)
  have hpw : (S : Set ℕ).Pairwise fun i j => ‖w i‖₊ ≠ ‖w j‖₊ := fun i hi j hj hij =>
    hnn _ _ (hwdist i (by simpa using hi) j (by simpa using hj) hij)
  have hkey : ‖∑ j ∈ S, v j‖₊ * ‖π‖₊ = ‖u - π‖₊ * ‖∑ j ∈ S, w j‖₊ := by
    refine nnnorm_sum_mul_eq_mul_nnnorm_sum S v w ‖π‖₊ ‖u - π‖₊ hpv hpw ?_
    intro j hj
    have h := hvw j hj
    have hc : ((‖v j‖₊ * ‖π‖₊ : NNReal) : ℝ) = ((‖u - π‖₊ * ‖w j‖₊ : NNReal) : ℝ) := by
      push_cast
      exact h
    exact_mod_cast hc
  rw [hsum1, hsum2]
  have hR : ((‖∑ j ∈ S, v j‖₊ * ‖π‖₊ : NNReal) : ℝ)
      = ((‖u - π‖₊ * ‖∑ j ∈ S, w j‖₊ : NNReal) : ℝ) := by rw [hkey]
  push_cast at hR
  exact hR

end Digits

/-! ## §3 `σ` と跳び `i`

★`σ` は **`K`-代数準同型**でよい(Galois性・全単射性・等長性のどれも要らない)。 -/

section Jump

variable {K L : Type*} [Field K] [NormedField L] [IsUltrametricDist L] [Algebra K L]

/-- ★★**`σ` 版の主結果** —— `‖σx − x‖ · ‖π‖ = ‖σπ − π‖ · ‖x − a_0‖`。

★`σ` が `π` に何をするかだけで `x` 全体への作用が決まる、というのが要点である。 -/
theorem norm_algHom_sub_mul_norm_eq {n : ℕ} {π : L} {a : ℕ → K} (σ : L →ₐ[K] L)
    (hn : 0 < n) (hπ1 : ‖π‖ < 1) (hu : ‖σ π - π‖ < ‖π‖)
    (hchar : ∀ j, 0 < j → j < n → ‖(j : L)‖ = 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    ‖σ (digitSum n π a) - digitSum n π a‖ * ‖π‖
      = ‖σ π - π‖ * ‖digitSum n π a - algebraMap K L (a 0)‖ := by
  rw [map_digitSum]
  exact norm_digitSum_sub_mul_norm_eq hn hπ1 hu hchar hval

/-- ★★**退化検査 / 非空虚性** —— `σ` が `π` を動かすなら、`σ` に固定される `x` は
`K` に入る(`x = a_0`)。

★これは「`L^σ = K`」の桁展開版であり、主結果が**空虚でない**ことの証拠である
(主結果の等式の右辺が `0` になるのは `x ∈ K` のときに限る、と言っている)。
★★分岐理論も Galois 理論も使わずに `L^σ = K` の片側が出る。 -/
theorem eq_algebraMap_digit_zero_of_algHom_eq {n : ℕ} {π : L} {a : ℕ → K} (σ : L →ₐ[K] L)
    (hn : 0 < n) (hπ1 : ‖π‖ < 1) (hu : ‖σ π - π‖ < ‖π‖) (hσπ : σ π ≠ π)
    (hchar : ∀ j, 0 < j → j < n → ‖(j : L)‖ = 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hfix : σ (digitSum n π a) = digitSum n π a) :
    digitSum n π a = algebraMap K L (a 0) := by
  have h := norm_algHom_sub_mul_norm_eq (a := a) σ hn hπ1 hu hchar hval
  rw [hfix, sub_self, norm_zero, zero_mul] at h
  have hne : ‖σ π - π‖ ≠ 0 := norm_ne_zero_iff.mpr (sub_ne_zero.mpr hσπ)
  have hz : ‖digitSum n π a - algebraMap K L (a 0)‖ = 0 := by
    rcases mul_eq_zero.mp h.symm with h1 | h2
    · exact absurd h1 hne
    · exact h2
  exact sub_eq_zero.mp (norm_eq_zero.mp hz)

/-- ★★★**持ち場の目標** —— 跳び `i`(`‖σπ − π‖ = ‖π‖^{i+1}`、すなわち
`σ ∈ G_i ∖ G_{i+1}`)のとき

`‖x − a_0‖ = ‖π‖^{−i} · ‖σx − x‖`.

★★**1 段の定数は `|π_L|^{−i}` である**(等号なので、これが最良)。
★`i ≥ 1`(= `σ ∈ G_1`、暴分岐)を仮定している。`i = 0` だと `‖σπ−π‖ = ‖π‖` で
`hu` が壊れる(そのときは `L/K` が順分岐で、木の tame 側が担当する)。 -/
theorem norm_sub_digit_zero_eq_zpow_mul {n i : ℕ} {π : L} {a : ℕ → K} (σ : L →ₐ[K] L)
    (hn : 0 < n) (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < n → ‖(j : L)‖ = 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    ‖digitSum n π a - algebraMap K L (a 0)‖
      = ‖π‖ ^ (-(i : ℤ)) * ‖σ (digitSum n π a) - digitSum n π a‖ := by
  have hu : ‖σ π - π‖ < ‖π‖ := by
    rw [hbreak]
    calc ‖π‖ ^ (i + 1) < ‖π‖ ^ 1 := by
          exact pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 (by omega)
      _ = ‖π‖ := pow_one _
  have hmain := norm_algHom_sub_mul_norm_eq (a := a) σ hn hπ1 hu hchar hval
  rw [hbreak] at hmain
  have hne : ‖π‖ ≠ 0 := ne_of_gt hπ0
  have hpow : ‖π‖ ^ (i + 1) = ‖π‖ ^ (i : ℤ) * ‖π‖ := by
    rw [zpow_natCast, ← pow_succ]
  rw [hpow] at hmain
  have hcancel : ‖σ (digitSum n π a) - digitSum n π a‖
      = ‖π‖ ^ (i : ℤ) * ‖digitSum n π a - algebraMap K L (a 0)‖ :=
    mul_right_cancel₀ hne (by rw [hmain]; ring)
  rw [hcancel, ← mul_assoc, ← zpow_add₀ hne]
  simp

/-- ★★★**`d(x, K) = ‖π‖^{−i} · ‖σx − x‖`** —— 最良近似込みの形。

`y = a_0 ∈ K` が実際に `‖π‖^{−i}·‖σx − x‖` を達成し(第 1 主張)、
かつ **どの `c ∈ K` もそれより近くならない**(第 2 主張)。
★★つまり `‖π‖^{−i}` は 1 段の**最良の定数**である。

★★★**次の波への警告**: これを塔に沿って繰り返すと定数は **掛かり算**で積み上がる
(超距離は和を `max` に潰すが積は潰さない)。★ファイル冒頭の
「1 段が一様なら Ax が出るは偽」を必ず読むこと。 -/
theorem exists_norm_sub_algebraMap_eq_zpow_mul {n i : ℕ} {π : L} {a : ℕ → K} (σ : L →ₐ[K] L)
    (hn : 0 < n) (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < n → ‖(j : L)‖ = 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    ∃ y : K, ‖digitSum n π a - algebraMap K L y‖
        = ‖π‖ ^ (-(i : ℤ)) * ‖σ (digitSum n π a) - digitSum n π a‖ ∧
      ∀ c : K, ‖π‖ ^ (-(i : ℤ)) * ‖σ (digitSum n π a) - digitSum n π a‖
        ≤ ‖digitSum n π a - algebraMap K L c‖ := by
  refine ⟨a 0, norm_sub_digit_zero_eq_zpow_mul σ hn hπ0 hπ1 hi hbreak hchar hval, fun c => ?_⟩
  rw [← norm_sub_digit_zero_eq_zpow_mul σ hn hπ0 hπ1 hi hbreak hchar hval]
  exact norm_digitSum_sub_digit_zero_le hn hπ0 hπ1 hval a c

end Jump

/-! ## §4 `.src`

★本ファイルは `AxSenTate` の入力であって、原典が独立に立てた項目ではない
(冒頭「逸脱の記録 (3)」)。`SenLemma.lean` と同じ項目を指す。 -/

def norm_pow_sub_pow_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def nnnorm_sum_digit_eq_sup.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_digitSum_sub_digit_zero_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_digitSum_sub_mul_norm_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_algHom_sub_mul_norm_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def eq_algebraMap_digit_zero_of_algHom_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_sub_digit_zero_eq_zpow_mul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_eq_zpow_mul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 使っている公理の一覧 -/

#print axioms norm_eq_of_norm_sub_lt
#print axioms norm_pow_succ_sub_pow_succ_le
#print axioms norm_pow_sub_pow_eq
#print axioms zpow_mul_pow_ne
#print axioms nnnorm_sum_mul_eq_mul_nnnorm_sum
#print axioms digitSum
#print axioms digitSum_sub_digit_zero
#print axioms map_digitSum
#print axioms nnnorm_sum_digit_eq_sup
#print axioms norm_digitSum_sub_digit_zero_le
#print axioms norm_digitSum_sub_mul_norm_eq
#print axioms norm_algHom_sub_mul_norm_eq
#print axioms eq_algebraMap_digit_zero_of_algHom_eq
#print axioms norm_sub_digit_zero_eq_zpow_mul
#print axioms exists_norm_sub_algebraMap_eq_zpow_mul

end ABC3.Found.PGC
