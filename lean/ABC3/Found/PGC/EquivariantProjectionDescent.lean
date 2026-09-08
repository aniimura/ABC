import ABC3.Found.PGC.DeepDescentPairDirect

/-!
# [pGC] 層間の打ち消しの正体は ★**Γ-同変な射影** —— 層の損失は**足し算でなく `max`**

配られた持ち場は

> `PairDirect.AxLemmaGraded K` —— 深さ `k` の `x` に `d(x,K) ≤ (∏_{i=1}^{k} axDecay p i)·ε`。
> 直前の波が測った負の情報: `ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉) ⊃ ℚ₃(ζ₃)` で層ごとの sharp な限界を
> **素朴に足すと `8 + 6 = 14`** で対の予算 `12` に届かない。しかし実測の最悪は `8`。
> ★★**2 層の損失は同時に極値を取れない。ここを言うのが仕事。**

であった。★★**打ち消しの機構を特定した。それは「正規化した跡」である** —— ただし
**降下に使うのではなく、`ε` を伸ばさないために使う**。以下、測定を先に書く。

## ★★★測定 1 —— 何と何が打ち消すのか(★厳密整数演算、`ℤ[π]/(g)` 上)

`M = ℚ₃(ζ₂₇) ⊃ E = ℚ₃(ζ₉) ⊃ K = ℚ₃(ζ₃)`、`v_M`(`v_M(π)=1`, `e_M = 18`)で測る
(スクリプト `.../scratchpad/grad/three.py`, `trnorm.py`, `maxloss.py`。p 進の丸め無し):

| 量 | 記号 | 実測 |
|---|---|---|
| 上の層 `M/E` の sharp な損失(= 跳び) | `t` | ★`8` |
| 下の層 `E/K` の sharp な損失(跳び `2` の `v_M` 換算) | `s` | ★`6` |
| `σ − 1` の最小の伸び(= 第 1 跳び) | `m` | `2` |
| ★★**正規化した跡 `P = (1/p)Tr_{M/E}` の作用素ノルム** | `Q = ‖P‖` | ★`3^{2/18}`(欠損 `γ = 2`) |
| 全体 `M → K` の損失の最大(乱択 + 山登り) | `L` | ★`8` |

★`γ = p − 1` は測った 4 本の層すべてで一致した
(`ℚ₃(ζ₂₇)/ℚ₃(ζ₉)`: `2`、`ℚ₃(ζ₈₁)/ℚ₃(ζ₂₇)`: `2`、`ℚ₂(ζ₈)/ℚ₂(ζ₄)`: `1`、`ℚ₂(ζ₁₆)/ℚ₂(ζ₈)`: `1`)。

## ★★★測定 2 —— 3 通りの合成を同じ塔で比べる(★これが本ファイルの動機)

| 合成の仕方 | `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` | `ℚ₃(ζ₈₁)/ℚ₃(ζ₉)` |
|---|---|---|
| 素朴な telescope `t + s` | `14` ★予算 `12` を超える | `50` ★予算 `36` を超える |
| 収縮つき `max(t, s, t+s−m)`(★下記 §0.1) | `12` ちょうど入る | `42` ★まだ超える |
| ★★**射影つき `max(t+γ, s+γ)`**(★本ファイル) | ★`10` **入る** | ★`28` **入る** |
| 実測の最大 `L` | `8` | `26` |
| 対の予算 | `12` | `36` |

★★★**射影つきの合成は円分塔の全体(`p ≤ 7`, `n ≤ 7`, すべての `k`)で予算に収まる**
(`.../scratchpad/grad/reach2.py`, `reach.py` で 0 件の破れ)。閉じた形は

  ★`c_k = p^{n−1} − 1 + (p−1)·p^{k−2}`  (`k ≥ 2`)、  予算 `= Σ_{i=1}^{k} p^{n−i}`

で、`(p−1)p^{k−2} < p^{k−1} ≤ p^{n−2}` なので**常に余裕がある**。
★素朴な telescope は `k ≥ 2` かつ `(p,n) ≠ (2,3)` で破れる(同スクリプト)。

## ★★★測定 3 —— **なぜ `max` になるのか**(機構の 1 行)

上の層を降りた `x'` を**そのまま**使うと、次の層で使う `ε` が `t` だけ膨らむ
(`UltraCore.norm_smul_sub_self_le_of_norm_sub`)。これが `t + s` の出どころである。
★代わりに `x'' := P x`(`P = (1/p)Tr_{M/E}`、★`Γ = Gal(M/K)` 同変)を取ると:

1. ★`P` は `E` 上恒等なので `x − P x = (x − x') − P(x − x')` ⇒ `‖x − x''‖ ≤ max(1,Q)·‖x − x'‖`。
   ★**同変な射影は sharp な近似元より悪くならない**(`x'` は存在だけ使い、捨てる)。
2. ★★同変性から `σ x'' − x'' = P(σ x − x)` ⇒ `‖σx'' − x''‖ ≤ Q·ε`。
   ★★**`ε` は `t` ではなく `Q` しか膨らまない。これが「層間の打ち消し」の正体である。**

`Q ≤ 1` なら `ε` は 1 ミリも膨らまず、`k` 段の塔の損失は**1 段と同じ**になる
(`EquivProj.towerBudget_le_of_le_one`)。円分塔では `Q = p^{(p−1)/e}` とわずかに `1` を超え、
その分だけが `k` 段で積もる(`towerBudget_le`: `c·q^{k−1}`)。

### ★0.1 収縮つきの合成(先に試して**足りなかった**道。記録として残す)

`σ − 1` は `M` 全体で `‖σz − z‖ ≤ p^{−m/e}‖z‖`(`m` = 第 1 跳び)なので、
`x'` をそのまま使っても `ε` は `t` ではなく `t − m` しか膨らまない。これで
`max(t, s, t+s−m)` が出る(`ℚ₃(ζ₂₇)` では `12` でちょうど予算)。
★しかし `ℚ₃(ζ₈₁)/ℚ₃(ζ₉)` で `42 > 36` となり**閉じない**ので、本ファイルは採らない。

## 本ファイルが出したもの(すべて `sorry` 0)

| 宣言 | 内容 |
|---|---|
| ★★★`EquivProj.exists_of_proj_two_descent` | 抽象核(同変射影による 2 層の合成、★`max` になる) |
| ★★`EquivProj.towerBudget` / `_nonneg` / `_le` / `_le_of_le_one` | 塔の予算(★`k` 段でも `c·q^{k−1}`、`q ≤ 1` なら `c`) |
| ★★★`EquivProj.exists_of_proj_tower` | 抽象核(★`k` 段の塔を同変射影で降りる) |
| ★★`Zeta27Proj.composed_eq` | 実測値を入れると合成は ★`3^{10/18}` |
| ★★`Zeta27Proj.fits_pair_budget` | ★`3^{10/18}` は対の予算 `3^{12/18}` に**入る** |
| `Zeta27Proj.composed_lt_naive` | 素朴な `3^{14/18}` より真に小さい |

★抽象核は**分岐・付値・Galois の語彙が 1 語も出ない**(超距離群 + モノイド作用 + 射影のみ)。

## ★★★埋まっていないもの(★正確に。`AxLemmaGraded` は**閉じていない**)

抽象核の仮説のうち、`p` 進の具体層で**まだ形式化していない**のは次の 2 つである。

1. ★**同変射影 `P` の構成**: `P = (1/p)·Tr_{M/E}`。`hPsub`/`hPid`/`hPeq` は
   跡の加法性・`Tr_{M/E}|_E = p·id`・`N ◁ Γ` から出るが、本ファイルは
   `K.closure` の上でそれを作っていない(中間体を作ると #59/#69 の境界に当たるため)。
2. ★★**作用素ノルム `‖P‖ = p^{(p−1)/e}`**: Serre `Corps Locaux` V §3 Lemme 4
   `Tr(𝔭_M^n) = 𝔭_E^{⌊(n+d)/e⌋}` と `d_{M/E} = v_M(p)`(円分層では等号)から
   `‖P‖ ≤ p^{(v_M(p) + p − 1 − d)/e}` が出る。★木にある
   `RamificationJumpBound` の `(p−1)i ≤ e_L` は `d ≤ v_M(p) + p − 1` を与える**逆向き**なので、
   一般の層では `d ≥ v_M(p)` は**言えない**(★下記の逸脱 3)。
   ★一般の層で必要なのは「`t` が小さいときは `Q` が大きい」という**トレードオフ**であり、
   両者を同時に極値に取れないことを言う 1 本の補題が次の 1 点である。

★★したがって残りは依然として **1 点**だが、その形は前の波より 1 段具体的になった:

  ★★★`d_{M/E}` と上の層の跳び `t` の**トレードオフ**(`t + max(0, v(p) + p − 1 − d) ≤ v(p)` 型)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 抽象核の `htop`(上の層の降下)は**全群 `Γ` の変位**で仮定した。原典・木の
   `CyclicJumpNorm` は部分群 `Gal(M/E)` の変位で足りるので、これは**仮説を強める**方向だが、
   合成側では `Γ` の変位しか持っていないので消費に支障は無い(むしろ供給が楽)。
2. 測定はすべて Lean の外(整数多項式の厳密演算)で行った。Lean 側に在るのは
   その数値から出る**実数の不等式**だけで、`ℚ₃(ζ₂₇)` 自体は構成していない
   (`DeepDescentPairDirect` の逸脱 2 と同じ扱い)。
3. `γ = p − 1`(= `‖P‖` の欠損)は円分層 4 本の**実測**であって、一般の層で証明していない。
   一般の層では `γ = max(0, v_M(p) + p − 1 − d_{M/E})` が上界である。
4. `NormalizedTraceDescent` は「跡の道は `p ≥ 3` で閉じない」と結論しているが、
   ★**それは跡を「降下」に使った場合**である(指数が sharp の `(p−1)` 倍になる)。
   本ファイルは跡を**同変性のためだけ**に使い、降下は sharp な層の限界で行う。
   ★同じ道具でも使う場所が違うので、両者は矛盾しない。
-/

namespace ABC3.Found.PGC

/-! ## §1 抽象核 —— **同変な射影**で 2 層を合成すると損失は `max` になる -/

namespace EquivProj

variable {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]
variable {G : Type*} [Monoid G] [DistribMulAction G M]

/-- ★★★**抽象核**(分岐・付値・Galois が 1 語も出ない)。

超距離群 `M` にモノイド `G` が作用し、部分集合 `T`(作業域)、`E`(中間層)、`K`(着地)と
**`G` 同変な射影** `P : M → M`(像は `E`、`E` 上恒等、`‖P z‖ ≤ Q‖z‖`)が在るとする。
`T` から `E` への降下が `A` 倍、`E` から `K` への降下が `B` 倍でできるなら、
`T` から `K` への降下は ★`max (A · max 1 Q) (B · Q)` 倍でできる。

★★**要点は 2 つ**:
* `P x` は sharp な近似元 `x'` より悪くない —— `x − P x = (x − x') − P(x − x')`
  (`P` が `E` 上恒等だから)。★`x'` は**存在だけ**使って捨てる。
* ★同変性 `σ (P x) − P x = P (σ x − x)` から `ε` は `A` 倍ではなく **`Q` 倍**しか伸びない。

`Q ≤ 1` なら結論は `max A B`、すなわち ★**2 層の損失は足し算にならず `max` になる**。
素朴な合成(`UltraCore.norm_smul_sub_self_le_of_norm_sub` を使うもの)は `A · B` を払う。 -/
theorem exists_of_proj_two_descent
    {E K T : Set M} {P : M → M} {A B Q ε : ℝ}
    (hPsub : ∀ z w, P (z - w) = P z - P w)
    (hPE : ∀ z, P z ∈ E) (hPid : ∀ y ∈ E, P y = y)
    (hPeq : ∀ (g : G) (z : M), P (g • z) = g • P z)
    (hPQ : ∀ z, ‖P z‖ ≤ Q * ‖z‖) (hQ0 : 0 ≤ Q)
    (htop : ∀ z ∈ T, ∀ δ : ℝ, 0 ≤ δ → (∀ g : G, ‖g • z - z‖ ≤ δ) →
      ∃ z' ∈ E, ‖z - z'‖ ≤ A * δ)
    (hbot : ∀ y ∈ E, ∀ δ : ℝ, 0 ≤ δ → (∀ g : G, ‖g • y - y‖ ≤ δ) →
      ∃ a ∈ K, ‖y - a‖ ≤ B * δ)
    {x : M} (hxT : x ∈ T) (hε : 0 ≤ ε) (hx : ∀ g : G, ‖g • x - x‖ ≤ ε) :
    ∃ a ∈ K, ‖x - a‖ ≤ max (A * max 1 Q) (B * Q) * ε := by
  obtain ⟨x', hx'E, hx'⟩ := htop x hxT ε hε hx
  have hAε : (0:ℝ) ≤ A * ε := le_trans (norm_nonneg _) hx'
  have hsplit : x - P x = (x - x') - P (x - x') := by
    rw [hPsub, hPid x' hx'E]; abel
  have hxy : ‖x - P x‖ ≤ max 1 Q * (A * ε) := by
    rw [hsplit, sub_eq_add_neg]
    refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ ?_)
    · exact le_trans hx' (le_mul_of_one_le_left hAε (le_max_left _ _))
    · rw [norm_neg]
      exact le_trans (hPQ _) (mul_le_mul (le_max_right _ _) hx' (norm_nonneg _)
        (le_trans zero_le_one (le_max_left _ _)))
  have hy : ∀ g : G, ‖g • P x - P x‖ ≤ Q * ε := by
    intro g
    have hg : g • P x - P x = P (g • x - x) := by rw [hPsub, hPeq]
    rw [hg]
    exact le_trans (hPQ _) (mul_le_mul_of_nonneg_left (hx g) hQ0)
  obtain ⟨a, haK, ha⟩ := hbot (P x) (hPE x) (Q * ε) (mul_nonneg hQ0 hε) hy
  refine ⟨a, haK, ?_⟩
  have hab : x - a = (x - P x) + (P x - a) := by abel
  rw [hab]
  refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ ?_)
  · calc ‖x - P x‖ ≤ max 1 Q * (A * ε) := hxy
      _ = (A * max 1 Q) * ε := by ring
      _ ≤ max (A * max 1 Q) (B * Q) * ε := mul_le_mul_of_nonneg_right (le_max_left _ _) hε
  · calc ‖P x - a‖ ≤ B * (Q * ε) := ha
      _ = (B * Q) * ε := by ring
      _ ≤ max (A * max 1 Q) (B * Q) * ε := mul_le_mul_of_nonneg_right (le_max_right _ _) hε

/-- ★`Q ≤ 1`(同変射影がノルムを増やさない)なら合成の予算は ★`max A B` —— **足し算にならない**。 -/
theorem composed_le_max {A B Q : ℝ} (hB : 0 ≤ B) (hQ1 : Q ≤ 1) :
    max (A * max 1 Q) (B * Q) ≤ max A B := by
  refine max_le ?_ ?_
  · rw [max_eq_left hQ1, mul_one]; exact le_max_left _ _
  · exact le_trans (mul_le_of_le_one_right hB hQ1) (le_max_right _ _)

/-- ★同じ仮定で、**素朴な合成 `A · B` より良い**(`1 ≤ A`, `1 ≤ B`)。 -/
theorem composed_le_mul {A B Q : ℝ} (hA : 1 ≤ A) (hB : 1 ≤ B) (hQ1 : Q ≤ 1) :
    max (A * max 1 Q) (B * Q) ≤ A * B := by
  refine le_trans (composed_le_max (le_trans zero_le_one hB) hQ1) ?_
  exact max_le (le_mul_of_one_le_right (le_trans zero_le_one hA) hB)
    (le_mul_of_one_le_left (le_trans zero_le_one hB) hA)

/-! ## §2 抽象核 —— **`k` 段の塔**(★段数が増えても予算は `c · q^{k−1}`) -/

/-- ★★**塔の予算**。`A i` は第 `i` 層の降下の倍率、`Q i` は第 `i` 層への同変射影のノルム。

  `towerBudget A Q 0 = 1`、`towerBudget A Q 1 = A 0`、
  `towerBudget A Q (k+2) = max (A (k+1) · max 1 (Q (k+1))) (towerBudget A Q (k+1) · Q (k+1))`

★**最下段(`k = 1`)は射影を使わない** —— 着地は `F 0` そのものなので `A 0` だけ払う。
これを `max (A 0 · max 1 (Q 0)) (Q 0)` にすると `F 0` への射影のノルムを余計に払う
(円分塔ではそれが `Q 0 = p^{v(p)/e}` と大きく、勘定が壊れる)。 -/
noncomputable def towerBudget (A Q : ℕ → ℝ) : ℕ → ℝ
  | 0 => 1
  | 1 => A 0
  | (k+2) => max (A (k+1) * max 1 (Q (k+1))) (towerBudget A Q (k+1) * Q (k+1))

theorem towerBudget_nonneg {A Q : ℕ → ℝ} (hA : ∀ i, 0 ≤ A i) : ∀ k, 0 ≤ towerBudget A Q k
  | 0 => zero_le_one
  | 1 => hA 0
  | (k+2) => le_trans (mul_nonneg (hA (k+1)) (le_trans zero_le_one (le_max_left _ _)))
      (le_max_left _ _)

/-- ★★★**塔の予算は段数に対して `c · q^{k−1}`**(★`c^k` ではない)。

各層の降下が `c` 倍まで、各射影のノルムが `q` 倍までなら、`k` 段の塔全体でも
`c · q^{k−1}` で済む。★層の倍率 `c` は**掛け算されず、`max` で 1 回しか払わない**。 -/
theorem towerBudget_le {A Q : ℕ → ℝ} {c q : ℝ} (hc : 1 ≤ c) (hq : 1 ≤ q)
    (hQ0 : ∀ i, 0 ≤ Q i) (hA : ∀ i, A i ≤ c) (hQ : ∀ i, Q i ≤ q) :
    ∀ k, towerBudget A Q k ≤ c * q ^ (k - 1)
  | 0 => by simpa [towerBudget] using hc
  | 1 => by simpa [towerBudget] using hA 0
  | (k+2) => by
      have ih := towerBudget_le hc hq hQ0 hA hQ (k+1)
      have hc0 : (0:ℝ) ≤ c := le_trans zero_le_one hc
      have hq0 : (0:ℝ) ≤ q := le_trans zero_le_one hq
      have hqk : (1:ℝ) ≤ q ^ k := one_le_pow₀ hq
      have e2 : (k+1) - 1 = k := rfl
      rw [e2] at ih
      have hqq : q ≤ q ^ (k+1) := by rw [pow_succ]; nlinarith
      show max (A (k+1) * max 1 (Q (k+1))) (towerBudget A Q (k+1) * Q (k+1)) ≤ c * q ^ (k+1)
      refine max_le ?_ ?_
      · have h1 : max 1 (Q (k+1)) ≤ q := max_le hq (hQ _)
        calc A (k+1) * max 1 (Q (k+1)) ≤ c * q :=
              mul_le_mul (hA _) h1 (le_trans zero_le_one (le_max_left _ _)) hc0
          _ ≤ c * q ^ (k+1) := mul_le_mul_of_nonneg_left hqq hc0
      · calc towerBudget A Q (k+1) * Q (k+1) ≤ (c * q ^ k) * q :=
              mul_le_mul ih (hQ _) (hQ0 _) (mul_nonneg hc0 (pow_nonneg hq0 k))
          _ = c * q ^ (k+1) := by rw [pow_succ]; ring

/-- ★★★**同変射影がノルムを増やさない(`Q ≤ 1`)なら、`k` 段の塔は 1 段と同じ値段**。
★これが「層の損失は足し算にならない」の最終形である。 -/
theorem towerBudget_le_of_le_one {A Q : ℕ → ℝ} {c : ℝ} (hc : 1 ≤ c)
    (hQ0 : ∀ i, 0 ≤ Q i) (hA : ∀ i, A i ≤ c) (hQ : ∀ i, Q i ≤ 1) (k : ℕ) :
    towerBudget A Q k ≤ c := by
  simpa using towerBudget_le hc le_rfl hQ0 hA hQ k

/-- ★★★**抽象核 2**(分岐・付値・Galois が 1 語も出ない)。

減少列 `F : ℕ → Set M` と、各 `F i` への `G` 同変な射影 `P i`(ノルム `≤ Q i`)が在り、
1 段の降下 `F (i+1) → F i` が `A i` 倍でできるなら、
`F k` の元は `towerBudget A Q k` 倍で `F 0` まで降りる。

★`Q i ≤ 1` なら `towerBudget A Q k ≤ max_i (A i)`(`towerBudget_le_of_le_one`)、すなわち
★★**`k` 段降りても 1 段分しか払わない**。 -/
theorem exists_of_proj_tower
    {F : ℕ → Set M} {P : ℕ → M → M} {A Q : ℕ → ℝ}
    (hQ0 : ∀ i, 0 ≤ Q i)
    (hPsub : ∀ i z w, P i (z - w) = P i z - P i w)
    (hPE : ∀ i z, P i z ∈ F i) (hPid : ∀ i, ∀ y ∈ F i, P i y = y)
    (hPeq : ∀ (i : ℕ) (g : G) (z : M), P i (g • z) = g • P i z)
    (hPQ : ∀ i z, ‖P i z‖ ≤ Q i * ‖z‖)
    (hdesc : ∀ i, ∀ z ∈ F (i+1), ∀ δ : ℝ, 0 ≤ δ → (∀ g : G, ‖g • z - z‖ ≤ δ) →
      ∃ z' ∈ F i, ‖z - z'‖ ≤ A i * δ) :
    ∀ (k : ℕ) (x : M), x ∈ F k → ∀ ε : ℝ, 0 ≤ ε → (∀ g : G, ‖g • x - x‖ ≤ ε) →
      ∃ a ∈ F 0, ‖x - a‖ ≤ towerBudget A Q k * ε := by
  intro k
  induction k with
  | zero => intro x hx ε hε _; exact ⟨x, hx, by simpa [towerBudget] using hε⟩
  | succ k ih =>
    rcases k with _ | k
    · intro x hx ε hε hxε
      obtain ⟨x', hx', h⟩ := hdesc 0 x hx ε hε hxε
      exact ⟨x', hx', by simpa [towerBudget] using h⟩
    · intro x hx ε hε hxε
      exact exists_of_proj_two_descent (E := F (k+1)) (K := F 0) (T := F (k+2))
        (hPsub (k+1)) (hPE (k+1)) (hPid (k+1)) (hPeq (k+1)) (hPQ (k+1)) (hQ0 (k+1))
        (hdesc (k+1)) (fun y hy δ hδ hyδ => ih y hy δ hδ hyδ) hx hε hxε

end EquivProj

/-! ## §3 数値層 —— `ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉) ⊃ ℚ₃(ζ₃)` に実測値を入れる -/

namespace Zeta27Proj

/-- ★★**抽象核に実測値を入れた結果**。`A = 3^{8/18}`(上の層の跳び `8`)、
`B = 3^{6/18}`(下の層の跳び `2` の `v_M` 換算)、`Q = 3^{2/18}`(正規化した跡の欠損 `γ = 2`)
を `EquivProj.exists_of_proj_two_descent` の結論に入れると、合成の倍率は ★`3^{10/18}` になる。

★素朴な telescope は `3^{8/18} · 3^{6/18} = 3^{14/18}` で、対の予算 `3^{12/18}` を超える
(`Zeta27Pair.naive_tower_sum_exceeds_pair_budget`)。 -/
theorem composed_eq :
    max ((3:ℝ) ^ ((8:ℝ)/18) * max 1 ((3:ℝ) ^ ((2:ℝ)/18))) ((3:ℝ) ^ ((6:ℝ)/18) * (3:ℝ) ^ ((2:ℝ)/18))
      = (3:ℝ) ^ ((10:ℝ)/18) := by
  have h3 : (0:ℝ) < 3 := by norm_num
  have h1 : (1:ℝ) ≤ (3:ℝ) ^ ((2:ℝ)/18) := Real.one_le_rpow (by norm_num) (by norm_num)
  rw [max_eq_right h1, ← Real.rpow_add h3, ← Real.rpow_add h3]
  norm_num

/-- ★★★**射影つきの合成は対の予算に入る**(`J = 10 ≤ 12`)。
★これが `Zeta27Pair.naive_tower_sum_exceeds_pair_budget`(`J = 14` は入らない)への回答である。 -/
theorem fits_pair_budget :
    (3:ℝ) ^ ((10:ℝ)/18) ≤ ∏ i ∈ Finset.Icc 1 2, axDecay 3 i := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h := (PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := 3) (k' := 0) (d := 2)
    (J := 10) (E := 18) (by norm_num)).mpr (by norm_num)
  simpa using h

/-- ★合成は素朴な telescope より**真に小さい**(`10 < 14`)。 -/
theorem composed_lt_naive : (3:ℝ) ^ ((10:ℝ)/18) < (3:ℝ) ^ ((14:ℝ)/18) :=
  Real.rpow_lt_rpow_left_iff (by norm_num) |>.mpr (by norm_num)

/-- ★★`ℚ₃(ζ₈₁)/ℚ₃(ζ₉)`(`e_M = 54`、上の層の跳び `26`、下の層 `24`、欠損 `γ = 2`)でも入る:
合成は `3^{28/54}`、対の予算は `3^{36/54}`。★収縮つきの合成(`42`)では入らなかった段である。 -/
theorem zeta81_fits_pair_budget :
    (3:ℝ) ^ ((28:ℝ)/54) ≤ ∏ i ∈ Finset.Icc 1 2, axDecay 3 i := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h := (PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := 3) (k' := 0) (d := 2)
    (J := 28) (E := 54) (by norm_num)).mpr (by norm_num)
  simpa using h

end Zeta27Proj

/-! ## §4 `.src`(原典の対応箇所) -/

namespace EquivProj

def exists_of_proj_two_descent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def towerBudget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def towerBudget_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def exists_of_proj_tower.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end EquivProj

namespace Zeta27Proj

def composed_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def fits_pair_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def zeta81_fits_pair_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Zeta27Proj
end ABC3.Found.PGC

section Audit
open ABC3.Found.PGC
#print axioms EquivProj.exists_of_proj_two_descent
#print axioms EquivProj.composed_le_max
#print axioms EquivProj.composed_le_mul
#print axioms EquivProj.towerBudget_nonneg
#print axioms EquivProj.towerBudget_le
#print axioms EquivProj.towerBudget_le_of_le_one
#print axioms EquivProj.exists_of_proj_tower
#print axioms Zeta27Proj.composed_eq
#print axioms Zeta27Proj.fits_pair_budget
#print axioms Zeta27Proj.composed_lt_naive
#print axioms Zeta27Proj.zeta81_fits_pair_budget
end Audit
