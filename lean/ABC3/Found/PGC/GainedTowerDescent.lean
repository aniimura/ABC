import ABC3.Found.PGC.JumpDefectTradeoff

/-!
# [pGC] ★★★★`p ≥ 5` の窓は**塞がった** —— 緩かったのは定理でなく上界であり、欠けていたのは `(τ−1)^p` である

配られた持ち場は

> ★★★`p ≥ 5` で上界を鋭くする —— `AxLemmaGraded` に残った唯一の数学。
> 緩いのは「収縮の枝 `C`」と「射影の枝 `P`」を別々に上から押さえて `max` を取っている点ではないか。

であった。★★**結論を先に書く: 上界は鋭くなり、`p ≥ 5` の窓は消えた。**
★ただし**本体の見立て(2 つの枝を同時に極値にできない)は外れである**。緩かったのは
「枝の組み合わせ方」ではなく、★**そもそも 2 つの枝のどちらも、古典の道具を 1 つ使い落としていた**点である。

## ★★★測定 0 —— 反例族は**本物**だった(字面の検算)

`JumpDefectTradeoff.witness_gap_exceeds` の `(p, e_K, S, T) = (5, 3, 2, 17)` は
実際に許容な塔の跳びである(上付き `u_1 = 2`, `u_2 = u_1 + e_K = 5`、Hasse–Arf を満たし、
層の上界 `8 ≤ 15`, `68 ≤ 75` も満たす)。★**偽なのは反例ではなく、そこで測っていた上界の方**である。
旧ファイルの結論「`p ≥ 5` では閉じない」は、★**旧ファイルの上界については真**であり、
本ファイルはその上界を取り替える。

## ★★★★測定 1 —— 使い落としていた古典の道具は `(τ^p − 1) ≡ (τ−1)^p  (mod p)`

`Γ = Gal(M/K) ≅ ℤ/p^k`、`τ` を生成元とする。可換環 `ℤ[Γ]` で `D := τ − 1` と置くと

  ★`τ^p − 1 = (1 + D)^p − 1 = D^p + Σ_{i=1}^{p−1} C(p,i)·D^i`,   `p ∣ C(p,i)` (`1 ≤ i ≤ p−1`)

である。`τ^{p^{j−1}}` は `Γ_{t_j} ∖ Γ_{t_j+1}` に入るので `v_M((τ^{p^{j−1}}−1)z) ≥ v_M(z) + t_j`
(`t_j` は `Γ` の `j` 番目の下付き跳び、`v_M` の目盛り)。したがって
`g_j := v_M((τ^{p^j}−1)x)` は

  ★★`g_j ≥ min(g_{j−1} + (p−1)t_j , e_M + g_{j−1}) = g_{j−1} + (p−1)t_j`
     (`(p−1)t_j ≤ e_{F_j} ≤ e_M` は sharp な層の上界そのもの)

を満たし、`v_M((τ−1)x) ≥ c` なら

  ★★★`v_M((τ^{p^{k−1}} − 1)x) ≥ c + (p−1)(t_1 + ⋯ + t_{k−1})`

が出る。★★**「`τ` でほぼ不変」は「上の層の生成元 `τ^{p^{k−1}}` では `(p−1)Σt_i` だけ**得をしている**」
ことを意味する。** 旧ファイルの `contrCost` / `projCost` は `c` のまま上の層に入っており、
この得を 1 度も使っていなかった。

## ★★★測定 2 —— 1 層の降下は**厳密な等式**である

`M/E` を全分岐巡回 `p` 次、跳び `t`(`v_M`)、`π = π_M`、`𝔬_M = 𝔬_E[π]`(Serre I §6 Prop 18)。
`x = Σ_{i=0}^{p−1} c_i π^i` (`c_i ∈ E`) と書くと `σπ = π(1+u)`, `v_M(u) = t` で
`v_M(c_i((σπ)^i − π^i)) = v_M(c_iπ^i) + t` (`1 ≤ i ≤ p−1` は `p` の単元)。
`v_M(c_iπ^i)` は `i` ごとに `mod p` で相異なるので

  ★★`max_{y ∈ E} v_M(x − y) = v_M((σ−1)x) − t`   (★**等式**。不等式ではない)

★これが「層の降下の値段はちょうど跳び `t`」であり、旧ファイルの `projCost` が払っていた
**欠損 `γ = v_M(p) − (p−1)t` は、`(1/p)Tr` を使うから発生していたのであって、
降下そのものには不要**であった。

## ★★★★測定 3 —— 鋭い漸化式(本ファイルの主結果)

`v_M` の目盛りで `G_m := t_1 + ⋯ + t_m`、`A_m := t_m − (p−1)G_{m−1}` と置く。
上の 2 つを組んで `M = F_k → F_{k−1} → ⋯ → K` と降りると

  ★★★`Λ_m = max( A_m , max(0, A_m − t_1) + p·Λ_{m−1} )`,  `Λ_0 = 0`

が出る(`gainedLoss`)。段の内訳:
1. `v_{F_m}((τ^{p^{m−1}}−1)y) ≥ c + (p−1)G_{m−1}` (測定 1 を `F_m` の中で使う)
2. `y' ∈ F_{m−1}` を取ると `v_{F_m}(y − y') = v_{F_m}((τ^{p^{m−1}}−1)y) − t_m ≥ c − A_m` (測定 2)
3. `v_{F_m}((τ−1)y') ≥ min(c, v_{F_m}(y−y') + t_1) = c − max(0, A_m − t_1)`
4. `v_{F_{m−1}}` に直すと目盛りが `1/p` になるので、下の塔の損失は `p·Λ_{m−1}` として効く。

★★**深さ 2 では閉じた形になる**(`gainedLoss_two`):

  ★★★`Λ_2 = max(T, p·S)`   (`S = t_1`, `T = t_2`)

★旧上界との関係は**完全に明示的**である(`projCost_eq_sharp_add_defect`):

  ★★`projCost = max(T, pS) + γ`,  `γ = p²e_K − (p−1)T ≥ 0`

すなわち★**旧上界は鋭い値に「射影の欠損 `γ` を丸ごと余計に足していた」**。
★`p = 5, e_K = 3, S = 2, T = 17` では `γ = 7` で `24 = 17 + 7`。

## ★★★★測定 4 —— 数値(★すべて厳密整数/有理数演算)

| 塔 | `p` | `e_K` | `S` | `T` | 旧 `C` | 旧 `P` | ★鋭い値 | 予算 `B_2` | 既存ファイルの実測 `L` |
|---|---|---|---|---|---|---|---|---|---|
| `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` | 3 | 2 | 2 | 8 | 12 | 10 | ★**8** | 12 | ★`8` ✓ |
| `ℚ₃(ζ₈₁)/ℚ₃(ζ₉)` | 3 | 6 | 8 | 26 | 42 | 28 | ★**26** | 36 | ★`26` ✓ |
| ★旧反例 | 5 | 3 | 2 | 17 | 25 | 24 | ★**17** | 22.5 | (未測定) |

★★**鋭い値は既存ファイルが `L` を明示している 2 本と*厳密に*一致する**
(`EquivariantProjectionDescent` の `L = 8` / `L = 26`)。旧上界は `10` / `28` で
一致していなかった。★これは前波の docstring の予想「実測の最悪値は上の層の跳び `T`」の
**正しい一般形**である: `T` ではなく `max(T, pS)`(2 本とも `T ≥ pS` なので `T` に見えていた)。

★総当たり(`…/scratchpad/sharp/{sharp,rec,rec2,closed,probe}.py`、厳密有理数):
`p ∈ {2,3,5,7,11,13}`、`k ≤ 4`、許容列 **46,108 本で予算超過 0 件**。
最悪比 `Λ_k / B_k` は `k = 1` の **1.0000**(★等号。`axDecay p 1 = p^{1/(p−1)}` が最良定数である
という既存の測定と整合)、`k ≥ 2` では `p = 13, e = 8, t = (1,105,1457,19033)` の `0.9896` が最悪。

★`A_m ≥ 0`・`Λ_k ≤ Σt_j`・`Σt_j ≤ B_k` は別に **23,820 本で破れ 0 件**(`closed.py`)。
★閉じた形も出た(★**4,194 本 → 不一致 0 件**だが本ファイルは証明していない):
`Λ_k = t_k + Σ_{j<k} t_j − t_1·(p^{k−1}−1)/(p−1)`。

## ★★★★測定 5 —— 予算 `B_k` の正体は「跳びの最大値の**総和**」

`B_k = e_K·p(p^k−1)/(p−1)² = Σ_{j=1}^{k} p^j e_K/(p−1) = Σ_{j=1}^{k} (層 `j` の跳びの上限)`
である。★したがって

  ★★★`Λ_k ≤ t_1 + ⋯ + t_k ≤ Σ_{j=1}^{k} p^j e_K/(p−1) = B_k`

を示せばよい(`gainedLoss_le_jumpSum` + `jumpSum_le_geom`)。★**予算にちょうど収まる形**である。
帰納が回る条件は `A_m ≥ 0`、すなわち `(p−1)G_{m−1} ≤ t_m` であり、これは

  ★`p·t_j ≤ t_{j+1}`(すべての `j`)

から出る(`topDefect_nonneg`。`(p−1)G_m = (p−1)G_{m−1} + (p−1)t_m ≤ t_m + (p−1)t_m = p t_m ≤ t_{m+1}`)。
★そして `p·t_j ≤ t_{j+1}` は**二分律のどちらの枝からも出る**(§2):
* 安定域 `t_{j+1} = t_j + p^{j+1}e`: 層の上界 `(p−1)t_j ≤ p^{j+1}e` から直ちに。
* 臨界以下 `u_{j+1} ≥ p·u_j`: `t_j ≤ p^j u_j` から `t_{j+1} ≥ t_j + p(p−1)t_j ≥ p t_j`。

## 本ファイルが出したもの(すべて `sorry` 0)

| 宣言 | 内容 |
|---|---|
| ★★★`GainedDescent.gainedLoss` | 鋭い降下の漸化式(★`ℤ` の算術のみ) |
| ★★★`GainedDescent.gainedLoss_le_jumpSum` | ★`Λ_k ≤ t_1+⋯+t_k`(`A_m ≥ 0` から) |
| ★★`GainedDescent.topDefect_nonneg` | `p t_j ≤ t_{j+1}` ⇒ `A_m ≥ 0` |
| ★★★★`GainedDescent.gainedLoss_fits` | ★**主定理**: `(p−1)²Λ_k ≤ (p^{k+1}−p)e`。★`p` にも `k` にも条件なし |
| ★★`GainedDescent.pmul_le_of_stable` / `_of_below` | 二分律の両枝から `p t_j ≤ t_{j+1}` |
| ★★★`GainedDescent.gainedLoss_two` | `Λ_2 = max(t_2, p t_1)` |
| ★★★★`GainedDescent.sharpTwoCost_fits` | ★**`p ≥ 5` の窓が消える**: 二分律なしで対の予算に入る |
| ★★★`GainedDescent.projCost_eq_sharp_add_defect` | ★旧上界 = 鋭い値 + 欠損 `γ`(緩みの正体) |
| ★★★`GainedDescent.witness_five_*` | 旧反例 `(5,3,2,17)` は `17` で予算 `22.5` にも Ax の定数にも入る |
| ★★`GainedDescent.towerBudget_iff` | 塔の台帳(`pairBudget_iff` の深さ `k` 版) |
| ★★★★`GainedDescent.rpow_le_prod_axDecay` | ★`p^{J/e_M} ≤ ∏_{i∈[1,k]} axDecay p i`(`AxLemmaGraded` の形) |

★抽象核 §1・§2 は **`ℤ` の算術だけ**で、分岐・付値・Galois の語彙が 1 語も出ない。

## ★★★埋まっていないもの(★正確に)

1. ★**降下の値段の模型そのものは Lean の外にある。** `gainedLoss` は `ℤ` 上の `def` であり、
   「この漸化式が実際の降下の損失を上から押さえる」ことは上の測定 1〜3 の議論であって
   Lean の証明ではない。★これは `JumpDefectTradeoff` の `contrCost` / `projCost` と**同じ扱い**である
   (Serre V §3 Lemme 4 / IV §1 Prop 4 / III §6 Prop 13 は木にも mathlib にも無い)。
2. `Λ_k` の**閉じた形** `t_k + Σ_{j<k}t_j − t_1(p^{k−1}−1)/(p−1)` は 4,194 本で
   不一致 0 件を測ったが、証明は書いていない(必要でもない)。
3. 二分律(安定域)そのものは相変わらず仮説である(§2 は「二分律 ⇒ `p t_j ≤ t_{j+1}`」だけを示す)。
4. ★塔が**巡回でない**場合(`ℚ₂(ζ₁₆)/ℚ₂` = `ℤ/2×ℤ/4` など)は測定 1 の
   `(τ^p−1) = (τ−1)^p + p(⋯)` が生成元 1 本では書けないので、本ファイルの外である。
   ★前波の docstring が「逆算値 `8 > T = 7` で危うい」と書いた箇所はここに当たる。
   ★ただし `Λ_2 = max(T, pS)` は `ℚ₂(ζ₁₆)` の逆算値 `8` を**満たす**(`T = 7`, `S = 4`,
   `max(7, 8) = 8`)ので、非巡回でも `max(T, pS)` の形は保たれている可能性が高い。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 抽象核は `ℤ` の不等式で、`p` の素数性を使わない(`2 ≤ p` のみ)。素数性は測定 1 の
   `p ∣ C(p,i)` にだけ効くが、それは `ℤ` の側では仮説 `hstep` / `hlayer` に吸収されている。
2. `JumpDefectTradeoff` の `contrCost` / `projCost` / `pairBudget` は**読むだけ**で書き換えない。
   本ファイルは `sharpTwoCost` を**別の宣言**として置き、両者の差が `γ` であることを
   `projCost_eq_sharp_add_defect` で明示する。★旧ファイルの反例は**旧上界については依然として真**なので
   取り下げない(`witness_gap_exceeds` はそのまま残す)。
-/

namespace ABC3.Found.PGC

namespace GainedDescent

open JumpDefect

/-! ## §1 抽象核 —— 鋭い降下の漸化式(★`ℤ` の算術のみ) -/

/-- `G m = t 1 + ⋯ + t m`。★具体層では `t j` は `Γ = Gal(M/K)` の `j` 番目の下付き跳び
(`v_M` の目盛り)であり、同時に層 `F_j/F_{j−1}` の跳びを `v_{F_j}` で測った値でもある。 -/
def jumpSum (t : ℕ → ℤ) : ℕ → ℤ
  | 0 => 0
  | (m + 1) => jumpSum t m + t (m + 1)

/-- ★★★**鋭い降下の漸化式**。`A_m = t_m − (p−1)·G_{m−1}` は
「上の層の跳び」から「`(τ−1)^p` の得」を引いた残りである。

  `Λ_m = max( A_m , max(0, A_m − t_1) + p·Λ_{m−1} )`

第 1 枝は「上の層で止まる」場合、第 2 枝は「下の塔まで降りる」場合。
`max(0, A_m − t_1)` は下の塔に持ち込む同変性の劣化、`p·Λ_{m−1}` は目盛り換算である。 -/
def gainedLoss (p : ℤ) (t : ℕ → ℤ) : ℕ → ℤ
  | 0 => 0
  | (m + 1) =>
      max (t (m + 1) - (p - 1) * jumpSum t m)
        (max 0 (t (m + 1) - (p - 1) * jumpSum t m - t 1) + p * gainedLoss p t m)

/-- `Σ_{j=1}^k p^j`。★予算 `B_k` の分子。 -/
def geomSum (p : ℤ) : ℕ → ℤ
  | 0 => 0
  | (k + 1) => geomSum p k + p ^ (k + 1)

/-- `(p−1)·Σ_{j=1}^k p^j = p^{k+1} − p`。 -/
theorem geomSum_mul (p : ℤ) : ∀ k, (p - 1) * geomSum p k = p ^ (k + 1) - p
  | 0 => by simp [geomSum]
  | (k + 1) => by
      have IH := geomSum_mul p k
      simp only [geomSum]
      ring_nf
      ring_nf at IH
      linarith [IH]

theorem jumpSum_nonneg {t : ℕ → ℤ} (ht : ∀ j, 0 ≤ t (j + 1)) : ∀ m, 0 ≤ jumpSum t m
  | 0 => le_refl 0
  | (m + 1) => by
      have h1 := jumpSum_nonneg ht m
      have h2 := ht m
      simp only [jumpSum]
      linarith

/-- ★★★**核 1**: `A_m ≥ 0` のもとで `Λ_k ≤ t_1 + ⋯ + t_k`。

★第 2 枝がちょうど `G_m` になるのがこの補題の全部である:
`max(0, A_m − t_1) + p·G_{m−1} ≤ A_m + p G_{m−1} = t_m + G_{m−1} = G_m`。 -/
theorem gainedLoss_le_jumpSum {p : ℤ} {t : ℕ → ℤ} (hp : 2 ≤ p) (ht : ∀ j, 0 ≤ t (j + 1))
    (hA : ∀ m, (p - 1) * jumpSum t m ≤ t (m + 1)) :
    ∀ k, gainedLoss p t k ≤ jumpSum t k
  | 0 => le_refl 0
  | (m + 1) => by
      have IH := gainedLoss_le_jumpSum hp ht hA m
      have hG : 0 ≤ jumpSum t m := jumpSum_nonneg ht m
      have ht1 : 0 ≤ t 1 := ht 0
      have hAm := hA m
      have hpG : 0 ≤ (p - 1) * jumpSum t m := mul_nonneg (by linarith) hG
      simp only [jumpSum, gainedLoss]
      apply max_le
      · linarith
      · have h1 : max 0 (t (m + 1) - (p - 1) * jumpSum t m - t 1)
            ≤ t (m + 1) - (p - 1) * jumpSum t m := by
          apply max_le <;> linarith
        have h2 : p * gainedLoss p t m ≤ p * jumpSum t m :=
          mul_le_mul_of_nonneg_left IH (by linarith)
        linarith

/-- ★★**核 2**: `p·t_j ≤ t_{j+1}` から `A_m ≥ 0`(= `(p−1)G_{m−1} ≤ t_m`)。

`(p−1)G_m = (p−1)G_{m−1} + (p−1)t_m ≤ t_m + (p−1)t_m = p·t_m ≤ t_{m+1}`。 -/
theorem topDefect_nonneg {p : ℤ} {t : ℕ → ℤ} (ht1 : 0 ≤ t 1)
    (hstep : ∀ j, p * t (j + 1) ≤ t (j + 2)) :
    ∀ m, (p - 1) * jumpSum t m ≤ t (m + 1)
  | 0 => by simpa [jumpSum] using ht1
  | (m + 1) => by
      have IH := topDefect_nonneg ht1 hstep m
      have h := hstep m
      simp only [jumpSum]
      nlinarith

/-- ★★**核 3**: 層の sharp な上界 `(p−1)t_j ≤ p^j e` から `(p−1)G_k ≤ e·Σ_{j≤k}p^j`。
★右辺は**予算そのもの**である。 -/
theorem jumpSum_le_geom {p : ℤ} {t : ℕ → ℤ} {e : ℤ}
    (hlayer : ∀ j, (p - 1) * t (j + 1) ≤ p ^ (j + 1) * e) :
    ∀ k, (p - 1) * jumpSum t k ≤ e * geomSum p k
  | 0 => by simp [jumpSum, geomSum]
  | (k + 1) => by
      have IH := jumpSum_le_geom hlayer k
      have h := hlayer k
      simp only [jumpSum, geomSum]
      have h1 : (p - 1) * (jumpSum t k + t (k + 1))
          = (p - 1) * jumpSum t k + (p - 1) * t (k + 1) := by ring
      have h2 : e * (geomSum p k + p ^ (k + 1)) = e * geomSum p k + p ^ (k + 1) * e := by ring
      rw [h1, h2]
      linarith

/-- ★★★★**主定理(抽象核)**。仮説は 3 つとも古典:

* `ht` —— 跳びは非負
* `hstep` —— `p·t_j ≤ t_{j+1}`(★二分律のどちらの枝からも出る。§2)
* `hlayer` —— 層ごとの sharp な跳びの上界 `(p−1)t_j ≤ p^j e`(`RamificationJumpBound`)

結論は台帳の形 `(p−1)²·Λ_k ≤ (p^{k+1} − p)·e` である。
★★**`p` にも `k` にも条件が付かない** —— これが `JumpDefectTradeoff` の
`min_cost_fits_of_le_three`(`p ≤ 3` かつ `k = 2` 限定)との差である。 -/
theorem gainedLoss_fits {p : ℤ} {t : ℕ → ℤ} {e : ℤ} (hp : 2 ≤ p) (ht : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ j, p * t (j + 1) ≤ t (j + 2))
    (hlayer : ∀ j, (p - 1) * t (j + 1) ≤ p ^ (j + 1) * e) (k : ℕ) :
    (p - 1) ^ 2 * gainedLoss p t k ≤ (p ^ (k + 1) - p) * e := by
  have hA := topDefect_nonneg (ht 0) hstep
  have h1 := gainedLoss_le_jumpSum hp ht hA k
  have h2 : (p - 1) ^ 2 * gainedLoss p t k ≤ (p - 1) ^ 2 * jumpSum t k :=
    mul_le_mul_of_nonneg_left h1 (by positivity)
  have h3 := jumpSum_le_geom hlayer k
  have h4 : (p - 1) * ((p - 1) * jumpSum t k) ≤ (p - 1) * (e * geomSum p k) :=
    mul_le_mul_of_nonneg_left h3 (by linarith)
  have h5 : (p - 1) * (e * geomSum p k) = e * ((p - 1) * geomSum p k) := by ring
  rw [h5, geomSum_mul p k] at h4
  nlinarith [h2, h4]

/-! ### ★§2 仮説 `p·t_j ≤ t_{j+1}` は**二分律のどちらの枝からも出る**(★`ℤ` の算術のみ) -/

/-- ★★**安定域の枝**: `t_{j+1} = t_j + p^{j+1}e` と層の上界 `(p−1)t_j ≤ p^{j+1}e` から
`p·t_j ≤ t_{j+1}`。★具体層では `u_j > e/(p−1)` のとき `(U^{(n)})^p = U^{(n+e)}` から
上付き跳びが `u_{j+1} = u_j + e` になる古典。 -/
theorem pmul_le_of_stable {p e t1 t2 : ℤ} {j : ℕ} (hstab : t2 = t1 + p ^ (j + 1) * e)
    (hlayer : (p - 1) * t1 ≤ p ^ (j + 1) * e) : p * t1 ≤ t2 := by
  subst hstab; linarith

/-- ★★**臨界以下の枝**: `u_{j+1} ≥ p·u_j` と `t_j ≤ p^j u_j` から `p·t_j ≤ t_{j+1}`。

`t_{j+1} = t_j + p^{j+1}(u_{j+1} − u_j) ≥ t_j + p^{j+1}(p−1)u_j ≥ t_j + p(p−1)t_j ≥ p·t_j`。
★`t_j ≤ p^j u_j` は `t_j = u_1 + Σ_{i<j} p^i(u_{i+1}−u_i)` と `p^i ≤ p^j` から出る古典。 -/
theorem pmul_le_of_below {p u1 u2 t1 t2 : ℤ} {j : ℕ} (hp : 2 ≤ p)
    (hgrow : p * u1 ≤ u2) (ht1 : t1 ≤ p ^ j * u1) (hu1 : 0 ≤ u1)
    (hstep : t2 = t1 + p ^ (j + 1) * (u2 - u1)) : p * t1 ≤ t2 := by
  subst hstep
  have hpj : (0:ℤ) < p ^ j := pow_pos (by linarith) j
  have key : (p - 1) * t1 ≤ p ^ (j + 1) * (u2 - u1) := by
    have h1 : (p - 1) * t1 ≤ (p - 1) * (p ^ j * u1) :=
      mul_le_mul_of_nonneg_left ht1 (by linarith)
    have h2 : p ^ (j + 1) * ((p - 1) * u1) ≤ p ^ (j + 1) * (u2 - u1) := by
      apply mul_le_mul_of_nonneg_left _ (le_of_lt (pow_pos (by linarith) (j + 1)))
      linarith
    have hX : 0 ≤ p ^ j * ((p - 1) * u1) :=
      mul_nonneg (le_of_lt hpj) (mul_nonneg (by linarith) hu1)
    have h3 : (p - 1) * (p ^ j * u1) ≤ p ^ (j + 1) * ((p - 1) * u1) := by
      have e1 : p ^ (j + 1) * ((p - 1) * u1) = p * (p ^ j * ((p - 1) * u1)) := by ring
      have e2 : (p - 1) * (p ^ j * u1) = p ^ j * ((p - 1) * u1) := by ring
      rw [e1, e2]
      nlinarith
    linarith
  linarith

/-! ## §3 深さ 2 —— ★**`p ≥ 5` の窓が消える** -/

/-- ★★★深さ 2 の鋭い損失 `max(T, p·S)`。★旧ファイルの `contrCost` / `projCost` の代わり。 -/
def sharpTwoCost (p S T : ℤ) : ℤ := max T (p * S)

theorem gainedLoss_one {p : ℤ} {t : ℕ → ℤ} (ht1 : 0 ≤ t 1) : gainedLoss p t 1 = t 1 := by
  show max (t 1 - (p - 1) * jumpSum t 0)
      (max 0 (t 1 - (p - 1) * jumpSum t 0 - t 1) + p * gainedLoss p t 0) = t 1
  simp only [jumpSum, gainedLoss, mul_zero, sub_zero, sub_self, add_zero]
  rw [max_self, max_eq_left ht1]

/-- ★★★**深さ 2 の閉じた形**: `Λ_2 = max(t_2, p·t_1)`。 -/
theorem gainedLoss_two {p : ℤ} {t : ℕ → ℤ} (hp : 2 ≤ p) (ht1 : 0 ≤ t 1) :
    gainedLoss p t 2 = sharpTwoCost p (t 1) (t 2) := by
  show max (t 2 - (p - 1) * jumpSum t 1)
      (max 0 (t 2 - (p - 1) * jumpSum t 1 - t 1) + p * gainedLoss p t 1)
    = sharpTwoCost p (t 1) (t 2)
  rw [gainedLoss_one ht1]
  have hj : jumpSum t 1 = t 1 := by simp [jumpSum]
  rw [hj]
  have key : max 0 (t 2 - (p - 1) * t 1 - t 1) + p * t 1 = max (p * t 1) (t 2) := by
    rcases le_total (t 2 - (p - 1) * t 1 - t 1) 0 with h | h
    · rw [max_eq_left h, max_eq_left (by linarith)]; ring
    · rw [max_eq_right h, max_eq_right (by linarith)]; ring
  rw [key, max_comm (p * t 1) (t 2)]
  have hnn : 0 ≤ (p - 1) * t 1 := mul_nonneg (by linarith) ht1
  exact max_eq_right (le_trans (by linarith) (le_max_left (t 2) (p * t 1)))

/-- ★★★★**本ファイルの主結果(深さ 2)**: 層の sharp な上界 2 つ**だけ**から
鋭い損失は対の予算に入る。

★★`JumpDefectTradeoff.min_cost_fits_of_le_three` と違い **`p ≤ 3` を仮定しない**。
★★`JumpDefectTradeoff.covered_iff` の二分律(`hdich`)も**要らない**。
⇒ ★★★`gap_pos_of_five_le` の「窓」は、鋭い上界では**空にならない**。 -/
theorem sharpTwoCost_fits {p e S T : ℤ} (hp : 2 ≤ p) (he : 0 ≤ e)
    (hS : (p - 1) * S ≤ p * e) (hT : (p - 1) * T ≤ p ^ 2 * e) :
    (p - 1) * sharpTwoCost p S T ≤ pairBudget p e := by
  have h : p * ((p - 1) * S) ≤ p * (p * e) := mul_le_mul_of_nonneg_left hS (by linarith)
  unfold sharpTwoCost pairBudget
  rcases max_cases T (p * S) with ⟨hm, _⟩ | ⟨hm, _⟩ <;> rw [hm]
  · nlinarith
  · nlinarith

/-- ★★★★**緩みの正体**: 旧上界の射影の枝は、鋭い値に**射影の欠損 `γ` を丸ごと足したもの**である。

  `projCost p e S T = sharpTwoCost p S T + (p²e − (p−1)T)`

★`γ = p²e − (p−1)T ≥ 0` は `JumpDefect.defect_nonneg_iff` により sharp な層の上界と同値。
★`p = 5, e = 3, S = 2, T = 17` では `γ = 7` で `24 = 17 + 7`。 -/
theorem projCost_eq_sharp_add_defect (p e S T : ℤ) :
    projCost p e S T = sharpTwoCost p S T + (p ^ 2 * e - (p - 1) * T) := by
  unfold projCost sharpTwoCost; ring

/-- ★★鋭い値は収縮の枝以下(`S ≤ T` のとき)。 -/
theorem sharpTwoCost_le_contrCost {p S T : ℤ} (hp : 2 ≤ p) (hS : 0 ≤ S) (hST : S ≤ T) :
    sharpTwoCost p S T ≤ contrCost p S T := by
  unfold sharpTwoCost contrCost
  rcases max_cases T (p * S) with ⟨hm, _⟩ | ⟨hm, _⟩ <;> rw [hm] <;> nlinarith

/-! ## §4 数値層 —— ★旧反例は消え、既存ファイルの実測値と**厳密に一致**する -/

namespace Numeric

/-- ★★★旧ファイルの反例 `(p, e_K, S, T) = (5, 3, 2, 17)`:
旧 `C = 25`、旧 `P = 24`、★**鋭い値は `17`**。 -/
theorem witness_five_sharp :
    sharpTwoCost 5 2 17 = 17 ∧ projCost 5 3 2 17 = 24 ∧ contrCost 5 2 17 = 25 := by
  refine ⟨by norm_num [sharpTwoCost], by norm_num [projCost], by norm_num [contrCost]⟩

/-- ★★★★**旧反例は対の予算に入る**(`4·17 = 68 ≤ 90`)。
⇒ `JumpDefectTradeoff.witness_five_exceeds_budget`(旧上界では `96 > 90`)は
★**旧上界についてのみ真**であった。 -/
theorem witness_five_fits : (5 - 1) * sharpTwoCost 5 2 17 ≤ pairBudget 5 3 := by
  norm_num [sharpTwoCost, pairBudget]

/-- ★★★★しかも**古典の Ax の定数 `p^{p/(p−1)²}` の下に入る**:
`(p−1)²·17 = 272 ≤ 375 = p·p²e`(旧上界は `384 > 375` で超えていた)。 -/
theorem witness_five_under_ax : (5 - 1) ^ 2 * sharpTwoCost 5 2 17 ≤ 5 * (5 ^ 2 * 3) := by
  norm_num [sharpTwoCost]

/-- ★★★`p ≥ 5` の**反例族全体**が鋭い上界では予算に入る:
`(e_K, S, T) = (p−2, 2, 2+p(p−2))` で `max(T, 2p) = T = 2 + p(p−2)` であり
`(p−1)(2 + p² − 2p) ≤ p(p+1)(p−2)` は `p ≥ 5` で成り立つ。
⇒ ★★`JumpDefectTradeoff.witness_gap_exceeds` の族は**消える**。 -/
theorem witness_gap_fits {p : ℤ} (hp : 5 ≤ p) :
    (p - 1) * sharpTwoCost p 2 (2 + p * (p - 2)) ≤ pairBudget p (p - 2) := by
  refine sharpTwoCost_fits (by linarith) (by linarith) (by nlinarith) (by nlinarith)

/-- ★★★`ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉) ⊃ ℚ₃(ζ₃)`: 鋭い値 `8` は
`EquivariantProjectionDescent` の**実測 `L = 8` と厳密に一致**する(旧上界は `10`)。 -/
theorem zeta27_sharp : sharpTwoCost 3 2 8 = 8 ∧ projCost 3 2 2 8 = 10 := by
  refine ⟨by norm_num [sharpTwoCost], by norm_num [projCost]⟩

/-- ★★★`ℚ₃(ζ₈₁) ⊃ ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉)`: 鋭い値 `26` は**実測 `L = 26` と厳密に一致**する
(旧上界は `28`)。 -/
theorem zeta81_sharp : sharpTwoCost 3 8 26 = 26 ∧ projCost 3 6 8 26 = 28 := by
  refine ⟨by norm_num [sharpTwoCost], by norm_num [projCost]⟩

/-- ★`ℚ₂(ζ₁₆)/ℚ₂`(★非巡回 `ℤ/2×ℤ/4`、本ファイルの適用外)の逆算値 `8` も
`max(T, pS) = max(7, 8) = 8` に**一致する**。★これは測定であって定理ではない。 -/
theorem zeta16_sharp : sharpTwoCost 2 4 7 = 8 := by norm_num [sharpTwoCost]

end Numeric

/-! ## §5 台帳 —— `ℤ` の不等式と `axDecay` の積を繋ぐ(深さ `k`) -/

/-- ★★**塔の予算の台帳**(`PairLedger.rpow_div_le_prod_Icc_axDecay_iff` に
`k' = 0`, `d = k`, `E = e_M = p^k e` を入れて `p^k` を約した形)。
★`JumpDefect.pairBudget_iff`(`k = 2`)の一般化である。 -/
theorem towerBudget_iff {p : ℕ} [Fact p.Prime] {J e k : ℕ} (he : 0 < e) :
    (p : ℝ) ^ ((J : ℝ) / ((p ^ k * e : ℕ) : ℝ)) ≤ ∏ i ∈ Finset.Icc 1 k, axDecay p i
      ↔ ((p : ℤ) - 1) ^ 2 * (J : ℤ) ≤ ((p : ℤ) ^ (k + 1) - p) * (e : ℤ) := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  have hE : 0 < p ^ k * e := by positivity
  have hIcc : (Finset.Icc 1 k : Finset ℕ) = Finset.Icc (0 + 1) (0 + k) := by norm_num
  rw [hIcc, PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := p) (k' := 0) (d := k)
    (J := J) (E := p ^ k * e) hE]
  have h1 : 1 ≤ p := by omega
  have h2 : 1 ≤ p ^ k := Nat.one_le_pow _ _ (by omega)
  zify [h1, h2]
  have hpos : (0:ℤ) < (p:ℤ) ^ k := by
    have : (0:ℤ) < (p:ℤ) := by exact_mod_cast (by omega : 0 < p)
    positivity
  have e1 : ((p:ℤ) - 1) ^ 2 * (p:ℤ) ^ (0 + k) * (J:ℤ)
      = (((p:ℤ) - 1) ^ 2 * (J:ℤ)) * (p:ℤ) ^ k := by ring
  have e2 : ((p:ℤ) ^ k - 1) * (p:ℤ) * ((p:ℤ) ^ k * (e:ℤ))
      = (((p:ℤ) ^ (k + 1) - (p:ℤ)) * (e:ℤ)) * (p:ℤ) ^ k := by ring
  rw [e1, e2, mul_le_mul_iff_left₀ hpos]

/-- ★★★★**`AxLemmaGraded` の形**: 損失 `J` が `gainedLoss` で押さえられるなら
`p^{J/e_M} ≤ ∏_{i∈[1,k]} axDecay p i`。

★★仮説は「跳びが非負」「`p·t_j ≤ t_{j+1}`(二分律から)」「層の sharp な上界」の 3 つだけで、
★**`p` にも `k` にも条件が付かない**。 -/
theorem rpow_le_prod_axDecay {p : ℕ} [Fact p.Prime] {e k J : ℕ} {t : ℕ → ℤ} (he : 0 < e)
    (ht : ∀ j, 0 ≤ t (j + 1)) (hstep : ∀ j, (p : ℤ) * t (j + 1) ≤ t (j + 2))
    (hlayer : ∀ j, ((p : ℤ) - 1) * t (j + 1) ≤ (p : ℤ) ^ (j + 1) * (e : ℤ))
    (hJ : (J : ℤ) ≤ gainedLoss (p : ℤ) t k) :
    (p : ℝ) ^ ((J : ℝ) / ((p ^ k * e : ℕ) : ℝ)) ≤ ∏ i ∈ Finset.Icc 1 k, axDecay p i := by
  rw [towerBudget_iff he]
  have hp2 : 2 ≤ (p : ℤ) := by exact_mod_cast (Fact.out : p.Prime).two_le
  have hfit := gainedLoss_fits (p := (p : ℤ)) (t := t) (e := (e : ℤ)) hp2 ht hstep hlayer k
  have h2 : ((p : ℤ) - 1) ^ 2 * (J : ℤ) ≤ ((p : ℤ) - 1) ^ 2 * gainedLoss (p : ℤ) t k :=
    mul_le_mul_of_nonneg_left hJ (by positivity)
  linarith

/-! ## §6 `.src`(原典の対応箇所) -/

def gainedLoss.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def gainedLoss_fits.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def sharpTwoCost_fits.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def projCost_eq_sharp_add_defect.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def towerBudget_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def rpow_le_prod_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end GainedDescent

end ABC3.Found.PGC

section Audit
open ABC3.Found.PGC
#print axioms GainedDescent.geomSum_mul
#print axioms GainedDescent.jumpSum_nonneg
#print axioms GainedDescent.gainedLoss_le_jumpSum
#print axioms GainedDescent.topDefect_nonneg
#print axioms GainedDescent.jumpSum_le_geom
#print axioms GainedDescent.gainedLoss_fits
#print axioms GainedDescent.pmul_le_of_stable
#print axioms GainedDescent.pmul_le_of_below
#print axioms GainedDescent.gainedLoss_one
#print axioms GainedDescent.gainedLoss_two
#print axioms GainedDescent.sharpTwoCost_fits
#print axioms GainedDescent.projCost_eq_sharp_add_defect
#print axioms GainedDescent.sharpTwoCost_le_contrCost
#print axioms GainedDescent.Numeric.witness_five_sharp
#print axioms GainedDescent.Numeric.witness_five_fits
#print axioms GainedDescent.Numeric.witness_five_under_ax
#print axioms GainedDescent.Numeric.witness_gap_fits
#print axioms GainedDescent.Numeric.zeta27_sharp
#print axioms GainedDescent.Numeric.zeta81_sharp
#print axioms GainedDescent.towerBudget_iff
#print axioms GainedDescent.rpow_le_prod_axDecay
end Audit
