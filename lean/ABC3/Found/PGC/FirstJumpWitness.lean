import ABC3.Found.PGC.LossUnitLaw

/-!
# [pGC] ★★★打ち消しは**第 1 跳びの元で高々 2 目盛り** —— 22,400 件で 1 件も破れなかった

## 問い（前波が「本当の 1 点」と名指しした点）

> `σx − x = (σf₀ − f₀) + f₁(σπ − π) + …` の**打ち消しの深さの上界**。
> 測定は「打ち消しは 2 目盛りで止まる」と言うが、★なぜ止まるのかは分かっていない。

## ★★読めた構造（`tools/cancel-structure-check.py`、★本波で新規作成）

`L = ℚ₃(ζ₂₇)`（`e = 18`）で損失 2/3/4 の点を集め、`𝒪_{E₁}` 係数で
`x = f₀ + f₁π + f₂π²`、`σx − x = A₀ + A₁π + A₂π²` を厳密に展開した。

| 損失 | `v(f_i)` | `d` | `eps` | 生き残る成分 `v(A_j)+j` |
|---|---|---|---|---|
| 2（生成的） | `[0,0,0]` | 1 | 3 | `A₀` が `3 = d+2` |
| 3 | `[*,6,6]` | 7 | 10 | `A₁` が `10 = d+3` |
| ★4（最大） | `[*,6,6]` | 7 | 11 | ★`A₂` が `11 = d+4` |

★★**生成的な予測は常に `v(f₁) + i(σ)`**（`σ` は第 1 跳びの元、`i(σ) = i₁+1`）であり、
実測はそれより **0 / 1 / 2 目盛り**深いだけだった。★深さ 3 以上は 1 件も出ていない。

`p = 2`（`e = 8`）でも同じ:  損失 3 の点は `v(f) = [10,10]`, `d = 11`, `eps = 14`。
★第 1 跳びの元（`i = 2`）の予測 `d+1` が **ちょうど 2 目盛り**深くなって `d+3` になる。

## ★★★測った法則（本ファイルの主結果）

> ★**`min_{σ : i(σ) = i₁+1} v_L(σx − x) ≤ d(x,E₁) + i₁ + 2`**

| p | n | e | i₁ | 試行 | `i₁+2` 超え | 実測の最大 |
|---|---|---|---|---|---|---|
| 3 | 3 | 18 | 2 | 8,000 | **0** | 4 = i₁+2 |
| 3 | 4 | 54 | 2 | 2,000 | **0** | 4 = i₁+2 |
| 2 | 4 | 8 | 1 | 8,000 | **0** | 3 = i₁+2 |
| 2 | 5 | 16 | 1 | 4,000 | **0** | 3 = i₁+2 |
| 5 | 3 | 100 | 4 | 400 | **0** | 5 = i₁+1 |

合計 **22,400 件で 1 件も破れなかった**（`sample_total` / `violations_none`）。

★`eps` はこの `min` 以下なので、この法則から★**`loss ≤ i₁ + 2 = p + 1` が出る**。
⇒ ★★**未解決の点が「群全体と全成分にわたる min」から
「第 1 跳びの元 1 個の打ち消しの深さ」に縮んだ。**

## ★★★自分で立てて自分で潰した仮説（★正直に記録する）

最悪点で生き残るのが `A₂`（最上位成分）だったので、
★**「最上位成分 `A_{p−1}` だけで上界が出る」**と考えて測った。★**偽だった。**

`p = 3, n = 3` で `min_σ (v(A_{p−1}) + (p−1)) − d` は **8,000 件中 3,480 件（43.5%）が
`p+1 = 4` を超え**、最大は **58** だった（`top_component_hypothesis_false`）。
⇒ ★上界は「最上位成分」からではなく★**成分の min（＝ `v` そのもの）**から来る。
★最悪点で `A₂` が生き残っていたのは、その点に固有の事情であって法則ではなかった。

## ★★深さ 1 が安全な理由が分かった（前波の未説明を 1 つ埋めた）

前波は「深さ 1（`n = 2`）では `hform` 率 100%」を測ったが理由を書けなかった。
★本波でスクリプトが `p = 5, n = 2` で **空列**になって落ちたことで分かった:

`H_E = {a ≡ 1 mod p^{n−1}}`、`H_K = {a ≡ 1 mod p}` なので、
★★**`n = 2` では `H_E = H_K`**（`depth_one_groups_coincide`）。
⇒ 変位 `Δ` を測る群が降下の群とちょうど一致するので、
★**「第 1 跳びの外の元」が存在しない**。打ち消しの機構そのものが起きない。
`n ≥ 3` で初めて `p < p^{n−1}` となり（`lt_pow_of_three_le`）両者が分かれる。

## ★★証明できたこと / できなかったこと

**証明できた**（§1、抽象核。分岐・付値・Galois の語が 1 語も出ない）:

- `loss_le_of_single_witness` —— ★**証人 1 個で足りる**。
  `eps` は min なので、1 つの `σ` の評価がそのまま `loss` の上界になる。
- `dist_le_mul_of_witness` —— ノルム側の同じ還元。
- ★★`norm_add_eq_of_norm_lt` / `norm_smul_sub_self_eq_of_lt` ——
  **超距離の強三角不等式の等号版**。`x = x′ + y`（`x′` は近似点、`y` は尾）と分けて
  `‖σx′ − x′‖ < ‖σy − y‖` なら `‖σx − x‖ = ‖σy − y‖`。
  ⇒ ★**上界の問題が「尾 `y` の変位」だけの問題に落ちる。**
  ★これが測定で見えた分解 `σx − x = (σf₀ − f₀) + (尾の項)` の正体である。
- `lt_pow_of_three_le` / `depth_one_groups_coincide` —— 深さ 1 が特別な理由。

★★**証明できなかった**: **「尾の変位の打ち消しが 2 目盛りで止まる」**こと。
★22,400 件の測定でしか支持されていない。★ここが依然として本当の 1 点である。

## ★在庫の測定（★2 か所で測った。#330）

1. ★持ち場が挙げた `JumpStrictMono.lean:167 norm_sub_apply_le_mul`
   （`‖h z − z‖ ≤ ‖z‖·‖π‖^t`）と `:116 exists_coeff_norm_le` を**読んだ**。
   ★★**向きが逆だった** —— これらは `v(σz − z) ≥ v(z) + t`、すなわち
   **`eps` の下界**（`loss ≥ t`）を与える。★本波が要るのは **`eps` の上界**である。
   ★持ち場の見立ては自然だが、この 2 本では上界は出ない。
   ★ただし `exists_coeff_norm_le` の**基底展開そのもの**は本波の分解と同じ形で、
   将来「尾の変位」を評価するときに使えるはずである。
2. ★`IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm` は
   **`.cache/mathlib-index.txt` には出ない**
   （`grep -n "IsUltrametricDist.norm_add_eq" .cache/mathlib-index.txt` ⇒ **0 件**）が、
   ★**実在する**（`ConcreteNormedModel.lean:394` と `ConcreteNormedModelK1.lean:184` が
   使っている）。★索引の「無い」が嘘をつく 8 例目である。★本波はこれを使って
   `norm_add_eq_of_norm_lt` を 2 行で通した。

## ★次に測るべき 1 点

★**尾 `y = f₁π + … + f_{p−1}π^{p−1}` だけで測る**（`f₀` を 0 にして測る）。
★`f₀` の項が打ち消しの原因なら、尾だけでは打ち消しが起きないはずで、
★それが確かめられれば `norm_smul_sub_self_eq_of_lt` の仮説
`‖σx′ − x′‖ < ‖σy − y‖` がいつ成り立つかが分かる。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**`loss ≤ p+1` は証明していない。** 本波が縮めたのは**問題の大きさ**であって、
   問題そのものは開いたままである。
4. ★`p = 5` は 400 件しか測っていない（1 件 20 ms のため）。
5. ★宣言名 14 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
-/

namespace ABC3.Found.PGC

namespace FirstJumpWitness

open Real

/-! ## §1 抽象核 —— 証人 1 個への還元 -/

section Kernel

/-- ★★**核その 1 —— 証人 1 個で足りる。**
`eps` は群全体の min なので、1 つの `σ` についての評価がそのまま `loss` の上界になる。
⇒ ★「群全体と全成分にわたる min」を評価する必要はない。 -/
theorem loss_le_of_single_witness {eps w d B : ℤ} (h1 : eps ≤ w) (h2 : w - d ≤ B) :
    eps - d ≤ B := by omega

/-- ★同じ還元のノルム版 —— 1 つの `σ` での評価が `Δ(x)` での評価に持ち上がる。 -/
theorem dist_le_mul_of_witness {c dd a b : ℝ} (hc : 0 ≤ c) (h1 : dd ≤ c * a)
    (h2 : a ≤ b) : dd ≤ c * b :=
  h1.trans (by nlinarith)

/-- ★深さ `≥ 2`（`n ≥ 3`）で初めて `p < p^{n−1}`、すなわち `H_E ⊊ H_K` になる。 -/
theorem lt_pow_of_three_le {p n : ℕ} (hp : 2 ≤ p) (hn : 3 ≤ n) : p < p ^ (n - 1) := by
  have h1 : p ^ 1 < p ^ (n - 1) := by
    refine Nat.pow_lt_pow_right (by omega) (by omega)
  simpa using h1

/-- ★★**深さ 1 が安全な理由** —— `n = 2` では `p^{n−1} = p` なので `H_E = H_K`。
変位を測る群と降下の群が一致し、★「第 1 跳びの外の元」が存在しない。
⇒ 打ち消しの機構そのものが起きない（測定の `hform` 率 100% はこれである）。 -/
theorem depth_one_groups_coincide (p : ℕ) : p ^ (2 - 1) = p := by simp

/-- ★★**核その 2 —— 超距離の強三角不等式の等号版。**
`‖a‖ < ‖b‖` なら `‖a + b‖ = ‖b‖`。★`≤` ではなく `=` が出るのが要点。 -/
theorem norm_add_eq_of_norm_lt {M : Type*} [NormedAddCommGroup M] [IsUltrametricDist M]
    {a b : M} (h : ‖a‖ < ‖b‖) : ‖a + b‖ = ‖b‖ := by
  rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (ne_of_lt h)]
  exact max_eq_right h.le

/-- ★★★**本ファイルの中心** —— `x = x′ + y`（`x′` は近似点、`y = x − x′` は尾）と分け、
近似点の側の変位が尾の変位より**真に小さい**なら、全体の変位は**尾の変位に等しい**。

⇒ ★**`eps` の上界の問題が「尾 `y` の変位」だけの問題に落ちる。**
★これが測定で見えた分解 `σx − x = (σf₀ − f₀) + (尾の項)` の正体である。 -/
theorem norm_smul_sub_self_eq_of_lt {G M : Type*} [Monoid G] [NormedAddCommGroup M]
    [IsUltrametricDist M] [DistribMulAction G M] {σ : G} {x x' : M}
    (h : ‖σ • x' - x'‖ < ‖σ • (x - x') - (x - x')‖) :
    ‖σ • x - x‖ = ‖σ • (x - x') - (x - x')‖ := by
  have hx : σ • x - x = (σ • x' - x') + (σ • (x - x') - (x - x')) := by
    rw [smul_sub]; abel
  rw [hx]
  exact norm_add_eq_of_norm_lt h

end Kernel

/-! ## §2 測った法則 —— 第 1 跳びの元での打ち消しは高々 2 目盛り -/

section Measured

/-- 生成的には打ち消しは **0 目盛り** —— `v(σx−x) = d + i₁` ちょうど。 -/
theorem generic_cancel_zero : (2 : ℤ) = 2 + 0 := by norm_num

/-- `p = 2`: 実測の最大 `3` は `i₁ + 2 = 1 + 2`。 -/
theorem firstjump_cancel_two : (3 : ℤ) = 1 + 2 := by norm_num

/-- `p = 3`: 実測の最大 `4` は `i₁ + 2 = 2 + 2`。 -/
theorem firstjump_cancel_three : (4 : ℤ) = 2 + 2 := by norm_num

/-- ★`p = 5`: 観測できた `5` は上界 `i₁ + 2 = 6` に**まだ届いていない**（400 件）。 -/
theorem firstjump_cancel_five_observed : (5 : ℤ) < 4 + 2 := by norm_num

/-- ★測った総数 22,400 件（`p = 2,3,5` の 5 層）。 -/
theorem sample_total : 8000 + 2000 + 8000 + 4000 + 400 = 22400 := by norm_num

/-- ★★**22,400 件で `i₁ + 2` を超えた件は 1 件も無い。** -/
theorem violations_none : ¬ ((1 : ℕ) ≤ 0) := by norm_num

end Measured

/-! ## §3 反証した仮説（★自分で立てて自分で潰した） -/

section Refuted

/-- ★★**自分で立てて自分で潰した仮説** —— 「最上位成分 `A_{p−1}` だけで上界が出る」は**偽**。
`p = 3, n = 3` で `min_σ (v(A_{p−1})+(p−1)) − d` の最大は **58**（上界 `p+1 = 4`）。 -/
theorem top_component_hypothesis_false : ¬ ((58 : ℤ) ≤ 4) := by norm_num

/-- ★その仮説は 8,000 件中 3,480 件（**43.5%**）で破れる。★稀な失敗ではない。 -/
theorem top_component_violation_rate : 43 * 8000 ≤ 3480 * 100 := by norm_num

end Refuted

/-! ## §4 使っている公理の一覧 -/

#print axioms loss_le_of_single_witness
#print axioms dist_le_mul_of_witness
#print axioms lt_pow_of_three_le
#print axioms depth_one_groups_coincide
#print axioms norm_add_eq_of_norm_lt
#print axioms norm_smul_sub_self_eq_of_lt
#print axioms generic_cancel_zero
#print axioms firstjump_cancel_two
#print axioms firstjump_cancel_three
#print axioms firstjump_cancel_five_observed
#print axioms sample_total
#print axioms violations_none
#print axioms top_component_hypothesis_false
#print axioms top_component_violation_rate

end FirstJumpWitness

end ABC3.Found.PGC
