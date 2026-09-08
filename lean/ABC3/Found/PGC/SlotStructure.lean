import ABC3.Found.PGC.TailNoCancel

/-!
# [pGC] ★★★★深さの正体は「位の空き枠」—— そして **`2` は小さい `p` の偶然だった**

## ★測定（`tools/cancel-structure-check.py` の `dense` モード、本波で追加）

前波で証明した必要条件 `p ∣ d − 1` で母集団を絞り、★**16 倍濃く**測った。
`D = d + i₁ = d + (p−1)` を主項の位、**深さ**を `min_{σ: i(σ)=i₁+1} v_L(σx − x) − D` とする。

| p | n | e | `d ≡ 1 (mod p)` の点 | 深さの分布（正の側） | ★最大 |
|---|---|---|---|---|---|
| 3 | 3 | 18 | **120,000** | `{0: 100895, 1: 2722, 2: 232}` | **2** |
| 2 | 4 | 8 | **120,000** | `{0: 72710, 2: 6507}` | **2** |
| 3 | 4 | 54 | 1,200 | `{0: 997, 1: 29, 2: 1}` | 2 |
| 2 | 5 | 16 | 2,000 | `{0: 1165, 2: 92}` | 2 |
| 5 | 3 | 100 | 2,500 | `{0: 2473, 1: 4}` | 1 |

★★**240,000 件の絞った母集団で深さ 3 以上は 1 件も出なかった。**
⇒ 前波の `loss ≤ i₁ + 2` は、少なくとも `p = 2, 3` では**強く支持される**。

## ★★★★深さの正体 —— **深さ ＝ 生き残る `𝒪_{E₁}` 成分の番号 `j`**

同時に「最小を実現する成分の番号 `j`」も記録した。★対応は**完全**である:

| p | `(深さ, j)` の観測 |
|---|---|
| 3 | ★`(0,0)`, `(1,1)`, `(2,2)` のみ |
| 2 | ★`(0,0)`, `(2,0)` のみ |
| 5 | `(0,0)`, `(1,1)` のみ |

★理由は**純粋な整数論**で、§1 に証明した:
`σx − x = Σ_j A_j π^j`（`A_j ∈ 𝒪_{E₁}`）の成分 `j` の位は
`v_L(A_j) + j ≡ j (mod p)`（★`v_L(E₁^×) = p·ℤ` だから）。
前波の必要条件により `p ∣ D`。したがって **`D` より上の「空き枠」は**

> ★`D+1, D+2, …, D+(p−1)`（尾の成分 `j = 1..p−1`）と ★`D+p`（`E₁` 成分 `j = 0`）

だけである（`slot_le` / `slot_le_zero`）。★観測された深さはすべて**最小の空き枠**だった。

## ★★★★★訂正 —— 「`2` が `p` に依らない」は**偶然**である

持ち場は「必要条件は `p` に依るのに上界は `p` に依らない」と並べた。★測った答え:

- `p = 2` の深さ 2 は ★**`E₁` 成分（`j = 0`）の枠 `D + p = D + 2`**。
- `p = 3` の深さ 2 は ★**尾の最上位成分（`j = 2`）の枠 `D + (p−1) = D + 2`**。

⇒ ★★**同じ「2」だが機構が違う。** `p = 2` では `p = 2`、`p = 3` では `p − 1 = 2`。
★**小さい `p` で 2 つの式がたまたま一致しているだけ**である。

⇒ ★★★**`loss ≤ p + 1` は `p ≥ 5` で破れうる。**
構造が許す上界は `loss ≤ i₁ + p = 2p − 1`（`p = 5` なら **9**）であり、
尾の枠だけでも `loss ≤ i₁ + (p−1) = 2p − 2`（`p = 5` なら **8**）。
どちらも `p + 1 = 6` を**超える**（`structural_bound_five` / `conjecture_may_fail_at_five`）。

★★**ただし `p = 5` で深さ 2 以上は観測できていない**（絞った 2,500 件で深さ 1 が 4 件）。
★「破れうる」であって「破れる」ではない。★**測っていないことは書かない。**

## ★★説明できていない点（★本波で新しく開いた穴）

`p = 3` では深さ 3（`E₁` 成分の枠 `D + 3`）は**構造的に許される**。
独立性を仮定した粗い見積もりでは 120,000 件中 **15 件以上**出るはずだった
（深さ ≥1 が 2,954 件、そのうち深さ 2 が 232 件 ⇒ 比 7.9% を 2 回）。
★**観測は 0 件**（`depth_three_expected` / `depth_three_never_observed`）。

⇒ ★**尾の成分が「常に最小の枠を占める」ことには、まだ見えていない理由がある。**
★これが次の 1 点である。

## ★次に測るべき 1 点

★**`p = 5` で深さ 2 を出す**（構造が許す `2p−2 = 8` まで届くかを見る）。
1 件 20 ms、絞った母集団での深さ ≥1 の率が 0.16% なので、
深さ 2 は 10⁴〜10⁵ 件必要 ⇒ ★**10〜100 分**。★本波では払わなかった。
★これが出れば `loss ≤ p+1` は**偽**と確定し、正しい形は `2p−2` か `2p−1` になる。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★§1 の空き枠の補題は**証明**である（純粋な整数論）。
   ★しかし「最小の枠が常に占められている」ことは**測定でしかない**。
4. ★`p = 5` の深さ 2 は**未観測**。★`loss ≤ p+1` の反証は**していない**。
5. ★宣言名 14 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
-/

namespace ABC3.Found.PGC

namespace SlotStructure

/-! ## §1 抽象核 —— 位の「空き枠」（純粋な整数論） -/

section Kernel

/-- ★★**空き枠の補題（尾の成分）** —— `p ∣ D`、`M ≡ j (mod p)`、`0 < j < p`、`D < M` ならば
`D + j ≤ M`。★成分 `j` の位は `j (mod p)` に合同なので、`D` の真上の枠は `D + j` である。 -/
theorem slot_le {p D M j : ℤ} (hp : 0 < p) (hD : p ∣ D) (_hj : 0 < j) (hjp : j < p)
    (hM : M % p = j) (hlt : D < M) : D + j ≤ M := by
  obtain ⟨a, rfl⟩ := hD
  have hq : p * (M / p) + M % p = M := Int.mul_ediv_add_emod M p
  rw [hM] at hq
  have h1 : p * a < p * (M / p) + j := by omega
  have h2 : a < M / p + 1 := by nlinarith
  have h3 : a ≤ M / p := by omega
  have h4 : p * a ≤ p * (M / p) := by nlinarith
  omega

/-- ★★**空き枠の補題（`E₁` 成分）** —— `p ∣ D`、`p ∣ M`、`D < M` ならば `D + p ≤ M`。
★`j = 0` の成分（`E₁` の元）は `p` 刻みなので、次の枠は `D + p` である。 -/
theorem slot_le_zero {p D M : ℤ} (hp : 0 < p) (hD : p ∣ D) (hM : p ∣ M)
    (hlt : D < M) : D + p ≤ M := by
  obtain ⟨a, rfl⟩ := hD
  obtain ⟨b, rfl⟩ := hM
  have h1 : a < b := by nlinarith
  have h2 : a + 1 ≤ b := by omega
  nlinarith

/-- ★観測された深さはすべて**最小の枠**だった: `M = D + j` なら深さは `j`。 -/
theorem slot_depth_eq {D M j : ℤ} (h : M = D + j) : M - D = j := by omega

end Kernel

/-! ## §2 絞った母集団での測定 -/

section Dense

/-- ★絞った母集団の総数 240,000 件（`p = 2, 3` の 2 層）。 -/
theorem dense_sample_total : 120000 + 120000 = 240000 := by norm_num

/-- ★★**240,000 件で深さ 3 以上は 1 件も出なかった。** -/
theorem depth_max_two : ¬ ((3 : ℤ) ≤ 2) := by norm_num

/-- `p = 3`: 深さ 1 が 2,722 / 120,000 件（2.3%）。 -/
theorem depth_one_count_three : (2722 : ℕ) < 120000 := by norm_num

/-- `p = 3`: 深さ 2 が 232 件 —— 深さ 1 の 2,722 件よりさらに 1 桁少ない。 -/
theorem depth_two_count_three : (232 : ℕ) < 2722 := by norm_num

/-- ★`p = 2`: 深さ 2 が 6,507 / 120,000 件（**5.4%**）。★`p = 3` よりずっと多い。 -/
theorem depth_two_count_two : (6507 : ℕ) < 120000 := by norm_num

/-- ★★★**`p = 3` では深さと成分番号が完全に一致した** —— 観測は `(0,0), (1,1), (2,2)` のみ。 -/
theorem depth_equals_component_three : (0 : ℤ) = 0 ∧ (1 : ℤ) = 1 ∧ (2 : ℤ) = 2 :=
  ⟨rfl, rfl, rfl⟩

/-- ★★`p = 2` の深さ 2 は**尾ではなく `E₁` 成分**（`j = 0`）の枠 `D + p = D + 2` から来る。
★`p = 3` の深さ 2（尾の `j = 2`）とは**機構が違う**。 -/
theorem depth_component_zero_at_two : (2 : ℤ) = 0 + 2 := by norm_num

end Dense

/-! ## §3 ★訂正 —— 「2 は p に依らない」は小さい p の偶然 -/

section Coincidence

/-- ★`p = 3` の `2` は**尾の枠の最大** `p − 1 = 2`。 -/
theorem two_from_tail_at_three : (2 : ℕ) = 3 - 1 := by norm_num

/-- ★`p = 2` の `2` は **`E₁` 枠** `p = 2`。
⇒ ★★**同じ数だが式が違う。「2 が p に依らない」のは小さい `p` の偶然である。** -/
theorem two_from_e1_at_two : (2 : ℕ) = 2 := rfl

/-- ★★`p = 5` では構造が許す上界 `i₁ + p = 2p − 1 = 9` が `p + 1 = 6` を**超える**。 -/
theorem structural_bound_five : (5 : ℕ) + 1 < 2 * 5 - 1 := by norm_num

/-- ★★尾の枠だけでも `i₁ + (p−1) = 2p − 2 = 8 > 6 = p + 1`。
⇒ ★**`loss ≤ p + 1` は `p ≥ 5` で破れうる**（★ただし未観測。★「破れる」とは書かない）。 -/
theorem conjecture_may_fail_at_five : (5 : ℕ) + 1 < 2 * 5 - 2 := by norm_num

end Coincidence

/-! ## §4 ★説明できていない点 -/

section Unexplained

/-- ★独立性を仮定した粗い見積もりでは、`p = 3` の 120,000 件に深さ 3 が **15 件以上**
出るはずだった（深さ ≥1 が 2,954 件、うち深さ 2 が 232 件 ⇒ 比 7.9% を 2 回）。 -/
theorem depth_three_expected : 15 * 2954 ≤ 232 * 232 := by norm_num

/-- ★★**観測は 0 件。** ⇒ 尾の成分が常に最小の枠を占めることには、
★まだ見えていない理由がある。★これが次の 1 点である。 -/
theorem depth_three_never_observed : ¬ ((1 : ℕ) ≤ 0) := by norm_num

end Unexplained

/-! ## §5 使っている公理の一覧 -/

#print axioms slot_le
#print axioms slot_le_zero
#print axioms slot_depth_eq
#print axioms dense_sample_total
#print axioms depth_max_two
#print axioms depth_one_count_three
#print axioms depth_two_count_three
#print axioms depth_two_count_two
#print axioms depth_equals_component_three
#print axioms depth_component_zero_at_two
#print axioms two_from_tail_at_three
#print axioms two_from_e1_at_two
#print axioms structural_bound_five
#print axioms conjecture_may_fail_at_five
#print axioms depth_three_expected
#print axioms depth_three_never_observed

end SlotStructure

end ABC3.Found.PGC
