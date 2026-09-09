import ABC3.Found.PGC.LossTwoPSubTwo

/-!
# [pGC] ★★★★★残る仮説が「**指数が違う**」だけになった —— 一意性は要らなかった

## ★★★★★開いて分かったこと（★一意性は**回り道**だった）

前波で私は「`hz` の供給には**展開の一意性**が要る」と書いた。★本波で開いて測った結果、
★★**一意性は要らない。**

`z = B·π^j + R` と分けたとき、★**`‖B·π^j‖ ≠ ‖R‖` でありさえすれば**
超距離の等号版（`norm_add_eq_max_of_norm_ne_norm`）で

`‖z‖ = max(‖B·π^j‖, ‖R‖) ≥ ‖B·π^j‖`

が出る（`norm_le_of_norm_ne`、★**2 行**）。★これが `hz` そのものである。
★★**「どの分解か」は問題ではなく、「2 つのノルムが違う」ことだけが要る。**

⇒ ★★**残る仮説が `hz`（存在＋一意性）から
`hsep : ‖B·π^j‖ ≠ ‖R‖`（分離条件）に落ちた。**

## ★★★さらに 1 段 —— 分離条件は「**指数が違う**」に落ちる

`loss_le_of_exp_ne` : `‖B‖ = ‖π‖^{d+(p−1)}`、`‖R‖ = ‖π‖^v` のとき

> ★**`d + (p−1) + j ≠ v` ならば `‖π‖^{d+(2p−2)} ≤ ‖B·π^j + R‖`**

★★**これが今日の到達点である。**残る仮説は★**ただ 1 つの自然数の不等式**
`d + (p−1) + j ≠ v` で、★分岐の言葉では
「成分 `j` の位と、残りの部分の位が**一致しない**」——
`v(B·π^j) ≡ j (mod p)` と `v(R) ≢ j (mod p)` から出る古典である
（★`SlotStructure.slot_le` で本日**証明済み**の枠の議論と同じ）。

## ★★5 段 + 供給の現在地

| 段 | Lean | 状態 |
|---|---|---|
| (0)〜(6) | `TailNoCancel` / `RhoFactorization` / `BinomialFirstOrder` / `TopIndexSurvives` / `ComponentWitness` | ★**定理** |
| (7) 全体 | `LossTwoPSubTwo.loss_le_two_p_sub_two_of_component` | ★**定理**（仮説 `hz`） |
| ★(8) `hz` の供給 | ★**本ファイル `norm_le_of_norm_ne`** | ★**定理**（2 行） |
| ★(9) 分離の供給 | ★**本ファイル `norm_ne_of_exp_ne`** | ★**定理** |
| ★★残り | `d+(p−1)+j ≠ v`（★位が `p` を法で違う） | ★**未着手** |

## ★★次の波が最短で入るための `file:line`（★持ち場の要請）

★**残る 1 点**: `v(R) ≢ j (mod p)` を示すこと。材料:

- `SlotStructure.lean` の `slot_le` / `slot_le_zero`（★**本日証明済み**。
  「成分 `j` の位は `j (mod p)` に合同」の議論）。
- `TotallyRamifiedLayer.lean` の `not_dvd_sub_of_lt`（★`π^i` の位が `n` を法で相異なる）。
- `JumpStrictMono.lean:116 exists_coeff_norm_le`（★もはや**要らない** ——
  ★本波で一意性の道が不要と分かったため）。

★5 段の証明の**本文**は `TopIndexSurvives.lean` の docstring にある（★4 波連続で書く）。

## ★母集団（★規律どおり）

★**`p ≥ 3`**。`p = 2` は上界が `2p−1 = 3`、理由も証明済み
（`TwoIsDegenerate.mul_sub_one_eq_self_iff`）。
実測は `p ≥ 3` で破れ **0**（4,187 件）、`p = 2` でのみ 5.7%。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**自己訂正**: 前波の「一意性が要る」は**回り道だった**。
   ★超距離の等号版で 2 行で済んだ。★**「要る」と書くときも、まず 2 行を試すこと。**
4. ★★残る仮説（位が `p` を法で違うこと）は**証明していない**。
   ★「`loss ≤ 2p−2` が（無条件に）定理になった」とは書かない。
5. ★宣言名 4 件は **実ファイル**で衝突 0 を確かめた（#348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace SeparatedComponent

/-! ## §1 抽象核 —— ノルムが違えば片方は全体以下 -/

section Kernel

/-- ★★★**`hz` の供給（2 行）** —— `‖a‖ ≠ ‖r‖` なら `‖a‖ ≤ ‖a + r‖`。

★超距離の**等号版**（`norm_add_eq_max_of_norm_ne_norm`）だけで出る。
★★前波で「展開の一意性が要る」と書いたのは**回り道**だった。 -/
theorem norm_le_of_norm_ne {M : Type*} [NormedAddCommGroup M] [IsUltrametricDist M]
    {a r : M} (h : ‖a‖ ≠ ‖r‖) : ‖a‖ ≤ ‖a + r‖ := by
  rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm h]
  exact le_max_left _ _

end Kernel

/-! ## §2 `loss ≤ 2p − 2` —— 仮説が**分離条件**になった -/

section Separated

/-- ★★★★`loss ≤ 2p−2`（★仮説が**分離条件**になった版）。

`z = B·π^j + R` で `‖B·π^j‖ ≠ ‖R‖` なら結論が出る。
★「どの分解か」は問題ではなく、★**2 つのノルムが違う**ことだけが要る。 -/
theorem loss_le_of_separated {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π B R : M} {p d j : ℕ}
    (hπ0 : 0 ≤ ‖π‖) (hπ1 : ‖π‖ ≤ 1) (hp : 2 ≤ p) (hj : j ≤ p - 1)
    (hB : ‖B‖ = ‖π‖ ^ (d + (p - 1)))
    (hsep : ‖B * π ^ j‖ ≠ ‖R‖) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖B * π ^ j + R‖ :=
  LossTwoPSubTwo.loss_le_two_p_sub_two_of_component hπ0 hπ1 hp hj hB
    (norm_le_of_norm_ne hsep)

/-- `0 < ‖π‖ < 1` で指数が違えば `‖π‖` の冪も違う。 -/
theorem norm_ne_of_exp_ne {M : Type*} [NormedField M] {π a r : M} {v1 v2 : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (ha : ‖a‖ = ‖π‖ ^ v1) (hr : ‖r‖ = ‖π‖ ^ v2)
    (hv : v1 ≠ v2) : ‖a‖ ≠ ‖r‖ := by
  rw [ha, hr]
  rcases lt_or_gt_of_ne hv with h | h
  · exact ne_of_gt (pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 h)
  · exact ne_of_lt (pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 h)

/-- ★★★★★**今日の到達点** —— `loss ≤ 2p−2` の仮説が
★**ただ 1 つの自然数の不等式 `d + (p−1) + j ≠ v`** になった。

分岐の言葉では「成分 `j` の位と残りの部分の位が一致しない」で、
`v(B·π^j) ≡ j (mod p)`、`v(R) ≢ j (mod p)` から出る古典。
★`SlotStructure.slot_le`（本日証明済み）と同じ枠の議論である。 -/
theorem loss_le_of_exp_ne {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π B R : M} {p d j v : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p) (hj : j ≤ p - 1)
    (hB : ‖B‖ = ‖π‖ ^ (d + (p - 1)))
    (hR : ‖R‖ = ‖π‖ ^ v) (hv : d + (p - 1) + j ≠ v) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖B * π ^ j + R‖ := by
  have hBj : ‖B * π ^ j‖ = ‖π‖ ^ (d + (p - 1) + j) :=
    ComponentWitness.norm_component_eq_pow hB
  exact loss_le_of_separated (le_of_lt hπ0) (le_of_lt hπ1) hp hj hB
    (norm_ne_of_exp_ne hπ0 hπ1 hBj hR hv)

end Separated

/-! ## §3 使っている公理の一覧 -/

#print axioms norm_le_of_norm_ne
#print axioms loss_le_of_separated
#print axioms norm_ne_of_exp_ne
#print axioms loss_le_of_exp_ne

end SeparatedComponent

end ABC3.Found.PGC
