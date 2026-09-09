import ABC3.Found.PGC.ComponentWitness

/-!
# [pGC] ★★★★★`loss ≤ 2p − 2` が **1 つの Lean の型**になった（仮説 1 本つき）

## ★★★★★次の波が最初に読む場所（★到達点はこの 1 本）

> `loss_le_two_p_sub_two_of_component` :
> `‖B‖ = ‖π‖^{d+(p−1)}`、`j ≤ p−1`、★**`‖B·π^j‖ ≤ ‖z‖`** ならば
> ★★**`‖π‖^{d+(2p−2)} ≤ ‖z‖`**（＝ `v(z) ≤ d + (2p−2)`、＝ **`loss ≤ 2p−2`**）

★★**残る仮説は `hz : ‖B·π^j‖ ≤ ‖z‖` ただ 1 本**である。
★他の 5 段はすべて定理になっている（下の表）。

## ★どう選んだか（★開いてから決めた。本日 8 波連続の形）

前波で名前をつけた「展開の一意性」を測った:

- `hz` を `JumpStrictMono.lean:116 exists_coeff_norm_le` から供給するには、
  ★同定理が**存在**しか言わないので、`σx − x` の展開の係数が
  代数から計算した `B_{j₀}` と**一致する**こと（＝一意性）が要る。
- 一意性の材料は `TotallyRamifiedLayer.linearIndependent_of_ne_mod` に在るが、
  ★`Fin n` 添字の `LinearIndependent` から 2 つの表現の一致を出す配管が要る。
  ★**費用は書かない**（開いたが試していない）。

⇒ ★★持ち場が挙げたもう 1 つの道 —— ★**`hz` を仮説のまま残して
`loss ≤ 2p−2` を 1 つの型として書く** —— を取った。
★**今日この鎖の到達点が、これで 1 本の定理の型に見えるようになる**からである。

## ★★★5 段の現在地（★全部 Lean。★仮説は 1 本だけ）

| 段 | Lean | 状態 |
|---|---|---|
| (0) 打ち消しには `d ≡ 1 (mod p)` | `TailNoCancel.dvd_add_sub_one_iff` | ★**定理** |
| (1) `ρ = w + wπ`、`v(w) = p` | `RhoFactorization` / `ZetaUnitFactor` / `ResidueUnitNorm` | ★**定理** |
| (2) 主項の同定・残りは 2 次 | `BinomialFirstOrder.norm_add_pow_sub_linear_le` | ★**定理** |
| (3) `j₀` は単数（`1 ≤ j₀ ≤ p−1`） | `TopIndexSurvives.not_dvd_of_pos_lt` | ★**定理** |
| (4) `δ` は深い（`p < p(p−1) ⟺ p ≥ 3`） | `TopIndexSurvives.delta_deeper_iff_three_le` | ★**定理** |
| (5) 上界の算術 | `TopIndexSurvives.loss_bound_from_index` | ★**定理** |
| (6) 成分から結論へ | `ComponentWitness.loss_le_of_component_witness` | ★**定理** |
| ★(7) 全体 | ★**本ファイル** | ★**定理（仮説 `hz` 1 本）** |
| ★★残り | `hz` の供給（展開の一意性） | ★**未着手** |

## ★★出口までの道（★すべて既存の定理で繋がる）

1. `hz` が供給できれば `loss ≤ 2p−2`（本ファイル）。
2. ⇒ 1 段の定数は `p^{(2p−2)/e_k}`。
3. ⇒ `SharpPrizeBound.sharp_exponent_sum_le` が指数の総和 `2p/(p−1)` を閉じる。
4. ⇒ `AxTowerDecay.axLemma_of_wildDescent` が **`AxLemma K (p^{2p/(p−1)})`**。
5. ⇒ `AxSenTate`、すなわち `cor_3_1` の入口。

★★**残っているのは 1 の `hz` だけ**である。

## ★母集団（★規律どおり添える）

★**`p ≥ 3`**。`p = 2` は上界が `2p−1 = 3` で、
理由も証明済み（`TwoIsDegenerate.mul_sub_one_eq_self_iff`：`p(p−1) = p ⟺ p = 2`）。
★実測: `p ≥ 3` で `loss ≤ 2p−2` の破れ **0**（4,187 件）、`p = 2` でのみ 5.7% 破れる。

## ★★今日この鎖が立てた規律（★7 つ。★次の波はここから始めてよい）

1. **費用を先に測ってから選ぶ**。
2. **開いていないものは「未測定」と書き、費用は書かない**。
3. **見積もらずに開いてから決める**。
4. **法則を立てたら 3 つ目の場合で測る**（#350）。
5. **「偽」と書くときは母集団を添える**。
6. ★**「無い」と書くときはどこを測ったかを添える**。
7. **宣言名は実ファイルで数える**（#348）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**`hz` は証明していない。**★「`loss ≤ 2p−2` が（無条件に）定理になった」とは書かない。
   ★本ファイルの定理は**仮説 1 本つき**である。
4. ★5 段の証明の**本文**は `TopIndexSurvives.lean` の docstring にある
   （★3 波連続で書いている）。
5. ★宣言名 2 件は **実ファイル**で衝突 0 を確かめた（#348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace LossTwoPSubTwo

/-! ## §1 `loss ≤ 2p − 2` —— 仮説 1 本つきの Lean 定理 -/

section Main

/-- ★★★★★**今日この鎖の到達点** —— `loss ≤ 2p − 2` の Lean 定理（★仮説 1 本つき）。

`z = σx − x`、`π` は `L` の素元、`B` は生き残る成分 `j` の係数、
`d = d(x, E₁)` の目盛りとする。

・`hB : ‖B‖ = ‖π‖^{d+(p−1)}` —— ★`j` が単数だから主項が消えない
  （`TopIndexSurvives.not_dvd_of_pos_lt`）ことの帰結。
・★`hz : ‖B·π^j‖ ≤ ‖z‖` —— ★★**唯一残っている仮説**。
  `JumpStrictMono.exists_coeff_norm_le` の結論の形で、供給には展開の一意性が要る。

⇒ `‖π‖^{d+(2p−2)} ≤ ‖z‖`、すなわち **`loss = v(z) − d ≤ 2p − 2`**。

★母集団は **`p ≥ 3`**（`p = 2` は `2p−1`）。 -/
theorem loss_le_two_p_sub_two_of_component {M : Type*} [NormedField M]
    {z π B : M} {p d j : ℕ}
    (hπ0 : 0 ≤ ‖π‖) (hπ1 : ‖π‖ ≤ 1) (hp : 2 ≤ p) (hj : j ≤ p - 1)
    (hB : ‖B‖ = ‖π‖ ^ (d + (p - 1)))
    (hz : ‖B * π ^ j‖ ≤ ‖z‖) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖z‖ := by
  have hrw : d + (2 * p - 2) = (d + (p - 1)) + (p - 1) := by omega
  rw [hrw]
  exact ComponentWitness.loss_le_of_component_witness hπ0 hπ1 hj hB hz

/-- 読み方 —— `v(z) ≤ d + (2p−2)` は `loss = v(z) − d ≤ 2p − 2` のこと。 -/
theorem loss_exponent_value {p d : ℕ} (hp : 2 ≤ p) :
    (d + (2 * p - 2)) - d = 2 * p - 2 := by omega

end Main

/-! ## §2 使っている公理の一覧 -/

#print axioms loss_le_two_p_sub_two_of_component
#print axioms loss_exponent_value

end LossTwoPSubTwo

end ABC3.Found.PGC
