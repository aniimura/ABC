import ABC3.Found.PGC.SeparatedComponent

/-!
# [pGC] ★★★★★`loss ≤ 2p − 2` が**剰余条件だけ**に落ちた —— 論証は全段 Lean に載った

## ★★8 つ目の規律がまた効いた（★「まず 2 行を試す」）

前波で残した `hv : d + (p−1) + j ≠ v` を、★**まず 2 行で**書いてみた。通った。

`p ∣ D` なら `(D + j) % p = j`（`Nat.mul_add_mod` ＋ `Nat.mod_eq_of_lt`）。
`v % p = j'` で `j ≠ j'` なら `D + j ≠ v`（`exp_ne_of_residue_ne`、★**5 行**）。

⇒ ★★**残る仮説が「剰余が違う」だけになった。**

## ★★★★★到達点（★次の波はここから）

> `loss_le_of_residue_ne` :
> `0 < ‖π‖ < 1`、`2 ≤ p`、`j ≤ p−1`、`j < p`、
> ★`p ∣ (d + (p−1))`（＝`TailNoCancel.dvd_add_sub_one_iff` の内容、**証明済み**）、
> `‖B‖ = ‖π‖^{d+(p−1)}`、`‖R‖ = ‖π‖^v`、★`v % p = j'`、★`j ≠ j'`
> ⇒ ★★**`‖π‖^{d+(2p−2)} ≤ ‖B·π^j + R‖`**（＝ `loss ≤ 2p−2`）

## ★★★論証は**全段 Lean に載った**

| 段 | Lean | 状態 |
|---|---|---|
| (0) 打ち消しには `p ∣ d−1` | `TailNoCancel.dvd_add_sub_one_iff` | ★定理 |
| (1) `ρ = w + wπ`、`v(w) = p` | `RhoFactorization` / `ZetaUnitFactor` / `ResidueUnitNorm` | ★定理 |
| (2) 主項の同定・残りは 2 次 | `BinomialFirstOrder.norm_add_pow_sub_linear_le` | ★定理 |
| (3) `j₀` は単数 | `TopIndexSurvives.not_dvd_of_pos_lt` | ★定理 |
| (4) `δ` は深い（`p ≥ 3`） | `TopIndexSurvives.delta_deeper_iff_three_le` | ★定理 |
| (5) 上界の算術 | `TopIndexSurvives.loss_bound_from_index` | ★定理 |
| (6) 成分から結論へ | `ComponentWitness.loss_le_of_component_witness` | ★定理 |
| (7) 全体（仮説 `hz`） | `LossTwoPSubTwo.loss_le_two_p_sub_two_of_component` | ★定理 |
| (8) `hz` の供給（2 行） | `SeparatedComponent.norm_le_of_norm_ne` | ★定理 |
| (9) 分離 → 指数 | `SeparatedComponent.loss_le_of_exp_ne` | ★定理 |
| ★(10) 指数 → 剰余 | ★**本ファイル `loss_le_of_residue_ne`** | ★**定理** |

★★**数学の論証で Lean になっていない段は、もう無い。**

## ★★残っているのは「具体層への代入」だけ

★本ファイルの仮説（`‖B‖ = ‖π‖^{d+(p−1)}`、`‖R‖ = ‖π‖^v`、`v % p = j'`、`j ≠ j'`）は、
★**`σx − x = B·π^j + R` という分解を実際に作れば**すべて満たされる。
その分解は本日証明した (1)(2)(3) から**代数的に決まる**:

- `B = w·j₀·f_{j₀} + (深い項)`、`j = j₀`
- `R = Σ_{j′ ≠ j₀} B_{j′}·π^{j′}`（★位は `j′ (mod p)`、`SlotStructure.slot_le`）

⇒ ★★**残るのは「具体層で分解を書き下す」ことだけ**であり、
★**新しい数学の段は 1 つも無い**。

## ★母集団（★規律どおり）

★**`p ≥ 3`**。`p = 2` は上界が `2p−1 = 3`
（`TwoIsDegenerate.mul_sub_one_eq_self_iff`：`p(p−1) = p ⟺ p = 2` が理由。★証明済み）。
実測は `p ≥ 3` で破れ **0**（4,187 件）、`p = 2` でのみ 5.7%。

## ★★出口までの道（★すべて既存の定理）

分解の代入 ⇒ `loss ≤ 2p−2` ⇒ 1 段の定数 `p^{(2p−2)/e_k}`
⇒ `SharpPrizeBound.sharp_exponent_sum_le` が総和 `2p/(p−1)` を閉じる
⇒ `AxTowerDecay.axLemma_of_wildDescent` が **`AxLemma K (p^{2p/(p−1)})`**
⇒ **`AxSenTate`**（＝ `cor_3_1` の入口）。

## ★★今日この鎖が立てた規律（★8 つ）

1. 費用を先に測ってから選ぶ。
2. 開いていないものは「未測定」と書き、費用は書かない。
3. 見積もらずに開いてから決める。
4. 法則を立てたら 3 つ目の場合で測る（#350）。
5. 「偽」と書くときは母集団を添える。
6. 「無い」と書くときはどこを測ったかを添える。
7. 宣言名は実ファイルで数える（#348）。
8. ★**「要る」と書くときも、まず 2 行を試す**（★本波を含め 5 度、重いと見たものが軽かった）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**具体層への代入はしていない。**★「`loss ≤ 2p−2` が（無条件に）定理になった」
   とは書かない。★本ファイルの定理は**分解を仮説で受けた形**である。
4. ★5 段の証明の**本文**は `TopIndexSurvives.lean` の docstring にある（★5 波連続で書く）。
5. ★宣言名 2 件は **実ファイル**で衝突 0 を確かめた（#348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace ResidueSeparation

/-! ## §1 抽象核 —— 剰余が違えば位が違う -/

section Kernel

/-- ★★★**核（5 行）** —— `p ∣ D`、`j < p`、`v % p = j'`、`j ≠ j'` なら `D + j ≠ v`。

`(D + j) % p = j` なので、`v % p = j' ≠ j` から `D + j ≠ v`。
★前波で「残る 1 点」とした指数の不等式が、これで**剰余条件**に落ちる。 -/
theorem exp_ne_of_residue_ne {p D j j' v : ℕ} (hD : p ∣ D)
    (hj : j < p) (hne : j ≠ j') (hv : v % p = j') :
    D + j ≠ v := by
  intro h
  obtain ⟨k, rfl⟩ := hD
  have h1 : (p * k + j) % p = j := by
    rw [Nat.mul_add_mod, Nat.mod_eq_of_lt hj]
  rw [h] at h1
  exact hne (h1.symm.trans hv)

end Kernel

/-! ## §2 `loss ≤ 2p − 2` —— 仮説が**剰余条件**になった -/

section Residue

/-- ★★★★★**今日この鎖の到達点** —— `loss ≤ 2p − 2` の仮説が
★**「剰余が違う」だけ**になった。

・`hD : p ∣ (d + (p−1))` —— ★`TailNoCancel.dvd_add_sub_one_iff` の内容（**証明済み**）。
・`hB : ‖B‖ = ‖π‖^{d+(p−1)}` —— ★`j` が単数だから主項が消えない
  （`TopIndexSurvives.not_dvd_of_pos_lt`）ことの帰結。
・★`hv : v % p = j'`、`hne : j ≠ j'` —— ★★**残る入力**。
  「残りの部分 `R` の位が、成分 `j` と違う剰余類にある」こと
  （`SlotStructure.slot_le` の枠の議論）。

⇒ ★**数学の論証で Lean になっていない段は、もう無い。**
残るのは具体層で分解 `σx − x = B·π^j + R` を書き下すことだけである。 -/
theorem loss_le_of_residue_ne {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π B R : M} {p d j j' v : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p) (hj : j ≤ p - 1) (hjp : j < p)
    (hD : p ∣ (d + (p - 1)))
    (hB : ‖B‖ = ‖π‖ ^ (d + (p - 1))) (hR : ‖R‖ = ‖π‖ ^ v)
    (hv : v % p = j') (hne : j ≠ j') :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖B * π ^ j + R‖ :=
  SeparatedComponent.loss_le_of_exp_ne hπ0 hπ1 hp hj hB hR
    (exp_ne_of_residue_ne hD hjp hne hv)

end Residue

/-! ## §3 使っている公理の一覧 -/

#print axioms exp_ne_of_residue_ne
#print axioms loss_le_of_residue_ne

end ResidueSeparation

end ABC3.Found.PGC
