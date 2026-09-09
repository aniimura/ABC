import ABC3.Found.PGC.BinomialFirstOrder

/-!
# [pGC] ★★★★成分抽出を**仮説 1 本**にした —— 5 段が Lean に載り、残りが名前つきの 1 点に

## ★★★訂正 —— 「相対の冪基底は無い」は**私の測り方が狭かった**

前波まで私は「`𝒪_L = 𝒪_{E₁}[π]` の相対の冪基底は mathlib に見つからない」と書いていた。
★本波で木を開いたら、★★**ノルムの言葉で木に在った**:

- `TotallyRamifiedLayer.lean:133 adjoin_eq_top_of_valK`
  —— `[Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]` で
  「底の値群が `‖π‖^{nℤ}` なら `M = K(π)`」。★**`PAdicLocalField` を使っていない。**
- 同ファイルの `linearIndependent_of_ne_mod` —— `π^0,…,π^{n−1}` の 1 次独立性。
- `JumpStrictMono.lean:116 exists_coeff_norm_le` —— ★**展開 `z = Σ c_l π^l` と
  `‖c_l π^l‖ ≤ ‖z‖`**。

★★**`exists_coeff_norm_le` の結論こそ、本波が要っていたもの**である
（★私は以前「向きが逆」と書いたが、それは同ファイル `:167` の
`norm_sub_apply_le_mul` についてであって、★`:116` ではなかった）。
★**「mathlib に無い」は正しかったが、「木にも無い」は測っていなかった。**

## ★★★本波の定理 —— 成分 1 つで上界が出る

> `loss_le_of_component_witness` :
> `‖B‖ = ‖π‖^D`、`j ≤ p−1`、★**`‖B·π^j‖ ≤ ‖z‖`** ならば `‖π‖^{D+(p−1)} ≤ ‖z‖`

★仮説 `hz : ‖B·π^j‖ ≤ ‖z‖` は **`exists_coeff_norm_le` の結論そのもの**である。
★`‖B‖ = ‖π‖^D` は「`j₀` が単数だから主項が消えない」（`TopIndexSurvives.not_dvd_of_pos_lt`）
の帰結である。

⇒ ★★★**前波の 5 段が、これで Lean の上に並んだ:**

| 段 | Lean |
|---|---|
| (0) `d ≡ 1 (mod p)` が必要 | `TailNoCancel.dvd_add_sub_one_iff`（**定理**） |
| (1) `ρ = w + wπ`、`v(w) = p` | `RhoFactorization` / `ZetaUnitFactor`（**定理**） |
| (2) 主項の同定・残りは 2 次 | `BinomialFirstOrder.norm_add_pow_sub_linear_le`（**定理**、前波） |
| (3) `j₀` は単数 | `TopIndexSurvives.not_dvd_of_pos_lt`（**定理**） |
| (4) `δ` は深い（`p ≥ 3`） | `TopIndexSurvives.delta_deeper_iff_three_le`（**定理**） |
| (5) 上界の算術 | `TopIndexSurvives.loss_bound_from_index`（**定理**） |
| ★成分から結論へ | ★**本波 `loss_le_of_component_witness`（定理、仮説 1 本）** |

## ★★残っているのはただ 1 つ（★名前がついた）

> ★**仮説 `hz` を `exists_coeff_norm_le` から実際に供給すること。**

そのために要るのは、★`σx − x` の `π`-展開の係数 `c_{j₀}` が
★**代数から計算した `B_{j₀} = w·j₀·f_{j₀} + (深い項)` と一致する**こと
（＝展開の**一意性**）。★`linearIndependent_of_ne_mod` がその材料である。
★★**この 1 点だけが未着手**である（★費用は書かない。★開いたが試していない）。

## ★★★次の波が見落とさないための注記（★2 波連続で書く）

★★**5 段の証明の本文は `TopIndexSurvives.lean` の docstring にある。**
★`loss ≤ 2p−2` を**そのまま述べた** Lean 定理はまだ無い ——
本波の `loss_le_of_component_witness` が**その最後の 1 段**を、
★仮説 `hz` つきで載せた形である。

★母集団: **`p ≥ 3`**。`p = 2` は `2p−1 = 3` が上界で、
理由も証明済み（`TwoIsDegenerate`、`delta_deeper_iff_three_le`）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**本波は仮説 `hz` を証明していない。**★「`loss ≤ 2p−2` が定理になった」とは書かない。
4. ★★**自己訂正**: 「相対の冪基底は無い」は**測り方が狭かった** ——
   mathlib には無いが、★**木にノルムの言葉で在る**（`TotallyRamifiedLayer` /
   `JumpStrictMono`）。★「無い」と書くときは**どこを測ったか**を添えること
   （★本日の「『偽』には母集団を添える」と同じ規律）。
5. ★宣言名 2 件は **実ファイル**で衝突 0 を確かめた（#348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace ComponentWitness

/-! ## §1 成分 1 つから損失の上界が出る -/

section Witness

/-- `‖B‖ = ‖π‖^D` なら `‖B·π^j‖ = ‖π‖^{D+j}`。 -/
theorem norm_component_eq_pow {M : Type*} [NormedField M] {π B : M} {j D : ℕ}
    (hB : ‖B‖ = ‖π‖ ^ D) : ‖B * π ^ j‖ = ‖π‖ ^ (D + j) := by
  rw [norm_mul, hB, norm_pow, ← pow_add]

/-- ★★★★**5 段の証明の最後の 1 段**（★仮説 `hz` つき）。

成分 `j`（`j ≤ p−1`）の係数 `B` が `‖B‖ = ‖π‖^D` を満たし、
★`‖B·π^j‖ ≤ ‖z‖`（＝ `JumpStrictMono.exists_coeff_norm_le` の結論）なら
`‖π‖^{D+(p−1)} ≤ ‖z‖`、すなわち ★**`v(z) ≤ D + (p−1)`**。

★`d = D − (p−1)` なので `loss = v(z) − d ≤ 2p−2`
（`TopIndexSurvives.loss_bound_from_index`）。

★★残るのは `hz` を実際に供給すること —— `σx − x` の展開の係数が
代数から計算した `B_{j₀}` と一致すること（展開の一意性）。★これだけが未着手。 -/
theorem loss_le_of_component_witness {M : Type*} [NormedField M]
    {z π B : M} {p j D : ℕ} (hπ0 : 0 ≤ ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    (hj : j ≤ p - 1) (hB : ‖B‖ = ‖π‖ ^ D) (hz : ‖B * π ^ j‖ ≤ ‖z‖) :
    ‖π‖ ^ (D + (p - 1)) ≤ ‖z‖ := by
  have h1 : ‖B * π ^ j‖ = ‖π‖ ^ (D + j) := norm_component_eq_pow hB
  have h2 : ‖π‖ ^ (D + (p - 1)) ≤ ‖π‖ ^ (D + j) :=
    pow_le_pow_of_le_one hπ0 hπ1 (by omega)
  rw [h1] at hz
  linarith

end Witness

/-! ## §2 使っている公理の一覧 -/

#print axioms norm_component_eq_pow
#print axioms loss_le_of_component_witness

end ComponentWitness

end ABC3.Found.PGC
