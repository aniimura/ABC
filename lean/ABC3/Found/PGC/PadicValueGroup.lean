import ABC3.Found.PGC.ValKFromBase
import Mathlib.NumberTheory.Padics.PadicNumbers

/-!
# [pGC] 最後の外部入力が定理になった —— `ℚ_p` の値群と絶対次数

## 持ち場（前波の表の 2 行）

★**唯一の外部入力**「`ℚ_p` の値群が `‖p‖^ℤ`」と、未実装の「`finrank ℚ_p L = n`」。

## ★結果 —— **両方とも定理になった**

| 項目 | 状態 | 中身 |
|---|---|---|
| `ℚ_p` の値群 | ★**定理**（§1 `padic_value_group`） | `Padic.norm_eq_zpow_neg_valuation` ＋ `Padic.norm_p` ＋ `inv_zpow`/`zpow_neg` の **1 行** |
| 像に移す | §2 `valgroup_image` | `[NormedAlgebra]` なら `norm_algebraMap'` で等長 |
| `hvalbase` | ★**定理**（§3 `valbase_padic`） | 要るのは `‖π‖^n = ‖p‖` だけ |
| `hvalK` | ★★**定理**（§3 `valK_layer_padic`） | 入力は「2 つの次数」と `‖π‖^n = ‖p‖` |
| `finrank ℚ_p L = n` | ★**定理**（§3 `finrank_padic_eq`） | `PowerSpanTop`（前々波）＋ `valbase_padic`。★循環しない |

★`ℚ_p` の値群は **1 行**だった。★「唯一の外部入力」と身構えたが、
`Padic.norm_eq_zpow_neg_valuation`（`x ≠ 0 → ‖x‖ = (p:ℝ)^{-x.valuation}`）と
`Padic.norm_p`（`‖(p:ℚ_p)‖ = (p:ℝ)⁻¹`）を繋ぐだけだった。★索引に `Padic.norm_p` は
**在った**（`grep "‖(p : ℚ_\[p\])‖" .cache/mathlib-index.txt` で当たった。
★名前ではなく**結論のリテラル**で引いた）。

## ★★★残っているもの（★数学の入力はもう無い）

`valK_layer_padic` / `finrank_padic_eq` に渡すものを全部並べる:

| 入力 | 種類 | 状態 |
|---|---|---|
| `‖π‖^n = ‖p‖` | 数学 | ★`ZetaSubOnePrime.norm_zeta_sub_one_pow`（自分）✓ |
| `‖μ‖ = ‖π‖^p` | 数学 | ★`ZetaStepRatio.norm_zeta_step`（自分）✓ |
| 次数 `n` のモニック関係式 | 数学 | 円分多項式（`ZetaSubOnePrime` が Eisenstein で使った材料）★未実装 |
| `Algebra.adjoin ℚ_p {π} = ⊤` | 定義 | 塔の定義 |
| `q ≠ 1`（`L ≠ E₁`） | 定義 | 塔の定義 |
| `NormedAlgebra` / `IsScalarTower` / `FiniteDimensional` ×3 | ★**インスタンス** | ★模型の構成 |

⇒ ★**残るのは「関係式 1 本」と「模型の構成（インスタンス）」だけ**である。

## ★配管の記録（`#print axioms` の読み方）

型検査が落ちた宣言でも `leanfile.mjs --full` は `#print axioms` の行を出す。★そのとき
**`sorryAx` が混じる**（エラー回復で `sorry` が入るため）。実際に出た文（逐語）:

```
'ABC3.Found.PGC.PadicValueGroup.valbase_padic' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
```

★`sorryAx` を見たら、★**その上のエラーを先に読む**（本当に `sorry` を書いたわけではない）。
本ファイルは修正後に 5 件すべて `[propext, Classical.choice, Quot.sound]` を確認した。

## 逸脱の記録

- §2 は `[NormedAlgebra K M]` ＋ `[NormOneClass M]` を仮定する。★具体層では
  `‖(1 : L)‖ = 1` は体のノルムから自動だが、★インスタンスとしては**模型の構成**に属する。
- §1 の `padic_value_group` は `{x}` 暗黙なので、`∀ a, a ≠ 0 → …` の形に渡すときは
  ★`fun _ ha => padic_value_group ha` と**η 展開**が要る（#356 の親戚）。
-/

namespace ABC3.Found.PGC

namespace PadicValueGroup

/-! ## §1 ★`ℚ_p` の値群は `‖p‖^ℤ`（唯一の外部入力） -/

section Padic

variable {p : ℕ} [Fact p.Prime]

/-- ★★**唯一の外部入力**だったもの。`Padic.norm_eq_zpow_neg_valuation` と
`Padic.norm_p` から 3 行で出る。 -/
theorem padic_value_group {x : ℚ_[p]} (hx : x ≠ 0) :
    ∃ m : ℤ, ‖x‖ = ‖(p : ℚ_[p])‖ ^ m := by
  refine ⟨x.valuation, ?_⟩
  rw [Padic.norm_eq_zpow_neg_valuation hx, Padic.norm_p, inv_zpow, zpow_neg]

end Padic

/-! ## §2 像に移す（`NormedAlgebra` なら等長） -/

section Image

/-- `[NormedAlgebra K M]` かつ `‖(1 : M)‖ = 1` なら `algebraMap` は等長なので、
`K` の値群の主張がそのまま像に移る。 -/
theorem valgroup_image {K M : Type*} [NormedField K] [NormedField M] [NormedAlgebra K M]
    [NormOneClass M] {P : K}
    (hK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖a‖ = ‖P‖ ^ m) :
    ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖algebraMap K M P‖ ^ m := by
  intro a ha
  obtain ⟨m, hm⟩ := hK a ha
  rw [norm_algebraMap', norm_algebraMap']
  exact ⟨m, hm⟩

end Image

/-! ## §3 ★`hvalbase` が定理になった -/

section Base

variable {p : ℕ} [Fact p.Prime]

/-- ★★底 `ℚ_p` の条件が**定理**になった。要るのは `‖π‖^n = ‖p‖` の 1 本だけ
（＝ `ZetaSubOnePrime.norm_zeta_sub_one_pow`、`n = φ(p^{m+1})`）。 -/
theorem valbase_padic {M : Type*} [NormedField M] [NormedAlgebra ℚ_[p] M] [NormOneClass M]
    {π : M} {n : ℕ} (hpn : ‖π‖ ^ n = ‖algebraMap ℚ_[p] M (p : ℚ_[p])‖) :
    ∀ a : ℚ_[p], a ≠ 0 → ∃ m : ℤ, ‖algebraMap ℚ_[p] M a‖ = ‖π‖ ^ ((n : ℤ) * m) :=
  ValKFromBase.valbase_of_uniformizer hpn
    (valgroup_image (fun _ ha => padic_value_group ha))

/-- ★★★**`hvalK` が円分の材料だけから出る**。
入力は「2 つの次数」と「`‖π‖^n = ‖p‖`」だけ。 -/
theorem valK_layer_padic {E M : Type*} [Field E] [NormedField M] [IsUltrametricDist M]
    [NormedAlgebra ℚ_[p] M] [NormOneClass M] [Algebra ℚ_[p] E] [Algebra E M]
    [IsScalarTower ℚ_[p] E M] [FiniteDimensional ℚ_[p] M] [FiniteDimensional E M]
    [FiniteDimensional ℚ_[p] E]
    {π : M} {n q : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hn : Module.finrank ℚ_[p] M = n) (hq : Module.finrank E M = q)
    (hpn : ‖π‖ ^ n = ‖algebraMap ℚ_[p] M (p : ℚ_[p])‖) :
    ∀ a : E, a ≠ 0 → ∃ m : ℤ, ‖algebraMap E M a‖ = ‖π‖ ^ ((q : ℤ) * m) :=
  ValKFromBase.valK_of_base_layer hπ0 hπ1 hn hq (valbase_padic hpn)

/-- ★★**絶対次数も定理になった** —— `‖π‖^n = ‖p‖` ＋ 次数 `n` のモニック関係式
＋ `L = ℚ_p(π)` から `finrank ℚ_p L = n`。★循環しない（`valbase_padic` は次数を要求しない）。 -/
theorem finrank_padic_eq {M : Type*} [NormedField M] [IsUltrametricDist M]
    [NormedAlgebra ℚ_[p] M] [NormOneClass M] [FiniteDimensional ℚ_[p] M]
    {π : M} {n : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hpn : ‖π‖ ^ n = ‖algebraMap ℚ_[p] M (p : ℚ_[p])‖)
    (c : ℕ → ℚ_[p])
    (hrel : π ^ n = ∑ j ∈ Finset.range n, algebraMap ℚ_[p] M (c j) * π ^ j)
    (hadj : Algebra.adjoin ℚ_[p] ({π} : Set M) = ⊤) :
    Module.finrank ℚ_[p] M = n :=
  PowerSpanTop.finrank_eq_of_relation hπ0 hπ1 (valbase_padic hpn) c hrel hadj

end Base

/-! ## §4 使っている公理の一覧 -/

#print axioms padic_value_group
#print axioms valgroup_image
#print axioms valbase_padic
#print axioms valK_layer_padic
#print axioms finrank_padic_eq

end PadicValueGroup

end ABC3.Found.PGC
