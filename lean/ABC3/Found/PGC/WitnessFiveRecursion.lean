import ABC3.Found.PGC.AdmissibleClosedForm

/-!
# [pGC] ★★(L) は**木からは測れない** ＋ ★★★本日の区切りの記録（棚卸し最新版）

## ★(L) を開いた結果 —— `(未測定)` は「まだ測っていない」ではなく「木からは測れない」

`GainedTowerDescent.lean:87` の表の最終列「既存ファイルの実測 `L`」が
旧反例の行だけ `(未測定)` である。★開いて分かったこと:

- 他の 2 行の `L`（`8` / `26`）は `EquivariantProjectionDescent` が
  **`ℚ₃(ζ₂₇)` と `ℚ₃(ζ₈₁)` という実在の体**で測った値である。
- ★旧反例 `(p, e_K, S, T) = (5, 3, 2, 17)` は
  `JumpDefectTradeoff` 由来の**パラメータの組**であって、★**対応する体が木に無い**。
  （`e_K = 3` は ℚ₅ 上の分岐指数 3 の底、`e_L = 15` の完全分岐拡大が要る。）

⇒ ★★**「実測 `L`」はその行については定義できない。**「未測定」ではなく
★**「木からは測れない」**が正しい書き方である。★2 波続けて「残り」に挙げていたが、
★**開いてみたら測る対象が存在しなかった。**

★★在庫の測定（#348 が 2 波連続で効いた）: 代わりに書こうとした
「タプルの許容性」は ★**既に木にあった** ——
`JumpDefectTradeoff.lean:304 witness_five_admissible` が
`(p−1)·2 ≤ 5·3`、`(p−1)·17 ≤ 5²·3`、★`17 = 2 + 5·3`（安定枝の関係）、
`¬((p−1)·2 ≤ 3)` を**すべて証明済み**である。★書かずに済んだ。

## ★本波が足した 1 点 —— 漸化式と深さ 2 の公式の突き合わせ

`witness_five_gainedLoss` : `gainedLoss 5 tFive 2 = 17`。
★`sharpTwoCost 5 2 17 = 17`（`GainedDescent.Numeric.witness_five_sharp`）と**一致**する。
⇒ ★**一般の漸化式が、深さ 2 の閉じた公式を旧反例で再現する**ことの検査。
★これは今まで木に無かった（`witness_five_*` は `sharpTwoCost` 側だけを見ていた）。

## ★★★★区切りの記録 —— 棚卸し（2026-09-09 現在。★`ZetaUnitFactor.lean` §0 の更新版）

### 穴の現状

| 穴 | 状態 |
|---|---|
| ①不分岐 | (b) **閉**／(a)「不分岐 ⇒ `∃a, Tr(a)=1`」残（★測定済み・重い） |
| ②′一様定数 | ★**測れた**（目盛りで測ると深さに依らない） |
| ③幾何減衰 | ★**支持された**（鋭い定数も公比 `1/p` の等比列） |
| ④測定点 | **消えた**（積では通る） |
| ⑤出口の限界 | ★**確定**（`axDecay` は 1 目盛り足りない） |

### `hform` の鎖（★**16 波**。§0 の 12 波から 4 波追加）

| # | ファイル | 結論 |
|---|---|---|
| 1–12 | （`ZetaUnitFactor.lean` §0 の表を参照） | `hform` の破れ → 空き枠 → `p=2` の退化 → 柱が仮定ゼロ |
| 13 | `ResidueUnitNorm` | ★`p ∤ b ⇒ ‖b‖ = 1` を **Bézout だけ**で（#69 に触れず）。★柱の残る仮定は `‖p‖<1` の**設定宣言**のみ |
| 14 | `GainedLossClosedForm` | ★閉じた形は**任意の正列では偽**（22,620 本中一致 0.1〜0.6%）。正しい仮説 `hbr` の下で**証明** |
| 15 | `AdmissibleClosedForm` | ★★許容列（`p·t_j ≤ t_{j+1}`）では `hbr` が**自動**。⇒ 木の閉じた形は**無条件に正しかった** |
| 16 | ★本ファイル | (L) は木からは測れない。漸化式が深さ 2 の公式を再現することを確認 |

### 掘らなくてよいと確定した道（4 本、変化なし）

`AxWildDescentDecay`（循環）／予算関数 `F` の載せ替え／`FirstJumpRoute` の再パラメータ化／B2。

### ★本日の方法論（本日この鎖で確立し、台帳に載ったもの）

1. ★**費用を先に測ってから選ぶ**（数学の費用は 6 波連続で当たった）。
2. ★★**開いていないものは「未測定」と書き、費用は書かない**
   （配管の見積もりを 3 度外したあとの規律）。
3. ★★**見積もらずに開いてから決める**（4 波連続で機能。★本波もこれで (L) が
   「測れない」と分かった）。
4. ★**法則を立てたら 3 つ目の場合で測る**（#350。本日の偶然の一致 4 件）。
5. ★★**「偽」と書くときは母集団を添える**（#14 波と #15 波の対で確立）。
6. ★宣言名は**実ファイル**で数える（#348。★本波を含め 3 回、既存宣言を見つけて
   書かずに済んだ／改名した）。

### ★残り（★見積もりは書かない）

- ①(a)「不分岐 ⇒ `∃a, Tr(a)=1`」 —— ★測定済み・重い。
- ②′の一様な証人 / ③ / ⑤ —— ★「塔の外を排除する」問題。`KrasnerCeiling` で
  Krasner では原理的に届かないと定理化済み。
- ★柱が仮定ゼロになったことが `AxWildDescent` の側に何を与えるか —— ★**未測定**。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
   ★`GainedTowerDescent.lean:87` の `(未測定)` も**直していない** ——
   「木からは測れない」は**本ファイルに書く**。
3. ★`tFive` は `t 1 = 2`, `t 2 = 17` だけを決める最小の関数である
   （`k = 2` の計算に他の値は使われない）。★許容列の仮説は**課していない**
   （`gainedLoss` を直接計算しているだけなので不要）。
4. ★宣言名 3 件は **実ファイル**で衝突 0 を確かめた（#348）。
   ★書こうとした `witness_five_admissible` は**既に木にあった**ので書かなかった。
5. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace WitnessFiveRecursion

open GainedDescent

/-! ## §1 旧反例 `(p, e_K, S, T) = (5, 3, 2, 17)` を漸化式で回す -/

section Recursion

/-- 旧反例の跳びの列（`t 1 = 2`, `t 2 = 17`）。★`k = 2` の計算には他の値は使われない。 -/
def tFive : ℕ → ℤ := fun j => if j = 1 then 2 else 17

/-- ★★旧反例 `(5,3,2,17)` を**一般の漸化式**で回すと `Λ_2 = 17`。 -/
theorem witness_five_gainedLoss : gainedLoss 5 tFive 2 = 17 := by
  norm_num [gainedLoss, jumpSum, tFive]

/-- ★★★**漸化式と深さ 2 の閉じた公式が一致する**（どちらも `17`）。

★`GainedDescent.Numeric.witness_five_sharp` は `sharpTwoCost` 側だけを見ていた。
★本波は `gainedLoss` 側から回して同じ値になることを確かめた。 -/
theorem witness_five_recursion_matches_sharp :
    gainedLoss 5 tFive 2 = sharpTwoCost 5 2 17 := by
  rw [witness_five_gainedLoss]
  norm_num [sharpTwoCost]

end Recursion

/-! ## §2 使っている公理の一覧 -/

#print axioms witness_five_gainedLoss
#print axioms witness_five_recursion_matches_sharp

end WitnessFiveRecursion

end ABC3.Found.PGC
