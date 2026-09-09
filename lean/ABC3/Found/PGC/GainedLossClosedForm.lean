import ABC3.Found.PGC.GainedTowerDescent

/-!
# [pGC] ★★★`gainedLoss` の閉じた形は**一般には偽** —— 正しい仮説を特定して証明した

## ★どう選んだか（★本波の形 —— 見積もらずに**開いてから決めた**）

前波で私は「`GainedTowerDescent.lean:94` は**未測定**（開いていません）」とだけ書いた。
★本波はまず**開いた**。同ファイルの測定 4 に**未証明の 1 点が明示**されていた:

> ★閉じた形も出た（★**4,194 本 → 不一致 0 件だが本ファイルは証明していない**）:
> `Λ_k = t_k + Σ_{j<k} t_j − t_1·(p^{k−1}−1)/(p−1)`

★**明示された未証明**なので、そこを測った。所要 3 往復。

## ★★測定（`tools/gainedloss-closed-check.py`、★本波で新規作成）

漸化式（`GainedTowerDescent.lean:184`）:
`Λ_{m+1} = max( A_{m+1} , max(0, A_{m+1} − t_1) + p·Λ_m )`、`A_{m+1} = t_{m+1} − (p−1)·J_m`。

`t_j ∈ {1,…,12}`、`k ≤ 4` を総当たり（各素数 **22,620 本**）:

| p | 一致 | 割合 | ★不一致 |
|---|---|---|---|
| 2 | 131 | 0.6% | **22,489** |
| 3 | 40 | 0.2% | **22,580** |
| 5 | 23 | 0.1% | **22,597** |

⇒ ★★**閉じた形は一般には偽である。**最小の反例は `p=3, t ≡ 1, k=2`:
`Λ_2 = 3` に対し閉じた形は `1`（`closed_form_needs_hypothesis`）。

★★**元の docstring は誤りではない** —— 「4,194 本」は**許容列に制限した母集団**だったはずである。
★ただし**無条件に読める書き方**なので、ここに条件を明示する
（CLAUDE.md「訂正は自分のファイルに書く」）。

## ★★★特定した条件（★#350 に従い 3 つ目の素数まで測った）

> ★**2 段目以降のすべての段で「第 2 枝が勝ち」かつ「`A_m − t_1 ≥ 0`」**

| p | 条件を満たすもの: 一致 | 不一致 |
|---|---|---|
| 2 | 131 | ★**0** |
| 3 | 40 | ★**0** |
| 5 | 23 | ★**0** |

★3 素数すべてで**例外なし**。⇒ この条件が正しい仮説である。

## ★★証明（§2）

条件を**等式の形**で仮説にした:

`hbr : ∀ m, Λ_{m+2} = (t_{m+2} − (p−1)·J_{m+1} − t_1) + p·Λ_{m+1}`

★これは「第 2 枝が勝ち、かつ `max(0,·)` が正」を 1 本の式で書いたもので、
★**測定でそのまま検査できる形**である。この下で

> `gainedLoss_closed` : `Λ_{k+1} = t_{k+1} + J_k − t_1·geomLow p k`

を **`k` についての帰納法**で証明した。★核は `geomLow p (j+1) = 1 + p·geomLow p j`
（閉じた形の `(p^{k−1}−1)/(p−1)` を**割り算を使わずに**書いた形）だけである。

帰納段の計算:
`Λ_{k+2} = (t_{k+2} − (p−1)J_{k+1} − t_1) + p(t_{k+1} + J_k − t_1·g_k)`
で `J_{k+1} = J_k + t_{k+1}` と `g_{k+1} = 1 + p·g_k` を入れると
`= t_{k+2} + J_{k+1} − t_1·g_{k+1}` になる（`ring` 1 行）。

## ★残り（★見積もりは書かない）

1. ★条件（`hbr`）が**許容列（Hasse–Arf の跳びの列）で自動的に成り立つか**は
   ★**未測定**である。★元の docstring の「4,194 本」がその母集団なら成り立つはずだが、
   ★私はその母集団の定義を**読んでいない**。
2. ★`GainedTowerDescent.lean:94` の表の「旧反例 `p=5, e_K=3, S=2, T=17`」の
   **既存ファイルの実測 `L`** は同ファイルで `(未測定)` のまま。★本波でも測っていない。
3. ①(a) は測定済みで重い。②′/③/⑤ は変化なし。
4. ★棚卸しは `ZetaUnitFactor.lean` の §0（3 波前）にある。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
   ★`GainedTowerDescent.lean` の docstring も**直していない** —— 他が読んでいるので、
   条件の明示は**本ファイルに書いた**。
3. ★`open GainedDescent` を書いた（`gainedLoss` / `jumpSum` はその名前空間にある）。
   ★最初 `Unknown identifier 'gainedLoss'` が出た —— #68 の「import していない」ではなく
   **名前空間を開いていない**形だった。★同じエラー文でも原因が 2 通りある。
4. ★宣言名 3 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（#348）。
5. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace GainedLossClosedForm

open GainedDescent

/-! ## §1 `1 + p + ⋯ + p^{j−1}` -/

section Geom

/-- `1 + p + ⋯ + p^{j−1}`。★閉じた形の `(p^{k−1}−1)/(p−1)` を**割り算を使わずに**書いた形。 -/
def geomLow (p : ℤ) : ℕ → ℤ
  | 0 => 0
  | (j + 1) => 1 + p * geomLow p j

/-- ★帰納法の核 —— `geomLow p (j+1) = 1 + p·geomLow p j`。 -/
theorem geomLow_succ (p : ℤ) (j : ℕ) : geomLow p (j + 1) = 1 + p * geomLow p j := rfl

end Geom

/-! ## §2 閉じた形（★仮説付き） -/

section Closed

/-- ★★★**`GainedTowerDescent.lean:94` が「証明していない」と書いた閉じた形。**

仮説 `hbr` は「2 段目以降のすべての段で**第 2 枝が勝ち**、かつ `A_m − t_1 ≥ 0`」を
1 本の等式で書いたもの。★測定でそのまま検査できる形である。

★この仮説の下で `Λ_{k+1} = t_{k+1} + J_k − t_1·geomLow p k` を**帰納法で証明**した。
★測定: 条件を満たす列は `p = 2,3,5` で一致 131/40/23、**不一致 0**。 -/
theorem gainedLoss_closed {p : ℤ} {t : ℕ → ℤ}
    (hbr : ∀ m, gainedLoss p t (m + 2)
      = (t (m + 2) - (p - 1) * jumpSum t (m + 1) - t 1) + p * gainedLoss p t (m + 1))
    (ht1 : 0 ≤ t 1) :
    ∀ k, gainedLoss p t (k + 1) = t (k + 1) + jumpSum t k - t 1 * geomLow p k
  | 0 => by
      rw [gainedLoss_one ht1]
      simp [jumpSum, geomLow]
  | (k + 1) => by
      have IH := gainedLoss_closed hbr ht1 k
      rw [hbr k, IH]
      simp only [jumpSum, geomLow]
      ring

end Closed

/-! ## §3 仮説を外すと偽 -/

section Counterexample

/-- ★★**仮説を外すと偽** —— `p = 3`, `t ≡ 1`, `k = 2` で `Λ_2 = 3` だが閉じた形は `1`。

★総当たり（`t_j ∈ {1..12}`, `k ≤ 4`, 各素数 22,620 本）で一致は
`p=2` で 0.6%、`p=3` で 0.2%、`p=5` で 0.1% しかない。
⇒ ★**閉じた形は一般には偽**であり、仮説が要る。 -/
theorem closed_form_needs_hypothesis :
    gainedLoss 3 (fun _ => 1) 2 ≠ (fun _ : ℕ => (1 : ℤ)) 2
      + jumpSum (fun _ => 1) 1 - (fun _ : ℕ => (1 : ℤ)) 1 * geomLow 3 1 := by
  norm_num [gainedLoss, jumpSum, geomLow]

end Counterexample

/-! ## §4 使っている公理の一覧 -/

#print axioms geomLow_succ
#print axioms gainedLoss_closed
#print axioms closed_form_needs_hypothesis

end GainedLossClosedForm

end ABC3.Found.PGC
