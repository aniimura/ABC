import ABC3.Found.PGC.GainedLossClosedForm

/-!
# [pGC] ★★★★閉じた形は**許容列では無条件に正しい** —— 木の主張は正しかった

## ★どう選んだか（★本波の形 —— 見積もらずに**開いてから**決めた）

前波で私は (H)「`hbr` が許容列で自動的に成り立つかは**未測定**（元の母集団の定義を
読んでいません）」と書いた。★本波はまず**その母集団の定義を開いた**。

`GainedTowerDescent.lean:295-303` の `pmul_le_of_stable` / `pmul_le_of_below` が
両枝から出しているのは ★**`p·t_j ≤ t_{j+1}`** —— これが「許容列」の中身である。
★開いて分かったので、それを仮説にして測った。所要 3 往復。

## ★★測定（`#350` に従い **4 素数**で測った）

`t_j ∈ {1..9}`、`k ≤ 4`、`p·t_j ≤ t_{j+1}` を満たす列だけを取る:

| p | 許容列 | `hbr` 成立 | ★不成立 | 閉じた形一致 |
|---|---|---|---|---|
| 2 | 45 | 45 | ★**0** | 45 |
| 3 | 22 | 22 | ★**0** | 22 |
| 5 | 14 | 14 | ★**0** | 14 |
| 7 | 12 | 12 | ★**0** | 12 |

⇒ ★**許容列では例外なく成り立つ。**（前波は任意の正列で 0.1〜0.6% しか一致しなかった。）

## ★★★証明（§1〜§3。★測定ではなく定理）

1. `t1_le_A` —— ★`A_m := t_m − (p−1)·J_{m−1}` は `A_{m+1} ≥ A_m` で単調、
   かつ `A_2 = t_2 − (p−1)t_1 ≥ p t_1 − (p−1)t_1 = t_1`。⇒ ★**`A_m ≥ t_1`**。
2. `t1_le_gainedLoss` —— `Λ_1 = t_1`、`Λ_{m+2} ≥ A_{m+2} ≥ t_1`。⇒ ★**`Λ_j ≥ t_1`**。
3. `hbr_of_admissible` —— 1 から `max(0, A−t_1) = A−t_1`、
   2 と `p ≥ 2`, `t_1 ≥ 0` から `p·Λ ≥ 2Λ ≥ 2t_1 ≥ t_1` なので
   `(A−t_1) + pΛ ≥ A`。⇒ ★**第 2 枝が必ず勝つ**。
4. `gainedLoss_closed_of_admissible` —— 前波の `gainedLoss_closed` に 3 を代入。

⇒ ★★★**`GainedTowerDescent.lean:94` の「4,194 本 → 不一致 0 件だが本ファイルは
証明していない」は、これで証明になった。**★木の主張は正しかった。

## ★★2 波を合わせた正しい姿

| 母集団 | 閉じた形 |
|---|---|
| 任意の正の列 | ★**偽**（前波 `GainedLossClosedForm.closed_form_needs_hypothesis`） |
| ★許容列（`p·t_j ≤ t_{j+1}`） | ★★**真**（本波、無条件） |

★前波で私は「一般には偽」と書いた。★**それは正しいが、木が語っていた母集団では真**である。
★★**「偽」と書くときは母集団を必ず添えること** —— 本波の自己訂正である。

## ★配管（`lean-idioms.md` #351 に登録した）

`simp only [gainedLoss]` と書いたら**再帰呼び出しまで開かれ**、
仮説 `hpL : t 1 ≤ p * gainedLoss p t (m + 1)` とゴールの展開形が
**別の原子**になって `linarith failed to find a contradiction` になった。
★外側 1 回だけ `rfl` で開いて直した。★逆に `t1_le_A` では
`simp only [jumpSum] at IH ⊢` と**両方**開いて直した。

## ★残り（★見積もりは書かない）

1. ★`GainedTowerDescent.lean:94` の表の「旧反例 `p=5, e_K=3, S=2, T=17`」の
   **既存ファイルの実測 `L`** は `(未測定)` のまま。★本波でも測っていない。
2. ①(a) は測定済みで重い。②′/③/⑤ は変化なし。
3. ★棚卸しは `ZetaUnitFactor.lean` の §0（5 波前）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
   ★`GainedTowerDescent.lean:94` の「証明していない」という記述も**直していない** ——
   証明できたことは**本ファイルに書く**（CLAUDE.md「訂正は自分のファイルに」）。
3. ★`open GainedDescent GainedLossClosedForm` を書いた（前波の #68 の別形）。
4. ★宣言名 4 件は **実ファイル**で衝突 0 を先に確かめた（#348）。
5. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace AdmissibleClosedForm

open GainedDescent GainedLossClosedForm

/-! ## §1 許容性から `A_m ≥ t_1` -/

section AboveT1

variable {p : ℤ} {t : ℕ → ℤ}

/-- ★★許容列（`p·t_j ≤ t_{j+1}`）では `A_m := t_m − (p−1)·J_{m−1}` が `t_1` 以上。
★`A_{m+1} − A_m = t_{m+1} − p·t_m ≥ 0` で単調、底は `A_2 ≥ t_1`。 -/
theorem t1_le_A (hadm : ∀ j, p * t (j + 1) ≤ t (j + 2)) :
    ∀ m, t 1 ≤ t (m + 2) - (p - 1) * jumpSum t (m + 1)
  | 0 => by
      have h := hadm 0
      simp only [jumpSum]
      linarith
  | (m + 1) => by
      have IH := t1_le_A hadm m
      have h := hadm (m + 1)
      simp only [jumpSum] at IH ⊢
      linarith

/-- ★`Λ_j ≥ t_1`（`Λ_1 = t_1`、`Λ_{m+2} ≥ A_{m+2} ≥ t_1`）。 -/
theorem t1_le_gainedLoss (_hp : 2 ≤ p) (ht1 : 0 ≤ t 1)
    (hadm : ∀ j, p * t (j + 1) ≤ t (j + 2)) :
    ∀ m, t 1 ≤ gainedLoss p t (m + 1)
  | 0 => by rw [gainedLoss_one ht1]
  | (m + 1) => by
      have hA := t1_le_A hadm m
      have : t (m + 2) - (p - 1) * jumpSum t (m + 1) ≤ gainedLoss p t (m + 2) := by
        simp only [gainedLoss]
        exact le_max_left _ _
      linarith

end AboveT1

/-! ## §2 許容列では第 2 枝が必ず勝つ -/

section Branch

variable {p : ℤ} {t : ℕ → ℤ}

/-- ★★★**許容列では第 2 枝が必ず勝つ。**

`max(0, A−t_1) = A−t_1`（`t1_le_A`）と `p·Λ ≥ 2Λ ≥ 2t_1 ≥ t_1`（`t1_le_gainedLoss`）から
`(A−t_1) + p·Λ ≥ A`。⇒ 前波が**仮説として置いた** `hbr` が**定理になる**。 -/
theorem hbr_of_admissible (hp : 2 ≤ p) (ht1 : 0 ≤ t 1)
    (hadm : ∀ j, p * t (j + 1) ≤ t (j + 2)) (m : ℕ) :
    gainedLoss p t (m + 2)
      = (t (m + 2) - (p - 1) * jumpSum t (m + 1) - t 1) + p * gainedLoss p t (m + 1) := by
  have hA := t1_le_A hadm m
  have hL := t1_le_gainedLoss hp ht1 hadm m
  have hpos : (0 : ℤ) ≤ t (m + 2) - (p - 1) * jumpSum t (m + 1) - t 1 := by linarith
  have hL0 : (0 : ℤ) ≤ gainedLoss p t (m + 1) := le_trans ht1 hL
  have h2L : 2 * gainedLoss p t (m + 1) ≤ p * gainedLoss p t (m + 1) :=
    mul_le_mul_of_nonneg_right hp hL0
  have hpL : t 1 ≤ p * gainedLoss p t (m + 1) := by linarith
  have hunfold : gainedLoss p t (m + 2)
      = max (t (m + 2) - (p - 1) * jumpSum t (m + 1))
        (max 0 (t (m + 2) - (p - 1) * jumpSum t (m + 1) - t 1)
          + p * gainedLoss p t (m + 1)) := rfl
  rw [hunfold, max_eq_right hpos, max_eq_right (by linarith)]

end Branch

/-! ## §3 許容列では閉じた形が**無条件に**成り立つ -/

section Closed

variable {p : ℤ} {t : ℕ → ℤ}

/-- ★★★★**`GainedTowerDescent.lean:94` の「証明していない」が証明になった。**

許容列（`p·t_j ≤ t_{j+1}`、`0 ≤ t 1`、`2 ≤ p`）では**無条件に**
`Λ_{k+1} = t_{k+1} + J_k − t_1·geomLow p k`。

★前波は「一般には偽」を示した。★本波は「木が語っていた母集団では真」を示した。
★★**「偽」と書くときは母集団を必ず添えること。** -/
theorem gainedLoss_closed_of_admissible (hp : 2 ≤ p) (ht1 : 0 ≤ t 1)
    (hadm : ∀ j, p * t (j + 1) ≤ t (j + 2)) (k : ℕ) :
    gainedLoss p t (k + 1) = t (k + 1) + jumpSum t k - t 1 * geomLow p k :=
  gainedLoss_closed (hbr_of_admissible hp ht1 hadm) ht1 k

end Closed

/-! ## §4 使っている公理の一覧 -/

#print axioms t1_le_A
#print axioms t1_le_gainedLoss
#print axioms hbr_of_admissible
#print axioms gainedLoss_closed_of_admissible

end AdmissibleClosedForm

end ABC3.Found.PGC
