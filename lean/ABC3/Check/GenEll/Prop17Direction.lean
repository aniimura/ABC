/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17

/-!
# ★★★★★★★★★★[GenEll] `Proposition 1.7` の右の `≲` —— **印字どおりの向きでは偽**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★これは `Gap/GenEll/BDDirection.lean` が待っていた `falsifier` である

`Gap/GenEll/BDDirection.lean` は、`≲` の**印字された定義**（物理 p.5）

> if there exists a ["constant"] C ∈ R such that β(x) − α(x) ≤ C (respectively,

と、**同じ論文の中の用法**（`Proposition 1.6`・`Theorem 2.1`・`Proposition 1.4, (ii)`）が
互いに逆を向いていることを記録し、こう書いていた:

> ★**`falsifier` を書けないうちは ③ を名乗ってはならない**

★★★★★★**本ファイルがその `falsifier` である**——`Proposition 1.7, (i)` の右の `≲` を
印字どおりの向き（`BDle`）で読むと、**主張が偽になる**。

## ★★★反例の中身（幾何の実例つき）

`Proposition 1.7` の設定で **`φ = id`**（`Y = Z`、`K = F`）を取る:

| 条件 | 満たされるか |
|---|---|
| (b) `D_ℚ = φ_ℚ^{-1}(E_ℚ)_red` | ★`D ≔ E_red` と置けばよい |
| (c) `φ_ℚ` は `(U_Y)_ℚ → (U_Z)_ℚ` 上有限エタール | ★恒等射だから |
| (d) 分岐指数が `e` を割る | ★すべて `e_P = 1`。任意の `e` を割る |

このとき

* `log-diff(K) − log-diff(F) = 0`（同じ体だから）
* `log-cond_D = log-cond_E`（`log-cond` は根基で定義されている——`§9-966`）
* ★したがって**我々が証明した左の等式は成り立つ**: `log-cond_E − log-cond_D = 0`

★★★ところが右の `≲` を**印字どおり**（`BDle diff ((1−1/e)·cond_E)`、
すなわち `(1−1/e)·log-cond_E(x) − 0 ≤ C`）と読むと、
`e = 2` で `log-cond_E(x)/2 ≤ C` を**すべての点で**要求することになる。
★導手は点にわたって非有界だから、そのような `C` は無い。

★★★★★★**通常の読み**（`log-diff の差 ≤ (1−1/e)·log-cond_E + C`、
本プロジェクトの `BDge`）なら `0 ≤ (1−1/e)·log-cond_E` で**成り立つ**
——`prop_1_7_ordinary_direction_holds` がそれを確かめる。

## ★★分類について —— それでも ③ を名乗らない

`Gap/GenEll/BDDirection.lean` の規律に従い、本ファイルは **`Check`** に置く。
★示したのは「印字どおりの読みでは `Proposition 1.7, (i)` が偽である」ことであって、
「原典に穴がある」ことではない。★★`≲` の印字された定義が誤植である可能性と、
我々が `≲` の流儀を読み違えている可能性の**両方**が残っている。
★本ファイルは事実（機械検証された反例）だけを置く。
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GenEll

/-! ## ★★★★★★★★★★印字どおりの向きでは偽 -/

/-- ★★★★★★★★★★**`Proposition 1.7, (i)` の右の `≲` は、印字どおりの向き
（`BDle`）では偽である**。

★仮定として置いてあるのは、我々が実際に証明したこと（左の**等式**）と、
導手が非負であることだけである。★★それでも印字どおりの右の `≲` は破れる。

★反例: `φ = id` の場合（`K = F`、`D = E_red`）。
`log-diff` の差は `0`、`log-cond_D = log-cond_E`、そして導手は非有界。 -/
theorem prop_1_7_printed_direction_false :
    ¬ ∀ (Pt : Type) (condE condD diffY : Pt → ℝ) (e : ℕ), 0 < e →
        (∀ x, condE x - condD x = diffY x) →
        (∀ x, 0 ≤ condE x) →
        (∀ x, 0 ≤ condD x) →
        BDle diffY (fun x => (1 - 1 / (e : ℝ)) * condE x) := by
  intro h
  obtain ⟨C, hC⟩ := h ℕ (fun n => (n : ℝ)) (fun n => (n : ℝ)) (fun _ => 0) 2 (by norm_num)
    (fun n => by ring) (fun n => by positivity) (fun n => by positivity)
  obtain ⟨n, hn⟩ := exists_nat_gt (2 * C)
  have hb := hC n
  norm_num at hb
  linarith

/-! ## ★★★★★通常の読みなら成り立つ -/

/-- ★★★★★**同じ設定で、通常の読み（`BDge`）なら成り立つ**。

    `log-diff の差 ≤ (1 − 1/e)·log-cond_E + C`

★これが `Found/GenEll/Prop17.lean` の `prop_1_7` が取っている向きである。
★★上の `prop_1_7_printed_direction_false` と合わせて、
**2 つの読みが同じ主張ではないこと**が機械で確かめられた。 -/
theorem prop_1_7_ordinary_direction_holds :
    ∀ (Pt : Type) (condE condD diffY : Pt → ℝ) (e : ℕ), 0 < e →
        (∀ x, condE x - condD x = diffY x) →
        (∀ x, 0 ≤ condE x) →
        (∀ x, (1 / (e : ℝ)) * condE x ≤ condD x) →
        BDge diffY (fun x => (1 - 1 / (e : ℝ)) * condE x) := by
  intro Pt condE condD diffY e _ heq _ hlow
  refine ⟨0, fun x => ?_⟩
  have h1 := heq x
  have h2 := hlow x
  nlinarith

/-- ★★★★**`Gap/GenEll/BDDirection.lean` の用例表に 4 例目が加わった**。

| 箇所 | 印字 | 印字どおりに読むと |
|---|---|---|
| `Proposition 1.6` | `log-cond_D ≲ ht_D` | 表題 "Conductor **Bounded by** the Height" と逆 |
| `Theorem 2.1`, (i)(ii) | `ht ≲ (1+ϵ)(…)` | abc 予想と逆 |
| `Proposition 1.4`, (ii) | `ht_L̄ ≳ 0` | 「下に有界」と逆 |
| ★`Proposition 1.7`, (i) 右 | `log-diff の差 ≲ (1−1/e)·log-cond_E` | ★★**偽になる**(本ファイル) |

★★★前の 3 例は「表題・文脈と逆」であった。★★★★★**本例は「偽になる」**
——これが `Gap` が待っていた `falsifier` である。 -/
theorem bd_direction_fourth_case : True := trivial

end ABC3.Check.GenEll
