import ABC3.Skeleton.GenEll.Section2

/-!
# ★[GenEll] `Theorem 2.1` の載せ替えが取りこぼしたもの(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.11。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★このファイルは「まだ取れていない」ことの記録である

2026-08-27(第 426 ブロック)に `theorem_2_1` を構成へ載せ替え、`sorry` を外した。
★**しかし取れたのは (i) ⟹ (ii) だけである**。

★★原文の `Theorem 2.1` は**同値**であり、**実質は (ii) ⟹ (i) の側**
(原文 p.11–p.13 の 3 ページ、noncritical Belyi maps と `Proposition 1.7` を使う)である。

★★★**逆向きが取れていないことを、型で残しておく。**
`theorem_2_1_converse_pending` は「まだ書けない」ことの記録であって、主張ではない。

## ★★結論が空虚でないこと

`theorem_2_1_nonvacuous` で、実際に (i) が成り立つ場合に (ii) が出ることを確かめる。
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GenEll

/-- ★**結論は空虚ではない** —— 全体で成り立つ場合に部分集合でも出る。 -/
theorem theorem_2_1_nonvacuous :
    BDle (fun x : ↥((Set.univ : Set ℕ) ∩ (Set.univ : Set ℕ)) => (1 + (0:ℝ)) * ((0:ℝ) + 0))
         (fun x : ↥((Set.univ : Set ℕ) ∩ (Set.univ : Set ℕ)) => (0 : ℝ)) :=
  ABC3.Skeleton.GenEll.theorem_2_1 (fun _ : ℕ => (0 : ℝ)) (fun _ => 0) (fun _ => 0) 0
    Set.univ Set.univ ⟨0, fun _ => by norm_num⟩

/-- ★★★★**逆向き((ii) ⟹ (i))は取れていない** —— その事実を型で残す。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

★★★★★**これは主張ではない。**「何が残っているか」を機械が読める形にしただけである。
残っているのは:

| 要るもの | 状態 |
|---|---|
| noncritical Belyi maps(`[NCBelyi] Theorem 2.5`) | 一般の曲線への帰着(Riemann–Roch)が未了 |
| `Proposition 1.7` | 局所から大域への組み立てが未了 |
| 双曲曲線の étale 基本群の構造 | `[Stacks] 58.6` に原典、大きさ未測定 |
| compactly bounded subset の support の理論 | `Example 1.3, (ii)`、posit のまま |

★`WaitingFor` として書くのが本来だが、`Interface` の外なので docstring で記録する。 -/
def theorem_2_1_converse_pending : String :=
  "(ii) ⟹ (i) は未了。noncritical Belyi maps・Proposition 1.7・étale 基本群の構造が要る"

end ABC3.Check.GenEll
