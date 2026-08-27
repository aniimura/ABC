import ABC3.Skeleton.GenEll.Section1

/-!
# ★[GenEll] `Proposition 1.7` の仮定は空虚ではない(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★なぜこの検査が要るか

2026-08-27(第 423 ブロック)に `prop_1_7` を構成へ載せ替えたとき、
原文の証明の **prime-to-`Σ` の部分を仮定 `hlow` / `hup` として受けた**。

★**仮定を置いて `sorry` を消すのは、それ自体では進捗ではない**
——満たされない仮定を置けば statement は空虚になる。

★★本ファイルは仮定が**実際に満たされる場合がある**ことを機械検証する。
`Σ = ∅`、すべての量を `0` に取れば `hlow` / `hup` は成り立ち、
`prop_1_7` の結論(2 つの `BDle`)が実際に得られる。

★★★これは「elementary theory of differents を実装した」ではない。
本当に要るのは **`hlow` / `hup` を局所の different から供給する大域の組み立て**であり、
局所の核は `TameRamification.lean` / `DifferentKummer.lean` / `TotallyRamified.lean`
に取ってある(第 374–411、馴分岐 6/6)。
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GenEll

/-- ★**`prop_1_7` の仮定は満たされうる** —— 退化した場合(すべて `0`、`Σ = ∅`)。 -/
theorem prop_1_7_hypotheses_satisfiable :
    (∀ _x : ℕ, (0 : ℝ) - 0 ≤ ((0 : ℝ) - 0) + 0)
  ∧ (∀ _x : ℕ, (1 - 1 / ((1 : ℕ) : ℝ)) * (0 : ℝ) ≤ ((0 : ℝ) - 0) + 0)
  ∧ (∀ _x : ℕ, (0 : ℝ) ≤ ∑ q ∈ (∅ : Finset ℕ), Real.log q) := by
  refine ⟨fun _ => by norm_num, fun _ => by norm_num, fun _ => by simp⟩

/-- ★★**結論が実際に出る** —— 満たされる仮定で `prop_1_7` を叩くと 2 つの `BDle` が得られる。

★これで「仮定が偽だから何でも言える」ではないことが確かめられた。 -/
theorem prop_1_7_conclusion_obtained :
    BDle (fun _ : ℕ => (0 : ℝ) - 0) (fun _ : ℕ => (0 : ℝ) - 0)
  ∧ BDle (fun _ : ℕ => (0 : ℝ) - 0)
         (fun _ : ℕ => (1 - 1 / ((1 : ℕ) : ℝ)) * (0 : ℝ)) :=
  ABC3.Skeleton.GenEll.prop_1_7
    (fun _ : ℕ => (0 : ℝ)) (fun _ => 0) (fun _ => 0) (fun _ => 0) (fun _ => 0) (fun _ => 0)
    1 Nat.one_pos ∅
    (fun _ => by norm_num) (fun _ => by norm_num)
    (fun _ => by simp) (fun _ => by simp)

end ABC3.Check.GenEll
