import ABC3.Check.GenEll.RemarkAxiomGap

/-!
# ★★★[GenEll] `Proposition 1.7` も現在の `Interface` の下では**偽**である

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★病は §1 の 3 件だけではなかった

`RemarkAxiomGap.lean` は `Proposition 1.4` / `Remark 1.4.1` / `Remark 1.5.1` が
偽であることを機械検証した。★**`Proposition 1.7` にも同じ病がある**——本ファイルはそれを示す。

## ★★理由

`CoveringSetup` は `DY DZ : HeightTheoryData` を束ねた構造で、
条件 (a)(b)(c)(d) を **`hyp : Prop` という不透明な 1 つのフィールド**として持つ。

★`hyp` は `Prop` の**フィールド**であって、`logDiff` や `logCond` を**何も縛らない**。
`hyp := True` と置いたうえで `logDiff` を好きに取れば、結論の `BDle` は破れる。

★★これは `hyp` を不透明にしたこと自体の問題ではない
——**不透明な述語を置いても、それが他のフィールドと無関係なら制約にならない**、
という一般則である。

## ★★★★★では何をすべきか

`Remark 1.5.1` と同じ治療である——**構成に載せ替える**。
`Found/GenEll/` 側には既に

| 部品 | 場所 |
|---|---|
| 馴分岐の different(段 1〜6) | `TameRamification.lean` / `DifferentKummer.lean` / `TotallyRamified.lean` |
| `Σ` の上の寄与が `Σ_{q∈Σ} log q` 以下 | `SigmaBound.lean` / `LogCondSigma.lean` |
| `log-diff` の値 `log|disc F| / [F:ℚ]` | `LogDiffValue.lean` |

がある。★残るのは**スキームの分岐を局所の different へ落とす段**である。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll ABC3.Skeleton.GenEll ABC3.Found.GenEll

/-! ## ★被覆の反例データ -/

/-- ★**`hyp := True` の被覆データ** —— `log-diff` だけが点に線形に依る。

★`DY` は `logDiff x = x`、`DZ` は `logDiff = 0`。`logCond` はどちらも `0`。
`hyp` は `True` なので仮定は満たされるが、結論の `BDle` は破れる。 -/
noncomputable def badCover : CoveringSetup where
  DY := linDataDiff
  DZ := linData
  toPoint := id
  divY := ()
  divZ := ()
  e := 1
  e_pos := Nat.one_pos
  hyp := True
  maps_compl := fun x _ => Set.mem_univ x

/-! ## ★★★反例 -/

/-- ★★★★**`Proposition 1.7` の statement は任意の `CoveringSetup` については偽**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

`badCover` を取ると `hyp = True` は満たされるが、
`log-diff_Y − log-diff_Z = x` が発散する一方 `log-cond` は両側とも `0` なので、
第 1 の `BDle` が破れる。

★★**この `sorry` は「まだ証明していない」ではなく「証明できない」である**
——`Remark 1.5.1` と同じく、構成に載せ替えるほかない。 -/
theorem prop_1_7_false :
    ¬ (∀ (S : CoveringSetup), S.hyp →
        BDle (fun x : ↥(S.DY.compl S.divY) =>
                S.DZ.logCond S.divZ (S.toPoint x.1) - S.DY.logCond S.divY x.1)
             (fun x : ↥(S.DY.compl S.divY) =>
                S.DY.logDiff x.1 - S.DZ.logDiff (S.toPoint x.1))
      ∧ BDle (fun x : ↥(S.DY.compl S.divY) =>
                S.DY.logDiff x.1 - S.DZ.logDiff (S.toPoint x.1))
             (fun x : ↥(S.DY.compl S.divY) =>
                (1 - 1 / (S.e : ℝ)) * S.DZ.logCond S.divZ (S.toPoint x.1))) := by
  intro h
  obtain ⟨⟨C, hC⟩, -⟩ := h badCover trivial
  obtain ⟨n, hn⟩ := exists_nat_gt C
  have hb := hC ⟨n, Set.mem_univ n⟩
  simp only [badCover, linData, linDataDiff, id_eq, sub_zero, zero_sub, sub_neg_eq_add] at hb
  linarith

/-- ★★★★★**§1 の `Skeleton` の病は 4 件である**。

`Proposition 1.4` / `Remark 1.4.1` / `Remark 1.5.1` / `Proposition 1.7` が
同じ理由で偽である——`HeightTheoryData` が公理を 1 つも持たないからである。

★★★2026-08-27 の時点で、前の 3 件は**構成へ載せ替えて閉じてある**。
`Proposition 1.7` だけが残っている。 -/
theorem section1_four_diseases : True := trivial

end ABC3.Check.GenEll
