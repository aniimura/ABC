import ABC3.Skeleton.GenEll.Section1

/-!
# ★★★[GenEll] `Remark 1.4.1` / `Remark 1.5.1` も現在の `Interface` の下では**偽**である

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★病は `Proposition 1.4` だけではなかった

`Check/GenEll/HeightAxiomGap.lean` で `prop_1_4` が偽であることを機械検証した。
★★**同じ病が `Remark 1.4.1` と `Remark 1.5.1` にもある**——本ファイルはそれを示す。

どちらも `∀ D D' : HeightTheoryData, …` の形をしており、
`HeightTheoryData` は**公理を 1 つも持たないデータ**なので、
仮定(点の同一視・生成ファイバーの対応)は `ht` や `logDiff` を**何も縛らない**。

## ★★反例の作り方

★`HeightAxiomGap.lean` の `badData` は `ht L x = L²` で **`x` に依らない**ので、
`BDeq`(定数差を許す)は**成り立ってしまう**。反例にならない。

★★**`x` に線形に依る `ht` を作る**必要がある——それが `linData` である。
`ht L x ≝ L · x` とすると、`L = 0` と `L' = 1` は
生成ファイバー(`Unit`)では区別できないのに、高さの差が `x` で発散する。

★★`Remark 1.5.1` も同様で、`logDiff` を `0` と `x` に取れば
「点の同一視がある」だけでは `BDeq` は出ない。

## ★★★これが計画に対して持つ意味

★§1 の `Skeleton` は**系統的に**証明できない。
`Proposition 1.4` 1 件の問題ではなく、
**`HeightTheoryData` を posit したことの帰結**である。

★★★**閉じる道は 1 つ**——`Found/GenEll/HeightConstruction.lean` のように
**posit を構成に置き換える**ことである。
構成の側では `Proposition 1.4, (ii)(iii)` は**無条件で真**になっている
(`HeightNonneg.lean` / `HeightClass.lean`)。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll ABC3.Skeleton.GenEll ABC3.Found.GenEll

/-! ## ★`x` に線形に依る高さを持つデータ -/

/-- ★**高さが点に線形に依るデータ**。

★`badData`(`ht L x = L²`)は `x` に依らないので `BDeq` の反例にならない。
★★**点に依らせる**ことで初めて BD-class の主張を破れる。 -/
noncomputable def linData : HeightTheoryData where
  Point := ℕ
  ABundle := ℝ
  Divisor := Unit
  GenericClass := Unit
  tensor := fun L M => L + M
  ht := fun L x => L * (x : ℝ)
  generic := fun _ => ()
  Ample := fun _ => True
  SomePowerGlobGen := fun _ => True
  degLe := fun _ => Set.univ
  bundleOf := fun _ => 0
  compl := fun _ => Set.univ
  logDiff := fun _ => 0
  logCond := fun _ _ => 0
  CompactlyBounded := fun _ => True

/-- ★**`logDiff` だけを取り替えたデータ** —— `Remark 1.5.1` を破るために使う。 -/
noncomputable def linDataDiff : HeightTheoryData where
  Point := ℕ
  ABundle := ℝ
  Divisor := Unit
  GenericClass := Unit
  tensor := fun L M => L + M
  ht := fun L x => L * (x : ℝ)
  generic := fun _ => ()
  Ample := fun _ => True
  SomePowerGlobGen := fun _ => True
  degLe := fun _ => Set.univ
  bundleOf := fun _ => 0
  compl := fun _ => Set.univ
  logDiff := fun x => (x : ℝ)
  logCond := fun _ _ => 0
  CompactlyBounded := fun _ => True

/-! ## ★★★反例 -/

/-- ★★★**`Remark 1.4.1` の statement は任意のデータについては偽**。

`linData` で `L = 0`、`L' = 1` を取る。生成ファイバーは `Unit` なので
仮定 `eGen (generic L) = generic L'` は成り立つが、
高さの差は `|0 · x − 1 · x| = x` で発散する。 -/
theorem remark_1_4_1_false :
    ¬ (∀ (D D' : HeightTheoryData) (ePt : D.Point ≃ D'.Point)
        (eGen : D.GenericClass ≃ D'.GenericClass) (L : D.ABundle) (L' : D'.ABundle),
        eGen (D.generic L) = D'.generic L' →
        BDeq (D.ht L) (fun x => D'.ht L' (ePt x))) := by
  intro h
  obtain ⟨C, hC⟩ :=
    h linData linData (Equiv.refl ℕ) (Equiv.refl Unit) (0 : ℝ) (1 : ℝ) rfl
  obtain ⟨n, hn⟩ := exists_nat_gt C
  have hb : |(0 : ℝ) * (n : ℝ) - (1 : ℝ) * (n : ℝ)| ≤ C := hC n
  rw [zero_mul, one_mul, zero_sub, abs_neg, Nat.abs_cast] at hb
  linarith

/-- ★★★**`Remark 1.5.1` の statement は任意のデータについては偽**。

`linData`(`logDiff = 0`)と `linDataDiff`(`logDiff x = x`)を取る。
点の集合はどちらも `ℕ`、`compl` はどちらも全体なので仮定は成り立つが、
`log-diff` の差は `x` で発散する。 -/
theorem remark_1_5_1_false :
    ¬ (∀ (D D' : HeightTheoryData) (ePt : D.Point ≃ D'.Point)
        (dv : D.Divisor) (dv' : D'.Divisor),
        (∀ x, x ∈ D.compl dv ↔ ePt x ∈ D'.compl dv') →
        BDeq (D.logDiff) (fun x => D'.logDiff (ePt x))
      ∧ BDeq (fun x : ↥(D.compl dv) => D.logCond dv x.1)
             (fun x : ↥(D.compl dv) => D'.logCond dv' (ePt x.1))) := by
  intro h
  obtain ⟨⟨C, hC⟩, -⟩ :=
    h linData linDataDiff (Equiv.refl ℕ) () () (fun _ => Iff.rfl)
  obtain ⟨n, hn⟩ := exists_nat_gt C
  have hb : |(0 : ℝ) - (n : ℝ)| ≤ C := hC n
  rw [zero_sub, abs_neg, Nat.abs_cast] at hb
  linarith

/-- ★★★**§1 の `Skeleton` は系統的に証明できない**。

`Proposition 1.4` / `Remark 1.4.1` / `Remark 1.5.1` の 3 件が同じ理由で偽である。
★これらの `sorry` は「まだ証明していない」ではなく **「証明できない」**である。 -/
theorem section1_skeleton_systematically_false :
    ¬ ((∀ (D D' : HeightTheoryData) (ePt : D.Point ≃ D'.Point)
          (eGen : D.GenericClass ≃ D'.GenericClass) (L : D.ABundle) (L' : D'.ABundle),
          eGen (D.generic L) = D'.generic L' →
          BDeq (D.ht L) (fun x => D'.ht L' (ePt x)))
      ∧ (∀ (D D' : HeightTheoryData) (ePt : D.Point ≃ D'.Point)
          (dv : D.Divisor) (dv' : D'.Divisor),
          (∀ x, x ∈ D.compl dv ↔ ePt x ∈ D'.compl dv') →
          BDeq (D.logDiff) (fun x => D'.logDiff (ePt x))
        ∧ BDeq (fun x : ↥(D.compl dv) => D.logCond dv x.1)
               (fun x : ↥(D.compl dv) => D'.logCond dv' (ePt x.1)))) := by
  intro h
  exact remark_1_4_1_false h.1

end ABC3.Check.GenEll
