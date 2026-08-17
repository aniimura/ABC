import ABC3.Skeleton.GenEll.Section1

/-!
# ★★★[GenEll] `Proposition 1.4` は現在の `Interface` の下では**偽**である

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★何が起きているか

`Skeleton/GenEll/Section1.lean` の `prop_1_4` は

> `∀ D : HeightTheoryData, (ht の加法性) ∧ … ∧ (Northcott)`

という形をしている。★しかし **`HeightTheoryData` は公理を 1 つも持たないデータ**である
(`Interface/GenEll/HeightTheory.lean` の docstring が明言している——
「本ファイルは**公理を 1 つも持たない**(データと述語だけ)」)。

★★**したがって `prop_1_4` は「任意のデータについて」を主張しており、
原文より強く、しかも偽である。**

★これは `sorry` が「まだ証明していない」のではなく
**「証明できない」**ことを意味する。★★**その区別を機械検証された事実にする**のが
本ファイルである。

## ★原文は偽ではない

原文の `Proposition 1.4` は**本物の算術直線束**についての主張であり、
そこでは加法性も Northcott も成り立つ。
★★**偽になったのは、我々が `Interface` で「語彙だけ」を posit し、
`Skeleton` でそれを全称量化したからである。**

★`tools/check.mjs` 冒頭 B5 は「条件を posit して `sorry` を消しても進捗ではない」
と警告する。★★**本件はその裏側**——**条件を posit しないまま全称量化すると
statement が偽になる**。B5 と対になる穴である。

## ★★これが計画に対して持つ意味

`§1` の 9 項目を閉じるには、次のどちらかが要る:

1. `HeightTheoryData` に**公理を足す**(加法性・Northcott・…)。
   ★しかしそれは B5 そのもの——`Proposition 1.4` が仮定の言い換えになる。
2. ★★**`HeightTheoryData` を posit ではなく構成に置き換える**
   ——本物の算術直線束を作る。

★**2 しか道がない。** これが「Arakelov を作る」ことの中身である。

## ★反例(下で構成する)

`Point ≝ ℕ`、`ABundle ≝ ℝ`、`tensor ≝ (+)`、`ht L x ≝ L²`、
`degLe d ≝ univ`、`Ample ≝ True` と置く。すると

- **(i) 加法性が破れる**: `ht(1⊗1) = 4` だが `ht(1) + ht(1) = 2`
- **(iv) Northcott が破れる**: `{x ∈ univ | ht 0 x ≤ 0} = ℕ` は無限

★★1 つのデータで **4 条のうち 2 条**が破れる。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll

/-- ★**反例のデータ**。

`Point ≝ ℕ`(無限)、`ABundle ≝ ℝ`、`tensor ≝ (+)`、`ht L x ≝ L²`。
他のフィールドは最も単純なもので埋める(反例には効かない)。 -/
noncomputable def badData : HeightTheoryData where
  Point := ℕ
  ABundle := ℝ
  Divisor := Unit
  GenericClass := Unit
  tensor := fun L M => L + M
  ht := fun L _ => L ^ 2
  generic := fun _ => ()
  Ample := fun _ => True
  SomePowerGlobGen := fun _ => True
  degLe := fun _ => Set.univ
  bundleOf := fun _ => 0
  compl := fun _ => Set.univ
  logDiff := fun _ => 0
  logCond := fun _ _ => 0
  CompactlyBounded := fun _ => True

/-- ★★**`Proposition 1.4, (i)`(高さの加法性)は任意のデータについては偽**。

`ht(1⊗1) = (1+1)² = 4` だが `ht(1) + ht(1) = 1 + 1 = 2`。 -/
theorem prop_1_4_i_false :
    ¬ (∀ D : HeightTheoryData, ∀ (L M : D.ABundle) (x : D.Point),
        D.ht (D.tensor L M) x = D.ht L x + D.ht M x) := by
  intro h
  have h1 : ((1 : ℝ) + 1) ^ 2 = (1 : ℝ) ^ 2 + (1 : ℝ) ^ 2 :=
    h badData (1 : ℝ) (1 : ℝ) (0 : ℕ)
  norm_num at h1

/-- ★★**`Proposition 1.4, (iv)`(Northcott)は任意のデータについては偽**。

`degLe d = univ`(= `ℕ`)で `ht 0 x = 0 ≤ 0` なので、集合は `ℕ` 全体になる。 -/
theorem prop_1_4_iv_false :
    ¬ (∀ D : HeightTheoryData, ∀ (L : D.ABundle) (d : ℕ) (C : ℝ), 0 < d →
        D.Ample (D.generic L) → {x ∈ D.degLe d | D.ht L x ≤ C}.Finite) := by
  intro h
  have hfin : {x : ℕ | x ∈ (Set.univ : Set ℕ) ∧ ((0 : ℝ)) ^ 2 ≤ 0}.Finite :=
    h badData (0 : ℝ) 1 0 Nat.one_pos trivial
  have huniv : {x : ℕ | x ∈ (Set.univ : Set ℕ) ∧ ((0 : ℝ)) ^ 2 ≤ 0}
      = (Set.univ : Set ℕ) := by
    ext x
    simp
  rw [huniv] at hfin
  exact Set.infinite_univ hfin

/-- ★★★**`prop_1_4` の statement 全体が偽**。

★したがって `Skeleton/GenEll/Section1.lean` の `prop_1_4` の `sorry` は
「まだ証明していない」ではなく **「証明できない」**である。 -/
theorem prop_1_4_statement_false :
    ¬ (∀ D : HeightTheoryData,
        (∀ (L M : D.ABundle) (x : D.Point),
            D.ht (D.tensor L M) x = D.ht L x + D.ht M x)
      ∧ (∀ L : D.ABundle, D.SomePowerGlobGen (D.generic L) →
            ABC3.Found.GenEll.BDge (D.ht L) (fun _ => (0 : ℝ)))
      ∧ (∀ L L' : D.ABundle, D.generic L = D.generic L' →
            ABC3.Found.GenEll.BDeq (D.ht L) (D.ht L'))
      ∧ (∀ (L : D.ABundle) (d : ℕ) (C : ℝ), 0 < d → D.Ample (D.generic L) →
            {x ∈ D.degLe d | D.ht L x ≤ C}.Finite)) := by
  intro h
  exact prop_1_4_i_false (fun D => (h D).1)

end ABC3.Check.GenEll
