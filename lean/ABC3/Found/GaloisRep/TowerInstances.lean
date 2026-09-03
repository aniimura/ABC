/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35Descend
import ABC3.Meta.Claim

/-!
# 第 1222 ブロック —— **中間体への塔のインスタンス**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——組み上げの最後の障害

第 1221（`lemma_3_5_height_ineq_descend`）は `L''` について

    [Algebra (𝓞 L) L''] [IsScalarTower (𝓞 L) L L''] [IsScalarTower (𝓞 L) (𝓞 L'') L'']
    [Module.Finite (𝓞 L) (𝓞 L'')] [Algebra.IsIntegral (𝓞 L) (𝓞 L'')]

を要求する（`DegInfBaseChange` 以来の慣行）。

☆**測ったこと（2026-09-02）**——`M : IntermediateField L L̄` が `L` 上有限次のとき:

| インスタンス | 出るか |
|---|---|
| `NumberField M` | ★mathlib の `NumberField.of_module_finite` |
| `IsScalarTower ℚ L M` | ★自動 |
| `Algebra (𝓞 L) M` | ★自動 |
| `Module.Finite (𝓞 L) (𝓞 M)` | ★自動 |
| `Algebra.IsIntegral (𝓞 L) (𝓞 M)` | ★自動 |
| `IsScalarTower (𝓞 L) L M` | ☆**自動では出ない**——`of_algebraMap_eq fun _ => rfl` |
| `IsScalarTower (𝓞 L) (𝓞 M) M` | ☆**自動では出ない**——`of_algebraMap_eq fun _ => rfl` |

★★★どちらも `rfl` で出る——`Algebra (𝓞 L) M` は
`𝓞 L ⊆ L → M` の制限だからである。
☆自動で出ないのは mathlib のインスタンス探索が
`Algebra (𝓞 L) M` を別経路で解決してしまうためで、**数学の穴ではない**。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain

/-- ☆`IsScalarTower (𝓞 L) L M`——`rfl` で出る（第 1222）。 -/
theorem isScalarTower_ringOfIntegers_base (L : Type) [Field L] [NumberField L]
    (M : IntermediateField L (AlgebraicClosure L)) :
    IsScalarTower (𝓞 L) L M :=
  IsScalarTower.of_algebraMap_eq fun _ => rfl

/-- ☆`IsScalarTower (𝓞 L) (𝓞 M) M`——`rfl` で出る（第 1222）。 -/
theorem isScalarTower_ringOfIntegers_top (L : Type) [Field L] [NumberField L]
    (M : IntermediateField L (AlgebraicClosure L)) :
    IsScalarTower (𝓞 L) (𝓞 M) M :=
  IsScalarTower.of_algebraMap_eq fun _ => rfl

/-- ★★★★★★★★★★★★
**第 1221 が要るインスタンスは中間体で全部揃う**——★**無条件**（第 1222）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`M` が `L` 上有限次なら `NumberField M` と 2 つの塔がそろい、
残りは mathlib が自動で出す。

★★★これで第 1219（`[M:L] ≤ l−1` の `M` と点 `Q'`）を
第 1221 に**そのまま入れられる**。 -/
theorem numberField_and_towers (L : Type) [Field L] [NumberField L]
    (M : IntermediateField L (AlgebraicClosure L)) [FiniteDimensional L M] :
    NumberField M ∧ IsScalarTower (𝓞 L) L M ∧ IsScalarTower (𝓞 L) (𝓞 M) M :=
  ⟨NumberField.of_module_finite L M,
   isScalarTower_ringOfIntegers_base L M,
   isScalarTower_ringOfIntegers_top L M⟩

/-! ## ★出典の紐付け(`.src`) -/

def isScalarTower_ringOfIntegers_base.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsScalarTower (𝓞 L) L M は rfl で出る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isScalarTower_ringOfIntegers_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsScalarTower (𝓞 L) (𝓞 M) M は rfl で出る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def numberField_and_towers.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(第 1221 が要るインスタンスは中間体で全部揃う。★無条件)",
    sectionId := "genell-lemma-3-5" }

def numberField_and_towers.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1222）**——第 1219 の `M` を第 1221 に入れる" ++
       "**最後の障害**であった。☆`IsScalarTower (𝓞 L) L M` と" ++
       "`IsScalarTower (𝓞 L) (𝓞 M) M` は自動では出ないが、どちらも `rfl` で出る" ++
       "——`Algebra (𝓞 L) M` が `𝓞 L ⊆ L → M` の制限だからである。" ++
       "★**数学の穴ではなくインスタンス探索の経路の問題**である。") 2 ]

end ABC3.Found.GaloisRep
