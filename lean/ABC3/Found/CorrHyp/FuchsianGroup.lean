/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.Complex.UpperHalfPlane.Basic
import Mathlib.Analysis.Complex.UpperHalfPlane.MoebiusAction
import Mathlib.Analysis.Complex.UpperHalfPlane.ProperAction
import Mathlib.NumberTheory.ModularForms.ProperlyDiscontinuous
import Mathlib.GroupTheory.Schreier

/-!
# [CorrHyp] Track B の第一歩 —— Fuchsian 群を mathlib 上で実装する

`ResearchPaper/corrhyp-goal.md` §3(Track B)の入口。`Interface/CorrHyp/HyperbolicCurve.lean`
が posit した `Fuchsian` を、実際の `SL(2,ℝ)` の離散部分群として実装する
(まだ `HyperbolicCurveData` の具体的インスタンスは組んでいない——これは
その最初の 1 ブロックである)。

## ★★逸脱(記録)—— `PSL₂(ℝ)⁰` ではなく `SL(2,ℝ)` の部分群として扱う

原文は `Γ ⊆ PSL₂(ℝ)⁰`(射影特殊線形群の単位元成分)と書く。mathlib の
`ℍ` 上の Möbius 作用(`UpperHalfPlane.instProperlyDiscontinuousSL2RSubgroup` 等)は
`SL(2,ℝ)` を経由して定義されており、中心 `{±I}` は `ℍ` 上自明に作用する。
`SL(2,ℝ)` の部分群として扱っても、`ℍ` 上で見える作用・固有不連続性・
(離散な)商の構造は同じである。★`PSL₂(ℝ)⁰`(単位元成分、向きを保つ側)への
限定も本ファイルでは行っていない——`SL(2,ℝ)` 自体は連結なので、単位元成分の
制約は自動的に満たされる(`PSL₂(ℝ)` 全体には向きを反転する側が無い)。 -/

namespace ABC3.Found.CorrHyp

open UpperHalfPlane Matrix

/-- Fuchsian 群: `SL(2,ℝ)` の離散部分群。 -/
structure FuchsianGroup where
  /-- 部分群本体。 -/
  Γ : Subgroup (SpecialLinearGroup (Fin 2) ℝ)
  /-- 離散性(原文 §2: 「Γ ⊆ PSL₂(ℝ)⁰ の部分群として扱う」に暗黙に伴う仮定——
  基本群の deck 変換群の像は必ず離散である)。 -/
  discrete : DiscreteTopology Γ

attribute [instance] FuchsianGroup.discrete

/-- `F.Γ` は `ℍ` に固有不連続に作用する。

★**sorry 無し**——mathlib の `UpperHalfPlane.instProperlyDiscontinuousSL2RSubgroup`
から無条件に従う。原文 §2 冒頭の「Π ↪ Aut(H) = PSL₂(ℝ)⁰」の deck 変換の
議論が要求する性質そのもの。 -/
instance properlyDiscontinuous (F : FuchsianGroup) :
    ProperlyDiscontinuousSMul F.Γ ℍ := inferInstance

/-- `Γ₁` が `Γ₂` の部分群(`Interface` の `Sub`、`Proposition 3.2` で使う)。 -/
def IsSub (F1 F2 : FuchsianGroup) : Prop := F1.Γ ≤ F2.Γ

/-- `Γ₁` が `Γ₂` の中で有限指数(`Interface` の `FiniteIndexIn`)。 -/
def IsFiniteIndexIn (F1 F2 : FuchsianGroup) : Prop :=
  IsSub F1 F2 ∧ (F1.Γ.subgroupOf F2.Γ).FiniteIndex

/-- `Group.FG` は群同型で運べる(mathlib に無かったので補った、20 行)。 -/
private theorem fg_transport {G H : Type*} [Group G] [Group H] (e : G ≃* H) [Group.FG G] :
    Group.FG H := by
  obtain ⟨S, hS, hSfin⟩ := Group.fg_iff.mp ‹Group.FG G›
  refine Group.fg_iff.mpr ⟨e.toMonoidHom '' S, ?_, hSfin.image _⟩
  rw [← MonoidHom.map_closure, hS, Subgroup.map_top_of_surjective _ e.surjective]

/-- `Γ` が有限生成なら、その中の有限指数部分群も有限生成
(Schreier の補題、mathlib `Subgroup.fg_of_index_ne_zero`)。

★**sorry 無し**。原文が非明示的に使っている「Γ_X は有限生成(位相的基本群)」
という事実の系——`Theorem 3.3`/`Lemma 5.5` の証明の核となる Schreier の議論の
最初のブロック。 -/
theorem fg_of_finiteIndexIn {F1 F2 : FuchsianGroup} (hfg : Group.FG F2.Γ)
    (h : IsFiniteIndexIn F1 F2) : Group.FG F1.Γ := by
  haveI : Group.FG F2.Γ := hfg
  haveI : (F1.Γ.subgroupOf F2.Γ).FiniteIndex := h.2
  haveI : Group.FG (F1.Γ.subgroupOf F2.Γ) := Subgroup.fg_of_index_ne_zero _
  exact fg_transport (Subgroup.subgroupOfEquivOfLe h.1)

end ABC3.Found.CorrHyp
