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

/-! ## `FEt` の圏構造 —— `Definition 1.1`-`1.4` の具体化に要る 4 点

`Space := FuchsianGroup`・`FEt A B := IsFiniteIndexIn A B` というモデルで
`Corr`(`Definition 1.1`)を実装するには、`comp`(合成)と `pullback`(ファイバー積、
`Definition 1.4`)が要る。ここで両方を実際に構成する——`Skeleton/CorrHyp/Section1.lean`
の `Corr`/`comp'` が要求する 4 点(`comp`・`pullback`・`pbFst`・`pbSnd`)そのもの。 -/

/-- `Γ₁ ⊓ Γ₂` の離散性(離散群の部分群は離散)。`pullback` の台。 -/
def inter (F1 F2 : FuchsianGroup) : FuchsianGroup where
  Γ := F1.Γ ⊓ F2.Γ
  discrete := by
    haveI := F1.discrete
    haveI : DiscreteTopology (↥((F1.Γ ⊓ F2.Γ).subgroupOf F1.Γ)) := inferInstance
    exact (Subgroup.subgroupOfContinuousMulEquivOfLe
      (inf_le_left : F1.Γ ⊓ F2.Γ ≤ F1.Γ)).toHomeomorph.discreteTopology

/-- `FEt` の合成(`Interface.CorrHyp.HyperbolicCurveData.comp`)。

★**sorry 無し**——有限相対指数の乗法性(`Subgroup.relIndex_mul_relIndex`)から。 -/
theorem isFiniteIndexIn_comp {A B C : FuchsianGroup} (hAB : IsFiniteIndexIn A B)
    (hBC : IsFiniteIndexIn B C) : IsFiniteIndexIn A C := by
  refine ⟨hAB.1.trans hBC.1, ?_⟩
  haveI hAB' : A.Γ.IsFiniteRelIndex B.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hAB.2
  haveI hBC' : B.Γ.IsFiniteRelIndex C.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hBC.2
  rw [Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero] at hAB' hBC'
  rw [← Subgroup.isFiniteRelIndex_iff_finiteIndex, Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero,
      ← Subgroup.relIndex_mul_relIndex A.Γ B.Γ C.Γ hAB.1 hBC.1]
  exact mul_ne_zero hAB' hBC'

/-- `pullback f g` から `A` への射影(`Interface.CorrHyp.HyperbolicCurveData.pbFst`)。

★**sorry 無し**——`(A ⊓ B).relIndex C = (A⊓B).relIndex A · A.relIndex C` かつ
左辺が有限(`A`・`B` とも `C` の中で有限指数だから、`relIndex_inf_ne_zero`)なので、
積が有限なら両因子とも有限。 -/
theorem isFiniteIndexIn_pbFst {A B C : FuchsianGroup} (hA : IsFiniteIndexIn A C)
    (hB : IsFiniteIndexIn B C) : IsFiniteIndexIn (inter A B) A := by
  refine ⟨inf_le_left, ?_⟩
  haveI hA' : A.Γ.IsFiniteRelIndex C.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hA.2
  haveI hB' : B.Γ.IsFiniteRelIndex C.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hB.2
  rw [Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero] at hA' hB'
  have hinf : (A.Γ ⊓ B.Γ).relIndex C.Γ ≠ 0 := Subgroup.relIndex_inf_ne_zero hA' hB'
  rw [← Subgroup.isFiniteRelIndex_iff_finiteIndex, Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero]
  show (A.Γ ⊓ B.Γ).relIndex A.Γ ≠ 0
  intro h0
  have heq : (A.Γ ⊓ B.Γ).relIndex A.Γ * A.Γ.relIndex C.Γ = (A.Γ ⊓ B.Γ).relIndex C.Γ :=
    Subgroup.relIndex_mul_relIndex (A.Γ ⊓ B.Γ) A.Γ C.Γ inf_le_left hA.1
  rw [h0, zero_mul] at heq
  exact hinf heq.symm

/-- `pullback f g` から `B` への射影(`Interface.CorrHyp.HyperbolicCurveData.pbSnd`)。
`inter A B = inter B A`(`⊓` の可換性)へ帰着して `pbFst` を使い回す。 -/
theorem isFiniteIndexIn_pbSnd {A B C : FuchsianGroup} (hA : IsFiniteIndexIn A C)
    (hB : IsFiniteIndexIn B C) : IsFiniteIndexIn (inter A B) B := by
  have : inter A B = inter B A := by
    unfold inter; congr 1; exact inf_comm _ _
  rw [this]
  exact isFiniteIndexIn_pbFst hB hA

end ABC3.Found.CorrHyp
