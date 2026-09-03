/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.Instance
import Mathlib.NumberTheory.ModularForms.CongruenceSubgroups
import Mathlib.NumberTheory.ModularForms.ArithmeticSubgroups

/-!
# [CorrHyp] §1 の非空虚性——モジュラー群 `SL(2,ℤ)` と合同部分群 `Γ(2)`

`corrHypInstance`(`Instance.lean`)の `Space := FuchsianGroup` は抽象的には
非空虚だが、それだけでは「原文の意味での双曲曲線として実際に出現するデータ」
であることを示さない(`corrhyp-goal.md` の `StackType.waiting` と同型の懸念)。
ここでは**教科書的に意味のある**具体例——モジュラー群 `SL(2,ℤ)`(`SL(2,ℝ)`
への像として離散、`Matrix.SpecialLinearGroup.discreteSpecialLinearGroupIntRangeSL`)
と、その主合同部分群 `Γ(2)`(`CongruenceSubgroup.Gamma 2`、mathlib に
`FiniteIndex` インスタンス済み)——を使って、`Definition 1.5`
(`IsIsogenous`)が **`corrHypInstance` において真に非空虚**であることを示す:
`SL(2,ℤ)` と `Γ(2)` は**相異なる** `Space`(`FG2_ne`)でありながら
**isogenous**(`isIsogenous_witness`)。

★逐語ではなく非空虚性の witness(G2 ゲート相当)なので `.src` は個々の
補題には付けない——`corrHypInstance` 自体の非空虚性を補強する位置づけ。 -/

namespace ABC3.Found.CorrHyp

open Matrix

/-- `SL(2,ℤ) → SL(2,ℝ)`(各成分を整数から実数へ)。 -/
noncomputable def φ₂ : Matrix.SpecialLinearGroup (Fin 2) ℤ →* Matrix.SpecialLinearGroup (Fin 2) ℝ :=
  Matrix.SpecialLinearGroup.map (Int.castRingHom ℝ)

theorem φ₂_injective : Function.Injective φ₂ := Matrix.SpecialLinearGroup.map_intCast_injective

/-- `SL(2,ℤ)` の `SL(2,ℝ)` における像(モジュラー群そのもの)。 -/
noncomputable def Γ_SL2Z : Subgroup (Matrix.SpecialLinearGroup (Fin 2) ℝ) := φ₂.range

instance : DiscreteTopology Γ_SL2Z :=
  Matrix.SpecialLinearGroup.discreteSpecialLinearGroupIntRangeSL

/-- `Γ(2)`(mod 2 で単位行列に合同、`SL(2,ℤ)` の主合同部分群)の `SL(2,ℝ)` への像。 -/
noncomputable def Γ_Gamma2 : Subgroup (Matrix.SpecialLinearGroup (Fin 2) ℝ) :=
  (CongruenceSubgroup.Gamma 2).map φ₂

theorem Γ_Gamma2_le : Γ_Gamma2 ≤ Γ_SL2Z := Subgroup.map_le_range φ₂ (CongruenceSubgroup.Gamma 2)

instance : DiscreteTopology Γ_Gamma2 :=
  (Subgroup.subgroupOfContinuousMulEquivOfLe Γ_Gamma2_le).toHomeomorph.discreteTopology

/-- `Γ(2)` は `SL(2,ℤ)` の中で有限指数(mathlib の `Gamma 2` の
`FiniteIndex` インスタンスを、単射な `φ₂` を通して像へ運んだもの——
`Subgroup.relIndex_map_map` で `ker = ⊥` を消し込むだけ)。

★**sorry 無し**。標準3公理のみ。 -/
theorem Γ_Gamma2_finiteIndex : (Γ_Gamma2.subgroupOf Γ_SL2Z).FiniteIndex := by
  rw [Subgroup.finiteIndex_iff]
  show Γ_Gamma2.relIndex Γ_SL2Z ≠ 0
  have hker : φ₂.ker = ⊥ := (MonoidHom.ker_eq_bot_iff φ₂).mpr φ₂_injective
  have hrange : Γ_SL2Z = Subgroup.map φ₂ ⊤ := MonoidHom.range_eq_map φ₂
  rw [hrange]
  unfold Γ_Gamma2
  rw [Subgroup.relIndex_map_map, hker, sup_bot_eq, sup_bot_eq, Subgroup.relIndex_top_right]
  exact (CongruenceSubgroup.instFiniteIndexGamma 2).index_ne_zero

/-- モジュラー群(`SL(2,ℤ)` の像)を `FuchsianGroup` として。 -/
noncomputable def FG_SL2Z : FuchsianGroup := ⟨Γ_SL2Z, inferInstance⟩

/-- `Γ(2)` を `FuchsianGroup` として。 -/
noncomputable def FG_Gamma2 : FuchsianGroup := ⟨Γ_Gamma2, inferInstance⟩

theorem FG_isFiniteIndexIn : IsFiniteIndexIn FG_Gamma2 FG_SL2Z :=
  ⟨Γ_Gamma2_le, Γ_Gamma2_finiteIndex⟩

/-- `T = [[1,1],[0,1]]`(標準生成元)は `Γ(2)` に入らない
(mod 2 で `(0,1)` 成分が `1 ≠ 0`)。 -/
theorem T_not_mem_Gamma2 : ModularGroup.T ∉ CongruenceSubgroup.Gamma 2 := by
  rw [CongruenceSubgroup.Gamma_mem]
  intro ⟨_, h, _, _⟩
  have hT : (ModularGroup.T : Matrix.SpecialLinearGroup (Fin 2) ℤ) 0 1 = 1 := by
    simp [ModularGroup.T]
  rw [hT] at h
  exact absurd h (by decide)

theorem φ₂T_not_mem_Γ_Gamma2 : φ₂ ModularGroup.T ∉ Γ_Gamma2 := by
  unfold Γ_Gamma2
  rw [Subgroup.mem_map]
  rintro ⟨y, hy, hyT⟩
  have hyeq : y = ModularGroup.T := φ₂_injective hyT
  rw [hyeq] at hy
  exact T_not_mem_Gamma2 hy

/-- `SL(2,ℤ)` と `Γ(2)` は(`FuchsianGroup` として)**相異なる**——
`T` が前者には入るが後者には入らない witness で確認。 -/
theorem FG_ne : FG_Gamma2 ≠ FG_SL2Z := by
  intro h
  have hΓ : Γ_Gamma2 = Γ_SL2Z := congrArg FuchsianGroup.Γ h
  exact φ₂T_not_mem_Γ_Gamma2 (hΓ ▸ (⟨ModularGroup.T, rfl⟩ : φ₂ ModularGroup.T ∈ Γ_SL2Z))

open ABC3.Skeleton.CorrHyp in
/-- **[CorrHyp] `Definition 1.1`(`Corr`)の `corrHypInstance` における実現**
——`SL(2,ℤ)` から `Γ(2)` への correspondence(`Γ(2)` 自身を経由、
`Γ(2) ↪ SL(2,ℤ)` と `Γ(2) ↪ Γ(2)`(恒等)の対)。`FEt` が本物の
有限指数条件(`IsFiniteIndexIn`)で埋まっているので、この `C` が
実際に `X`・`Y` 両方に有限指数で入るという非自明な内容を持つ。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corr_witness : Corr corrHypInstance FG_SL2Z FG_Gamma2 :=
  ⟨FG_Gamma2, PLift.up FG_isFiniteIndexIn, PLift.up (isFiniteIndexIn_refl FG_Gamma2)⟩

def corr_witness.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 3, item := "Definition 1.1", sectionId := "corrhyp-def-1-1" }

open ABC3.Skeleton.CorrHyp in
/-- **[CorrHyp] `Definition 1.3`(`Corr.transpose`)の `corrHypInstance` に
おける実現**——`corr_witness`(`SL(2,ℤ) → Γ(2)`)の転置(`Γ(2) → SL(2,ℤ)`)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corr_witness_transpose : Corr corrHypInstance FG_Gamma2 FG_SL2Z :=
  Corr.transpose corrHypInstance corr_witness

def corr_witness_transpose.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.3", sectionId := "corrhyp-def-1-3" }

open ABC3.Skeleton.CorrHyp in
/-- **[CorrHyp] `Definition 1.4`(`Corr.comp'`)の `corrHypInstance` における
実現**——`corr_witness` とその転置の合成(`SL(2,ℤ) → Γ(2) → SL(2,ℤ)`、
`pullback`/`pbFst`/`pbSnd` を実際に経由する)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corr_witness_comp : Corr corrHypInstance FG_SL2Z FG_SL2Z :=
  Corr.comp' corrHypInstance corr_witness corr_witness_transpose

def corr_witness_comp.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.4", sectionId := "corrhyp-def-1-4" }

open ABC3.Skeleton.CorrHyp in
/-- **[CorrHyp] `Definition 1.5`(`IsIsogenous`)の `corrHypInstance` における
非空虚性**——`SL(2,ℤ)` と `Γ(2)` という**相異なる**(`FG_ne`)`Space` が
isogenous であることを、モジュラー群と主合同部分群という教科書的な例で
具体的に示した。

★**sorry 無し**。標準3公理のみ。 -/
theorem isIsogenous_witness : IsIsogenous corrHypInstance FG_SL2Z FG_Gamma2 := ⟨corr_witness⟩

def isIsogenous_witness.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.5", sectionId := "corrhyp-def-1-5" }

/-!
## ★`Definition 1.2`(`IsTrivial`)は claim しない——正直な限界

`FEt X Y := PLift (IsFiniteIndexIn X Y)` は `Prop` を包んだだけの型なので、
**証明無関係性により同じ `X Y` の間の `FEt` は常に唯一**である。ゆえに
`Corr.IsTrivial`(「ある `γ` で `α = β ∘ γ` と書ける」)は
`corrHypInstance` では**任意の correspondence について自動的に真**になる
(`γ := c.α` を取れば `cases` だけで閉じる——下の反例で確認済み)。
これは原文の「trivial」(複数の互いに異なる有限被覆が存在しうる中で、
特別にこの形に書けるもの)を**捉えていない**——`FEt` を `Prop` ではなく
本物の `Type`(例えば被覆変換群を込めた型)にしない限り、`Definition 1.2`
はこの instance で意味のある形では実装できない。★`.src` は付けない。 -/

open ABC3.Skeleton.CorrHyp in
theorem corr_witness_isTrivial_degenerate : Corr.IsTrivial corrHypInstance corr_witness := by
  refine ⟨corr_witness.α, ?_⟩
  cases corr_witness.α
  cases (corrHypInstance.comp corr_witness.β corr_witness.α)
  rfl

end ABC3.Found.CorrHyp
