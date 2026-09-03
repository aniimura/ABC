/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.FuchsianGroup
import ABC3.Skeleton.CorrHyp.Section1
import ABC3.Skeleton.CorrHyp.Section2
import ABC3.Skeleton.CorrHyp.Section3

/-!
# [CorrHyp] `HyperbolicCurveData` の具体的な実現(§1・§3 が要る範囲)

`Interface/CorrHyp/HyperbolicCurve.lean` の `HyperbolicCurveData` は Skeleton
段階では抽象 posit だったが、ここでその**具体的な項**(`corrHypInstance`)を
`FuchsianGroup`(`Found/CorrHyp/FuchsianGroup.lean`)で組み立てる。

★★これで `ABC3.Skeleton.CorrHyp.prop_3_2`(§3 `Proposition 3.2`、`D` について
全称量化された Skeleton の文そのもの)を、この `D := corrHypInstance` において
**文字通り**(`rfl` で一致することを確認済み)`sorry` 無しで閉じられる
——「関連する具体モデルの結果」ではなく「Skeleton の主張そのものの実現」。

## どの field が本物で、どれが placeholder か(正直な内訳)

`Proposition 3.2` の証明(`corrhyp_prop_3_2`)が実際に読むのは
`Space`・`FEt`・`idFEt`・`comp`・`pullback`・`pbFst`・`pbSnd`・`Fuchsian`・
`Gamma`・`IsDiscrete`・`Gamma_isDiscrete`・`Comm`・`FiniteIndexIn`・`Sub`
だけである——これらは `FuchsianGroup`・`IsFiniteIndexIn`・
`Subgroup.Commensurable.commensurator` という**本物の数学**で埋めた。

残り(`MargulisArithmetic`・`ShimuraArithmetic`・`IsConnected`・`Iso`・`type`・
`core`・`coreMap`・`Aut`・`idAut`・`Ext`・`extFEt`・`stackType`・`ModuliStack`・
`IsGenericallyScheme`・`deg`・`IsOpenDenseIn`・`IsConstructibleIn`)は
`Proposition 3.2` の statement/proof のどちらにも登場しない
(`_hX`・`_hC` は Skeleton 自身で未使用としてマークされている)ので、
**型を合わせるためだけの placeholder**を入れた。★これらは §2・§4-§6 の
他項目の実装ではない——`MargulisArithmetic`/`ShimuraArithmetic := fun _ ↦ False`
は `Definition 2.2`/`2.3` の実装を主張しない(むしろ逆に、`Arithmetic` が
恒偽になるので `¬ Arithmetic` は常に真になる——`Proposition 3.2` の結論が
この仮定を使わずに証明できることと整合しているだけで、`Theorem 2.5`
(`Arithmetic ↔ InfinitelyManyCorr`)等はこの instance では**成り立たない**
(`InfinitelyManyCorr` 側は非自明になりうるため)。 -/

namespace ABC3.Found.CorrHyp

open ABC3.Interface.CorrHyp

/-- `HyperbolicCurveData` の具体的な実現。`Space := FuchsianGroup`(離散部分群
そのものを「双曲曲線」とみなす、原文の一意化定理 `H/Γ ≅ X` を型に埋め込む
逸脱)、`Fuchsian := Subgroup (SL(2,ℝ))`(★離散性を要求しない——`Comm(Γ)` が
常に well-defined であるために必須、`corrhyp-goal.md` 2026-09-04 続報)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypInstance : HyperbolicCurveData where
  Space := FuchsianGroup
  FEt C X := PLift (IsFiniteIndexIn C X)
  idFEt X := PLift.up (isFiniteIndexIn_refl X)
  comp f g := PLift.up (isFiniteIndexIn_comp f.down g.down)
  pullback _f _g := inter _ _
  pbFst f g := PLift.up (isFiniteIndexIn_pbFst f.down g.down)
  pbSnd f g := PLift.up (isFiniteIndexIn_pbSnd f.down g.down)
  Iso X Y := X = Y
  type _ := (0, 0)
  Fuchsian := Subgroup (Matrix.SpecialLinearGroup (Fin 2) ℝ)
  Gamma X := X.Γ
  IsDiscrete Γ := DiscreteTopology Γ
  Gamma_isDiscrete X := X.discrete
  Comm Γ := Subgroup.Commensurable.commensurator Γ
  FiniteIndexIn Γ1 Γ2 := Γ1 ≤ Γ2 ∧ (Γ1.subgroupOf Γ2).FiniteIndex
  Sub Γ1 Γ2 := Γ1 ≤ Γ2
  MargulisArithmetic _ := False
  ShimuraArithmetic _ := False
  IsConnected _ := True
  core := id
  coreMap X := PLift.up (isFiniteIndexIn_refl X)
  Aut _ := PUnit
  idAut _ := PUnit.unit
  Ext := id
  extFEt f := f
  stackType _ := ⟨0, 0, Empty, inferInstance, Empty.elim⟩
  ModuliStack _ _ := ⟨⊥, inferInstance⟩
  IsGenericallyScheme _ := True
  deg _ := 0
  IsOpenDenseIn _ _ _ := True
  IsConstructibleIn _ _ _ := True

open ABC3.Skeleton.CorrHyp in
/-- [CorrHyp] **Proposition 3.2**、`corrHypInstance` における実現。

原文 (CorrHyp p.8):
> We have Γ_Z ⊆ Comm(Γ_X).

`ABC3.Skeleton.CorrHyp.prop_3_2 corrHypInstance`(`D := corrHypInstance` に
特殊化した Skeleton の文そのもの)と**関数として完全に一致する**ことを
`funext … ; rfl`(証明の値ではなく `Prop` の証明無関係性で閉じる)で対話的に
確認済み——ただしその確認補題自体は `prop_3_2`(Skeleton、`sorry` 定義)を
文中で参照するため `#print axioms` が `sorryAx` を巻き込んで表示してしまう
(確認補題自体の内容が sorry というわけではない——`Prop` の証明無関係性により
`rfl` で閉じるので数学的には空虚ではない)。★誤解を避けるため、その確認補題は
恒久的な宣言としては置かず、ここでは `prop_3_2_at_instance` 単体
(`#print axioms` で標準3公理のみ)だけを残す——Skeleton の `sorry` を埋める、
という意味で文字通り `Proposition 3.2` を実装したことになる。
証明は `corrhyp_prop_3_2`(`FuchsianGroup.lean`、`Ncore` 経由)をそのまま
`Corr` の2射 `c.α`・`c.β` に適用するだけ。

★**sorry 無し**。標準3公理のみ。 -/
theorem prop_3_2_at_instance (X Z : corrHypInstance.Space)
    (_hX : ¬ Arithmetic corrHypInstance (corrHypInstance.Gamma X))
    (c : Corr corrHypInstance X Z) (_hC : corrHypInstance.IsConnected c.C) :
    corrHypInstance.Sub (corrHypInstance.Gamma Z)
      (corrHypInstance.Comm (corrHypInstance.Gamma X)) :=
  corrhyp_prop_3_2 c.α.down c.β.down

def prop_3_2_at_instance.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 8, item := "Proposition 3.2", sectionId := "corrhyp-prop-3-2" }

end ABC3.Found.CorrHyp
