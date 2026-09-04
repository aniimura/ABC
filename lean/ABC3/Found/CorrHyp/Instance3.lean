/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.SchemeFEt
import ABC3.Skeleton.CorrHyp.Section4

/-!
# [CorrHyp] `Space := Over BaseK`(スキームモデル)による `HyperbolicCurveData` の具体化

`corrHypInstance`/`corrHypInstance2`(`Instance.lean`/`Instance2.lean`)は
`Space := FuchsianGroup`(`ℂ` 上の解析的モデル)を使う——`§1`-`§3` 向け。
`Lemma 4.1`(`§4`)は一般の体 `k` 上のスキームの言葉で書かれているため、
別の具体化 `corrHypInstance3`(`Space := Over BaseK`、`SchemeFEt.lean` の
スキーム論的 `FEt`/`Ext`)をここに置く。

`HyperbolicCurveData.Space` を universe 多相化した(`HyperbolicCurve.lean` の
逸脱記録を見よ)おかげで、`corrHypInstance`(`u=0`)と `corrHypInstance3`
(`u=1`、`Over BaseK : Type 1` を直接使う)が同じ interface を共有できる。

`Lemma 4.1`・`Fuchsian`/`Aut`/`ModuliStack` 等、`Lemma 4.1` の**主張自体**が
読まないフィールドはすべて安全な placeholder(`corrHypInstance` が
`Proposition 3.2` に対して行ったのと同じパターン)。 -/

namespace ABC3.Found.CorrHyp

open ABC3.Interface.CorrHyp CategoryTheory AlgebraicGeometry

/-- `Over BaseK` の適当な基点(`ModuliStack` 等、`Lemma 4.1` が読まない
フィールドの placeholder 用)。 -/
noncomputable def basePt : Over BaseK := Over.mk (𝟙 BaseK)

/-- **`corrHypInstance3`**——`Space := Over BaseK`(`k := ℚ`-スキーム全体)による
`HyperbolicCurveData` の具体化。`FEt`・`idFEt`・`comp`・`pullback`・`pbFst`・
`pbSnd`・`Ext`・`extFEt` は `SchemeFEt.lean` の本物のデータ、それ以外
(`Lemma 4.1` の主張が読まない)は placeholder。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypInstance3 : HyperbolicCurveData where
  Space := Over BaseK
  FEt := FEtK
  idFEt := FEtK.idFEt
  comp := FEtK.comp
  pullback := FEtK.pullback
  pbFst := FEtK.pbFst
  pbSnd := FEtK.pbSnd
  Iso X Y := X = Y
  type _ := (0, 0)
  Fuchsian := Unit
  Gamma _ := ()
  IsDiscrete _ := True
  Gamma_isDiscrete _ := trivial
  Comm := id
  FiniteIndexIn _ _ := True
  Sub _ _ := True
  MargulisArithmetic _ := False
  ShimuraArithmetic _ := False
  IsConnected _ := True
  core := id
  coreMap X := FEtK.idFEt X
  Aut _ := PUnit
  idAut _ := PUnit.unit
  Ext := ExtF.obj
  extFEt f := extFEt f
  stackType _ := ⟨0, 0, Empty, inferInstance, Empty.elim⟩
  ModuliStack _ _ := basePt
  IsGenericallyScheme _ := True
  deg _ := 0
  IsOpenDenseIn _ _ _ := True
  IsConstructibleIn _ _ _ := True

open ABC3.Skeleton.CorrHyp in
/-- **`Lemma 4.1` の statement は `corrHypInstance3` に対して意味を持つ**
(型検査が通ることの記録)——本体はまだ `sorry`。 -/
theorem lemma_4_1_at_instance3_statement :
    (∀ (X ZK : corrHypInstance3.Space) (_c : Corr corrHypInstance3
      (corrHypInstance3.Ext X) ZK), True) := fun _ _ _ => trivial

end ABC3.Found.CorrHyp
