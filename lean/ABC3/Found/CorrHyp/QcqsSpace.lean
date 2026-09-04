/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.SchemeFEt

/-!
# [CorrHyp] `Lemma 4.1` へ向けた第四歩 —— `Space` を qcqs スキームへ絞る

`ExtLimit.lean` で発見した設計上の要求(`AffineTransitionLimit.lean` の
被覆補題群が `CompactSpace`/`QuasiSeparatedSpace` を要求し、これは一般の
`X : Over BaseK` では保証されない)を受けて、`Space` を「準コンパクト・
準分離(qcqs、= 有限型の反映)なスキーム」に絞った部分型 `QcqsSpace` を導入し、
`FEtK`/`ExtF` の圏構造がこの部分型上でも(`SchemeFEt.lean` の
`FEtK_pullback_compactSpace` 等により)閉じることを確認する。

双曲曲線は原文の定義から常に有限型なので、この絞り込みは原文の意図に忠実
——`corrHypInstance3`(`Space := Over BaseK`、一般のスキームを許す)を
置き換えるのではなく、`Lemma 4.1` を正しく述べるための**新しい**具体化
(`corrHypInstance4`、まだ未着手)の土台として使う。 -/

namespace ABC3.Found.CorrHyp

open CategoryTheory AlgebraicGeometry Limits

/-- `Space` を有限型(qcqs)スキームに絞った部分型——`Lemma 4.1` の被覆補題
適用に要る前提(`ExtLimit.lean` の `extDiagram_map_isAffineHom` の docstring
を見よ)。

★`Type 1`(`Over BaseK` 自身が `Type 1` に住むので、その部分型も
`Type 1`)。 -/
def QcqsSpace : Type 1 := {X : Over BaseK // CompactSpace X.left ∧ QuasiSeparatedSpace X.left}

/-- `QcqsSpace` 上の有限エタール射——`FEtK` をそのまま流用。 -/
noncomputable def QcqsFEt (A B : QcqsSpace) : Type _ := FEtK A.1 B.1

/-- 恒等射。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def QcqsFEt.idFEt (A : QcqsSpace) : QcqsFEt A A := FEtK.idFEt A.1

/-- 合成。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def QcqsFEt.comp {A B C : QcqsSpace} (f : QcqsFEt A B) (g : QcqsFEt B C) :
    QcqsFEt A C := FEtK.comp f g

/-- ファイバー積——`FEtK.pullback` の台が `CompactSpace`/`QuasiSeparatedSpace`
であることを `FEtK_pullback_compactSpace`/`_quasiSeparatedSpace`
(`SchemeFEt.lean`、`B` 側がqcqsであれば十分)で保証し、`QcqsSpace` の
元として返す。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def QcqsFEt.pullback {A B C : QcqsSpace} (f : QcqsFEt A C) (g : QcqsFEt B C) :
    QcqsSpace :=
  haveI := B.2.1
  haveI := B.2.2
  ⟨FEtK.pullback f g, FEtK_pullback_compactSpace f g, FEtK_pullback_quasiSeparatedSpace f g⟩

/-- 第一射影。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def QcqsFEt.pbFst {A B C : QcqsSpace} (f : QcqsFEt A C) (g : QcqsFEt B C) :
    QcqsFEt (QcqsFEt.pullback f g) A := FEtK.pbFst f g

/-- 第二射影。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def QcqsFEt.pbSnd {A B C : QcqsSpace} (f : QcqsFEt A C) (g : QcqsFEt B C) :
    QcqsFEt (QcqsFEt.pullback f g) B := FEtK.pbSnd f g

/-- 係数拡大 `(-)_K`——`ExtF.obj` の台が `CompactSpace`/`QuasiSeparatedSpace`
であることを `ExtF_compactSpace`/`_quasiSeparatedSpace`(`SchemeFEt.lean`)
で保証する。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def QcqsExt (X : QcqsSpace) : QcqsSpace :=
  haveI := X.2.1
  haveI := X.2.2
  ⟨ExtF.obj X.1, ExtF_compactSpace X.1, ExtF_quasiSeparatedSpace X.1⟩

/-- `Ext` の有限エタール射への作用——`extFEt` をそのまま流用。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def QcqsExtFEt {A B : QcqsSpace} (f : QcqsFEt A B) :
    QcqsFEt (QcqsExt A) (QcqsExt B) := extFEt f

/- ★★次の一手(未着手): `corrHypInstance4`(`Space := QcqsSpace`、`FEt :=
QcqsFEt`、`Ext := QcqsExt` 等、`corrHypInstance3` と同じ placeholder パターン)
を構成し、`Lemma 4.1` の statement を型検査する。その上で、`ExtLimit.lean` の
`isLimit_extCone`(`X.left` の qcqs 性を仮定に追加すれば `AffineTransitionLimit.
lean` の被覆補題群が使える状態になる)と `FieldLimit.lean` の
`standardEtalePairRingBaseChange` を組み合わせ、`Z_K` のアフィン開被覆を
étale-locus のレベルで細分してから各片に降下を適用し、遷移射の一致の降下
(`Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`)で貼り合わせる、
という最終組み立てに進む。`corrhyp-goal.md` §4 に記録した計画の続き。 -/

end ABC3.Found.CorrHyp
