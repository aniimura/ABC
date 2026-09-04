/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.Morphisms.Etale
import Mathlib.AlgebraicGeometry.Morphisms.Finite
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Algebra.Rat

/-!
# [CorrHyp] `Lemma 4.1` へ向けた第二歩 —— `Space := Over BaseK` 上の `FEt`・`Ext`

`FieldLimit.lean` は「`K` を有限生成 `k`-部分環の余極限として見る」という
環論側の足場だった。ここでは `HyperbolicCurveData` の `Space`・`FEt`・`Ext` に
実際に使える**スキーム論側**の足場を構成する。

## 設計(記録)

`Space` を「`k`-スキームと`K`-スキームを区別する Sum 型」にする案は、
`corrhyp-goal.md` §4 に記録したとおり複雑さで頓挫した。代わりに
**「係数の制限」**(restriction of scalars)を使う: `K`-スキーム `Y → Spec K` は
`Spec K → Spec k` と合成すれば常に `k`-スキーム `Y → Spec k` として見られる。
ゆえに `Space := Over (Spec k)` 一本(`Over BaseK`)に固定し、`Ext X := X ×_{Spec k}
Spec K`(のち `Spec K → Spec k` と合成して `Over BaseK` に戻す)とすればよい。

`k := ℚ`・`K := ℝ` を選んだ(`AdjoinRoot` 経由の `ℚ(i)` は `Field` instance に
`Fact (Irreducible p)` を要求し、余計な複雑さを生むため回避——`Lemma 4.1` は
`K/k` の代数性を要求しない)。

`FEt`(有限 étale 射)は `AlgebraicGeometry.Etale`・`AlgebraicGeometry.IsFinite`
という mathlib の `MorphismProperty` 基盤(合成・恒等・base change 安定性が
すべて `instance` として自動で手に入る)にそのまま乗る——`idFEt`・`comp`・
`pullback`・`pbFst`・`pbSnd` はいずれも `infer_instance` と
`MorphismProperty.pullback_fst`/`pullback_snd` だけで閉じた。

`Ext`(射側 `extFEt` まで込み)は、手作りの `pullback.map`(`Limits.pullback.map`)
ではなく **mathlib 自身の `Over.pullback f ⋙ Over.map f`**(`f : specK ⟶ BaseK`)
という関手の合成として構成すると、`CategoryTheory.MorphismProperty.overPullbackMap`
が base change 安定性からの遺伝を丸ごと与えてくれる——`.map` の `Over.map` 側は
`.left` を変えない(post-compose するだけ)ので、`overPullbackMap` 一発で
`IsFinite`/`Etale` の遺伝が閉じる。手作りの `pullback.map`+`IsPullback` 貼り合わせ
(`IsPullback.of_right` 等)を試みたが、`Over.pullback`/`overPullbackMap` の組が
遥かに短く済んだ(配管の教訓、`tools/lean-idioms.md` に追記予定)。 -/

namespace ABC3.Found.CorrHyp

open CategoryTheory AlgebraicGeometry Limits

/-- `Spec ℚ`(基礎体 `k` のスキーム)。 -/
noncomputable def BaseK : Scheme := Spec (CommRingCat.of ℚ)

/-- `Space := Over BaseK`(`k`-スキーム全体)における有限 étale 射
(`Definition 1.1` の `FEt` に対応)。 -/
noncomputable def FEtK (X Y : Over BaseK) : Type _ :=
  {f : X ⟶ Y // IsFinite f.left ∧ Etale f.left}

/-- 恒等射は有限 étale。 -/
noncomputable def FEtK.idFEt (X : Over BaseK) : FEtK X X :=
  ⟨𝟙 X, by simp; infer_instance, by simp; infer_instance⟩

/-- 有限 étale 射の合成は有限 étale(mathlib の合成 instance から自動)。 -/
noncomputable def FEtK.comp {A B C : Over BaseK} (f : FEtK A B) (g : FEtK B C) : FEtK A C :=
  ⟨f.1 ≫ g.1, by
    have := f.2.1; have := f.2.2; have := g.2.1; have := g.2.2
    simp; infer_instance, by
    have := f.2.1; have := f.2.2; have := g.2.1; have := g.2.2
    simp; infer_instance⟩

/-- `f : FEtK A C`・`g : FEtK B C` のファイバー積(`Definition 1.4` の `C₁ ×_Y C₂`)
——`Over` 圏側の pullback 機構は経由せず、台となるスキームの `Limits.pullback` を
直接 `Over.mk` で包む。 -/
noncomputable def FEtK.pullback {A B C : Over BaseK} (f : FEtK A C) (g : FEtK B C) :
    Over BaseK :=
  Over.mk (Limits.pullback.fst f.1.left g.1.left ≫ A.hom)

/-- ファイバー積からの第一射影は有限 étale(`MorphismProperty.pullback_fst` が
`g` の性質から遺伝させる)。 -/
noncomputable def FEtK.pbFst {A B C : Over BaseK} (f : FEtK A C) (g : FEtK B C) :
    FEtK (FEtK.pullback f g) A :=
  ⟨Over.homMk (Limits.pullback.fst f.1.left g.1.left) rfl,
    MorphismProperty.pullback_fst f.1.left g.1.left g.2.1,
    MorphismProperty.pullback_fst f.1.left g.1.left g.2.2⟩

/-- ファイバー積からの第二射影は有限 étale(`MorphismProperty.pullback_snd` が
`f` の性質から遺伝させる)。 -/
noncomputable def FEtK.pbSnd {A B C : Over BaseK} (f : FEtK A C) (g : FEtK B C) :
    FEtK (FEtK.pullback f g) B :=
  ⟨Over.homMk (Limits.pullback.snd f.1.left g.1.left) (by
      show Limits.pullback.snd f.1.left g.1.left ≫ B.hom =
        Limits.pullback.fst f.1.left g.1.left ≫ A.hom
      rw [← Over.w f.1, ← Over.w g.1, ← Category.assoc, ← Category.assoc, pullback.condition]),
    MorphismProperty.pullback_snd f.1.left g.1.left f.2.1,
    MorphismProperty.pullback_snd f.1.left g.1.left f.2.2⟩

/-- `Spec ℝ`(拡大体 `K` のスキーム)。 -/
noncomputable def specK : Scheme := Spec (CommRingCat.of ℝ)

/-- `Spec K → Spec k`(`ℚ → ℝ` の構造射、係数拡大 `(−)_K` の相方)。 -/
noncomputable def toBaseK : specK ⟶ BaseK :=
  Spec.map (CommRingCat.ofHom (algebraMap ℚ ℝ))

/-- **`k`-スキーム `X` を `K` へ係数拡大し、また `k`-スキームとして見る関手**
(`(−)_K`、`Definition`/`§4` 冒頭)——mathlib 自身の `Over.pullback toBaseK`
(`Over BaseK ⥤ Over specK`、base change)と `Over.map toBaseK`
(`Over specK ⥤ Over BaseK`、`Spec K → Spec k` との後合成)の合成。 -/
noncomputable def ExtF : Over BaseK ⥤ Over BaseK :=
  Over.pullback toBaseK ⋙ Over.map toBaseK

/-- `Over.map` は `.left`(台となるスキームの射)を変えない——`ExtF.map f` の
有限 étale 性を `(Over.pullback toBaseK).map f` の有限 étale性に帰着させる鍵。 -/
theorem ExtF.map_left_eq {A B : Over BaseK} (f : A ⟶ B) :
    (ExtF.map f).left = ((Over.pullback toBaseK).map f).left := by
  simp [ExtF, Over.map]

/-- `ExtF` は有限性を保つ(`overPullbackMap` が base change 安定性から与える)。 -/
noncomputable instance ExtF.instFinite {A B : Over BaseK} (f : A ⟶ B) [IsFinite f.left] :
    IsFinite (ExtF.map f).left := by
  rw [ExtF.map_left_eq]
  exact MorphismProperty.overPullbackMap toBaseK f ‹IsFinite f.left›

/-- `ExtF` は étale 性を保つ(`overPullbackMap` が base change 安定性から与える)。 -/
noncomputable instance ExtF.instEtale {A B : Over BaseK} (f : A ⟶ B) [Etale f.left] :
    Etale (ExtF.map f).left := by
  rw [ExtF.map_left_eq]
  exact MorphismProperty.overPullbackMap toBaseK f ‹Etale f.left›

/-- **`extFEt`**——`Ext` の有限 étale 射への作用(`HyperbolicCurveData.extFEt`
そのもの)。`ExtF` の関手性(合成・恒等の保存)から `FEtK` の圏構造との両立は
自動。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def extFEt {A B : Over BaseK} (f : FEtK A B) : FEtK (ExtF.obj A) (ExtF.obj B) :=
  have := f.2.1
  have := f.2.2
  ⟨ExtF.map f.1, ExtF.instFinite f.1, ExtF.instEtale f.1⟩

end ABC3.Found.CorrHyp
