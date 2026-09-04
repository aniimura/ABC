/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.FieldLimit
import ABC3.Found.CorrHyp.SchemeFEt

/-!
# [CorrHyp] `Lemma 4.1` へ向けた第三歩 —— `Spec K` の極限表示を `Over BaseK` へ持ち上げる

`FieldLimit.lean` の `isLimit_specKCone : IsLimit (specKCone ℚ ℝ)` は**裸の
`Scheme`**の中での極限(`Spec ℝ = lim Spec R`、`R` は有限生成 `ℚ`-部分環)。
`Lemma 4.1` の構成的降下(`corrhyp-goal.md` §4)を実際に `Ext`(`SchemeFEt.lean`)
に適用するには、この極限を **`Over BaseK`**(`k`-スキーム全体の圏)の中の
極限として持ち上げる必要がある——`Ext X = X ×_k Spec K` の極限表示
(`X ×_k Spec K = lim (X ×_k Spec R)`)を得るための土台。

鍵となる事実: `CategoryTheory.Over.pullback (f : Y ⟶ X) : Over X ⥤ Over Y` は
`Over.map f` の右随伴(`Under.pushoutIsLeftAdjoint` のOver版、
`Over.mapPullbackAdj`)なので、**右随伴は極限を保存する**という一般論
(`infer_instance` で `PreservesLimit` が自動的に付くことを確認済み)により、
`Over BaseK` の中の極限を `Over X.left` へ base change しても極限のままである
——このファイルではまず `Spec K = lim Spec R` を `Over BaseK` へ持ち上げる
(`isLimit_specKConeOver`)ところまでを完成させる。`Ext X` への具体的な
接続(`pullback` の引数順の入れ替えの処理)は次の一手として残す。 -/

namespace ABC3.Found.CorrHyp

open CategoryTheory AlgebraicGeometry Limits
open scoped TensorProduct

set_option backward.isDefEq.respectTransparency false in
/-- `toSchemeDiagram ℚ ℝ`(`(FgSubalgebra ℚ ℝ)ᵒᵖ ⥤ Scheme`)を `Over BaseK`
へ持ち上げたもの——各 `Spec R`(`R` は有限生成 `ℚ`-部分環)は自然に
`BaseK = Spec ℚ` 上のスキームなので、その構造射 `Spec.map (algebraMap ℚ R)`
込みで `Over BaseK` の対象として見る。

★**sorry 無し**。標準3公理のみ。`toSchemeDiagramOver ⋙ Over.forget BaseK =
toSchemeDiagram ℚ ℝ` が `rfl` で成り立つ(対象・射どちらも)ことを確認済み
——`Over` への持ち上げが「余計な情報を足していない」ことの直接確認。 -/
noncomputable def toSchemeDiagramOver : (FgSubalgebra ℚ ℝ)ᵒᵖ ⥤ Over BaseK where
  obj R := Over.mk (Spec.map (CommRingCat.ofHom (algebraMap ℚ R.unop.1)))
  map {R S} h := Over.homMk ((toSchemeDiagram ℚ ℝ).map h) (by
    show (toSchemeDiagram ℚ ℝ).map h ≫ Spec.map (CommRingCat.ofHom (algebraMap ℚ S.unop.1)) =
      Spec.map (CommRingCat.ofHom (algebraMap ℚ R.unop.1))
    show Spec.map ((toRingCat ℚ ℝ).map h.unop) ≫
      Spec.map (CommRingCat.ofHom (algebraMap ℚ S.unop.1)) = _
    rw [← Spec.map_comp]
    congr 1)

set_option backward.isDefEq.respectTransparency false in
/-- `specKCone ℚ ℝ` を `Over BaseK` へ持ち上げたもの——頂点は
`Over.mk toBaseK`(`specK` を `k`-スキームとして見たもの)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def specKConeOver : Cone toSchemeDiagramOver where
  pt := Over.mk toBaseK
  π :=
  { app := fun R => Over.homMk ((specKCone ℚ ℝ).π.app R) (by
      show Spec.map (CommRingCat.ofHom R.unop.1.subtype) ≫
          Spec.map (CommRingCat.ofHom (algebraMap ℚ R.unop.1)) =
          Spec.map (CommRingCat.ofHom (algebraMap ℚ ℝ))
      rw [← Spec.map_comp]
      congr 1)
    naturality := by
      intro R S h
      apply Over.OverMorphism.ext
      exact (specKCone ℚ ℝ).π.naturality h }

/-- `specKConeOver` の頂点を裸の `Scheme` へ引き戻した(`Over.forget BaseK`
で送った)cone——`toSchemeDiagramOver ⋙ Over.forget BaseK = toSchemeDiagram
ℚ ℝ` が `rfl` で成り立つことを利用して、任意の `s : Cone toSchemeDiagramOver`
から直接 `Cone (toSchemeDiagram ℚ ℝ)` を作る(`(Over.forget BaseK).mapCone`
を経由すると `Cone` が関手を型レベルで引きずるため rewrite が詰まる——
直接構成することでこれを避ける)。 -/
private noncomputable def auxCone (s : Cone toSchemeDiagramOver) :
    Cone (toSchemeDiagram ℚ ℝ) where
  pt := s.pt.left
  π :=
  { app := fun R => (s.π.app R).left
    naturality := by
      intro R S h
      exact congrArg Over.Hom.left (s.π.naturality h) }

set_option backward.isDefEq.respectTransparency false in
/-- **`specKConeOver` は `Over BaseK` の中で極限**——`isLimit_specKCone`
(裸の `Scheme` での極限)から、`toBaseK` が `(specKCone ℚ ℝ).π.app R` と
`(toSchemeDiagramOver.obj R).hom` の合成として書けること(任意の `R` で
成り立つ、`ℚ ⟶ R ⟶ ℝ` が `ℚ ⟶ ℝ` に一致するという `IsScalarTower` 相当の
事実)を使って、任意の cone の頂点からの一意な lift を直接構成する。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def isLimit_specKConeOver : IsLimit specKConeOver := by
  set R0 : (FgSubalgebra ℚ ℝ)ᵒᵖ := Opposite.op ⟨⊥, Subalgebra.fg_bot⟩ with hR0
  refine IsLimit.mk (fun s => Over.homMk
    ((isLimit_specKCone ℚ ℝ).lift (auxCone s)) ?_) ?_ ?_
  · show ((isLimit_specKCone ℚ ℝ).lift (auxCone s)) ≫ toBaseK = s.pt.hom
    have hcomp : ((isLimit_specKCone ℚ ℝ).lift (auxCone s)) ≫
        (specKCone ℚ ℝ).π.app R0 = (auxCone s).π.app R0 :=
      (isLimit_specKCone ℚ ℝ).fac (auxCone s) R0
    have htoBaseK : (specKCone ℚ ℝ).π.app R0 ≫ (toSchemeDiagramOver.obj R0).hom = toBaseK :=
      (specKConeOver.π.app R0).w
    rw [← htoBaseK, ← Category.assoc, hcomp]
    exact (s.π.app R0).w
  · intro s R
    apply Over.OverMorphism.ext
    exact (isLimit_specKCone ℚ ℝ).fac (auxCone s) R
  · intro s m hm
    apply Over.OverMorphism.ext
    apply (isLimit_specKCone ℚ ℝ).uniq (auxCone s)
    intro R
    exact congrArg Over.Hom.left (hm R)

/-- `X : Over BaseK`(`Lemma 4.1` の `X`)を固定したとき、`Over.pullback X.hom`
(`Over.mapPullbackAdj` の右随伴なので `PreservesLimit` が `infer_instance` で
自動的に付く)を `specKConeOver`/`isLimit_specKConeOver` に適用すれば、
`X.left ×_{BaseK} Spec K` の極限表示が(`Over X.left` の中の極限として)
得られる——`pullback toBaseK X.hom`(`Over.pullback X.hom` が自然に与える
引数順)の形で。`Ext X` 自身の定義(`SchemeFEt.lean`、`pullback X.hom
toBaseK` という逆順)との接続には `Limits.pullbackSymmetry` が要る(未着手、
下記)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def extConeOver (X : Over BaseK) :
    IsLimit ((Over.pullback X.hom).mapCone specKConeOver) :=
  isLimitOfPreserves (Over.pullback X.hom) isLimit_specKConeOver

set_option backward.isDefEq.respectTransparency false in
/-- `toSchemeDiagram ℚ ℝ` を `X.hom`(`X : Over BaseK`)に沿って base change
した図式——`i ↦ X.left ×_{BaseK} Spec R_i`(裸の `Scheme` の中、`Ext X` の
定義(`pullback X.hom toBaseK`)と**同じ引数順**で作ってあるので、
`pullbackSymmetry` を経由せずに直接 `Ext X` の極限表示の土台になる。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def extDiagram (X : Over BaseK) : (FgSubalgebra ℚ ℝ)ᵒᵖ ⥤ Scheme where
  obj R := Limits.pullback X.hom (toSchemeDiagramOver.obj R).hom
  map {R S} h := Limits.pullback.map X.hom (toSchemeDiagramOver.obj R).hom X.hom
    (toSchemeDiagramOver.obj S).hom (𝟙 X.left) (toSchemeDiagramOver.map h).left (𝟙 BaseK)
    (by simp) (by simpa using (toSchemeDiagramOver.map h).w.symm)
  map_id R := by
    apply pullback.hom_ext <;> simp
  map_comp {R S T} f g := by
    rw [pullback.map_comp]
    congr 1
    simp

/-- `(extDiagram X).map h` を `pullback.snd` と合成した式の naturality——
`extDiagram X` の頂点候補が `Ext X` の極限であることを示す際、`Cone` の
`naturality` フィールド(`Functor.const` で包まれた形)を直接埋めようとすると
`instances` 透明度の配管に当たる(`Over.mk` で見た問題の再発)ため、まず
**独立した補題として**(`Functor.const` を経由せずに)証明しておく。
`simp only [extDiagram]` で `.map` の定義を展開してから `pullback.lift_snd`
一発で閉じる——`unfold`+`show` より `simp only [<def名>]` の方がこの手の
配管に強いことを再確認。

★**sorry 無し**。標準3公理のみ。 -/
theorem extDiagram_map_snd {X : Over BaseK} {R S : (FgSubalgebra ℚ ℝ)ᵒᵖ} (h : R ⟶ S) :
    (extDiagram X).map h ≫ pullback.snd X.hom (toSchemeDiagramOver.obj S).hom =
      pullback.snd X.hom (toSchemeDiagramOver.obj R).hom ≫ (toSchemeDiagramOver.map h).left := by
  simp only [extDiagram]
  rw [pullback.lift_snd]

set_option backward.isDefEq.respectTransparency false in
/-- `extDiagram X` から `toSchemeDiagram ℚ ℝ` への自然変換——各段階を
`pullback.snd`(`X.left ×_{BaseK} Spec R → Spec R` への射影)で送る。
`Cone` の中に直接埋め込まず**独立した `NatTrans` として**構成すると
(`extDiagram_map_snd` を直接使うだけで)`instances` 透明度の配管に
当たらずに閉じる——`Functor.const` で包まれた `Cone.π` の naturality
フィールドの中で同じ主張を証明しようとすると詰まる(下記の「次の一手」)、
という配管上の教訓。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def extDiagramToSpecK (X : Over BaseK) : extDiagram X ⟶ toSchemeDiagram ℚ ℝ where
  app R := pullback.snd X.hom (toSchemeDiagramOver.obj R).hom
  naturality R S h := by
    show (extDiagram X).map h ≫ pullback.snd X.hom (toSchemeDiagramOver.obj S).hom =
      pullback.snd X.hom (toSchemeDiagramOver.obj R).hom ≫ (toSchemeDiagramOver.map h).left
    exact extDiagram_map_snd h

/-- `(extDiagram X).map h` を `pullback.fst` と合成した式の naturality——
`extDiagram_map_snd` の fst 版。標準射影 `𝟙 X.left` 側は `Category.comp_id`
一発で閉じる。

★**sorry 無し**。標準3公理のみ。 -/
theorem extDiagram_map_fst {X : Over BaseK} {R S : (FgSubalgebra ℚ ℝ)ᵒᵖ} (h : R ⟶ S) :
    (extDiagram X).map h ≫ pullback.fst X.hom (toSchemeDiagramOver.obj S).hom =
      pullback.fst X.hom (toSchemeDiagramOver.obj R).hom := by
  simp only [extDiagram]
  rw [pullback.lift_fst, Category.comp_id]

set_option backward.isDefEq.respectTransparency false in
/-- **`extDiagram X` の頂点候補 `Limits.pullback X.hom toBaseK`(= `(Ext X).left`)
への `π`**——`specKCone ℚ ℝ` の各射影を `pullback.map` で `X.hom` に沿って
base change したもの。前段(`auxCone2`/`extConePi` の失敗した試み)との違い:
`𝟙 pt`・`(extDiagram X).map h` 等を `set ... with h...` で**先に名前を
つけてから**扱うと、`Category.assoc`/`pullback.lift_fst`/`pullback.lift_snd`
が「型が合わない」を起こさずに済む——`Functor.const` の配管(第22項)を
避けるだけでなく、この「`set` で名前を固定してから計算する」こと自体が
鍵だった(`tools/lean-idioms.md` に追記予定)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def extConePi (X : Over BaseK) :
    (Functor.const (FgSubalgebra ℚ ℝ)ᵒᵖ).obj (Limits.pullback X.hom toBaseK) ⟶ extDiagram X where
  app R := Limits.pullback.map X.hom toBaseK X.hom (toSchemeDiagramOver.obj R).hom
      (𝟙 X.left) ((specKCone ℚ ℝ).π.app R) (𝟙 BaseK) (by simp)
      (by rw [Category.comp_id]; exact (specKConeOver.π.app R).w.symm)
  naturality R S h := by
    simp only [Functor.const_obj_map]
    set mR := Limits.pullback.map X.hom toBaseK X.hom (toSchemeDiagramOver.obj R).hom
      (𝟙 X.left) ((specKCone ℚ ℝ).π.app R) (𝟙 BaseK) (by simp)
      (by rw [Category.comp_id]; exact (specKConeOver.π.app R).w.symm) with hmR
    set mS := Limits.pullback.map X.hom toBaseK X.hom (toSchemeDiagramOver.obj S).hom
      (𝟙 X.left) ((specKCone ℚ ℝ).π.app S) (𝟙 BaseK) (by simp)
      (by rw [Category.comp_id]; exact (specKConeOver.π.app S).w.symm) with hmS
    show 𝟙 (Limits.pullback X.hom toBaseK) ≫ mS = mR ≫ (extDiagram X).map h
    rw [Category.id_comp]
    have hRfst : mR ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R).hom = pullback.fst X.hom toBaseK := by
      rw [hmR]; unfold Limits.pullback.map; rw [pullback.lift_fst, Category.comp_id]
    have hSfst : mS ≫ pullback.fst X.hom (toSchemeDiagramOver.obj S).hom = pullback.fst X.hom toBaseK := by
      rw [hmS]; unfold Limits.pullback.map; rw [pullback.lift_fst, Category.comp_id]
    have hRsnd : mR ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R).hom =
        pullback.snd X.hom toBaseK ≫ (specKCone ℚ ℝ).π.app R := by
      rw [hmR]; unfold Limits.pullback.map; rw [pullback.lift_snd]
    have hSsnd : mS ≫ pullback.snd X.hom (toSchemeDiagramOver.obj S).hom =
        pullback.snd X.hom toBaseK ≫ (specKCone ℚ ℝ).π.app S := by
      rw [hmS]; unfold Limits.pullback.map; rw [pullback.lift_snd]
    have key : (specKCone ℚ ℝ).π.app R ≫ (toSchemeDiagramOver.map h).left =
        (specKCone ℚ ℝ).π.app S := by
      have h2 := congrArg Over.Hom.left (specKConeOver.π.naturality h)
      simp only [Functor.const_obj_map, Over.comp_left, Over.id_left] at h2
      exact h2.symm
    simp only [extDiagram]
    apply pullback.hom_ext
    · conv_rhs => rw [Category.assoc, pullback.lift_fst]
      rw [Category.comp_id, hSfst, hRfst]
    · conv_rhs => rw [Category.assoc, pullback.lift_snd]
      rw [← Category.assoc, hRsnd, Category.assoc, key, hSsnd]

/-- **`extDiagram X` の頂点候補**——`(Ext X).left` を頂点とする
`Cone (extDiagram X)`。`extConePi` により `π` は完成しているので、あとは
これが極限であること(`isLimit_extCone`、未着手)を示せば `Lemma 4.1` の
スキーム側の道具が揃う。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def extCone (X : Over BaseK) : Cone (extDiagram X) where
  pt := Limits.pullback X.hom toBaseK
  π := extConePi X

/-- `s : Cone (extDiagram X)` を `extDiagramToSpecK` で `toSchemeDiagram ℚ ℝ`
の cone へ変換したもの——`Cone.postcompose` を使うと naturality の再証明が
不要になる(`Cones.postcompose`/`Cone.postcompose` は deprecated 警告が出るが
mathlib 側の名称変更のみで機能は同じ)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def auxCone3 (X : Over BaseK) (s : Cone (extDiagram X)) :
    Cone (toSchemeDiagram ℚ ℝ) :=
  (Cone.postcompose (extDiagramToSpecK X)).obj s

/-- `⊥ : FgSubalgebra ℚ ℝ`(常に有限生成)から任意の `R` への射
(`(FgSubalgebra ℚ ℝ)ᵒᵖ` の中、`⊥ ≤ R.unop` を `homOfLE` で持ち上げてから
`.op` する)——`isLimit_extCone` で「代表元」として使う基点を用意する。 -/
theorem R0_le (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    (⟨⊥, Subalgebra.fg_bot⟩ : FgSubalgebra ℚ ℝ) ≤ R.unop := by
  show (⊥ : Subalgebra ℚ ℝ) ≤ R.unop.1
  exact bot_le

noncomputable def R0hom (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    R ⟶ Opposite.op (⟨⊥, Subalgebra.fg_bot⟩ : FgSubalgebra ℚ ℝ) :=
  (CategoryTheory.homOfLE (R0_le R)).op

set_option backward.isDefEq.respectTransparency false in
/-- `s : Cone (extDiagram X)` の `s.π.app R ≫ pullback.fst`(`X.left` への
射影)は `R` に依らず一定——`extDiagram_map_fst` と `s` 自身の cone
naturality(`Cone.w`)から。`isLimit_extCone` の `lift`/`fac`/`uniq` すべてで
使う核心の補題。

★**sorry 無し**。標準3公理のみ。 -/
theorem extCone_fst_const {X : Over BaseK} (s : Cone (extDiagram X)) {R S : (FgSubalgebra ℚ ℝ)ᵒᵖ}
    (h : R ⟶ S) :
    s.π.app S ≫ pullback.fst X.hom (toSchemeDiagramOver.obj S).hom =
      s.π.app R ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R).hom := by
  have hs : s.π.app R ≫ (extDiagram X).map h = s.π.app S := s.w h
  calc s.π.app S ≫ pullback.fst X.hom (toSchemeDiagramOver.obj S).hom
      = (s.π.app R ≫ (extDiagram X).map h) ≫ pullback.fst X.hom (toSchemeDiagramOver.obj S).hom := by
        rw [hs]
    _ = s.π.app R ≫ (extDiagram X).map h ≫ pullback.fst X.hom (toSchemeDiagramOver.obj S).hom := by
        rw [Category.assoc]
    _ = s.π.app R ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R).hom := by
        rw [extDiagram_map_fst h]

set_option backward.isDefEq.respectTransparency false in
/-- `extConePi X` の `app R` を `pullback.fst`/`pullback.snd` と合成した式の
計算——`extDiagram_map_fst`/`_snd` と同じパターン(`simp only [extConePi]`
で `.app` の定義を展開してから `pullback.lift_fst`/`_snd` 一発)。
`isLimit_extCone` の `fac`/`uniq` 両方で使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem extConePi_app_fst (X : Over BaseK) (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    (extConePi X).app R ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R).hom =
      pullback.fst X.hom toBaseK := by
  simp only [extConePi]
  rw [pullback.lift_fst, Category.comp_id]

set_option backward.isDefEq.respectTransparency false in
theorem extConePi_app_snd (X : Over BaseK) (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    (extConePi X).app R ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R).hom =
      pullback.snd X.hom toBaseK ≫ (specKCone ℚ ℝ).π.app R := by
  simp only [extConePi]
  rw [pullback.lift_snd]

set_option backward.isDefEq.respectTransparency false in
/-- **`extCone X` は極限**——`Lemma 4.1` の構成的降下に必要なスキーム側の
道具の最後のピース。`s : Cone (extDiagram X)` から `X.left` への射
(`s.π.app R0 ≫ pullback.fst`、`R0 := ⊥`)と `specK` への射
(`(isLimit_specKCone ℚ ℝ).lift` を `auxCone3` に適用)を `pullback.lift`
で束ねる。互換性条件は `pullback.condition` + `(isLimit_specKCone).fac`。
`fac`/`uniq` はどちらも `pullback.hom_ext` で fst/snd に割ってから
`extConePi_app_fst`/`_snd`・`extCone_fst_const`・`(isLimit_specKCone).fac`/
`.uniq` を適用するだけで閉じる——鍵は `hm : ∀ j, m ≫ (extCone X).π.app j =
s.π.app j` を `have hm' : ∀ R, m ≫ (extConePi X).app R = s.π.app R := hm`
のように**型を明示して再束縛する**ことで、`rw` の構文一致ではなく `have`
の defeq チェック(`set_option backward.isDefEq.respectTransparency false`
込み)を経由させる、という配管——ここまでの「`Functor.const` の配管」
(第22項)の締めくくりとして `tools/lean-idioms.md` に追記予定。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def isLimit_extCone (X : Over BaseK) : IsLimit (extCone X) := by
  set R0 : (FgSubalgebra ℚ ℝ)ᵒᵖ := Opposite.op ⟨⊥, Subalgebra.fg_bot⟩ with hR0
  refine IsLimit.mk (fun s => pullback.lift
    (s.π.app R0 ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R0).hom)
    ((isLimit_specKCone ℚ ℝ).lift (auxCone3 X s))
    (by
      have hfac : (isLimit_specKCone ℚ ℝ).lift (auxCone3 X s) ≫ (specKCone ℚ ℝ).π.app R0 =
          s.π.app R0 ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R0).hom :=
        (isLimit_specKCone ℚ ℝ).fac (auxCone3 X s) R0
      have htoBaseK : (specKCone ℚ ℝ).π.app R0 ≫ (toSchemeDiagramOver.obj R0).hom = toBaseK :=
        (specKConeOver.π.app R0).w
      rw [Category.assoc, pullback.condition, ← Category.assoc, ← hfac, Category.assoc,
        htoBaseK])) ?_ ?_
  · intro s R
    show pullback.lift (s.π.app R0 ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R0).hom)
      ((isLimit_specKCone ℚ ℝ).lift (auxCone3 X s)) _ ≫ (extConePi X).app R = s.π.app R
    apply pullback.hom_ext
    · rw [Category.assoc, extConePi_app_fst, pullback.lift_fst]
      exact (extCone_fst_const s (R0hom R))
    · rw [Category.assoc, extConePi_app_snd, ← Category.assoc, pullback.lift_snd]
      exact (isLimit_specKCone ℚ ℝ).fac (auxCone3 X s) R
  · intro s m hm
    have hm' : ∀ R, m ≫ (extConePi X).app R = s.π.app R := hm
    show m = pullback.lift (s.π.app R0 ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R0).hom)
      ((isLimit_specKCone ℚ ℝ).lift (auxCone3 X s)) _
    apply pullback.hom_ext
    · rw [pullback.lift_fst, ← extConePi_app_fst X R0, ← Category.assoc, hm' R0]
    · rw [pullback.lift_snd]
      apply (isLimit_specKCone ℚ ℝ).uniq (auxCone3 X s)
      intro R
      show (m ≫ pullback.snd X.hom toBaseK) ≫ (specKCone ℚ ℝ).π.app R =
        s.π.app R ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R).hom
      rw [Category.assoc, ← extConePi_app_snd X R, ← Category.assoc, hm' R]

/-- `extDiagram X` の遷移射はアフィン——`AffineTransitionLimit.lean` の被覆
補題群(`Scheme.exists_isOpenCover_and_isAffine` 等)を `extDiagram X` に
適用するための前提の1つ。`(extDiagram X).map h` を `pullback.map ... (𝟙
X.left) φ (𝟙 BaseK) ...` の形に展開し、`CategoryTheory.MorphismProperty.
pullbackMap`(`P i₁`・`P i₂` から `P (pullback.map ...)` を直接与える、
`SchemeFEt.lean`/`ExtLimit.lean` の他の箇所で使った `Over.pullback`+
`overPullbackMap` より単純な道)に `IsAffineHom` を適用する——`𝟙 X.left`
は自明にアフィン、`φ`(`toSchemeDiagramOver` 自身の遷移射)は
`toSchemeDiagram_isAffineHom`(`FieldLimit.lean`)から。

★**sorry 無し**。標準3公理のみ。

★★**残る前提(未確認・記録)**: 同じ被覆補題群は `∀ i, CompactSpace
(D.obj i)`・`∀ i, QuasiSeparatedSpace (D.obj i)` も要求するが、
`(extDiagram X).obj R = pullback X.hom (D R).hom` の `CompactSpace` は
`infer_instance` では自動的に付かない(`(D R).hom` はアフィンだが `X.hom`
は一般には何の制約も無いため、`X.left` 自身の準コンパクト性を要求する
可能性が高い)。双曲曲線は常に有限型(⟹ qcqs)なので原文の意図とは合致するが、
`corrHypInstance3` の `Space := Over BaseK` は一般のスキームを許すため、
`Lemma 4.1` を厳密に述べるには `Space` を qcqs スキームに絞るか、
`lemma_4_1` の統合先で `X` に有限型の仮定を追加する設計変更が要る——
`corrhyp-goal.md` に記録。 -/
theorem extDiagram_map_isAffineHom (X : Over BaseK) {R S : (FgSubalgebra ℚ ℝ)ᵒᵖ} (h : R ⟶ S) :
    IsAffineHom ((extDiagram X).map h) := by
  simp only [extDiagram]
  have hi1 : IsAffineHom (𝟙 X.left) := inferInstance
  have hi2 : IsAffineHom (toSchemeDiagramOver.map h).left := by
    show IsAffineHom ((toSchemeDiagram ℚ ℝ).map h)
    infer_instance
  exact MorphismProperty.pullbackMap hi1 hi2 (by simp)
    (by simpa using (toSchemeDiagramOver.map h).w.symm)

/-!
## 準コンパクトなスキームは有限アフィン開被覆を持つ

`Lemma 4.1` の被覆の組み立て(`corrhyp-goal.md` §4)の最初の一手——
`AffineTransitionLimit.lean` の `Scheme.exists_isOpenCover_and_isAffine`
に渡す入力(有限アフィン開被覆)を用意する。`Scheme.affineCover`
(mathlib既存、任意のスキームがアフィンスキームによる開被覆を持つ)の
台の開集合が `CompactSpace` により有限部分被覆へ落ちることを、一般
位相の `IsCompact.elim_finite_subcover` で示す——CorrHyp 固有の内容を
一切使わない、任意のスキームに対する一般的な事実。 -/

/-- **準コンパクトなスキームは(`affineCover` の中に)有限アフィン開被覆を
持つ**——`isAffineOpen_opensRange`(アフィンからの開埋め込みの像はアフィン
開)と一般位相の有限部分被覆の存在(`IsCompact.elim_finite_subcover`)を
合成する。

★**sorry 無し**。標準3公理のみ。CorrHyp に依存しない一般的な事実。 -/
theorem Scheme.exists_finite_affineOpenCover (X : Scheme) [CompactSpace X] :
    ∃ s : Finset X.affineCover.I₀,
      TopologicalSpace.IsOpenCover (fun i : s => (X.affineCover.f i).opensRange) ∧
      ∀ i : s, IsAffineOpen (X.affineCover.f i).opensRange := by
  have hcov : ⋃ i, Set.range (X.affineCover.f i) = Set.univ := X.affineCover.iUnion_range
  have hopen : ∀ i, IsOpen (Set.range (X.affineCover.f i)) :=
    fun i => (X.affineCover.f i).isOpenEmbedding.isOpen_range
  obtain ⟨s, hs⟩ := (isCompact_univ (X := X)).elim_finite_subcover
    (fun i => Set.range (X.affineCover.f i)) hopen (by rw [hcov])
  refine ⟨s, ?_, fun i => isAffineOpen_opensRange (X.affineCover.f i)⟩
  show (⨆ i : s, (X.affineCover.f i).opensRange) = ⊤
  apply TopologicalSpace.Opens.ext
  rw [TopologicalSpace.Opens.coe_iSup, TopologicalSpace.Opens.coe_top]
  apply Set.eq_univ_of_univ_subset
  intro x _
  have hx : x ∈ (⋃ i ∈ s, Set.range (X.affineCover.f i)) := hs (Set.mem_univ x)
  simp only [Set.mem_iUnion] at hx ⊢
  obtain ⟨i, hi, hxi⟩ := hx
  exact ⟨⟨i, hi⟩, hxi⟩

/-!
## `Ext X` の有限アフィン開被覆を有限段階へ降ろす

`Scheme.exists_finite_affineOpenCover`(`Ext X` は `CompactSpace` なので
有限アフィン開被覆を持つ)と `isLimit_extCone`(`Ext X` の台は `extDiagram X`
の極限)を `AffineTransitionLimit.lean` の `Scheme.exists_isOpenCover_and_
isAffine` に渡し、この被覆がある有限段階 `R`(`X.left ×_{BaseK} Spec R`)の
有限アフィン開被覆から来ることを示す——`Lemma 4.1` の構成的降下で最も
中心的な一手。 -/

/-- `(toSchemeDiagramOver.obj R).hom`(`Spec R → BaseK`)はアフィン射。 -/
theorem toSchemeDiagramOver_hom_isAffineHom (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    IsAffineHom (toSchemeDiagramOver.obj R).hom := by
  show IsAffineHom (Spec.map (CommRingCat.ofHom (algebraMap ℚ R.unop.1)))
  infer_instance

/-- `(toSchemeDiagramOver.obj R).left`(有限段階 `R` の台、`Spec R`)は
アフィン——「1アフィン片の降下」で `arrowIsoSpecΓOfIsAffine`(mathlib)を
使うための前提。 -/
theorem toSchemeDiagramOver_obj_isAffine (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    IsAffine (toSchemeDiagramOver.obj R).left := by
  show IsAffine (Spec (CommRingCat.of R.unop.1)); infer_instance

/-- `extDiagram X` の各段は `CompactSpace`——`(D R).hom` がアフィンなので
その base change `pullback.fst` もアフィン、よって `X.left` が
`CompactSpace` なら各段も `CompactSpace`。

★**sorry 無し**。標準3公理のみ。 -/
theorem extDiagram_obj_compactSpace (X : Over BaseK) [CompactSpace X.left]
    (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) : CompactSpace ((extDiagram X).obj R) := by
  show CompactSpace (Limits.pullback X.hom (toSchemeDiagramOver.obj R).hom : Scheme)
  have : IsAffineHom (Limits.pullback.fst X.hom (toSchemeDiagramOver.obj R).hom) :=
    MorphismProperty.pullback_fst X.hom (toSchemeDiagramOver.obj R).hom
      (toSchemeDiagramOver_hom_isAffineHom R)
  exact QuasiCompact.compactSpace_of_compactSpace
    (Limits.pullback.fst X.hom (toSchemeDiagramOver.obj R).hom)

/-- `extDiagram X` の各段は `QuasiSeparatedSpace`(`extDiagram_obj_
compactSpace` と同じ議論)。

★**sorry 無し**。標準3公理のみ。 -/
theorem extDiagram_obj_quasiSeparatedSpace (X : Over BaseK) [QuasiSeparatedSpace X.left]
    (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) : QuasiSeparatedSpace ((extDiagram X).obj R) := by
  show QuasiSeparatedSpace (Limits.pullback X.hom (toSchemeDiagramOver.obj R).hom : Scheme)
  have : IsAffineHom (Limits.pullback.fst X.hom (toSchemeDiagramOver.obj R).hom) :=
    MorphismProperty.pullback_fst X.hom (toSchemeDiagramOver.obj R).hom
      (toSchemeDiagramOver_hom_isAffineHom R)
  exact quasiSeparatedSpace_of_quasiSeparated
    (Limits.pullback.fst X.hom (toSchemeDiagramOver.obj R).hom)

/-- **`Ext X` の有限アフィン開被覆はある有限段階から来る**——`Lemma 4.1`
の構成的降下の核心部品。`ExtF.obj X = (Ext X)`(qcqs、`X` が qcqs のとき)
の有限アフィン開被覆(`Scheme.exists_finite_affineOpenCover`)を
`isLimit_extCone`(`Ext X` の台は `extDiagram X` の極限)経由で
`Scheme.exists_isOpenCover_and_isAffine` に渡し、`X.left ×_{BaseK} Spec R`
(ある有限生成 `k`-部分環 `R`)の有限アフィン開被覆に降ろす——降ろした先の
各片 `V j` が、元の `Ext X` 側の片(`(ExtF.obj X).left.affineCover.f (j:s)`
の像)の**引き戻しそのもの**であることも保持する(`hVprop j`)——`c.α`
(相関の有限エタール射)をこの片へ制限する際に使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_extDiagram_finite_affine_descent (X : Over BaseK) [CompactSpace X.left]
    [QuasiSeparatedSpace X.left] :
    ∃ (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) (s : Finset (ExtF.obj X).left.affineCover.I₀)
      (t : Finset s) (V : t → ((extDiagram X).obj R).Opens),
      TopologicalSpace.IsOpenCover V ∧ ∀ j : t, IsAffineOpen (V j) ∧
        (((ExtF.obj X).left).affineCover.f (j : s)).opensRange = (extCone X).π.app R ⁻¹ᵁ V j := by
  haveI hcompact : CompactSpace (ExtF.obj X).left := ExtF_compactSpace X
  haveI hqs : QuasiSeparatedSpace (ExtF.obj X).left := ExtF_quasiSeparatedSpace X
  obtain ⟨s, hcov, haff⟩ := Scheme.exists_finite_affineOpenCover (ExtF.obj X).left
  haveI hRcompact : ∀ R, CompactSpace ((extDiagram X).obj R) := extDiagram_obj_compactSpace X
  haveI hRqs : ∀ R, QuasiSeparatedSpace ((extDiagram X).obj R) :=
    extDiagram_obj_quasiSeparatedSpace X
  haveI hRaff : ∀ {R S : (FgSubalgebra ℚ ℝ)ᵒᵖ} (f : R ⟶ S), IsAffineHom ((extDiagram X).map f) :=
    fun f => extDiagram_map_isAffineHom X f
  obtain ⟨i, t, V, hVcov, hVprop⟩ := Scheme.exists_isOpenCover_and_isAffine
    (extDiagram X) (extCone X) (isLimit_extCone X)
    (fun j : s => (((ExtF.obj X).left).affineCover.f j).opensRange) hcov (fun j => haff j)
  exact ⟨i, s, t, V, hVcov, hVprop⟩

/-- **スキームレベルの `Etale` から環レベルの `Algebra.Etale` への橋渡し**——
`Lemma 4.1` の構成的降下の最後の接続部品。`[Etale α]` なスキーム射を
アフィン開 `U`(target 側)へ制限すると、誘導される環準同型
`α.appLE U (α ⁻¹ᵁ U) le_rfl` による代数構造で `Algebra.Etale` が成り立つ
——`α ⁻¹ᵁ U` がアフィンなのは `[IsFinite α]`(→ `IsAffineHom`)から
(`IsAffineHom.isAffine_preimage`)。これで各アフィン片の上で
`exists_finite_standardEtaleCover`(`FieldLimit.lean`)が使える。

★**sorry 無し**。標準3公理のみ。 -/
theorem Etale.algebraEtale_appLE {C Y : Scheme} (α : C ⟶ Y) [IsFinite α] [Etale α]
    (U : Y.Opens) (hU : IsAffineOpen U) :
    letI : Algebra Γ(Y, U) Γ(C, α ⁻¹ᵁ U) :=
      (Scheme.Hom.appLE α U (α ⁻¹ᵁ U) le_rfl).hom.toAlgebra
    Algebra.Etale Γ(Y, U) Γ(C, α ⁻¹ᵁ U) := by
  have hV : IsAffineOpen (α ⁻¹ᵁ U) := IsAffineHom.isAffine_preimage U hU
  exact Etale.etale_appLE α hU hV le_rfl

/-! ## `Ext X` の台を、有限段階でのbase changeとして特定する

`Lemma 4.1` の「1アフィン片の降下」に要る最後の鍵——`Ext X`(=`extCone X`
の頂点)の有限段階 `R` への射影 `(extConePi X).app R` が、**まさに
`Spec K → Spec R` に沿った base change の射影そのもの**であることを示す。
これにより `V ⊆ (extDiagram X).obj R` がアフィンなら、`Γ(preimage of V,
preimage of V)` が `pullbackSpecIso`(mathlib既存)経由で `Γ(V,V) ⊗[R] K`
と書けるようになる——`standardEtalePairRingBaseChange`(`FieldLimit.lean`)
を実際に適用するための土台。 -/

/-- `pullback.congrHom` の `.hom` と `pullback.fst`/`pullback.snd` の関係
(`.inv` 版の `congrHom_inv` から)——mathlib に `.hom` 版が無かったので補う。 -/
theorem pullback_congrHom_hom_fst {C : Type*} [CategoryTheory.Category C]
    [CategoryTheory.Limits.HasPullbacks C] {X Y Z : C} {f₁ f₂ : X ⟶ Z} {g₁ g₂ : Y ⟶ Z}
    (h₁ : f₁ = f₂) (h₂ : g₁ = g₂) :
    (Limits.pullback.congrHom h₁ h₂).hom ≫ Limits.pullback.fst f₂ g₂ =
      Limits.pullback.fst f₁ g₁ := by
  rw [← Iso.eq_inv_comp, Limits.pullback.congrHom_inv, Limits.pullback.lift_fst, Category.comp_id]

/-- `pullback_congrHom_hom_fst` の snd 版。 -/
theorem pullback_congrHom_hom_snd {C : Type*} [CategoryTheory.Category C]
    [CategoryTheory.Limits.HasPullbacks C] {X Y Z : C} {f₁ f₂ : X ⟶ Z} {g₁ g₂ : Y ⟶ Z}
    (h₁ : f₁ = f₂) (h₂ : g₁ = g₂) :
    (Limits.pullback.congrHom h₁ h₂).hom ≫ Limits.pullback.snd f₂ g₂ =
      Limits.pullback.snd f₁ g₁ := by
  rw [← Iso.eq_inv_comp, Limits.pullback.congrHom_inv, Limits.pullback.lift_snd, Category.comp_id]

/-- `Spec K → Spec R`(有限段階 `R` への遷移射)——`specKCone ℚ ℝ` の
射影そのもの。 -/
noncomputable def phiR (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) : specK ⟶ (toSchemeDiagramOver.obj R).left :=
  (specKCone ℚ ℝ).π.app R

/-- `phiR R` は `toBaseK`(`Spec K → Spec ℚ`)を `Spec R` 経由で分解する
——`specKConeOver`(`Over BaseK` の中のcone)の `.w` そのもの、任意の `R`
で成り立つ(`isLimit_specKConeOver` の証明が `R0` という特定の代表元
だけで使っていたのと同じ事実の一般形)。 -/
theorem phiR_comp (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    phiR R ≫ (toSchemeDiagramOver.obj R).hom = toBaseK :=
  (specKConeOver.π.app R).w

/-- **`Ext X` の台は、有限段階 `R` のファイバー積 `(extDiagram X).obj R`
を `Spec K → Spec R` に沿って base change したものに同型**——
`pullbackLeftPullbackSndIso`(mathlib、pullback の pasting)を
`phiR_comp` で `toBaseK` に付け替えるだけ。 -/
noncomputable def extConeIso (X : Over BaseK) (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    Limits.pullback (pullback.snd X.hom (toSchemeDiagramOver.obj R).hom) (phiR R) ≅
      Limits.pullback X.hom toBaseK :=
  (pullbackLeftPullbackSndIso X.hom (toSchemeDiagramOver.obj R).hom (phiR R)) ≪≫
    pullback.congrHom rfl (phiR_comp R)

/-- **`extConePi X .app R`(`Ext X` の `R` への射影)は、`extConeIso` の
下で「`P_R` への射影」そのもの**——`Lemma 4.1` の「1アフィン片の降下」
構成の核となる事実。`extConePi_app_fst`/`_snd`(定義そのもの)を
`pullback.hom_ext` で照合するだけだが、`m`/`n` の型を `set` で明示的に
`pullback X.hom toBaseK ⟶ pullback X.hom (toSchemeDiagramOver.obj R).hom`
に固定しておかないと `Functor.const`/`extDiagram.obj` の非簡約な展開で
`pullback.hom_ext` が「instances 透明度で型が合わない」を起こす
(`tools/lean-idioms.md` 第22項の教訓、`have` ではなく `set` で型注釈つきの
束縛をするのが鍵——`have` は非 `Prop` の値を消してしまうので後段で `rw`
できなくなる、という追加の教訓も得た)。

★**sorry 無し**。標準3公理のみ。 -/
theorem extConePi_app_eq (X : Over BaseK) (R : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    (extConePi X).app R = (extConeIso X R).inv ≫
      Limits.pullback.fst (pullback.snd X.hom (toSchemeDiagramOver.obj R).hom) (phiR R) := by
  set m : Limits.pullback X.hom toBaseK ⟶ Limits.pullback X.hom (toSchemeDiagramOver.obj R).hom :=
    (extConePi X).app R with hm
  set n : Limits.pullback X.hom toBaseK ⟶ Limits.pullback X.hom (toSchemeDiagramOver.obj R).hom :=
    (extConeIso X R).inv ≫
      Limits.pullback.fst (pullback.snd X.hom (toSchemeDiagramOver.obj R).hom) (phiR R) with hn
  have hm' : m ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R).hom = pullback.fst X.hom toBaseK := by
    rw [hm]; exact extConePi_app_fst X R
  have hm'' : m ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R).hom =
      pullback.snd X.hom toBaseK ≫ phiR R := by
    rw [hm]; exact extConePi_app_snd X R
  have hn' : n ≫ pullback.fst X.hom (toSchemeDiagramOver.obj R).hom = pullback.fst X.hom toBaseK := by
    rw [hn, Category.assoc]
    have h1 : (extConeIso X R).hom ≫ pullback.fst X.hom toBaseK =
        Limits.pullback.fst (pullback.snd X.hom (toSchemeDiagramOver.obj R).hom) (phiR R) ≫
          pullback.fst X.hom (toSchemeDiagramOver.obj R).hom := by
      show ((pullbackLeftPullbackSndIso X.hom (toSchemeDiagramOver.obj R).hom (phiR R)) ≪≫
          pullback.congrHom rfl (phiR_comp R)).hom ≫ pullback.fst X.hom toBaseK = _
      rw [Iso.trans_hom, Category.assoc, pullback_congrHom_hom_fst]
      exact pullbackLeftPullbackSndIso_hom_fst X.hom (toSchemeDiagramOver.obj R).hom (phiR R)
    rw [← h1, ← Category.assoc, Iso.inv_hom_id, Category.id_comp]
  have hn'' : n ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R).hom =
      pullback.snd X.hom toBaseK ≫ phiR R := by
    rw [hn, Category.assoc]
    have hcond : Limits.pullback.fst (pullback.snd X.hom (toSchemeDiagramOver.obj R).hom) (phiR R) ≫
        pullback.snd X.hom (toSchemeDiagramOver.obj R).hom =
        Limits.pullback.snd (pullback.snd X.hom (toSchemeDiagramOver.obj R).hom) (phiR R) ≫ phiR R :=
      Limits.pullback.condition
    rw [hcond]
    have h2 : (extConeIso X R).hom ≫ pullback.snd X.hom toBaseK =
        Limits.pullback.snd (pullback.snd X.hom (toSchemeDiagramOver.obj R).hom) (phiR R) := by
      show ((pullbackLeftPullbackSndIso X.hom (toSchemeDiagramOver.obj R).hom (phiR R)) ≪≫
          pullback.congrHom rfl (phiR_comp R)).hom ≫ pullback.snd X.hom toBaseK = _
      rw [Iso.trans_hom, Category.assoc, pullback_congrHom_hom_snd]
      exact pullbackLeftPullbackSndIso_hom_snd X.hom (toSchemeDiagramOver.obj R).hom (phiR R)
    rw [← h2, ← Category.assoc, ← Category.assoc, Iso.inv_hom_id, Category.id_comp]
  have heq : m = n := by
    apply pullback.hom_ext
    · rw [hm', hn']
    · rw [hm'', hn'']
  rw [hm, hn]; exact heq

/-! ## `Ext X` の `X.left` 由来のアフィン片は `Spec(Γ(U,U) ⊗[ℚ] ℝ)`
——generic flatness を要らなくする戦略

`X.left` の(**`R` に依らない、`X` 自身の**)アフィン開被覆 `{U_i}` を
そのまま `Ext X = X.left ×_k Spec K` へ base change すれば、各片は
`Spec(Γ(U_i,U_i) ⊗[ℚ] ℝ)` になる——`P_R`(有限段階のファイバー積)の
**任意の**アフィン開`V_j`から出発する(`exists_extDiagram_finite_affine_
descent`の戦略)と、後段で `Γ(V_j,V_j)` が `R` 上平坦とは限らないため
generic flatness(EGA IV、mathlibに現状無い)が要る——`X.left`自身の
アフィン開から出発すればテンソルが**常に `ℚ`(体)上**になり、
`Γ(U_i,U_i) ⊗[ℚ] -` は自動的に完全関手(平坦)になるので、この問題が
構造的に発生しない。`corrhyp-goal.md` §4 参照。 -/

/-- `X.left` の affine open `U` 上の切断環を `ℚ`-代数にする標準的な
環準同型(`BaseK = Spec ℚ` から `U.ι ≫ X.hom : U ⟶ BaseK` を
`Spec.preimage` で環準同型に直したもの)。 -/
noncomputable def pieceRingHom (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U) :
    CommRingCat.of ℚ ⟶ Γ(X.left, U) :=
  Spec.preimage (hU.isoSpec.inv ≫ (U.ι ≫ X.hom))

/-- `pieceRingHom` の定義方程式——`hU.isoSpec.hom` で移せば `U.ι ≫ X.hom`
と一致する。 -/
theorem pieceRingHom_spec (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U) :
    hU.isoSpec.hom ≫ Spec.map (pieceRingHom X U hU) = U.ι ≫ X.hom := by
  unfold pieceRingHom
  have heq : Spec.map (Spec.preimage (hU.isoSpec.inv ≫ (U.ι ≫ X.hom))) =
      hU.isoSpec.inv ≫ (U.ι ≫ X.hom) := Spec.map_preimage _
  rw [heq]
  have hgoal : hU.isoSpec.hom ≫ hU.isoSpec.inv ≫ (U.ι ≫ X.hom) = U.ι ≫ X.hom := by
    rw [Iso.hom_inv_id_assoc]
  exact hgoal

/-- `pullback (i.hom ≫ g) f ≅ pullback g f`(`i` が同型なら、片方の脚を
同型で付け替えても pullback は変わらない)——mathlib に無かったので
`pullback.map_isIso`(両側 `IsIso`)経由で補う。 -/
noncomputable def pullbackHomIsoLeft {X Y Z W : Scheme} (i : X ≅ Y) (g : Y ⟶ Z) (f : W ⟶ Z) :
    (Limits.pullback (i.hom ≫ g) f : Scheme) ≅ (Limits.pullback g f : Scheme) :=
  asIso (pullback.map (i.hom ≫ g) f g f i.hom (𝟙 _) (𝟙 _) (by simp) (by simp))

/-- `pieceRingHom` による `Γ(X.left,U)` の `ℚ`-代数構造。 -/
noncomputable def pieceAlgebra (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U) :
    Algebra ℚ Γ(X.left, U) :=
  (pieceRingHom X U hU).hom.toAlgebra

/-- **`Ext X` の `U`(`X.left` のアフィン開、`R` に依らない)上のアフィン片は
`Spec(Γ(U,U) ⊗[ℚ] ℝ)`**——`Lemma 4.1` の構成的降下で generic flatness を
回避する鍵となる構成(ファイル冒頭の節を見よ)。`pullbackRestrictIsoRestrict`
+`pullbackSymmetry`+`pullbackRightPullbackFstIso`(pullback の pasting)で
`(pullback.fst X.hom toBaseK) ⁻¹ᵁ U` を `pullback (U.ι ≫ X.hom) toBaseK`
と同一視し、`pieceRingHom_spec`+`pullbackHomIsoLeft` で `U.ι ≫ X.hom` を
`Spec.map (pieceRingHom X U hU)` に付け替えてから `pullbackSpecIso` を
適用する。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piecePullbackIso (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U) :
    letI := pieceAlgebra X U hU
    ((pullback.fst X.hom toBaseK) ⁻¹ᵁ U : Scheme) ≅
      Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] ℝ)) := by
  letI := pieceAlgebra X U hU
  calc ((pullback.fst X.hom toBaseK) ⁻¹ᵁ U : Scheme)
      ≅ (Limits.pullback (pullback.fst X.hom toBaseK) U.ι : Scheme) :=
        (pullbackRestrictIsoRestrict _ U).symm
    _ ≅ (Limits.pullback U.ι (pullback.fst X.hom toBaseK) : Scheme) := pullbackSymmetry _ _
    _ ≅ (Limits.pullback (U.ι ≫ X.hom) toBaseK : Scheme) := pullbackRightPullbackFstIso X.hom toBaseK U.ι
    _ ≅ (Limits.pullback (hU.isoSpec.hom ≫ Spec.map (pieceRingHom X U hU)) toBaseK : Scheme) :=
        (pieceRingHom_spec X U hU) ▸ Iso.refl _
    _ ≅ (Limits.pullback (Spec.map (pieceRingHom X U hU)) toBaseK : Scheme) :=
        pullbackHomIsoLeft hU.isoSpec (Spec.map (pieceRingHom X U hU)) toBaseK
    _ ≅ Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] ℝ)) := pullbackSpecIso ℚ Γ(X.left, U) ℝ

/-- `Ext X` の `U`(`X.left` のアフィン開)由来のアフィン片はアフィン
——`piecePullbackIso`の頂点が`IsIso`であることから。 -/
theorem piece_isAffineOpen (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U) :
    IsAffineOpen (pullback.fst X.hom toBaseK ⁻¹ᵁ U : ((ExtF.obj X).left).Opens) := by
  show IsAffine (pullback.fst X.hom toBaseK ⁻¹ᵁ U : Scheme)
  haveI : IsIso (piecePullbackIso X U hU).hom := (piecePullbackIso X U hU).isIso_hom
  exact IsAffine.of_isIso (piecePullbackIso X U hU).hom

/-- **`Ext X` へ有限エタールに写る `C` を、`X.left` のアフィン開 `U`
由来のアフィン片へ制限すると、環レベルの `Algebra.Etale` になる**——
`Etale.algebraEtale_appLE` を `piece_isAffineOpen` と組み合わせるだけ。
`Lemma 4.1` の「1アフィン片の降下」で `exists_finite_standardEtaleCover`
(`FieldLimit.lean`)を呼び出す直前の、スキーム→環の最後の接続点。
残る一手(未着手): `Γ((ExtF.obj X).left, V)` を `piecePullbackIso` 経由
で `Γ(U,U) ⊗[ℚ] ℝ`(`exists_fg_subalgebra_tensor_standardEtalePair_
baseChange` が読む形)へ変換すること——`Scheme.ΓSpecIso` を介した
環レベルの同型の輸送が要る。

★**sorry 無し**。標準3公理のみ。 -/
theorem piece_algebraEtale (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI V := (pullback.fst X.hom toBaseK ⁻¹ᵁ U : ((ExtF.obj X).left).Opens)
    letI : Algebra Γ((ExtF.obj X).left, V) Γ(C, α ⁻¹ᵁ V) :=
      (Scheme.Hom.appLE α V (α ⁻¹ᵁ V) le_rfl).hom.toAlgebra
    Algebra.Etale Γ((ExtF.obj X).left, V) Γ(C, α ⁻¹ᵁ V) :=
  Etale.algebraEtale_appLE α _ (piece_isAffineOpen X U hU)

/-- **`Algebra.Etale` は底環の同型に沿って輸送できる**——`RingHom.Etale.
respectsIso`(mathlib、pre-compose で保たれる)+`RingHom.etale_
algebraMap`(mathlib、`Algebra.Etale` との往復)を組み合わせるだけ。
`piece_algebraEtale`(環が`Γ((ExtF.obj X).left,V)`のまま)を
`pieceRingEquiv`が与える同型で`Γ(U,U)⊗[ℚ]ℝ`へ書き換えるのに使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem algebraEtale_transport {R R' S : Type} [CommRing R] [CommRing R'] [CommRing S]
    [Algebra R S] (e : R' ≃+* R) (h : Algebra.Etale R S) :
    letI : Algebra R' S := ((algebraMap R S).comp e.toRingHom).toAlgebra
    Algebra.Etale R' S := by
  letI : Algebra R' S := ((algebraMap R S).comp e.toRingHom).toAlgebra
  have h1 : (algebraMap R S).Etale := RingHom.etale_algebraMap.mpr h
  have h2 : ((algebraMap R S).comp e.toRingHom).Etale :=
    RingHom.Etale.respectsIso.2 (algebraMap R S) e h1
  have h3 : (algebraMap R' S) = (algebraMap R S).comp e.toRingHom := rfl
  rw [← RingHom.etale_algebraMap, h3]
  exact h2

/-- **`Ext X` の `U`(`X.left`のアフィン開)由来のアフィン片は
`Γ(U,U) ⊗[ℚ] ℝ` に(環として)同型**——`piecePullbackIso`(スキームの
同型)を `Scheme.Opens.topIso`(`Γ(↥U,⊤) ≅ Γ(X,U)`)・`Scheme.Γ.mapIso`
(スキーム圏の同型を`Γ`関手で送る)・`Scheme.ΓSpecIso`(`Γ(Spec R,⊤)≅R`)
で繋いで環レベルの同型に変換したもの。`piece_algebraEtale_tensor`の
核となる部品。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def pieceRingEquiv (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U) :
    letI := pieceAlgebra X U hU
    Γ((ExtF.obj X).left, pullback.fst X.hom toBaseK ⁻¹ᵁ U) ≃+* Γ(X.left, U) ⊗[ℚ] ℝ := by
  letI := pieceAlgebra X U hU
  set V := (pullback.fst X.hom toBaseK ⁻¹ᵁ U : ((ExtF.obj X).left).Opens) with hVdef
  have e1 : Γ((ExtF.obj X).left, V) ≃+* Γ((V : Scheme), ⊤) :=
    (Scheme.Opens.topIso V).symm.commRingCatIsoToRingEquiv
  have e2 : Γ((V : Scheme), ⊤) ≃+* Γ(Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] ℝ)), ⊤) :=
    (Scheme.Γ.mapIso (piecePullbackIso X U hU).symm.op).commRingCatIsoToRingEquiv
  have e3 : Γ(Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] ℝ)), ⊤) ≃+* Γ(X.left, U) ⊗[ℚ] ℝ :=
    (Scheme.ΓSpecIso (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] ℝ))).commRingCatIsoToRingEquiv
  exact e1.trans (e2.trans e3)

/-- **`Lemma 4.1`「1アフィン片の降下」の、スキーム→環の橋渡しの完成形**
——`C` を `Ext X` の `X.left` 由来のアフィン片へ制限すると、
`exists_finite_standardEtaleCover`・`exists_fg_subalgebra_tensor_
standardEtalePair_baseChange`(`FieldLimit.lean`)が直接読める形
(`Γ(U,U) ⊗[ℚ] ℝ` 上の `Algebra.Etale`)になる。`piece_algebraEtale`
(環はまだ `Γ((ExtF.obj X).left,V)`)を `pieceRingEquiv`+
`algebraEtale_transport` で輸送しただけ。

★**sorry 無し**。標準3公理のみ。 -/
theorem piece_algebraEtale_tensor (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
        (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
        (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
    Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) := by
  letI := pieceAlgebra X U hU
  letI : Algebra Γ((ExtF.obj X).left, pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    (Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.toAlgebra
  exact algebraEtale_transport (pieceRingEquiv X U hU).symm (piece_algebraEtale X U hU C α)

/-! ## `Γ(C,piece)`全体を一度に`R`レベルへ降ろす(設計の再考、2026-09-04)

`corrHypGlueData`の`Z i`(`descendPiece`)は個別の標準的エタール元
`f i`ごとに独立した`R_i`から`ℝ`へbase changeしたものであり、`R`
レベルの候補片そのもの(`Spec(P₀.Ring)`)は使い捨てられている。`C`は
`X`の`extDiagram`のような`R`段階近似の塔を持たないため、個別の`R_i`
を`⊔`で合流させる(`piece_descends_iso_promote`/`_family_promote`、
既存)よりも、`Γ(C,piece)`**全体**(`piece_algebraEtale_tensor`により
`(A⊗ℝ)`上`Etale`=`FinitePresentation`)を、`Algebra.Presentation.
ofFinitePresentation`(mathlib、有限表示から明示的な生成元・関係式を
取り出す)+`exists_fg_subalgebra_tensor_mvPolynomial_finset`(既存、
関係式の係数を単一の共通`R`へ降ろす)で**一度に**降ろす方が筋が良い
かもしれない、という代替案を検証している。 -/

open scoped Classical in
/-- **`Γ(C,piece)`の有限表示が持つ関係式を単一の共通`R`へ降ろす**——
`Algebra.Presentation.ofFinitePresentation`が与える関係式の族
(`Fin m`個)に`exists_fg_subalgebra_tensor_mvPolynomial_finset`を
適用するだけ。`Γ(C,piece)`全体を一度に`R`レベルへ降ろす計画の
中核部品——次の一手は、この`R`のもとで`MvPolynomial(Fin n)(A⊗R.1)
⧸ Ideal.span{降ろした関係式}`という`R`レベルの代数`S_0`を構成し、
`S_0 ⊗_{A⊗R.1}(A⊗ℝ) ≅ Γ(C,piece)`を示すこと(`Algebra.TensorProduct.
quotientTensorEquiv`(mathlib)+`Algebra.Presentation.baseChange`
(mathlib)を組み合わせる見込み、まだ未着手)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def pieceAlgebra_relation_descend_R (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
        (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
        (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
    haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      piece_algebraEtale_tensor X U hU C α
    haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      inferInstance
    FgSubalgebra ℚ ℝ := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  exact (exists_fg_subalgebra_tensor_mvPolynomial_finset (Γ(X.left, U))
    (Finset.image (Algebra.Presentation.ofFinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))).relation Finset.univ)).choose

open scoped Classical in
/-- `pieceAlgebra_relation_descend_R`が与える`R`のもとで、`Γ(C,piece)`
の関係式`k`番目(`P.relation k`)を実際に`R`レベルへ降ろした多項式
——`exists_fg_subalgebra_tensor_mvPolynomial_finset`の`choose_spec`
から`k`番目の成分を取り出すだけ。 -/
noncomputable def pieceAlgebra_relation_descend_q₀ (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
        (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
        (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
    haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      piece_algebraEtale_tensor X U hU C α
    haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      inferInstance
    Fin (Algebra.Presentation.ofFinitePresentationRels (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))) →
    MvPolynomial (Fin (Algebra.Presentation.ofFinitePresentationVars (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))))
      (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  intro k
  exact ((exists_fg_subalgebra_tensor_mvPolynomial_finset (Γ(X.left, U))
    (Finset.image (Algebra.Presentation.ofFinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))).relation Finset.univ)).choose_spec
    ((Algebra.Presentation.ofFinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))).relation k)
    (Finset.mem_image_of_mem _ (Finset.mem_univ k))).choose

open scoped Classical in
/-- **`Γ(C,piece)`の`R`レベルモデル`S_0`**——`pieceAlgebra_relation_
descend_q₀`が与える降ろした関係式による`MvPolynomial`商。`Γ(C,piece)`
全体を一度に`R`レベルへ降ろす計画の目標そのもの——次の一手は
`S_0 ⊗_{A⊗R.1}(A⊗ℝ) ≅ Γ(C,piece)`を示すこと(`MvPolynomial.
algebraTensorAlgEquiv`(mathlib、`A⊗MvPolynomial σ R ≃ₐ[A] MvPolynomial
σ A`)+`Algebra.TensorProduct.quotientTensorEquiv`(mathlib)+`Algebra.
Generators.ker_eq_ker_aeval_val`(mathlib、`P.ker`が`aeval`の核である
こと)を組み合わせる見込み)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def pieceAlgebra_R_model (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
        (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
        (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
    haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      piece_algebraEtale_tensor X U hU C α
    haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      inferInstance
    Type := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  letI hCR : CommRing (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) :=
    inferInstance
  exact MvPolynomial (Fin (Algebra.Presentation.ofFinitePresentationVars (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))))
      (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) ⧸
    Ideal.span (Set.range (pieceAlgebra_relation_descend_q₀ X U hU C α))

set_option maxHeartbeats 1000000 in
open scoped Classical in
/-- **`pieceAlgebra_R_model`(`S_0`)の base change が`Γ(C,piece)`を
復元すること**——`quotient_mvPolynomial_baseChange`(`FieldLimit.lean`、
汎用)を`R := Γ(U,U)⊗[ℚ]R.1`、`A := Γ(U,U)⊗[ℚ]ℝ`へ specialize すると、
残る作業は「降ろした関係式`q₀ k`の像が元の関係式`P.relation k`に一致
すること」(`pieceAlgebra_relation_descend_q₀`の定義そのもの、`choose_
spec`で直接取り出せる)だけ。そこから`Ideal.map_span`で像の集合として
のイデアルの一致に落とし、`Algebra.Presentation.span_range_relation_
eq_ker`+`Algebra.Generators.ker_eq_ker_aeval_val`(いずれもmathlib)で
`P.ker`(`= aeval P.val`の核)まで書き換え、最後に`RingHom.
quotientKerEquivOfSurjective`(第一同型定理、mathlib)で`Γ(C,piece)`
そのものに一致させる。`Γ(C,piece)`を`R`レベルへ一度に降ろす計画の
「base changeで元へ戻ること」の完成——次の一手は、`S_0`を候補片
(`descendPiece`の`R`レベル版)として使い、`transitionElem`/`gdT`/
`cocycle`の`R`レベル版を組み立てること。

★**sorry 無し**。標準3公理のみ。 -/
theorem pieceAlgebra_R_model_baseChange (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
        (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
        (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
    haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      piece_algebraEtale_tensor X U hU C α
    haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      inferInstance
    letI hCR : CommRing (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) :=
      inferInstance
    letI : Algebra (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
        (Γ(X.left, U) ⊗[ℚ] ℝ) :=
      (Algebra.TensorProduct.map (AlgHom.id ℚ Γ(X.left, U))
        (Subalgebra.val (pieceAlgebra_relation_descend_R X U hU C α).1)).toRingHom.toAlgebra
    Nonempty (
      (MvPolynomial (Fin (Algebra.Presentation.ofFinitePresentationVars (Γ(X.left, U) ⊗[ℚ] ℝ)
          Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))))
          (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) ⧸
        Ideal.span (Set.range (pieceAlgebra_relation_descend_q₀ X U hU C α)))
        ⊗[Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1] (Γ(X.left, U) ⊗[ℚ] ℝ)
      ≃+* Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))) := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  letI hCR : CommRing (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) :=
    inferInstance
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
      (Γ(X.left, U) ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ Γ(X.left, U))
      (Subalgebra.val (pieceAlgebra_relation_descend_R X U hU C α).1)).toRingHom.toAlgebra
  set P := Algebra.Presentation.ofFinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ)
    Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) with hP
  refine ⟨(quotient_mvPolynomial_baseChange (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
    (Γ(X.left, U) ⊗[ℚ] ℝ) (Fin (Algebra.Presentation.ofFinitePresentationVars (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))))
    (Ideal.span (Set.range (pieceAlgebra_relation_descend_q₀ X U hU C α)))).trans ?_⟩
  rw [Ideal.map_span]
  have hq₀ : ∀ k, MvPolynomial.map (algebraMap (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
      (Γ(X.left, U) ⊗[ℚ] ℝ)) (pieceAlgebra_relation_descend_q₀ X U hU C α k) = P.relation k := by
    intro k
    exact ((exists_fg_subalgebra_tensor_mvPolynomial_finset (Γ(X.left, U))
      (Finset.image (Algebra.Presentation.ofFinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ)
        Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))).relation Finset.univ)).choose_spec
      (P.relation k) (Finset.mem_image_of_mem _ (Finset.mem_univ k))).choose_spec
  have himg : (MvPolynomial.map (algebraMap (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
      (Γ(X.left, U) ⊗[ℚ] ℝ))) '' (Set.range (pieceAlgebra_relation_descend_q₀ X U hU C α))
      = Set.range P.relation := by
    rw [← Set.range_comp]
    exact congrArg Set.range (funext hq₀)
  rw [himg, P.span_range_relation_eq_ker, P.ker_eq_ker_aeval_val]
  exact RingHom.quotientKerEquivOfSurjective P.aeval_val_surjective

/-- **`piecePullbackIso` の有限段階版**——`Ext X` の代わりに有限段階
`(extDiagram X).obj R'` を使うと、`U`(`X.left`のアフィン開)由来のアフィン
片は `Spec(Γ(U,U) ⊗[ℚ] R'.1)` になる。証明の骨格は `piecePullbackIso`
と同一(`toBaseK` を `(toSchemeDiagramOver.obj R').hom` に置き換えた
だけ)——`(toSchemeDiagramOver.obj R').hom` が最初から `Spec.map` の形
なので、`piecePullbackIso` で必要だった`arrowIsoSpecΓOfIsAffine`相当の
変換が不要になる分だけ簡潔。`Lemma 4.1`の「`c.C`(既存のスキーム)から
有限段階への射を`glueMorphisms`で直接構成する」戦略の核となる部品。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piecePullbackIsoStage (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (R' : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    letI := pieceAlgebra X U hU
    ((pullback.fst X.hom (toSchemeDiagramOver.obj R').hom ⁻¹ᵁ U : Scheme)) ≅
      Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] R'.unop.1)) := by
  letI := pieceAlgebra X U hU
  have step1 : ((pullback.fst X.hom (toSchemeDiagramOver.obj R').hom) ⁻¹ᵁ U : Scheme) ≅
      (Limits.pullback (Spec.map (pieceRingHom X U hU)) (toSchemeDiagramOver.obj R').hom : Scheme) := by
    calc ((pullback.fst X.hom (toSchemeDiagramOver.obj R').hom) ⁻¹ᵁ U : Scheme)
        ≅ (Limits.pullback (pullback.fst X.hom (toSchemeDiagramOver.obj R').hom) U.ι : Scheme) :=
          (pullbackRestrictIsoRestrict _ U).symm
      _ ≅ (Limits.pullback U.ι (pullback.fst X.hom (toSchemeDiagramOver.obj R').hom) : Scheme) :=
          pullbackSymmetry _ _
      _ ≅ (Limits.pullback (U.ι ≫ X.hom) (toSchemeDiagramOver.obj R').hom : Scheme) :=
          pullbackRightPullbackFstIso X.hom (toSchemeDiagramOver.obj R').hom U.ι
      _ ≅ (Limits.pullback (hU.isoSpec.hom ≫ Spec.map (pieceRingHom X U hU))
            (toSchemeDiagramOver.obj R').hom : Scheme) :=
          (pieceRingHom_spec X U hU) ▸ Iso.refl _
      _ ≅ (Limits.pullback (Spec.map (pieceRingHom X U hU)) (toSchemeDiagramOver.obj R').hom : Scheme) :=
          pullbackHomIsoLeft hU.isoSpec (Spec.map (pieceRingHom X U hU)) (toSchemeDiagramOver.obj R').hom
  refine step1.trans ?_
  show (Limits.pullback (Spec.map (pieceRingHom X U hU))
    (Spec.map (CommRingCat.ofHom (algebraMap ℚ R'.unop.1))) : Scheme) ≅ _
  exact pullbackSpecIso ℚ Γ(X.left, U) R'.unop.1

/-- **`Ext X` 側の`U`のアフィン片は、`R'`側の同じアフィン片の
`extConePi X .app R'` による preimage そのもの**——`piecePullbackIso`
(`Ext X`側)と`piecePullbackIsoStage`(`R'`側)が「同じ`U`」を指して
いることの位相的な裏付け。`extConePi_app_fst`(定義)+`Scheme.Hom.
comp_preimage`(preimageの合成則)だけで閉じる。`C`の局所片を
`R'`側へ降ろす際、`Ext X`側の記述と`R'`側の記述が同じ場所を指している
ことを保証するのに使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem piece_preimage_eq (X : Over BaseK) (U : X.left.Opens) (R' : (FgSubalgebra ℚ ℝ)ᵒᵖ) :
    (extConePi X).app R' ⁻¹ᵁ (pullback.fst X.hom (toSchemeDiagramOver.obj R').hom ⁻¹ᵁ U) =
      pullback.fst X.hom toBaseK ⁻¹ᵁ U := by
  rw [← Scheme.Hom.comp_preimage, extConePi_app_fst]
  rfl

/-! ## `StandardEtalePair.Ring` を実際のスキームの局所片として実現する
——作業単位3(GlueData組み立て)の第一歩

`FieldLimit.lean`の`exists_fg_subalgebra_tensor_standardEtale_elem`で
「有限段階の候補局所片`P₀.Ring`」を環として構成できるようになった。
ここでは`P₀.Ring`を実際の**スキーム**として実現し(`Spec P₀.Ring`)、
それが有限段階の空間`Spec R`上で étale であること、かつ`K`へbase change
すると元の(`Ext X`側の)局所片`Spec P.Ring`にちょうど一致することを示す
——`c'.C`をGlueDataから組み立てる際の各ピースの実体になる。 -/

/-- `StandardEtalePair.Ring`の`Spec`への持ち上げ——`Spec P.Ring ⟶ Spec R`
という自然な射。`c'.C`の候補となる局所片の実体。 -/
noncomputable def standardEtalePairSpecMap {R : Type} [CommRing R] (P : StandardEtalePair R) :
    Spec (CommRingCat.of P.Ring) ⟶ Spec (CommRingCat.of R) :=
  Spec.map (CommRingCat.ofHom (algebraMap R P.Ring))

/-- **`standardEtalePairSpecMap`は常にétale**——`StandardEtalePair.Ring`に
mathlibが自動で与える`Algebra.Etale R P.Ring`インスタンスを、
`HasRingHomProperty.Spec_iff`(`AlgebraicGeometry.Etale`が`RingHom.Etale`
の`Spec`への持ち上げとして特徴づけられること)+`RingHom.etale_algebraMap`
で運ぶだけ。

★**sorry 無し**。標準3公理のみ。 -/
instance standardEtalePairSpecMap_etale {R : Type} [CommRing R] (P : StandardEtalePair R) :
    Etale (standardEtalePairSpecMap P) := by
  show Etale (Spec.map (CommRingCat.ofHom (algebraMap R P.Ring)))
  rw [HasRingHomProperty.Spec_iff (P := @Etale)]
  exact RingHom.etale_algebraMap.mpr inferInstance

/-- **候補局所片`Spec P₀.Ring`を`K`へbase changeすると、実際に元の
（`K`段階の）局所片`Spec P.Ring`(`P := P₀.map(algebraMap R K)`)に一致
する**——`pullbackSpecIso`(pullbackの計算)と`standardEtalePairRingBaseChange`
(環側のbase change)を`Algebra.TensorProduct.comm`で順序を合わせてから
合成するだけ。`c'.C`の各ピースが、実際に`c.C`の対応するピースへbase
changeで戻ることを保証する、GlueData組み立ての核心部品。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def standardEtalePairPullbackIso {R K : Type} [CommRing R] [CommRing K] [Algebra R K]
    (P₀ : StandardEtalePair R) :
    pullback (standardEtalePairSpecMap P₀) (Spec.map (CommRingCat.ofHom (algebraMap R K))) ≅
      Spec (CommRingCat.of (P₀.map (algebraMap R K)).Ring) := by
  have e : TensorProduct R P₀.Ring K ≃+* (P₀.map (algebraMap R K)).Ring :=
    (Algebra.TensorProduct.comm R P₀.Ring K).toRingEquiv.trans
      (standardEtalePairRingBaseChange P₀).toRingEquiv
  exact (pullbackSpecIso R P₀.Ring K).trans (Scheme.Spec.mapIso e.symm.toCommRingCatIso.op)

/-! ## `C` 自体の局所片がアフィンであること——`piece_algebraEtale_tensor` の
土台を確認する

`piece_algebraEtale_tensor`は「`Γ(C, α⁻¹(piece))`が環として
`Algebra.Etale`である」ことを述べるが、これが`α⁻¹(piece)`自体の
スキーム構造を正しく捉えるには、この開集合が**アフィン**であること
(`Γ`が全データを持つこと)を確認する必要がある——`α`が有限(ゆえに
アフィン射)であることと、`piece_isAffineOpen`(`Ext X`側のアフィン片が
アフィン)から`IsAffineOpen.preimage`で直ちに従う。 -/

/-- **`α`が有限のとき、`Ext X`のアフィン片`piece`の`α`による逆像は
アフィン**——`IsFinite α → IsAffineHom α`(`IsFinite`の定義に含まれる)
と`piece_isAffineOpen`(`Ext X`側)から`IsAffineOpen.preimage`で従う。
`exists_finite_standardEtaleCover`(`FieldLimit.lean`)を`Γ(C, α⁻¹(piece))`
へ適用する前に、この開集合が実際に`Spec Γ(C, α⁻¹(piece))`そのもの
であることを保証する。

★**sorry 無し**。標準3公理のみ。 -/
theorem piece_preimage_isAffineOpen (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] :
    IsAffineOpen (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U) : C.Opens) :=
  (piece_isAffineOpen X U hU).preimage α

/-- **`C`の局所片`α⁻¹(piece)`は`Spec Γ(C, α⁻¹(piece))`そのものと同一視
できる**——`piece_preimage_isAffineOpen`の`isoSpec`。`exists_finite_
standardEtaleCover`が与える基本開被覆`D(f_i)`を、実際に`C`内の開集合
として解釈するための橋渡し。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piecePreimageIso (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] :
    (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U) : Scheme) ≅
      Spec (CommRingCat.of Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))) :=
  (piece_preimage_isAffineOpen X U hU C α).isoSpec

set_option maxHeartbeats 1000000 in
/-- **`Γ(C,piece)`の`R`レベル候補片そのもの**——`pieceAlgebra_R_model`
(`S_0`)を`Spec`に渡すだけ。`descendPiece`(`StandardEtalePair`経由、
`standardEtalePairSpecMap`のpullbackとして構成、まだℝレベル)とは
異なり、こちらは**正真正銘`R`レベル**のアフィンスキーム(`Spec(A⊗R.1)`
上有限型)である。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def descendPieceR (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] : Scheme := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  letI hCR : CommRing (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) :=
    inferInstance
  exact Spec (CommRingCat.of (MvPolynomial (Fin (Algebra.Presentation.ofFinitePresentationVars
      (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))))
      (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) ⧸
    Ideal.span (Set.range (pieceAlgebra_relation_descend_q₀ X U hU C α))))

set_option maxHeartbeats 1000000 in
/-- `descendPieceR`から候補片の`R`レベル底`Spec(A⊗R.1)`への構造射——
`algebraMap`をここで1回だけ名前を付けて確定させる。配管の注意:
これを独立した`def`にせず`descendPieceR_iso`の証明中に生の`algebraMap`
式を直接書くと、`Ideal.Quotient`商の`CommRing`インスタンスが場所ごとに
別々に(非`defeq`に)導出されてしまい`pullbackSpecIso`との単一化に
失敗する——`standardEtalePairSpecMap`と同じ「名前を1回だけ確定させる」
パターンで解消した(新しい失敗形、`tools/lean-idioms.md`と同系統)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def descendPieceR_toBase (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    descendPieceR X U hU C α ⟶
      Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)) := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  letI hCR : CommRing (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) :=
    inferInstance
  unfold descendPieceR
  exact Spec.map (CommRingCat.ofHom (algebraMap (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
    (MvPolynomial (Fin (Algebra.Presentation.ofFinitePresentationVars (Γ(X.left, U) ⊗[ℚ] ℝ)
        Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))))
        (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) ⧸
      Ideal.span (Set.range (pieceAlgebra_relation_descend_q₀ X U hU C α)))))

set_option maxHeartbeats 1000000 in
/-- `Spec(A⊗ℝ)`から候補片の`R`レベル底`Spec(A⊗R.1)`への構造射
(`algebraMap (A⊗R.1)(A⊗ℝ)`のSpec)——`descendPieceR_toBase`と対にして
`pullback`を取ると、`descendPieceR`が実際に`α⁻¹(piece)`のbase change先
であることが言える(`descendPieceR_iso`)。 -/
noncomputable def descendPieceR_reBaseMap (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] ℝ)) ⟶
      Spec (CommRingCat.of (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)) := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  letI hCR : CommRing (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) :=
    inferInstance
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
      (Γ(X.left, U) ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ Γ(X.left, U))
      (Subalgebra.val (pieceAlgebra_relation_descend_R X U hU C α).1)).toRingHom.toAlgebra
  exact Spec.map (CommRingCat.ofHom (algebraMap (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
    (Γ(X.left, U) ⊗[ℚ] ℝ)))

set_option maxHeartbeats 1000000 in
/-- **`descendPieceR`は実際に`α⁻¹(piece)`のbase change先である**——
`pullbackSpecIso`(mathlib、`pullback(Spec.map algebraMap R S)(Spec.map
algebraMap R T) ≅ Spec(S⊗[R]T)`)と`pieceAlgebra_R_model_baseChange`
(`S_0⊗[R](A⊗ℝ)≅Γ(C,piece)`)・`piecePreimageIso`(`Spec Γ(C,piece)≅
α⁻¹(piece)`)を合成するだけ。`standardEtalePairPullbackIso`と全く同じ
形の主張だが、こちらは**正真正銘`R`レベル**の`descendPieceR`について
成り立つ——`Γ(C,piece)`を`R`レベルへ一度に降ろす計画の、スキーム
レベルでの完成形。

★**sorry 無し**。標準3公理のみ。 -/
theorem descendPieceR_iso (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    Nonempty (pullback (descendPieceR_toBase X U hU C α) (descendPieceR_reBaseMap X U hU C α) ≅
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U) : Scheme)) := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  haveI : Algebra.FinitePresentation (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    inferInstance
  letI hCR : CommRing (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) :=
    inferInstance
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
      (Γ(X.left, U) ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ Γ(X.left, U))
      (Subalgebra.val (pieceAlgebra_relation_descend_R X U hU C α).1)).toRingHom.toAlgebra
  unfold descendPieceR_toBase descendPieceR_reBaseMap descendPieceR
  refine ⟨(pullbackSpecIso (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1)
    (MvPolynomial (Fin (Algebra.Presentation.ofFinitePresentationVars (Γ(X.left, U) ⊗[ℚ] ℝ)
        Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))))
        (Γ(X.left, U) ⊗[ℚ] (pieceAlgebra_relation_descend_R X U hU C α).1) ⧸
      Ideal.span (Set.range (pieceAlgebra_relation_descend_q₀ X U hU C α)))
    (Γ(X.left, U) ⊗[ℚ] ℝ)).trans ?_⟩
  refine (Scheme.Spec.mapIso (pieceAlgebra_R_model_baseChange X U hU C α).some.symm.toCommRingCatIso.op).trans ?_
  exact (piecePreimageIso X U hU C α).symm

/-! ## 「項目(d)の第二段」実データ接続への第一歩——`X.left`の異なる
アフィン開`U`・`V`(`V≤U`)にまたがる`c.C`側の片同士の制限写像

`2026-09-05(続き11)`の精査で判明した通り、ℝレベルでは`α⁻¹(U)`・`α⁻¹(V)`
は共に`c.C`自身の開集合なので、両者の重なりの貼り合わせデータは
`c.C`の層構造が既に自動的に与える(非自明な遷移同型は不要)——非自明に
なるのは、これを`descendPieceR`(Rレベル、抽象的なスキーム)側で
比較するときだけである。この制限写像`piece_restrict_hom`は、
`exists_mvPolynomial_quotient_specIso_descend`の`ψ`(2つの`Presentation`
の生成元同士の対応)を実際に構成するための第一部品——`Algebra.
Presentation.ofFinitePresentation`の生成元(`P.val`)をこの写像で送り、
overlap上の元として実現するのに使う。 -/

open CategoryTheory AlgebraicGeometry Limits in
/-- `V≤U`のとき、`c.C`側の対応する片も同じ向きの包含になること
(`α⁻¹`・`pullback.fst⁻¹`という2段の逆像がどちらも単調であることから
直ちに従う)。 -/
theorem piece_le_of_le (X : Over BaseK) (U V : X.left.Opens) (hV : V ≤ U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) :
    (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ V) : C.Opens) ≤
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U) : C.Opens) :=
  fun _ hx => hV hx

open CategoryTheory AlgebraicGeometry Limits in
/-- **`C`側の片`α⁻¹(piece)`の、`X.left`側の開集合の包含`V≤U`に沿った
制限写像**——標準的な層の制限(`C.presheaf.map (homOfLE ...).op`)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piece_restrict_hom (X : Over BaseK) (U V : X.left.Opens) (hV : V ≤ U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) :
    Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) ⟶ Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ V)) :=
  C.presheaf.map (homOfLE (piece_le_of_le X U V hV C α)).op

/-! ### `piece_restrict_hom`を基本開の被覆(`2026-09-05続き13`の簡略化)へ
特殊化する

「共通の細分`W`」問題は、`U_i`・`U_j`を任意のアフィン開のペアとして
扱うのではなく、単一のアフィン開`U`の基本開被覆`X.basicOpen f_i`で
覆う設計に切り替えれば、`Scheme.basicOpen_mul`(mathlib、既存)により
`W := X.basicOpen (f*g) = X.basicOpen f ⊓ X.basicOpen g`が自動的に
アフィンになり(基本開同士の交わりは再び基本開)、分離性の仮定が
一切不要になる——`corrHypGlueData`(既存のℝレベル版)が最初から
採用していた設計そのもの。 -/

open CategoryTheory AlgebraicGeometry Limits in
/-- **`piece_restrict_hom`の基本開特殊化(左側)**——`D(f)`側の片から
`D(f*g) = D(f)⊓D(g)`側の片への制限写像。`Scheme.basicOpen_mul`
(mathlib)+`inf_le_left`で`D(f*g) ≤ D(f)`を得るだけ。 -/
noncomputable def piece_restrict_hom_basicOpen_left (X : Over BaseK) (U : X.left.Opens) (f g : Γ(X.left, U))
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) :
    Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen f))) ⟶
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen (f * g)))) :=
  piece_restrict_hom X (X.left.basicOpen f) (X.left.basicOpen (f * g))
    (by rw [Scheme.basicOpen_mul]; exact inf_le_left) C α

open CategoryTheory AlgebraicGeometry Limits in
/-- **`piece_restrict_hom`の基本開特殊化(右側)**——`D(g)`側の片から
`D(f*g)`側の片への制限写像(左側の対)。 -/
noncomputable def piece_restrict_hom_basicOpen_right (X : Over BaseK) (U : X.left.Opens) (f g : Γ(X.left, U))
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) :
    Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen g))) ⟶
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen (f * g)))) :=
  piece_restrict_hom X (X.left.basicOpen g) (X.left.basicOpen (f * g))
    (by rw [Scheme.basicOpen_mul]; exact inf_le_right) C α

/-! ### `C`側の片`piece(D(f*g))`を`piece(D(f))`の基本開として実現する
(`2026-09-05続き16`)

`descendPieceR`(`R`レベル)を実際の局所化として貼り合わせるための
足がかり——`piece(D(f*g))`(`α⁻¹(pullback.fst⁻¹(D(f*g)))`)が、`C`の
なかで**`piece(D(f))`のある元による基本開そのもの**であることを示す。
`X.left`レベルの事実`X.basicOpen(f*g)=X.basicOpen(g|_{D(f)})`
(`Scheme.basicOpen_mul`+`Scheme.basicOpen_res`、`rw`の自動`rfl`閉じで
1行で証明できることを`2026-09-05続き15`で確認済み)を、`mathlib`の
`Scheme.preimage_basicOpen`(`f⁻¹ᵁ(Y.basicOpen r)=X.basicOpen(f.app r)`)
で`pullback.fst`・`α`の2段の逆像へ**別々の`have`として**押し出す
——`rw`を2回連鎖させると「motive is not type correct」(`(ExtF.obj X).left`
と`pullback X.hom toBaseK`が`rw`の構文一致チェックでは同一視されない)
になったため、中間結果を`(ExtF.obj X).left.basicOpen`の形で明示的に
型注釈した`have`として確定させ、`exact`(defeq判定)で個別に閉じる
という2段構成にして解消した(新しい失敗形、配管の一種)。 -/

open CategoryTheory AlgebraicGeometry Limits in
/-- **`piece(D(f*g))=piece(D(f)).basicOpen(...)`の局所化パラメータ**——
`g`を`D(f)`へ制限してから`pullback.fst`・`α`の2段で`C`側へ押し出した元。
`Γ(C,piece(D(f)))`の元として、次の`piece_basicOpen_mul_eq`のRHSに現れる。 -/
noncomputable def piece_basicOpen_localizationElem (X : Over BaseK) (U : X.left.Opens) (f g : Γ(X.left, U))
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) :
    Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen f))) :=
  α.app (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen f))
    ((pullback.fst X.hom toBaseK).app (X.left.basicOpen f)
      (X.left.presheaf.map (homOfLE (X.left.basicOpen_le f)).op g))

open CategoryTheory AlgebraicGeometry Limits in
/-- **`C`側の片`piece(D(f*g))`は`piece(D(f))`の基本開そのものである**——
`piece_basicOpen_localizationElem`を局所化パラメータとする。`R`レベルの
`descendPieceR`側で、`D(f*g)`の`S_0`を`D(f)`の`S_0`のこの元(の`R`レベル
持ち上げ)による`Localization.Away`として実際に構成するための、
スキームレベルでの根拠——`GlueData.f i j`が開埋め込みであることを
`Spec`の左では言えたが、`descendPieceR`(`Spec`の右、`R`レベルの抽象
モデル)側でも同じ開埋め込みが実現できることを、この事実が保証する
土台になる(まだ`R`レベルへの引き上げ自体は未着手)。

★**sorry 無し**。標準3公理のみ。 -/
theorem piece_basicOpen_mul_eq (X : Over BaseK) (U : X.left.Opens) (f g : Γ(X.left, U))
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) :
    (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen (f * g))) : C.Opens) =
      C.basicOpen (piece_basicOpen_localizationElem X U f g C α) := by
  have hcore : X.left.basicOpen (f * g) = X.left.basicOpen
      (X.left.presheaf.map (homOfLE (X.left.basicOpen_le f)).op g) := by
    rw [Scheme.basicOpen_mul, Scheme.basicOpen_res]
  have h1 : (pullback.fst X.hom toBaseK ⁻¹ᵁ (X.left.basicOpen (f * g)) : (ExtF.obj X).left.Opens) =
      (ExtF.obj X).left.basicOpen ((pullback.fst X.hom toBaseK).app (X.left.basicOpen f)
        (X.left.presheaf.map (homOfLE (X.left.basicOpen_le f)).op g)) := by
    rw [hcore]; exact Scheme.preimage_basicOpen _ _
  rw [h1]
  exact Scheme.preimage_basicOpen _ _

/-- **アフィン開`U`上の基本開`X.basicOpen f`は`Spec(Localization.Away f)`
そのものと同一視できる**——`mathlib`の`basicOpenIsoSpecAway`は`X := Spec R`
の場合限定だったので、一般のアフィン開`U`(`X`自体はアフィンでなくてよい)
版として一般化した。`hU.isLocalization_basicOpen`(`Γ(X,X.basicOpen f)`が
`Localization.Away f`の実現になっていること)+`(hU.basicOpen f).isoSpec`
を合成するだけ。`exists_finite_standardEtaleCover`が返す基本開被覆
`D(f_i)`(`C`内、`piece_preimage_isAffineOpen`で得たアフィン開の上)を、
実際に`Spec(Localization.Away f_i)`として認識するのに使う——CorrHyp
に依存しない一般的な事実。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def IsAffineOpen.basicOpenIsoSpecAway {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U)
    (f : Γ(X, U)) :
    (X.basicOpen f : Scheme) ≅ Spec (CommRingCat.of (Localization.Away f)) := by
  haveI := hU.isLocalization_basicOpen (f := f)
  have e : Γ(X, X.basicOpen f) ≃ₐ[Γ(X, U)] Localization.Away f :=
    IsLocalization.algEquiv (Submonoid.powers f) _ _
  exact (hU.basicOpen f).isoSpec.trans (Scheme.Spec.mapIso e.symm.toRingEquiv.toCommRingCatIso.op)

/-! ## 1ピース分の完全な連結——`standardEtalePairPullbackIso`・
`IsAffineOpen.basicOpenIsoSpecAway`・`exists_fg_subalgebra_tensor_
standardEtalePair_mapEq`の合成

前回の詰まり(`algebraMap (A⊗R.1)(A⊗ℝ)`と`letI`で導入した`Algebra`
インスタンスの構文不一致)は、`exists_fg_subalgebra_tensor_
standardEtalePair_mapEq`が返す等式を`letI`のスコープ内で`algebraMap`形へ
**defeqの代入だけで**(`▸`や`show`を使わず、ただの`:=`で)通せることに
気づいて解消した——`standardEtalePairPullbackIso`を**再証明せず**
そのまま呼び、変換した等式で`▸`するだけにしたのが鍵(前回は環の同型を
インラインで再導出しようとして`TensorProduct`の引数順を取り違えていた)。 -/

open scoped TensorProduct in
/-- **標準エタールな元`f`の局所化`Localization.Away f`は、ある有限段階の
候補スキーム(`standardEtalePairSpecMap P₀`)を`A⊗ℝ`へbase changeした
ものと同型**——`exists_fg_subalgebra_tensor_standardEtalePair_mapEq`
(環の降下)+`standardEtalePairPullbackIso`(base changeのスキーム版)を
`letI`のスコープ内でdefeqの代入だけで繋ぐ。`Algebra.IsStandardEtale`
から得られる`StandardEtalePresentation`(`Nonempty`)を`Classical.choice`
的に取り出して使う。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def onePieceSchemeIso {A B : Type} [CommRing A] [Algebra ℚ A]
    [CommRing B] [Algebra (A ⊗[ℚ] ℝ) B] (f : B)
    [Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)] :
    ∃ (R : FgSubalgebra ℚ ℝ) (P₀ : StandardEtalePair (A ⊗[ℚ] R.1)),
      letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
      Nonempty ((Spec (CommRingCat.of (Localization.Away f)) : Scheme) ≅
        pullback (standardEtalePairSpecMap P₀)
          (Spec.map (CommRingCat.ofHom (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ))))) := by
  obtain ⟨Pres⟩ :=
    ‹Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)›.nonempty_standardEtalePresentation
  obtain ⟨R, P₀, hP₀⟩ := exists_fg_subalgebra_tensor_standardEtalePair_mapEq A Pres.P
  letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
  have hP₀' : P₀.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ)) = Pres.P := hP₀
  refine ⟨R, P₀, ⟨?_⟩⟩
  have i2 : (pullback (standardEtalePairSpecMap P₀)
      (Spec.map (CommRingCat.ofHom (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ)))) : Scheme) ≅
      Spec (CommRingCat.of Pres.P.Ring) := hP₀' ▸ standardEtalePairPullbackIso P₀
  have i1 : (Spec (CommRingCat.of (Localization.Away f)) : Scheme) ≅ Spec (CommRingCat.of Pres.P.Ring) :=
    Scheme.Spec.mapIso Pres.equivRing.toRingEquiv.toCommRingCatIso.symm.op
  exact i1.trans i2.symm

open scoped TensorProduct in
/-- **作業単位1・3の核心の合流点——`C`(または任意のスキーム`X`)の
標準エタール局所片`X.basicOpen f`は、有限段階の候補局所片のbase change
にちょうど一致する**——`IsAffineOpen.basicOpenIsoSpecAway`(位相・
アフィン性)と`onePieceSchemeIso`(環→スキームの降下)を合成するだけ。
`Lemma 4.1`の「1アフィン片の降下」において、`exists_finite_
standardEtaleCover`が返す各`f_i`をこの補題に適用すれば、`C`の対応する
開集合が有限段階の候補スキームのbase changeと一致することが**直接**
得られる——GlueDataの各ピースの構成に使う核心部品。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piece_descends_iso {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    (f : Γ(X, U)) [Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)] :
    ∃ (R : FgSubalgebra ℚ ℝ) (P₀ : StandardEtalePair (A ⊗[ℚ] R.1)),
      letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
      Nonempty ((X.basicOpen f : Scheme) ≅
        pullback (standardEtalePairSpecMap P₀)
          (Spec.map (CommRingCat.ofHom (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ))))) := by
  obtain ⟨R, P₀, ⟨e⟩⟩ := onePieceSchemeIso (A := A) f
  exact ⟨R, P₀, ⟨(IsAffineOpen.basicOpenIsoSpecAway hU f).trans e⟩⟩

/-! ## `R`レベルの貼り合わせへ向けた第一歩——`piece_descends_iso`の
降下先`R`を、より粗い共通段階`R'`へ昇格する

`Lemma 4.1`の真の`k`(=`ℚ`)への降下には、`piece_descends_iso`が与える
`R`レベルの候補片`Spec(P₀.Ring)`自体を、`R`レベル(`FgSubalgebra ℚ ℝ`
の圏)で貼り合わせる必要がある(2026-09-04、`corrhyp-goal.md`に記録
した訂正)。異なる添字`i`ごとに得られる`R_i`は一般に異なるため、
比較の前に共通の粗い段階`R'`(`R_i ⊔ R_j`のような)へ揃える必要がある
——ここではその**第一歩**(1つの`piece_descends_iso`の結果を任意の
より粗い`R'`へ昇格する)を用意する。 -/

/-- **`piece_descends_iso`が与える降下先`R`を、より粗い共通段階`R'`
(`R≤R'`)へ昇格しても、`X.basicOpen f`との比較同型が保たれる**——
`standardEtalePairPullbackIso`(`R`レベルの候補片は、より粗い段階への
base changeでも`.map`の`Spec`として実現できる、既存の一般的事実)を
2回+`exists_fg_subalgebra_tensor_standardEtalePair_promote`(既存、
`FieldLimit.lean`)を組み合わせるだけ——**新しい`mathlib`の降下定理
(cofiltered極限からの降下)は不要**だった。`R`レベルの複数添字を
共通`R'`へ合流させる際の核心部品。

★**sorry 無し**。標準3公理のみ。 -/
@[reducible]
noncomputable def piece_descends_iso_promote {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    (f : Γ(X, U)) [Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)]
    (R' : FgSubalgebra ℚ ℝ) (hle : (piece_descends_iso (A := A) hU f).choose ≤ R') :
    letI : Algebra (A ⊗[ℚ] R'.1) (A ⊗[ℚ] ℝ) :=
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom.toAlgebra
    Nonempty ((X.basicOpen f : Scheme) ≅
      pullback (standardEtalePairSpecMap
          ((piece_descends_iso (A := A) hU f).choose_spec.choose.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle)).toRingHom))
        (Spec.map (CommRingCat.ofHom (algebraMap (A ⊗[ℚ] R'.1) (A ⊗[ℚ] ℝ))))) := by
  set R := (piece_descends_iso (A := A) hU f).choose with hR
  set P₀ := (piece_descends_iso (A := A) hU f).choose_spec.choose with hP₀
  letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
  letI : Algebra (A ⊗[ℚ] R'.1) (A ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom.toAlgebra
  obtain ⟨e⟩ := (piece_descends_iso (A := A) hU f).choose_spec.choose_spec
  refine ⟨e.trans ((standardEtalePairPullbackIso P₀).trans ?_)⟩
  have hpromote : (P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle)).toRingHom).map
      (algebraMap (A ⊗[ℚ] R'.1) (A ⊗[ℚ] ℝ)) = P₀.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ)) :=
    exists_fg_subalgebra_tensor_standardEtalePair_promote A R R' hle P₀
      (P₀.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ))) rfl
  rw [← hpromote]
  exact (standardEtalePairPullbackIso
    (P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle)).toRingHom)).symm

/-- `piece_descends_iso`の降下先`R`を、`IsStandardEtale`を明示的な証明
として受け取る形で取り出す——`Finset`の各元ごとに別々の証明を渡して
族を組み立てるのに使う(`descendPieceOfProof`と同じ理由)。 -/
noncomputable def piece_descends_iso_R_of_proof {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    (f : Γ(X, U)) (hf : Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)) :
    FgSubalgebra ℚ ℝ :=
  letI := hf
  (piece_descends_iso (A := A) hU f).choose

open scoped Classical in
/-- **有限個の`piece_descends_iso`の降下先`R_i`の共通上界`R'`**——
`exists_fgSubalgebra_upperBound`(既存、`FieldLimit.lean`)を
`piece_descends_iso_R_of_proof`の族へ適用するだけ。`R`レベルの複数
添字を単一の共通段階へ合流させる際の第一歩。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piece_descends_iso_R_upperBound {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    {ι : Type} (t : Finset ι) (f : ι → Γ(X, U))
    (hf : ∀ i ∈ t, Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away (f i))) :
    FgSubalgebra ℚ ℝ :=
  (exists_fgSubalgebra_upperBound t
    (fun i => if h : i ∈ t then piece_descends_iso_R_of_proof (A := A) hU (f i) (hf i h)
      else ⟨⊥, Subalgebra.fg_bot⟩)).choose

open scoped Classical in
/-- `piece_descends_iso_R_upperBound`が実際に各`R_i`の上界であること
——`exists_fgSubalgebra_upperBound`の`.choose_spec`をそのまま呼ぶだけ。
`piece_descends_iso_promote`の`hle`引数として使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem piece_descends_iso_R_upperBound_spec {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    {ι : Type} (t : Finset ι) (f : ι → Γ(X, U))
    (hf : ∀ i ∈ t, Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away (f i)))
    (i : ι) (hi : i ∈ t) :
    piece_descends_iso_R_of_proof (A := A) hU (f i) (hf i hi) ≤
      piece_descends_iso_R_upperBound (A := A) hU t f hf := by
  unfold piece_descends_iso_R_upperBound
  have h := (exists_fgSubalgebra_upperBound t
    (fun i => if h : i ∈ t then piece_descends_iso_R_of_proof (A := A) hU (f i) (hf i h)
      else ⟨⊥, Subalgebra.fg_bot⟩)).choose_spec i hi
  simpa [hi] using h

open scoped Classical in
/-- **有限個の`piece_descends_iso`の族を、単一の共通段階`R'`
(`piece_descends_iso_R_upperBound`)へ実際に揃える**——各`i∈t`ごとに
`letI := hf i hi`の下で`piece_descends_iso_promote`を
`piece_descends_iso_R_upperBound_spec`が与える`hle`とともに適用する
だけ(`piece_descends_iso_R_of_proof`と`piece_descends_iso_promote`
内部の`(piece_descends_iso ...).choose`は`rfl`で一致することを確認
済み)。`R`レベルの複数添字を単一の共通段階`R'`へ合流させる最終形——
これで`R`レベルの候補片`Spec(P₀_i'.Ring)`(`P₀_i'`は`R'`レベルへ
昇格した後の`StandardEtalePair`)がすべて**同じ`R'`**の上に揃った。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piece_descends_iso_family_promote {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    {ι : Type} (t : Finset ι) (f : ι → Γ(X, U))
    (hf : ∀ i ∈ t, Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away (f i)))
    (i : ι) (hi : i ∈ t) :
    letI := hf i hi
    letI : Algebra (A ⊗[ℚ] (piece_descends_iso_R_upperBound (A := A) hU t f hf).1) (A ⊗[ℚ] ℝ) :=
      (Algebra.TensorProduct.map (AlgHom.id ℚ A)
        (Subalgebra.val (piece_descends_iso_R_upperBound (A := A) hU t f hf).1)).toRingHom.toAlgebra
    Nonempty ((X.basicOpen (f i) : Scheme) ≅
      pullback (standardEtalePairSpecMap
          ((piece_descends_iso (A := A) hU (f i)).choose_spec.choose.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A)
              (Subalgebra.inclusion (piece_descends_iso_R_upperBound_spec (A := A) hU t f hf i hi))).toRingHom))
        (Spec.map (CommRingCat.ofHom (algebraMap (A ⊗[ℚ] (piece_descends_iso_R_upperBound (A := A) hU t f hf).1)
          (A ⊗[ℚ] ℝ))))) := by
  letI := hf i hi
  exact piece_descends_iso_promote (A := A) hU (f i) (piece_descends_iso_R_upperBound (A := A) hU t f hf)
    (piece_descends_iso_R_upperBound_spec (A := A) hU t f hf i hi)

/- ★★次の一手(未着手): `piece_descends_iso_family_promote`で全員が
揃った`R'`レベルの候補片`Spec(P₀_i'.Ring)`同士の**重なり**(`transitionElem`
/`gdT`/`cocycle`の`R`レベル版)を構築する——ただし`ℝ`レベルの
`transitionElem`は`X.presheaf.map`(層の制限写像)を使っており、`R`
レベルでは対応する「層」が`extDiagram X`の段階`R'`におけるアンビエント
スキームの層(`(extDiagram X).obj (op R')`)になる見込みで、この対応を
精密に詰める必要がある(まだ未着手)。`corrhyp-goal.md`に記録。 -/

/-! ## GlueDataの遷移射——`piece_descends_iso`を重なりへ制限する

前回、遷移射(`D(f_i)`と`D(f_j)`の重なりの比較)を`Localization`の
局所化の局所化という抽象的な環同型で構成しようとして「捩れ」の壁に
再度当たったが、`piece_descends_iso`(既に完成)を重なり
`X.basicOpen(f_i·f_j)`へ**制限**するだけで済むと判明した——抽象的な
環レベルの独立検証は一切不要。ここではその制限を実際に構成する。 -/

/-- スキームの同型`e : X ≅ Z`について、`e.hom`による像は`e.inv`による
逆像に一致する——`e.hom ⁻¹ᵁ`の単射性(`preimage_image_eq`+`e.hom≫e.inv
=𝟙`)から。同型の下での開集合の像を「逆側の射の逆像」として計算する
ための一般的な橋渡し(CorrHyp非依存)。

★**sorry 無し**。標準3公理のみ。 -/
theorem Scheme.hom_image_iso_eq_inv_preimage {X Z : Scheme} (e : X ≅ Z) (W : X.Opens) :
    e.hom ''ᵁ W = e.inv ⁻¹ᵁ W := by
  have h1 : e.hom ⁻¹ᵁ (e.hom ''ᵁ W) = W := e.hom.preimage_image_eq W
  have h2 : e.hom ⁻¹ᵁ (e.inv ⁻¹ᵁ W) = W := by
    rw [← Scheme.Hom.comp_preimage, Iso.hom_inv_id]; simp
  have hinj : Function.Injective (fun (V : Z.Opens) => e.hom ⁻¹ᵁ V) := by
    intro V₁ V₂ hV
    have := congrArg (fun (V : X.Opens) => e.inv ⁻¹ᵁ V) hV
    simp only [← Scheme.Hom.comp_preimage] at this
    simpa using this
  exact hinj (h1.trans h2.symm)

/-- **開埋め込みの像への制限`isoImage`は部分開集合の包含と両立する**——
`V ≤ W`のとき、`X.homOfLE(V≤W)`で`W`から`V`へ制限してから`isoImage W`で
`Y`側へ渡すのと、先に`isoImage V`で`Y`側へ渡してから像の包含
`Y.homOfLE(p''ᵁV≤p''ᵁW)`で制限するのとが一致する——`isoImage`の自然性。
`t_fac`・`cocycle`(`transitionElemIso`が積の分解と両立すること)の核心
部品。証明は`(p''ᵁW).ι`で両辺をキャンセルするだけ(`isoImage_hom_ι`+
`homOfLE_ι`)。CorrHyp非依存の一般的事実。

★**sorry 無し**。標準3公理のみ。 -/
theorem Scheme.Hom.isoImage_naturality {X Y : Scheme} (p : X ⟶ Y) [IsOpenImmersion p]
    {V W : X.Opens} (h : V ≤ W) :
    X.homOfLE h ≫ (p.isoImage W).hom = (p.isoImage V).hom ≫ Y.homOfLE (Scheme.Hom.image_mono p h) := by
  rw [← cancel_mono (p ''ᵁ W).ι]
  simp [Scheme.Hom.isoImage_hom_ι, Scheme.homOfLE_ι, Scheme.homOfLE_ι_assoc]

/-- `Scheme.Hom.isoImage_naturality`の逆向き(`.inv`版)——`isoImage`の
逆同型についての同じ自然性。`isoImage_naturality`から
`cancel_epi (p.isoImage V).hom`で導く。

★**sorry 無し**。標準3公理のみ。 -/
theorem Scheme.Hom.isoImage_inv_naturality {X Y : Scheme} (p : X ⟶ Y) [IsOpenImmersion p]
    {V W : X.Opens} (h : V ≤ W) :
    (p.isoImage V).inv ≫ X.homOfLE h = Y.homOfLE (Scheme.Hom.image_mono p h) ≫ (p.isoImage W).inv := by
  rw [← cancel_epi (p.isoImage V).hom, Iso.hom_inv_id_assoc, ← Category.assoc,
    ← Scheme.Hom.isoImage_naturality p h, Category.assoc, Iso.hom_inv_id, Category.comp_id]

attribute [reassoc] Scheme.Hom.isoImage_naturality Scheme.Hom.isoImage_inv_naturality

/-- **`eqToIso`は`homOfLE`と両立する**(`.inv`版)——`A=A'`・`B=B'`のとき、
「`B'≤A'`に沿って制限してから`A'`を`A`へ同一視する」のと「先に`B'`を`B`へ
同一視してから`B≤A`に沿って制限する」のが一致する。`eqToIso`の引数は
`Eq`(`Prop`)なので証明の中身は問わない(証明無関係性)——`subst`+`simp`で
機械的に示せる。`t_fac`(`transitionElemIso`の自然性)で、精密な等式
(`transitionElem_basicOpen_eq`)を経由した2つの異なる`homOfLE`を
橋渡しするために使う。CorrHyp非依存の一般的事実。

★**sorry 無し**。標準3公理のみ。 -/
theorem eqToIso_homOfLE_comm {X : Scheme} {A A' B B' : X.Opens} (hA : A = A') (hB : B = B')
    (h : B ≤ A) (h' : B' ≤ A')
    (eA : (A : Scheme) = (A' : Scheme)) (eB : (B : Scheme) = (B' : Scheme)) :
    X.homOfLE h' ≫ (eqToIso eA).inv = (eqToIso eB).inv ≫ X.homOfLE h := by
  subst hA; subst hB
  simp

/-- `eqToIso_homOfLE_comm`の`.hom`版(向きが逆)。

★**sorry 無し**。標準3公理のみ。 -/
theorem eqToIso_homOfLE_comm' {X : Scheme} {A A' B B' : X.Opens} (hA : A = A') (hB : B = B')
    (h : B ≤ A) (h' : B' ≤ A')
    (eA : (A : Scheme) = (A' : Scheme)) (eB : (B : Scheme) = (B' : Scheme)) :
    X.homOfLE h ≫ (eqToIso eA).hom = (eqToIso eB).hom ≫ X.homOfLE h' := by
  subst hA; subst hB
  simp

/-- **重なり`X.basicOpen(f₁·f₂)`を、`X.basicOpen f₁`とその「候補片」`Z`
との任意の同型`e`(`piece_descends_iso`/`onePieceSchemeIso`が与える)の
下で`Z`の基本開集合として実現する**——`W`(`X.basicOpen f₁`自身のスキーム
の中の対応する開集合)を経由して、(a) `(X.basicOpen f₁).ι`で`X`側の
`X.basicOpen(f₁·f₂)`に戻ることと、(b) `e.hom`で`Z`側の基本開集合
`Z.basicOpen s`に写ることの両方を保証する。`s`は具体的な切断として
与えられる。

2つの標準エタール片`D(f_i)`・`D(f_j)`をそれぞれの候補片`Z_i`・`Z_j`へ
写したとき、両方とも`X.basicOpen(f_i·f_j)`という**同じ`X`内の開集合**
を経由するので、遷移射(`Z_i`の対応する開集合と`Z_j`の対応する開集合の
同型)は`e_i.hom`・`e_j.hom`のこの記述を合成するだけで得られる——
抽象的な環レベルの独立検証(`Localization`の局所化の局所化の比較)は
不要になる。GlueDataの遷移射構成の核心。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_transitionOpen_eq_basicOpen {X : Scheme} {U : X.Opens} (f₁ f₂ : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    ∃ s : Γ(Z, ⊤),
      (X.basicOpen f₁).ι ''ᵁ ((X.basicOpen f₁).toScheme.basicOpen
          ((X.basicOpen f₁).topIso.inv
            (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op f₂))) = X.basicOpen (f₁ * f₂) ∧
      e.hom ''ᵁ ((X.basicOpen f₁).toScheme.basicOpen
          ((X.basicOpen f₁).topIso.inv
            (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op f₂))) = Z.basicOpen s := by
  refine ⟨e.inv.app ⊤ ((X.basicOpen f₁).topIso.inv
    (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op f₂)), ?_, ?_⟩
  · rw [Scheme.Opens.ι_image_basicOpen_topIso_inv, Scheme.basicOpen_res, Scheme.basicOpen_mul]
  · rw [Scheme.hom_image_iso_eq_inv_preimage, Scheme.preimage_basicOpen]
    rfl

/-- **`Spec Γ(X,U)`側の基本開被覆(環レベル、`PrimeSpectrum.basicOpen`)は
`X`側の基本開被覆(`X.basicOpen`)にそのまま運べる**——`hU.fromSpec`
(`Spec Γ(X,U) ⟶ X`、値域がちょうど`U`)の`preimage_basicOpen`との可換性
+`image_preimage_eq_opensRange_inf`(開埋め込みの像と逆像の関係)を
組み合わせるだけ。`exists_finite_standardEtaleCover`(`FieldLimit.lean`)
が与える環レベルの被覆を、実際に`C`(または任意のスキーム)内の開被覆
として認識するための橋渡し——`Scheme.openCoverOfIsOpenCover`と組み合わ
せれば`C`自身の`OpenCover`が作れ、`Scheme.Cover.gluedCover`(cocycle
自動証明済み)を適用する道が開ける。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_scheme_basicOpen_cover_of_ring {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U)
    {ι : Type} (t : Finset ι) (f : ι → Γ(X, U))
    (hcov : (⨆ i ∈ t, PrimeSpectrum.basicOpen (f i)) = ⊤) :
    (⨆ i ∈ t, X.basicOpen (f i)) = U := by
  have hpre : hU.fromSpec ⁻¹ᵁ (⨆ i ∈ t, X.basicOpen (f i)) = ⊤ := by
    simp_rw [Scheme.Hom.preimage_iSup, hU.fromSpec_preimage_basicOpen]
    exact hcov
  have hrange : hU.fromSpec.opensRange = U := by
    apply TopologicalSpace.Opens.ext; exact hU.range_fromSpec
  have himg := Scheme.Hom.image_preimage_eq_opensRange_inf hU.fromSpec (⨆ i ∈ t, X.basicOpen (f i))
  rw [hpre, Scheme.Hom.image_top_eq_opensRange, hrange] at himg
  apply le_antisymm
  · exact iSup₂_le (fun i _ => X.basicOpen_le (f i))
  · exact himg.le.trans inf_le_right

/-- **`exists_transitionOpen_eq_basicOpen`の帰結を、開集合の等式ではなく
実際のスキーム同型として取り出す**——重なり`X.basicOpen(f₁·f₂)`が、
候補片`Z`の対応する基本開集合`Z.basicOpen s`と**直接同型**であることを
示す。`Scheme.Hom.isoImage`(開埋め込みの、その開集合への制限が同型で
あること)+`eqToIso`(`exists_transitionOpen_eq_basicOpen`が与える開集合
の等式)を合成するだけ。

2つの標準エタール片`i`・`j`にこれを適用すると、どちらも同じ
`X.basicOpen(f_i·f_j)`と同型になるので、それらを合成すれば
`Z_i.basicOpen(s_i) ≅ Z_j.basicOpen(s_j)`という**遷移射そのもの**が
直接得られる——GlueDataの`t`フィールドの構成材料。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_transitionIso {X : Scheme} {U : X.Opens} (f₁ f₂ : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    ∃ s : Γ(Z, ⊤), Nonempty
      ((X.basicOpen (f₁ * f₂) : Scheme) ≅ (Z.basicOpen s : Scheme)) := by
  obtain ⟨s, h1, h2⟩ := exists_transitionOpen_eq_basicOpen f₁ f₂ e
  refine ⟨s, ⟨?_⟩⟩
  set W := (X.basicOpen f₁).toScheme.basicOpen
    ((X.basicOpen f₁).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op f₂)) with hW
  have iso1 : (W : Scheme) ≅ (X.basicOpen (f₁ * f₂) : Scheme) :=
    ((X.basicOpen f₁).ι.isoImage W).trans (eqToIso (by rw [h1]))
  have iso2 : (W : Scheme) ≅ (Z.basicOpen s : Scheme) :=
    (e.hom.isoImage W).trans (eqToIso (by rw [h2]))
  exact iso1.symm.trans iso2

/-- **`exists_transitionIso`は3重(以上)の重なりも同じ1回の適用でカバー
する**——`f₂`に任意の元(たとえば`f_j*f_k`という積)を渡せるので、GlueDataの
`t'`(3重の重なりでの整合性)構成に新しい数学的内容は不要と分かる。
`s.prod g`(有限個の元の積)を渡すだけ。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_transitionIso_finset {X : Scheme} {U : X.Opens} (f₁ : Γ(X, U))
    {κ : Type} (s : Finset κ) (g : κ → Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    ∃ s' : Γ(Z, ⊤), Nonempty
      ((X.basicOpen (f₁ * s.prod g) : Scheme) ≅ (Z.basicOpen s' : Scheme)) :=
  exists_transitionIso f₁ (s.prod g) e

/-- **遷移射の座標を`Classical.choice`を経由せず直接の式として定義する**
——`e.inv`・`topIso.inv`・制限写像という3つの環準同型の合成を`f₂`に
適用するだけ。`Classical.choice`の不透明性を経由しないので、乗法性等の
代数的性質を自由に示せるようになる(`t'`の構成に必須)。`GlueData`の族
`f`・`Z`・`e`を導入する前に、単体の`X`・`U`・`Z`・`e`について自己完結的に
定義しておく(`exists_transitionIso`と同じ流儀)——族の`f`・`Z`・`e`との
名前衝突を避けるため。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def transitionElem {X : Scheme} {U : X.Opens} (f₁ f₂ : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) : Γ(Z, ⊤) :=
  e.inv.app ⊤ ((X.basicOpen f₁).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op f₂))

/-- **`transitionElem`の乗法性**——3つの環準同型の合成が乗法的である
ことから`map_mul`を3回適用するだけ。`t'`(3重の重なり)の構成の核心。

★**sorry 無し**。標準3公理のみ。 -/
theorem transitionElem_mul {X : Scheme} {U : X.Opens} (f₁ f₂ f₃ : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    transitionElem f₁ (f₂ * f₃) e = transitionElem f₁ f₂ e * transitionElem f₁ f₃ e := by
  unfold transitionElem
  rw [map_mul, map_mul, map_mul]
  rfl

/-- **`transitionElem`が満たす精密な等式**——`exists_transitionOpen_eq_
basicOpen`の証明をそのまま`transitionElem`について直接示したもの(存在
命題を経由しない)。

★**sorry 無し**。標準3公理のみ。 -/
theorem transitionElem_basicOpen_eq {X : Scheme} {U : X.Opens} (f₁ f₂ : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    (X.basicOpen f₁).ι ''ᵁ ((X.basicOpen f₁).toScheme.basicOpen
        ((X.basicOpen f₁).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op f₂)))
      = X.basicOpen (f₁ * f₂) ∧
    e.hom ''ᵁ ((X.basicOpen f₁).toScheme.basicOpen
        ((X.basicOpen f₁).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op f₂)))
      = Z.basicOpen (transitionElem f₁ f₂ e) := by
  unfold transitionElem
  refine ⟨?_, ?_⟩
  · rw [Scheme.Opens.ι_image_basicOpen_topIso_inv, Scheme.basicOpen_res, Scheme.basicOpen_mul]
  · rw [Scheme.hom_image_iso_eq_inv_preimage, Scheme.preimage_basicOpen]
    rfl

/-- **`transitionElem`版の遷移同型**——`X.basicOpen(f₁·g) ≅ Z.basicOpen
(transitionElem f₁ g e)`。`exists_transitionIso`の`transitionElem`版。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def transitionElemIso {X : Scheme} {U : X.Opens} (f₁ g : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    (X.basicOpen (f₁ * g) : Scheme) ≅ (Z.basicOpen (transitionElem f₁ g e) : Scheme) := by
  have h1 := (transitionElem_basicOpen_eq f₁ g e).1
  have h2 := (transitionElem_basicOpen_eq f₁ g e).2
  set W := (X.basicOpen f₁).toScheme.basicOpen
    ((X.basicOpen f₁).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op g))
    with hW
  have iso1 : (W : Scheme) ≅ (X.basicOpen (f₁ * g) : Scheme) :=
    ((X.basicOpen f₁).ι.isoImage W).trans (eqToIso (by rw [h1]))
  have iso2 : (W : Scheme) ≅ (Z.basicOpen (transitionElem f₁ g e) : Scheme) :=
    (e.hom.isoImage W).trans (eqToIso (by rw [h2]))
  exact iso1.symm.trans iso2

/-- `eqToIso`を、ある開集合の族の等式(`Opens`のまま)を`congrArg`で
コード化されたスキームレベルの等式へ持ち上げたものへ適用すると、`.ι`
と可換——`subst`すれば両辺とも`Iso.refl`になるだけの軽い事実だが、
`transitionElemIso`(`eqToIso`を直接使う構成)の自然性を調べるのに
必要になる(`Scheme.isoOfEq_hom_ι`の`eqToIso`版に相当)。 -/
theorem eqToIso_congrArg_scheme_hom_ι {Z : Scheme} {V V' : Z.Opens} (hOpens : V = V') :
    (eqToIso (congrArg (fun (W : Z.Opens) => (W : Scheme)) hOpens)).hom ≫ V'.ι = V.ι := by
  subst hOpens
  simp

/-- **`transitionElemIso`の`.ι`との可換性**——`X.basicOpen(f₁*g)`から
`X.basicOpen f₁`への制限(`X.homOfLE`)を`e`で転送すると、`Z.basicOpen
(transitionElem f₁ g e)`から`Z`自身への包含(`gdF`そのもの)に一致する。
`iso1`(`X`側、`.ι`)・`iso2`(`Z`側、`e.hom`の`isoImage`)それぞれを
`Scheme.Hom.isoImage_hom_ι`+`eqToIso_congrArg_scheme_hom_ι`で開き、
`(X.basicOpen f₁).ι`が単射(open immersionでmono)であることによる
簡約(`cancel_mono`)で`iso1`側を`X.homOfLE`の形にまとめる。

`piecesGluedCoverVIso`(`V`成分比較同型)の`gdF`(`corrHypGlueData`の
`fst`)との自然性の核心部品——項目(b)のNatIso構成に向けた最後の
ピース。

★**sorry 無し**。標準3公理のみ。 -/
theorem transitionElemIso_hom_ι {X : Scheme} {U : X.Opens} (f₁ g : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    (transitionElemIso f₁ g e).hom ≫ (Z.basicOpen (transitionElem f₁ g e)).ι =
      X.homOfLE (show X.basicOpen (f₁ * g) ≤ X.basicOpen f₁ by
          rw [X.basicOpen_mul f₁ g]; exact inf_le_left) ≫ e.hom := by
  set hle : X.basicOpen (f₁ * g) ≤ X.basicOpen f₁ := by
    rw [X.basicOpen_mul f₁ g]; exact inf_le_left with hle_def
  unfold transitionElemIso
  simp only [Iso.trans_hom, Iso.symm_hom]
  set h1 := (transitionElem_basicOpen_eq f₁ g e).1
  set h2 := (transitionElem_basicOpen_eq f₁ g e).2
  set W := (X.basicOpen f₁).toScheme.basicOpen
    ((X.basicOpen f₁).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op g)) with hW
  set iso1 : (W : Scheme) ≅ (X.basicOpen (f₁ * g) : Scheme) :=
    ((X.basicOpen f₁).ι.isoImage W).trans (eqToIso (by rw [h1])) with hiso1
  have key1 : iso1.hom ≫ X.homOfLE hle = W.ι := by
    rw [← cancel_mono (X.basicOpen f₁).ι]
    rw [Category.assoc, Scheme.homOfLE_ι]
    show ((X.basicOpen f₁).ι.isoImage W).hom ≫ (eqToIso (congrArg (fun (V:X.Opens) => (V:Scheme)) h1)).hom ≫
      (X.basicOpen (f₁ * g)).ι = W.ι ≫ (X.basicOpen f₁).ι
    rw [eqToIso_congrArg_scheme_hom_ι h1]
    exact Scheme.Hom.isoImage_hom_ι _ _
  have key2 : iso1.inv ≫ W.ι = X.homOfLE hle := by
    rw [← key1, Iso.inv_hom_id_assoc]
  have key3 : (Scheme.Hom.isoImage e.hom W).hom ≫ (eqToIso (congrArg (fun (V:Z.Opens) => (V:Scheme)) h2)).hom ≫
      (Z.basicOpen (transitionElem f₁ g e)).ι = W.ι ≫ e.hom := by
    rw [eqToIso_congrArg_scheme_hom_ι h2]
    exact Scheme.Hom.isoImage_hom_ι _ _
  show iso1.inv ≫ (Scheme.Hom.isoImage e.hom W).hom ≫ (eqToIso (congrArg (fun (V:Z.Opens) => (V:Scheme)) h2)).hom ≫
    (Z.basicOpen (transitionElem f₁ g e)).ι = X.homOfLE hle ≫ e.hom
  rw [key3, ← key2, Category.assoc]

/-- `transitionElem`の座標を与える制限元は、積の分解と両立する
(`X.presheaf.map`の乗法性+`basicOpen_mul`)——`t_fac`で`transitionElemIso`
の自然性を組み立てるための土台となる開集合の包含。 -/
theorem transitionElem_restrict_mul_le {X : Scheme} {U : X.Opens} (f₁ g h : Γ(X, U)) :
    (X.basicOpen f₁ : Scheme).basicOpen
      ((ConcreteCategory.hom (X.basicOpen f₁).topIso.inv)
        ((ConcreteCategory.hom (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op)) (g * h)))
    ≤ (X.basicOpen f₁ : Scheme).basicOpen
      ((ConcreteCategory.hom (X.basicOpen f₁).topIso.inv)
        ((ConcreteCategory.hom (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op)) g)) := by
  rw [map_mul, map_mul, Scheme.basicOpen_mul]; exact inf_le_left

/-- `t_fac`の核心である`transitionElemIso_inv_naturality`を組み立てる
5つの部品のうち1つ目——`eqToIso_homOfLE_comm`を`transitionElem`の
`Z`側の精密な等式に適用しただけ。個別の`theorem`として切り出す理由:
本体の`transitionElemIso_inv_naturality`が1つの巨大な`have`ブロックとして
書かれていると、`X.presheaf`絡みの型検査が`whnf`のheartbeat上限に
達してしまう(配管の教訓、`tools/lean-idioms.md`参照)。 -/
theorem transitionElemIso_step1 {X : Scheme} {U : X.Opens} (f₁ g h : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    Z.homOfLE (show Z.basicOpen (transitionElem f₁ (g*h) e) ≤ Z.basicOpen (transitionElem f₁ g e) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_left)
      ≫ (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).2)).inv
    = (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
      Z.homOfLE (Scheme.Hom.image_mono e.hom (transitionElem_restrict_mul_le f₁ g h)) :=
  eqToIso_homOfLE_comm (transitionElem_basicOpen_eq f₁ g e).2 (transitionElem_basicOpen_eq f₁ (g*h) e).2
    (Scheme.Hom.image_mono e.hom (transitionElem_restrict_mul_le f₁ g h))
    (show Z.basicOpen (transitionElem f₁ (g*h) e) ≤ Z.basicOpen (transitionElem f₁ g e) by
      rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_left)
    (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).2)
    (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)

/-- 部品2つ目——`isoImage_inv_naturality`を`e.hom`に適用しただけ。 -/
theorem transitionElemIso_step2 {X : Scheme} {U : X.Opens} (f₁ g h : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    (Scheme.Hom.isoImage e.hom _).inv ≫ (X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le f₁ g h)
    = Z.homOfLE (Scheme.Hom.image_mono e.hom (transitionElem_restrict_mul_le f₁ g h)) ≫
        (Scheme.Hom.isoImage e.hom _).inv :=
  Scheme.Hom.isoImage_inv_naturality e.hom (transitionElem_restrict_mul_le f₁ g h)

/-- 部品3つ目——`isoImage_naturality`を`(X.basicOpen f₁).ι`に適用しただけ。 -/
theorem transitionElemIso_step3 {X : Scheme} {U : X.Opens} (f₁ g h : Γ(X, U)) :
    (X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le f₁ g h) ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom
    = (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le f₁ g h)) :=
  Scheme.Hom.isoImage_naturality (X.basicOpen f₁).ι (transitionElem_restrict_mul_le f₁ g h)

/-- 部品4つ目——`eqToIso_homOfLE_comm'`を`transitionElem`の`X`側の
精密な等式に適用しただけ。 -/
theorem transitionElemIso_step4 {X : Scheme} {U : X.Opens} (f₁ g h : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le f₁ g h)) ≫
      (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom
    = (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).1)).hom ≫
      X.homOfLE (show X.basicOpen (f₁*(g*h)) ≤ X.basicOpen (f₁*g) by
        rw [show f₁*(g*h) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left) :=
  eqToIso_homOfLE_comm' (transitionElem_basicOpen_eq f₁ g e).1 (transitionElem_basicOpen_eq f₁ (g*h) e).1
    (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le f₁ g h))
    (show X.basicOpen (f₁*(g*h)) ≤ X.basicOpen (f₁*g) by
      rw [show f₁*(g*h) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left)
    (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)
    (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).1)

/-- 部品4つ目と5つ目(`X.homOfLE≫ι = ι`という`homOfLE_ι`)を1つにまとめた
もの——`transitionElemIso_inv_naturality`本体の`calc`の最終段でこの形の
まま使う(`step4`単体を使い、あとから`homOfLE_ι`を別の`calc`段として
足すと、その段だけ`whnf`のheartbeat上限に達してしまうため)。 -/
theorem transitionElemIso_step45 {X : Scheme} {U : X.Opens} (f₁ g h : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le f₁ g h)) ≫
      (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
      (X.basicOpen (f₁*g)).ι
    = (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).1)).hom ≫
      (X.basicOpen (f₁*(g*h))).ι := by
  rw [← Category.assoc, transitionElemIso_step4 f₁ g h e, Category.assoc,
    Scheme.homOfLE_ι X (show X.basicOpen (f₁*(g*h)) ≤ X.basicOpen (f₁*g) by
      rw [show f₁*(g*h) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left)]

/-- **`transitionElemIso`の自然性——`t_fac`の核心**。`transitionElemIso
f₁ (g·h) e`の`.inv`を、`transitionElemIso f₁ g e`の`.inv`に、両側の
基本開集合の包含(`Z.basicOpen(t_{f₁(gh)})≤Z.basicOpen(t_{f₁g})`・
`X.basicOpen(f₁(gh))≤X.basicOpen(f₁g)`)を挟んで橋渡しできる——
`transitionElemIso f₁ (g·h) e`を`X.basicOpen(f₁g)`へ制限したものが、
`transitionElemIso f₁ g e`と一致するという主張。

**配管の教訓(新項目、`tools/lean-idioms.md`へ追記予定)**: この証明は
`unfold transitionElemIso`+`simp`で得られる巨大な`isoImage`/`eqToIso`
の塔に対して、`rw`(や`simp`単体のさらなる適用)で部品を差し込もうと
すると、**常に**「`instances`透明度で`X.presheaf`/`X.presheaf.obj`の
型が合わない」というエラーに当たる——`unfold`+`simp`で作られる
`Scheme.basicOpen`絡みの巨大な項に対して、`rw`・`simp`・`conv`は
(たとえ挿入する事実自体がその場で`have`として正しく型検査できても)
congruence motiveの構築に失敗する。**唯一有効だったのは、`rw`を一切
使わず、`calc`+`congrArg`+`(Category.assoc _ _ _).symm`という
**term-mode**の構成だけで結合子を組み立てること**——`congrArg`は
ゴールに対して`rw`のようなmotive探索(`kabstract`)を行わず、与えられた
関数と等式から直接新しい項を構成するので、この壁を経由しない。
さらに、各部品(`step1`〜`step5`)を`transitionElemIso_inv_naturality`
の中の`have`として書くと、それだけで`whnf`のheartbeat上限(400万でも
不足)に達する——**個別の独立した`theorem`として先に証明し切ってから、
`calc`の中では名前を参照するだけにする**ことで初めて軽くなった
(依存関係が閉じた小さな項を`theorem`として確定させると、以降の
`whnf`はその項を不透明な定数として扱えるため)。

★**sorry 無し**。標準3公理のみ。 -/
theorem transitionElemIso_inv_naturality {X : Scheme} {U : X.Opens} (f₁ g h : Γ(X, U))
    {Z : Scheme} (e : (X.basicOpen f₁ : Scheme) ≅ Z) :
    Z.homOfLE (show Z.basicOpen (transitionElem f₁ (g*h) e) ≤ Z.basicOpen (transitionElem f₁ g e) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_left)
      ≫ (transitionElemIso f₁ g e).inv
    = (transitionElemIso f₁ (g*h) e).inv ≫
      X.homOfLE (show X.basicOpen (f₁*(g*h)) ≤ X.basicOpen (f₁*g) by
        rw [show f₁*(g*h) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left) := by
  rw [← cancel_mono (X.basicOpen (f₁ * g)).ι]
  unfold transitionElemIso
  simp only [Iso.trans_hom, Iso.trans_inv, Iso.symm_hom, Iso.symm_inv, Category.assoc,
    Scheme.Hom.isoImage_hom_ι, Scheme.Hom.isoImage_hom_ι_assoc,
    Scheme.isoOfEq_hom, Scheme.homOfLE_ι, Scheme.homOfLE_ι_assoc]
  calc _ = ((eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        Z.homOfLE (Scheme.Hom.image_mono e.hom (transitionElem_restrict_mul_le f₁ g h))) ≫
        (Scheme.Hom.isoImage e.hom _).inv ≫ (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      congrArg (· ≫ (Scheme.Hom.isoImage e.hom _).inv ≫ (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι) (transitionElemIso_step1 f₁ g h e)
    _ = (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        (Z.homOfLE (Scheme.Hom.image_mono e.hom (transitionElem_restrict_mul_le f₁ g h)) ≫
          (Scheme.Hom.isoImage e.hom _).inv) ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      (Category.assoc _ _ _).symm
    _ = (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        ((Scheme.Hom.isoImage e.hom _).inv ≫ (X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le f₁ g h)) ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      congrArg (fun x =>
        (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫ x ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι) (transitionElemIso_step2 f₁ g h e).symm
    _ = (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        (Scheme.Hom.isoImage e.hom _).inv ≫
        ((X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le f₁ g h) ≫
          (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom) ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      (Category.assoc _ _ _).symm
    _ = (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        (Scheme.Hom.isoImage e.hom _).inv ≫
        ((Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
          X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le f₁ g h))) ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      congrArg (fun x =>
        (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        (Scheme.Hom.isoImage e.hom _).inv ≫ x ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
        (X.basicOpen (f₁*g)).ι) (transitionElemIso_step3 f₁ g h)
    _ = (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        (Scheme.Hom.isoImage e.hom _).inv ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le f₁ g h)) ≫
          (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e).1)).hom ≫
          (X.basicOpen (f₁*g)).ι) :=
      (Category.assoc _ _ _).symm
    _ = (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        (Scheme.Hom.isoImage e.hom _).inv ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        ((eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).1)).hom ≫
          (X.basicOpen (f₁*(g*h))).ι) :=
      congrArg (fun x =>
        (eqToIso (congrArg (fun (V : Z.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (g*h) e).2)).inv ≫
        (Scheme.Hom.isoImage e.hom _).inv ≫ (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫ x)
        (transitionElemIso_step45 f₁ g h e)

/-! ## `Scheme.GlueData`の配線——作業単位3の最終段

これまで積み上げた部品(`piece_descends_iso`・`exists_transitionIso`・
`transitionElem`)を使って、`CategoryTheory.GlueData`の12フィールドの
うち、まず**対象レベルのデータ**(`V`・`f`・`t`・`t'`)を組み立てる。
数学的な核心はすべて完成しているので、残りは純粋な配線(`f_mono`・
`f_hasPullback`・`f_id`は`basicOpen`の`.ι`が標準的に持つ性質から自動、
`t_id`は`gdT`が対角成分で打ち消し合うことから、`t'`は`transitionElem`の
乗法性から構成できる)。残る`t_fac`・`cocycle`は次の一手として残す。 -/

variable {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U)) (Z : J → Scheme)
  (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)

/-- 重なりの候補片(`V(i,j)`)——`i`側の候補片`Z i`の中の、`j`との重なりに
対応する基本開集合。`transitionElem`をそのまま使う。

`@[reducible]`——`gdF`の`instances`透明度での展開に必要(`t_fac`・`cocycle`で
`pullback.fst/snd (gdF …) (gdF …)`と`isPullback_opens_inf`由来の
`pullback.fst/snd U.ι V.ι`を`simp`で結びつけるとき、`gdF`が不透明だと
「`instances`透明度で型が合わない」という一見謎のエラーになる——`gdF`が
`.ι`へ展開できないせいで、同じ`pullback`対象なのに別物として扱われて
しまうため)。 -/
@[reducible]
noncomputable def gdV (p : J × J) : Scheme :=
  (Z p.1).basicOpen (transitionElem (f p.1) (f p.2) (e p.1))

/-- `V(i,j) ⟶ U i`——基本開集合の標準的な開埋め込み`.ι`。`@[reducible]`の理由は
`gdV`と同じ。 -/
@[reducible]
noncomputable def gdF (i j : J) : gdV f Z e (i, j) ⟶ Z i :=
  ((Z i).basicOpen (transitionElem (f i) (f j) (e i))).ι

/-- **遷移射`t(i,j) : V(i,j) ≅ V(j,i)`**——`transitionElem_basicOpen_eq`
(精密な等式)を`i`側・`j`側それぞれで使い、`isoImage`(開埋め込みの制限が
同型)+`eqToIso`で組み立てた2つの同型を、`X.basicOpen(f_i·f_j) =
X.basicOpen(f_j·f_i)`(可換環の乗法の可換性)で繋ぐ。GlueDataの遷移射
そのもの。

★**sorry 無し**。標準3公理のみ。`let`(`obtain`ではなく)で`.choose_spec`
相当の等式を束縛しているのが鍵——`obtain`(内部で`And.rec`を生成)だと
後で`unfold`しても`rec`が簡約されずスタックするが、`let`+`.1`/`.2`射影
なら`unfold`後に`simp`で素直に簡約できる(`tools/lean-idioms.md` #29)。 -/
noncomputable def gdT (i j : J) : gdV f Z e (i, j) ≅ gdV f Z e (j, i) :=
  let h1i := (transitionElem_basicOpen_eq (f i) (f j) (e i)).1
  let h2i := (transitionElem_basicOpen_eq (f i) (f j) (e i)).2
  let h1j := (transitionElem_basicOpen_eq (f j) (f i) (e j)).1
  let h2j := (transitionElem_basicOpen_eq (f j) (f i) (e j)).2
  let Wi := (X.basicOpen (f i)).toScheme.basicOpen
    ((X.basicOpen (f i)).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le (f i))).op (f j)))
  let Wj := (X.basicOpen (f j)).toScheme.basicOpen
    ((X.basicOpen (f j)).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le (f j))).op (f i)))
  let iso1 : (Wi : Scheme) ≅ gdV f Z e (i, j) :=
    ((e i).hom.isoImage Wi).trans (eqToIso (by rw [h2i]))
  let iso2 : (Wi : Scheme) ≅ (X.basicOpen (f i * f j) : Scheme) :=
    ((X.basicOpen (f i)).ι.isoImage Wi).trans (eqToIso (by rw [h1i]))
  let iso3 : (Wj : Scheme) ≅ gdV f Z e (j, i) :=
    ((e j).hom.isoImage Wj).trans (eqToIso (by rw [h2j]))
  let iso4 : (Wj : Scheme) ≅ (X.basicOpen (f j * f i) : Scheme) :=
    ((X.basicOpen (f j)).ι.isoImage Wj).trans (eqToIso (by rw [h1j]))
  let hcomm : X.basicOpen (f i * f j) = X.basicOpen (f j * f i) := by rw [mul_comm]
  iso1.symm.trans (iso2.trans ((X.isoOfEq hcomm).trans (iso4.symm.trans iso3)))

/-- **`t_id`の完成**——`i=j`のとき`gdT`は自明な同型(恒等射)になる。
`gdT`を定義する4つの部分同型が(`i=j`のとき)ペアで打ち消し合う。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdT_id (i : J) : gdT f Z e i i = Iso.refl _ := by
  apply Iso.ext
  unfold gdT
  simp only [show (X.isoOfEq (by rw [mul_comm] :
    X.basicOpen (f i * f i) = X.basicOpen (f i * f i))) = Iso.refl _ from by simp,
    Iso.refl_hom, Iso.trans_hom, Iso.symm_hom, Iso.trans_inv, Category.assoc, Category.id_comp]
  simp only [← Category.assoc]
  simp

/-- `gdF`は常にmono(基本開集合の標準的な開埋め込み`.ι`の一般的性質)——
GlueDataの`f_mono`フィールド。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdF_mono (i j : J) : Mono (gdF f Z e i j) :=
  inferInstanceAs (Mono (((Z i).basicOpen (transitionElem (f i) (f j) (e i))).ι))

/-- `gdF`は常に開埋め込み——GlueDataの`f_open`フィールド。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdF_isOpenImmersion (i j : J) : IsOpenImmersion (gdF f Z e i j) :=
  inferInstanceAs (IsOpenImmersion (((Z i).basicOpen (transitionElem (f i) (f j) (e i))).ι))

/-- `gdF i j`と`gdF i k`はpullbackを持つ(開埋め込みの一般的性質)——
GlueDataの`f_hasPullback`フィールド。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdF_hasPullback (i j k : J) : HasPullback (gdF f Z e i j) (gdF f Z e i k) := by
  unfold gdF; infer_instance

/-- スキームの同型`e:X≅Z`の下で`⊤`の像は`⊤`——CorrHyp非依存の一般的事実。 -/
theorem Scheme.hom_image_top_eq_top {X Z : Scheme} (e : X ≅ Z) :
    e.hom ''ᵁ (⊤ : X.Opens) = ⊤ := by
  apply TopologicalSpace.Opens.ext
  simp only [Scheme.Hom.coe_image, TopologicalSpace.Opens.coe_top]
  exact Set.image_univ_of_surjective e.hom.homeomorph.surjective

/-- 開集合`V`が`⊤`なら`V.ι`は同型——`X.topIso`に帰着するだけ
(CorrHyp非依存の一般的事実)。 -/
theorem isIso_ι_of_eq_top {X : Scheme} (V : X.Opens) (h : V = ⊤) : IsIso V.ι := by
  subst h
  rw [show (⊤ : X.Opens).ι = X.topIso.hom from Scheme.Hom.ext' rfl]
  infer_instance

/-- **`gdV`の対角成分(`i=j`)は`⊤`**——選ばれた座標が実際に`X.basicOpen f_i`
自身の定義元(`RingedSpace.isUnit_res_basicOpen`で単元と分かる)から来る
ことと、`e i`が同型であることから、重なり`X.basicOpen(f_i·f_i)`は
`Z i`全体に対応する。`f_id`(対角成分が同型)の核心部品。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdV_diag_eq_top (i : J) :
    (Z i).basicOpen (transitionElem (f i) (f i) (e i)) = ⊤ := by
  have h2 := (transitionElem_basicOpen_eq (f i) (f i) (e i)).2
  have hunit : IsUnit (X.presheaf.map (homOfLE (X.basicOpen_le (f i))).op (f i)) :=
    AlgebraicGeometry.RingedSpace.isUnit_res_basicOpen X.toRingedSpace (f i)
  have hWtop : (X.basicOpen (f i)).toScheme.basicOpen
      ((X.basicOpen (f i)).topIso.inv (X.presheaf.map (homOfLE (X.basicOpen_le (f i))).op (f i)))
      = ⊤ :=
    Scheme.basicOpen_of_isUnit (X.basicOpen (f i)) (hunit.map (X.basicOpen (f i)).topIso.inv.hom)
  rw [← h2, hWtop, Scheme.hom_image_top_eq_top]

/-- **`f_id`の完成**——`gdF`の対角成分(`i=i`)は同型。`gdV_diag_eq_top`
(候補片の対角成分が`⊤`であること)+`isIso_ι_of_eq_top`から直ちに従う。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdF_id (i : J) : IsIso (gdF f Z e i i) := by
  show IsIso ((Z i).basicOpen (transitionElem (f i) (f i) (e i))).ι
  exact isIso_ι_of_eq_top _ (gdV_diag_eq_top f Z e i)

/-- **候補片内の重なり(`pullback(gdF i j, gdF i k)`)は基本開集合として
実現できる**——`pullback U.ι V.ι ≅ U⊓V`(`isPullback_opens_inf`)+
`transitionElem_mul`(逆向き、`basicOpen_mul`と組み合わせて`⊓`を単一の
`basicOpen`にまとめる)。`t'`の構成の核心部品。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def gdVpullbackIso (i j k : J) :
    pullback (gdF f Z e i j) (gdF f Z e i k) ≅
      (Z i).basicOpen (transitionElem (f i) (f j * f k) (e i)) := by
  have hpb : pullback (gdF f Z e i j) (gdF f Z e i k) ≅
      (((Z i).basicOpen (transitionElem (f i) (f j) (e i))) ⊓
        ((Z i).basicOpen (transitionElem (f i) (f k) (e i))) : (Z i).Opens) := by
    unfold gdF
    exact (isPullback_opens_inf _ _).isoPullback.symm
  have heq : (((Z i).basicOpen (transitionElem (f i) (f j) (e i))) ⊓
        ((Z i).basicOpen (transitionElem (f i) (f k) (e i))) : (Z i).Opens) =
      (Z i).basicOpen (transitionElem (f i) (f j * f k) (e i)) := by
    rw [transitionElem_mul, Scheme.basicOpen_mul]
  exact hpb.trans ((Z i).isoOfEq heq)

/-- **`t'(i,j,k) : pullback(f i j, f i k) ≅ pullback(f j k, f j i)`
——GlueDataの`t'`フィールドそのもの**。`gdVpullbackIso`(候補片内の重なり
を基本開集合として実現)+`transitionElemIso`(基本開集合を`X`の対応する
基本開集合に戻す)+乗法の可換性・結合性(`ring`)で`i`側・`j`側を橋渡し
する。3重の重なりに新しい数学的内容は不要——`transitionElem`が任意の
積を受け付けることから`exists_transitionIso_finset`で確認したとおり。

これで`Scheme.GlueData`の12フィールド中11個
(`J・U・V・f・t・t_id・f_mono・f_open・f_hasPullback・f_id・t'`)が
完成した。残るのは`t_fac`・`cocycle`(`t'`が満たす2つの可換性の等式)
のみ。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def gdT' (i j k : J) :
    pullback (gdF f Z e i j) (gdF f Z e i k) ≅ pullback (gdF f Z e j k) (gdF f Z e j i) :=
  (gdVpullbackIso f Z e i j k).trans <|
    (transitionElemIso (f i) (f j * f k) (e i)).symm.trans <|
    (X.isoOfEq (by rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).trans <|
    (transitionElemIso (f j) (f k * f i) (e j)).trans (gdVpullbackIso f Z e j k i).symm

/-- **`gdT`を`transitionElemIso`だけで書き直した形**——`gdT`本来の定義
(`isoImage`4段の塔)と、`transitionElemIso`(こちらも内部は同じ塔だが、
自然性補題`transitionElemIso_inv/hom_naturality`をすでに持っている)が
実は同じ射であることを`unfold`+`rfl`で直接確認する。`t_fac`で`gdT`の
内部にもう一度立ち入らずに済むための橋渡し。

★**sorry 無し**。標準3公理のみ。`unfold gdT transitionElemIso; rfl`が
そのまま通る——`gdT`の`iso1`〜`iso4`は`transitionElemIso`内部の
`iso1_TEI`/`iso2_TEI`と定義上完全に同じ組み立て方をしているため。 -/
theorem gdT_eq_transitionElemIso (i j : J) :
    (gdT f Z e i j).hom = (transitionElemIso (f i) (f j) (e i)).inv ≫
      (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom ≫
      (transitionElemIso (f j) (f i) (e j)).hom := by
  unfold gdT transitionElemIso
  rfl

/-- `X.isoOfEq`の`.hom`と`X.homOfLE`の合成が単一の`X.homOfLE`にまとまる
——`Scheme.isoOfEq_hom`+`Scheme.homOfLE_homOfLE`を繋ぐだけ。`@[reassoc]`
を付け、末尾に何か続く形の書き換えにも使えるようにする
(`gdT_hom_comp_gdF`で使う)。CorrHyp非依存の一般的事実。 -/
@[reassoc]
theorem homOfLE_isoOfEq_comp' {X : Scheme} {A B C : X.Opens} (heq : A = B) (hBC : B ≤ C) (hAC : A ≤ C) :
    (X.isoOfEq heq).hom ≫ X.homOfLE hBC = X.homOfLE hAC := by
  rw [Scheme.isoOfEq_hom, Scheme.homOfLE_homOfLE]

/-- **`gdT`(項目(b)のNatIso自然性の`snd`側で必要)と`gdF`(反対側)の
合成——`transitionElemIso`だけの式へ帰着**——`gdT_eq_transitionElemIso`
+`transitionElemIso_hom_ι`(`j`側)+`homOfLE_isoOfEq_comp'`で
「`(gdT i j).hom ≫ gdF(j,i)`」を「`transitionElemIso(f i)(f j)(e i)`の
`.inv`(`i`側はそのまま)+`X.homOfLE`+`(e j).hom`」の形へ書き直す。

★配管の注意: `rw[gdT_eq_transitionElemIso]`直後に`show`で式全体を
再掲示すると、埋め込まれた`(X.isoOfEq ...)`の証明項同士の(命題として
等しいが構文的には別の)照合で`instances`透明度の壁に当たり
`whnf`が終わらない(`#31`の新しい現れ方)——`show`を避け、代わりに
`rw[Category.assoc, Category.assoc]`で結合を先に揃えてから
`transitionElemIso_hom_ι`・`homOfLE_isoOfEq_comp'_assoc`を適用する
ことで回避した(元の証明項に一切触れないので壁に当たらない)。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdT_hom_comp_gdF (i j : J) :
    (gdT f Z e i j).hom ≫ gdF f Z e j i =
      (transitionElemIso (f i) (f j) (e i)).inv ≫
      X.homOfLE (show X.basicOpen (f i * f j) ≤ X.basicOpen (f j) by
        rw [X.basicOpen_mul (f i) (f j)]; exact inf_le_right) ≫ (e j).hom := by
  show (gdT f Z e i j).hom ≫ ((Z j).basicOpen (transitionElem (f j) (f i) (e j))).ι = _
  rw [gdT_eq_transitionElemIso, Category.assoc, Category.assoc, transitionElemIso_hom_ι,
    homOfLE_isoOfEq_comp'_assoc]

/-- `gdVpullbackIso`の`.hom`と`pullback.fst`の関係——`isPullback_opens_
inf`が与える`isoPullback_inv_fst`を、`gdVpullbackIso`自身の構成
(`hpb.trans(isoOfEq heq)`)と組み合わせて`(Z i).homOfLE`1つにまとめる。
`t_fac`で`pullback.fst`を`gdVpullbackIso`側の言葉へ変換するための部品。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdVpullbackIso_hom_comp_homOfLE_fst (i j k : J) :
    (gdVpullbackIso f Z e i j k).hom ≫
      (Z i).homOfLE (show (Z i).basicOpen (transitionElem (f i) (f j * f k) (e i)) ≤
        (Z i).basicOpen (transitionElem (f i) (f j) (e i)) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_left)
      = pullback.fst (gdF f Z e i j) (gdF f Z e i k) := by
  unfold gdVpullbackIso gdF
  simp only [Iso.trans_hom, Category.assoc, Scheme.isoOfEq_hom, Scheme.homOfLE_homOfLE, id_eq,
    Iso.symm_hom, IsPullback.isoPullback_inv_fst]

/-- `gdVpullbackIso_hom_comp_homOfLE_fst`の`snd`版。 -/
theorem gdVpullbackIso_hom_comp_homOfLE_snd (i j k : J) :
    (gdVpullbackIso f Z e i j k).hom ≫
      (Z i).homOfLE (show (Z i).basicOpen (transitionElem (f i) (f j * f k) (e i)) ≤
        (Z i).basicOpen (transitionElem (f i) (f k) (e i)) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
      = pullback.snd (gdF f Z e i j) (gdF f Z e i k) := by
  unfold gdVpullbackIso gdF
  simp only [Iso.trans_hom, Category.assoc, Scheme.isoOfEq_hom, Scheme.homOfLE_homOfLE, id_eq,
    Iso.symm_hom, IsPullback.isoPullback_inv_snd]

/-- `gdVpullbackIso_hom_comp_homOfLE_snd`の`.inv`版(`.hom`を消して
`pullback.snd`と`(Z i).homOfLE`だけの関係にする)——`gdT'`の末尾
`(gdVpullbackIso …).inv ≫ pullback.snd(…)`をそのまま`(Z i).homOfLE`へ
変換するための部品。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdVpullbackIso_inv_comp_pullback_snd (i j k : J) :
    (gdVpullbackIso f Z e i j k).inv ≫ pullback.snd (gdF f Z e i j) (gdF f Z e i k)
      = (Z i).homOfLE (show (Z i).basicOpen (transitionElem (f i) (f j * f k) (e i)) ≤
        (Z i).basicOpen (transitionElem (f i) (f k) (e i)) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right) := by
  have h := gdVpullbackIso_hom_comp_homOfLE_snd f Z e i j k
  have step := congrArg ((gdVpullbackIso f Z e i j k).inv ≫ ·) h
  simp only [← Category.assoc, Iso.inv_hom_id, Category.id_comp] at step
  exact step.symm

/-! ### `transitionElemIso`の自然性、積の順序を入れ替えた版(`h·g`)

`t_fac`の`j`側(`gdT'`が`transitionElemIso (f j) (f k * f i) (e j)`を
使うのに対し、`gdT`は`transitionElemIso (f j) (f i) (e j)`を使う)では、
生き残る元(`f i`)が積の**後ろ**に来る(`f k * f i`)。前に証明した
`transitionElemIso_inv/hom_naturality`は生き残る元が**前**に来る形
(`g * h`、`g`が生き残る)だったので、同じ補題は直接使えない——`h * g`
(`h`が消え、`g`が生き残る)の場合を独立に証明し直す必要がある
(`inf_le_left`/`inf_le_right`が入れ替わるだけで、証明の骨格は完全に
同じ)。 -/

theorem transitionElem_restrict_mul_le2 (f₁ h g : Γ(X, U)) :
    (X.basicOpen f₁ : Scheme).basicOpen
      ((ConcreteCategory.hom (X.basicOpen f₁).topIso.inv)
        ((ConcreteCategory.hom (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op)) (h * g)))
    ≤ (X.basicOpen f₁ : Scheme).basicOpen
      ((ConcreteCategory.hom (X.basicOpen f₁).topIso.inv)
        ((ConcreteCategory.hom (X.presheaf.map (homOfLE (X.basicOpen_le f₁)).op)) g)) := by
  rw [map_mul, map_mul, Scheme.basicOpen_mul]; exact inf_le_right

theorem transitionElemIso_step1' (f₁ h g : Γ(X, U)) {Z' : Scheme} (e' : (X.basicOpen f₁ : Scheme) ≅ Z') :
    Z'.homOfLE (show Z'.basicOpen (transitionElem f₁ (h*g) e') ≤ Z'.basicOpen (transitionElem f₁ g e') by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
      ≫ (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').2)).inv
    = (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
      Z'.homOfLE (Scheme.Hom.image_mono e'.hom (transitionElem_restrict_mul_le2 f₁ h g)) :=
  eqToIso_homOfLE_comm (transitionElem_basicOpen_eq f₁ g e').2 (transitionElem_basicOpen_eq f₁ (h*g) e').2
    (Scheme.Hom.image_mono e'.hom (transitionElem_restrict_mul_le2 f₁ h g))
    (show Z'.basicOpen (transitionElem f₁ (h*g) e') ≤ Z'.basicOpen (transitionElem f₁ g e') by
      rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
    (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').2)
    (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)

theorem transitionElemIso_step2' (f₁ h g : Γ(X, U)) {Z' : Scheme} (e' : (X.basicOpen f₁ : Scheme) ≅ Z') :
    (Scheme.Hom.isoImage e'.hom _).inv ≫ (X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le2 f₁ h g)
    = Z'.homOfLE (Scheme.Hom.image_mono e'.hom (transitionElem_restrict_mul_le2 f₁ h g)) ≫
        (Scheme.Hom.isoImage e'.hom _).inv :=
  Scheme.Hom.isoImage_inv_naturality e'.hom (transitionElem_restrict_mul_le2 f₁ h g)

theorem transitionElemIso_step3' (f₁ h g : Γ(X, U)) :
    (X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le2 f₁ h g) ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom
    = (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le2 f₁ h g)) :=
  Scheme.Hom.isoImage_naturality (X.basicOpen f₁).ι (transitionElem_restrict_mul_le2 f₁ h g)

theorem transitionElemIso_step4' (f₁ h g : Γ(X, U)) {Z' : Scheme} (e' : (X.basicOpen f₁ : Scheme) ≅ Z') :
    X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le2 f₁ h g)) ≫
      (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom
    = (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').1)).hom ≫
      X.homOfLE (show X.basicOpen (f₁*(h*g)) ≤ X.basicOpen (f₁*g) by
        rw [show f₁*(h*g) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left) :=
  eqToIso_homOfLE_comm' (transitionElem_basicOpen_eq f₁ g e').1 (transitionElem_basicOpen_eq f₁ (h*g) e').1
    (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le2 f₁ h g))
    (show X.basicOpen (f₁*(h*g)) ≤ X.basicOpen (f₁*g) by
      rw [show f₁*(h*g) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left)
    (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)
    (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').1)

theorem transitionElemIso_step45' (f₁ h g : Γ(X, U)) {Z' : Scheme} (e' : (X.basicOpen f₁ : Scheme) ≅ Z') :
    X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le2 f₁ h g)) ≫
      (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
      (X.basicOpen (f₁*g)).ι
    = (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').1)).hom ≫
      (X.basicOpen (f₁*(h*g))).ι := by
  rw [← Category.assoc, transitionElemIso_step4' f₁ h g e', Category.assoc,
    Scheme.homOfLE_ι X (show X.basicOpen (f₁*(h*g)) ≤ X.basicOpen (f₁*g) by
      rw [show f₁*(h*g) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left)]

set_option maxHeartbeats 1000000 in
/-- **`transitionElemIso`の自然性(積の順序`h·g`版)**——`transitionElemIso_
inv_naturality`と同じ主張を、積の順序を入れ替えて示す。証明は完全に同じ
term-mode の`calc`+`congrArg`+`(Category.assoc _ _ _).symm`パターン
(`tools/lean-idioms.md` #31)。

★**sorry 無し**。標準3公理のみ。 -/
theorem transitionElemIso_inv_naturality2 (f₁ h g : Γ(X, U)) {Z' : Scheme} (e' : (X.basicOpen f₁ : Scheme) ≅ Z') :
    Z'.homOfLE (show Z'.basicOpen (transitionElem f₁ (h*g) e') ≤ Z'.basicOpen (transitionElem f₁ g e') by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
      ≫ (transitionElemIso f₁ g e').inv
    = (transitionElemIso f₁ (h*g) e').inv ≫
      X.homOfLE (show X.basicOpen (f₁*(h*g)) ≤ X.basicOpen (f₁*g) by
        rw [show f₁*(h*g) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left) := by
  rw [← cancel_mono (X.basicOpen (f₁ * g)).ι]
  unfold transitionElemIso
  simp only [Iso.trans_hom, Iso.trans_inv, Iso.symm_hom, Iso.symm_inv, Category.assoc,
    Scheme.Hom.isoImage_hom_ι, Scheme.Hom.isoImage_hom_ι_assoc,
    Scheme.isoOfEq_hom, Scheme.homOfLE_ι, Scheme.homOfLE_ι_assoc]
  calc _ = ((eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        Z'.homOfLE (Scheme.Hom.image_mono e'.hom (transitionElem_restrict_mul_le2 f₁ h g))) ≫
        (Scheme.Hom.isoImage e'.hom _).inv ≫ (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      congrArg (· ≫ (Scheme.Hom.isoImage e'.hom _).inv ≫ (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι) (transitionElemIso_step1' f₁ h g e')
    _ = (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        (Z'.homOfLE (Scheme.Hom.image_mono e'.hom (transitionElem_restrict_mul_le2 f₁ h g)) ≫
          (Scheme.Hom.isoImage e'.hom _).inv) ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      (Category.assoc _ _ _).symm
    _ = (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        ((Scheme.Hom.isoImage e'.hom _).inv ≫ (X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le2 f₁ h g)) ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      congrArg (fun x =>
        (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫ x ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι) (transitionElemIso_step2' f₁ h g e').symm
    _ = (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        (Scheme.Hom.isoImage e'.hom _).inv ≫
        ((X.basicOpen f₁ : Scheme).homOfLE (transitionElem_restrict_mul_le2 f₁ h g) ≫
          (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom) ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      (Category.assoc _ _ _).symm
    _ = (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        (Scheme.Hom.isoImage e'.hom _).inv ≫
        ((Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
          X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le2 f₁ h g))) ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι :=
      congrArg (fun x =>
        (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        (Scheme.Hom.isoImage e'.hom _).inv ≫ x ≫
        (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
        (X.basicOpen (f₁*g)).ι) (transitionElemIso_step3' f₁ h g)
    _ = (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        (Scheme.Hom.isoImage e'.hom _).inv ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        (X.homOfLE (Scheme.Hom.image_mono (X.basicOpen f₁).ι (transitionElem_restrict_mul_le2 f₁ h g)) ≫
          (eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ g e').1)).hom ≫
          (X.basicOpen (f₁*g)).ι) :=
      (Category.assoc _ _ _).symm
    _ = (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        (Scheme.Hom.isoImage e'.hom _).inv ≫
        (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫
        ((eqToIso (congrArg (fun (V : X.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').1)).hom ≫
          (X.basicOpen (f₁*(h*g))).ι) :=
      congrArg (fun x =>
        (eqToIso (congrArg (fun (V : Z'.Opens) => (V : Scheme)) (transitionElem_basicOpen_eq f₁ (h*g) e').2)).inv ≫
        (Scheme.Hom.isoImage e'.hom _).inv ≫ (Scheme.Hom.isoImage (X.basicOpen f₁).ι _).hom ≫ x)
        (transitionElemIso_step45' f₁ h g e')

/-- `transitionElemIso_inv_naturality2`の`.hom`版(向きが逆)——`.inv`版
から`cancel_epi`/`cancel_mono`相当の代数操作(`Iso.hom_inv_id`等)で導く。
`gdT'_key`で`gdT'`のj側(`transitionElemIso (f j) (f k*f i) (e j)`)を
`gdT`のj側(`transitionElemIso (f j) (f i) (e j)`)へ橋渡しするのに使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem transitionElemIso_hom_naturality2 (f₁ h g : Γ(X, U)) {Z' : Scheme} (e' : (X.basicOpen f₁ : Scheme) ≅ Z') :
    (transitionElemIso f₁ (h*g) e').hom ≫
      Z'.homOfLE (show Z'.basicOpen (transitionElem f₁ (h*g) e') ≤ Z'.basicOpen (transitionElem f₁ g e') by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
    = X.homOfLE (show X.basicOpen (f₁*(h*g)) ≤ X.basicOpen (f₁*g) by
        rw [show f₁*(h*g) = (f₁*g)*h by ring, Scheme.basicOpen_mul]; exact inf_le_left)
      ≫ (transitionElemIso f₁ g e').hom := by
  have hN := transitionElemIso_inv_naturality2 f₁ h g e'
  have step1 := congrArg ((transitionElemIso f₁ (h*g) e').hom ≫ ·) hN
  simp only [← Category.assoc, Iso.hom_inv_id, Category.id_comp] at step1
  have step2 := congrArg (· ≫ (transitionElemIso f₁ g e').hom) step1
  simp only [Category.assoc, Iso.inv_hom_id, Category.comp_id] at step2
  exact step2

/-- **`t_fac`の核心の等式**——`gdT'`のi側・j側それぞれで`transitionElemIso`
の自然性(`_inv_naturality`・`_hom_naturality2`)を適用し、残った純粋な
`X`レベルの結合律・可換律の等式(`isoOfEq`+`homOfLE`の合成、証明無関係性
で自動的に閉じる)で橋渡しする。`gdVpullbackIso`の`.hom`をキャンセルした
あとに残る本体。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdT'_key (i j k : J) :
    (transitionElemIso (f i) (f j * f k) (e i)).inv ≫
        (X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
          rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).hom ≫
        (transitionElemIso (f j) (f k * f i) (e j)).hom ≫
        (Z j).homOfLE (show (Z j).basicOpen (transitionElem (f j) (f k * f i) (e j)) ≤
          (Z j).basicOpen (transitionElem (f j) (f i) (e j)) by
          rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
      = (Z i).homOfLE (show (Z i).basicOpen (transitionElem (f i) (f j * f k) (e i)) ≤
          (Z i).basicOpen (transitionElem (f i) (f j) (e i)) by
          rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_left) ≫
        (transitionElemIso (f i) (f j) (e i)).inv ≫
        (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom ≫
        (transitionElemIso (f j) (f i) (e j)).hom := by
  have key3 : (X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
        rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).hom ≫
      X.homOfLE (show X.basicOpen (f j * (f k * f i)) ≤ X.basicOpen (f j * f i) by
        rw [show f j * (f k * f i) = (f j * f i) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left)
    = X.homOfLE (show X.basicOpen (f i * (f j * f k)) ≤ X.basicOpen (f i * f j) by
        rw [show f i * (f j * f k) = (f i * f j) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left) ≫
      (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom := by
    simp only [Scheme.isoOfEq_hom, Scheme.homOfLE_homOfLE]
  have key2 : (X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
        rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).hom ≫
      (transitionElemIso (f j) (f k * f i) (e j)).hom ≫
      (Z j).homOfLE (show (Z j).basicOpen (transitionElem (f j) (f k * f i) (e j)) ≤
        (Z j).basicOpen (transitionElem (f j) (f i) (e j)) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
    = X.homOfLE (show X.basicOpen (f i * (f j * f k)) ≤ X.basicOpen (f i * f j) by
        rw [show f i * (f j * f k) = (f i * f j) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left) ≫
      (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom ≫
      (transitionElemIso (f j) (f i) (e j)).hom := by
    have hbridge := transitionElemIso_hom_naturality2 (f j) (f k) (f i) (e j)
    calc _ = (X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
          rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).hom ≫
        (X.homOfLE (show X.basicOpen (f j * (f k * f i)) ≤ X.basicOpen (f j * f i) by
          rw [show f j * (f k * f i) = (f j * f i) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left) ≫
        (transitionElemIso (f j) (f i) (e j)).hom) :=
      congrArg (fun x => (X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
          rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).hom ≫ x) hbridge
    _ = ((X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
          rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).hom ≫
        X.homOfLE (show X.basicOpen (f j * (f k * f i)) ≤ X.basicOpen (f j * f i) by
          rw [show f j * (f k * f i) = (f j * f i) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left)) ≫
        (transitionElemIso (f j) (f i) (e j)).hom :=
      (Category.assoc _ _ _).symm
    _ = (X.homOfLE (show X.basicOpen (f i * (f j * f k)) ≤ X.basicOpen (f i * f j) by
        rw [show f i * (f j * f k) = (f i * f j) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left) ≫
        (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom) ≫
        (transitionElemIso (f j) (f i) (e j)).hom :=
      congrArg (· ≫ (transitionElemIso (f j) (f i) (e j)).hom) key3
    _ = _ := Category.assoc _ _ _
  have hbridge2 := transitionElemIso_inv_naturality (f i) (f j) (f k) (e i)
  calc _ = (transitionElemIso (f i) (f j * f k) (e i)).inv ≫
      (X.homOfLE (show X.basicOpen (f i * (f j * f k)) ≤ X.basicOpen (f i * f j) by
        rw [show f i * (f j * f k) = (f i * f j) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left) ≫
      (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom ≫
      (transitionElemIso (f j) (f i) (e j)).hom) :=
    congrArg (fun x => (transitionElemIso (f i) (f j * f k) (e i)).inv ≫ x) key2
  _ = ((transitionElemIso (f i) (f j * f k) (e i)).inv ≫
      X.homOfLE (show X.basicOpen (f i * (f j * f k)) ≤ X.basicOpen (f i * f j) by
        rw [show f i * (f j * f k) = (f i * f j) * f k by ring, Scheme.basicOpen_mul]; exact inf_le_left)) ≫
      (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom ≫
      (transitionElemIso (f j) (f i) (e j)).hom :=
    (Category.assoc _ _ _).symm
  _ = ((Z i).homOfLE (show (Z i).basicOpen (transitionElem (f i) (f j * f k) (e i)) ≤
        (Z i).basicOpen (transitionElem (f i) (f j) (e i)) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_left) ≫
      (transitionElemIso (f i) (f j) (e i)).inv) ≫
      (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom ≫
      (transitionElemIso (f j) (f i) (e j)).hom :=
    congrArg (· ≫ (X.isoOfEq (show X.basicOpen (f i * f j) = X.basicOpen (f j * f i) by rw [mul_comm])).hom ≫
      (transitionElemIso (f j) (f i) (e j)).hom) hbridge2.symm
  _ = _ := Category.assoc _ _ _

/-- **`t_fac`——`Scheme.GlueData`の11個目のフィールド**。`gdT'`(i側・j側
それぞれで`gdVpullbackIso`+`transitionElemIso`を経由する)と`gdT`
(`transitionElemIso`2つを`isoOfEq(可換律)`で橋渡し)が、`pullback.fst`・
`pullback.snd`を経由して両立することを示す。`gdT'`を`gdT_eq_
transitionElemIso`で書き換え、`pullback.fst`/`pullback.snd`を
`gdVpullbackIso_hom_comp_homOfLE_fst`/`gdVpullbackIso_inv_comp_pullback_
snd`で`(Z i).homOfLE`の言葉に変換し、共通の`(gdVpullbackIso i j k).hom`
を`congrArg`でキャンセルすれば`gdT'_key`にちょうど一致する。

これで`Scheme.GlueData`の12フィールド中11個
(`J・U・V・f・t・t_id・f_mono・f_open・f_hasPullback・f_id・t'・t_fac`)
が完成した。残るのは`cocycle`(`t' i j k ≫ t' j k i ≫ t' k i j = 𝟙`)
のみ。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdT'_t_fac (i j k : J) :
    (gdT' f Z e i j k).hom ≫ pullback.snd (gdF f Z e j k) (gdF f Z e j i)
      = pullback.fst (gdF f Z e i j) (gdF f Z e i k) ≫ (gdT f Z e i j).hom := by
  unfold gdT'
  simp only [Iso.trans_hom, Iso.symm_hom, Category.assoc]
  rw [gdT_eq_transitionElemIso, ← gdVpullbackIso_hom_comp_homOfLE_fst f Z e i j k]
  rw [show (gdVpullbackIso f Z e j k i).inv ≫ pullback.snd (gdF f Z e j k) (gdF f Z e j i)
      = (Z j).homOfLE (show (Z j).basicOpen (transitionElem (f j) (f k * f i) (e j)) ≤
          (Z j).basicOpen (transitionElem (f j) (f i) (e j)) by
        rw [transitionElem_mul, Scheme.basicOpen_mul]; exact inf_le_right)
      from gdVpullbackIso_inv_comp_pullback_snd f Z e j k i]
  rw [← Category.assoc]
  exact congrArg ((gdVpullbackIso f Z e i j k).hom ≫ ·) (gdT'_key f Z e i j k)

/-! ### `cocycle`——`Scheme.GlueData`最後のフィールド

**配管の教訓(`tools/lean-idioms.md`へ追記予定)**: `cocycle`は当初
`t_fac`と同じ「`unfold gdT'`してから`calc`/`congrArg`で組み立てる」
方法を試みたが、**`unfold gdT'`をゴールに含まれる`gdT'`の`2つ以上の
出現に同時に適用する**と、それ以降のあらゆる型検査(`rw`・`simp`・
`exact`・`refine`・`show`すべて)が極端に重くなる(`maxHeartbeats`を
2000万まで上げても完走しない)ことが分かった——`t_fac`は`gdT'`の
出現が常に1つだけだったので問題にならなかった。**回避策は`gdT'`を
1つずつ`gdT'_hom_eq`/`gdT'_inv_eq`という名前付きの事実として先に
確定させ、以降は`unfold`せず`rw`でこれらの名前を参照するだけにする**
こと——これなら`gdT'`が何度出現しても軽い。

もう1つの教訓: `congrArg`で式の一部を書き換えるとき、**書き換え対象を
2つの合成の"間に挟む"**(`congrArg (fun x => A ≫ x ≫ B) h`)と極端に
重くなることがある——`A`と`B`の型がジェネリックな`x`に対して整合する
ことを確認する型検査が高くつくと見られる。**常に「前だけ」
(`congrArg (· ≫ K) h`)か「後ろだけ」(`congrArg (K ≫ ·) h`)の
`congrArg`を順番に適用する**(`t_fac`の`gdT'_key`と同じ流儀)ことで
劇的に軽くなる。 -/

/-- `(gdT' i j k).hom`を`unfold`せずに参照できる、明示的な形。`gdT'`は
`gdVpullbackIso`+`transitionElemIso`+`isoOfEq`の`.trans`/`.symm`だけで
組み立てられているので、`unfold`+`Iso.trans_hom`/`Iso.symm_hom`だけで
安全に取り出せる(1つの`gdT'`だけを`unfold`するので軽い)。 -/
theorem gdT'_hom_eq (i j k : J) :
    (gdT' f Z e i j k).hom = (gdVpullbackIso f Z e i j k).hom ≫
      (transitionElemIso (f i) (f j * f k) (e i)).inv ≫
      (X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
        rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).hom ≫
      (transitionElemIso (f j) (f k * f i) (e j)).hom ≫ (gdVpullbackIso f Z e j k i).inv := by
  unfold gdT'
  simp only [Iso.trans_hom, Iso.symm_hom, Category.assoc]

/-- `gdT'_hom_eq`の`.inv`版。 -/
theorem gdT'_inv_eq (i j k : J) :
    (gdT' f Z e i j k).inv = (gdVpullbackIso f Z e j k i).hom ≫
      (transitionElemIso (f j) (f k * f i) (e j)).inv ≫
      (X.isoOfEq (show X.basicOpen (f i * (f j * f k)) = X.basicOpen (f j * (f k * f i)) by
        rw [show f i * (f j * f k) = f j * (f k * f i) by ring])).inv ≫
      (transitionElemIso (f i) (f j * f k) (e i)).hom ≫ (gdVpullbackIso f Z e i j k).inv := by
  unfold gdT'
  simp only [Iso.trans_inv, Iso.symm_inv, Iso.symm_hom, Category.assoc]

/-- **`cocycle`の核心の等式**——残る2つの`gdT'`(`j k i`・`k i j`)の合成が
1つ目(`i j k`)の逆射に一致する。`gdT'_hom_eq`/`gdT'_inv_eq`で3つの
`gdT'`をそれぞれ独立に(`unfold`を経由せず)明示形へ変換し、共通の
`gdVpullbackIso`・`transitionElemIso`の対をキャンセルすると、残るのは
3つの結合律の`homOfLE`が1つに合成される等式(`Scheme.homOfLE_homOfLE`、
証明無関係性で自動的に閉じる)だけになる。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdT'_pair (i j k : J) :
    (gdT' f Z e j k i).hom ≫ (gdT' f Z e k i j).hom = (gdT' f Z e i j k).inv := by
  rw [gdT'_hom_eq f Z e j k i, gdT'_hom_eq f Z e k i j, gdT'_inv_eq f Z e i j k]
  simp only [Category.assoc, Iso.inv_hom_id_assoc, Iso.hom_inv_id_assoc,
    Scheme.isoOfEq_hom, Scheme.isoOfEq_inv]
  have hcombine : X.homOfLE (show X.basicOpen (f j * (f k * f i)) ≤ X.basicOpen (f k * (f i * f j)) by
        rw [show f j * (f k * f i) = f k * (f i * f j) by ring]) ≫
      X.homOfLE (show X.basicOpen (f k * (f i * f j)) ≤ X.basicOpen (f i * (f j * f k)) by
        rw [show f k * (f i * f j) = f i * (f j * f k) by ring])
    = X.homOfLE (show X.basicOpen (f j * (f k * f i)) ≤ X.basicOpen (f i * (f j * f k)) by
        rw [show f j * (f k * f i) = f i * (f j * f k) by ring]) :=
    Scheme.homOfLE_homOfLE X _ _
  have step1 := congrArg (· ≫ (transitionElemIso (f i) (f j * f k) (e i)).hom ≫
      (gdVpullbackIso f Z e i j k).inv) hcombine
  have step2 := congrArg ((gdVpullbackIso f Z e j k i).hom ≫
      (transitionElemIso (f j) (f k * f i) (e j)).inv ≫ ·) step1
  simpa [Category.assoc] using step2

/-- **`cocycle`——`Scheme.GlueData`の12個目、最後のフィールド**。
`gdT'_pair`(残る2つの`gdT'`の合成が1つ目の逆射に一致する)と
`Iso.hom_inv_id`から直ちに従う。

これで`Scheme.GlueData`の**12フィールドすべて**
(`J・U・V・f・f_mono・f_hasPullback・f_id・t・t_id・t'・t_fac・cocycle`)
が完成した。

★**sorry 無し**。標準3公理のみ。 -/
theorem gdT'_cocycle (i j k : J) :
    (gdT' f Z e i j k).hom ≫ (gdT' f Z e j k i).hom ≫ (gdT' f Z e k i j).hom = 𝟙 _ := by
  have h := congrArg ((gdT' f Z e i j k).hom ≫ ·) (gdT'_pair f Z e i j k)
  simpa [Category.assoc] using h

/-- **`corrHypGlueData`——`Scheme.GlueData`の構造体インスタンス本体**。
ここまで積み上げた12個の部品(`gdV`・`gdF`・`gdF_mono`・
`gdF_hasPullback`・`gdF_id`・`gdT`・`gdT_id`・`gdT'`・`gdT'_t_fac`・
`gdT'_cocycle`・`gdF_isOpenImmersion`)をそのままフィールドとして渡す
だけで組み立てられる——`t`・`t'`は`Iso`(`gdT`・`gdT'`)の`.hom`を取る
だけ、`t_id`は`gdT_id`(`gdT i i = Iso.refl _`)から`.hom`を取って
`rfl`で閉じるだけ。

これで**`Scheme.GlueData`(候補片`Z`の族から貼り合わせデータへ)が
完全に組み立てられた**——`Lemma 4.1`の構成的降下のうち、GlueDataの
配線という最大の技術的部分が完了した。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypGlueData : Scheme.GlueData where
  toGlueData := {
    J := J
    U := Z
    V := gdV f Z e
    f := gdF f Z e
    f_mono := gdF_mono f Z e
    f_hasPullback := gdF_hasPullback f Z e
    f_id := gdF_id f Z e
    t := fun i j => (gdT f Z e i j).hom
    t_id := fun i => by rw [gdT_id]; rfl
    t' := fun i j k => (gdT' f Z e i j k).hom
    t_fac := gdT'_t_fac f Z e
    cocycle := gdT'_cocycle f Z e
  }
  f_open := gdF_isOpenImmersion f Z e

/-! ### ロードマップ項目(b)への第一歩——`U`自身の`OpenCover`

`corrHypGlueData.glued`(貼り合わせてできるスキーム)が`U`に同型である
ことを示すには、比較対象として「`Z`の族が実際に`U`を覆っている」ことを
体現する`Scheme.OpenCover`が要る。`f i`が`U`を(基本開集合として)覆う
という前提(`⨆ i, X.basicOpen (f i) = U`)から、`mathlib`の
`Scheme.Cover.mkOfCovers`(全射性の条件だけからOpenCoverを組み立てる
スマートコンストラクタ)で直接構成する——各`i`ごとの写像は
`(e i).inv ≫ X.homOfLE (…) : Z i ⟶ U`(`X.basicOpen (f i) ⊆ U`への
包含を`e i`で`Z i`側へ転送したもの)。

配管の注意: 全射性の証明に出てくる被覆条件の不等式
(`X.basicOpen (f i) ≤ U`)を`have`でインライン展開すると
(`#31`と同系統の)`whnf`のheartbeat上限に達して詰まる——
`basicOpen_le_of_iSup_eq`という独立した`theorem`に先出しして名前で
参照するだけで解消する(`#31`の「部品は先に独立した`theorem`にする」
という教訓の再確認)。 -/

/-- `f i`による基本開集合は常に`U`に含まれる——`⨆ i, X.basicOpen (f i) = U`
という被覆条件から。`piecesOpenCover`の全射性の証明をインラインで書くと
`whnf`のheartbeat上限に達するので、先に独立した`theorem`として確定させた
もの(`#31`と同系統の教訓)。 -/
theorem basicOpen_le_of_iSup_eq {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (hcover : ⨆ i, X.basicOpen (f i) = U) (i : J) : X.basicOpen (f i) ≤ U :=
  (le_iSup (fun i => X.basicOpen (f i)) i).trans_eq hcover

set_option maxHeartbeats 1000000 in
/-- **`Z`の族が`U`を覆うという被覆条件から、`U`自身の`OpenCover`を直接
構成する**——`Scheme.Cover.mkOfCovers`(全射性の条件だけからOpenCoverを
組み立てるスマートコンストラクタ)に、各`i`ごとの写像
`(e i).inv ≫ X.homOfLE (…) : Z i ⟶ U`と、`TopologicalSpace.Opens.
mem_iSup`(点がある`X.basicOpen (f i)`に属することの言い換え)から
得られる全射性を渡すだけ。`corrHypGlueData.glued ≅ U`
(ロードマップ項目(b)、未完成)を示す比較対象として使う——`mathlib`の
`Scheme.Cover.gluedCover`+`Scheme.Cover.fromGlued`(`IsIso`)がこの
`OpenCover`から`U`への同型を既製で与える。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piecesOpenCover {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) : (U : Scheme).OpenCover :=
  Scheme.Cover.mkOfCovers J Z
    (fun i => (e i).inv ≫ X.homOfLE (basicOpen_le_of_iSup_eq f hcover i))
    (fun x => by
      have hx : x.1 ∈ ⨆ i, X.basicOpen (f i) := by rw [hcover]; exact x.2
      obtain ⟨i, hi⟩ := TopologicalSpace.Opens.mem_iSup.mp hx
      refine ⟨i, (e i).hom ⟨x.1, hi⟩, ?_⟩
      show ((e i).hom ≫ (e i).inv ≫ X.homOfLE (basicOpen_le_of_iSup_eq f hcover i)) ⟨x.1, hi⟩ = x
      rw [← Category.assoc, Iso.hom_inv_id, Category.id_comp]
      apply Subtype.ext
      simp)

/-! ### 項目(b)の核心部品——`V`成分の比較同型

`piecesOpenCover f Z e hcover`の`gluedCover`(mathlib既製、`.glued ≅ U`
が`Scheme.Cover.fromGlued`の`IsIso`性からすでに得られる)と
`corrHypGlueData f Z e`(独自配線)の`.diagram`(`MultispanIndex`)同士の
`NatIso`を構成するには、`U`成分(恒等——両者とも`Z`をそのままobjectに
使う)と`V`成分(`pullback (𝒰.f i)(𝒰.f j) ≅ gdV f Z e (i,j)`)の比較が
要る。ここでは`V`成分の比較同型を組み立てる——3段の`calc`:
`pullback (𝒰.f i)(𝒰.f j) ≅ pullback (X.homOfLE hi)(X.homOfLE hj)`
(`e i`・`e j`をpullbackの脚から追い出す、`pullbackHomIsoLeft`+
`pullbackSymmetry`)→`≅ (X.basicOpen (f i) ⊓ X.basicOpen (f j) :
X.Opens)`(`U`への埋め込みを経由するpullbackは`⊓`に一致する——
`isPullback_opens_inf`の`U`版、`pullbackIsPullbackOfCompMono`
(mathlib、mono後合成してもpullbackが変わらない)+一意性で得る)
→`≅ X.basicOpen (f i * f j)`(`Scheme.basicOpen_mul`)→
`≅ (Z i).basicOpen (transitionElem (f i)(f j)(e i)) = gdV f Z e (i,j)`
(`transitionElemIso`、既存)。 -/

/-- `A,B ≤ U`のとき、`X.homOfLE`同士のpullbackは`A⊓B`に一致する——
`isPullback_opens_inf`の`U`相対版。`pullbackIsPullbackOfCompMono`
(mathlib、`Mono i`のとき`pullback f g`は`pullback (f≫i)(g≫i)`の
極限錐でもある)を`U.ι`(常にmono、開埋め込み)に適用し、
`Scheme.homOfLE_ι`(`X.homOfLE e ≫ V.ι = U.ι`)で`.ι`へ書き換えてから
`isPullback_opens_inf`と一意性(`IsPullback.isoIsPullback`)で比較する。
CorrHyp非依存の一般的事実。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def pullbackHomOfLEIso {X : Scheme} {U A B : X.Opens} (hA : A ≤ U) (hB : B ≤ U) :
    (pullback (X.homOfLE hA) (X.homOfLE hB) : Scheme) ≅ (A ⊓ B : X.Opens) :=
  have h1 : IsPullback (pullback.fst (X.homOfLE hA) (X.homOfLE hB)) (pullback.snd (X.homOfLE hA) (X.homOfLE hB))
      (X.homOfLE hA ≫ U.ι) (X.homOfLE hB ≫ U.ι) :=
    IsPullback.of_isLimit (pullbackIsPullbackOfCompMono (X.homOfLE hA) (X.homOfLE hB) U.ι)
  have h2 : IsPullback (pullback.fst (X.homOfLE hA) (X.homOfLE hB)) (pullback.snd (X.homOfLE hA) (X.homOfLE hB))
      A.ι B.ι := (Scheme.homOfLE_ι X hA) ▸ (Scheme.homOfLE_ι X hB) ▸ h1
  h2.isoIsPullback (A : Scheme) (B : Scheme) (isPullback_opens_inf A B)

/-- `pullbackHomOfLEIso`の`fst`との可換性——`IsPullback.isoIsPullback_hom_fst`
をそのまま呼ぶだけ。 -/
theorem pullbackHomOfLEIso_hom_fst {X : Scheme} {U A B : X.Opens} (hA : A ≤ U) (hB : B ≤ U) :
    (pullbackHomOfLEIso hA hB).hom ≫ X.homOfLE (inf_le_left : A ⊓ B ≤ A) =
      pullback.fst (X.homOfLE hA) (X.homOfLE hB) := by
  unfold pullbackHomOfLEIso
  exact IsPullback.isoIsPullback_hom_fst _ _ _ _

/-- `pullbackHomOfLEIso`の`snd`との可換性——同様。 -/
theorem pullbackHomOfLEIso_hom_snd {X : Scheme} {U A B : X.Opens} (hA : A ≤ U) (hB : B ≤ U) :
    (pullbackHomOfLEIso hA hB).hom ≫ X.homOfLE (inf_le_right : A ⊓ B ≤ B) =
      pullback.snd (X.homOfLE hA) (X.homOfLE hB) := by
  unfold pullbackHomOfLEIso
  exact IsPullback.isoIsPullback_hom_snd _ _ _ _

/-- `pullbackHomOfLEIso`を`Scheme.basicOpen_mul`で単一の基本開集合の形
(`X.basicOpen (f₁*f₂)`)にまとめたもの。 -/
noncomputable def pullbackHomOfLEIsoBasicOpen {X : Scheme} {U : X.Opens} (f₁ f₂ : Γ(X, U))
    (h1 : X.basicOpen f₁ ≤ U) (h2 : X.basicOpen f₂ ≤ U) :
    (pullback (X.homOfLE h1) (X.homOfLE h2) : Scheme) ≅ (X.basicOpen (f₁ * f₂) : Scheme) :=
  (pullbackHomOfLEIso h1 h2).trans (X.isoOfEq (X.basicOpen_mul f₁ f₂).symm)

/-- `pullbackHomOfLEIsoBasicOpen`の`fst`との可換性——`Scheme.isoOfEq_hom`
(`X.isoOfEq e`の`.hom`は`X.homOfLE`そのもの)+`Scheme.homOfLE_homOfLE`
(合成が単一の`homOfLE`にまとまる)で`pullbackHomOfLEIso_hom_fst`へ帰着。 -/
theorem pullbackHomOfLEIsoBasicOpen_hom_fst {X : Scheme} {U : X.Opens} (f₁ f₂ : Γ(X, U))
    (h1 : X.basicOpen f₁ ≤ U) (h2 : X.basicOpen f₂ ≤ U) :
    (pullbackHomOfLEIsoBasicOpen f₁ f₂ h1 h2).hom ≫
      X.homOfLE (show X.basicOpen (f₁ * f₂) ≤ X.basicOpen f₁ by
        rw [X.basicOpen_mul f₁ f₂]; exact inf_le_left) =
      pullback.fst (X.homOfLE h1) (X.homOfLE h2) := by
  unfold pullbackHomOfLEIsoBasicOpen
  rw [Iso.trans_hom, Category.assoc, Scheme.isoOfEq_hom, Scheme.homOfLE_homOfLE]
  exact pullbackHomOfLEIso_hom_fst h1 h2

/-- `pullbackHomOfLEIsoBasicOpen`の`snd`との可換性——同様。 -/
theorem pullbackHomOfLEIsoBasicOpen_hom_snd {X : Scheme} {U : X.Opens} (f₁ f₂ : Γ(X, U))
    (h1 : X.basicOpen f₁ ≤ U) (h2 : X.basicOpen f₂ ≤ U) :
    (pullbackHomOfLEIsoBasicOpen f₁ f₂ h1 h2).hom ≫
      X.homOfLE (show X.basicOpen (f₁ * f₂) ≤ X.basicOpen f₂ by
        rw [X.basicOpen_mul f₁ f₂]; exact inf_le_right) =
      pullback.snd (X.homOfLE h1) (X.homOfLE h2) := by
  unfold pullbackHomOfLEIsoBasicOpen
  rw [Iso.trans_hom, Category.assoc, Scheme.isoOfEq_hom, Scheme.homOfLE_homOfLE]
  exact pullbackHomOfLEIso_hom_snd h1 h2

/-! ### 項目(b)の残り——`φ(i,j)`の自然性(`fst`/`snd`との可換性)

`piecesGluedCoverVIso`(`V`成分比較同型)がNatIsoの1成分として使える
ためには、`fst`/`snd`(`corrHypGlueData`側の`gdF`・`gdT≫gdF`、
`piecesOpenCover.gluedCover`側の`pullback.fst`・
`pullbackSymmetry.hom≫pullback.fst`)と可換であることを示す必要がある。
`φ(i,j)`は`pullbackEInvIso`(`e`を追い出す)→`pullbackHomOfLEIsoBasicOpen`
(`⊓`を`basicOpen`にまとめる)→`transitionElemIso`(`X`側から`Z i`側へ
転送)の3段構成なので、各段の自然性を個別に示してから繋ぐ。ここでは
まず`pullbackHomIsoLeft`(mathlibに無かったので既存で補った一般的事実)
自身の`fst`/`snd`自然性(`pullback.map`の定義から`pullback.lift_fst`/
`_snd`で直ちに従う)を確立する。 -/

/-- `pullbackHomIsoLeft`の`fst`との可換性——`pullback.map`の定義
(`pullback.lift`)から`pullback.lift_fst`で直ちに従う。`@[reassoc]`
を付け、`pullbackEInvIso`の自然性の証明で「末尾に何かが続く」形の
書き換えにも使えるようにする。 -/
@[reassoc]
theorem pullbackHomIsoLeft_hom_fst' {X Y Z W : Scheme} (i : X ≅ Y) (g : Y ⟶ Z) (f : W ⟶ Z) :
    (pullbackHomIsoLeft i g f).hom ≫ pullback.fst g f = pullback.fst (i.hom ≫ g) f ≫ i.hom := by
  unfold pullbackHomIsoLeft
  exact pullback.lift_fst _ _ _

/-- `pullbackHomIsoLeft`の`snd`との可換性——同様。 -/
@[reassoc]
theorem pullbackHomIsoLeft_hom_snd' {X Y Z W : Scheme} (i : X ≅ Y) (g : Y ⟶ Z) (f : W ⟶ Z) :
    (pullbackHomIsoLeft i g f).hom ≫ pullback.snd g f = pullback.snd (i.hom ≫ g) f := by
  unfold pullbackHomIsoLeft pullback.map
  exact pullback.lift_snd _ _ _

/-- `piecesOpenCover`の脚`(e i).inv ≫ X.homOfLE (h i)`同士のpullbackは、
`e i`・`e j`をpullbackの脚から追い出す(`pullbackHomIsoLeft`+
`pullbackSymmetry`、いずれも既存の一般的事実)ことで、`X.homOfLE`同士
のpullbackへ帰着する。 -/
noncomputable def pullbackEInvIso {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (h : ∀ i, X.basicOpen (f i) ≤ U) (i j : J) :
    (pullback ((e i).inv ≫ X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j)) : Scheme) ≅
      (pullback (X.homOfLE (h i)) (X.homOfLE (h j)) : Scheme) :=
  (pullbackHomIsoLeft (e i).symm (X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j))).trans
    ((pullbackSymmetry (X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j))).trans
      ((pullbackHomIsoLeft (e j).symm (X.homOfLE (h j)) (X.homOfLE (h i))).trans
        (pullbackSymmetry (X.homOfLE (h j)) (X.homOfLE (h i)))))

/-- **`pullbackEInvIso`の`fst`との可換性**——`(e i).symm.hom`(`show`で
明示形へ)に対し、内側から順に`pullbackSymmetry_hom_comp_fst`・
`pullbackHomIsoLeft_hom_snd'`・`pullbackSymmetry_hom_comp_snd`・
`pullbackHomIsoLeft_hom_fst'`(すべて非`_assoc`版、末尾に何も続かない
ので)を`rw`で辿るだけ。

★配管の注意(`#31`の新しい現れ方): `(e i).inv`をそのまま`show`に
書くと、`pullbackHomIsoLeft`が内部で使う`(e i).symm.hom`との型不一致
(`instances`透明度では見分けられない)で`rw`が失敗する——`show`の
中では`(e i).inv`ではなく`(e i).symm.hom`で統一して書くことで解消。

★**sorry 無し**。標準3公理のみ。 -/
theorem pullbackEInvIso_hom_fst {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (h : ∀ i, X.basicOpen (f i) ≤ U) (i j : J) :
    (pullbackEInvIso f Z e h i j).hom ≫ pullback.fst (X.homOfLE (h i)) (X.homOfLE (h j)) =
      pullback.fst ((e i).inv ≫ X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j)) ≫ (e i).inv := by
  show (pullbackHomIsoLeft (e i).symm (X.homOfLE (h i)) ((e j).symm.hom ≫ X.homOfLE (h j))).hom ≫
      (pullbackSymmetry (X.homOfLE (h i)) ((e j).symm.hom ≫ X.homOfLE (h j))).hom ≫
      (pullbackHomIsoLeft (e j).symm (X.homOfLE (h j)) (X.homOfLE (h i))).hom ≫
      (pullbackSymmetry (X.homOfLE (h j)) (X.homOfLE (h i))).hom ≫
      pullback.fst (X.homOfLE (h i)) (X.homOfLE (h j)) =
    pullback.fst ((e i).symm.hom ≫ X.homOfLE (h i)) ((e j).symm.hom ≫ X.homOfLE (h j)) ≫ (e i).symm.hom
  rw [pullbackSymmetry_hom_comp_fst, pullbackHomIsoLeft_hom_snd',
    pullbackSymmetry_hom_comp_snd, pullbackHomIsoLeft_hom_fst']

/-- **`pullbackEInvIso`の`snd`との可換性**——`fst`版と対称、ただし
「末尾に`(e j).symm.hom`が続く」形になる2段は`_assoc`版を使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem pullbackEInvIso_hom_snd {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (h : ∀ i, X.basicOpen (f i) ≤ U) (i j : J) :
    (pullbackEInvIso f Z e h i j).hom ≫ pullback.snd (X.homOfLE (h i)) (X.homOfLE (h j)) =
      pullback.snd ((e i).inv ≫ X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j)) ≫ (e j).inv := by
  show (pullbackHomIsoLeft (e i).symm (X.homOfLE (h i)) ((e j).symm.hom ≫ X.homOfLE (h j))).hom ≫
      (pullbackSymmetry (X.homOfLE (h i)) ((e j).symm.hom ≫ X.homOfLE (h j))).hom ≫
      (pullbackHomIsoLeft (e j).symm (X.homOfLE (h j)) (X.homOfLE (h i))).hom ≫
      (pullbackSymmetry (X.homOfLE (h j)) (X.homOfLE (h i))).hom ≫
      pullback.snd (X.homOfLE (h i)) (X.homOfLE (h j)) =
    pullback.snd ((e i).symm.hom ≫ X.homOfLE (h i)) ((e j).symm.hom ≫ X.homOfLE (h j)) ≫ (e j).symm.hom
  rw [pullbackSymmetry_hom_comp_snd, pullbackHomIsoLeft_hom_fst',
    pullbackSymmetry_hom_comp_fst_assoc, pullbackHomIsoLeft_hom_snd'_assoc]

/-- **`V`成分の比較同型(`h`を任意に取った、`piecesOpenCover`非依存の
中間形)**——`pullbackEInvIso`(`e`を追い出す)+`pullbackHomOfLEIsoBasicOpen`
(`⊓`を単一の`basicOpen`にまとめる)+`transitionElemIso`(`X`側から`Z i`
側への転送、既存)を繋ぐだけ。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def pullbackHomOfLE_gdV_iso {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (h : ∀ i, X.basicOpen (f i) ≤ U) (i j : J) :
    (pullback ((e i).inv ≫ X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j)) : Scheme) ≅
      gdV f Z e (i, j) :=
  (pullbackEInvIso f Z e h i j).trans
    ((pullbackHomOfLEIsoBasicOpen (f i) (f j) (h i) (h j)).trans (transitionElemIso (f i) (f j) (e i)))

/-- **`pullbackHomOfLE_gdV_iso`の`fst`との可換性**——3段
(`pullbackEInvIso`・`pullbackHomOfLEIsoBasicOpen`・`transitionElemIso`)
それぞれの`fst`/`.ι`自然性(`pullbackEInvIso_hom_fst`・
`pullbackHomOfLEIsoBasicOpen_hom_fst`・`transitionElemIso_hom_ι`)を
内側から順に繋ぐ——「末尾に何も続かない」箇所は非`_assoc`版で
書き換えてから、最後に`(e i).inv≫(e i).hom=𝟙`で閉じる。

★**sorry 無し**。標準3公理のみ。 -/
theorem pullbackHomOfLE_gdV_iso_hom_fst {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (h : ∀ i, X.basicOpen (f i) ≤ U) (i j : J) :
    (pullbackHomOfLE_gdV_iso f Z e h i j).hom ≫ gdF f Z e i j =
      pullback.fst ((e i).inv ≫ X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j)) := by
  show (pullbackEInvIso f Z e h i j).hom ≫ (pullbackHomOfLEIsoBasicOpen (f i) (f j) (h i) (h j)).hom ≫
      (transitionElemIso (f i) (f j) (e i)).hom ≫ gdF f Z e i j = _
  unfold gdF
  rw [transitionElemIso_hom_ι]
  rw [show (pullbackHomOfLEIsoBasicOpen (f i) (f j) (h i) (h j)).hom ≫
      (X.homOfLE (show X.basicOpen (f i * f j) ≤ X.basicOpen (f i) by
        rw [X.basicOpen_mul (f i) (f j)]; exact inf_le_left) ≫ (e i).hom)
    = pullback.fst (X.homOfLE (h i)) (X.homOfLE (h j)) ≫ (e i).hom from by
      rw [← Category.assoc, pullbackHomOfLEIsoBasicOpen_hom_fst]]
  rw [← Category.assoc, pullbackEInvIso_hom_fst, Category.assoc, Iso.inv_hom_id, Category.comp_id]

/-- **`pullbackHomOfLE_gdV_iso`の`snd`との可換性**——`fst`版と対称、
`gdF f Z e i j`の代わりに`(gdT f Z e i j).hom ≫ gdF f Z e j i`
(`GlueData.diagram`の`snd`フィールドそのもの)を使う。`gdT_hom_comp_gdF`
で`transitionElemIso(f i)(f j)(e i))`の`.inv`まで帰着させてから、
同じ`transitionElemIso(f i)(f j)(e i))`の`.hom`(`φ`の最終段)と
`Iso.hom_inv_id_assoc`で相殺するのが`fst`版との違い。

★**sorry 無し**。標準3公理のみ。 -/
theorem pullbackHomOfLE_gdV_iso_hom_snd {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (h : ∀ i, X.basicOpen (f i) ≤ U) (i j : J) :
    (pullbackHomOfLE_gdV_iso f Z e h i j).hom ≫ (gdT f Z e i j).hom ≫ gdF f Z e j i =
      pullback.snd ((e i).inv ≫ X.homOfLE (h i)) ((e j).inv ≫ X.homOfLE (h j)) := by
  show (pullbackEInvIso f Z e h i j).hom ≫ (pullbackHomOfLEIsoBasicOpen (f i) (f j) (h i) (h j)).hom ≫
      (transitionElemIso (f i) (f j) (e i)).hom ≫ (gdT f Z e i j).hom ≫ gdF f Z e j i = _
  rw [gdT_hom_comp_gdF, Iso.hom_inv_id_assoc]
  rw [show (pullbackHomOfLEIsoBasicOpen (f i) (f j) (h i) (h j)).hom ≫
      (X.homOfLE (show X.basicOpen (f i * f j) ≤ X.basicOpen (f j) by
        rw [X.basicOpen_mul (f i) (f j)]; exact inf_le_right) ≫ (e j).hom)
    = pullback.snd (X.homOfLE (h i)) (X.homOfLE (h j)) ≫ (e j).hom from by
      rw [← Category.assoc, pullbackHomOfLEIsoBasicOpen_hom_snd]]
  rw [← Category.assoc, pullbackEInvIso_hom_snd, Category.assoc, Iso.inv_hom_id, Category.comp_id]

/-- `piecesOpenCover`の`.f i`が定義どおりの明示形であること(単一
出現の`unfold`なので軽い、`#31`/`#32`の教訓)。 -/
theorem piecesOpenCover_f_eq {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) (i : J) :
    (piecesOpenCover f Z e hcover).f i = (e i).inv ≫ X.homOfLE (basicOpen_le_of_iSup_eq f hcover i) := by
  unfold piecesOpenCover
  rfl

/-- **`piecesOpenCover`の`gluedCover`(mathlib)の`V`成分と
`corrHypGlueData`の`V`成分(`gdV`)を比較する同型**——項目(b)のNatIso
構成に必要な2つの成分(`U`成分は恒等、`V`成分がこれ)のうち、より
本質的な方。`piecesOpenCover_f_eq`で`.f i`を明示形へ書き換えてから
`pullbackHomOfLE_gdV_iso`を適用するだけ。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def piecesGluedCoverVIso {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) (i j : J) :
    (pullback ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) : Scheme) ≅
      gdV f Z e (i, j) := by
  rw [piecesOpenCover_f_eq f Z e hcover i, piecesOpenCover_f_eq f Z e hcover j]
  exact pullbackHomOfLE_gdV_iso f Z e (basicOpen_le_of_iSup_eq f hcover) i j

/-- **`piecesGluedCoverVIso`の`fst`との可換性(NatIso自然性の`fst`半分、
完成)**——`pullbackHomOfLE_gdV_iso_hom_fst`と`piecesOpenCover_f_eq`
(`rfl`による定義的一致なので`rw`を経由せず`show`+`exact`だけで届く、
`#31`の「`rw`を避けdefeqで直接繋ぐ」教訓の別の現れ方)を合成するだけ。

★**sorry 無し**。標準3公理のみ。 -/
theorem piecesGluedCoverVIso_hom_fst {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) (i j : J) :
    (piecesGluedCoverVIso f Z e hcover i j).hom ≫ gdF f Z e i j =
      pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) := by
  show (pullbackHomOfLE_gdV_iso f Z e (basicOpen_le_of_iSup_eq f hcover) i j).hom ≫ gdF f Z e i j = _
  exact pullbackHomOfLE_gdV_iso_hom_fst f Z e (basicOpen_le_of_iSup_eq f hcover) i j

/-- **`piecesGluedCoverVIso`の`snd`との可換性(NatIso自然性の`snd`半分、
完成)**——`pullbackHomOfLE_gdV_iso_hom_snd`+`pullbackSymmetry_hom_comp_
fst`(mathlib既製、`(pullbackSymmetry f g).hom≫pullback.fst g f=
pullback.snd f g`)を合成、`piecesOpenCover_f_eq`への接続は最後の`rfl`
(定義的一致)で閉じる。**これで`fst`・`snd`ともに揃い、項目(b)のNatIso
自然性が完全に確立した**。

★**sorry 無し**。標準3公理のみ。 -/
theorem piecesGluedCoverVIso_hom_snd {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) (i j : J) :
    (piecesGluedCoverVIso f Z e hcover i j).hom ≫ (gdT f Z e i j).hom ≫ gdF f Z e j i =
      (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
        pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i) := by
  show (pullbackHomOfLE_gdV_iso f Z e (basicOpen_le_of_iSup_eq f hcover) i j).hom ≫
      (gdT f Z e i j).hom ≫ gdF f Z e j i = _
  rw [pullbackHomOfLE_gdV_iso_hom_snd, pullbackSymmetry_hom_comp_fst]
  rfl

/-- **`corrHypGlueData`の`t≫f`(`GlueData.diagram`の`snd`)コクイザライザ
条件が`piecesOpenCover`の`gluedCover`のそれと両立する**——
`piecesGluedCoverVIso_hom_fst`/`_hom_snd`(`fst`/`snd`自然性)と
`Scheme.Cover.gluedCover`自身の`GlueData.glue_condition`(mathlib既製、
コクイザライザの定義関係式)を組み合わせるだけ。`Multicoequalizer.desc`
で`corrHypGlueData.glued`から`gluedCover.glued`への射を作るのに必要な
唯一の条件。

★配管の注意(`#31`の新しい現れ方、これまでで最も執拗): `𝒢.ι i`・
`𝒢.f i j`のような、`gluedCover`(`piecesOpenCover`から作った
`Scheme.GlueData`)の`J`成分が`piecesOpenCover`自身の`I₀`(`=J`、
定義的には等しい)と`instances`透明度で一致しないため、`rw`は
おろか`simp only […] at h`のような**ローカル仮定への書き換えでも**
同じ壁に当たる(`set 𝒢 := …`で名前を付けても、`letI`で型注釈しても
解消しない場合があった)。**唯一有効だった方法**: 中間の`have`を
一切`rw`/`simp`で後から書き換えず、**最初から`calc`+`congrArg`+
`Category.assoc`(項として)だけで組み立てる**——`piecesGluedCoverVIso_
hom_fst`/`_hom_snd`から得た事実を`Iso.inv_hom_id_assoc`と`.trans`で
直接合成し、最終的な`calc`チェーンも`congrArg (· ≫ …)`/
`Category.assoc _ _ _`(`.symm`込み)だけで積み上げる。これで一度も
`rw`/`simp`を使わずに完走できた。

★**sorry 無し**。標準3公理のみ。 -/
theorem corrHypGlueData_compat {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) (i j : J) :
    gdF f Z e i j ≫ (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i =
      ((gdT f Z e i j).hom ≫ gdF f Z e j i) ≫
        (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι j := by
  have hglue :
      (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
        pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i) ≫
        (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι j =
      pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) ≫
        (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i :=
    (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).glue_condition i j
  have hgdF : gdF f Z e i j = (piecesGluedCoverVIso f Z e hcover i j).inv ≫
      pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) :=
    (Iso.inv_hom_id_assoc (piecesGluedCoverVIso f Z e hcover i j) (gdF f Z e i j)).symm.trans
      (congrArg ((piecesGluedCoverVIso f Z e hcover i j).inv ≫ ·) (piecesGluedCoverVIso_hom_fst f Z e hcover i j))
  have hgdTF : (gdT f Z e i j).hom ≫ gdF f Z e j i =
      (piecesGluedCoverVIso f Z e hcover i j).inv ≫
      (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
      pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i) :=
    (Iso.inv_hom_id_assoc (piecesGluedCoverVIso f Z e hcover i j)
      ((gdT f Z e i j).hom ≫ gdF f Z e j i)).symm.trans
      (congrArg ((piecesGluedCoverVIso f Z e hcover i j).inv ≫ ·) (piecesGluedCoverVIso_hom_snd f Z e hcover i j))
  calc gdF f Z e i j ≫ (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i
      = ((piecesGluedCoverVIso f Z e hcover i j).inv ≫
          pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)) ≫
          (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i :=
        congrArg (· ≫ (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i) hgdF
    _ = (piecesGluedCoverVIso f Z e hcover i j).inv ≫
          pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) ≫
          (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i :=
        Category.assoc _ _ _
    _ = (piecesGluedCoverVIso f Z e hcover i j).inv ≫
          (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
          pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i) ≫
          (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι j :=
        congrArg ((piecesGluedCoverVIso f Z e hcover i j).inv ≫ ·) hglue.symm
    _ = ((piecesGluedCoverVIso f Z e hcover i j).inv ≫
          (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
          pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i)) ≫
          (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι j :=
        (Category.assoc _ _ _).symm
    _ = ((gdT f Z e i j).hom ≫ gdF f Z e j i) ≫
          (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι j :=
        congrArg (· ≫ (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι j) hgdTF.symm

/-- **`corrHypGlueData`から`piecesOpenCover`の`gluedCover`への射**——
`Multicoequalizer.desc`(コクイザライザの普遍性)を`corrHypGlueData_
compat`(唯一必要な整合条件)とともに適用するだけ。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypGlueData_toGluedCover {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) :
    (corrHypGlueData f Z e).glued ⟶ (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).glued :=
  Multicoequalizer.desc (corrHypGlueData f Z e).diagram _
    (fun i => (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i)
    (fun ⟨i, j⟩ => corrHypGlueData_compat f Z e hcover i j)

/-- `corrHypGlueData_toGluedCover`の逆向きの整合条件——
`corrHypGlueData_compat`と同じ`piecesGluedCoverVIso_hom_fst`/`_hom_snd`
から、今度は`corrHypGlueData`自身の`GlueData.glue_condition`(mathlib
既製)を経由して導く。`φ.inv`を消して`gdF`/`gdT`絡みの言葉から
`pullback.fst`/`pullbackSymmetry`絡みの言葉へ変換したのち、`φ`が
iso であること(`Iso.hom_inv_id`)で`φ`自体を打ち消す、という
`corrHypGlueData_compat`とは逆方向のキャンセルが必要になる。

★配管の注意: `corrHypGlueData_compat`と同じく、一度も`rw`/`simp`を
使わず`calc`+`congrArg`+`Category.assoc`(項として)だけで組み立てた
(`#31`の壁の回避法、再確認)。

★**sorry 無し**。標準3公理のみ。 -/
theorem gluedCover_compat {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) (i j : J) :
    pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) ≫
        (corrHypGlueData f Z e).ι i =
      ((pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
        pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i)) ≫
        (corrHypGlueData f Z e).ι j := by
  have hglue2 : (gdT f Z e i j).hom ≫ gdF f Z e j i ≫ (corrHypGlueData f Z e).ι j =
      gdF f Z e i j ≫ (corrHypGlueData f Z e).ι i :=
    (corrHypGlueData f Z e).glue_condition i j
  have hgdF : gdF f Z e i j = (piecesGluedCoverVIso f Z e hcover i j).inv ≫
      pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) :=
    (Iso.inv_hom_id_assoc (piecesGluedCoverVIso f Z e hcover i j) (gdF f Z e i j)).symm.trans
      (congrArg ((piecesGluedCoverVIso f Z e hcover i j).inv ≫ ·) (piecesGluedCoverVIso_hom_fst f Z e hcover i j))
  have hgdTF : (gdT f Z e i j).hom ≫ gdF f Z e j i =
      (piecesGluedCoverVIso f Z e hcover i j).inv ≫
      (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
      pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i) :=
    (Iso.inv_hom_id_assoc (piecesGluedCoverVIso f Z e hcover i j)
      ((gdT f Z e i j).hom ≫ gdF f Z e j i)).symm.trans
      (congrArg ((piecesGluedCoverVIso f Z e hcover i j).inv ≫ ·) (piecesGluedCoverVIso_hom_snd f Z e hcover i j))
  have hcancel : (piecesGluedCoverVIso f Z e hcover i j).inv ≫
      (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
      pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i) ≫
      (corrHypGlueData f Z e).ι j =
      (piecesGluedCoverVIso f Z e hcover i j).inv ≫
      pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) ≫
      (corrHypGlueData f Z e).ι i :=
    calc (piecesGluedCoverVIso f Z e hcover i j).inv ≫
        (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
        pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i) ≫
        (corrHypGlueData f Z e).ι j
        = ((piecesGluedCoverVIso f Z e hcover i j).inv ≫
            (pullbackSymmetry ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)).hom ≫
            pullback.fst ((piecesOpenCover f Z e hcover).f j) ((piecesOpenCover f Z e hcover).f i)) ≫
            (corrHypGlueData f Z e).ι j := (Category.assoc _ _ _).symm
      _ = ((gdT f Z e i j).hom ≫ gdF f Z e j i) ≫ (corrHypGlueData f Z e).ι j :=
          congrArg (· ≫ (corrHypGlueData f Z e).ι j) hgdTF.symm
      _ = (gdT f Z e i j).hom ≫ gdF f Z e j i ≫ (corrHypGlueData f Z e).ι j := Category.assoc _ _ _
      _ = gdF f Z e i j ≫ (corrHypGlueData f Z e).ι i := hglue2
      _ = ((piecesGluedCoverVIso f Z e hcover i j).inv ≫
            pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j)) ≫
            (corrHypGlueData f Z e).ι i :=
          congrArg (· ≫ (corrHypGlueData f Z e).ι i) hgdF
      _ = (piecesGluedCoverVIso f Z e hcover i j).inv ≫
            pullback.fst ((piecesOpenCover f Z e hcover).f i) ((piecesOpenCover f Z e hcover).f j) ≫
            (corrHypGlueData f Z e).ι i := Category.assoc _ _ _
  have hfinal := congrArg ((piecesGluedCoverVIso f Z e hcover i j).hom ≫ ·) hcancel
  rw [← Category.assoc, ← Category.assoc, Iso.hom_inv_id, Category.id_comp,
    ← Category.assoc, ← Category.assoc, Iso.hom_inv_id, Category.id_comp] at hfinal
  exact hfinal.symm

/-- **`piecesOpenCover`の`gluedCover`から`corrHypGlueData`への射(逆
向き)**——`gluedCover_compat`とともに`Multicoequalizer.desc`を適用
するだけ。 -/
noncomputable def gluedCover_toCorrHypGlueData {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) :
    (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).glued ⟶ (corrHypGlueData f Z e).glued :=
  Multicoequalizer.desc (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).diagram _
    (fun i => (corrHypGlueData f Z e).ι i)
    (fun ⟨i, j⟩ => gluedCover_compat f Z e hcover i j)

/-- `corrHypGlueData_toGluedCover`・`gluedCover_toCorrHypGlueData`が
互いに逆であること(1本目)——`Multicoequalizer.hom_ext`(コクイザライザ
からの射は各ピースへの制限で決まる)+両方の`Multicoequalizer.π_desc`
(定義そのもの)を組み合わせるだけ。 -/
theorem corrHypGlueData_toGluedCover_comp {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) :
    corrHypGlueData_toGluedCover f Z e hcover ≫ gluedCover_toCorrHypGlueData f Z e hcover = 𝟙 _ := by
  apply Multicoequalizer.hom_ext
  intro i
  have step1 : (Multicoequalizer.π (corrHypGlueData f Z e).diagram i ≫ corrHypGlueData_toGluedCover f Z e hcover) ≫
      gluedCover_toCorrHypGlueData f Z e hcover =
      (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i ≫
        gluedCover_toCorrHypGlueData f Z e hcover :=
    congrArg (· ≫ gluedCover_toCorrHypGlueData f Z e hcover)
      (Multicoequalizer.π_desc (corrHypGlueData f Z e).diagram _
        (fun i => (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i)
        (fun ⟨i, j⟩ => corrHypGlueData_compat f Z e hcover i j) i)
  have step2 : (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i ≫
      gluedCover_toCorrHypGlueData f Z e hcover = (corrHypGlueData f Z e).ι i :=
    Multicoequalizer.π_desc (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).diagram _
      (fun i => (corrHypGlueData f Z e).ι i) (fun ⟨i, j⟩ => gluedCover_compat f Z e hcover i j) i
  have hgoal : Multicoequalizer.π (corrHypGlueData f Z e).diagram i ≫
      (corrHypGlueData_toGluedCover f Z e hcover ≫ gluedCover_toCorrHypGlueData f Z e hcover) =
      (corrHypGlueData f Z e).ι i :=
    ((Category.assoc _ _ _).symm.trans step1).trans step2
  exact hgoal.trans (Category.comp_id _).symm

/-- 互いに逆であること(2本目、反対向き)——同じ構成の対称形。 -/
theorem gluedCover_toCorrHypGlueData_comp {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) :
    gluedCover_toCorrHypGlueData f Z e hcover ≫ corrHypGlueData_toGluedCover f Z e hcover = 𝟙 _ := by
  apply Multicoequalizer.hom_ext
  intro i
  have step1 : (Multicoequalizer.π (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).diagram i ≫
        gluedCover_toCorrHypGlueData f Z e hcover) ≫ corrHypGlueData_toGluedCover f Z e hcover =
      (corrHypGlueData f Z e).ι i ≫ corrHypGlueData_toGluedCover f Z e hcover :=
    congrArg (· ≫ corrHypGlueData_toGluedCover f Z e hcover)
      (Multicoequalizer.π_desc (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).diagram _
        (fun i => (corrHypGlueData f Z e).ι i) (fun ⟨i, j⟩ => gluedCover_compat f Z e hcover i j) i)
  have step2 : (corrHypGlueData f Z e).ι i ≫ corrHypGlueData_toGluedCover f Z e hcover =
      (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i :=
    Multicoequalizer.π_desc (corrHypGlueData f Z e).diagram _
      (fun i => (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i)
      (fun ⟨i, j⟩ => corrHypGlueData_compat f Z e hcover i j) i
  have hgoal : Multicoequalizer.π (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).diagram i ≫
      (gluedCover_toCorrHypGlueData f Z e hcover ≫ corrHypGlueData_toGluedCover f Z e hcover) =
      (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).ι i :=
    ((Category.assoc _ _ _).symm.trans step1).trans step2
  exact hgoal.trans (Category.comp_id _).symm

/-- **`corrHypGlueData.glued ≅ (piecesOpenCover ...).gluedCover.glued`**
——`corrHypGlueData_toGluedCover`/`gluedCover_toCorrHypGlueData`が
互いに逆であることから`Iso.mk`でまとめるだけ。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypGlueData_gluedIso {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) :
    (corrHypGlueData f Z e).glued ≅ (Scheme.Cover.gluedCover (piecesOpenCover f Z e hcover)).glued where
  hom := corrHypGlueData_toGluedCover f Z e hcover
  inv := gluedCover_toCorrHypGlueData f Z e hcover
  hom_inv_id := corrHypGlueData_toGluedCover_comp f Z e hcover
  inv_hom_id := gluedCover_toCorrHypGlueData_comp f Z e hcover

/-- **★★★ロードマップ項目(b)の完成★★★: `corrHypGlueData.glued ≅ U`**
——`corrHypGlueData_gluedIso`(`corrHypGlueData.glued ≅ gluedCover.glued`、
このファイルで構成)と`Scheme.Cover.fromGlued`(mathlib既製、
`gluedCover.glued ≅ U`、`IsIso`性を`asIso`で取り出す)を`.trans`で
繋ぐだけ。`Z`の族(標準的な有限étale候補片)から作った`corrHypGlueData`
の貼り合わせ結果が、実際に元の`U`に同型であることを示す——`Lemma 4.1`
の構成的降下で「候補片の族が実際に`C`(の対応する開集合)を再構成する」
ことを保証する、GlueData配線の最終目標。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypGlueData_glued_iso {X : Scheme} {U : X.Opens} {J : Type} (f : J → Γ(X, U))
    (Z : J → Scheme) (e : ∀ i, (X.basicOpen (f i) : Scheme) ≅ Z i)
    (hcover : ⨆ i, X.basicOpen (f i) = U) :
    (corrHypGlueData f Z e).glued ≅ (U : Scheme) :=
  (corrHypGlueData_gluedIso f Z e hcover).trans (asIso (piecesOpenCover f Z e hcover).fromGlued)

/-! ### `corrHypGlueData`の具体化(ロードマップ項目(a)の第一歩)

`corrHypGlueData`はここまで完全に抽象的(`X,U,J,f,Z,e`は任意)だった。
ここでは`piece_descends_iso`(既に存在する、単一のstandard-étale元
`f`に対する有限段階スキームの存在命題)から`.choose`/`.choose_spec`で
具体的な`Z i`・`e i`の族を取り出し、`corrHypGlueData`を実際に
呼び出せる形に落とし込む。

★重要な発見: `corrHypGlueData`の12個のフィールド(特に`t_fac`・
`cocycle`)は`Z i`・`Z j`(異なる添字)を互いに比較する必要が一切ない
——各`e i`は添字`i`ごとに独立に使われるだけなので、`.choose`由来の
「アルゴリズム的に不透明な」`Z i`同士でも問題は起きない。これは
以前`transitionElem`で「複数の`.choose`呼び出しが代数的に無関係に
なってしまう」問題に苦しんだのとは対照的(`corrhyp-goal.md`参照)。 -/

/-- `piece_descends_iso`の存在証明から、単一のstandard-étale元`f`に
対応する有限段階スキームを`.choose`で取り出したもの。 -/
noncomputable def descendPiece {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    (f : Γ(X, U)) [Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)] : Scheme :=
  letI R := (piece_descends_iso hU f).choose
  letI P₀ := (piece_descends_iso hU f).choose_spec.choose
  letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
  pullback (standardEtalePairSpecMap P₀)
    (Spec.map (CommRingCat.ofHom (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ))))

/-- `descendPiece`が`X.basicOpen f`と同型であることの比較同型
(`piece_descends_iso`の`Nonempty`成分から`.some`で取り出す)。 -/
noncomputable def descendPieceIso {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    (f : Γ(X, U)) [Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)] :
    (X.basicOpen f : Scheme) ≅ descendPiece (A := A) hU f := by
  unfold descendPiece
  exact (piece_descends_iso hU f).choose_spec.choose_spec.some

/-- `descendPiece`の instance 版(`[Algebra.IsStandardEtale ...]`を
typeclass ではなく明示的な証明として受け取る)——`Finset`の添字ごとに
`IsStandardEtale`証明を個別に渡して族を組み立てるのに使う。 -/
noncomputable def descendPieceOfProof {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    (f : Γ(X, U)) (hf : Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)) : Scheme :=
  letI := hf
  descendPiece (A := A) hU f

/-- `descendPieceIso`の instance 版。 -/
noncomputable def descendPieceIsoOfProof {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    (f : Γ(X, U)) (hf : Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away f)) :
    (X.basicOpen f : Scheme) ≅ descendPieceOfProof (A := A) hU f hf := by
  unfold descendPieceOfProof
  exact descendPieceIso (A := A) hU f

/-- **`corrHypGlueData`をCorrHypの実際のデータへ具体化する**: `Γ(X,U)`
上のstandard-étale元の有限族`f : ι → Γ(X,U)`(`Finset t`で添字づけ)が
与えられたとき、各`f i`ごとに`descendPieceOfProof`で候補片を作り、
`corrHypGlueData`をその族に適用する。

★注意: これはまだ被覆条件`⨆ i∈t, X.basicOpen (f i) = U`を使っていない
——それは`corrHypGlueData.glued ≅ U`(ロードマップ項目(b))で必要になる。
また`A`・`X`・`U`と実際の`corrHypInstance4`・`Ext`・`C`との接続も未着手
(ロードマップ項目(c))。ここは「有限standard-étale被覆→GlueData」という
CorrHyp非依存の再利用可能な部品。

★**sorry 無し**。 -/
noncomputable def corrHypGlueDataOfCover {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    {ι : Type} (t : Finset ι) (f : ι → Γ(X, U))
    (hf : ∀ i ∈ t, Algebra.IsStandardEtale (A ⊗[ℚ] ℝ) (Localization.Away (f i))) :
    Scheme.GlueData :=
  corrHypGlueData (X := X) (U := U) (J := {i // i ∈ t}) (fun i => f i.1)
    (fun i => descendPieceOfProof (A := A) hU (f i.1) (hf i.1 i.2))
    (fun i => descendPieceIsoOfProof (A := A) hU (f i.1) (hf i.1 i.2))

/-- **`corrHypGlueDataOfCover`をさらに`Algebra.Etale`という1個の仮定
だけから自動的に得る**: `Γ(X,U)`が`A⊗[ℚ]ℝ`上étaleであれば
(`exists_finite_standardEtaleCover`で)有限standard-étale被覆
`ι,t,f`が自動的に存在するので、それを`corrHypGlueDataOfCover`に
渡すだけでよい。`ι,t,f`自体は`Exists.choose`で取り出した不透明な値
のまま外へは出さない(`corrHypGlueDataOfEtale_cover`が使う場面でも
同じ式`exists_finite_standardEtaleCover (A ⊗[ℚ] ℝ) Γ(X, U)`を再度
書くことで、定義的に同じ`.choose`結果を指すようにしている)。

★**sorry 無し**。 -/
noncomputable def corrHypGlueDataOfEtale {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    [Algebra.Etale (A ⊗[ℚ] ℝ) Γ(X, U)] : Scheme.GlueData :=
  letI h := exists_finite_standardEtaleCover (A ⊗[ℚ] ℝ) Γ(X, U)
  corrHypGlueDataOfCover (A := A) hU h.choose_spec.choose h.choose_spec.choose_spec.choose
    h.choose_spec.choose_spec.choose_spec.2

/-- `corrHypGlueDataOfEtale`が使う族`(ι,t,f)`は実際に`U`を覆う——
`exists_finite_standardEtaleCover`の被覆条件(環レベル、
`PrimeSpectrum.basicOpen`)を`exists_scheme_basicOpen_cover_of_ring`で
スキームレベル(`X.basicOpen`)へ変換するだけ。ロードマップ項目(b)
(`corrHypGlueData.glued ≅ U`)で必要になる被覆条件そのもの。 -/
theorem corrHypGlueDataOfEtale_cover {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    [Algebra.Etale (A ⊗[ℚ] ℝ) Γ(X, U)] :
    letI h := exists_finite_standardEtaleCover (A ⊗[ℚ] ℝ) Γ(X, U)
    (⨆ i ∈ h.choose_spec.choose, X.basicOpen (h.choose_spec.choose_spec.choose i)) = U := by
  exact exists_scheme_basicOpen_cover_of_ring hU _ _
    (exists_finite_standardEtaleCover (A ⊗[ℚ] ℝ) Γ(X, U)).choose_spec.choose_spec.choose_spec.1

/-- **`corrHypGlueDataOfEtale`の`.glued`が実際に`U`に同型である**——
`corrHypGlueData_glued_iso`(項目(b)、完成済み)を、`corrHypGlueDataOfEtale`
自身の内部構成(`corrHypGlueDataOfCover`経由で`corrHypGlueData`を
`J:={i//i∈t}`として呼ぶ)へそのまま適用するだけ。`corrHypGlueDataOfEtale_
cover`が与える`Finset`添字の被覆条件(`⨆i∈t,...`)を`iSup_subtype'`
(mathlib、`⨆i∈t,f i = ⨆i:{x//x∈t},f i.1`)で`corrHypGlueData_glued_iso`
が要求する部分型添字の形へ変換する。

これで「étaleな環拡大から作ったGlueDataの貼り合わせ結果が実際に元の
`U`を再構成する」ことが、抽象的な`corrHypGlueDataOfEtale`のレベルで
確立された——項目(b)で示した`IsColimit`比較を、この特殊化された文脈
でもう一度示す必要は無い(まさに一般化しておいた価値がここで出る)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypGlueDataOfEtale_glued_iso {A : Type} [CommRing A] [Algebra ℚ A]
    {X : Scheme} {U : X.Opens} (hU : IsAffineOpen U) [Algebra (A ⊗[ℚ] ℝ) Γ(X, U)]
    [Algebra.Etale (A ⊗[ℚ] ℝ) Γ(X, U)] :
    (corrHypGlueDataOfEtale (A := A) hU).glued ≅ (U : Scheme) := by
  unfold corrHypGlueDataOfEtale corrHypGlueDataOfCover
  exact corrHypGlueData_glued_iso _ _ _ (iSup_subtype'.symm.trans (corrHypGlueDataOfEtale_cover (A := A) hU))

/-! ### ロードマップ項目(c)への配線: `Ext X`へ有限étaleに写る実際の`C`

`piece_algebraEtale_tensor`(既存、`Lemma 4.1`「1アフィン片の降下」の
スキーム→環の橋渡しの完成形として以前に構築済み)の結論は、
`corrHypGlueDataOfEtale`が要求する`[Algebra.Etale (A⊗[ℚ]ℝ) Γ(X,U)]`
の形にちょうど一致する——`A := Γ(X.left,U)`(`X:Over BaseK`のアフィン
開`U`由来)、`corrHypGlueDataOfEtale`側の`X`には`C`(`Ext X`へ有限étale
に写る任意のスキーム)、`U`には`α ⁻¹ᵁ (piece)`を代入するだけでよい。
`hU`(アフィン性)は`piece_preimage_isAffineOpen`(既存)がちょうど
与える。これで「étaleな環拡大→GlueData」というCorrHyp非依存の部品
(`corrHypGlueDataOfEtale`)が、`Ext`/`ExtF`という実際のCorrHyp構成
データへ初めて直接接続された。 -/

/-- **`corrHypGlueDataOfEtale`を`Ext X`へ有限étaleに写る実際の`C`へ
接続する**: `piece_algebraEtale_tensor`の結論をそのまま
`corrHypGlueDataOfEtale`の`Algebra.Etale`引数へ流し込むだけの配線。
`X:Over BaseK`のアフィン開`U`由来の片`α⁻¹(piece)`(`C`側)に対応する
GlueDataを直接与える——`Lemma 4.1`の`c.α : c.C ⟶ D.Ext X`を
`X:=X.1`・`C:=c.C.1.left`・`α:=c.α.1.left`として代入すれば
(`corrHypInstance4`の下で`[IsFinite α][Etale α]`は`c.α`の証明成分
からそのまま得られる)、`Corr`の実データへの適用になる——ただし
そこへの実際の代入(`corrHypInstance4`・`QcqsSpace`・`Corr`の射影)は
まだ行っていない。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrPieceGlueData (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] : Scheme.GlueData :=
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  corrHypGlueDataOfEtale (A := Γ(X.left, U)) (piece_preimage_isAffineOpen X U hU C α)

/-- `corrPieceGlueData`が使う族は、実際に`α⁻¹(piece)`(`C`側で
`X.left`のアフィン開`U`に対応する片)を覆う——`corrHypGlueDataOfEtale_cover`
をそのまま呼ぶだけ。 -/
theorem corrPieceGlueData_cover (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    letI := pieceAlgebra X U hU
    letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
        (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
        (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
    haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
      piece_algebraEtale_tensor X U hU C α
    letI h := exists_finite_standardEtaleCover (Γ(X.left, U) ⊗[ℚ] ℝ)
      Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))
    (⨆ i ∈ h.choose_spec.choose,
      C.basicOpen (h.choose_spec.choose_spec.choose i)) = α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U) := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  exact corrHypGlueDataOfEtale_cover (A := Γ(X.left, U)) (piece_preimage_isAffineOpen X U hU C α)

/-- **`corrPieceGlueData ... |>.glued`が`α⁻¹(piece)`に同型である**——
`corrHypGlueDataOfEtale_glued_iso`(既に一般に確立済み)を、
`corrPieceGlueData = corrHypGlueDataOfEtale (A:=Γ(X.left,U)) ...`という
定義の一致(`show`で明示)からそのまま特殊化するだけ。以前の
「(b)IsColimit同士の比較、genuinely new」という見積もりは、項目(b)を
一般形(`corrHypGlueData_glued_iso`)で片付けておいたことで、この
特殊化ではもう新しい証明を要さない。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrPieceGlueData_glued_iso (X : Over BaseK) (U : X.left.Opens) (hU : IsAffineOpen U)
    (C : Scheme) (α : C ⟶ (ExtF.obj X).left) [IsFinite α] [Etale α] :
    (corrPieceGlueData X U hU C α).glued ≅ (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U) : Scheme) := by
  letI := pieceAlgebra X U hU
  letI : Algebra (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    ((Scheme.Hom.appLE α (pullback.fst X.hom toBaseK ⁻¹ᵁ U)
      (α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
      (pieceRingEquiv X U hU).symm.toRingHom).toAlgebra
  haveI : Algebra.Etale (Γ(X.left, U) ⊗[ℚ] ℝ) Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U)) :=
    piece_algebraEtale_tensor X U hU C α
  show (corrHypGlueDataOfEtale (A := Γ(X.left, U)) (piece_preimage_isAffineOpen X U hU C α)).glued ≅ _
  exact corrHypGlueDataOfEtale_glued_iso (piece_preimage_isAffineOpen X U hU C α)

/- ★★次の一手(未着手、2026-09-04に理解を訂正): 「`c.C.1.left`の外側の
貼り合わせ」ではない——`α⁻¹(piece)`は`C`自身の中に既に開集合として
存在する(トートロジー、貼り合わせ不要)。真に必要なのは、
`descendPiece`(`corrHypGlueData`の`Z i`)がまだ`ℝ`レベルのまま
(`piece_descends_iso`が与える`Spec(P₀.Ring)`という**`R`レベル**の
候補片が`descendPiece`内部で使い捨てられ、外へ取り出されていない)
ことを踏まえ、この`R`レベルの候補片自体を`R`レベル(`FgSubalgebra
ℚ ℝ`の圏)で貼り合わせる新しい層——異なる添字・異なるアフィン片
ごとに異なりうる`R_i`を共通精密化`R_i⊔R_j`へ持ち上げて比較する新しい
議論が要り、`ℝ`レベルで構築した`transitionElem`/`gdT`/`cocycle`一式
に相当する困難をもう一段繰り返すことを意味する(新しい独立した規模
の数学的内容、既存部品の配線だけでは済まない)。

さらに`α・β`脚と整合性の等式`h▸extCorr D c'=c`の構成——ただし
`h:ZK=D.Ext Z`という文字通りの命題的等号は、`Corr`定義の`Nonempty C`
欠落および`QcqsSpace`が同型類の商でないことと組み合わさると証明不可能
(あるいは偽)になりうるという構造的懸念が判明している(`corrhyp-
goal.md`2026-09-04の該当エントリに詳細記録、拙速な`Corr`修正は既に
完成済みの§1 5/5を後退させかねないため見送り中)。 -/

end ABC3.Found.CorrHyp
