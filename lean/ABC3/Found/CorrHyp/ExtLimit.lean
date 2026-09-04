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

end ABC3.Found.CorrHyp
