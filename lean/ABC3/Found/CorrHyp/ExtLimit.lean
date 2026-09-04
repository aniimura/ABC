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

/- ★★次の一手(未着手): `extCone X` が極限であることを示す
(`isLimit_extCone : IsLimit (extCone X)`)——`s : Cone (extDiagram X)` から
`pullback.lift (s.π.app R0 ≫ pullback.fst ...) ((isLimit_specKCone ℚ ℝ).lift
(auxCone3 X s)) _` を構成する。互換性条件(`pullback.lift` の第三引数)は
`pullback.condition`(`fst ≫ X.hom = snd ≫ (D R0).hom`)・`extCone_fst_const`
非依存の `hfac`/`htoBaseK` の `calc` チェーンで**完全に閉じることを確認済み**
(この commit の直前まで動作した)。残るのは `fac`(`lift ≫ π.app R = π_R`、
`pullback.hom_ext` で fst/snd に割ってから `extCone_fst_const`/`(isLimit_
specKCone).fac` で閉じる見込み)と `uniq`——`pullback.map`/`pullback.lift` が
入れ子になった式で `Category.assoc`/`unfold Limits.pullback.map` の適用順を
間違えると `(specKCone ℚ ℝ).pt` と `specK` の型不一致に化ける、という配管が
最後まで残った。閉じれば `FieldLimit.lean` の
`standardEtalePairRingBaseChange`(環側の base change)と合わせて
`Lemma 4.1` の構成的降下に必要なスキーム側・環側の道具が完全に揃う。 -/

end ABC3.Found.CorrHyp
