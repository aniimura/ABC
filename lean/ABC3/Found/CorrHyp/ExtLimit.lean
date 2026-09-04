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

/- ★★次の一手(未着手): `X : Over BaseK`(`Lemma 4.1` の `X`)を固定し、
`Over.pullback X.hom : Over BaseK ⥤ Over X.left`(`X.hom` に沿った base
change の関手、`Over.mapPullbackAdj` の右随伴なので `PreservesLimit` が
`infer_instance` で自動的に付くことを確認済み)を `specKConeOver`/
`isLimit_specKConeOver` に適用すれば、`X.left ×_{BaseK} Spec K` の極限表示
(`X.left ×_{BaseK} Spec K = lim (X.left ×_{BaseK} Spec R)`)が得られる。
唯一の詰まりは `Limits.pullback X.hom toBaseK`(`Ext X` の定義に使う順序)と
`Limits.pullback toBaseK X.hom`(`Over.pullback X.hom` が自然に与える順序)の
**引数の順序が入れ替わっている**こと——`pullbackSymmetry`(mathlib既存の
`pullback f g ≅ pullback g f`)で埋める必要がある。埋まれば
`FieldLimit.lean` の `standardEtalePairRingBaseChange` と組み合わせて
`Lemma 4.1` の構成的降下の全体(スキーム側の極限表示+環側の base change)が
揃う。 -/

end ABC3.Found.CorrHyp
