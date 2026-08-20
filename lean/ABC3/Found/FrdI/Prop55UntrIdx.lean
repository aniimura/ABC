/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop33Classes
import ABC3.Found.FrdI.Prop55Untr

/-!
# [FrdI] Proposition 5.5, (ii) —— 添字圏の同一視(un-tr の側を閉じる)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> erality that C is of isotropic type. Then it follows immediately from the definition of

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C

## ★★何が残っていたか

`Prop55Untr.lean` は「`Hom^pf(A,B)` を unit-equivalence で割ったもの」と
「各段の `Hom^un-tr` の**帰納極限**」が一致することを示した。
★しかし `(𝒞^un-tr)^pf` の射の集合は、**`𝒞^un-tr` 自身の添字圏**の上の帰納極限である。
★★**添字圏が違う** —— `𝒞^istr` の Frobenius 型射の対 と `𝒞^un-tr` のそれ。

★本ファイルはその 2 つの添字圏を突き合わせる。

## ★★★鍵は 1 本の一般補題

**`final_of_essSurj_of_thin`** —— filtered な源から**細い**的への
本質的全射な関手は**終尾(final)**である。

★証明は 2 行:
* comma 圏が空でない —— 本質的全射性
* comma 圏が連結 —— 源の filtered 性で共通上界を取れば、
  ★**的が細いので 2 本の射が自動的に一致する**

★★**細さ(`idx_hom_ext`)は `𝒞` が totally epimorphic であることから来る** ——
すなわち **pre-Frobenioid の定義そのもの**である。

## ★仮定について

`idxToUnTr` が本質的全射であることは、`istrToUnTr` の**充満性**(在庫)と
**Frobenius 型・`degFr` の両向きの保存**
(`unTr_isFrobeniusType_iff'` / `unTr_degFr_iff`、在庫)から出る。
★★**新しい仮定は 1 つも要らない。**

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `final_of_essSurj_of_thin` | ★filtered から細い圏への本質的全射な関手は終尾 |
| `biFrToUnTr` / `idxToUnTr` | 添字圏のあいだの関手 |
| `idxToUnTr_essSurj` | 本質的全射(`istrToUnTr` の充満性から) |
| `idxToUnTr_final` | ★★**終尾** |
| `istrToUnTr_frobTransport` | 遷移写像が対応すること(`frobTransport` の一意性 1 本) |
| `stageNatIso` | 図式の同一視 |
| `prop_5_5_ii_untr` | ★★★★**`Hom^pf(A,B)/∼ ≃ Hom^pf_{𝒞^un-tr}(A,B)`** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe uJ vJ uK vK

/-! ## ★1. 一般補題 —— filtered から細い圏への本質的全射な関手は終尾 -/

/-- ★★★**filtered な源から細い的への本質的全射な関手は終尾(final)**。

★comma 圏が空でないのは本質的全射性から。連結性は、源の filtered 性で
共通上界を取れば、★**的が細いので 2 本の射が自動的に一致する**ことから。 -/
theorem final_of_essSurj_of_thin {J : Type uJ} [Category.{vJ} J] [IsFiltered J]
    {K : Type uK} [Category.{vK} K] (G : J ⥤ K)
    (hthin : ∀ {X Y : K} (f g : X ⟶ Y), f = g) [G.EssSurj] : G.Final := by
  refine ⟨fun d => ?_⟩
  haveI : Nonempty (StructuredArrow d G) :=
    ⟨StructuredArrow.mk (G.objObjPreimageIso d).inv⟩
  refine zigzag_isConnected (fun X Y => ?_)
  let c := IsFiltered.max X.right Y.right
  let uX : X.right ⟶ c := IsFiltered.leftToMax _ _
  let uY : Y.right ⟶ c := IsFiltered.rightToMax _ _
  have hXY : X.hom ≫ G.map uX = Y.hom ≫ G.map uY := hthin _ _
  refine Zigzag.trans (Zigzag.of_hom (StructuredArrow.homMk uX rfl :
      X ⟶ StructuredArrow.mk (X.hom ≫ G.map uX))) ?_
  refine Zigzag.of_inv (StructuredArrow.homMk uY ?_ :
      Y ⟶ StructuredArrow.mk (X.hom ≫ G.map uX))
  exact hXY.symm

/-! ## ★2. 添字圏のあいだの関手 -/

section UntrIdx

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)
  (F₁ : FrobenioidCore (istrPre P Fc)) (F₂ : FrobenioidCore (unTrPre P Fc))

/-- ★**`𝒞^{bi-Fr}` のあいだの関手** —— 各成分を `istrToUnTr` で送る。

★射のクラスが保たれるのは `unTr_isFrobeniusType_iff'` と `degFr` の一致から。 -/
def biFrToUnTr : BiFr (istrPre P Fc) F₁ ⥤ BiFr (unTrPre P Fc) F₂ where
  obj Z := ⟨((istrToUnTr P).obj Z.obj.1, (istrToUnTr P).obj Z.obj.2)⟩
  map {Z W} u :=
    ⟨((istrToUnTr P).map u.hom.1, (istrToUnTr P).map u.hom.2),
      (unTr_isFrobeniusType_iff' P Fc _).mpr u.property.1,
      (unTr_isFrobeniusType_iff' P Fc _).mpr u.property.2.1,
      u.property.2.2⟩
  map_id Z := by
    apply InducedWideCategory.Hom.ext
    exact Prod.ext ((istrToUnTr P).map_id _) ((istrToUnTr P).map_id _)
  map_comp u u' := by
    apply InducedWideCategory.Hom.ext
    exact Prod.ext ((istrToUnTr P).map_comp _ _) ((istrToUnTr P).map_comp _ _)

/-- ★底の対象は動かない。 -/
theorem biFrToUnTr_obj_biFrObj (A B : Istr P) :
    (biFrToUnTr P Fc F₁ F₂).obj (biFrObj (istrPre P Fc) F₁ A B)
      = biFrObj (unTrPre P Fc) F₂ A B := rfl

/-- ★★**添字圏のあいだの関手** —— コスライスへ持ち上げただけ。 -/
def idxToUnTr (A B : Istr P) :
    IdxPf (istrPre P Fc) F₁ A B ⥤ IdxPf (unTrPre P Fc) F₂ A B :=
  Under.post (biFrToUnTr P Fc F₁ F₂)

/-- ★★**本質的全射** —— `istrToUnTr` の充満性で対を持ち上げ、
Frobenius 型と `degFr` を逆向きに引き戻す。 -/
instance idxToUnTr_essSurj (A B : Istr P) :
    (idxToUnTr P Fc F₁ F₂ A B).EssSurj where
  mem_essImage d := by
    obtain ⟨a, ha⟩ := (istrToUnTr P).map_surjective d.hom.hom.1
    obtain ⟨b, hb⟩ := (istrToUnTr P).map_surjective d.hom.hom.2
    have hFa : IsFrobeniusType (istrPre P Fc) a := by
      refine (unTr_isFrobeniusType_iff' P Fc a.hom).mp ?_
      have h1 : IsFrobeniusType (unTrPre P Fc) d.hom.hom.1 := d.hom.property.1
      rw [← ha] at h1
      exact h1
    have hFb : IsFrobeniusType (istrPre P Fc) b := by
      refine (unTr_isFrobeniusType_iff' P Fc b.hom).mp ?_
      have h1 : IsFrobeniusType (unTrPre P Fc) d.hom.hom.2 := d.hom.property.2.1
      rw [← hb] at h1
      exact h1
    have hdeg : (istrPre P Fc).degFr a = (istrPre P Fc).degFr b := by
      have h : (unTrPre P Fc).degFr d.hom.hom.1 = (unTrPre P Fc).degFr d.hom.hom.2 :=
        d.hom.property.2.2
      rw [← ha, ← hb] at h
      exact h
    let W : BiFr (istrPre P Fc) F₁ := ⟨d.right.obj⟩
    let u : biFrObj (istrPre P Fc) F₁ A B ⟶ W := ⟨(a, b), hFa, hFb, hdeg⟩
    refine ⟨Under.mk u, ⟨Under.isoMk (Iso.refl _) ?_⟩⟩
    show (biFrToUnTr P Fc F₁ F₂).map u ≫ 𝟙 _ = d.hom
    rw [Category.comp_id]
    apply InducedWideCategory.Hom.ext
    exact Prod.ext ha hb

/-- ★★★**終尾(final)** —— 一般補題を当てるだけ。
細さは `idx_hom_ext`(`𝒞^un-tr` が totally epimorphic であることの帰結)。 -/
instance idxToUnTr_final (A B : Istr P) :
    (idxToUnTr P Fc F₁ F₂ A B).Final :=
  final_of_essSurj_of_thin _ (fun f g => idx_hom_ext f g)

/-! ## ★3. 図式の同一視 -/

/-- ★★**遷移写像が対応する** —— `frobTransport` の**一意性** 1 本で出る。

★`istrToUnTr` を `Proposition 1.10, (i)` の四角形に当てれば、
その像が `𝒞^un-tr` 側の四角形になる。 -/
theorem istrToUnTr_frobTransport {A' B' A'' B'' : Istr P}
    (a : A' ⟶ A'') (ha : IsFrobeniusType (istrPre P Fc) a)
    (b : B' ⟶ B'') (hb : IsFrobeniusType (istrPre P Fc) b)
    (hd : (istrPre P Fc).degFr a = (istrPre P Fc).degFr b) (φ : A' ⟶ B')
    (ha' : IsFrobeniusType (unTrPre P Fc) ((istrToUnTr P).map a))
    (hb' : IsFrobeniusType (unTrPre P Fc) ((istrToUnTr P).map b))
    (hd' : (unTrPre P Fc).degFr ((istrToUnTr P).map a)
      = (unTrPre P Fc).degFr ((istrToUnTr P).map b)) :
    (istrToUnTr P).map (frobTransport (F := F₁) a ha b hb hd φ)
      = frobTransport (F := F₂) _ ha' _ hb' hd' ((istrToUnTr P).map φ) :=
  (frobTransport_eq (F := F₂) _ ha' _ hb' hd' _ _ (by
    rw [← (istrToUnTr P).map_comp, ← (istrToUnTr P).map_comp,
      frobTransport_spec (F := F₁) a ha b hb hd φ])).symm

/-- ★各段の同一視 —— 商と `Hom_{𝒞^un-tr}` は同じもの。 -/
noncomputable def stageIso (A B : Istr P) (Z : IdxPf (istrPre P Fc) F₁ A B) :
    (unTrQuotPf (istrPre P Fc) F₁ A B).functor.obj Z
      ≃ ((idxToUnTr P Fc F₁ F₂ A B) ⋙ homFunctorPf (unTrPre P Fc) F₂ A B).obj Z where
  toFun := Quotient.lift (fun φ => ULift.up ((istrToUnTr P).map φ.down))
    (fun _ _ h => congrArg ULift.up (Quotient.sound h))
  invFun y := Quotient.liftOn y.down
    (fun α => Quotient.mk ((unTrQuotPf (istrPre P Fc) F₁ A B).setoid Z)
      (ULift.up (ObjectProperty.homMk α)))
    (fun _ _ h => Quotient.sound h)
  left_inv z := by
    refine Quotient.inductionOn z (fun φ => ?_)
    obtain ⟨φ⟩ := φ
    rfl
  right_inv y := by
    obtain ⟨y⟩ := y
    refine Quotient.inductionOn y (fun α => ?_)
    rfl

/-- ★★**図式の同一視** —— 自然性は `istrToUnTr_frobTransport` そのもの。 -/
noncomputable def stageNatIso (A B : Istr P) :
    (unTrQuotPf (istrPre P Fc) F₁ A B).functor
      ≅ (idxToUnTr P Fc F₁ F₂ A B) ⋙ homFunctorPf (unTrPre P Fc) F₂ A B :=
  NatIso.ofComponents (fun Z => (stageIso P Fc F₁ F₂ A B Z).toIso) (by
    intro Z W u
    ext z
    refine Quotient.inductionOn z (fun φ => ?_)
    exact congrArg ULift.up (istrToUnTr_frobTransport P Fc F₁ F₂
      u.right.hom.1 u.right.property.1 u.right.hom.2 u.right.property.2.1
      u.right.property.2.2 φ.down
      ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.1
      ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.2.1
      ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.2.2))

/-! ## ★4. 主定理 -/

/-- ★★★**帰納極限の一致** —— 図式の同一視と終尾性を継ぐ。 -/
noncomputable def homPfUnTrColimIso (A B : Istr P) :
    HomColim ((unTrQuotPf (istrPre P Fc) F₁ A B).functor)
      ≅ HomPf (unTrPre P Fc) F₂ A B :=
  HasColimit.isoOfNatIso (stageNatIso P Fc F₁ F₂ A B) ≪≫
    Functor.Final.colimitIso (idxToUnTr P Fc F₁ F₂ A B) (homFunctorPf (unTrPre P Fc) F₂ A B)

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)**(un-tr の側)——
`(𝒞^pf)^un-tr` の射の集合と `(𝒞^un-tr)^pf` の射の集合の**自然な全単射**。

★左辺は「`Hom^pf` を unit-equivalence で割ったもの」、
右辺は「`𝒞^un-tr` の `Hom^pf`」である。

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C -/
noncomputable def prop_5_5_ii_untr (A B : Istr P) :
    Quotient (HomColim.quotKer (unTrQuotPf (istrPre P Fc) F₁ A B))
      ≃ HomPf (unTrPre P Fc) F₂ A B :=
  (HomColim.quotEquiv _).trans (homPfUnTrColimIso P Fc F₁ F₂ A B).toEquiv

end UntrIdx

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (ii)`(un-tr の側、添字圏の同一視まで込み)。 -/
def prop_5_5_ii_untr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — (𝒞^pf)^un-tr と (𝒞^un-tr)^pf の射の集合の自然な全単射(添字圏の同一視)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
