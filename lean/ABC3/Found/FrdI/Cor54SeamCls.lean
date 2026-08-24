/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54Birat
import ABC3.Found.FrdI.Thm52Model
import ABC3.Found.FrdI.Cor54Seam

/-!
# [FrdI] Corollary 5.4 の縦の矢印の継ぎ目 —— 基準切断を揃えれば閉じる

原文 (FrdI p.104):
> with the 1-commutative diagram of Proposition 5.3 [involving perfections, unit-

★★`Cor54SeamUnTr.lean` の測定で、継ぎ目に残る 1 条が
「**基準切断 `S₁` / `S₂` を揃えれば消える**」ものであることが分かった。
本ファイルはその**揃えた場合の計算**を閉じる。

## ★本ファイルが閉じること

| 定理 | 中身 |
|---|---|
| `FPPath.mapAlong` | ★`F-𝒫-path` を `Ψ` で送る(`Ψ` が pre-step と `𝒫` を保てばよい) |
| `FPPath.gpMap_cls_mapAlong` | ★★★★★**運んだ path の類は、もとの類を `η` で送ったものと一致** |
| `pathSeamIso` | ★★★★★★`PathCat` の層の 1-可換図式(**条なし**) |
| `conjIsoOfSquare` | ★★★★四角形を擬逆で共役する(一般補題) |
| `cSeamIso` | ★★★★★★**`𝒞` の層の 1-可換図式(条なし)** —— 継ぎ目そのもの |

## ★★中身(記録)

`Theorem 5.2, (iv)` の類は

```
FPPath.cls p = -spanCls p.toObj _ p.toRef
             = -Φ.gpMapOn (inv (Base toObj)) (toGp (Div toRef) - toGp (Div toObj))
```

なので、`Div` の対応(`Corollary 4.11, (iv)` の `hdivc`)と
底の 1-可換図式 `sq` の自然性だけで運べる。★要になるのは

```
ΨB.map (inv (Base₁ toObj)) ≫ sq.app vertex = sq.app X ≫ inv (Base₂ (Ψ toObj))
```

という 1 本で、これは `sq` の `toObj` での自然性の**両側を逆射でずらしたもの**である。

★★★これで「基準を揃えれば `γ(d) = 0`」が**計算として**確かめられ、
そのまま `pathSeamIso`(`PathCat` の層)→ `cSeamIso`(`𝒞` の層)と組み上がった。

★★★★★**結論**: `Corollary 5.4` の縦の矢印の継ぎ目は、
`modelType_equiv`(`Classical.choice` が切断を選ぶ)ではなく
**`thm_5_2_iv` を使って切断を明示し、`S₂ := Ψ_*(S₁)`(在庫 `BaseSection.map`、
`Corollary 5.7, (i)`)と揃えれば、仮定なしで 1-可換**である。
`Cor54SeamUnTr.lean` の `hex` は、この形では**そもそも要らない**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 uA vA uB vB uA2 vA2 uB2 vB2

section ClsMap

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★**`F-𝒫-path` を `Ψ` で送る**。

★`Ψ` が pre-step を pre-step へ送り(`hPS`)、`𝒫` の対象を `𝒫` の対象へ送れば(`hsec`)、
path はそのまま送れる。★`hsec` は `S₂ := S₁.map Ψ`(`Corollary 5.7, (i)`)なら**定義から**成り立つ。 -/
def FPPath.mapAlong (Ψ : C₁ ⥤ C₂) {S₁ : BaseSection P₁} {S₂ : BaseSection P₂}
    (hsec : ∀ {A : C₁}, S₁.objP A → S₂.objP (Ψ.obj A))
    (hPS : ∀ {A B : C₁} (f : A ⟶ B), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    {X : C₁} (p : FPPath S₁ X) : FPPath S₂ (Ψ.obj X) where
  ref := Ψ.obj p.ref
  ref_mem := hsec p.ref_mem
  vertex := Ψ.obj p.vertex
  toRef := Ψ.map p.toRef
  toObj := Ψ.map p.toObj
  toRef_preStep := hPS _ p.toRef_preStep
  toObj_preStep := hPS _ p.toObj_preStep

/-- ★★★★**要の 1 本** —— `sq` の自然性を逆射でずらした形。 -/
theorem seam_base_key (Ψ : C₁ ⥤ C₂) (ΨB : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    {V X : C₁} (φ : V ⟶ X) (h₁ : IsIso (P₁.Base φ)) (h₂ : IsIso (P₂.Base (Ψ.map φ))) :
    ΨB.map (@inv _ _ _ _ (P₁.Base φ) h₁) ≫ sq.hom.app V
      = sq.hom.app X ≫ @inv _ _ _ _ (P₂.Base (Ψ.map φ)) h₂ := by
  haveI i₁ := h₁
  haveI i₂ := h₂
  have hnat : ΨB.map (P₁.Base φ) ≫ sq.hom.app X
      = sq.hom.app V ≫ P₂.Base (Ψ.map φ) := sq.hom.naturality φ
  have e1 : ΨB.map (P₁.Base φ) ≫ (ΨB.map (@inv _ _ _ _ (P₁.Base φ) h₁) ≫ sq.hom.app V)
      = sq.hom.app V := by
    rw [← Category.assoc, ← ΨB.map_comp, IsIso.hom_inv_id, ΨB.map_id, Category.id_comp]
  have e2 : ΨB.map (P₁.Base φ) ≫ (sq.hom.app X ≫ @inv _ _ _ _ (P₂.Base (Ψ.map φ)) h₂)
      = sq.hom.app V := by
    have s1 : ΨB.map (P₁.Base φ) ≫ (sq.hom.app X ≫ @inv _ _ _ _ (P₂.Base (Ψ.map φ)) h₂)
        = (ΨB.map (P₁.Base φ) ≫ sq.hom.app X) ≫ @inv _ _ _ _ (P₂.Base (Ψ.map φ)) h₂ :=
      (Category.assoc _ _ _).symm
    have s2 : (sq.hom.app V ≫ P₂.Base (Ψ.map φ)) ≫ @inv _ _ _ _ (P₂.Base (Ψ.map φ)) h₂
        = sq.hom.app V ≫ (P₂.Base (Ψ.map φ) ≫ @inv _ _ _ _ (P₂.Base (Ψ.map φ)) h₂) :=
      Category.assoc _ _ _
    have s3 : P₂.Base (Ψ.map φ) ≫ @inv _ _ _ _ (P₂.Base (Ψ.map φ)) h₂ = 𝟙 _ :=
      @IsIso.hom_inv_id _ _ _ _ (P₂.Base (Ψ.map φ)) h₂
    rw [s1, hnat]
    exact s2.trans ((congrArg (fun t => sq.hom.app V ≫ t) s3).trans (Category.comp_id _))
  exact (cancel_epi (ΨB.map (P₁.Base φ))).mp (e1.trans e2.symm)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 5.4 の継ぎ目の本体** ——
**運んだ path の類は、もとの類を `η` で送ったものと一致する。**

★入力は `Corollary 4.11, (iii)(iv)` が与えるもの(`ΨB` / `η` / `sq` / `hdivc`)と
`Ψ` が pre-step・`𝒫` を保つことだけである。 -/
theorem FPPath.gpMap_cls_mapAlong (Ψ : C₁ ⥤ C₂) (ΨB : D₁ ⥤ D₂)
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    (hdivc : ∀ {A B : C₁} (φ : A ⟶ B),
      phiIsoApp ΨB η ((P₁.toElem.obj A).base) (P₁.Div φ)
        = Φ₂.map (sq.hom.app A) (P₂.Div (Ψ.map φ)))
    {S₁ : BaseSection P₁} {S₂ : BaseSection P₂}
    (hsec : ∀ {A : C₁}, S₁.objP A → S₂.objP (Ψ.obj A))
    (hPS : ∀ {A B : C₁} (f : A ⟶ B), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    {X : C₁} (p : FPPath S₁ X) :
    gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj X).base)) p.cls
      = gpMap _ (Φ₂.map (sq.hom.app X)) (p.mapAlong Ψ hsec hPS).cls := by
  haveI i₁ : IsIso (P₁.Base p.toObj) := p.toObj_preStep.2
  haveI i₂ : IsIso (P₂.Base (Ψ.map p.toObj)) := (hPS _ p.toObj_preStep).2
  have hkey := seam_base_key Ψ ΨB sq p.toObj i₁ i₂
  -- ★`Div` の対応を `Gp` の層に持ち上げる
  have hdivgp : ∀ {A B : C₁} (φ : A ⟶ B),
      gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj A).base)) (toGp _ (P₁.Div φ))
        = Φ₂.gpMapOn (sq.hom.app A) (toGp _ (P₂.Div (Ψ.map φ))) := by
    intro A B φ
    rw [gpMap_toGp, hdivc]
    exact (gpMap_toGp _ (Φ₂.map (sq.hom.app A)) _).symm
  have hsub : gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj p.vertex).base))
        (toGp _ (P₁.Div p.toRef) - toGp _ (P₁.Div p.toObj))
      = Φ₂.gpMapOn (sq.hom.app p.vertex)
        (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj))) := by
    have hr : Φ₂.gpMapOn (sq.hom.app p.vertex)
          (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj)))
        = Φ₂.gpMapOn (sq.hom.app p.vertex) (toGp _ (P₂.Div (Ψ.map p.toRef)))
          - Φ₂.gpMapOn (sq.hom.app p.vertex) (toGp _ (P₂.Div (Ψ.map p.toObj))) :=
      map_sub _ _ _
    rw [map_sub, hr, hdivgp, hdivgp]
    rfl
  -- ★符号を外した形で計算する(型の見た目を揃えるため `gpMap` の第 1 引数を明示する)
  have hspan : gpMap (Φ₁.val (P₁.toElem.obj X).base)
        (phiIsoApp ΨB η ((P₁.toElem.obj X).base)) (spanCls p.toObj i₁ p.toRef)
      = gpMap (Φ₂.val (P₂.toElem.obj (Ψ.obj X)).base) (Φ₂.map (sq.hom.app X))
        (spanCls (Ψ.map p.toObj) i₂ (Ψ.map p.toRef)) := by
    have e1 : gpMap (Φ₁.val (P₁.toElem.obj X).base)
          (phiIsoApp ΨB η ((P₁.toElem.obj X).base)) (spanCls p.toObj i₁ p.toRef)
        = Φ₂.gpMapOn (ΨB.map (@inv _ _ _ _ (P₁.Base p.toObj) i₁))
            (gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj p.vertex).base))
              (toGp _ (P₁.Div p.toRef) - toGp _ (P₁.Div p.toObj))) := by
      rw [spanCls_eq]
      exact gpMap_phiIsoApp_nat ΨB η (@inv _ _ _ _ (P₁.Base p.toObj) i₁) _
    rw [e1, hsub, spanCls_eq]
    show Φ₂.gpMapOn (ΨB.map (@inv _ _ _ _ (P₁.Base p.toObj) i₁))
        (Φ₂.gpMapOn (sq.hom.app p.vertex)
          (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj))))
      = Φ₂.gpMapOn (sq.hom.app X)
        (Φ₂.gpMapOn (@inv _ _ _ _ (P₂.Base (Ψ.map p.toObj)) i₂)
          (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj))))
    calc Φ₂.gpMapOn (ΨB.map (@inv _ _ _ _ (P₁.Base p.toObj) i₁))
          (Φ₂.gpMapOn (sq.hom.app p.vertex)
            (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj))))
        = Φ₂.gpMapOn (ΨB.map (@inv _ _ _ _ (P₁.Base p.toObj) i₁) ≫ sq.hom.app p.vertex)
            (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj))) :=
          (Φ₂.gpMapOn_comp _ _ _).symm
      _ = Φ₂.gpMapOn (sq.hom.app X ≫ (@inv _ _ _ _ (P₂.Base (Ψ.map p.toObj)) i₂))
            (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj))) := by
          rw [hkey]
          rfl
      _ = Φ₂.gpMapOn (sq.hom.app X)
            (Φ₂.gpMapOn (@inv _ _ _ _ (P₂.Base (Ψ.map p.toObj)) i₂)
              (toGp _ (P₂.Div (Ψ.map p.toRef)) - toGp _ (P₂.Div (Ψ.map p.toObj)))) :=
          Φ₂.gpMapOn_comp _ _ _
  show gpMap (Φ₁.val (P₁.toElem.obj X).base)
        (phiIsoApp ΨB η ((P₁.toElem.obj X).base)) (-(spanCls p.toObj i₁ p.toRef))
    = gpMap (Φ₂.val (P₂.toElem.obj (Ψ.obj X)).base) (Φ₂.map (sq.hom.app X))
        (-(spanCls (Ψ.map p.toObj) i₂ (Ψ.map p.toRef)))
  rw [map_neg, map_neg, hspan]
  rfl

/-! ## ★3. `PathCat` の層で継ぎ目を組む -/

section PathSeam

variable [IsConnected D₁] [IsConnected D₂]

/-- ★★**`Ψ` が誘導する `𝒞̃₁ ⥤ 𝒞̃₂`** —— 対象は path ごと送る。 -/
def pathMapAlong (Ψ : C₁ ⥤ C₂) {S₁ : BaseSection P₁} {S₂ : BaseSection P₂}
    (hsec : ∀ {A : C₁}, S₁.objP A → S₂.objP (Ψ.obj A))
    (hPS : ∀ {A B : C₁} (f : A ⟶ B), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f)) :
    PathCat S₁ ⥤ PathCat S₂ where
  obj X := ⟨Ψ.obj X.obj, X.path.mapAlong Ψ hsec hPS⟩
  map f := Ψ.map f
  map_id X := Ψ.map_id X.obj
  map_comp f g := Ψ.map_comp f g

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 5.4 の縦の矢印の継ぎ目**(`PathCat` の層、**条なし**)——

```
𝒞̃₁ ⥤ Model₁ --Ψ^model--> Model₂
 |                          ‖
 Ψ̃                          ‖
 v                          ‖
𝒞̃₂ ⥤ Model₂ ═══════════════╝
```

★★**基準切断を `S₂ := Ψ_*(S₁)` と揃えれば、この四角形は仮定なしで 1-可換**である
(`hsec` がその「揃える」の中身)。

★入力は `Corollary 4.11, (iii)(iv)` が与えるもの ——
底の 1-可換図式 `sq`、因子の単系の同型 `η`、`Div` の対応 `hdivc`、
Frobenius 次数の保存 `hdegc` —— と、
`Ψ` が pre-step・`𝒫` を保つこと(`hPS` / `hsec`)だけである。
★★★`u` の帳尻は `Div_B` の単射性(`hinj`)から**自動**であり、
対象の類の帳尻は `FPPath.gpMap_cls_mapAlong` により **`u = 0`** で足りる。 -/
noncomputable def pathSeamIso (Ψ : C₁ ⥤ C₂) (ΨB : D₁ ⥤ D₂)
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    (hdivc : ∀ {A B : C₁} (φ : A ⟶ B),
      phiIsoApp ΨB η ((P₁.toElem.obj A).base) (P₁.Div φ)
        = Φ₂.map (sq.hom.app A) (P₂.Div (Ψ.map φ)))
    (hdegc : ∀ {A B : C₁} (φ : A ⟶ B), P₁.degFr φ = P₂.degFr (Ψ.map φ))
    {G₁ : Frobenioid P₁} (R₁ : RatFnData P₁ G₁)
    (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ Z : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) Z)
    {S₁ : BaseSection P₁} (Fs₁ : ℕ+ →* SectionEnd S₁) (hFs₁ : IsFrobeniusSection S₁ Fs₁)
    {G₂ : Frobenioid P₂} (R₂ : RatFnData P₂ G₂)
    (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ Z : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) Z)
    {S₂ : BaseSection P₂} (Fs₂ : ℕ+ →* SectionEnd S₂) (hFs₂ : IsFrobeniusSection S₂ Fs₂)
    (hsec : ∀ {A : C₁}, S₁.objP A → S₂.objP (Ψ.obj A))
    (hPS : ∀ {A B : C₁} (f : A ⟶ B), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (Fo : ModelDataHomOver ΨB R₁.model R₂.model)
    (hphi : ∀ d : D₁, Fo.phiHom d = phiIsoApp ΨB η d)
    (hinj : ∀ d : D₂, Function.Injective (R₂.divB d)) :
    pathToModel R₁ hiso₁ hfn₁ S₁ Fs₁ hFs₁ ⋙ Fo.functor
      ≅ pathMapAlong Ψ hsec hPS ⋙ pathToModel R₂ hiso₂ hfn₂ S₂ Fs₂ hFs₂ :=
  ModelData.squareOfBaseOfInj hinj R₂.bneg R₂.bneg_add
    (fun X => sq.app X.obj)
    (fun {_ _} f => sq.hom.naturality f)
    (fun {X _} f => by
      show Fo.phiHom _ (P₁.Div f) = Φ₂.map (sq.hom.app X.obj) (P₂.Div (Ψ.map f))
      rw [hphi]
      exact hdivc f)
    (fun {_ _} f => hdegc f)
    (fun _ => 0)
    (fun X => by
      simp only [map_zero, add_zero]
      show gpMap _ (Fo.phiHom _) X.path.cls
        = Φ₂.gpMapOn (sq.hom.app X.obj) ((X.path.mapAlong Ψ hsec hPS).cls)
      rw [hphi]
      exact FPPath.gpMap_cls_mapAlong Ψ ΨB η sq hdivc hsec hPS X.path)

end PathSeam

end ClsMap

/-! ## ★4. 圏同値で `𝒞` の層へ降ろす -/

section Conj

variable {A : Type uA} [Category.{vA} A] {B : Type uB} [Category.{vB} B]
  {A' : Type uA2} [Category.{vA2} A'] {B' : Type uB2} [Category.{vB2} B']

/-- ★★★★**四角形を擬逆で共役する** ——
`K ⋙ F₂ = F₁ ⋙ H` かつ `F₁` `F₂` が圏同値なら `F₁⁻¹ ⋙ K ≅ H ⋙ F₂⁻¹`。

★`Theorem 5.2, (iv)` の圏同値は `𝒞̃ ⥤ 𝒞` の擬逆を経由するので、
`PathCat` の層で作った 1-可換図式を `𝒞` の層へ降ろすのにこれが要る。 -/
noncomputable def conjIsoOfSquare (F₁ : A ⥤ B) [F₁.IsEquivalence]
    (F₂ : A' ⥤ B') [F₂.IsEquivalence] (K : A ⥤ A') (H : B ⥤ B')
    (h : K ⋙ F₂ = F₁ ⋙ H) :
    F₁.asEquivalence.inverse ⋙ K ≅ H ⋙ F₂.asEquivalence.inverse :=
  Functor.isoWhiskerLeft F₁.asEquivalence.inverse
      ((Functor.rightUnitor K).symm
        ≪≫ Functor.isoWhiskerLeft K F₂.asEquivalence.unitIso
        ≪≫ (Functor.associator K F₂ F₂.asEquivalence.inverse).symm
        ≪≫ Functor.isoWhiskerRight (eqToIso h) F₂.asEquivalence.inverse
        ≪≫ Functor.associator F₁ H F₂.asEquivalence.inverse)
    ≪≫ (Functor.associator F₁.asEquivalence.inverse F₁
          (H ⋙ F₂.asEquivalence.inverse)).symm
    ≪≫ Functor.isoWhiskerRight F₁.asEquivalence.counitIso
        (H ⋙ F₂.asEquivalence.inverse)
    ≪≫ Functor.leftUnitor _

end Conj

section CSeam

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  [IsConnected D₁] [IsConnected D₂]

/-- ★★`Ψ̃` は忘却関手と**厳密に**可換(`rfl`)。 -/
theorem pathMapAlong_pathForget (Ψ : C₁ ⥤ C₂) {S₁ : BaseSection P₁} {S₂ : BaseSection P₂}
    (hsec : ∀ {A : C₁}, S₁.objP A → S₂.objP (Ψ.obj A))
    (hPS : ∀ {A B : C₁} (f : A ⟶ B), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f)) :
    pathMapAlong Ψ hsec hPS ⋙ pathForget S₂ = pathForget S₁ ⋙ Ψ := rfl

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 5.4 の縦の矢印の継ぎ目**(`𝒞` の層、**条なし**)——

```
𝒞₁ ≌ Model₁ --Ψ^model--> Model₂
 |                         ‖
 Ψ                         ‖
 v                         ‖
𝒞₂ ≌ Model₂ ═══════════════╝
```

★★`Theorem 5.2, (iv)` の圏同値(`thm_5_2_iv`)を**切断を明示して**取り、
`S₂` を `Ψ` と揃えて取れば(`hsec`)、この四角形は**仮定なしで 1-可換**である。

★★★`Cor54SeamUnTr.lean` の測定「基準切断を揃えないと帳尻が合わない」への回答が
これである —— `modelType_equiv`(`Classical.choice` が切断を選ぶ)ではなく
`thm_5_2_iv` を使い、`S₂ := Ψ_*(S₁)`(在庫 `BaseSection.map`、`Corollary 5.7, (i)`)
と取れば、残っていた帳尻は**消える**。 -/
noncomputable def cSeamIso (Ψ : C₁ ⥤ C₂) (ΨB : D₁ ⥤ D₂)
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    (hdivc : ∀ {A B : C₁} (φ : A ⟶ B),
      phiIsoApp ΨB η ((P₁.toElem.obj A).base) (P₁.Div φ)
        = Φ₂.map (sq.hom.app A) (P₂.Div (Ψ.map φ)))
    (hdegc : ∀ {A B : C₁} (φ : A ⟶ B), P₁.degFr φ = P₂.degFr (Ψ.map φ))
    {G₁ : Frobenioid P₁} (R₁ : RatFnData P₁ G₁)
    (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ Z : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) Z)
    {S₁ : BaseSection P₁} (Fs₁ : ℕ+ →* SectionEnd S₁) (hFs₁ : IsFrobeniusSection S₁ Fs₁)
    {G₂ : Frobenioid P₂} (R₂ : RatFnData P₂ G₂)
    (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ Z : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) Z)
    {S₂ : BaseSection P₂} (Fs₂ : ℕ+ →* SectionEnd S₂) (hFs₂ : IsFrobeniusSection S₂ Fs₂)
    (hsec : ∀ {A : C₁}, S₁.objP A → S₂.objP (Ψ.obj A))
    (hPS : ∀ {A B : C₁} (f : A ⟶ B), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (Fo : ModelDataHomOver ΨB R₁.model R₂.model)
    (hphi : ∀ d : D₁, Fo.phiHom d = phiIsoApp ΨB η d)
    (hinj : ∀ d : D₂, Function.Injective (R₂.divB d)) :
    (thm_5_2_iv R₁ hiso₁ hfn₁ S₁ Fs₁ hFs₁).functor ⋙ Fo.functor
      ≅ Ψ ⋙ (thm_5_2_iv R₂ hiso₂ hfn₂ S₂ Fs₂ hFs₂).functor := by
  haveI := pathForget_isEquivalence G₁ S₁
  haveI := pathForget_isEquivalence G₂ S₂
  haveI := pathToModel_isEquivalence R₁ hiso₁ hfn₁ S₁ Fs₁ hFs₁
  haveI := pathToModel_isEquivalence R₂ hiso₂ hfn₂ S₂ Fs₂ hFs₂
  have hκ := conjIsoOfSquare (pathForget S₁) (pathForget S₂)
    (pathMapAlong Ψ hsec hPS) Ψ (pathMapAlong_pathForget Ψ hsec hPS)
  exact Functor.associator _ _ _
    ≪≫ Functor.isoWhiskerLeft (pathForget S₁).asEquivalence.inverse
        (pathSeamIso Ψ ΨB η sq hdivc hdegc R₁ hiso₁ hfn₁ Fs₁ hFs₁
          R₂ hiso₂ hfn₂ Fs₂ hFs₂ hsec hPS Fo hphi hinj)
    ≪≫ (Functor.associator _ _ _).symm
    ≪≫ Functor.isoWhiskerRight hκ (pathToModel R₂ hiso₂ hfn₂ S₂ Fs₂ hFs₂)
    ≪≫ Functor.associator _ _ _

end CSeam

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Corollary 5.4` の継ぎ目の本体
(`F-𝒫-path` の類は `Ψ` で運べる。★基準切断を揃えた場合の計算)。 -/
def FPPath.gpMap_cls_mapAlong.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 継ぎ目の本体(F-𝒫-path の類は Ψ で運べる)",
    sectionId := "frdi-cor-5-4" }

/-- ★★★★★★locator —— `Corollary 5.4` の縦の矢印の継ぎ目(`𝒞` の層、**条なし**)。

★`Theorem 5.2, (iv)` の圏同値を切断を明示して取り、`S₂` を `Ψ` と揃えれば
四角形は仮定なしで 1-可換になる。 -/
def cSeamIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印の継ぎ目(𝒞 の層、基準切断を揃えた形)",
    sectionId := "frdi-cor-5-4" }

/-- ★★★★locator —— 四角形を擬逆で共役する一般補題
(`Theorem 5.2, (iv)` の圏同値を `𝒞` の層へ降ろすのに要る)。 -/
def conjIsoOfSquare.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (iv) — 四角形を擬逆で共役する(𝒞̃ から 𝒞 へ降ろす)",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
