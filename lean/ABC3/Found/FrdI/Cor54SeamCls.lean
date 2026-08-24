/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54Birat
import ABC3.Found.FrdI.Thm52Model

/-!
# [FrdI] Corollary 5.4 の継ぎ目の本体 —— `F-𝒫-path` の類は `Ψ` で運べる

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

★★★これで「基準を揃えれば `γ(d) = 0`」が**計算として**確かめられた ——
`Cor54SeamUnTr.lean` の `hex` は、`S₂ := S₁.map Ψ`(在庫 `BaseSection.map`、
`Corollary 5.7, (i)`)と取れば `u = 0` で満たせる、というのが本ファイルの帰結である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

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

end ClsMap

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Corollary 5.4` の継ぎ目の本体
(`F-𝒫-path` の類は `Ψ` で運べる。★基準切断を揃えた場合の計算)。 -/
def FPPath.gpMap_cls_mapAlong.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 継ぎ目の本体(F-𝒫-path の類は Ψ で運べる)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
