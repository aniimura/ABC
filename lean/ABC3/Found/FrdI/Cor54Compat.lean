/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54
import ABC3.Found.FrdI.Cor54Birat
import ABC3.Found.FrdI.Prop53Base
import ABC3.Found.FrdI.Prop53Diag

/-!
# [FrdI] Corollary 5.4 の最後の文 —— `Ψ^rlf` の作り方は完全化と 1-可換

原文 (FrdI p.104):
> with the 1-commutative diagram of Proposition 5.3 [involving perfections, unit-

★`Corollary 5.4` の最後の文は「`Ψ` から `Ψ^rlf` を作る操作は、
`Proposition 5.3` の 1-可換図式(完全化・単位自明化などが出てくるもの)と
**1-可換である**」と言う。

★★我々の書き方では、`Proposition 5.3` の図式の下の行の右の矢印は
**係数の拡大 `σ : S ⟶ S'`**(`Prop53Base.lean` の `scBaseFunctor`)であり、
`Ψ^rlf` は**因子の単系の同型 `η` を `⊗_S` したもの**(`Cor54.lean` の `psiSc`)である。
したがって主張は次の四角形の可換性になる:

```
(𝒞₁^un-tr)^pf --Ψ^pf--> (𝒞₂^un-tr)^pf
     |                        |
   係数の拡大               係数の拡大
     v                        v
   𝒞₁^rlf     --Ψ^rlf-->    𝒞₂^rlf
```

★★★**中身は `scBase σ` と `scMap θ` が可換であること 1 本だけ**である ——
どちらもテンソル積の別々の因子に効くのだから当たり前で、
それが `scBase_scMap`(`Prop53Base.lean`)である。

## ★本ファイルが閉じるもの

1. 最後の文(**完全化・実化との 1-可換性**)—— `scBaseFunctor_comp_psiSc`。
2. `Proposition 5.3` の**縦の矢印の `𝒞^un-tr` の段**との 1-可換性 ——
   `untrToSc_comp_psiSc`。★中身は `scMap θ (1 ⊗ m) = 1 ⊗ θ m` だけである。
3. 縦の矢印の**いちばん上の段** `𝒞 ⥤ 𝒞^un-tr` との 1-可換性 ——
   `cToUnTr_comp_psiUnTr`。★`psiUnTr_square`(`Theorem 3.4, (iv)`)に
   `toIstrOfIsotropic_comp_psiIstr`(両辺が**同じ項**になるので `rfl`)を継ぐだけ。

## ★★まだ実装していない条(記録)

`Corollary 5.4` の残りは

* **1-一意性**(`Ψ^rlf` が同型を除いて一意)
* **rigidity**(図式の合成関手がすべて rigid)
* ★**縦の矢印を 1 本に繋ぐ最後の継ぎ目** ——
  `𝒞^un-tr ≌ Model(Φ, Φ^birat)`(`unTr_modelFrobenioid`)が
  `Ψ^un-tr` と `untrModelHomOver` を**取り違えない**こと、すなわち

    `(unTr_modelFrobenioid Fc₁ G₁).functor ⋙ (untrModelHomOver …).functor`
      `= psiUnTr Ψ … ⋙ (unTr_modelFrobenioid Fc₂ G₂).functor`

  これは `Corollary 4.11, (iii)(iv)` の内容そのものであり、本ファイルの範囲外。

  ★★★**2026-08-24 の分析(次に着手する人へ)**: 両辺は
  `UnTr P₁ ⥤ ModelData.Obj (unTr_ratFnData Fc₂ G₂ …)` である。
  model Frobenioid の射は **4 成分 `(base, div, deg, u)`** で、
  `toElem` が見るのは前 3 つだけ、`u` は有理関数の成分である。
  ★したがって継ぎ目は 2 つに割れる:
  * **前 3 成分** —— `Corollary 4.11, (iv)` の 1-可換図式(在庫 `cor_4_11_iv_square`)と
    model 型の圏同値が `𝔽` と 1-可換であること(在庫 `modelTypeEquiv_comp_toElem`)。
  * **`u` の成分** —— `Φ^birat` の輸送(在庫 `phiBiratOn_transport_of_cor411`)。
  ★★要るのは「**model Frobenioid への 2 つの関手が
  `toElem` の合成と `B` 成分で一致すれば等しい**」という 1 本である。
  ★1-一意性の側(`Cor54Uniq.lean`)は
  `natTrans_ext_of_essSurj` / `phiIso_ext` で `(Ψ_𝒟, η)` へ帰着済み。

★`hbirat` は `Cor54Birat.lean` の `phiBiratOn_transport_of_cor411` で
**すでに `Corollary 4.10`/`4.11` から導けている**(`Cor54.lean` 冒頭の記述は古い)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

section Compat

variable {S S' : Type} [CommSemiring S] [CommSemiring S']

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★★**`Corollary 5.4` の最後の文**(`ModelData` の射の段)——
「`Ψ` から `Ψ^rlf` を作る操作」と「係数の拡大」は可換。

★中身は `scBase σ (scMap θ x) = scMap θ (scBase σ x)`(`scBase_scMap`)だけである。 -/
theorem scBaseHom_compOver_scModelHomOver (σ : S →+* S')
    (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (G₁ : Frobenioid P₁) (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
    (hcharInj₁ : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₁.map α)))
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (ScT S (Φ₁.val A)))
    (hcharInj₁' : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₁.map α)))
    (hint₁' : ∀ A : D₁, IsIntegralMonoid (ScT S' (Φ₁.val A)))
    (G₂ : Frobenioid P₂) (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
    (hcharInj₂ : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₂.map α)))
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (ScT S (Φ₂.val A)))
    (hcharInj₂' : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₂.map α)))
    (hint₂' : ∀ A : D₂, IsIntegralMonoid (ScT S' (Φ₂.val A)))
    (hfsmD₁ : IsOfFSMType D₁) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)), y ∈ phiBiratOn G₁ d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d)) :
    (scBaseHom σ G₁ hiso₁ hfn₁ hcharInj₁ hint₁ hcharInj₁' hint₁' hfsmD₁).compOver
        (scModelHomOver S' ΨB η G₁ hiso₁ hfn₁ hcharInj₁' hint₁'
          G₂ hiso₂ hfn₂ hcharInj₂' hint₂' hfsmD₁ hfsmD₂ hbirat)
      = (scModelHomOver S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁
          G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hfsmD₁ hfsmD₂ hbirat).compHom
        (scBaseHom σ G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hcharInj₂' hint₂' hfsmD₂) := by
  refine ModelDataHomOver.ext' (fun d => ?_) (fun d => ?_)
  · exact AddMonoidHom.ext (fun x => (scBase_scMap σ (phiIsoApp ΨB η d) x).symm)
  · exact AddMonoidHom.ext (fun x =>
      Subtype.ext (gpMap_scBase_scMap σ (phiIsoApp ΨB η d) _).symm)

/-- ★★★★★★**`Corollary 5.4` の最後の文**(関手の等式として)。

```
(𝒞₁^un-tr)^pf --Ψ^pf--> (𝒞₂^un-tr)^pf --係数の拡大--> 𝒞₂^rlf
        ＝
(𝒞₁^un-tr)^pf --係数の拡大--> 𝒞₁^rlf --Ψ^rlf--> 𝒞₂^rlf
```
-/
theorem scBaseFunctor_comp_psiSc (σ : S →+* S')
    (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (G₁ : Frobenioid P₁) (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
    (hcharInj₁ : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₁.map α)))
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (ScT S (Φ₁.val A)))
    (hcharInj₁' : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₁.map α)))
    (hint₁' : ∀ A : D₁, IsIntegralMonoid (ScT S' (Φ₁.val A)))
    (G₂ : Frobenioid P₂) (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
    (hcharInj₂ : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₂.map α)))
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (ScT S (Φ₂.val A)))
    (hcharInj₂' : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₂.map α)))
    (hint₂' : ∀ A : D₂, IsIntegralMonoid (ScT S' (Φ₂.val A)))
    (hfsmD₁ : IsOfFSMType D₁) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)), y ∈ phiBiratOn G₁ d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d)) :
    scBaseFunctor σ G₁ hiso₁ hfn₁ hcharInj₁ hint₁ hcharInj₁' hint₁' hfsmD₁
        ⋙ psiSc S' ΨB η G₁ hiso₁ hfn₁ hcharInj₁' hint₁'
          G₂ hiso₂ hfn₂ hcharInj₂' hint₂' hfsmD₁ hfsmD₂ hbirat
      = psiSc S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁
          G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hfsmD₁ hfsmD₂ hbirat
        ⋙ scBaseFunctor σ G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hcharInj₂' hint₂' hfsmD₂ := by
  unfold scBaseFunctor psiSc
  rw [← ModelDataHom.compOver_functor, ← ModelDataHomOver.compHom_functor,
    scBaseHom_compOver_scModelHomOver]

end Compat

/-! ## ★`Proposition 5.3` の縦の矢印(`𝒞^un-tr` の段)との 1-可換性 -/

section VerticalUnTr

variable {S : Type} [CommSemiring S]

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  [IsConnected D₁] [IsConnected D₂]

/-- ★★★★**`Ψ` が `𝒞^un-tr` の model data のあいだに誘導する射**
`(Φ₁, Φ₁^birat) ⟶_{Ψ_𝒟} (Φ₂, Φ₂^birat)`。

★`scModelHomOver` の「係数をつけない版」である ——
`Corollary 4.11` から来る `η`(因子の単系の同型)と `hbirat`(`Φ^birat` を保つこと)
だけで組める。 -/
noncomputable def untrModelHomOver (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (Φ₁.val A)) (hfsmD₁ : IsOfFSMType D₁)
    (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (Φ₂.val A)) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)),
      y ∈ phiBiratOn (unTr_frobenioid P₁ Fc₁ G₁) d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn (unTr_frobenioid P₂ Fc₂ G₂) (ΨB.obj d)) :
    ModelDataHomOver ΨB (unTr_ratFnData Fc₁ G₁ hint₁ hfsmD₁).model
      (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).model where
  phiHom d := phiIsoApp ΨB η d
  phiNat f x := phiIsoApp_nat ΨB η f x
  bmonHom d :=
    AddMonoidHom.codRestrict
      ((gpMap _ (phiIsoApp ΨB η d)).comp
        (phiBiratOn (unTr_frobenioid P₁ Fc₁ G₁) d).subtype) _
      (fun x => hbirat d x x.2)
  bmonNat := fun f x => Subtype.ext (gpMap_phiIsoApp_nat ΨB η f _)
  divCompat _ _ := rfl

/-- ★★★★★★**`Proposition 5.3` の縦の矢印(`𝒞^un-tr ⥤ 𝒞^rlf`)は `Ψ` と 1-可換**
(`ModelData` の射の段)。

★中身は `scMap θ (1 ⊗ m) = 1 ⊗ θ m`(`scMap_toSc`)だけである。 -/
theorem untrToScHom_compOver_scModelHomOver
    (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (Φ₁.val A))
    (hcharInj₁ : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₁.map α)))
    (hintS₁ : ∀ A : D₁, IsIntegralMonoid (ScT S (Φ₁.val A))) (hfsmD₁ : IsOfFSMType D₁)
    (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (Φ₂.val A))
    (hcharInj₂ : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₂.map α)))
    (hintS₂ : ∀ A : D₂, IsIntegralMonoid (ScT S (Φ₂.val A))) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)),
      y ∈ phiBiratOn (unTr_frobenioid P₁ Fc₁ G₁) d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn (unTr_frobenioid P₂ Fc₂ G₂) (ΨB.obj d)) :
    (untrToScHom S Fc₁ G₁ hint₁ hcharInj₁ hintS₁ hfsmD₁).compOver
        (scModelHomOver S ΨB η
          (unTr_frobenioid P₁ Fc₁ G₁) (unTr_isotropic P₁ Fc₁)
            (fun Z => (unTr_isOfModelType Fc₁ G₁).2 Z) hcharInj₁ hintS₁
          (unTr_frobenioid P₂ Fc₂ G₂) (unTr_isotropic P₂ Fc₂)
            (fun Z => (unTr_isOfModelType Fc₂ G₂).2 Z) hcharInj₂ hintS₂
          hfsmD₁ hfsmD₂ hbirat)
      = (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat).compHom
        (untrToScHom S Fc₂ G₂ hint₂ hcharInj₂ hintS₂ hfsmD₂) := by
  refine ModelDataHomOver.ext' (fun d => ?_) (fun d => ?_)
  · exact AddMonoidHom.ext (fun x => scMap_toSc (phiIsoApp ΨB η d) x)
  · exact AddMonoidHom.ext (fun x =>
      Subtype.ext (gpMap_scMap_toScGp (phiIsoApp ΨB η d) _))

/-- ★★★★★★**同上、関手の等式として**。

```
𝒞₁^un-tr --Ψ^un-tr--> 𝒞₂^un-tr --> 𝒞₂^rlf
      ＝
𝒞₁^un-tr --> 𝒞₁^rlf --Ψ^rlf--> 𝒞₂^rlf
```
-/
theorem untrToSc_comp_psiSc (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (Φ₁.val A))
    (hcharInj₁ : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₁.map α)))
    (hintS₁ : ∀ A : D₁, IsIntegralMonoid (ScT S (Φ₁.val A))) (hfsmD₁ : IsOfFSMType D₁)
    (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (Φ₂.val A))
    (hcharInj₂ : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₂.map α)))
    (hintS₂ : ∀ A : D₂, IsIntegralMonoid (ScT S (Φ₂.val A))) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)),
      y ∈ phiBiratOn (unTr_frobenioid P₁ Fc₁ G₁) d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn (unTr_frobenioid P₂ Fc₂ G₂) (ΨB.obj d)) :
    (untrToScHom S Fc₁ G₁ hint₁ hcharInj₁ hintS₁ hfsmD₁).functor
        ⋙ psiSc S ΨB η
          (unTr_frobenioid P₁ Fc₁ G₁) (unTr_isotropic P₁ Fc₁)
            (fun Z => (unTr_isOfModelType Fc₁ G₁).2 Z) hcharInj₁ hintS₁
          (unTr_frobenioid P₂ Fc₂ G₂) (unTr_isotropic P₂ Fc₂)
            (fun Z => (unTr_isOfModelType Fc₂ G₂).2 Z) hcharInj₂ hintS₂
          hfsmD₁ hfsmD₂ hbirat
      = (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor
        ⋙ (untrToScHom S Fc₂ G₂ hint₂ hcharInj₂ hintS₂ hfsmD₂).functor := by
  unfold psiSc
  rw [← ModelDataHom.compOver_functor, ← ModelDataHomOver.compHom_functor,
    untrToScHom_compOver_scModelHomOver]

end VerticalUnTr

/-! ## ★縦の矢印のいちばん上の段 `𝒞 ⥤ 𝒞^un-tr` -/

section VerticalTop

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence]

/-- ★★isotropic 型なら `𝒞 ⥤ 𝒞^istr` は `Ψ^istr` と**厳密に**可換。

★対象は `⟨Ψ.obj A, _⟩` どうし(所属は `Prop` なので証明無関係)、
射は `Ψ.map` どうしで、両辺は**同じ項**になる。 -/
theorem toIstrOfIsotropic_comp_psiIstr (h₁ : IsOfQuasiIsotropicType C₁ P₁)
    (h₂ : IsOfQuasiIsotropicType C₂ P₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) :
    toIstrOfIsotropic hiso₁ ⋙ psiIstr Ψ P₁ P₂ h₁ h₂
      = Ψ ⋙ toIstrOfIsotropic hiso₂ := rfl

/-- ★★★★★**`Proposition 5.3` の縦の矢印のいちばん上の段 `𝒞 ⥤ 𝒞^un-tr` は
`Ψ` と 1-可換**(厳密に)。

★`psiUnTr_square`(`Theorem 3.4, (iv)` の四角形)に
`toIstrOfIsotropic_comp_psiIstr` を継ぐだけ。 -/
theorem cToUnTr_comp_psiUnTr (h₁ : IsOfQuasiIsotropicType C₁ P₁)
    (h₂ : IsOfQuasiIsotropicType C₂ P₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hUE : ∀ {A B : C₁} (α₁ α₂ : A ⟶ B),
      P₁.toElem.map α₁ = P₁.toElem.map α₂ →
        P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂)) :
    cToUnTr hiso₁ ⋙ psiUnTr Ψ h₁ h₂ hUE = Ψ ⋙ cToUnTr hiso₂ := by
  show toIstrOfIsotropic hiso₁ ⋙ (istrToUnTr P₁ ⋙ psiUnTr Ψ h₁ h₂ hUE) = _
  rw [← psiUnTr_square Ψ h₁ h₂ hUE]
  rfl

end VerticalTop

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Corollary 5.4` の最後の文
「`Ψ^rlf` の作り方は完全化・実化と 1-可換」。 -/
def scBaseFunctor_comp_psiSc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — Ψ^rlf の作り方は完全化・実化と 1-可換",
    sectionId := "frdi-cor-5-4" }

/-- ★★★★locator —— `Corollary 5.4` の 1-可換図式のうち
`Proposition 5.3` の縦の矢印の `𝒞^un-tr` の段。 -/
def untrToSc_comp_psiSc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印(𝒞^un-tr ⥤ 𝒞^rlf)は Ψ と 1-可換",
    sectionId := "frdi-cor-5-4" }

/-- ★★★★locator —— `Corollary 5.4` の 1-可換図式のうち
`Proposition 5.3` の縦の矢印のいちばん上の段 `𝒞 ⥤ 𝒞^un-tr`。 -/
def cToUnTr_comp_psiUnTr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印(𝒞 ⥤ 𝒞^un-tr)は Ψ と 1-可換",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
