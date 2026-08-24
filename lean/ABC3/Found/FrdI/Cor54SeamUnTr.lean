/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54Seam
import ABC3.Found.FrdI.Cor54Compat
import ABC3.Found.FrdI.Cor54Rigid

/-!
# [FrdI] Corollary 5.4 の縦の矢印の継ぎ目 —— `𝒞^un-tr` の段を `𝔽_Φ` の上で繋ぐ

原文 (FrdI p.104):
> with the 1-commutative diagram of Proposition 5.3 [involving perfections, unit-

★`Cor54Seam.lean` の `ModelData.squareOfElemIso` は
「**`𝔽_Φ` の上の自然同型 `α` ＋ `(cls, u)` の帳尻**」から
model Frobenioid への 2 関手の同型を作る。本ファイルはその **`α` を実際に組み**、
継ぎ目を `(cls, u)` の帳尻 3 条だけに落とす。

## ★★`α` の組み方(4 段)

```
(𝒞₁^un-tr ⥤ Model₁ ⥤ Model₂) ⋙ 𝔽
  = 𝒞₁^un-tr ⥤ Model₁ ⋙ (Ψ^model ⋙ 𝔽)      -- 結合則
  ≅ 𝒞₁^un-tr ⥤ Model₁ ⋙ (𝔽 ⋙ Ψ_𝔽)          -- ★modelDataHomOverElemIso
  ≅ (𝒞₁^un-tr ⋙ 𝔽) ⋙ Ψ_𝔽                    -- ★modelTypeEquiv_comp_toElem
  ≅ Ψ^un-tr ⋙ (𝒞₂^un-tr ⋙ 𝔽)                -- ★★cor_4_11_iv_square
  = (Ψ^un-tr ⋙ 𝒞₂^un-tr ⥤ Model₂) ⋙ 𝔽       -- ★modelTypeEquiv_comp_toElem の逆
```

★★**4 段のうち 3 段は在庫**であり、新しいのは 1 段目
(`modelDataHomOverElemIso`、`Cor54Seam.lean`)だけである。
`untrModelHomOver` の `phiHom` は `phiIsoApp ΨB η d = (η.hom.app (op d)).hom` そのものなので、
その段の仮定 `hη` は **`rfl`** で埋まる。

## ★★★継ぎ目に残るもの(記録)

これで `base` / `div` / `deg` の側は**完全に閉じた**(`untrSeamElemIso`)。
`untrSeamIso` が示すとおり、残るのは **`uu` / `hc` / `hu` の 3 条**だけ ——
すなわち「`𝒞₁^un-tr` の対象の類 `cls` と射の `u` が、`Ψ` を通した後で
**`Φ^birat` の元 1 つ分だけずれる**」という**有理関数の側の帳尻**である。

## ★★★★★★測って分かったこと(2026-08-24)——**基準切断を揃えないと `hc` は偽**

★★残り 3 条を実際に潰そうとして、**この形のままでは `hc` は成り立たない**ことが分かった。
次に着手する人はここから始めてほしい。

`Theorem 5.2, (iv)` の対象の類は `Thm52Model.lean` の

```
pathObj R X = ⟨(P.toElem.obj X.obj).base, X.path.cls⟩
FPPath.cls  = Φ.gpMapOn (inv (Base toObj)) (toGp (Div toObj) - toGp (Div toRef))
```

であり、★**`toRef`(基準対象 `ref ∈ 𝒫` への pre-step)に依存する**。
したがって `cls` は **`BaseSection S` の取り方に依存**する。

* `H₁ ⋙ Fo.functor` の類 …… `S₁` から測ったもの(を `η` で送ったもの)
* `Ψu ⋙ H₂` の類 …………… `S₂` から測ったもの

`hdivc` は `Div` の対応を与えるので、差は
**`Ψu (S₁ d)` と `S₂ (ΨB d)` の間の span の類 `γ(d)`** ちょうどになる。
★`hc` は「`γ(d) ∈ Φ₂^birat`」を要求するが、
`Φ₂^birat ⊊ Gp(Φ₂)`(算術の例では類群が非自明)である以上、
★★**`S₁` と `S₂` を無関係に選べば `γ(d)` は `Φ₂^birat` の外に出る。**

★★★そして `unTr_modelFrobenioid` は `modelType_equiv` 経由で
**`Classical.choice hm.1` が切断を選ぶ**(`Cor54Rigid.lean` の
`modelTypeEquiv_comp_toElem` を見よ)ので、両側の切断は**無関係**である。

### ★直し方(在庫あり)

**`modelType_equiv` ではなく `thm_5_2_iv` を使い、切断を明示的に渡す。**
そのうえで `S₂ := S₁.map Ψu …`(★在庫 `BaseSection.map`、`Cor57Base.lean`、
`Corollary 5.7, (i)`)と取れば、`Ψu (S₁ d)` は `S₂` の `𝒫` の対象そのものになり、
基準が揃って `γ(d)` が消える(あるいは span の類として明示的に `uu` を与える)。

★`BaseSection.map` の仮定は原文 `Corollary 5.7, (i)` のもの ——
底の 1-可換図式 `BaseSquare`、pull-back 射の保存、Frobenius-trivial 性の保存。
★★これは `Corollary 5.4` の入力(`Corollary 4.11, (ii)(iii)(iv)`)と同じ層にあり、
**新しい数学は要らない**。要るのは `Theorem 5.2, (iv)` の切断を
`modelType_equiv` の内側から外へ出す作業である。

★なお `divB` は包含そのものなので、`uu` ＋ `hc` は
「差が `Φ₂^birat` に入る」という **membership 1 条**に畳める。
★★`cor_4_11_ii_uniq` は**底の関手 `Gf : 𝒟₁ ⥤ 𝒟₂` の** 1-一意性であって
`𝒞^un-tr` の関手の一意性ではないので、そこからは降りてこない(測定済み)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section SeamUnTr

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  [IsConnected D₁] [IsConnected D₂]

variable (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
  (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
  (hint₁ : ∀ A : D₁, IsIntegralMonoid (Φ₁.val A)) (hfsmD₁ : IsOfFSMType D₁)
  (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
  (hint₂ : ∀ A : D₂, IsIntegralMonoid (Φ₂.val A)) (hfsmD₂ : IsOfFSMType D₂)
  (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)),
    y ∈ phiBiratOn (unTr_frobenioid P₁ Fc₁ G₁) d →
    gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn (unTr_frobenioid P₂ Fc₂ G₂) (ΨB.obj d))
  (Ψu : UnTr P₁ ⥤ UnTr P₂)
  (sq : (unTrPre P₁ Fc₁).proj ⋙ ΨB ≅ Ψu ⋙ (unTrPre P₂ Fc₂).proj)
  (hdivc : ∀ {A B : UnTr P₁} (φ : A ⟶ B),
    (η.hom.app (Opposite.op (((unTrPre P₁ Fc₁).toElem.obj A).base))).hom
        ((unTrPre P₁ Fc₁).Div φ)
      = Φ₂.map (sq.hom.app A) ((unTrPre P₂ Fc₂).Div (Ψu.map φ)))
  (hdeg : ∀ {A B : UnTr P₁} (φ : A ⟶ B),
    (unTrPre P₂ Fc₂).degFr (Ψu.map φ) = (unTrPre P₁ Fc₁).degFr φ)

/-- ★★★★★★**[FrdI] Corollary 5.4 の縦の矢印の継ぎ目の `𝔽_Φ` の側** ——

`𝒞^un-tr` の段の 2 つの関手は、**`𝔽_Φ` へ落とせば自然同型**である。

★入力は `Corollary 4.11, (iii)(iv)` が与えるもの(`ΨB` / `η` / `sq` / `hdivc` / `hdeg`)だけ。 -/
noncomputable def untrSeamElemIso :
    ((unTr_modelFrobenioid Fc₁ G₁ hint₁ hfsmD₁).functor
        ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor)
        ⋙ ModelData.toElem
      ≅ (Ψu ⋙ (unTr_modelFrobenioid Fc₂ G₂ hint₂ hfsmD₂).functor) ⋙ ModelData.toElem :=
  Functor.associator _ _ _
    ≪≫ Functor.isoWhiskerLeft _
        (modelDataHomOverElemIso
          (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat)
          η.hom (fun _ _ => rfl))
    ≪≫ (Functor.associator _ _ _).symm
    ≪≫ Functor.isoWhiskerRight
        (modelTypeEquiv_comp_toElem (unTr_ratFnData Fc₁ G₁ hint₁ hfsmD₁)
          (unTr_isotropic P₁ Fc₁) (unTr_isOfModelType Fc₁ G₁))
        (elemFrobMapOver ΨB η.hom)
    ≪≫ cor_4_11_iv_square Ψu ΨB sq η.hom hdivc hdeg
    ≪≫ Functor.isoWhiskerLeft Ψu
        (modelTypeEquiv_comp_toElem (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂)
          (unTr_isotropic P₂ Fc₂) (unTr_isOfModelType Fc₂ G₂)).symm
    ≪≫ (Functor.associator _ _ _).symm

/-- ★★`untrSeamElemIso` が与える**底の同型**。 -/
noncomputable abbrev untrSeamBase (X : UnTr P₁) :
    ((((unTr_modelFrobenioid Fc₁ G₁ hint₁ hfsmD₁).functor
        ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁
            Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor)).obj X).base
      ≅ (((Ψu ⋙ (unTr_modelFrobenioid Fc₂ G₂ hint₂ hfsmD₂).functor)).obj X).base :=
  ModelData.elemIsoBase
    (untrSeamElemIso ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat Ψu sq hdivc hdeg) X

/-- ★★★★★★**[FrdI] Corollary 5.4 の縦の矢印の継ぎ目** ——

**残るのは有理関数の側の帳尻 `uu` / `hc` / `hu` の 3 条だけ**である。

★`base` / `div` / `deg` の側は `untrSeamElemIso` が閉じている
(`Corollary 4.11, (iii)(iv)` ＋ `Theorem 5.2, (iv)` の圏同値の 1-可換性)。
★★`neg` は `RatFnData` の `bneg`(`B` が group-like であること)、
`hdivOn` は `P₂` の divisorial 性からそのまま取れる。 -/
noncomputable def untrSeamIso
    (uu : ∀ X : UnTr P₁,
      (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).model.bmon.val
        (((((unTr_modelFrobenioid Fc₁ G₁ hint₁ hfsmD₁).functor
          ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁
              Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor)).obj X).base))
    (hc : ∀ X : UnTr P₁,
      ((((unTr_modelFrobenioid Fc₁ G₁ hint₁ hfsmD₁).functor
          ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁
              Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor)).obj X).cls
        = (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).model.phi.gpMapOn
            (untrSeamBase ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat
              Ψu sq hdivc hdeg X).hom
            (((Ψu ⋙ (unTr_modelFrobenioid Fc₂ G₂ hint₂ hfsmD₂).functor).obj X).cls)
          + (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).model.divB _ (uu X))
    (hu : ∀ {X Y : UnTr P₁} (f : X ⟶ Y),
      (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).model.bmon.map
          ((((unTr_modelFrobenioid Fc₁ G₁ hint₁ hfsmD₁).functor
            ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁
                Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor)).map f).base (uu Y)
        + (((unTr_modelFrobenioid Fc₁ G₁ hint₁ hfsmD₁).functor
            ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁
                Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor).map f).u
        = (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).model.bmon.map
            (untrSeamBase ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat
              Ψu sq hdivc hdeg X).hom
            (((Ψu ⋙ (unTr_modelFrobenioid Fc₂ G₂ hint₂ hfsmD₂).functor).map f).u)
          + (((((Ψu ⋙ (unTr_modelFrobenioid Fc₂ G₂ hint₂ hfsmD₂).functor).map f).deg
              : ℕ+) : ℕ)) • (uu X)) :
    (unTr_modelFrobenioid Fc₁ G₁ hint₁ hfsmD₁).functor
        ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor
      ≅ Ψu ⋙ (unTr_modelFrobenioid Fc₂ G₂ hint₂ hfsmD₂).functor :=
  ModelData.squareOfElemIso (fun A => P₂.divisorial A)
    (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).bneg
    (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂).bneg_add
    (untrSeamElemIso ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat Ψu sq hdivc hdeg)
    uu hc hu

end SeamUnTr

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Corollary 5.4` の縦の矢印の継ぎ目の `𝔽_Φ` の側
(`base` / `div` / `deg` の 3 成分。★**条つき**: `Corollary 4.11, (iii)(iv)` の
入力 `sq` / `hdivc` / `hdeg` を仮定として受け取っている)。 -/
def untrSeamElemIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印の継ぎ目(𝒞^un-tr の段を 𝔽_Φ の上で繋ぐ)",
    sectionId := "frdi-cor-5-4" }

/-- ★★★★★★locator —— `Corollary 5.4` の縦の矢印の継ぎ目
(★**条つき**: 残るのは有理関数の側の帳尻 `uu` / `hc` / `hu` の 3 条だけ)。 -/
def untrSeamIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印の継ぎ目(有理関数の帳尻 3 条へ還元)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
