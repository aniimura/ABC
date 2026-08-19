/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55Rlf

/-!
# [FrdI] Proposition 5.3 の 1-可換図式 —— 縦の矢印を通す

原文 (FrdI p.103):
> [where the functor C →Cistr is the isotropification functor of Proposition 1.9, (v);

★`Proposition 5.3` の末尾の図式

```
𝒞       ⟶ 𝒞^istr       ⟶ 𝒞^pf
↓                          ↓
𝒞^un-tr ⟶ (𝒞^un-tr)^pf ⟶ 𝒞^rlf
```

のうち、**左の縦の矢印 `𝒞 ⥤ 𝒞^un-tr` と、そこから `𝒞^rlf` へ抜ける道**を通す。

★`𝒞` が isotropic 型なら `𝒞 ⥤ 𝒞^istr` は恒等的(`toIstrOfIsotropic`)、
`𝒞^istr ⥤ 𝒞^un-tr` は在庫にある(`istrToUnTr`、`Proposition 3.3, (iii)`)。
そこへ `untrToSc`(`Prop53Rlf.lean`)を継ぐ。

## ★★まだ実装していない条(記録)

* 上の行(`𝒞^istr ⟶ 𝒞^pf`)と右の縦の矢印 —— perfection がらみで、
  鎖 `prop55` の `p55i` / `p55ii-*` 待ちである。
* 図式の **1-可換性**そのもの(`𝒞 ⥤ 𝒞^rlf` の底への射影が `P.proj` と一致すること)。
  ★`pathToModel_toElem` が `rfl` なので出る見込みだが、まだ書いていない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {S : Type} [CommSemiring S]

/-- ★isotropic 型なら `𝒞 ⥤ 𝒞^istr` は恒等的
(`Proposition 1.9, (v)` の isotropification が恒等になる場合)。 -/
noncomputable def toIstrOfIsotropic (hiso : ∀ X : C, IsIsotropic P X) : C ⥤ Istr P where
  obj A := ⟨A, hiso A⟩
  map f := ObjectProperty.homMk f
  map_id _ := rfl
  map_comp _ _ := rfl

/-- ★★**`𝒞 ⥤ 𝒞^un-tr`** —— `Proposition 5.3` の図式の左の縦の矢印。 -/
noncomputable def cToUnTr (hiso : ∀ X : C, IsIsotropic P X) : C ⥤ UnTr P :=
  toIstrOfIsotropic hiso ⋙ istrToUnTr P

variable [IsConnected D]

/-- ★★★★★**`𝒞 ⥤ 𝒞^rlf`** —— `Proposition 5.3` の図式の
**左の縦の矢印 ＋ 下の行**の合成。

★`S = ℝ≥0` なら実化まで、`S = ℚ≥0` なら `(𝒞^un-tr)^pf` の行き先まで。 -/
noncomputable def cToSc (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D) :
    C ⥤ ScModelObj S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD :=
  cToUnTr hiso ⋙ untrToSc S Fc G hint hcharInj hintS hfsmD

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.3` の 1-可換図式の左の縦の矢印と、そこから実化へ抜ける道
(★**条つき**: 上の行(perfection)と 1-可換性そのものは未実装)。 -/
def cToSc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 1-可換図式の縦の矢印 𝒞 ⥤ 𝒞^un-tr ⥤ 𝒞^rlf",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
