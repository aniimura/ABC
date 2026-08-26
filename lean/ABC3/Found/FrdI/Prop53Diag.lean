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

/-! ## ★1-可換性と上の行 -/

omit [IsConnected D] in
/-- ★★★**左の縦の矢印は `𝔽_Φ` への射影と 1-可換**。

★★これが `Proposition 5.3` の図式の**可換性の中身**である ——
`ƒ → 𝒞^un-tr → 𝔽_Φ` と `𝒞 → 𝔽_Φ` が一致する。
★`istrToUnTr_comp_unTrToElem`(`UnTr.lean`)と、
`toIstrOfIsotropic` が `ι` で戻ること(定義から)の 2 本で `rfl` になる。 -/
theorem cToUnTr_comp_unTrToElem (hiso : ∀ X : C, IsIsotropic P X) :
    cToUnTr hiso ⋙ unTrToElem P = P.toElem := rfl

omit [IsConnected D] in
/-- ★★**底への射影とも 1-可換**。 -/
theorem cToUnTr_comp_proj (hiso : ∀ X : C, IsIsotropic P X) :
    cToUnTr hiso ⋙ (unTrToElem P ⋙ ElemFrobCat.proj) = P.proj := rfl

/-! ## ★図式の上の行 `𝒞 ⟶ 𝒞^istr ⟶ 𝒞^pf` -/

/-- ★★★**図式の上の行** —— `𝒞 ⟶ 𝒞^istr ⟶ 𝒞^pf`。

★isotropic 型なら `𝒞 ⟶ 𝒞^istr` は恒等的(`toIstrOfIsotropic`)、
`𝒞^istr ⟶ 𝒞^pf` は在庫(`toPfCat`、`Definition 3.1, (iii)`)。 -/
noncomputable def cToPf (Fc : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    (F₁ : FrobenioidCore (istrPre P Fc)) : C ⥤ PfCat (istrPre P Fc) F₁ :=
  toIstrOfIsotropic hiso ⋙ toPfCat (istrPre P Fc) F₁

omit [IsConnected D] in
/-- ★★**上の行も `𝔽_{Φ^pf}` への射影と 1-可換**
(`Proposition 3.2, (i)` の 1-可換図式を `𝒞` から見たもの)。 -/
theorem cToPf_comp_pfToElem (Fc : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    (F₁ : FrobenioidCore (istrPre P Fc)) :
    cToPf Fc hiso F₁ ⋙ pfToElem (istrPre P Fc) F₁
      = toIstrOfIsotropic hiso ⋙ (istrPre P Fc).toElem ⋙ elemToPfElem (istrPre P Fc) := by
  rw [cToPf, Functor.assoc]
  exact congrArg (fun G => toIstrOfIsotropic hiso ⋙ G)
    (Functor.ext_of_iso (pfSquare (istrPre P Fc) F₁) (fun _ => rfl) (fun _ => rfl))

/-! ## ★`Proposition 5.5, (iii)` の observation

原文 (FrdI p.105):
> morphisms of Cun-tr, Crlf are precisely the linear isometries [cf. Proposition 1.4,

★原文が挙げる 2 つの observation の 1 つ。★`𝒞^un-tr` / model Frobenioid については
すでに在庫がある(`unTr_isPullBack_iff` / `ModelData.model_isPullBack_iff`)ので、
ここでは**isotropic 型の Frobenioid 一般**の形に束ねておく。 -/

omit [IsConnected D] in
/-- ★★★**isotropic 型では「pull-back 射 ⟺ linear な等長射」**。

★`⟹` は `Definition 1.3, (iv)` の `pullBackLB`、
`⟸` は `Proposition 1.4, (ii)`(co-angular 性は `Proposition 1.4, (i)` から自動)。 -/
theorem isPullBack_iff_linear_isometric (Fc : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (φ : A ⟶ B) :
    IsPullBack P φ ↔ (IsLinear P φ ∧ IsIsometric P φ) := by
  constructor
  · intro h
    obtain ⟨hlb, hlin⟩ := Fc.pullBackLB φ h
    exact ⟨hlin, hlb.2⟩
  · rintro ⟨hlin, hisom⟩
    exact prop_1_4_ii_mpr P Fc φ ⟨prop_1_4_i P φ (fun Y _ => hiso Y), hisom⟩ hlin

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.3` の 1-可換図式の左の縦の矢印と、そこから実化へ抜ける道
(★**条つき**: 上の行(perfection)と 1-可換性そのものは未実装)。 -/
def cToSc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 1-可換図式の縦の矢印 𝒞 ⥤ 𝒞^un-tr ⥤ 𝒞^rlf",
    sectionId := "frdi-prop-5-3" }

/-- ★locator —— `Proposition 5.5, (iii)` の証明が挙げる observation
「pull-back 射はちょうど linear な等長射」。 -/
def isPullBack_iff_linear_isometric.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — pull-back 射はちょうど linear な等長射",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
