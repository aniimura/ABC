/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55UntrFun
import ABC3.Found.FrdI.Prop53Diag

/-!
# [FrdI] Proposition 5.3 の図式 —— **四角形の 1-可換性**

原文 (FrdI p.103):
> [where the functor C → Cistr is the isotropification functor of Proposition 1.9,

★★`Proposition 5.3` の 1-可換図式

```
𝒞  ──→ 𝒞^istr ──→ 𝒞^pf
│                    │
↓                    ↓
𝒞^un-tr ─→ (𝒞^un-tr)^pf ─→ 𝒞^rlf
```

のうち、**四角形**(右の縦の矢印までの部分)の 1-可換性を示す。

## ★右の縦の矢印は何か

原文は「the remaining functors are the functors that arise naturally from the
construction of the 'unit-trivialization', 'perfection', and 'realification'」と言う。
★すなわち `𝒞^pf ⟶ 𝒞^rlf` は

  `𝒞^pf ⟶ (𝒞^pf)^un-tr ≌ (𝒞^un-tr)^pf ⟶ 𝒞^rlf`

である。★真ん中の圏同値が `Proposition 5.5, (ii)`(`unTrPfEquivalence`)。

## ★★四角形の中身は 1 本の補題

`𝒞` の射 `φ` を 2 通りに送ると、どちらも

  `toHomPf ((istrToUnTr P).map φ)`

になる。★これは `homPfToUnTrPf_toHomPf`(`Prop55UntrFun.lean`)そのものである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section Square

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)
  (F₁ : FrobenioidCore (istrPre P Fc)) (F₂ : FrobenioidCore (unTrPre P Fc))

/-- ★**`𝒞^pf ⟶ (𝒞^pf)^un-tr`** —— 図式の右の縦の矢印の第 1 段。 -/
noncomputable def pfToUnTrPf
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X) :
    PfCat (istrPre P Fc) F₁ ⥤ UnTr (pfPre (istrPre P Fc) F₁) :=
  cToUnTr (P := pfPre (istrPre P Fc) F₁) hisoPf

/-- ★★**`𝒞^pf ⟶ (𝒞^un-tr)^pf`** —— 図式の右の縦の矢印(`𝒞^rlf` の 1 つ手前まで)。

★`𝒞^pf ⟶ (𝒞^pf)^un-tr` と `Proposition 5.5, (ii)` の圏同値の合成である。 -/
noncomputable def pfToUnTrPfCat
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X) :
    PfCat (istrPre P Fc) F₁ ⥤ PfCat (unTrPre P Fc) F₂ :=
  pfToUnTrPf P Fc F₁ hisoPf ⋙ unTrPfFunctor P Fc F₁ F₂

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の図式の四角形は 1-可換**。

```
𝒞  ──→ 𝒞^istr ──→ 𝒞^pf
│                    │
↓                    ↓
𝒞^un-tr ────────→ (𝒞^un-tr)^pf
```

★中身は `homPfToUnTrPf_toHomPf` 1 本 —— どちらの道も
`φ ↦ toHomPf ((istrToUnTr P).map φ)` になる。 -/
theorem cToPf_comp_pfToUnTrPfCat_obj (hiso : ∀ X : C, IsIsotropic P X)
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X)
    (X : C) :
    (cToPf Fc hiso F₁ ⋙ pfToUnTrPfCat P Fc F₁ F₂ hisoPf).obj X
      = (cToUnTr hiso ⋙ toPfCat (unTrPre P Fc) F₂).obj X := by
  rfl

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の図式の四角形は 1-可換**(射の側)。

★中身は `homPfToUnTrPf_toHomPf` 1 本 —— どちらの道も
`φ ↦ toHomPf ((istrToUnTr P).map φ)` になる。 -/
theorem cToPf_comp_pfToUnTrPfCat_map (hiso : ∀ X : C, IsIsotropic P X)
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X)
    {X Y : C} (φ : X ⟶ Y) :
    (cToPf Fc hiso F₁ ⋙ pfToUnTrPfCat P Fc F₁ F₂ hisoPf).map φ
      = (cToUnTr hiso ⋙ toPfCat (unTrPre P Fc) F₂).map φ :=
  homPfToUnTrPf_toHomPf (Fc := Fc) (F₁ := F₁) (F₂ := F₂)
    ((toIstrOfIsotropic hiso).map φ)

/-- ★★★★★★★**[FrdI] Proposition 5.3 の 1-可換図式**(任意の尻尾つき)——

```
𝒞  ──→ 𝒞^istr ──→ 𝒞^pf ────────┐
│                               ↓
↓                          (𝒞^un-tr)^pf ──T──→ ?
𝒞^un-tr ──────────────────────┘
```

★上の道と下の道は、`(𝒞^un-tr)^pf` から先の関手 `T` が何であっても一致する。
★★`T := untrPf_modelFrobenioid ⋙ pfRootToRlfFunctor_of_arb` と取れば
**原文の図式そのもの**(尻尾は `𝒞^rlf`)になる。 -/
theorem cToPf_comp_pfToUnTrPfCat_map_comp {E : Type u2} [Category.{v2} E]
    (T : PfCat (unTrPre P Fc) F₂ ⥤ E) (hiso : ∀ X : C, IsIsotropic P X)
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X)
    {X Y : C} (φ : X ⟶ Y) :
    (cToPf Fc hiso F₁ ⋙ pfToUnTrPfCat P Fc F₁ F₂ hisoPf ⋙ T).map φ
      = (cToUnTr hiso ⋙ toPfCat (unTrPre P Fc) F₂ ⋙ T).map φ :=
  congrArg T.map (cToPf_comp_pfToUnTrPfCat_map P Fc F₁ F₂ hiso hisoPf φ)

end Square

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.3` の 1-可換図式の**四角形**。 -/
def cToPf_comp_pfToUnTrPfCat_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 1-可換図式の四角形(𝒞 ⥤ 𝒞^pf ⥤ (𝒞^un-tr)^pf と 𝒞 ⥤ 𝒞^un-tr ⥤ (𝒞^un-tr)^pf)",
    sectionId := "frdi-prop-5-3" }

/-- ★★★★locator —— 図式の右の縦の矢印(`𝒞^pf ⥤ (𝒞^un-tr)^pf`、
`Proposition 5.5, (ii)` の圏同値を経由する形)。 -/
def pfToUnTrPfCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 図式の右の縦の矢印 𝒞^pf ⥤ (𝒞^pf)^un-tr ≌ (𝒞^un-tr)^pf",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
