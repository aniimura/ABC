/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55UntrFun
import ABC3.Found.FrdI.Prop53Diag
import ABC3.Found.FrdI.Prop53UntrPfModel
import ABC3.Found.FrdI.Prop53PfCatRoot

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

universe v u w u2 v2 ue ve

section Square

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)
  (F₁ : FrobenioidCore (istrPre P Fc)) (F₂ : FrobenioidCore (unTrPre P Fc))

/-- ★★★**`𝒞^pf`(根 1 の部分)は isotropic 型** —— 仮引数で受ける必要はない。

★`𝒞^istr` は isotropic 型(`istr_isOfIsotropicType`)なので Frobenius-isotropic 型でもあり
(`isOfFrobeniusIsotropicType_of_isotropic`)、`Proposition 3.2, (iii)` を
根 1 の部分へ降ろした `pf_isOfIsotropicType` がそのまま当たる。 -/
theorem pfIstr_isotropic (X : PfCat (istrPre P Fc) F₁) :
    IsIsotropic (pfPre (istrPre P Fc) F₁) X :=
  pf_isOfIsotropicType
    (isOfFrobeniusIsotropicType_of_isotropic (istr_isOfIsotropicType (F := Fc))) X

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
theorem cToPf_comp_pfToUnTrPfCat_map_comp {E : Type ue} [Category.{ve} E]
    (T : PfCat (unTrPre P Fc) F₂ ⥤ E) (hiso : ∀ X : C, IsIsotropic P X)
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X)
    {X Y : C} (φ : X ⟶ Y) :
    (cToPf Fc hiso F₁ ⋙ pfToUnTrPfCat P Fc F₁ F₂ hisoPf ⋙ T).map φ
      = (cToUnTr hiso ⋙ toPfCat (unTrPre P Fc) F₂ ⋙ T).map φ :=
  congrArg T.map (cToPf_comp_pfToUnTrPfCat_map P Fc F₁ F₂ hiso hisoPf φ)

/-! ## ★2. 右の縦の矢印を `𝒞^rlf` まで束ねる -/

section ToRlf

variable [IsConnected D] (G : Frobenioid P)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の図式の右の縦の矢印(全部)** `𝒞^pf ⟶ 𝒞^rlf`。

```
𝒞^pf ⟶ (𝒞^pf)^un-tr ≌ (𝒞^un-tr)^pf ⟶ model ⟶ 𝒞^rlf
```

★4 段の合成である:
1. `pfToUnTrPfCat` —— 上の 2 本(`Proposition 5.5, (ii)` の圏同値を含む)
2. `pfCatToRoot` —— 根 1 の部分を原文の `𝒞^pf` へ
3. `untrPf_modelFrobenioid` —— `(𝒞^un-tr)^pf ≌ model Frobenioid`
4. `pfRootToRlfFunctor_full` —— 係数を `ℚ` から `ℝ` へ -/
noncomputable def pfToRlfRight
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D)
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A))) :
    PfCat (istrPre P Fc) F₁ ⥤
      Crlf (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
        (birat_frobNormalized_of_unitTrivial (unTr_frobenioid P Fc G)
          (unTr_isotropic P Fc) (unTr_unitTrivial P Fc))
        hcharInj' hint' hfsmD :=
  pfToUnTrPfCat P Fc F₁ (unTr_frobenioidCore P Fc) hisoPf
    ⋙ pfCatToRoot (unTrPre P Fc) (unTr_frobenioidCore P Fc)
    ⋙ (untrPf_modelFrobenioid P Fc G hint hfsmD).functor
    ⋙ pfRootToRlfFunctor_full (unTr_frobIsotropicType Fc) (unTr_isotropic P Fc)
        (unTr_frobNormalizedType Fc G hint hfsmD) (unTr_unitTrivial P Fc)
        (unTr_frobenioid P Fc G) (untrPf_frobenioid P Fc G) hfsmD hcharInj' hint'

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の 1-可換図式は可換**(射の側、`𝒞^rlf` まで)。

```
𝒞  ──→ 𝒞^istr ──→ 𝒞^pf
│                    │
↓                    ↓
𝒞^un-tr ─→ (𝒞^un-tr)^pf ─→ 𝒞^rlf
```
-/
theorem cToPf_comp_pfToRlfRight_map (hiso : ∀ X : C, IsIsotropic P X)
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D)
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A)))
    {X Y : C} (φ : X ⟶ Y) :
    (cToPf Fc hiso F₁ ⋙ pfToRlfRight P Fc F₁ G hisoPf hint hfsmD hcharInj' hint').map φ
      = (cToUnTr hiso ⋙ toPfCat (unTrPre P Fc) (unTr_frobenioidCore P Fc)
          ⋙ (pfCatToRoot (unTrPre P Fc) (unTr_frobenioidCore P Fc)
            ⋙ (untrPf_modelFrobenioid P Fc G hint hfsmD).functor
            ⋙ pfRootToRlfFunctor_full (unTr_frobIsotropicType Fc) (unTr_isotropic P Fc)
              (unTr_frobNormalizedType Fc G hint hfsmD) (unTr_unitTrivial P Fc)
              (unTr_frobenioid P Fc G) (untrPf_frobenioid P Fc G) hfsmD
              hcharInj' hint')).map φ :=
  cToPf_comp_pfToUnTrPfCat_map_comp P Fc F₁ (unTr_frobenioidCore P Fc) _ hiso hisoPf φ

end ToRlf

/-! ## ★3. 図式の最初の矢印を **isotropification** に取り替える(条を落とす)

★★`Prop53Diag.lean` の `cToPf` / `cToUnTr` は `𝒞` が isotropic 型であること
(`hiso`)を使って `𝒞 ⟶ 𝒞^istr` を恒等的に作っていた。
★原文の第 1 の矢印は **`Proposition 1.9, (v)` の isotropification 関手**であり、
在庫(`isotropification`)にあるので、そちらに取り替えれば `hiso` は要らない。

★★★四角形の 1-可換性の中身(`homPfToUnTrPf_toHomPf`)は
**`Istr P` の任意の射について**述べてあるので、取り替えは無料である。 -/

section Isotropification

variable (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X)

/-- ★★**図式の上の行**(条なし)`𝒞 ⟶ 𝒞^istr ⟶ 𝒞^pf`。 -/
noncomputable def cToPfGen : C ⥤ PfCat (istrPre P Fc) F₁ :=
  isotropification P Fc ⋙ toPfCat (istrPre P Fc) F₁

/-- ★★**図式の左の縦の矢印**(条なし)`𝒞 ⟶ 𝒞^un-tr`。 -/
noncomputable def cToUnTrGen : C ⥤ UnTr P :=
  isotropification P Fc ⋙ istrToUnTr P

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の図式の四角形は 1-可換**(★条なし)。

```
𝒞  ──→ 𝒞^istr ──→ 𝒞^pf
│                    │
↓                    ↓
𝒞^un-tr ────────→ (𝒞^un-tr)^pf
```
-/
theorem cToPfGen_comp_pfToUnTrPfCat_map {X Y : C} (φ : X ⟶ Y) :
    (cToPfGen P Fc F₁ ⋙ pfToUnTrPfCat P Fc F₁ F₂ hisoPf).map φ
      = (cToUnTrGen P Fc ⋙ toPfCat (unTrPre P Fc) F₂).map φ :=
  homPfToUnTrPf_toHomPf (Fc := Fc) (F₁ := F₁) (F₂ := F₂)
    ((isotropification P Fc).map φ)

/-- ★★★★★★★**同上、`(𝒞^un-tr)^pf` から先の関手 `T` が何であっても**。

★★`T := pfCatToRoot ⋙ untrPf_modelFrobenioid ⋙ pfRootToRlfFunctor_full` と取れば
**原文の図式そのもの**(尻尾は `𝒞^rlf`)になる。 -/
theorem cToPfGen_comp_pfToUnTrPfCat_map_comp {E : Type ue} [Category.{ve} E]
    (T : PfCat (unTrPre P Fc) F₂ ⥤ E) {X Y : C} (φ : X ⟶ Y) :
    (cToPfGen P Fc F₁ ⋙ pfToUnTrPfCat P Fc F₁ F₂ hisoPf ⋙ T).map φ
      = (cToUnTrGen P Fc ⋙ toPfCat (unTrPre P Fc) F₂ ⋙ T).map φ :=
  congrArg T.map (cToPfGen_comp_pfToUnTrPfCat_map P Fc F₁ F₂ hisoPf φ)

end Isotropification

/-! ## ★4. `Proposition 5.3` 全体 -/

section Whole

variable [IsConnected D] (G : Frobenioid P)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の 1-可換図式**(★条なし)——

```
𝒞  ──→ 𝒞^istr ──→ 𝒞^pf
│                    │
↓                    ↓
𝒞^un-tr ─→ (𝒞^un-tr)^pf ─→ 𝒞^rlf
```

★上の道と下の道は一致する。残る仮引数は原文の常備仮定だけである:
`hint`(`Φ` が整域単系)、`hfsmD`(`𝒟` が FSM 型)、
`hcharInj'` / `hint'`(`ℝ≥0` 係数、`Φ` の perf-factorial 性から出るもの)。 -/
theorem prop_5_3_diagram
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D)
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A)))
    {X Y : C} (φ : X ⟶ Y) :
    (cToPfGen P Fc F₁ ⋙ pfToRlfRight P Fc F₁ G (pfIstr_isotropic P Fc F₁)
        hint hfsmD hcharInj' hint').map φ
      = (cToUnTrGen P Fc ⋙ toPfCat (unTrPre P Fc) (unTr_frobenioidCore P Fc)
          ⋙ (pfCatToRoot (unTrPre P Fc) (unTr_frobenioidCore P Fc)
            ⋙ (untrPf_modelFrobenioid P Fc G hint hfsmD).functor
            ⋙ pfRootToRlfFunctor_full (unTr_frobIsotropicType Fc) (unTr_isotropic P Fc)
              (unTr_frobNormalizedType Fc G hint hfsmD) (unTr_unitTrivial P Fc)
              (unTr_frobenioid P Fc G) (untrPf_frobenioid P Fc G) hfsmD
              hcharInj' hint')).map φ :=
  cToPfGen_comp_pfToUnTrPfCat_map_comp P Fc F₁ (unTr_frobenioidCore P Fc)
    (pfIstr_isotropic P Fc F₁) _ φ

end Whole

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

/-- ★★★★★★★**[FrdI] Proposition 5.3(命題全体)** の locator。

原文 (FrdI p.103):
> Proposition 5.3. (Realifications of Frobenioids) Suppose that Φ is perf-factorial.

★★**3 つの文の対応**:

| 原文 | 実装 |
|---|---|
| 第 1 文 `𝒞^rlf` の定義(`Φ^rlf` と `ℝ·Φ^birat`) | `Crlf` / `scModel NNReal`(`Prop53Rlf.lean`) |
| 第 2 文 `𝒞^un-tr` は `(Φ, Φ^birat)` の model | `unTr_isOfModelType` / `unTr_modelFrobenioid` / `unTr_ratFnData_bmon_val_of_C` |
| 第 2 文 `(𝒞^un-tr)^pf` は `(Φ^pf, ℚ·Φ^birat)` の model | `untrPf_isOfModelType` / `untrPf_modelFrobenioid` / `untrPf_ratFnData_bmon_val_of_C` |
| `ℚ·Φ^birat = Φ^birat ⊗_ℤ ℚ = (Φ^birat)^pf` | `qPhiBiratOn_map_eq` / `phiBiratOn_pf_eq_qPhiBiratOn_full` |
| 第 3 文 図式の上の行・左の縦 | `cToPfGen` / `cToUnTrGen`(isotropification 経由) |
| 第 3 文 図式の下の行 | `untrToSc` / `scBaseFunctor` / `untrToSc_comp_scBaseFunctor` |
| 第 3 文 図式の右の縦 | `pfToUnTrPfCat` / `pfToRlfRight` |
| 第 3 文 図式の 1-可換性 | `prop_5_3_diagram` |

★★**残る仮引数は原文の常備仮定だけ**である:
`hint`(`Φ` が整域単系)、`hfsmD`(`𝒟` が FSM 型)、
`hcharInj'` / `hint'`(`ℝ≥0` 係数、`Φ` の perf-factorial 性から出るもの)。 -/
def prop_5_3.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
