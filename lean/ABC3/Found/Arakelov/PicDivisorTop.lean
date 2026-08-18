import ABC3.Found.Arakelov.PicDivisor

/-!
# Arakelov (B2) 第 187 ブロック —— **空因子は単位元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★`ofDivisor_top` が出た

    𝒪_X(⊤) = 1  すなわち  idealSheaf ⊤ ≅ 𝒪_X

★機構は「`idealSections ⊤ U = ⊤`」の 1 行と、包含射 `idealSheaf D ⟶ 𝒪_X` が
各点で全単射であることだけである。

## ★★包含射も `homMk` で組む

第 178 と同じく、構造体を直に組むと `Module` の綴りが 2 通りになる。
★`PresheafOfModules.homMk`(アーベル群の射 + 線型性)なら
**部分加群の包含 `s ↦ ↑s` を書くだけ**で済む——加法も線型性も `rfl`。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `idealIncl` | ★イデアル層から構造層への包含 |
| `idealSections_top` | ★`⊤` の切断は全体 |
| `idealTopIso` | ★★★`idealPresheaf ⊤ ≅ 𝟙_` |
| `idealSheafTopIso` | ★★★層加群としての同型 |
| `ofDivisorSheaf_top` | ★★★★**空因子は単位元** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★**イデアル層から構造層への包含**。 -/
noncomputable def idealIncl (D : X.IdealSheafData) :
    idealPresheaf D ⟶ 𝟙_ (X.PresheafOfModules) :=
  PresheafOfModules.homMk
    { app := fun U => AddCommGrpCat.ofHom
        (AddMonoidHom.mk' (fun s : (idealSections D U.unop) => (s : (Γ(X, U.unop) : Type u)))
          (fun _ _ => rfl))
      naturality := by intro U V f; ext s; rfl }
    (fun U c s => rfl)

/-- ★**`⊤` の切断は全体である**。 -/
theorem idealSections_top (U : X.Opens) :
    idealSections (⊤ : X.IdealSheafData) U = ⊤ := by
  ext s
  simp [idealSections]

theorem bijective_idealIncl_top (U : (X.Opens)ᵒᵖ) :
    Function.Bijective (((idealIncl (⊤ : X.IdealSheafData)).app U).hom) := by
  refine ⟨fun a b h => Subtype.ext h, fun y => ?_⟩
  exact ⟨⟨y, by rw [idealSections_top]; trivial⟩, rfl⟩

/-- ★★★**`idealPresheaf ⊤ ≅ 𝟙_`**。 -/
noncomputable def idealTopIso :
    idealPresheaf (⊤ : X.IdealSheafData) ≅ 𝟙_ (X.PresheafOfModules) :=
  haveI happ : ∀ U, IsIso ((idealIncl (⊤ : X.IdealSheafData)).app U) := fun U =>
    (ConcreteCategory.isIso_iff_bijective _).2 (bijective_idealIncl_top U)
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map (idealIncl (⊤ : X.IdealSheafData))) := by
    rw [NatTrans.isIso_iff_isIso_app]
    intro U
    haveI := happ U
    exact inferInstanceAs
      (IsIso ((forget₂ _ AddCommGrpCat).map ((idealIncl (⊤ : X.IdealSheafData)).app U)))
  haveI : IsIso (idealIncl (⊤ : X.IdealSheafData)) :=
    isIso_of_reflects_iso _ (PresheafOfModules.toPresheaf _)
  asIso (idealIncl (⊤ : X.IdealSheafData))

/-- ★★★**層加群としての同型** `idealSheaf ⊤ ≅ 𝒪_X`。 -/
noncomputable def idealSheafTopIso :
    idealSheaf (⊤ : X.IdealSheafData) ≅ unitModules X where
  hom := ⟨idealTopIso.hom⟩
  inv := ⟨idealTopIso.inv⟩
  hom_inv_id := SheafOfModules.hom_ext idealTopIso.hom_inv_id
  inv_hom_id := SheafOfModules.hom_ext idealTopIso.inv_hom_id

/-- ★★★★**空因子は単位元である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `CartierPicData.ofDivisor_top` である。 -/
theorem ofDivisorSheaf_top (X : Scheme.{u}) :
    ofDivisorSheaf (⊤ : X.IdealSheafData) = 1 := by
  classical
  refine PicSheaf.mk_eq_mk ?_
  rw [divisorInvSheaf, dif_pos (isCartier_top X)]
  exact idealSheafTopIso

/-! ## ★出典の紐付け(`.src`) -/

def ofDivisorSheaf_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——空因子は単位元)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
