import ABC3.Found.Arakelov.PicUnitRing
import ABC3.Found.Arakelov.PicTrivialNoSheaf

/-!
# Arakelov (B2) 第 206 ブロック —— **単位の自己射は全射なら単元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-225 の補題 A の核

    φ : 𝟙_ ⟶ 𝟙_ が終対象で全射  ⟹  unitEndEquiv φ は単元

★★これがあると「局所自明どうしの局所全射は同型」が出る——共通の自明化開の上では
射は**係数の掛け算**(第 166 の `unitEndEquiv`)であり、全射なら係数は単元だからである。

## ★★機構は 3 行

`φ(y) = 1` なる `y` を取ると

    1 = unitVal (φ y) = unitVal (y • 1 の係数) * unitEndEquiv φ

なので `unitEndEquiv φ` は単元である。★`y` を**係数**として読み替える橋
`unitScal`(型付き恒等写像)を置くと、`y • 1 = y` が第 164 の `smul_unitOne` で出る。

★★★[[typed-identity-bridge]] の 4 例目である——単位対象の切断は
「値」としても「係数」としても読めるので、**両方の恒等写像**を用意しておくと
`rw` が素直に通る。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `unitScal` | ★切断を係数と見る橋 |
| `unitScal_smul_unitOne` / `unitVal_unitScal` | ★橋の 2 法則 |
| `isUnit_unitEndEquiv_of_surjective` | ★★★★**全射なら単元** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (U : X.Opens)

/-- ★単位対象の切断を係数環の元と見る(型付き恒等写像)。 -/
noncomputable def unitScal
    (x : ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u)) :
    ((((Over.forget U).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj
      (op (Over.mk (𝟙 U))) : Type u) := x

theorem unitScal_smul_unitOne
    (x : ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u)) :
    unitScal U x • unitOne U = x := smul_unitOne U (unitScal U x)

theorem unitVal_unitScal
    (x : ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u)) :
    unitVal U x = (unitScal U x : (Γ(X, U) : Type u)) := rfl

/-- ★★単位の自己準同型が終対象で全射なら、対応する係数は単元である。 -/
theorem isUnit_unitEndEquiv_of_surjective (φ : End (𝟙_ (PresheafModulesOn X U)))
    (hsurj : Function.Surjective ((φ.app (op (Over.mk (𝟙 U)))).hom)) :
    IsUnit (unitEndEquiv U φ) := by
  obtain ⟨y, hy⟩ := hsurj (unitOne U)
  refine IsUnit.of_mul_eq_one (unitVal U y) ?_
  have h2 : (φ.app (op (Over.mk (𝟙 U)))).hom y
      = unitScal U y • (φ.app (op (Over.mk (𝟙 U)))).hom (unitOne U) := by
    rw [← map_smul, unitScal_smul_unitOne]
  have h4 : unitVal U ((φ.app (op (Over.mk (𝟙 U)))).hom y)
      = unitVal U y * unitEndEquiv U φ := by
    rw [h2, unitVal_smul, unitVal_unitScal]
    rfl
  have h5 : unitVal U ((φ.app (op (Over.mk (𝟙 U)))).hom y) = 1 := by
    rw [hy]; exact unitVal_one U
  rw [h5] at h4
  exact (mul_comm _ _).trans h4.symm


/-! ## ★出典の紐付け(`.src`) -/

def isUnit_unitEndEquiv_of_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——単位の自己射は全射なら単元)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
