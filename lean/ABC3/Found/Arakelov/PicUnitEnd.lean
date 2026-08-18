import ABC3.Found.Arakelov.PicUnitMul

/-!
# Arakelov (B2) 第 164 ブロック —— ★★★★★★**単位対象の自己射環**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★双対の加群構造の鍵

mathlib の `Preadditive.moduleEndRight : Module (End Y) (X ⟶ Y)` により、
`Hom(F|_U, 𝟙_)` は `End(𝟙_)` 加群である。
★あとは環準同型 `End(𝟙_) →+* Γ(X,U)` があればよい。

## ★★★★逃げ道——**型付き恒等写像を橋にする**

`(𝟙_).obj t` の台は `Γ(X,U)` だが、`One` も `Mul` も `SMul` も
**綴りが違うと見つからない**。★5 通り試して全滅した後、

    def unitVal (x : ((𝟙_ …).obj t : Type u)) : (Γ(X,U) : Type u) := x

という**型付き恒等写像**を置いたら、

    unitVal (a • b) = unitVal a * unitVal b     ★★**`rfl`**
    unitVal (unitOne) = 1                       ★★**`rfl`**

が両方 `rfl` で通った。★★★これで作用と掛け算が**同じものだと言える**。

## ★★★試して全滅した手(記録)

| 手 | 結果 |
|---|---|
| `LinearMap.mul` / `LinearMap.lsmul` | `HMul` / `HSMul` が無い |
| **型注釈**で綴りを合わせる | ★注釈は推論型を**変えない**(本 session 3 例目) |
| `inferInstanceAs` で `CommRing` を足す | ★★**`SMul` の経路が 2 本になり `rw` が当たらなくなった** |
| RingCat 綴りに統一 | `One` が見つからない |
| **型付き恒等写像**(`unitVal`) | ★★★**通った** |

★★★★教訓: **instance を足すより、型付き恒等関数で橋を架ける方が安全**
——instance は既存の経路と競合しうるが、恒等関数は競合しない。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (U : X.Opens)

/-- ★単位対象の切断としての `1`。 -/
noncomputable def unitOne :
    ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u) :=
  (1 : (Γ(X, U) : Type u))

/-- ★単位対象の切断を `Γ(X,U)` の元と見る。 -/
noncomputable def unitVal
    (x : ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u)) :
    (Γ(X, U) : Type u) := x

/-- ★★作用は掛け算(`rfl` か?)。 -/
theorem unitVal_smul
    (a : ((((Over.forget U).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj
      (op (Over.mk (𝟙 U))) : Type u))
    (b : ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u)) :
    unitVal U (a • b) = unitVal U a * unitVal U b := rfl

/-- ★★`unitVal` は `1` を `1` に送る。 -/
theorem unitVal_one : unitVal U (unitOne U) = 1 := rfl

/-- ★`unitVal` は単射(恒等写像)。 -/
theorem unitVal_injective : Function.Injective (unitVal U) := fun _ _ h => h

/-- ★★`a • 1 = a`。 -/
theorem smul_unitOne
    (a : ((((Over.forget U).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj
      (op (Over.mk (𝟙 U))) : Type u)) :
    a • (unitOne U) = a :=
  unitVal_injective U (by rw [unitVal_smul, unitVal_one, mul_one])

/-- ★★★単位対象の自己射環から `Γ(X,U)` への評価。 -/
noncomputable def unitEnd : End (𝟙_ (PresheafModulesOn X U)) →+* (Γ(X, U) : Type u) where
  toFun φ := unitVal U ((φ.app (op (Over.mk (𝟙 U)))).hom (unitOne U))
  map_one' := rfl
  map_zero' := rfl
  map_add' := fun _ _ => rfl
  map_mul' := fun f g => by
    have hsm := (f.app (op (Over.mk (𝟙 U)))).hom.map_smul
      ((g.app (op (Over.mk (𝟙 U)))).hom (unitOne U)) (unitOne U)
    rw [smul_unitOne] at hsm
    show unitVal U ((f.app (op (Over.mk (𝟙 U)))).hom
      ((g.app (op (Over.mk (𝟙 U)))).hom (unitOne U))) = _
    rw [hsm, unitVal_smul, mul_comm]


/-! ## ★出典の紐付け(`.src`) -/

def unitEnd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——単位対象の自己射環)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
