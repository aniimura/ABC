import ABC3.Found.Arakelov.PicUnitRes

/-!
# Arakelov (B2) 第 169 ブロック —— ★★★★★★★★**双対の前層**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★mathlib に無い「層加群の双対」の第 1 歩

    dualPresheaf F : X.PresheafOfModules
      obj U = ((restrictPresheafFunctor X U).obj F ⟶ 𝟙_)

★`Mathlib/Algebra/Category/ModuleCat/Sheaf/` に `Monoidal.lean` は無く、
内部 Hom も双対も**存在しない**(2026-08-18 実測)。

## ★★組み上げに要った 4 つの部品

| 部品 | ブロック |
|---|---|
| `Hom` は `End(𝟙_)` 加群 | mathlib `Preadditive.moduleEndRight` |
| `End(𝟙_) ≃+* Γ(X,U)` | ★第 164–166 |
| `c • φ = φ ≫ unitMul c` | ★第 167 |
| 制限の半線型性 | ★第 168 |

## ★★★逃げ道 4 つ(すべて綴りの問題)

| 症状 | 逃げ道 |
|---|---|
| `Module ((…).obj V) (双対 at V)` が見つからない | ★**索引を反対圏にした instance をもう 1 本**(`dualModuleOp`) |
| `Functor.map_add` の `Additive` が無い | ★**`rfl`** で通る |
| `show` で `(restrictOnFunctor _).map φ` の型が合わない | ★`have h1` で先に書き換える |
| 最後の `•` が `Module.compHom` 経由で `rw` が当たらない | ★★**`exact (dual_smul_eq …).symm`**(defeq で通る) |

★★★★本 session で**5 度目**の綴り問題であり、毎回
「下流の要求に合わせる」「項で手渡す」で解決している。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules)

/-- ★★双対の前層。 -/
noncomputable def dualPresheaf : X.PresheafOfModules where
  obj U := ModuleCat.of ((X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u)
    ((restrictPresheafFunctor X U.unop).obj F ⟶ 𝟙_ (PresheafModulesOn X U.unop))
  map {U V} f :=
    letI : Module (((X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U) : Type u)
        ((restrictPresheafFunctor X V.unop).obj F ⟶ 𝟙_ (PresheafModulesOn X V.unop)) :=
      Module.compHom _ ((X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}).map f).hom
    ModuleCat.ofHom
      { toFun := fun φ => (restrictOnFunctor (leOfHom f.unop)).map φ
        map_add' := fun a b => rfl
        map_smul' := fun c φ => by
          have h1 : (restrictOnFunctor (leOfHom f.unop)).map (c • φ)
              = (restrictOnFunctor (leOfHom f.unop)).map (φ ≫ unitMul U.unop c) := by
            rw [dual_smul_eq]
          show (restrictOnFunctor (leOfHom f.unop)).map (c • φ) = _
          rw [h1, Functor.map_comp, unitMul_res]
          exact (dual_smul_eq F V.unop _ _).symm }


/-! ## ★出典の紐付け(`.src`) -/

def dualPresheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——双対の前層)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
