import ABC3.Found.Arakelov.PicTensorSection

/-!
# Arakelov (B1) 第 83 ブロック —— **前層射 `tilde M ⊗ tilde N ⟶ tilde (M ⊗ N)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★切断の写像を前層射に組む

第 82 ブロックの `tensorSectionMap` を全ての `U` について集めると
**前層加群の射**になる。

## ★★詰まった 1 点 —— 加群インスタンスの二経路

`ModuleCat.ofHom` は

    Module 𝒪(U) ((structurePresheafInModuleCat R (M⊗N) ⋙ forget₂ _ Ab).obj U)

を要求するが、**instance 検索が拾わない**([[ring-instance-two-paths]])。
★これは mathlib の `moduleStructurePresheaf` の定義の中で `letI` として
与えられているもので、外からは見えない。

★★**回避法**: `letI` で**エラーが要求した型そのまま**に書き、
`inferInstanceAs` で `structureSheafInType` 側の instance を渡す。
★★★`(tilde ..).val.obj U` と書いても**通らない**——展開形で書く必要がある(実測)。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tildeTensorApp` | ★各開集合での写像 |
| `tildeTensorPre` | ★★★★**前層射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace TensorProduct
open StructureSheaf

variable (R : CommRingCat.{u}) (M N : ModuleCat.{u} (R : Type u))

/-- ★**各開集合での写像**。 -/
noncomputable def tildeTensorApp (U : (Opens ↥(Spec R))ᵒᵖ) :
    ((tilde M).val ⊗ (tilde N).val).obj U
      ⟶ (tilde (ModuleCat.of (R : Type u) (M ⊗[(R : Type u)] N))).val.obj U :=
  letI : Module (((Spec R).presheaf ⋙ forget₂ CommRingCat RingCat).obj U : Type u)
      ((structurePresheafInModuleCat (R : Type u)
        (ModuleCat.of (R : Type u) (M ⊗[(R : Type u)] N))
        ⋙ forget₂ (ModuleCat (R : Type u)) Ab).obj U : Type u) :=
    inferInstanceAs (Module ((structureSheafInType (R : Type u) (R : Type u)).1.obj U : Type u)
      ((structureSheafInType (R : Type u) (M ⊗[(R : Type u)] N)).1.obj U : Type u))
  ModuleCat.ofHom (tensorSectionMap (R : Type u) M N U.unop)

/-- ★**`letI` を剥がす `rfl` 補題**——これが無いと `map_add` が発火しない
(第 25 ブロックで効いたのと同じ手)。 -/
@[simp] theorem tildeTensorApp_apply (U : (Opens ↥(Spec R))ᵒᵖ)
    (z : ((tilde M).val ⊗ (tilde N).val).obj U) :
    (tildeTensorApp R M N U).hom z = tensorSectionMap (R : Type u) M N U.unop z := rfl

/-- ★★★★**前層射** `tilde M ⊗ tilde N ⟶ tilde (M ⊗_R N)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★純テンソルでは `rfl`、加法性は上の `rfl` 補題で `map_add` を発火させる。 -/
noncomputable def tildeTensorPre :
    (tilde M).val ⊗ (tilde N).val
      ⟶ (tilde (ModuleCat.of (R : Type u) (M ⊗[(R : Type u)] N))).val where
  app := tildeTensorApp R M N
  naturality {X Y} f := by
    ext s
    induction s using TensorProduct.induction_on with
    | zero => exact (map_zero _).trans (map_zero _).symm
    | tmul a b => rfl
    | add x y hx hy => exact ((map_add _ _ _).trans (by rw [hx, hy])).trans (map_add _ _ _).symm

/-! ## ★出典の紐付け(`.src`) -/

def tildeTensorApp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断のテンソル写像の各開集合成分)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
