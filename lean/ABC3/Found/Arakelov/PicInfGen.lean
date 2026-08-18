import ABC3.Found.Arakelov.PicFreeNat

/-!
# Arakelov (B1) 第 37 ブロック —— **交わりの同型の生成元での値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★Y 側と X 側の両方に効く

最終段の両辺には、同じ形の計算が 1 回ずつ現れる:

| 側 | 何を計算するか |
|---|---|
| Y | `freeYonedaEquiv (iotaY V W)`(第 26 の同型の `inv`) |
| X | `freeYonedaEquiv (targetIso.inv)`(同じもの、`f⁻¹V, f⁻¹W` で) |

★★**汎用の束の圏で 1 回書けば両方に効く**——本ブロックがそれである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `free_map_app_freeMk` | ★`free` の `map` は生成元を生成元に送る |
| `freeYonedaEquiv_freeYonedaInfIso_inv` | ★★★★**交わりの同型の生成元での値** |

## ★★★実装の要点(2026-08-18 実測)

★`(freeTensorIso ..).hom` を `freeTensorHom` に**先に書き換える**こと。
そうしないと `erw` が `whnf` で**タイムアウト**する(第 25 ブロックに
`freeTensorIso_hom` を追加した理由)。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite Limits

/-! ## ★`free` の `map` は生成元を生成元に送る -/

/-- ★**`free` の `map` は生成元を生成元に送る**。 -/
theorem free_map_app_freeMk {C : Type u} [SmallCategory C] {R₀ : Cᵒᵖ ⥤ RingCat.{u}}
    {F G : Cᵒᵖ ⥤ Type u} (φ : F ⟶ G) (Z : Cᵒᵖ) (x : F.obj Z) :
    ((PresheafOfModules.free R₀).map φ).app Z (ModuleCat.freeMk x)
      = ModuleCat.freeMk (φ.app Z x) :=
  ModuleCat.free_map_apply _ _

/-! ## ★★★★交わりの同型の生成元での値 -/

variable {α : Type u} [Lattice α] {R : αᵒᵖ ⥤ CommRingCat.{u}}

/-- ★★★★★**交わりの同型の生成元での値**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 26 ブロックの `freeYonedaInfIso` の `inv` を生成元に当てると、
**2 つの包含射の生成元の純テンソル**になる。

★★これは最終段の**両辺**に現れる計算である
(Y 側は `V, W`、X 側は `f⁻¹V, f⁻¹W`)。 -/
theorem freeYonedaEquiv_freeYonedaInfIso_inv (V W : α) :
    PresheafOfModules.freeYonedaEquiv ((freeYonedaInfIso (R := R) V W).inv)
      = (ModuleCat.freeMk (homOfLE (inf_le_left : V ⊓ W ≤ V))
          : (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
              (yoneda.obj V)).obj (op (V ⊓ W)))
        ⊗ₜ (ModuleCat.freeMk (homOfLE (inf_le_right : V ⊓ W ≤ W))
          : (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
              (yoneda.obj W)).obj (op (V ⊓ W))) := by
  rw [freeYonedaEquiv_apply]
  show ((PresheafOfModules.free _).map (yonedaInfIso V W).inv
      ≫ (freeTensorIso (R := R) (yoneda.obj V) (yoneda.obj W)).hom).app _ _ = _
  erw [CategoryTheory.comp_apply]
  rw [freeTensorIso_hom, free_map_app_freeMk]
  erw [freeTensorHom_app_freeMk]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def freeYonedaEquiv_freeYonedaInfIso_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——交わりの同型の生成元での値)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
