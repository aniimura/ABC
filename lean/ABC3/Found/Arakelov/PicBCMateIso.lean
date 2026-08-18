import ABC3.Found.Arakelov.PicResGenVal

/-!
# Arakelov (B1) 第 53 ブロック —— ★★★★★★★★**Beck–Chevalley の mate は生成元で同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★取れたもの

    (f|)^*_pre ((free (yoneda W))|_V)  ⟶  (f^*_pre (free (yoneda W)))|_{f⁻¹V}

が**同型の合成に等しい**——したがって**生成元の上で mate は同型**である。

★★★第 29 ブロックの器具(生成元で同型なら全対象で同型)を当てれば、
**Beck–Chevalley が全対象で同型**になる。

## ★★★★★核心は `δ` のときと同じ

> **両辺とも `freeYonedaEquiv (unit)` の生成元での値に落ちる。**

★左辺は第 38 ブロック(`unit` の値 = 同型の `inv` の値)、
★★右辺は第 51 ブロック(制限した site でも `unit` は生成元を生成元に送る)。

## ★★実装の要点(2026-08-18 実測)

★★★**順序が命である**:

1. `hpm` / `comp_apply` で**先に剥がす**
2. **その後で** `isoHomUnitGenOn`(第 51)を当てる——2 段が 1 度に済む
3. 先に `freeYonedaEquivUnitGen` を当てると剥がす対象が消えて詰まる

★`erw` は補題の内部構造を展開してしまうので、
**`app` 版の補題**(第 52 に追加)を切っておく必要があった。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★**Beck–Chevalley の mate は生成元の上で同型の合成に等しい**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★核心は `δ`(第 40 ブロック)のときと同じ——
**両辺とも `freeYonedaEquiv (unit)` の生成元での値に落ちる**。 -/
theorem bcMate_free_eq (V W : Y.Opens) :
    (pullbackPreOn f V).map (restrictFreeYonedaIso V W).inv
        ≫ (bcMate f V).app (freeY W)
      = (pullbackOnFreeYonedaIso f V (objOn V W)).hom
        ≫ (restrictFreeYonedaIso ((Opens.map f.base).obj V)
            ((Opens.map f.base).obj W)).inv
        ≫ (restrictPresheafFunctor X ((Opens.map f.base).obj V)).map
            (pullbackFreeYonedaIso f W).inv := by
  refine ((adjOn f V).homEquiv _ _).injective ?_
  erw [Adjunction.homEquiv_naturality_left]
  show _ ≫ (adjOn f V).homEquiv _ _ (((adjOn f V).homEquiv _ _).symm _) = _
  rw [Equiv.apply_symm_apply]
  erw [Adjunction.homEquiv_naturality_right]
  rw [Adjunction.homEquiv_unit]
  refine PresheafOfModules.freeYonedaEquiv.injective ?_
  simp only [PresheafOfModules.freeYonedaEquiv_comp]
  have hL := freeYonedaEquiv_restrictFreeYonedaIso_inv V W
  erw [hL]
  have hpm : ∀ {A B : PresheafModulesOn X ((Opens.map f.base).obj V)} (g : A ⟶ B)
      (Z : Over V) (x : ((PresheafOfModules.pushforward (phiOn f V)).obj A).obj (op Z)),
      ((PresheafOfModules.pushforward (phiOn f V)).map g).app (op Z) x
        = g.app (op ((overPost f V).obj Z)) x := fun _ _ _ => rfl
  erw [hpm, CategoryTheory.comp_apply, hpm, CategoryTheory.comp_apply]
  erw [isoHomUnitGenOn f V (objOn V W)]
  have hres : ∀ {A B : Y.PresheafOfModules} (g : A ⟶ B) (Z : Over V)
      (x : ((restrictPresheafFunctor Y V).obj A).obj (op Z)),
      ((restrictPresheafFunctor Y V).map g).app (op Z) x = g.app (op Z.left) x :=
    fun _ _ _ => rfl
  erw [hres]
  have h38 := unit_app_eq_inv_app f W (homOfLE (inf_le_left : W ⊓ V ≤ W))
  erw [h38]
  have hresX : ∀ {A B : X.PresheafOfModules} (g : A ⟶ B)
      (Z : Over ((Opens.map f.base).obj V))
      (x : ((restrictPresheafFunctor X ((Opens.map f.base).obj V)).obj A).obj (op Z)),
      ((restrictPresheafFunctor X ((Opens.map f.base).obj V)).map g).app (op Z) x
        = g.app (op Z.left) x := fun _ _ _ => rfl
  erw [hresX]
  erw [Functor.mapIso_inv, free_map_app_freeMk]
  rfl


/-! ## ★出典の紐付け(`.src`) -/

def bcMate_free_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Beck–Chevalley の mate が生成元で同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
