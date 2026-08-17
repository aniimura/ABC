import ABC3.Found.Arakelov.PicPullTensorFree

/-!
# Arakelov (B1) 第 28 ブロック —— **余極限で同型を持ち上げる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★これが「生成元 → 全体」の器具である

第 27 ブロックで**生成元の上での同型**が出た。
★★それを全対象へ持ち上げる器具が本ブロックである:

> **両側が余極限を保ち、図式の上で同型なら、余極限の点でも同型。**

★★★mathlib に該当する補題は**無かった**(2026-08-18 実測)ので書いた。
機構は `IsColimit.coconePointsIsoOfNatIso` と自然性だけである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isIso_app_of_isColimit` | ★★★★★**余極限で同型を持ち上げる**(汎用) |
| `isIso_app_of_isColimit_of_hasColimit` | ★同上、`colimit` の形で |

## ★★★使い方

前層加群 `M` は **`free (yoneda X)` の余積の余核**である
(mathlib `isColimitFreeYonedaCoproductsCokernelCofork`)。
★したがって本補題を **2 回**当てればよい——余積で 1 回、余核で 1 回。
-/

universe u v u' v' u'' v''

namespace ABC3.Found.Arakelov

open CategoryTheory Limits

/-! ## ★★★★★余極限で同型を持ち上げる -/

/-- ★★★★★**余極限で同型を持ち上げる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★両側の関手が図式 `K` の余極限を保ち、`τ` が `K` の各対象で同型ならば、
`τ` は余極限の点でも同型である。

★★★証明は `IsColimit.coconePointsIsoOfNatIso` と `τ` の自然性だけである——
`τ.app c.pt` が「同型として作った射」と一致することを `hom_ext` で見る。 -/
theorem isIso_app_of_isColimit {C : Type u} [Category.{v} C] {D : Type u'} [Category.{v'} D]
    {J : Type u''} [Category.{v''} J]
    {A B : C ⥤ D} (τ : A ⟶ B) (K : J ⥤ C) (c : Cocone K) (hc : IsColimit c)
    [PreservesColimit K A] [PreservesColimit K B]
    (h : ∀ j, IsIso (τ.app (K.obj j))) : IsIso (τ.app c.pt) := by
  haveI : ∀ j, IsIso ((Functor.whiskerLeft K τ).app j) := h
  haveI : IsIso (Functor.whiskerLeft K τ) := NatIso.isIso_of_isIso_app _
  have key : τ.app c.pt
      = (IsColimit.coconePointsIsoOfNatIso (isColimitOfPreserves A hc)
          (isColimitOfPreserves B hc) (asIso (Functor.whiskerLeft K τ))).hom := by
    refine (isColimitOfPreserves A hc).hom_ext ?_
    intro j
    erw [IsColimit.comp_coconePointsIsoOfNatIso_hom]
    exact τ.naturality (c.ι.app j)
  rw [key]
  infer_instance

/-- ★同じことを `colimit` の形で述べたもの。 -/
theorem isIso_app_of_colimit {C : Type u} [Category.{v} C] {D : Type u'} [Category.{v'} D]
    {J : Type u''} [Category.{v''} J]
    {A B : C ⥤ D} (τ : A ⟶ B) (K : J ⥤ C) [HasColimit K]
    [PreservesColimit K A] [PreservesColimit K B]
    (h : ∀ j, IsIso (τ.app (K.obj j))) : IsIso (τ.app (colimit K)) :=
  isIso_app_of_isColimit τ K _ (colimit.isColimit K) h

/-! ## ★出典の紐付け(`.src`) -/

def isIso_app_of_isColimit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——余極限で同型を持ち上げる器具)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
