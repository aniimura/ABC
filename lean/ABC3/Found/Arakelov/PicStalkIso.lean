import ABC3.Found.Arakelov.PicAwayFree

/-!
# Arakelov (B1) 第 77 ブロック —— **茎で同型なら同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★層化を跨ぐ器具

§9-76 で測った通り、**層化は茎を変えない**(mathlib
`TopCat.Presheaf.stalkFunctor_map_unit_toSheafify_isIso`)。
★したがって `tensorModules`(= 前層テンソルの層化)の茎は
**前層テンソルの茎**そのものである。

★★これを使うには「層加群の射が全ての茎で同型なら同型」という器具が要る。
mathlib は `TopCat.Sheaf` の版を持つ(`app_isIso_of_stalkFunctor_map_iso`)ので、
★★★`SheafOfModules.toSheaf` を通して持ち込み、
`toPresheaf` と `forget` の 2 段で反射させる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isIso_of_stalkIso` | ★★★★**茎で同型なら層加群の射は同型** |

## ★★★これで比較射の証明が茎の計算に落ちる

    比較射が同型  ⟸  ∀ p, (M ⊗_R N)_p ≅ M_p ⊗_{R_p} N_p
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X : Scheme.{u}} {A B : X.Modules}

/-- ★★★★**層加群の射が全ての茎で同型なら、その射は同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「層化を跨ぐ」器具である——層化は茎を変えないので、
層化を含む対象どうしの比較も茎だけで決まる。 -/
theorem isIso_of_stalkIso (φ : A ⟶ B)
    (h : ∀ x : X, IsIso ((TopCat.Presheaf.stalkFunctor AddCommGrpCat.{u} x).map
      ((SheafOfModules.toSheaf X.ringCatSheaf).map φ).1)) :
    IsIso φ := by
  haveI := h
  haveI happ : ∀ U : (Opens X)ᵒᵖ,
      IsIso (((SheafOfModules.toSheaf X.ringCatSheaf).map φ).1.app U) := by
    intro U
    exact TopCat.Presheaf.app_isIso_of_stalkFunctor_map_iso _ U.unop
  haveI : IsIso ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map
      ((SheafOfModules.forget X.ringCatSheaf).map φ)) := by
    rw [NatTrans.isIso_iff_isIso_app]
    intro U
    exact happ U
  haveI : IsIso ((SheafOfModules.forget X.ringCatSheaf).map φ) :=
    isIso_of_reflects_iso _ (PresheafOfModules.toPresheaf X.ringCatSheaf.obj)
  exact (SheafOfModules.fullyFaithfulForget X.ringCatSheaf).isIso_of_isIso_map φ

/-! ## ★出典の紐付け(`.src`) -/

def isIso_of_stalkIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——茎で同型なら同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
