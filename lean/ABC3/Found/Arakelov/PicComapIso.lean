import ABC3.Found.Arakelov.PicBaseChange
import Mathlib.Algebra.Category.Ring.Basic

/-!
# Arakelov (B2) 第 199 ブロック —— **同型に沿った引き戻しでも可逆**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★開埋め込みで `comap` を運ぶための橋

mathlib の `ideal_comap_of_isOpenImmersion` は

    (I.comap f).ideal U = (I.ideal (f ''ᵁ U)).comap ((f.appIso U).inv.hom)

を与える——**環の同型に沿った引き戻し**である。したがって

    可逆性が同型に沿った引き戻しで保たれる

を示せば、開埋め込みに沿って可逆性を運べる。★これで一般の `X` の
アフィン開へ `Spec` の議論(第 198)を持ち込める。

## ★★第 197 を**両向き**に使う

第 197 は「`algebraMap R S` が全射 ⟹ `Module.Invertible R M` から `S` へ」。
★同型なら**役割を入れ替えて**逆向きにも使える(`algebraMap S R = e.inv` も全射)。
そのあと `↥(I.comap e.hom) ≃ₗ[R] ↥I` を手で組む。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `invertible_comap_of_isIso` | ★★★★**同型に沿った引き戻しでも可逆** |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory TensorProduct

variable {R S : CommRingCat.{u}} (e : R ≅ S) (I : Ideal (S : Type u))

/-- ★★★★**環の同型に沿った引き戻しでも可逆である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで開埋め込みに沿って可逆性を運べる。 -/
theorem invertible_comap_of_isIso [Module.Invertible (S : Type u) I] :
    Module.Invertible (R : Type u) (I.comap e.hom.hom) := by
  letI : Algebra (R : Type u) (S : Type u) := e.hom.hom.toAlgebra
  letI : Algebra (S : Type u) (R : Type u) := e.inv.hom.toAlgebra
  have heh : ∀ r : (R : Type u), algebraMap (R : Type u) (S : Type u) r = e.hom.hom r := fun _ => rfl
  have hhe : ∀ s : (S : Type u), algebraMap (S : Type u) (R : Type u) s = e.inv.hom s := fun _ => rfl
  have hinv : ∀ s : (S : Type u), e.hom.hom (e.inv.hom s) = s := fun s =>
    congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) s) e.inv_hom_id
  haveI : IsScalarTower (S : Type u) (R : Type u) (I : Type u) := by
    refine ⟨fun s r x => ?_⟩
    show (algebraMap (R : Type u) (S : Type u) (algebraMap (S : Type u) (R : Type u) s * r)) • x
      = s • ((algebraMap (R : Type u) (S : Type u) r) • x)
    rw [heh, hhe, map_mul, hinv, mul_smul, heh]
  have hsurj : Function.Surjective (algebraMap (S : Type u) (R : Type u)) := by
    show Function.Surjective e.inv.hom
    exact e.commRingCatIsoToRingEquiv.symm.surjective
  haveI : Module.Invertible (R : Type u) (I : Type u) :=
    invertible_of_surjective_algebraMap (R := (S : Type u)) (S := (R : Type u))
      (M := (I : Type u)) hsurj
  refine Module.Invertible.congr (M := (I : Type u)) ?_
  refine LinearEquiv.symm
    { toFun := fun x => ⟨e.hom.hom x.1, x.2⟩
      map_add' := fun a b => Subtype.ext (map_add _ _ _)
      map_smul' := fun r x => Subtype.ext (by
        show e.hom.hom (r * x.1) = (algebraMap (R : Type u) (S : Type u) r) * e.hom.hom x.1
        rw [map_mul, heh])
      invFun := fun y => ⟨e.inv.hom y.1, by
        show e.hom.hom (e.inv.hom y.1) ∈ I
        rw [hinv]; exact y.2⟩
      left_inv := fun x => Subtype.ext (congrArg (fun (m : _ ⟶ _) =>
        (CommRingCat.Hom.hom m) x.1) e.hom_inv_id)
      right_inv := fun y => Subtype.ext (hinv y.1) }


/-! ## ★出典の紐付け(`.src`) -/

def invertible_comap_of_isIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——同型に沿った引き戻しでも可逆)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
