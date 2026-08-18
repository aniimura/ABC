import ABC3.Found.Arakelov.PicQCoh
import ABC3.Found.Arakelov.PicInvOfModule

/-!
# Arakelov (B1) 第 144 ブロック —— ★★★★★★**可逆層の大域切断は可逆加群**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`equivPicRing` の**逆向き**が繋がった

第 133 で「可逆加群 ⟹ 可逆層」、本ブロックでその**逆**である:

    L : InvSheaf (Spec R)  ⟹  Module.Invertible R (Γ L.carrier)

## ★★筋——`tilde` を通して `R` に戻す

    Γ(Lc) ⊗ Γ(Li)  ≅  Γ(tilde (Γ Lc ⊗ Γ Li))       (tildeGammaIso)
                   ≅  Γ(tilde (Γ Lc) ⊗ tilde (Γ Li)) (第 91、逆向き)
                   ≅  Γ(Lc ⊗ Li)                    (★第 143)
                   ≅  Γ(𝒪)  ≅  R                    (L.isInv、第 82)

★★第 143(局所自明 ⟹ `F ≅ (Γ F)~`)がここで効く——
**可逆層は `tilde` の本質像に入る**という一段である。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `gammaTensorIso` | ★★★`Γ(Lc) ⊗ Γ(Li) ≅ R` |
| `invertible_gammaCarrier` | ★★★★★★**大域切断は可逆加群** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open scoped TensorProduct

variable (R : CommRingCat.{u})


/-- ★可逆層の大域切断のテンソル積は `R` である。 -/
noncomputable def gammaTensorIso (L : InvSheaf (Spec R)) :
    ModuleCat.of (R : Type u)
        ((AlgebraicGeometry.moduleSpecΓFunctor.obj L.carrier : Type u) ⊗[(R : Type u)]
          (AlgebraicGeometry.moduleSpecΓFunctor.obj L.inv : Type u))
      ≅ ModuleCat.of (R : Type u) (R : Type u) :=
  (tildeGammaIso R).app _ ≪≫
    AlgebraicGeometry.moduleSpecΓFunctor.mapIso
      ((tildeTensorIso R _ _).symm
        ≪≫ tensorModulesIso (tildeGammaIsoOfTrivial R L.carrier L.trivial)
          (tildeGammaIsoOfTrivial R L.inv L.invTrivial)
        ≪≫ L.isInv.some ≪≫ (tildeUnit R).symm)
    ≪≫ ((tildeGammaIso R).app _).symm

/-- ★★★可逆層の大域切断は可逆加群である。 -/
theorem invertible_gammaCarrier (L : InvSheaf (Spec R)) :
    Module.Invertible (R : Type u)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj L.carrier : Type u) :=
  Module.Invertible.left (gammaTensorIso R L).toLinearEquiv


/-! ## ★出典の紐付け(`.src`) -/

def invertible_gammaCarrier.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆層の大域切断は可逆加群)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
