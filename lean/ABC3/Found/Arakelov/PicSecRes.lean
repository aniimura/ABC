import ABC3.Found.Arakelov.PicSpecRes
import ABC3.Found.Arakelov.PicSectionTrivial
import ABC3.Found.Arakelov.PicBasicTrivial

/-!
# Arakelov (B1) 第 137 ブロック —— **自明な開集合での制限は局所化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★自明化の四角形で `𝒪` 側から `F` 側へ運ぶ

    Γ(F, D g) --res--> Γ(F, D(t·g))
        ≅ e₁               ≅ e₂
    𝒪(D g)   --res--> 𝒪(D(t·g))

★下段は第 136(`isLocalizedModule_specRes`)。
★★可換性は **`e.hom.naturality`** そのもの
——`Over.homMk (homOfLE h) : Over.mk (homOfLE h) ⟶ Over.mk (𝟙 (D g))` に沿って取る。

★★★上下の縦線が全単射なので
`IsLocalizedModule.comp_iff_of_bijective_left/right` で上段へ移る。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `secModOn` / `secTowerOn` | ★切断の `𝒪(U)` 加群構造と足場 |
| `secLinOfTrivial` | ★★自明化から `R` 線型写像 |
| `bijective_secLinOfTrivial` | ★全単射 |
| `secEquivOfTrivial` | ★`R` 線型同型 |
| `secLin_comm` | ★★★**四角形の可換性**(`naturality` 一発) |
| `isLocalizedModule_secRes` | ★★★★★**`F` 側の局所化** |

## ★★★型の 2 経路(記録、[[ring-instance-two-paths]] の 6・7 例目)

| 症状 | 逃げ道 |
|---|---|
| `(modulesSpecToSheaf.obj F).obj.obj (op U)` に `𝒪(U)` 作用が無い | ★`inferInstanceAs` |
| `ModuleCat.of Γ Γ` に `ringCatSheaf` 側の作用が無い | ★★`letI` + `inferInstanceAs`(**証明の中**で) |
| `r • y = algebraMap r * y` が `rfl` で**ない** | ★`Algebra.smul_def` |
| `IsScalarTower` の `smul_assoc` が `rfl` で**ない** | ★`IsScalarTower.of_algebraMap_smul` なら `rfl` |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (F : (Spec R).Modules)

noncomputable instance secModOn (U : (Spec R).Opens) :
    Module (Γ(Spec R, U) : Type u) (((modulesSpecToSheaf.obj F).obj.obj (op U)) : Type u) :=
  inferInstanceAs (Module (Γ(Spec R, U) : Type u) ((F.val.obj (op U)) : Type u))

instance secTowerOn (U : (Spec R).Opens) :
    IsScalarTower (R : Type u) (Γ(Spec R, U) : Type u)
      (((modulesSpecToSheaf.obj F).obj.obj (op U)) : Type u) := by
  refine IsScalarTower.of_algebraMap_smul fun r x => ?_
  rfl

/-- ★自明化から切断の `R` 線型写像。 -/
noncomputable def secLinOfTrivial (U : (Spec R).Opens)
    (e : (restrictPresheafFunctor (Spec R) U).obj F.val ≅ 𝟙_ (PresheafModulesOn (Spec R) U)) :
    (((modulesSpecToSheaf.obj F).obj.obj (op U)) : Type u) →ₗ[(R : Type u)]
      (Γ(Spec R, U) : Type u) where
  toFun := fun x => (sectionIsoOfTrivial U F.val e).hom.hom x
  map_add' := fun x y => map_add _ _ _
  map_smul' := fun r x => by
    show (sectionIsoOfTrivial U F.val e).hom.hom (r • x) = r • _
    rw [Algebra.smul_def]
    letI : Module (((Spec R).ringCatSheaf.obj.obj (op U)) : Type u)
        ((Γ(Spec R, U)) : Type u) :=
      inferInstanceAs (Module ((Γ(Spec R, U)) : Type u) ((Γ(Spec R, U)) : Type u))
    exact (sectionIsoOfTrivial U F.val e).hom.hom.map_smul _ x

theorem bijective_secLinOfTrivial (U : (Spec R).Opens)
    (e : (restrictPresheafFunctor (Spec R) U).obj F.val ≅ 𝟙_ (PresheafModulesOn (Spec R) U)) :
    Function.Bijective (secLinOfTrivial R F U e) :=
  ConcreteCategory.bijective_of_isIso (sectionIsoOfTrivial U F.val e).hom

noncomputable def secEquivOfTrivial (U : (Spec R).Opens)
    (e : (restrictPresheafFunctor (Spec R) U).obj F.val ≅ 𝟙_ (PresheafModulesOn (Spec R) U)) :
    (((modulesSpecToSheaf.obj F).obj.obj (op U)) : Type u) ≃ₗ[(R : Type u)]
      (Γ(Spec R, U) : Type u) :=
  LinearEquiv.ofBijective (secLinOfTrivial R F U e) (bijective_secLinOfTrivial R F U e)

theorem secLin_comm (g t : (R : Type u))
    (e : (restrictPresheafFunctor (Spec R) (specD R g)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R g))) :
    (secLinOfTrivial R F (specD R (t * g)) (trivialOfLe F.val (specDle R g t) e)).comp
        ((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle R g t)).op).hom
      = (specResAlgHom R g t).toLinearMap.comp (secLinOfTrivial R F (specD R g) e) := by
  refine LinearMap.ext fun x => ?_
  have hnat := e.hom.naturality
    (Over.homMk (homOfLE (specDle R g t))
      : (Over.mk (homOfLE (specDle R g t)) : Over (specD R g)) ⟶ Over.mk (𝟙 (specD R g))).op
  exact DFunLike.congr_fun (congrArg ModuleCat.Hom.hom hnat) x

theorem isLocalizedModule_secRes (g t : (R : Type u))
    (e : (restrictPresheafFunctor (Spec R) (specD R g)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R g))) :
    IsLocalizedModule (Submonoid.powers t)
      ((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle R g t)).op).hom := by
  refine (IsLocalizedModule.comp_iff_of_bijective_left (Submonoid.powers t)
    (secLinOfTrivial R F (specD R (t * g)) (trivialOfLe F.val (specDle R g t) e))
    (bijective_secLinOfTrivial R F _ _)).1 ?_
  rw [secLin_comm R F g t e]
  exact (IsLocalizedModule.comp_iff_of_bijective_right (Submonoid.powers t) _
    (bijective_secLinOfTrivial R F _ e)).2 (isLocalizedModule_specRes R g t)


/-! ## ★出典の紐付け(`.src`) -/

def isLocalizedModule_secRes.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自明な開集合での制限は局所化)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
