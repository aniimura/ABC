import ABC3.Found.Arakelov.PicLocEq

/-!
# Arakelov (B1) 第 141 ブロック —— **制限写像の短縮記法と等式版局所化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★`surj'` の証明は長い——先に道具を揃える

残る 1 欄 `surj'`(「`fⁿ y` は `⊤` から来る」)の証明は
**貼り合わせ + 2 段の最大値**で長くなる。★先に 4 つの道具を置く。

| 道具 | 用途 |
|---|---|
| `isLocalizedModule_secRes_eq` | ★開集合が**等式でしか**一致しないときの第 137 |
| `inf_specD` | ★`D(a) ⊓ D(b) = D(a·b)` |
| `secMap` | ★制限写像の短縮記法 |
| `secMap_trans` | ★制限の推移律 |

★★`IsCompatible` は `Opens.infLELeft` を使うので、
交わりを `D(a·b)` と**同定**する必要がある——そこで `inf_specD` と等式版が要る。
-/

universe u v

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (F : (Spec R).Modules)


/-- ★開集合を等式で指定した局所化。 -/
theorem isLocalizedModule_secRes_eq (V W : (Spec R).Opens) (g t : (R : Type u))
    (hV : V = specD R g) (hW : W = specD R (t * g)) (hle : W ≤ V)
    (e : (restrictPresheafFunctor (Spec R) (specD R g)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R g))) :
    IsLocalizedModule (Submonoid.powers t)
      ((modulesSpecToSheaf.obj F).obj.map (homOfLE hle).op).hom := by
  subst hV; subst hW
  exact isLocalizedModule_secRes R F g t e

/-- ★交わりは積の基本開集合。 -/
theorem inf_specD (a b : (R : Type u)) :
    specD R a ⊓ specD R b = specD R (a * b) := by
  show _ = PrimeSpectrum.basicOpen (a * b)
  rw [PrimeSpectrum.basicOpen_mul]
  rfl

/-- ★制限写像(短縮記法)。 -/
noncomputable abbrev secMap (U V : (Spec R).Opens) (h : V ≤ U) :
    (((modulesSpecToSheaf.obj F).obj.obj (op U)) : Type u) →ₗ[(R : Type u)]
      (((modulesSpecToSheaf.obj F).obj.obj (op V)) : Type u) :=
  ((modulesSpecToSheaf.obj F).obj.map (homOfLE h).op).hom

/-- ★制限の推移律。 -/
theorem secMap_trans (U V W : (Spec R).Opens) (h : V ≤ U) (h' : W ≤ V)
    (x : (((modulesSpecToSheaf.obj F).obj.obj (op U)) : Type u)) :
    secMap R F V W h' (secMap R F U V h x) = secMap R F U W (le_trans h' h) x := by
  show ((modulesSpecToSheaf.obj F).obj.map _).hom
    (((modulesSpecToSheaf.obj F).obj.map _).hom x) = _
  rw [← ConcreteCategory.comp_apply, ← Functor.map_comp]
  rfl


/-! ## ★出典の紐付け(`.src`) -/

def secMap_trans.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限写像の短縮記法)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
