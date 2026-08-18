import ABC3.Found.Arakelov.PicOverPresieve
import Mathlib.Algebra.Module.LocalizedModule.IsLocalization

/-!
# Arakelov (B1) 第 123 ブロック —— **局所化した生成元の乗法は全単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の最後の材料

`M_g ≅ R_g`(生成元の乗法)を `D(t·g)` へ局所化すると、
`M_{t·g} ≅ R_{t·g}`(制限した生成元の乗法)になる。

★★機構は第 112 ブロック(局所化は全単射を保つ)そのものである。

## ★★詰まり —— **statement の elaboration に instance が要る**

★`IsLocalizedModule.map` は**言明の中で**インスタンスを要求するので、
`haveI` を証明の中に置いても**遅い**(2026-08-24 実測)。
★★`instance` として**先に**宣言する必要がある。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `isLocalizationClosureMapInst` | ★係数側の局所化(instance) |
| `isLocalizedModuleLiftAwayInst` | ★加群側の局所化(instance) |
| `bijective_localizedGenerator` | ★★★★★**局所化した生成元の乗法は全単射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : Type u) [CommRing R] (g t : R)
  (M : Type u) [AddCommGroup M] [Module R M]

/-- ★**係数側の局所化**(instance 化)。 -/
instance isLocalizationClosureMapInst :
    IsLocalization ((Submonoid.closure ({t, g} : Set R)).map
        (algebraMap R (Localization (Submonoid.powers g))))
      (Localization (Submonoid.powers (t * g))) :=
  isLocalization_closure_map R g t

/-- ★**加群側の局所化**(instance 化)。 -/
instance isLocalizedModuleLiftAwayInst :
    IsLocalizedModule ((Submonoid.closure ({t, g} : Set R)).map
      (algebraMap R (Localization (Submonoid.powers g)))) (liftAwayMapA R g t M) :=
  (isLocalizedModule_iff_isBaseChange _
    (Localization (Submonoid.powers (t * g))) _).2 (isBaseChange_liftAwayMapA R g t M)

/-- ★★★★★**局所化した生成元の乗法は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `M_{t·g} ≅ R_{t·g}` が**制限した生成元で**与えられる。 -/
theorem bijective_localizedGenerator
    (e : LocalizedModule (Submonoid.powers g) M
      ≃ₗ[Localization (Submonoid.powers g)] Localization (Submonoid.powers g)) :
    Function.Bijective (IsLocalizedModule.map
      ((Submonoid.closure ({t, g} : Set R)).map
        (algebraMap R (Localization (Submonoid.powers g))))
      (Algebra.linearMap (Localization (Submonoid.powers g))
        (Localization (Submonoid.powers (t * g))))
      (liftAwayMapA R g t M) e.symm.toLinearMap) :=
  bijective_localizedMap _ _ _ _ e.symm.bijective

/-! ## ★出典の紐付け(`.src`) -/

def bijective_localizedGenerator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所化した生成元の乗法は全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
