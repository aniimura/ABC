import ABC3.Found.Arakelov.PicAwayUnit

/-!
# Arakelov (B1) 第 98 ブロック —— **`M_{t·g}` を `R_g` 加群と見る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★基底変換に載せるための土台

§9-101 の 5 段のうち、`IsBaseChange` を使うには
**`M_{t·g}` が `R_g` 加群でなければならない**。

★★`Algebra R_g R_{t·g}`(第 95 ブロック)と `Module R_{t·g} M_{t·g}` から
`Module.compHom` で作れる——★**instance では出ない**(実測)ので手で置く。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `awayMulModule` | ★`M_{t·g}` は `R_g` 加群 |
| `awayMulModuleTower` | ★★**`R → R_g → M_{t·g}` は足場** |
| `liftAwayMap` | ★★★**`M_g →ₗ M_{t·g}`**(第 97 の可逆性から) |
| `liftAwayMap_comp` | ★合成は `mk` に戻る |

## ★★★次

★`IsBaseChange.comp_iff` で `IsBaseChange R_{t·g} (liftAwayMap)` を出し、
`isLocalizedModule_iff_isBaseChange` で `IsLocalizedModule` に戻し、
`Module.free_of_isLocalizedModule` で自由性を運ぶ。
-/

universe u

namespace ABC3.Found.Arakelov

variable (R : Type u) [CommRing R] (g t : R)
  (M : Type u) [AddCommGroup M] [Module R M]

/-- ★**`M_{t·g}` は `R_g` 加群である**。 -/
noncomputable instance awayMulModule : Module (Localization (Submonoid.powers g))
    (LocalizedModule (Submonoid.powers (t * g)) M) :=
  Module.compHom _ (algebraMap (Localization (Submonoid.powers g))
    (Localization (Submonoid.powers (t * g))))

/-- ★★**`R → R_g → M_{t·g}` は足場をなす**。 -/
instance awayMulModuleTower : IsScalarTower R (Localization (Submonoid.powers g))
    (LocalizedModule (Submonoid.powers (t * g)) M) := by
  refine IsScalarTower.of_algebraMap_smul fun r x => ?_
  show (algebraMap (Localization (Submonoid.powers g))
    (Localization (Submonoid.powers (t * g)))
      (algebraMap R (Localization (Submonoid.powers g)) r)) • x = r • x
  rw [← IsScalarTower.algebraMap_apply, ← IsScalarTower.algebraMap_smul
    (Localization (Submonoid.powers (t * g))) r x]

/-- ★★★**`M_g →ₗ M_{t·g}`**——第 97 ブロックの可逆性から普遍性で作る。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
noncomputable def liftAwayMap : LocalizedModule (Submonoid.powers g) M →ₗ[R]
    LocalizedModule (Submonoid.powers (t * g)) M :=
  IsLocalizedModule.lift (Submonoid.powers g)
    (LocalizedModule.mkLinearMap (Submonoid.powers g) M)
    (LocalizedModule.mkLinearMap (Submonoid.powers (t * g)) M)
    (isUnit_end_powers R g t M)

/-- ★**合成は `mk` に戻る**。 -/
theorem liftAwayMap_comp :
    liftAwayMap R g t M ∘ₗ LocalizedModule.mkLinearMap (Submonoid.powers g) M
      = LocalizedModule.mkLinearMap (Submonoid.powers (t * g)) M :=
  IsLocalizedModule.lift_comp _ _ _ _

/-! ## ★出典の紐付け(`.src`) -/

def liftAwayMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——M_g から M_{t·g} への写像)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
