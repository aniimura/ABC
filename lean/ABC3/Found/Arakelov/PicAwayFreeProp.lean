import ABC3.Found.Arakelov.PicAwayModule
import Mathlib.RingTheory.LocalProperties.Projective

/-!
# Arakelov (B1) 第 99 ブロック —— **自由性は `D(g)` から `D(t·g)` へ運べる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★§9-101 の 5 段が閉じた

    M_g が R_g 上自由  ⟹  M_{t·g} が R_{t·g} 上自由

★★これで「可逆加群は基本開集合の上で自由」(第 76)が、**その部分開集合でも**成り立つ。

## ★★機構 —— 基底変換の言葉に移す

| 段 | 出典 |
|---|---|
| 1 | `M_g →ₗ M_{t·g}`(`IsLocalizedModule.lift`) | ★第 98 |
| 2 | `R_g` 線型化(`extendScalarsOfIsLocalization`) | ★本ブロック |
| 3 | `IsBaseChange R_g mk_g` / `IsBaseChange R_{t·g} mk_{t·g}` | ★mathlib |
| 4 | `IsBaseChange.comp_iff` | ★mathlib |
| 5 | `Module.free_of_isLocalizedModule` | ★mathlib |

★★★**「加群版の局所化推移」は mathlib に無いが、基底変換の言葉なら在る**
——§9-101 で測った通りである。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `awayMulTower'` | ★`R_g → R_{t·g} → M_{t·g}` は足場 |
| `liftAwayMapA` | ★`R_g` 線型版 |
| `isBaseChange_liftAwayMapA` | ★★基底変換であること |
| `free_away_mul` | ★★★★★★**自由性の伝播** |
-/

universe u

namespace ABC3.Found.Arakelov

variable (R : Type u) [CommRing R] (g t : R)
  (M : Type u) [AddCommGroup M] [Module R M]

/-- ★**`R_g → R_{t·g} → M_{t·g}` は足場をなす**。 -/
instance awayMulTower' : IsScalarTower (Localization (Submonoid.powers g))
    (Localization (Submonoid.powers (t * g)))
    (LocalizedModule (Submonoid.powers (t * g)) M) :=
  IsScalarTower.of_algebraMap_smul fun _ _ => rfl

/-- ★**`M_g →ₗ M_{t·g}` の `R_g` 線型版**。 -/
noncomputable def liftAwayMapA : LocalizedModule (Submonoid.powers g) M
    →ₗ[Localization (Submonoid.powers g)] LocalizedModule (Submonoid.powers (t * g)) M :=
  LinearMap.extendScalarsOfIsLocalization (Submonoid.powers g)
    (Localization (Submonoid.powers g)) (liftAwayMap R g t M)

/-- ★★**`liftAwayMapA` は基底変換である**。 -/
theorem isBaseChange_liftAwayMapA :
    IsBaseChange (Localization (Submonoid.powers (t * g))) (liftAwayMapA R g t M) := by
  have hg : IsBaseChange (Localization (Submonoid.powers g))
      (LocalizedModule.mkLinearMap (Submonoid.powers g) M) :=
    (isLocalizedModule_iff_isBaseChange (Submonoid.powers g) _ _).1 inferInstance
  refine (IsBaseChange.comp_iff hg).1 ?_
  have h : (liftAwayMapA R g t M : LocalizedModule (Submonoid.powers g) M →ₗ[R] _)
      ∘ₗ LocalizedModule.mkLinearMap (Submonoid.powers g) M
      = LocalizedModule.mkLinearMap (Submonoid.powers (t * g)) M := liftAwayMap_comp R g t M
  rw [h]
  exact (isLocalizedModule_iff_isBaseChange (Submonoid.powers (t * g)) _ _).1 inferInstance

/-- ★★★★★★**自由性は `D(g)` から `D(t·g)` へ運べる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで「可逆加群は基本開集合の上で自由」が、**その部分開集合でも**成り立つ。 -/
theorem free_away_mul
    [Module.Free (Localization (Submonoid.powers g)) (LocalizedModule (Submonoid.powers g) M)] :
    Module.Free (Localization (Submonoid.powers (t * g)))
      (LocalizedModule (Submonoid.powers (t * g)) M) := by
  haveI : IsLocalizedModule ((Submonoid.closure ({t, g} : Set R)).map
      (algebraMap R (Localization (Submonoid.powers g)))) (liftAwayMapA R g t M) := by
    haveI := isLocalization_closure_map R g t
    exact (isLocalizedModule_iff_isBaseChange _
      (Localization (Submonoid.powers (t * g))) _).2 (isBaseChange_liftAwayMapA R g t M)
  haveI := isLocalization_closure_map R g t
  exact Module.free_of_isLocalizedModule
    (Rₛ := Localization (Submonoid.powers (t * g)))
    ((Submonoid.closure ({t, g} : Set R)).map
      (algebraMap R (Localization (Submonoid.powers g))))
    (liftAwayMapA R g t M)

/-! ## ★出典の紐付け(`.src`) -/

def free_away_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自由性の伝播)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
