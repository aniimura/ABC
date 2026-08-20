import ABC3.Found.Arakelov.ArcTrivMetric

/-!
# Arakelov (C3) 第 254 ブロック —— **開埋め込みに沿ってノルムを運ぶ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★貼り合わせの前段 —— 「制限した層の自明化」から `X` 側のノルムへ

第 253 で**自明な層**の連続な計量が取れた。一般の `L` はアフィン開 `U` の上でしか
自明にならないので、`U` の側で作ったノルムを `X` の側へ運ぶ必要がある。

★★機構は 2 つの同型を繋ぐだけである:

    arcFiber (p ≫ j) L ≅ arcFiber p ((pullback j).obj L) ≅ arcFiber p (restrict L j)

| 段 | 使うもの |
|---|---|
| 引き戻しの合成 | ★第 65 ブロック `pullbackCompApp` |
| 制限 = 引き戻し | ★★mathlib `Scheme.Modules.restrictFunctorIsoPullback` |

## ★★★在庫の測定(2026-08-19)—— mathlib が半分持っていた

`Mathlib/AlgebraicGeometry/Modules/Sheaf.lean` に:

| 在庫 | 内容 |
|---|---|
| `Scheme.Modules.restrict` | ★開埋め込みに沿った制限 |
| `restrictAppIso` | ★★切断は `Γ(M, f ''ᵁ U)`(**`Iso.refl`**) |
| `restrictFunctorIsoPullback` | ★★★**制限 ≅ 引き戻し** |
| `restrictAdjunction` | ★押し出しの右随伴 |

★★★★したがって残るのは「`IsInvertibleSheaf` の(`Over V` 上の)自明化から
`restrict L V.ι ≅ 𝒪_V` を作る」段だけである——**`V.Opens` と `Over V` の対応**を書く仕事。

| 定義・定理 | 内容 |
|---|---|
| `arcFiberFactor` | ★★★ファイバーの分解 |
| `restrNorm` | ★★運んだノルム |
| `restrNorm_nonneg` / `_smul` / `_eq_zero_iff` | ★3 法則(同型で移すだけ) |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {U X : Scheme.{0}} (j : U ⟶ X) [IsOpenImmersion j] (L : X.Modules)

/-- ★★★開埋め込みで分解した点のファイバーは、制限した層のファイバーである。 -/
noncomputable def arcFiberFactor (p : Spec (CommRingCat.of ℂ) ⟶ U) :
    arcFiber (p ≫ j) L ≅ arcFiber p (Scheme.Modules.restrict L j) :=
  (moduleSpecΓFunctor (R := CommRingCat.of ℂ)).mapIso
    (pullbackCompApp p j L ≪≫ (Scheme.Modules.pullback p).mapIso
      ((Scheme.Modules.restrictFunctorIsoPullback j).app L).symm)

variable (t : Scheme.Modules.restrict L j ≅ unitModules U)

/-- ★★制限の自明化から、`X` 側のファイバーのノルムを作る。 -/
noncomputable def restrNorm (p : Spec (CommRingCat.of ℂ) ⟶ U)
    (w : ↥(arcFiber (p ≫ j) L)) : ℝ :=
  trivNorm t p ((arcFiberFactor j L p).hom.hom w)

theorem restrNorm_nonneg (p : Spec (CommRingCat.of ℂ) ⟶ U)
    (w : ↥(arcFiber (p ≫ j) L)) : 0 ≤ restrNorm j L t p w :=
  trivNorm_nonneg t p _

theorem restrNorm_smul (p : Spec (CommRingCat.of ℂ) ⟶ U)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber (p ≫ j) L)) :
    restrNorm j L t p (c • w) = ‖c‖ * restrNorm j L t p w := by
  show trivNorm t p ((arcFiberFactor j L p).hom.hom (c • w)) = _
  rw [map_smul]
  exact trivNorm_smul t p c _

theorem restrNorm_eq_zero_iff (p : Spec (CommRingCat.of ℂ) ⟶ U)
    (w : ↥(arcFiber (p ≫ j) L)) : restrNorm j L t p w = 0 ↔ w = 0 := by
  have hb : Function.Bijective (((arcFiberFactor j L p).hom).hom) :=
    ConcreteCategory.bijective_of_isIso _
  show trivNorm t p ((arcFiberFactor j L p).hom.hom w) = 0 ↔ w = 0
  rw [trivNorm_eq_zero_iff]
  constructor
  · intro h
    refine hb.1 ?_
    rw [h]
    exact (map_zero _).symm
  · intro h
    rw [h]
    exact map_zero _


/-! ## ★出典の紐付け(`.src`) -/

def restrNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みに沿ってノルムを運ぶ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
