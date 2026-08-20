import ABC3.Found.Arakelov.ArcLocalDirect
import ABC3.Found.Arakelov.ArcLiftOpen

/-!
# Arakelov (C3) 第 271 ブロック —— ★★★★**`X` の点での局所ノルム**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★貼り合わせのために `X` の点へ再指標付けする

第 270 の `localNorm` は **`V.toScheme` の点**で添字づけられている。
★計量は `X` のすべての点で要るので、`p ⁻¹ᵁ V = ⊤` なる `p` に対し
`liftToOpenOfTop`(第 259)で降ろす。

★★ファイバーの同型は 2 段:

    arcFiber p F  ≅(eqToIso, liftToOpen_fac)  arcFiber (q ≫ V.ι) F
                  ≅(第 254 arcFiberFactor)     arcFiber q (restrict F V.ι)

★★★3 法則は同型で移すだけである。

| 定義・定理 | 内容 |
|---|---|
| `arcFiberAt` | ★ファイバーを `V.toScheme` 側へ |
| `xNorm` | ★★`X` の点での局所ノルム |
| `xNorm_nonneg` / `_smul` / `_eq_zero_iff` | ★★★3 法則 |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★`X` の点のファイバーを `V.toScheme` の側へ移す。 -/
noncomputable def arcFiberAt (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤) :
    arcFiber p F ≅ arcFiber (liftToOpenOfTop V p h) (Scheme.Modules.restrict F V.ι) :=
  eqToIso (congrArg (fun r => arcFiber r F) (liftToOpen_fac V p h).symm)
    ≪≫ arcFiberFactor V.ι F (liftToOpenOfTop V p h)

/-- ★★`X` の点での局所ノルム。 -/
noncomputable def xNorm (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (w : ↥(arcFiber p F)) : ℝ :=
  localNorm V F e (liftToOpenOfTop V p h) ((arcFiberAt V F p h).hom.hom w)

theorem xNorm_nonneg (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (w : ↥(arcFiber p F)) : 0 ≤ xNorm V F e p h w :=
  localNorm_nonneg V F e _ _

theorem xNorm_smul (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber p F)) :
    xNorm V F e p h (c • w) = ‖c‖ * xNorm V F e p h w := by
  show localNorm V F e _ ((arcFiberAt V F p h).hom.hom (c • w)) = _
  rw [map_smul]
  exact localNorm_smul V F e _ c _

theorem xNorm_eq_zero_iff (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (w : ↥(arcFiber p F)) : xNorm V F e p h w = 0 ↔ w = 0 := by
  have hb : Function.Bijective (((arcFiberAt V F p h).hom).hom) :=
    ConcreteCategory.bijective_of_isIso _
  show localNorm V F e _ ((arcFiberAt V F p h).hom.hom w) = 0 ↔ w = 0
  rw [localNorm_eq_zero_iff]
  constructor
  · intro hz
    refine hb.1 ?_
    rw [hz]
    exact (map_zero _).symm
  · intro hz
    rw [hz]
    exact map_zero _


/-! ## ★出典の紐付け(`.src`) -/

def xNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——X の点での局所ノルム)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
