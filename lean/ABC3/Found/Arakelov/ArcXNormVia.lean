import ABC3.Found.Arakelov.ArcXNorm
import ABC3.Found.Arakelov.ArcLocalDirect
import ABC3.Found.Arakelov.ArcContTransport
import ABC3.Found.Arakelov.ArcOpenImmersion

/-!
# Arakelov (C3) 第 277 ブロック —— ★★★★**依存書き換えを回避する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-320 の残り(2)を解いた

`xNorm V F e (q ≫ V.ι) h w = localNorm V F e q (…)` を示したいが、
`rw [liftToOpenOfTop_comp]` は **motive is not type correct** で落ちる
——`arcFiber p F` が `p` に依存する型だからである。

★★★逃げ道は**選択を明示的な引数にする**こと:

    xNormVia V F e p (r) (hr : r ≫ V.ι = p) w := localNorm V F e r (…)

| 段 | 内容 |
|---|---|
| `xNormVia_indep` | ★★選択に依らない(★`subst` で `r₁ = r₂` を潰す) |
| `xNormVia_self` | ★`r := q`, `hr := rfl` では `eqToIso rfl` が消える(`rfl`) |
| `xNorm_eq_via` | ★`xNorm` は `xNormVia` の特別な場合(`rfl`) |
| `xNorm_comp` | ★★★★合成則 |

★★★★**依存書き換えを「選択の独立性」に置き換える**——
`rw` は型を書き換えられないが、`subst` は**束縛変数**なら書き換えられる。
★選択を引数にすれば束縛変数になる。

## ★これは定石として書ける

**「`f (choice p)` を書き換えたい」ときは、`choice` を引数に出して
「選択に依らない」を `subst` で示し、都合のよい選択で instantiate する。**
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★持ち上げを明示的に受けたノルム。 -/
noncomputable def xNormVia (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (r : Spec (CommRingCat.of ℂ) ⟶ V.toScheme) (hr : r ≫ V.ι = p)
    (w : ↥(arcFiber p F)) : ℝ :=
  localNorm V F e r ((eqToIso (congrArg (fun t => arcFiber t F) hr.symm)
    ≪≫ arcFiberFactor V.ι F r).hom.hom w)

/-- ★★持ち上げの取り方に依らない。 -/
theorem xNormVia_indep (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (r₁ r₂ : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (hr₁ : r₁ ≫ V.ι = p) (hr₂ : r₂ ≫ V.ι = p) (w : ↥(arcFiber p F)) :
    xNormVia V F e p r₁ hr₁ w = xNormVia V F e p r₂ hr₂ w := by
  have hr : r₁ = r₂ := comp_openImmersion_injective V.ι (hr₁.trans hr₂.symm)
  subst hr
  rfl

/-- ★★★`q ≫ V.ι` での値。 -/
theorem xNormVia_self (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (w : ↥(arcFiber (q ≫ V.ι) F)) :
    xNormVia V F e (q ≫ V.ι) q rfl w
      = localNorm V F e q ((arcFiberFactor V.ι F q).hom.hom w) :=
  rfl


/-- ★`xNorm` は `xNormVia` である。 -/
theorem xNorm_eq_via (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (w : ↥(arcFiber p F)) :
    xNorm V F e p h w
      = xNormVia V F e p (liftToOpenOfTop V p h) (liftToOpen_fac V p h) w :=
  rfl

/-- ★★★★依存書き換えを回避した合成則。 -/
theorem xNorm_comp (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (w : ↥(arcFiber (q ≫ V.ι) F)) :
    xNorm V F e (q ≫ V.ι) (comp_preimage_eq_top V q) w
      = localNorm V F e q ((arcFiberFactor V.ι F q).hom.hom w) := by
  rw [xNorm_eq_via]
  rw [xNormVia_indep V F e (q ≫ V.ι) _ q _ rfl]
  exact xNormVia_self V F e q w


/-! ## ★出典の紐付け(`.src`) -/

def xNorm_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——依存書き換えを選択の独立性に置き換える)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
