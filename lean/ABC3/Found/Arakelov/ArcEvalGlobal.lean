import ABC3.Found.Arakelov.ArcUnitEval
import ABC3.Found.Arakelov.ArcTopologyAffine
import ABC3.Found.Arakelov.ArcContCriterion

/-!
# Arakelov (C3) 第 252 ブロック —— ★★★★★**大域の正則関数は `X^arc` 上連続**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★これが「連続関数 `|s|_L`」の内実である

自明化 `t : L|_V ≅ 𝒪_V` を使うとノルムは `‖s‖(p) = |(t s)(p)|` になる。
★第 251 で右辺が `p^♯(t s)` であることを示した。
★★本ブロックで **`p ↦ p^♯(g)` が連続**であることを示す——これで連続性が閉じる。

## ★★機構は 3 段

| 段 | 内容 |
|---|---|
| `evalGlobal_spec` | ★アフィンでは `evalAffine` と一致(`ΓSpecIso_naturality` + `Spec.map_preimage`) |
| `continuous_evalGlobal_affine` | ★★アフィンでの連続性(第 5 ブロックの `continuous_evalAffine`) |
| `continuous_evalGlobal` | ★★★★★**一般の `X`**(第 249 の chart 判定) |

★★★アフィンの段は **mathlib の 2 行**で出た:

    (Spec.map f).appTop ≫ (ΓSpecIso S).hom = (ΓSpecIso R).hom ≫ f     (`ΓSpecIso_naturality`)
    Spec.map (Spec.preimage q) = q                                       (`Spec.map_preimage`)

★★★★一般の段は chart ごとに `isoSpec` で輸送し、`isoSpec.inv ≫ isoSpec.hom = 𝟙` を使うだけである。

| 定義・定理 | 内容 |
|---|---|
| `evalGlobal` | ★大域切断の複素点での値 |
| `evalGlobal_spec` | ★★アフィンでは `evalAffine` |
| `continuous_evalGlobal` | ★★★★★**大域の正則関数は連続** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}}

/-- ★大域切断の複素点での値。 -/
noncomputable def evalGlobal (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (g : ((X.presheaf.obj (op ⊤)) : Type)) : ℂ :=
  (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((Scheme.Hom.appTop p).hom g)

/-- ★★アフィンでは `evalAffine` と一致する。 -/
theorem evalGlobal_spec (A : CommRingCat.{0}) (q : Spec (CommRingCat.of ℂ) ⟶ Spec A)
    (b : (((Spec A).presheaf.obj (op ⊤)) : Type)) :
    evalGlobal q b = evalAffine A q ((Scheme.ΓSpecIso A).hom.hom b) := by
  have h := Scheme.ΓSpecIso_naturality (Spec.preimage q)
  rw [Spec.map_preimage] at h
  exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) b) h

/-- ★★★アフィンでは連続である。 -/
theorem continuous_evalGlobal_affine (A : CommRingCat.{0})
    (b : (((Spec A).presheaf.obj (op ⊤)) : Type)) :
    @Continuous _ _ (arcTopologyAffine A) _ (fun q => evalGlobal q b) := by
  have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ Spec A => evalGlobal q b)
      = (fun q => evalAffine A q ((Scheme.ΓSpecIso A).hom.hom b)) :=
    funext fun q => evalGlobal_spec A q b
  rw [he]
  exact continuous_evalAffine A _

/-- ★★★★★**大域の正則関数は `X^arc` 上連続である**。 -/
theorem continuous_evalGlobal (g : ((X.presheaf.obj (op ⊤)) : Type)) :
    @Continuous _ ℂ (arcTopology X) _ (fun p => evalGlobal p g) := by
  refine continuous_of_charts_factor _
    (fun U q => evalGlobal q ((Scheme.Hom.appTop U.2.isoSpec.inv).hom
      ((Scheme.Hom.appTop U.1.ι).hom g)))
    (fun U => continuous_evalGlobal_affine _ _) (fun U p => ?_)
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((Scheme.Hom.appTop (p ≫ U.1.ι)).hom g)
    = (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((Scheme.Hom.appTop (p ≫ U.2.isoSpec.hom)).hom
        ((Scheme.Hom.appTop U.2.isoSpec.inv).hom ((Scheme.Hom.appTop U.1.ι).hom g)))
  have hid : (Scheme.Hom.appTop U.2.isoSpec.inv) ≫ (Scheme.Hom.appTop U.2.isoSpec.hom)
      = 𝟙 _ := by
    rw [← Scheme.Hom.comp_appTop, U.2.isoSpec.hom_inv_id]
    rfl
  exact congrArg (fun (m : _ ⟶ _) => (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
    ((Scheme.Hom.appTop p).hom ((CommRingCat.Hom.hom m) ((Scheme.Hom.appTop U.1.ι).hom g))))
    hid.symm


/-! ## ★出典の紐付け(`.src`) -/

def continuous_evalGlobal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——大域の正則関数は X^arc 上連続)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
