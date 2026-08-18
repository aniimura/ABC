import ABC3.Found.Arakelov.PicUnitSurj
import ABC3.Found.Arakelov.PicLocalTrivial
import ABC3.Found.Arakelov.PicEvalBij

/-!
# Arakelov (B2) 第 207 ブロック —— ★★★★**自明化開の上では全射な射は同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-225 の補題 A

    L|_V ≅ 𝟙_、M|_V ≅ 𝟙_、φ.app V が全射  ⟹  φ|_V は同型

★★これは**可逆層を扱う後段すべてで効く**再利用可能な結果である
——「直線束の間の全射は同型」の局所版だからである。

## ★★機構は共役 1 回

    𝟙_ --eL.inv--> L|_V --φ|_V--> M|_V --eM.hom--> 𝟙_

★この共役 `conjToUnit` は単位の自己射であり、全射なら第 206 で**係数が単元**である。
`unitEndEquiv` は環同型なので、単元は同型に戻る。あとは

    φ|_V = eL.hom ≫ conjToUnit ≫ eM.inv

と分解すれば `IsIso` が合成で出る。

## ★★型の橋が 2 本要った([[typed-identity-bridge]])

`φ.app (op V)` の値は `M.obj (op V)` だが、`eM` が食うのは
`((restrictPresheafFunctor X V).obj M).obj (op (Over.mk (𝟙 V)))` である。
★第 173 の `fVal`(値の橋)と第 181 の `gVal`(逆向き)が**そのまま**当たった。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `conjToUnit` | ★共役して単位の自己射にする |
| `conjToUnit_surjective` | ★★全射性は共役で保たれる |
| `isIso_restrict_of_surjective` | ★★★★**自明化開の上では全射な射は同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} {L M : X.PresheafOfModules} (φ : L ⟶ M) (V : X.Opens)

/-- ★自明化を両側に噛ませて単位の自己射にする。 -/
noncomputable def conjToUnit
    (eL : (restrictPresheafFunctor X V).obj L ≅ 𝟙_ (PresheafModulesOn X V))
    (eM : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    End (𝟙_ (PresheafModulesOn X V)) :=
  eL.inv ≫ (restrictPresheafFunctor X V).map φ ≫ eM.hom

theorem conjToUnit_surjective
    (eL : (restrictPresheafFunctor X V).obj L ≅ 𝟙_ (PresheafModulesOn X V))
    (eM : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (hs : Function.Surjective ((φ.app (op V)).hom)) :
    Function.Surjective (((conjToUnit φ V eL eM).app (op (Over.mk (𝟙 V)))).hom) := by
  intro y
  obtain ⟨x, hx⟩ := hs (gVal M (op V) ((eM.inv.app (op (Over.mk (𝟙 V)))).hom y))
  refine ⟨(eL.hom.app (op (Over.mk (𝟙 V)))).hom (fVal L (op V) x), ?_⟩
  have h1 : (eL.inv.app (op (Over.mk (𝟙 V)))).hom
      ((eL.hom.app (op (Over.mk (𝟙 V)))).hom (fVal L (op V) x)) = fVal L (op V) x :=
    congrArg (fun (m : _ ⟶ _) => ((m.app (op (Over.mk (𝟙 V)))).hom (fVal L (op V) x))) eL.hom_inv_id
  show (eM.hom.app _).hom ((((restrictPresheafFunctor X V).map φ).app _).hom
    ((eL.inv.app _).hom ((eL.hom.app _).hom (fVal L (op V) x)))) = y
  rw [h1]
  have h2 : (((restrictPresheafFunctor X V).map φ).app (op (Over.mk (𝟙 V)))).hom
      (fVal L (op V) x) = fVal M (op V) ((φ.app (op V)).hom x) := rfl
  rw [h2, hx]
  show (eM.hom.app _).hom ((eM.inv.app _).hom y) = y
  exact congrArg (fun (m : _ ⟶ _) => ((m.app (op (Over.mk (𝟙 V)))).hom y)) eM.inv_hom_id




/-- ★★★★自明化開の上では、全射な射は同型である。 -/
theorem isIso_restrict_of_surjective
    (eL : (restrictPresheafFunctor X V).obj L ≅ 𝟙_ (PresheafModulesOn X V))
    (eM : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (hs : Function.Surjective ((φ.app (op V)).hom)) :
    IsIso ((restrictPresheafFunctor X V).map φ) := by
  have hu := isUnit_unitEndEquiv_of_surjective V (conjToUnit φ V eL eM)
  have hu := isUnit_unitEndEquiv_of_surjective V (conjToUnit φ V eL eM)
    (conjToUnit_surjective φ V eL eM hs)
  haveI : IsIso (conjToUnit φ V eL eM) := by
    obtain ⟨c, hc⟩ := hu
    have hcv : conjToUnit φ V eL eM = (unitEndEquiv V).symm (c : (Γ(X, V) : Type u)) := by
      rw [hc]; exact ((unitEndEquiv V).symm_apply_apply _).symm
    refine ⟨(unitEndEquiv V).symm ((c⁻¹ : (Γ(X, V) : Type u)ˣ) : (Γ(X, V) : Type u)), ?_, ?_⟩
    · rw [hcv]
      show (unitEndEquiv V).symm _ * (unitEndEquiv V).symm _ = 1
      rw [← map_mul, c.inv_mul, map_one]
    · rw [hcv]
      show (unitEndEquiv V).symm _ * (unitEndEquiv V).symm _ = 1
      rw [← map_mul, c.mul_inv, map_one]
  have hdecomp : (restrictPresheafFunctor X V).map φ
      = eL.hom ≫ conjToUnit φ V eL eM ≫ eM.inv := by
    rw [conjToUnit]
    simp
  rw [hdecomp]
  infer_instance

/-! ## ★出典の紐付け(`.src`) -/

def isIso_restrict_of_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自明化開の上では全射な射は同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
