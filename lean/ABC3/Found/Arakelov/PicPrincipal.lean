import ABC3.Found.Arakelov.PicDivisorMul
import ABC3.Found.Arakelov.PicPrincipalRing

/-!
# Arakelov (B2) 第 192 ブロック —— ★★★★★★**主因子と自明な直線束**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`ofDivisor_eq_one_iff` は**定義から出た**

    IsPrincipalDivisor X D := Nonempty (idealSheaf D ≅ 𝒪_X)

と定義すれば、`Pic` が同型類の商であることから

    ofDivisor X D = 1  ⟺  IsPrincipalDivisor X D      (Cartier のとき)

は `Quotient.exact` / `mk_eq_mk` の 2 行で出る。★★★**設計で消せる仕事**であった。

## ★★「主因子 ⟹ Cartier」も出る

`idealSheaf D ≅ 𝒪_X` なら各アフィン開 `A` で `D.ideal A ≅ Γ(A)`、
すなわち**自由階数 1** だから可逆である。★これで `isPrincipalDivisor_mul` に
第 190 の `mulModulesIso`(Cartier を要求する)を当てられる。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `modulesIsoApp` | ★層加群の同型は各開集合で加群の同型 |
| `IsPrincipalDivisor` | ★★★**主因子** |
| `isPrincipalDivisor_isCartier` | ★★主因子は Cartier |
| `ofDivisorSheaf_eq_one_iff` | ★★★★★★**`𝒪(D)` が自明 ⟺ `D` が主因子** |
| `isPrincipalDivisor_mul` | ★★★主因子は積で閉じる |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★層加群の同型は各開集合で加群の同型を与える。 -/
noncomputable def modulesIsoApp {M N : X.Modules} (e : M ≅ N) (U : (X.Opens)ᵒᵖ) :
    M.val.obj U ≅ N.val.obj U :=
  (PresheafOfModules.evaluation _ U).mapIso
    { hom := e.hom.val, inv := e.inv.val,
      hom_inv_id := congrArg SheafOfModules.Hom.val e.hom_inv_id,
      inv_hom_id := congrArg SheafOfModules.Hom.val e.inv_hom_id }

/-- ★★★主因子。 -/
def IsPrincipalDivisor (X : Scheme.{u}) (D : X.IdealSheafData) : Prop :=
  Nonempty (idealSheaf D ≅ unitModules X)

theorem isPrincipalDivisor_isCartier (D : X.IdealSheafData)
    (h : IsPrincipalDivisor X D) : IsCartier X D := by
  obtain ⟨e⟩ := h
  intro A
  haveI : Module.Invertible (Γ(X, A.1) : Type u)
      (((unitModules X).val.obj (op A.1)) : Type u) :=
    inferInstanceAs (Module.Invertible (Γ(X, A.1) : Type u) (Γ(X, A.1) : Type u))
  have hli := (modulesIsoApp e (op A.1)).toLinearEquiv
  refine cast (congrArg (fun (J : Ideal (Γ(X, A.1) : Type u)) =>
    Module.Invertible (Γ(X, A.1) : Type u) J) (idealSections_affine D A)) ?_
  exact Module.Invertible.congr hli.symm




theorem ofDivisorSheaf_eq_one_iff (D : X.IdealSheafData) (hD : IsCartier X D) :
    ofDivisorSheaf D = 1 ↔ IsPrincipalDivisor X D := by
  classical
  constructor
  · intro h
    obtain ⟨e⟩ := Quotient.exact h
    rw [divisorInvSheaf_carrier D hD] at e
    exact ⟨e⟩
  · rintro ⟨e⟩
    refine PicSheaf.mk_eq_mk ?_
    rw [divisorInvSheaf_carrier D hD]
    exact e

theorem isPrincipalDivisor_mul (D E : X.IdealSheafData)
    (hD : IsPrincipalDivisor X D) (hE : IsPrincipalDivisor X E) :
    IsPrincipalDivisor X (D * E) := by
  have hcD := isPrincipalDivisor_isCartier D hD
  obtain ⟨eD⟩ := hD
  obtain ⟨eE⟩ := hE
  exact ⟨(mulModulesIso D E hcD).symm ≪≫ tensorModulesIso eD eE
    ≪≫ tensorUnitLeft (unitModules X)⟩


/-! ## ★出典の紐付け(`.src`) -/

def ofDivisorSheaf_eq_one_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——𝒪(D) が自明なのは D が主因子のとき)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
