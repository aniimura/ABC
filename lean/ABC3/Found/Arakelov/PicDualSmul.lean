import ABC3.Found.Arakelov.PicUnitRing
import ABC3.Found.Arakelov.PicResTrans

/-!
# Arakelov (B2) 第 167 ブロック —— **双対の作用は合成**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★双対を前層に組む前の 2 補題

第 166 の加群構造は `Module.compHom` 経由なので、
`c • φ` の**具体形**を出しておく必要がある:

    c • φ = φ ≫ unitMul U c

★`RingEquiv.ofBijective` の `symm` は `Function.surjInv` なので
**定義からは `unitMul` と一致しない**——単射性で示す。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `unitEndEquiv_symm_apply` | ★環同型の逆は `unitMul` |
| `dual_smul_eq` | ★★**作用は合成**(`c • φ = φ ≫ unitMul c`) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules)

/-- ★環同型の逆は `unitMul`。 -/
theorem unitEndEquiv_symm_apply (U : X.Opens) (c : (Γ(X, U) : Type u)) :
    (unitEndEquiv U).symm c = unitMul U c := by
  apply (unitEndEquiv U).injective
  rw [RingEquiv.apply_symm_apply]
  exact (unitEnd_unitMul U c).symm

/-- ★★作用は合成である。 -/
theorem dual_smul_eq (U : X.Opens) (c : (Γ(X, U) : Type u))
    (φ : (restrictPresheafFunctor X U).obj F ⟶ 𝟙_ (PresheafModulesOn X U)) :
    c • φ = φ ≫ unitMul U c := by
  show ((unitEndEquiv U).symm c) • φ = _
  rw [unitEndEquiv_symm_apply]
  rfl


/-! ## ★出典の紐付け(`.src`) -/

def dual_smul_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——双対の作用は合成)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
