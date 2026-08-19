import ABC3.Found.Arakelov.PicResIso

/-!
# Arakelov (B2) 第 224 ブロック —— ★★★★★★**`hcompat` の核**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-254 の逃げ道(a)が効いた

§9-254 で「`rw` が証明部分まで照合するので `Subsingleton.elim` では足りない」と測り、
逃げ道の本命を「**第 222 を証明を引数に取る形で書き直す**」とした。★その通りであった。

    gammaSpec_appIso_inv' (h : fromSpec ''ᵁ (fromSpec ⁻¹ᵁ U) ≤ U) : … = X.presheaf.map (homOfLE h).op

★★証明 `h` を**外から渡す**ので、呼び出し側の項と**構文的に一致する**。
★★★これで `rw` が食い、あとは `simp` が同型を打ち消す。

## ★★摩擦の 6 例目の逃げ道が確定した

| # | 症状 | 逃げ道 |
|---|---|---|
| 6 | `rw` が証明部分まで照合する | ★**証明を引数に出す**(本ブロック) |

★★★これで [[ring-instance-two-paths]] の 6 例すべてに逃げ道がついた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gammaSpec_appIso_inv'` | ★★★第 222 の証明引数版 |
| `appIso_gammaSpec_comp` | ★★★★★★**`hcompat` の核** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} {U : X.Opens} (hU : IsAffineOpen U)

/-- ★★★**第 222 の証明引数版**——`rw` が食う形にする。 -/
theorem gammaSpec_appIso_inv' (h : (hU.fromSpec ''ᵁ (hU.fromSpec ⁻¹ᵁ U)) ≤ U) :
    (Scheme.ΓSpecIso (Γ(X, U))).inv
      ≫ (Spec Γ(X, U)).presheaf.map (eqToHom hU.fromSpec_preimage_self).op
      ≫ (hU.fromSpec.appIso (hU.fromSpec ⁻¹ᵁ U)).inv
      = X.presheaf.map (homOfLE h).op := by
  rw [Subsingleton.elim h (Set.image_preimage_subset hU.fromSpec U.1)]
  exact gammaSpec_appIso_inv hU


/-- ★★★★★★**`hcompat` の核**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★§9-254 で特定した逃げ道(証明を引数に出す)がそのまま効いた。 -/
theorem appIso_gammaSpec_comp (h : (hU.fromSpec ''ᵁ (hU.fromSpec ⁻¹ᵁ U)) ≤ U) :
    (hU.fromSpec.appIso (hU.fromSpec ⁻¹ᵁ U)).hom
        ≫ (Spec Γ(X, U)).presheaf.map (eqToHom hU.fromSpec_preimage_self.symm).op
        ≫ (Scheme.ΓSpecIso (Γ(X, U))).hom
      = @inv _ _ _ _ _ (isIso_res_of_eq h (image_preimage_self hU).symm) := by
  haveI := isIso_res_of_eq h (image_preimage_self hU).symm
  refine (IsIso.inv_eq_of_hom_inv_id ?_).symm
  rw [← gammaSpec_appIso_inv' hU h]
  simp [← Scheme.Hom.appIso_inv_naturality]


/-! ## ★出典の紐付け(`.src`) -/

def appIso_gammaSpec_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——hcompat の核)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
