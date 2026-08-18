import ABC3.Found.Arakelov.PicSurjPair

/-!
# Arakelov (B2) 第 222 ブロック —— ★★★**`ΓSpecIso` と `appIso` の関係**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-250 で特定した筋がそのまま当たった

    ΓSpecIso.inv ≫ 転送 ≫ (appIso (fromSpec⁻¹U)).inv = 制限 Γ(X,U) → Γ(X, fromSpec''ᵁ(fromSpec⁻¹U))

★★mathlib の 2 本を繋ぐだけである:

| 補題 | 場所 |
|---|---|
| ★`Scheme.Hom.app_appIso_inv` | `OpenImmersion.lean:210` |
| ★`IsAffineOpen.fromSpec_app_self` | `AffineScheme.lean:569` |

## ★★§9-249 で「未解決」と書いた所が 1 ブロックで済んだ

★§9-249 では「依存位置の転送で、逃げ道は未解決」と書いた。
★★§9-250 で `external/_refs` を grep し直したら**2 本とも在った**——
そして本ブロックは `rw` 2 回と `simp` 1 回である。

★★★**「無い」と思ったら grep し直す**——この区間 3 度目の教訓である
(§9-235 の `complEDS₂`、§9-237 の写像可換性、本節)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `app_appIso_top_inv` | ★`app_appIso_inv` を `fromSpec` に当てた形 |
| `gammaSpec_appIso_inv` | ★★★**`ΓSpecIso` と `appIso` の関係** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} {U : X.Opens} (hU : IsAffineOpen U)

/-- ★`app_appIso_inv` を `fromSpec` に当てた形。 -/
theorem app_appIso_top_inv :
    hU.fromSpec.app U ≫ (hU.fromSpec.appIso (hU.fromSpec ⁻¹ᵁ U)).inv
      = X.presheaf.map (homOfLE (Set.image_preimage_subset hU.fromSpec U.1)).op :=
  Scheme.Hom.app_appIso_inv _ _

/-- ★★★**`ΓSpecIso.inv ≫ 転送 ≫ appIso.inv` は制限である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これが `hcompat` の核である。 -/
theorem gammaSpec_appIso_inv :
    (Scheme.ΓSpecIso (Γ(X, U))).inv
      ≫ (Spec Γ(X, U)).presheaf.map (eqToHom hU.fromSpec_preimage_self).op
      ≫ (hU.fromSpec.appIso (hU.fromSpec ⁻¹ᵁ U)).inv
      = X.presheaf.map (homOfLE (Set.image_preimage_subset hU.fromSpec U.1)).op := by
  rw [← app_appIso_top_inv hU, IsAffineOpen.fromSpec_app_self]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def gammaSpec_appIso_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——ΓSpecIso と appIso の関係)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
