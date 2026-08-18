import ABC3.Found.Arakelov.PicAppTopSplit

/-!
# Arakelov (B2) 第 216 ブロック —— **`fromSpec` の切断レベルの値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★生成元の場合に要る 2 本

第 215 で `appTop` の分解が出たので、あとは `fromSpec` そのものの
切断レベルの姿を押さえればよい。

| 定理 | 内容 |
|---|---|
| `fromSpec_app_val` | ★★`fromSpec.app U` は `ΓSpecIso.inv` に制限を継いだもの |
| `fromSpec_image_top` | ★★★`fromSpec ''ᵁ ⊤ = U` |

★どちらも mathlib の `fromSpec_app_self` / `opensRange_fromSpec` から出る。

## ★★★測って分かった —— 最後の同定は **`eqToHom` の輸送**が重い

    (fromSpec.appIso ⊤).hom = (制限 fromSpec''ᵁ⊤ → U) ≫ ΓSpecIso.inv

を書こうとすると、`appIso_hom'` を展開した先で `simp` が
`Spec.map (X.presheaf.map (eqToHom …))` の形に正規化してしまい、
★**`eqToHom` の輸送を手で追う**必要がある(2026-08-19 実測)。

★★これは §9-232 で測った「転送を外した形は `whnf` で落ちる」と**同じ系統**の摩擦である
——`fromSpec ''ᵁ ⊤` と `U` は**等しいが綴りが違う**。
★★★見積もりを **1 → 2–4 ブロック**に更新する(`eqToHom` の輸送を丁寧に追う分)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `fromSpec_app_val` | ★★元のレベルの値 |
| `fromSpec_image_top` | ★★★`fromSpec ''ᵁ ⊤ = U` |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} {U : X.Opens} (hU : IsAffineOpen U)

/-- ★★`fromSpec` の `U` での app は `ΓSpecIso` の逆(制限つき)。 -/
theorem fromSpec_app_val (s : (Γ(X, U) : Type u)) :
    (hU.fromSpec.app U).hom s
      = ((Spec Γ(X, U)).presheaf.map (eqToHom hU.fromSpec_preimage_self).op).hom
        ((Scheme.ΓSpecIso (Γ(X, U))).inv.hom s) := by
  rw [IsAffineOpen.fromSpec_app_self]
  rfl




/-- ★★★`fromSpec ''ᵁ ⊤ = U`(第 200 で使った等式の再掲)。 -/
theorem fromSpec_image_top : hU.fromSpec ''ᵁ (⊤ : (Spec Γ(X, U)).Opens) = U := by
  simp [Scheme.Hom.image_top_eq_opensRange, hU.opensRange_fromSpec]

/-! ## ★出典の紐付け(`.src`) -/

def fromSpec_image_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——fromSpec の切断レベルの値)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
