import ABC3.Found.Arakelov.PicUnitInv

/-!
# Arakelov (B1) 第 39 ブロック —— **`(u ⊗ₘ v) ≫ μ` は成分ごとである**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`δ` の adjunct の左辺そのもの

第 33 ブロックで `δ` の adjunct は

    (unit P ⊗ₘ unit Q) ≫ μ

と書けた。★本ブロックはこれを**純テンソルの上で成分ごとに**書き下す:

    ((u ⊗ₘ v) ≫ μ).app (op V) (a ⊗ₜ b) = u.app a ⊗ₜ v.app b

## ★★★実装の順序が要る(2026-08-18 実測)

★★`erw` で一気に書き換えると **`μ` 自身が展開されて**しまい、
`pushforward_μ_app_tmul` が当たらなくなる(あるいは `whnf` でタイムアウトする)。

★★★**順序**:

    erw [comp_apply, tensorHom_app]      -- ★合成と ⊗ₘ を app に落とす
    simp only [tensorHom_tmul]           -- ★★`simp only` なら μ を展開しない
    erw [pushforward_μ_app_tmul]         -- ★第 32 ブロック
    rfl

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tensor_mu_app_tmul` | ★★★★**`((u ⊗ₘ v) ≫ μ).app (a ⊗ₜ b) = u.app a ⊗ₜ v.app b`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`(u ⊗ₘ v) ≫ μ` は純テンソルの上で成分ごとである**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `δ` の adjunct(第 33 ブロック)を要素の言葉にする最後の道具である。 -/
theorem tensor_mu_app_tmul {A B : X.PresheafOfModules} {P Q : Y.PresheafOfModules}
    (u : P ⟶ (PresheafOfModules.pushforward (pullbackPhi f)).obj A)
    (v : Q ⟶ (PresheafOfModules.pushforward (pullbackPhi f)).obj B)
    (V : Y.Opens) (a : P.obj (op V)) (b : Q.obj (op V)) :
    ((u ⊗ₘ v) ≫ Functor.LaxMonoidal.μ (PresheafOfModules.pushforward (pullbackPhi f)) A B).app
        (op V) (a ⊗ₜ b)
      = (show ((PresheafOfModules.pushforward (pullbackPhi f)).obj (A ⊗ B)).obj (op V) from
          u.app (op V) a ⊗ₜ v.app (op V) b) := by
  erw [CategoryTheory.comp_apply, PresheafOfModules.Monoidal.tensorHom_app]
  erw [pushforward_μ_app_tmul]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def tensor_mu_app_tmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——(u ⊗ₘ v) ≫ μ が成分ごとであること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
