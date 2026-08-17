import ABC3.Found.Arakelov.PicDeltaCalc

/-!
# Arakelov (B1) 第 34 ブロック —— **`unit` を要素の言葉に翻訳する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★最終段に要る 2 本

第 33 ブロックまでで `δ` の adjunct は

    ι ≫ (unit P ⊗ₘ unit Q) ≫ μ

と書けた。★これを「同型の合成」と突き合わせるには、右辺側の adjunct も要る:

| 本ブロック | 内容 |
|---|---|
| `homEquiv_pullbackFreeYonedaIso_hom` | ★★同型の `hom` の adjunct は `unit ≫ f_*(hom)` |
| `freeYonedaEquiv_unit_free` | ★★★`unit` の `freeYonedaEquiv` は `inv` のそれに等しい |

★★★2 本目が効く理由: `inv ≫ hom = 𝟙` なので

    hom.app (freeYonedaEquiv inv) = freeYonedaEquiv (inv ≫ hom) = freeMk (𝟙)

となり、**右辺は生成元そのもの**になる。

## ★★これで最終段は「要素の等式 1 本」である

    左辺: unit_P(freeMk i_V) ⊗ₜ unit_Q(freeMk i_W)   (第 31・32・本ブロック)
    右辺: freeMk (𝟙 U₀) の像                          (本ブロック)

★`U₀ = f⁻¹V ⊓ f⁻¹W = f⁻¹(V ⊓ W)`(第 23 ブロック、`rfl`)。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-! ## ★★同型の `hom` の adjunct -/

/-- ★★**同型の `hom` の adjunct は `unit ≫ f_*(hom)` である**。

★機構は `Adjunction.homEquiv_unit` そのものである。 -/
theorem homEquiv_pullbackFreeYonedaIso_hom (V : Y.Opens) :
    (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv _ _
        ((pullbackFreeYonedaIso f V).hom)
      = (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app (freeY V)
        ≫ (PresheafOfModules.pushforward (pullbackPhi f)).map
          ((pullbackFreeYonedaIso f V).hom) := by
  rw [Adjunction.homEquiv_unit]
  rfl

/-! ## ★★★`unit` を要素にする -/

/-- ★★★★**`unit` の `freeYonedaEquiv` は、第 24 ブロックの同型の `inv` のそれに等しい**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「抽象な `unit`」を**具体的な元**に変える最後の一手である。
★機構は第 31 ブロック(`unit` の書き下し)+ 第 33 ブロック(余表現の `homEquiv`)。

★★これがあると `inv ≫ hom = 𝟙` から

    hom.app (freeYonedaEquiv (unit)) = freeMk (𝟙)

が出る——すなわち **`unit` は生成元を生成元に送る**。 -/
theorem freeYonedaEquiv_unit_free (V : Y.Opens) :
    PresheafOfModules.freeYonedaEquiv
        ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app (freeY V))
      = PresheafOfModules.freeYonedaEquiv
        (M := (PresheafOfModules.pullback (pullbackPhi f)).obj (freeY V))
        (X := (Opens.map f.base).obj V)
        ((pullbackFreeYonedaIso f V).inv) := by
  rw [pullbackUnit_app_free, mathlibCorep_homEquiv]
  exact Equiv.apply_symm_apply _ _

/-! ## ★出典の紐付け(`.src`) -/

def freeYonedaEquiv_unit_free.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unit を要素の言葉に翻訳すること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
