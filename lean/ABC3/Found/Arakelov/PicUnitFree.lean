import ABC3.Found.Arakelov.PicDeltaLift

/-!
# Arakelov (B1) 第 31 ブロック —— **随伴の unit を生成元の上で書き下す**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★抽象な随伴の unit に、計算の手がかりを付ける

`PresheafOfModules.pullback` は随伴関手定理で作られているので、
`unit` も `counit` も**具体的な記述を持たない**。
★★しかし `δ` の定義は

    δ P Q = homEquiv.symm ((unit P ⊗ₘ unit Q) ≫ μ)

なので、`unit` が計算できなければ `δ` は計算できない。

## ★★★生成元の上でだけは計算できる

第 24 ブロックの `pullbackFreeYonedaIso` は
`CorepresentableBy.uniqueUpToIso` で作った。その `inv` は定義により

    (mathlib 側の homEquiv).symm (随伴の homEquiv (𝟙))

であり、`随伴の homEquiv (𝟙) = unit` である。★★★したがって

    unit (free (yoneda V)) = (mathlib 側の homEquiv) ((pullbackFreeYonedaIso f V).inv)

と**書き下せる**。mathlib 側の `homEquiv` は `freeYonedaEquiv` で明示的である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pullbackFreeYonedaIso_inv` | ★`inv` は 2 つの `homEquiv` の合成(`rfl`) |
| `pullbackUnit_app_free` | ★★★★**unit を生成元の上で書き下したもの** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★生成元(`free (yoneda V)`)の略記。 -/
noncomputable abbrev freeY (V : Y.Opens) : Y.PresheafOfModules :=
  (PresheafOfModules.free (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj (yoneda.obj V)

/-- ★mathlib 側の余表現(`pushforward ⋙ coyoneda` を `free (yoneda (f⁻¹V))` が余表現する)。 -/
noncomputable abbrev mathlibCorep (V : Y.Opens) :=
  PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy (pullbackPhi f) V

/-! ## ★同型の逆射は 2 つの `homEquiv` の合成である -/

/-- ★**`CorepresentableBy.uniqueUpToIso` の定義そのもの**(`rfl` で出る)。 -/
theorem pullbackFreeYonedaIso_inv (V : Y.Opens) :
    (pullbackFreeYonedaIso f V).inv
      = (mathlibCorep f V).homEquiv.symm
        ((pullbackCorepresentableBy f (freeY V)).homEquiv (𝟙 _)) := rfl

/-! ## ★★★★unit を生成元の上で書き下す -/

/-- ★★★★★**随伴の unit を生成元の上で書き下したもの**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `δ` を計算するための**唯一の手がかり**である——
`pullback` は随伴関手定理で作られているので、`unit` に他の記述は無い。

★機構は 3 段:
1. `pullbackCorepresentableBy` の `homEquiv` は随伴の `homEquiv` そのもの
2. 随伴の `homEquiv (𝟙)` は `unit`(`Adjunction.homEquiv_unit` で `𝟙` を入れる)
3. `uniqueUpToIso` の `inv` は 2 つの `homEquiv` の合成(上の `rfl`) -/
theorem pullbackUnit_app_free (V : Y.Opens) :
    (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app (freeY V)
      = (mathlibCorep f V).homEquiv ((pullbackFreeYonedaIso f V).inv) := by
  rw [pullbackFreeYonedaIso_inv, Equiv.apply_symm_apply]
  show _ = (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv _ _ (𝟙 _)
  rw [Adjunction.homEquiv_unit]
  simp only [Functor.id_obj, CategoryTheory.Functor.map_id]
  exact (Category.comp_id _).symm

/-! ## ★出典の紐付け(`.src`) -/

def pullbackUnit_app_free.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——随伴の unit を生成元の上で書き下すこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
