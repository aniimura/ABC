import ABC3.Found.Arakelov.PicGammaStruct

/-!
# Arakelov (B1) 第 71 ブロック —— **局所自明性から切断の階数 1 自由性へ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`Interface` の `IsInvertibleSheaf` に橋を架ける

`Interface/Arakelov/LineBundle.lean` の `IsInvertibleSheaf` は
**「各開集合が、階数 1 自由な開集合で覆われる」**を

    ∀ U, ∃ S 被覆篩, ∀ V i, S i → Nonempty (𝒪(V) ≃ₗ[𝒪(V)] F.val(V))

と `≃ₗ` で述べている(`Interface` は `Found` を import できないので
mathlib の語彙だけで書いてある)。

★★一方 `Found` 側の `IsLocallyTrivial` は**制限が単位に同型**と述べている:

    ∀ U, ∃ S 被覆篩, ∀ V i, S i → Nonempty ((restrict V).obj M ≅ 𝟙_)

★★★**この 2 つを繋ぐ**のが本ブロックである。

## ★★機構 —— `Over V` の終対象で評価する

`(restrict V).obj M ≅ 𝟙_ (PresheafModulesOn X V)` を
`Over V` の終対象 `Over.mk (𝟙 V)` で評価すると

    M.obj (op V) ≅ 𝒪(V)          (`𝒪(V)` 加群として)

★`PresheafOfModules` には `evaluation` 関手が**無い**(2026-08-18 実測)ので、
`hom`/`inv` を手で置き、両端の等式は `congrArg` で運ぶ。

★★`rw [← PresheafOfModules.comp_app]` は**通らない**——環インスタンスの
二経路([[ring-instance-two-paths]])に当たる。`congrArg` なら当たらない。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `sectionIsoOfTrivial` | ★★終対象で評価した切断の同型 |
| `isRankOneFree_of_isLocallyTrivial` | ★★★★**`Interface` の階数 1 条件** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X : Scheme.{u}}

/-- ★★**制限が単位に同型なら、その開集合での切断は `𝒪(V)` に同型**。

★`Over V` の終対象 `Over.mk (𝟙 V)` で評価するだけである。 -/
noncomputable def sectionIsoOfTrivial (V : X.Opens) (M : X.PresheafOfModules)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    M.obj (op V) ≅ ModuleCat.of (X.presheaf.obj (op V) : Type u)
      (X.presheaf.obj (op V) : Type u) where
  hom := e.hom.app (op (Over.mk (𝟙 V)))
  inv := e.inv.app (op (Over.mk (𝟙 V)))
  hom_inv_id :=
    congrArg (fun f : (restrictPresheafFunctor X V).obj M ⟶ _ =>
      f.app (op (Over.mk (𝟙 V)))) e.hom_inv_id
  inv_hom_id :=
    congrArg (fun f : 𝟙_ (PresheafModulesOn X V) ⟶ _ =>
      f.app (op (Over.mk (𝟙 V)))) e.inv_hom_id

/-- ★★★★**局所自明な層加群は `Interface` の階数 1 自由性を満たす**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `Found` の `InvSheaf` を `Interface` の `IsInvertibleSheaf` に
翻訳する段が繋がる。 -/
theorem isRankOneFree_of_isLocallyTrivial {M : X.PresheafOfModules}
    (h : IsLocallyTrivial X M) (U : X.Opens) :
    ∃ S : Sieve U, S ∈ (Opens.grothendieckTopology X) U ∧
      ∀ ⦃V : X.Opens⦄ (i : V ⟶ U), S.arrows i →
        Nonempty ((X.presheaf.obj (op V) : Type u)
          ≃ₗ[(X.presheaf.obj (op V) : Type u)] (M.obj (op V) : Type u)) := by
  obtain ⟨S, hS, hV⟩ := h U
  refine ⟨S, hS, fun V i hi => ?_⟩
  exact (hV i hi).map fun e =>
    (sectionIsoOfTrivial V M e).symm.toLinearEquiv

/-! ## ★出典の紐付け(`.src`) -/

def isRankOneFree_of_isLocallyTrivial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所自明性から階数 1 自由性へ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
