import ABC3.Found.Arakelov.PicFreeDescBij

/-!
# Arakelov (B1) 第 107 ブロック —— **単位からの射の全単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`Inhabited` の二経路が正体だった

§9-110〜§9-119 で **20 手以上**詰まった所である。★正体は

    @default (W ⟶ T) (@Unique.instInhabited _ _)      ← ゴール側
    @default (W ⟶ T) (@Unique.toInhabited _ _)        ← 補助定義側

という **`Unique` から `Inhabited` を取る道の二重性**だった。

★★★**解決**: 補助定義 `restrictSec` を
**`[Unique ((yoneda.obj T).obj (op W))]` で引数化**する
——ゴール側(第 106 ブロック)と**同じ class** なので `default` が一致する。

★★★★これは [[ring-instance-two-paths]] の **3 例目**である
(環・加群に続いて `Inhabited`)。

## ★★証明は 4 行になった

    rw [freeYonedaEquiv_symm_eq_desc, freeObjDesc_app]
    refine bijective_freeDesc_of_unique _ ?_
    convert hb using 2 with c
    congr 1

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `restrictSec` | ★型を固定した切断 |
| `bijective_freeYonedaEquiv_symm_app` | ★★★★★★**単位からの射の全単射性** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (V : X.Opens) (P : PresheafModulesOn X V)

/-- ★**型を固定した切断**——`P.map d.op s` を `P.obj (op W)` の元と見る。

★★`[Unique ((yoneda.obj T).obj (op W))]` で引数化するのが要点である
——これでゴール側と `default` が一致する。 -/
noncomputable def restrictSec (s : P.obj (op (Over.mk (𝟙 V)))) (W : Over V)
    [Unique ((yoneda.obj (Over.mk (𝟙 V))).obj (op W))] : (P.obj (op W) : Type u) :=
  P.map (default : (yoneda.obj (Over.mk (𝟙 V))).obj (op W)).op s

/-- ★★★★★★**単位からの射は、生成元の乗法が全単射なら全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで局所自明性の同型が「生成元の乗法」に帰着する。 -/
theorem bijective_freeYonedaEquiv_symm_app (s : P.obj (op (Over.mk (𝟙 V)))) (W : Over V)
    [Unique ((yoneda.obj (Over.mk (𝟙 V))).obj (op W))]
    (hb : Function.Bijective (fun c : ((((Over.forget V).op ⋙ X.presheaf)
        ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op W) : Type u) =>
      c • restrictSec V P s W)) :
    Function.Bijective ((PresheafOfModules.freeYonedaEquiv.symm s).app (op W)) := by
  rw [freeYonedaEquiv_symm_eq_desc, PresheafOfModules.freeObjDesc_app]
  refine bijective_freeDesc_of_unique _ ?_
  convert hb using 2 with c
  congr 1

/-! ## ★出典の紐付け(`.src`) -/

def bijective_freeYonedaEquiv_symm_app.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——単位からの射の全単射性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
