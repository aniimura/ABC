import ABC3.Found.Arakelov.PicResFree

/-!
# Arakelov (B1) 第 48 ブロック —— **制限した site でも生成元の引き戻しは生成元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 24 ブロックがそのまま移る

第 24 ブロック(`f^*(free (yoneda V)) ≅ free (yoneda f⁻¹V)`)の証明は
**`φ` に何も仮定していない**——`CorepresentableBy` の一意性だけである。

★★したがって**制限した site**(`Over V → Over (f⁻¹V)`)にもそのまま移る:

    (f|)^*_pre (free (yoneda Z)) ≅ free (yoneda (overPost Z))

★★★これで Beck–Chevalley の mate の**両端の対象**が計算できる。

## ★★対象の側は一致する(第 23 ブロックによる)

    (f|)^*_pre ((free (yoneda W))|_V) ≅ free (yoneda (overPost (W ⊓ V)))   ★第 47 + 本ブロック
    (f^*_pre (free (yoneda W)))|_{f⁻¹V} ≅ free (yoneda (f⁻¹W ⊓ f⁻¹V))      ★第 24 + 第 47

★`f⁻¹(W ⊓ V) = f⁻¹W ⊓ f⁻¹V`(第 23 `opensMap_inf`)なので**同じ対象**である。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `phiOn` | ★制限した site の環前層の射 |
| `pullbackOnCorepresentableBy` | ★★随伴からの余表現可能性(制限版) |
| `pullbackOnFreeYonedaIso` | ★★★★**制限した引き戻しの生成元上の具体形** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (V : Y.Opens)

/-- ★制限した site の環前層の射(第 20 ブロックの `alphaR` を当てたもの)。 -/
noncomputable abbrev phiOn :=
  alphaR (overPost f V)
    ((Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf)
    ((Over.forget V).op ⋙ Y.presheaf)
    (restrictedC f V)

/-- ★★**随伴からの余表現可能性**(制限版)——第 24 ブロックと同じ 4 行。 -/
noncomputable def pullbackOnCorepresentableBy (M : PresheafModulesOn Y V) :
    (PresheafOfModules.pushforward (phiOn f V) ⋙ coyoneda.obj (op M)).CorepresentableBy
      ((PresheafOfModules.pullback (phiOn f V)).obj M) where
  homEquiv := (PresheafOfModules.pullbackPushforwardAdjunction (phiOn f V)).homEquiv _ _
  homEquiv_comp g h :=
    (PresheafOfModules.pullbackPushforwardAdjunction (phiOn f V)).homEquiv_naturality_right h g

/-- ★★★★★**制限した引き戻しも生成元を生成元に送る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 24 ブロックの証明は `φ` に何も仮定していないので、
**制限した site にもそのまま移る**。 -/
noncomputable def pullbackOnFreeYonedaIso (Z : Over V) :
    (pullbackPreOn f V).obj
        ((PresheafOfModules.free
          (((Over.forget V).op ⋙ Y.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
          (yoneda.obj Z))
      ≅ (PresheafOfModules.free
          (((Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf)
            ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
          (yoneda.obj ((overPost f V).obj Z)) :=
  (pullbackOnCorepresentableBy f V _).uniqueUpToIso
    (PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy (phiOn f V) Z)

/-! ## ★出典の紐付け(`.src`) -/

def pullbackOnFreeYonedaIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限した引き戻しも生成元を生成元に送ること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
