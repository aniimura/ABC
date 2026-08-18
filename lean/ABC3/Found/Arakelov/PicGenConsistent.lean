import ABC3.Found.Arakelov.PicUnitGen

/-!
# Arakelov (B1) 第 51 ブロック —— **一般形は両方の site に効く**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★リファクタが正しいことを機械で確かめる

第 50 ブロックで `unit` の計算を `φ` について一般化した。
★★**一般化が元のものと一致すること**を確かめる:

    pullbackFreeIsoGen (pullbackPhi f) W = pullbackFreeYonedaIso f W    (`rfl`)

★★★これが `rfl` で出る以上、**第 24–40 ブロックの結果はそのまま保たれる**。

## ★★そして制限した site にも効く

    pullbackFreeIsoGen (phiOn f V) Z : (f|)^*_pre (free (yoneda Z)) ≅ free (yoneda (overPost Z))
    isoHomUnitGenGen (phiOn f V) Z   : unit は生成元を生成元に送る

★★★★これで Beck–Chevalley の mate を `δ` と**同じ手口**で計算できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pullbackFreeIsoGen_eq` | ★★一般形はスキームの場合と `rfl` で一致 |
| `pullbackOnFreeIsoGen` | ★制限した site での生成元の引き戻し |
| `isoHomUnitGenOn` | ★★★制限した site でも `unit` は生成元を生成元に送る |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (V : Y.Opens)

/-! ## ★★リファクタの一致 -/

/-- ★★**一般形は第 24 ブロックと `rfl` で一致する**——リファクタが結果を変えていない。 -/
theorem pullbackFreeIsoGen_eq (W : Y.Opens) :
    pullbackFreeIsoGen (pullbackPhi f) W = pullbackFreeYonedaIso f W := rfl

/-! ## ★制限した site での適用 -/

/-- ★**制限した site での生成元の引き戻し**(一般形の適用)。 -/
noncomputable abbrev pullbackOnFreeIsoGen (Z : Over V) :
    (PresheafOfModules.pullback (phiOn f V)).obj
        ((PresheafOfModules.free
          (((Over.forget V).op ⋙ Y.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
          (yoneda.obj Z))
      ≅ (PresheafOfModules.free
          (((Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf)
            ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
          (yoneda.obj ((overPost f V).obj Z)) :=
  pullbackFreeIsoGen (phiOn f V) Z

/-- ★★★**制限した site でも `unit` は生成元を生成元に送る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで Beck–Chevalley の mate を `δ`(第 40 ブロック)と**同じ手口**で計算できる。 -/
theorem isoHomUnitGenOn (Z : Over V) :
    ((pullbackOnFreeIsoGen f V Z).hom).app (op ((overPost f V).obj Z))
        (PresheafOfModules.freeYonedaEquiv
          (M := (PresheafOfModules.pushforward (phiOn f V)).obj
            ((PresheafOfModules.pullback (phiOn f V)).obj
              ((PresheafOfModules.free _).obj (yoneda.obj Z)))) (X := Z)
          ((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f V)).unit.app
            ((PresheafOfModules.free _).obj (yoneda.obj Z))))
      = ModuleCat.freeMk (𝟙 ((overPost f V).obj Z)) :=
  isoHomUnitGenGen (phiOn f V) Z

/-! ## ★出典の紐付け(`.src`) -/

def pullbackFreeIsoGen_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——一般形が元の形と一致すること)",
    sectionId := "genell-def-1-1-i" }

def isoHomUnitGenOn.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限した site でも unit が生成元を生成元に送ること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
