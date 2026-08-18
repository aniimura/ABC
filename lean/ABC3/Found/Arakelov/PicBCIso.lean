import ABC3.Found.Arakelov.PicBCMateIso
import ABC3.Found.Arakelov.PicFreeYonedaLift

/-!
# Arakelov (B1) 第 54 ブロック —— ★★★★★★★★★★**Beck–Chevalley が同型になった**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★取れたもの

    (f^*_pre L)|_{f⁻¹V}  ≅  (f|)^*_pre (L|_V)      すべての `L` について

★★★これが §9-52 で「**我々が作るしかない**」と測った Beck–Chevalley である。
**mathlib には無い**(2026-08-18 実測)。

## ★★★★道のり(第 43–54、12 ブロック)

| 段 | ブロック | 内容 |
|---|---|---|
| 枠 | 43–47 | 余極限保存・制限した引き戻し・四角形・mate・生成元の制限 |
| 一般化 | 48–52 | 制限した site の生成元・`unit` の一般化・制限同型の値 |
| 完成 | 53–54 | mate が生成元で同型 → **全対象へ持ち上げ** |

## ★★★★★器具が 2 度目の出番を得た

★第 28・29 ブロック(余極限で同型を持ち上げる)は `δ` のために作った。
★★**Beck–Chevalley でもそのまま効いた**——1 行で持ち上がる。

## ★★これで第 41 ブロックの仮定が落とせる

`L` が局所自明なら、`L|_V ≅ 𝒪_V` を持つ被覆があり、

    (f^*_pre L)|_{f⁻¹V} ≅ (f|)^*_pre (L|_V) ≅ (f|)^*_pre 𝒪_V ≅ 𝒪_{f⁻¹V}

(最後は第 42 ブロック)。★★★すなわち **`f^*_pre L` は局所自明**である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (V : Y.Opens)

/-! ## ★★★★★★生成元の上で同型 -/

/-- ★★★★★★**mate は生成元の上で同型である**。

★第 53 ブロックの等式の右辺は同型の合成であり、左端の
`(pullbackPreOn).map (iso.inv)` も同型だから。 -/
theorem isIso_bcMate_free (W : Y.Opens) :
    IsIso ((bcMate f V).app (freeY W)) := by
  haveI : IsIso ((pullbackPreOn f V).map (restrictFreeYonedaIso V W).inv) :=
    (pullbackPreOn f V).map_isIso _
  haveI : IsIso ((pullbackPreOn f V).map (restrictFreeYonedaIso V W).inv
      ≫ (bcMate f V).app (freeY W)) := by
    rw [bcMate_free_eq]
    haveI h1 : IsIso ((restrictPresheafFunctor X ((Opens.map f.base).obj V)).map
        (pullbackFreeYonedaIso f W).inv) :=
      (restrictPresheafFunctor X ((Opens.map f.base).obj V)).map_isIso _
    haveI h2 : IsIso ((pullbackOnFreeYonedaIso f V (objOn V W)).hom) := Iso.isIso_hom _
    haveI h3 : IsIso ((restrictFreeYonedaIso ((Opens.map f.base).obj V)
        ((Opens.map f.base).obj W)).inv) := Iso.isIso_inv _
    exact @IsIso.comp_isIso _ _ _ _ _ _ _ h2 (@IsIso.comp_isIso _ _ _ _ _ _ _ h3 h1)
  exact IsIso.of_isIso_comp_left ((pullbackPreOn f V).map (restrictFreeYonedaIso V W).inv) _

/-! ## ★★★★★★★★★★全対象で同型 -/

/-- ★★★★★★★★★★**Beck–Chevalley は全対象で同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★§9-52 で「**mathlib に無いので我々が作るしかない**」と測ったものである。

★機構は第 29 ブロック(生成元で同型なら全対象で同型)——
`δ` のために作った器具が**そのまま 2 度目の出番を得た**。 -/
theorem isIso_bcMate (P : Y.PresheafOfModules) :
    IsIso ((bcMate f V).app P) := by
  haveI hpull : PreservesColimitsOfSize.{u, u} (pullbackPre f) := pullbackPre_preservesColimits f
  exact isIso_app_of_freeYoneda (bcMate f V)
    (preservesColimitsOfSize_comp (restrictPresheafFunctor Y V) (pullbackPreOn f V)
      (restrictPresheafFunctor_preservesColimits Y V)
      ((PresheafOfModules.pullbackPushforwardAdjunction
        (phiOn f V)).leftAdjoint_preservesColimits))
    (preservesColimitsOfSize_comp (pullbackPre f)
      (restrictPresheafFunctor X ((Opens.map f.base).obj V)) hpull
      (restrictPresheafFunctor_preservesColimits X ((Opens.map f.base).obj V)))
    (fun W => isIso_bcMate_free f V W) P

/-- ★★★★★★★★★★**Beck–Chevalley の同型** `(f|)^*_pre (L|_V) ≅ (f^*_pre L)|_{f⁻¹V}`。 -/
noncomputable def bcIso (P : Y.PresheafOfModules) :
    (pullbackPreOn f V).obj ((restrictPresheafFunctor Y V).obj P)
      ≅ (restrictPresheafFunctor X ((Opens.map f.base).obj V)).obj ((pullbackPre f).obj P) := by
  exact @asIso _ _ _ _ ((bcMate f V).app P) (isIso_bcMate f V P)

/-! ## ★出典の紐付け(`.src`) -/

def isIso_bcMate.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Beck–Chevalley が同型であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
