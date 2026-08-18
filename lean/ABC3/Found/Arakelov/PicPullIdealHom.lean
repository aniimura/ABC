import ABC3.Found.Arakelov.PicSectionPull
import ABC3.Found.Arakelov.PicSchemeDelta

/-!
# Arakelov (B2) 第 205 ブロック —— **引き戻しの比較射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`ofDivisor_pullback` の第 3 歩

    (pullbackPre f).obj (idealPresheaf D)  ⟶  idealPresheaf (D.comap f)

★これを**随伴**で作る。`f^*` の具体形(テンソル積)は mathlib に無い(§9-221 で実測)ので、
**押し出し側**で射を書いてから随伴で移す。

## ★★押し出し側は「切断を引き戻すだけ」

    idealPresheaf D  ⟶  f_* (idealPresheaf (D.comap f))
    s ↦ f^# s

★第 204(切断は引き戻せる)がちょうどこの写像の**行き先の証明**である。
加法性・線型性・自然性はすべて `f.app` のそれに帰着する
(`map_add` / `map_mul` / `f.naturality`)。

★★`PresheafOfModules.homMk`(第 178 で見つけた逃げ道)を使うので
`Module` の綴りが 2 通りになる問題を踏まない。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `pushHom` | ★★押し出し側の射(切断の引き戻し) |
| `pullIdealHom` | ★★★**引き戻しの比較射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★**押し出し側の射**——切断を引き戻すだけ。 -/
noncomputable def pushHom :
    idealPresheaf D ⟶ (PresheafOfModules.pushforward (pullbackPhi f)).obj
      (idealPresheaf (D.comap f)) :=
  PresheafOfModules.homMk
    { app := fun V => AddCommGrpCat.ofHom
        (AddMonoidHom.mk' (fun s : (idealSections D V.unop) =>
          (⟨(f.app V.unop).hom s.1,
            app_mem_idealSections f D V.unop s.1 s.2⟩ :
              (idealSections (D.comap f) (f ⁻¹ᵁ V.unop)))) (fun _ _ => Subtype.ext (map_add _ _ _)))
      naturality := by
        intro V W g
        ext s
        exact Subtype.ext (congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) s.1)
          (f.naturality g)) }
    (fun V c s => Subtype.ext (map_mul (f.app V.unop).hom c s.1))




/-- ★★★随伴で移した前層の射。 -/
noncomputable def pullIdealHom :
    (pullbackPre f).obj (idealPresheaf D) ⟶ idealPresheaf (D.comap f) :=
  ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv
    (idealPresheaf D) (idealPresheaf (D.comap f))).symm (pushHom f D)


/-! ## ★出典の紐付け(`.src`) -/

def pullIdealHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻しの比較射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
