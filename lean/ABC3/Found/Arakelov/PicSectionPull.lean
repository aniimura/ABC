import ABC3.Found.Arakelov.PicSubschemeVanish

/-!
# Arakelov (B2) 第 204 ブロック —— **切断は引き戻せる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`ofDivisor_pullback` の第 2 歩

    s ∈ idealSections D V  ⟹  f^# s ∈ idealSections (D.comap f) (f⁻¹V)

★これで前層の射 `(pullbackPre f).obj (idealPresheaf D) ⟶ idealPresheaf (D.comap f)` を
随伴で作る材料が揃う。

## ★★機構は引き戻しの可換正方形 1 枚

`comap = (pullback.fst f ι).ker` なので、示すべきは
`(pullback.fst).app A (f^# s |_A) = 0` である。★`appLE` の合成則で左辺は
`(pullback.fst ≫ f).appLE V W` に等しく、可換正方形 `fst ≫ f = snd ≫ ι` で
`(pullback.snd ≫ ι).appLE V W` に移る。これは `ι.app V s` を経由するので
**第 203** で `0` である。

★★`appLE` の合成則(`appLE_comp_appLE` / `comp_appLE`)は mathlib に在り、
`app_eq_appLE` / `appLE_map` と合わせて**制限つきの射**を綺麗に扱える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `app_mem_idealSections_comap` | ★★★アフィン開ごとの所属 |
| `app_mem_idealSections` | ★★★★**切断は引き戻せる** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★★**引き戻した切断はアフィン開ごとに `D.comap f` のイデアルに入る**。 -/
theorem app_mem_idealSections_comap (V : Y.Opens) (s : (Γ(Y, V) : Type u))
    (hs : s ∈ idealSections D V) (A : X.affineOpens) (hA : A.1 ≤ f ⁻¹ᵁ V) :
    (X.presheaf.map (homOfLE hA).op).hom ((f.app V).hom s)
      ∈ (D.comap f).ideal A := by
  set g := Limits.pullback.fst f D.subschemeι with hg
  set h := Limits.pullback.snd f D.subschemeι with hh
  have hsq : g ≫ f = h ≫ D.subschemeι := Limits.pullback.condition
  set W : (Limits.pullback f D.subschemeι).Opens := g ⁻¹ᵁ A.1 with hW
  have hgf : W ≤ (g ≫ f) ⁻¹ᵁ V := by
    intro z hz
    exact hA hz
  have hhi : W ≤ (h ≫ D.subschemeι) ⁻¹ᵁ V := by rw [← hsq]; exact hgf
  -- 左辺は (g ≫ f).appLE
  have hleft : (g.app A.1).hom ((X.presheaf.map (homOfLE hA).op).hom ((f.app V).hom s))
      = ((g ≫ f).appLE V W hgf).hom s := by
    show ((f.app V ≫ X.presheaf.map (homOfLE hA).op) ≫ g.app A.1).hom s = _
    congr 1
    rw [Scheme.Hom.app_eq_appLE, Scheme.Hom.appLE_map, Scheme.Hom.app_eq_appLE,
      Scheme.Hom.appLE_comp_appLE]
  rw [Scheme.IdealSheafData.comap, Scheme.Hom.ker_apply _ A, RingHom.mem_ker, hleft]
  have heq : (g ≫ f).appLE V W hgf = (h ≫ D.subschemeι).appLE V W hhi := by
    congr 1
  rw [heq, Scheme.Hom.comp_appLE]
  show (h.appLE _ W _).hom ((D.subschemeι.app V).hom s) = 0
  rw [subschemeι_app_eq_zero D V s hs, map_zero]


/-- ★★★★**切断は引き戻せる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが前層の射を作る材料である。 -/
theorem app_mem_idealSections (V : Y.Opens) (s : (Γ(Y, V) : Type u))
    (hs : s ∈ idealSections D V) :
    (f.app V).hom s ∈ idealSections (D.comap f) (f ⁻¹ᵁ V) :=
  fun A hA => app_mem_idealSections_comap f D V s hs A hA

/-! ## ★出典の紐付け(`.src`) -/

def app_mem_idealSections.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断は引き戻せる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
