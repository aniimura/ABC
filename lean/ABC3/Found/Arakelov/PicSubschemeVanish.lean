import ABC3.Found.Arakelov.PicCartierComap

/-!
# Arakelov (B2) 第 203 ブロック —— **イデアル層の切断は部分スキームで消える**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★最後の 1 欄(`ofDivisor_pullback`)の第 1 歩

    s ∈ idealSections D V  ⟹  D.subschemeι.app V s = 0

★これがあると、引き戻した切断 `f^# s` が `D.comap f` の切断であることが
**引き戻しの可換正方形**から出る:

    pullback.fst ≫ f = pullback.snd ≫ D.subschemeι
      ⟹ (pullback.fst).app (f^# s) = (pullback.snd).app (ι.app s) = 0

## ★★機構は層の貼り合わせ 1 回

`idealSections D V` の定義は「**アフィン開ごと**に `D.ideal` に入る」であり、
`D.ideal B = ker (ι.app B)`(mathlib `ker_subschemeι_app`)なので、
`ι.app V s` は `ι⁻¹B` たちの上で消える。★`ι⁻¹B` たちは `ι⁻¹V` を覆うので
`D.subscheme` の構造層の**貼り合わせの一意性**から `0` である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `subschemeι_app_eq_zero` | ★★★**イデアル層の切断は部分スキームで消える** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {Y : Scheme.{u}} (D : Y.IdealSheafData)

/-- ★★★**イデアル層の切断は部分スキームの上で消える**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★アフィン開ごとの条件を、`D.subscheme` の構造層の貼り合わせで大域に上げる。 -/
theorem subschemeι_app_eq_zero (V : Y.Opens) (s : (Γ(Y, V) : Type u))
    (hs : s ∈ idealSections D V) :
    (D.subschemeι.app V).hom s = 0 := by
  refine D.subscheme.sheaf.eq_of_locally_eq'
    (fun B : { B : Y.affineOpens // B.1 ≤ V } => D.subschemeι ⁻¹ᵁ B.1.1)
    (D.subschemeι ⁻¹ᵁ V) (fun B => homOfLE (by exact Scheme.Hom.preimage_mono _ B.2))
    ?_ _ _ ?_
  · intro y hy
    obtain ⟨B, hB, hyB, hBV⟩ := Opens.isBasis_iff_nbhd.1 Y.isBasis_affineOpens
      (U := V) (x := D.subschemeι.base y) hy
    exact Opens.mem_iSup.2 ⟨⟨⟨B, hB⟩, hBV⟩, hyB⟩
  · intro B
    have hnat := D.subschemeι.naturality (homOfLE B.2).op
    have h1 : D.subscheme.sheaf.obj.map (homOfLE (Scheme.Hom.preimage_mono _ B.2)).op
        ((D.subschemeι.app V).hom s)
        = (D.subschemeι.app B.1.1).hom ((Y.presheaf.map (homOfLE B.2).op).hom s) :=
      congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) s) hnat.symm
    rw [h1]
    have h2 := hs B.1 B.2
    rw [← Scheme.IdealSheafData.ker_subschemeι_app] at h2
    exact (RingHom.mem_ker.1 h2).trans (map_zero _).symm


/-! ## ★出典の紐付け(`.src`) -/

def subschemeι_app_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——イデアル層の切断は部分スキームで消える)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
