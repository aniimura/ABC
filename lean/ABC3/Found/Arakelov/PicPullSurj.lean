import ABC3.Found.Arakelov.PicAppIsoVal
import ABC3.Found.Arakelov.PicAppLEApply
import ABC3.Found.Arakelov.PicSectionPull
import ABC3.Found.Arakelov.PicPullImage
import ABC3.Found.Arakelov.PicImgClosure

/-!
# Arakelov (B2) 第 218 ブロック —— ★★★★★★**比較射はアフィン開で全射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★§9-225 の補題 B が出た

    (D.comap f).ideal A' = (D.ideal B').map appLE
      ⟹ pullIdealHom.app A' は全射

★★★**`f^*` の具体形を一度も使っていない**——§9-221 で「具体形は塞がっている」と
測った所を、**像が部分加群であることだけ**で迂回した。

## ★★機構は `span_induction` 1 回

| 場合 | 根拠 |
|---|---|
| 生成元 | ★第 208(引き戻した切断は像に入る)+ 自然性 |
| 零・加法 | ★第 212 |
| スカラー倍 | ★第 209 |

★`Ideal.map` は生成元の span なので、この 4 つで全体が覆える。

## ★★イデアルの等式は**仮定**として受ける

第 210–217 で部品は揃ったが、`fromSpec` の転送を通した形なので、
★本ブロックは**等式を仮定に取る**——後で `A' = fromSpec ''ᵁ ⊤` の場合に当てる。
★★こうすると本ブロックが**転送に依存しない**再利用可能な形になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gen_in_img` | ★★★生成元の像は `imgIdeal` に入る |
| `surjective_pullIdealHom` | ★★★★★★**比較射はアフィン開で全射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★★**生成元の像は `imgIdeal` に入る**。 -/
theorem gen_in_img {B' : Y.affineOpens} {A' : X.affineOpens} (hle : A'.1 ≤ f ⁻¹ᵁ B'.1)
    (t : (Γ(Y, B'.1) : Type u)) (ht : t ∈ D.ideal B') :
    (X.presheaf.map (homOfLE hle).op).hom ((f.app B'.1).hom t) ∈ imgIdeal f D A'.1 := by
  have hts : t ∈ idealSections D B'.1 := (idealSections_affine D B') ▸ ht
  have h1 := mem_imgSet_of_pushHom f D B'.1 (⟨t, hts⟩ : ((idealPresheaf D).obj (op B'.1) : Type u))
  have h2 := imgSet_res f D hle _ h1
  exact ⟨_, h2, rfl⟩


/-- ★★★★★★**比較射はアフィン開で全射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`f^*` の具体形を一度も使わず、**像が部分加群であること**だけで出る。 -/
theorem surjective_pullIdealHom {B' : Y.affineOpens} {A' : X.affineOpens}
    (hle : A'.1 ≤ f ⁻¹ᵁ B'.1)
    (hideal : (D.comap f).ideal A' = (D.ideal B').map (f.appLE B'.1 A'.1 hle).hom) :
    Function.Surjective (((pullIdealHom f D).app (op A'.1)).hom) := by
  intro z
  have hz : sVal z ∈ (D.comap f).ideal A' := (idealSections_affine (D.comap f) A') ▸ z.2
  rw [hideal] at hz
  have hgen : ∀ w ∈ (D.ideal B').map (f.appLE B'.1 A'.1 hle).hom, w ∈ imgIdeal f D A'.1 := by
    intro w hw
    refine Submodule.span_induction ?_ (imgIdeal_zero f D A'.1)
      (fun _ _ _ _ ha hb => imgIdeal_add f D A'.1 ha hb)
      (fun c _ _ ha => by rw [smul_eq_mul]; exact imgIdeal_isSubmodule f D A'.1 c ha) hw
    rintro _ ⟨t, ht, rfl⟩
    exact gen_in_img f D hle t ht
  obtain ⟨x, ⟨y, rfl⟩, hval⟩ := hgen _ hz
  exact ⟨y, Subtype.ext hval⟩


/-! ## ★出典の紐付け(`.src`) -/

def surjective_pullIdealHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——比較射はアフィン開で全射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
