import ABC3.Found.Arakelov.PicPushAdj
import ABC3.Found.Arakelov.PicIdealMulHom
import ABC3.Found.Arakelov.PicCartierComap
import ABC3.Found.Arakelov.PicSurjIso

/-!
# Arakelov (B2) 第 209 ブロック —— **比較射の像**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-225 の補題 B の第 2 歩

比較射の像を扱う道具を揃える:

| 事実 | 根拠 |
|---|---|
| 像は引き戻した切断を含む | ★第 208 |
| 像は制限で運べる | ★自然性 |
| 像は掛け算で閉じる | ★線型性 |

★★★これで「像 ⊇ 生成元 かつ 像は部分加群 ⟹ 像 = 全体」が言える。
**`f^*` の具体形は最後まで使わない**。

## ★★`Submodule` でなく `Set` で扱う

`LinearMap.range` を使おうとすると `RingHomSurjective (RingHom.id …)` の
instance 解決に失敗する(係数環の綴りが `X.ringCatSheaf.obj.obj` 経由になるため)。
★`Set.range` と `Set.image` で扱い、必要な閉性だけを**手で**示すほうが短い
——[[ring-instance-two-paths]] の逃げ道のひとつである。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `imgSet` | ★比較射の像 |
| `mem_imgSet_of_pushHom` | ★★引き戻した切断は像に入る |
| `imgSet_res` | ★★像は制限で運べる |
| `imgIdeal` / `imgIdeal_isSubmodule` | ★★★像は掛け算で閉じる |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★比較射の像(集合として)。 -/
def imgSet (A : X.Opens) : Set ((idealPresheaf (D.comap f)).obj (op A) : Type u) :=
  Set.range ((pullIdealHom f D).app (op A)).hom

theorem mem_imgSet_of_pushHom (B : Y.Opens) (s : ((idealPresheaf D).obj (op B) : Type u)) :
    ((pushHom f D).app (op B)).hom s ∈ imgSet f D (f ⁻¹ᵁ B) :=
  ⟨_, (pushHom_app f D B s).symm⟩

/-- ★像は制限で運べる(自然性)。 -/
theorem imgSet_res {A A' : X.Opens} (h : A' ≤ A)
    (x : ((idealPresheaf (D.comap f)).obj (op A) : Type u)) (hx : x ∈ imgSet f D A) :
    ((idealPresheaf (D.comap f)).map (homOfLE h).op).hom x ∈ imgSet f D A' := by
  obtain ⟨y, rfl⟩ := hx
  refine ⟨(((pullbackPre f).obj (idealPresheaf D)).map (homOfLE h).op).hom y, ?_⟩
  exact congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) y)
    ((pullIdealHom f D).naturality (homOfLE h).op)




/-- ★像を `Γ(X,A)` の部分集合として見る(切断の値)。 -/
def imgIdeal (A : X.Opens) : Set (Γ(X, A) : Type u) :=
  sVal '' (imgSet f D A)

theorem imgIdeal_isSubmodule (A : X.Opens) :
    ∀ (c : (Γ(X, A) : Type u)) {a : (Γ(X, A) : Type u)}, a ∈ imgIdeal f D A →
      c * a ∈ imgIdeal f D A := by
  rintro c a ⟨x, ⟨y, rfl⟩, rfl⟩
  refine ⟨c • ((pullIdealHom f D).app (op A)).hom y, ⟨c • y, ?_⟩, ?_⟩
  · exact map_smul _ _ _
  · rfl

/-! ## ★出典の紐付け(`.src`) -/

def imgIdeal_isSubmodule.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——比較射の像)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
