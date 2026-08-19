import ABC3.Found.Arakelov.PicTrivRestrict

/-!
# Arakelov (B2) 第 238 ブロック —— ★★★★★**良い開集合の上で比較射は全単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★単射性は「全射性 + 直線束」から出る

第 235 ブロックで**全射性**が取れた。★単射性を別に作る必要は無い——
第 207 ブロック(自明化開の上では全射な射は同型)を噛ませればよい。

    V はアフィン ∧ V ≤ f⁻¹ᵁ B ∧ V で両辺が自明
      ⟹ 第 235 で全射 ⟹ 第 207 で同型 ⟹ 全単射

★★これが「良い開集合」の定義である。

## ★摩擦 #7 がもう一度出た

第 235 の全射性は座標 `hA.fromSpec ''ᵁ ⊤` で述べてある。
★`A` と `hA.fromSpec ''ᵁ ⊤` は**等しいが綴りが違う**ので、
`image_top_eq_opensRange` と `opensRange_fromSpec` で橋を架ける(`img_top`)。
★★今回は `Opens` の等式を `rw` するだけで通った——`appLE` の側に
依存位置を追い出してあったからである。

| 定理 | 内容 |
|---|---|
| `img_top` | ★`hA.fromSpec ''ᵁ ⊤ = A` |
| `bij_of_good` | ★★★★★**良い開集合の上で比較射は全単射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

/-- ★**`fromSpec` の像は元のアフィン開である**。 -/
theorem img_top {Z : Scheme.{u}} {A : Z.Opens} (hA : IsAffineOpen A) :
    hA.fromSpec ''ᵁ (⊤ : (Spec Γ(Z, A)).Opens) = A := by
  rw [Scheme.Hom.image_top_eq_opensRange, hA.opensRange_fromSpec]

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★★★★**良い開集合の上で比較射は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★全射性は第 235、単射性は第 207(全射 ⟹ 同型)から出る。 -/
theorem bij_of_good {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (eL : (restrictPresheafFunctor X A).obj ((pullbackPre f).obj (idealPresheaf D))
        ≅ 𝟙_ (PresheafModulesOn X A))
    (eM : (restrictPresheafFunctor X A).obj (idealPresheaf (D.comap f))
        ≅ 𝟙_ (PresheafModulesOn X A)) :
    Function.Bijective (((pullIdealHom f D).app (op A)).hom) := by
  have hle : hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens)
      ≤ f ⁻¹ᵁ (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens)) := by
    rw [img_top hA, img_top hB]; exact i
  have hs := surjective_affine f D hA hB i hle
  rw [img_top hA] at hs
  haveI := isIso_restrict_of_surjective (pullIdealHom f D) A eL eM hs
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map
      ((restrictPresheafFunctor X A).map (pullIdealHom f D))) := inferInstance
  rw [NatTrans.isIso_iff_isIso_app] at this
  haveI := this (op (Over.mk (𝟙 A)))
  have hb : Function.Bijective (((PresheafOfModules.toPresheaf _).map
      ((restrictPresheafFunctor X A).map (pullIdealHom f D))).app (op (Over.mk (𝟙 A)))) :=
    ConcreteCategory.bijective_of_isIso _
  exact hb

/-! ## ★出典の紐付け(`.src`) -/

def bij_of_good.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——良い開集合の上で比較射は全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
