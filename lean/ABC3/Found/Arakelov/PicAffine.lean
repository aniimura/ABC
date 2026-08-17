import ABC3.Found.Arakelov.PicInterface
import Mathlib.RingTheory.PicardGroup

/-!
# Arakelov (B1) 第 20 ブロック —— **アフィンでの大域切断**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★前層テンソルは各開集合ごとなので、大域切断は**自動的にモノイダル**

★★★これが第 16 ブロック(前層の段で組む)の**2 つ目の配当**である:

    (F ⊗ G)(U) = F(U) ⊗_{𝒪(U)} G(U)      -- 定義そのもの(`tensorObj_obj`)

★★したがって `Γ = (·)(⊤)` は乗法を保つ。
★★★**`equivPicRing` の乗法性が無料で出る**——`tilde` のモノイダル性は要らない。

## ★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `invertible_sections_top` | 可逆前層の大域切断は**可逆加群** |
| `picToRingPic` | ★`Pic (Spec R) →* CommRing.Pic Γ(Spec R, ⊤)` |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

/-- ★★★★**可逆前層の大域切断は可逆加群である**。

★機構は `L.isInv` を `⊤` で評価するだけ(前層テンソルは各開集合ごとだから)。 -/
theorem invertible_sections_top {X : Scheme.{u}} (L : InvertiblePresheaf X) :
    Module.Invertible ((X.presheaf.obj (op ⊤) : Type u))
      ((L.carrier.obj (op ⊤) : Type u)) := by
  obtain ⟨e⟩ := L.isInv
  exact Module.Invertible.left
    (((PresheafOfModules.evaluation _ (op (⊤ : X.Opens))).mapIso e).toLinearEquiv)

/-- ★大域切断の同型類(`CommRing.Pic` の元)。 -/
noncomputable def sectionsClass {X : Scheme.{u}} (L : InvertiblePresheaf X) :
    CommRing.Pic ((X.presheaf.obj (op ⊤) : Type u)) :=
  haveI := invertible_sections_top L
  CommRing.Pic.mk _ ((L.carrier.obj (op ⊤) : Type u))

/-- ★★★★★**大域切断は乗法を保つ**。

★★★これが「前層の段で組む」判断の配当である——前層テンソルは各開集合ごとなので、
`(L ⊗ M)(⊤) = L(⊤) ⊗ M(⊤)` が**定義そのもの**である。 -/
theorem sectionsClass_mul {X : Scheme.{u}} (L M : InvertiblePresheaf X) :
    sectionsClass (L.mul M) = sectionsClass L * sectionsClass M := by
  haveI := invertible_sections_top L
  haveI := invertible_sections_top M
  haveI := invertible_sections_top (L.mul M)
  have h1 : sectionsClass L * sectionsClass M
      = CommRing.Pic.mk ((X.presheaf.obj (op ⊤) : Type u))
        (TensorProduct ((X.presheaf.obj (op ⊤) : Type u))
          ((L.carrier.obj (op ⊤) : Type u)) ((M.carrier.obj (op ⊤) : Type u))) :=
    CommRing.Pic.mk_tensor.symm
  rw [h1]
  exact CommRing.Pic.mk_eq_mk_iff.2 ⟨LinearEquiv.refl _ _⟩

/-- ★同型な可逆前層は同じ類を与える。 -/
theorem sectionsClass_congr {X : Scheme.{u}} {L M : InvertiblePresheaf X}
    (e : L.carrier ≅ M.carrier) : sectionsClass L = sectionsClass M := by
  haveI := invertible_sections_top L
  haveI := invertible_sections_top M
  exact CommRing.Pic.mk_eq_mk_iff.2
    ⟨((PresheafOfModules.evaluation _ (op (⊤ : X.Opens))).mapIso e).toLinearEquiv⟩

/-- ★★★★★**`Pic X` から大域切断の Picard 群への群準同型**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★アフィン `X = Spec R` ではこれが `equivPicRing` の主役になる。 -/
noncomputable def picToSectionsPic (X : Scheme.{u}) :
    PicPre X →* CommRing.Pic ((X.presheaf.obj (op ⊤) : Type u)) where
  toFun := Quotient.lift sectionsClass (fun _ _ ⟨e⟩ => sectionsClass_congr e)
  map_one' := by
    show sectionsClass (InvertiblePresheaf.one X) = 1
    haveI := invertible_sections_top (InvertiblePresheaf.one X)
    exact CommRing.Pic.mk_eq_one_iff.2 ⟨LinearEquiv.refl _ _⟩
  map_mul' := by
    rintro ⟨L⟩ ⟨M⟩
    exact sectionsClass_mul L M

/-! ## ★出典の紐付け(`.src`) -/

def picToSectionsPic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Pic X から大域切断の Picard 群への群準同型)",
    sectionId := "genell-def-1-1-i" }

def invertible_sections_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆前層の大域切断は可逆加群)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
