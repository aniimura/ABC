/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop33UnTr

/-!
# [FrdI] Proposition 4.8, (iii) —— Frobenius-compact 対象の移送(下ごしらえ)

原文 (FrdI p.88):
> (iii) If C is of rationally standard type, then (Cistr)birat is of standard

★★`Frobenius-compact` は `𝒪^×` についての 3 条である。
★`Remark 4.5.1` が与えるのは `((𝒞^istr)^un-tr)^birat` の compact 対象で、
`Proposition 4.8, (iii)` が要求するのは `(𝒞^istr)^birat` のもの ——
**向きが逆**なのが `iii-b-compact` の障害だった。

## ★★★測って分かったこと(2026-08-19)

第 2 条(無限位数の単元が存在する)は**準同型の逆向きに落ちる**
(`infiniteOrder_of_map`、`Prop48Sec.lean`)。
したがって **`𝒪^×` の写像が全射でありさえすれば第 2 条は渡る**。

★その全射性は `Proposition 4.4, (ii)` の記述
(`𝒪^×(A^birat)` は `𝒪^▷(A)^gp` の像)を通し、
**`𝒪^▷(𝒞^istr) ↠ 𝒪^▷((𝒞^istr)^un-tr)`** に帰着する。
★本ファイルはその**最初の 1 本**を取る。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★★★★**`𝒪^▷` は `𝒞^istr → 𝒞^un-tr` で全射**。

★★`𝒞^un-tr` は射を `𝔽_Φ` の像で同一視した商なので、
`Base` も `degFr` も**定義的に一致する** ——
したがって `istrToUnTr` の充満性(在庫)で代表元を取れば、
それは自動的に `𝒪^▷` に入る。 -/
theorem otri_unTr_surjective (Fc : FrobenioidCore P) (A : UnTr P) (β : End A)
    (hβ : β ∈ OTri (unTrPre P Fc) A) :
    ∃ β₀ : End (show Istr P from A),
      β₀ ∈ OTri (istrPre P Fc) (show Istr P from A) ∧ (istrToUnTr P).map β₀ = β := by
  obtain ⟨β₀, hβ₀⟩ := (istrToUnTr P).map_surjective β
  refine ⟨β₀, ⟨?_, ?_⟩, hβ₀⟩
  · have h : (unTrPre P Fc).Base ((istrToUnTr P).map β₀)
        = (unTrPre P Fc).Base (𝟙 A) := by
      rw [hβ₀]
      exact hβ.1
    exact h
  · have h : (unTrPre P Fc).degFr ((istrToUnTr P).map β₀) = 1 := by
      rw [hβ₀]
      exact hβ.2
    exact h

def otri_unTr_surjective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — 𝒪^▷ は 𝒞^istr → 𝒞^un-tr で全射",
    sectionId := "frdi-prop-4-8" }

/-! ## ★★★`Remark 4.5.1` —— standard 型は `𝒞^istr` で保たれる

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is

★★`IsOfStandardType` は**フィールド 6 個**の構造体。
`𝒟` も `Φ` も変わらないので 2 つは**そのまま**、
`quasiIsotropic` は在庫、`groupLikeCompact` は
`Istr (istrPre P F)` が `Istr P` の全対象であることから元の条が効く。
★残る 2 つを本節で取る。 -/

/-- ★★★★**`𝒞^istr` は Frobenius-isotropic 型** —— **自明に成り立つ**。

★★`𝒞^istr` の対象はすべて isotropic なので、`Dd := A`、`φ := 𝟙 A` で足りる。
★元の `𝒞` の条件は要らない。 -/
theorem istr_isOfFrobeniusIsotropicType (F : FrobenioidCore P) :
    IsOfFrobeniusIsotropicType (istrPre P F) := fun A =>
  ⟨A, 𝟙 A, isFrobeniusType_of_isIso (istrPre P F) (𝟙 A), istr_isotropic P F A⟩

/-- ★★★★**`Frobenius-normalized` は `𝒞^istr` へ移る**。

★`istr_compat_Base` / `istr_compat_degFr` が `rfl` なので、
条件は `𝒞` のものと**定義的に同じ形**になる。 -/
theorem istr_isFrobeniusNormalized (F : FrobenioidCore P) (A : Istr P)
    (h : IsFrobeniusNormalized P A.obj) : IsFrobeniusNormalized (istrPre P F) A := by
  intro φ hφ α hα
  refine InducedCategory.hom_ext ?_
  have hb : IsBaseIdentity P (Functor.mapEnd A ((isotropicProp P).ι) φ) := hφ
  have hmem : (Functor.mapEnd A ((isotropicProp P).ι) α) ∈ OTri P A.obj := ⟨hα.1, hα.2⟩
  have hkey := h (Functor.mapEnd A ((isotropicProp P).ι) φ) hb _ hmem
  show ((isotropicProp P).ι).map (φ ≫ _) = ((isotropicProp P).ι).map (_ ≫ φ)
  refine (((isotropicProp P).ι).map_comp _ _).trans ?_
  refine Eq.trans ?_ ((((isotropicProp P).ι).map_comp _ _).symm)
  refine Eq.trans ?_ hkey
  exact congrArg (fun t => ((isotropicProp P).ι).map φ ≫ t)
    ((Functor.mapEnd A ((isotropicProp P).ι)).map_pow α _)

def istr_isOfFrobeniusIsotropicType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — 𝒞^istr は Frobenius-isotropic 型",
    sectionId := "frdi-remark-4-5-1" }

end ABC3.Found.FrdI
