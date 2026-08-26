/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.CategoryVocabulary

/-!
# [FrdI] Theorem 6.2, (iii) —— mono がすべて同型なら FSM 型

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> (Geometric Frobenioids) For i = 1, 2, let Vi be a proper

## ★★原文が「immediately」で畳んだ所

原文は `Theorem 6.2, (iii)` で「`𝒟` の mono はすべて同型なので `𝒟` は FSM 型、
したがって FSMFF 型」と述べる。★中身は**定義を開くだけ**である ——
`FSM-morphism` の定義(`§0`)が `fiberwise surjective ∧ mono` だから、
「mono ⟹ 同型」があれば `IsOfFSMType` はそのまま出る。

★★`FSM 型 ⟹ FSMFF 型` は在庫にある(`isOfFSMFFType_of_isOfFSMType`)。

## ★具体例では成り立っている

`Sec6GaloisCat.lean` の `finSubOp_isOfFSMType` が
`𝒟 = B(G)⁰ = (FinSub K K̄)ᵒᵖ` について**実際に**これを与えている
(中身は「`FinSub` の epi は同型」で、無限 Galois 理論の
`fixedField (fixingSubgroup L) = L` 1 本)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `isOfFSMType_of_mono_isIso` | ★mono がすべて同型なら FSM 型 |
| `isOfFSMFFType_of_mono_isIso` | ★FSMFF 型まで |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u

variable (C : Type u) [Category.{v} C]

/-- ★★★**mono がすべて同型なら `𝒞` は FSM 型**。

★`FSM-morphism = fiberwise surjective ∧ mono` なので、定義を開くだけ。 -/
theorem isOfFSMType_of_mono_isIso (h : ∀ {A B : C} (f : A ⟶ B), Mono f → IsIso f) :
    IsOfFSMType C :=
  fun _ _ β hβ => h β hβ.2

/-- ★★**mono がすべて同型なら `𝒞` は FSMFF 型**。 -/
theorem isOfFSMFFType_of_mono_isIso (h : ∀ {A B : C} (f : A ⟶ B), Mono f → IsIso f) :
    IsOfFSMFFType C :=
  isOfFSMFFType_of_isOfFSMType (isOfFSMType_of_mono_isIso C h)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Theorem 6.2, (iii)` の「`𝒟` の mono はすべて同型 ⟹ FSM 型」。 -/
def isOfFSMType_of_mono_isIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (iii) — 𝒟 の mono がすべて同型なら FSM 型(ゆえに FSMFF 型)",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.FrdI
