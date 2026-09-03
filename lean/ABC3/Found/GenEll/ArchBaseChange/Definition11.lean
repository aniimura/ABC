/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArchConj

/-!
# ArchBaseChange —— `[GenEll] Definition 1.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField
variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]

/-! ## ★★制限した素点の埋め込みは共役を除いて決まる -/

/-- ★**制限した素点は、埋め込みの合成が定める素点である**。 -/
theorem comap_eq_mk {k : Type*} [Field k] (w : InfinitePlace K) (f : k →+* K) :
    w.comap f = InfinitePlace.mk (w.embedding.comp f) := by
  conv_lhs => rw [← InfinitePlace.mk_embedding w]
  exact InfinitePlace.comap_mk _ _

/-- ★★**制限した素点の埋め込みは、合成した埋め込み**または**その共役**である。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**これが「計量が `ι_X` と両立する」ことを要求する理由である。**
共役の分だけ ℂ-点がずれるので、Green 関数が共役不変でないと
高さが定義体の取り方に依ってしまう。 -/
theorem embedding_comap_dichotomy {k : Type*} [Field k] (w : InfinitePlace K)
    (f : k →+* K) :
    (w.comap f).embedding = w.embedding.comp f
      ∨ (w.comap f).embedding
          = NumberField.ComplexEmbedding.conjugate (w.embedding.comp f) := by
  rw [comap_eq_mk]
  exact InfinitePlace.embedding_mk_eq _

def embedding_comap_dichotomy.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(計量と ι_X の両立を要求する理由——共役の場合分け)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
