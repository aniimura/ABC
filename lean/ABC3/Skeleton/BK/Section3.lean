import ABC3.Meta.Claim
import ABC3.Interface.BK.DeRhamComparison

/-!
# [BK] Lemma 3.8.1・Example 3.11

原典: S. Bloch, K. Kato, *L-Functions and Tamagawa Numbers of Motives*
(1990) pp.333-400。物理 p.23(Lemma 3.8.1)・物理 p.29(Example 3.11)。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

`[LocProP] Lemma 4.1` が直接引用する箇所。
-/

namespace ABC3.Skeleton.BK

open ABC3.Interface.BK

/-- **`[BK] Lemma 3.8.1`**。

内容 (BK p.23、260dpi 目視、OCR 層は文字化けのため逐語照合対象外——
Interface/BK/DeRhamComparison.lean 冒頭の注記参照): Let V be a de Rham
representation of Gal(K̄/K). Then H¹(K, B⁺_dR ⊗ V) → H¹(K, B_dR ⊗ V)
is injective. -/
def lemma_3_8_1 (E : DeRhamComparisonSetup) := E.lem381_injective

def lemma_3_8_1.nonvacuous := lemma_3_8_1 DeRhamComparisonSetup.example

def lemma_3_8_1.src : ABC3.Meta.Source :=
  { paper := "BK", pdfPage := 23, item := "Lemma 3.8.1", sectionId := "bk-lemma-3-8-1" }

/-- **`[BK] Example 3.11`**、式 (3.11.2)。

内容 (BK p.29、260dpi 目視、OCR 層は文字化けのため逐語照合対象外——
Interface/BK/DeRhamComparison.lean 冒頭の注記参照): Let A be an abelian
variety over K, T its total Tate module. Then H¹_e(K,T) = H¹_f(K,T)
= H¹_g(K,T). -/
def example_3_11 (E : DeRhamComparisonSetup) :=
  And.intro (Nonempty.intro E.ex311_ef) (Nonempty.intro E.ex311_fg)

def example_3_11.nonvacuous := example_3_11 DeRhamComparisonSetup.example

def example_3_11.src : ABC3.Meta.Source :=
  { paper := "BK", pdfPage := 29, item := "Example 3.11", sectionId := "bk-example-3-11" }

end ABC3.Skeleton.BK
