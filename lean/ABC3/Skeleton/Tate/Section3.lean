import ABC3.Meta.Claim
import ABC3.Interface.Tate.PDivisibleGroups

/-!
# [Tate] §3.3, Theorem 1・Theorem 2

原典: J. Tate, *p-Divisible Groups*(1967)pp.158-183(冒頭13頁のみ入手)。
物理 p.10(印字 p.176-177)。**目視確認 2026-09-04**(逐語の食い違い無し)。

`[LocProP] Lemma 2.2` が直接引用する箇所。
-/

namespace ABC3.Skeleton.Tate

open ABC3.Interface.Tate

/-- **`[Tate] §3.3, Theorem 1`**。

内容 (Tate p.10、260dpi 目視、OCR 層は文字化けのため逐語照合対象外——
Interface/Tate/PDivisibleGroups.lean 冒頭の注記参照): Theorem 1. We have
H0(𝒢, C) = K, and H1(𝒢, C) is a one-dimensional vector space over K. -/
def theorem_1 (E : TateCohomologySetup) :=
  And.intro (Nonempty.intro E.thm1i) (Nonempty.intro E.thm1ii)

def theorem_1.nonvacuous := theorem_1 TateCohomologySetup.example

def theorem_1.src : ABC3.Meta.Source :=
  { paper := "Tate", pdfPage := 10, item := "Theorem 1", sectionId := "tate-thm-1" }

/-- **`[Tate] §3.3, Theorem 2`**。

内容 (Tate p.10、260dpi 目視、OCR 層は文字化けのため逐語照合対象外——
Interface/Tate/PDivisibleGroups.lean 冒頭の注記参照): Theorem 2. Suppose
that there is a finite extension K0 of K contained in K∞ such that
K∞/K0 is totally ramified and Gal(K∞/K0) ≃ Zp. Then H0(𝒢, C(χ)) = 0
and H1(𝒢, C(χ)) = 0. -/
def theorem_2 (E : TateCohomologySetup) := E.thm2

def theorem_2.nonvacuous := theorem_2 TateCohomologySetup.example

def theorem_2.src : ABC3.Meta.Source :=
  { paper := "Tate", pdfPage := 10, item := "Theorem 2", sectionId := "tate-thm-2" }

end ABC3.Skeleton.Tate
