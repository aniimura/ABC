import ABC3.Meta.Claim
import ABC3.Interface.Falt1.GoodReduction

/-!
# [Falt1] Chapter I §3 —— Good reduction(2 項目)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)。物理 p.10・12(印字 p.263・265)。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。★§3 は 2/2。
-/

namespace ABC3.Skeleton.Falt1

open ABC3.Interface.Falt1

/-- **`[Falt1] Theorem 3.1`**。

内容 (Falt1 p.10、260dpi 目視): Suppose S is a normal finite R-algebra
such that S[1/p] is étale over R[1/p]. Let S_∞ denote the normalization
of S⊗_R R_∞. Then S_∞ is an almost étale covering of R_∞. -/
def theorem_3_1 (E : GoodReductionSetup) := E.thm31

def theorem_3_1.nonvacuous := theorem_3_1 GoodReductionSetup.example

def theorem_3_1.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 10, item := "Theorem 3.1", sectionId := "falt1-thm-3-1" }

/-- **`[Falt1] Theorem 3.2`**。

内容 (Falt1 p.12、260dpi 目視): Suppose S is a finite normal torsionfree
R-algebra such that S[1/(pu_1...u_d)] is étale over R[1/(pu_1...u_d)].
If S_n denotes the normalization of S⊗_R R_n and S_∞ the union of all
S_n, then S_∞ is almost étale over R_∞. -/
def theorem_3_2 (E : GoodReductionSetup) := E.thm32

def theorem_3_2.nonvacuous := theorem_3_2 GoodReductionSetup.example

def theorem_3_2.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 12, item := "Theorem 3.2", sectionId := "falt1-thm-3-2" }

end ABC3.Skeleton.Falt1
