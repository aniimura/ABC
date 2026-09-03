import ABC3.Meta.Claim
import ABC3.Interface.Falt1.AlmostEtale

/-!
# [Falt1] Chapter I §2 —— Almost unramified extensions(4 項目)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)。物理 p.6-8(印字 p.259-261)。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。★§2 は 4/4。
-/

namespace ABC3.Skeleton.Falt1

open ABC3.Interface.Falt1

/-- **`[Falt1] Definition 2.1`** —— almost étale covering の定義。

内容 (Falt1 p.6、260dpi 目視。OCR 層は数式記号を激しく壊すため地の文で
写す): Suppose A is a ring, B an A-algebra. B is called an almost étale
covering of A if (i) B[1/p] is a projective A[1/p]-module of finite rank
and an étale A[1/p]-algebra; (ii) the trace map tr_{B/A}: B[1/p]→A[1/p]
maps B into A; (iii) for the idempotent e_{B/A} corresponding to the
diagonal, p^ε e_{B/A} lies in the image of B⊗_A B for any ε>0. -/
def isAlmostEtaleCovering (E : AlmostEtaleSetup) := E.isAlmostEtale

def isAlmostEtaleCovering.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 6, item := "Definition 2.1",
    sectionId := "falt1-def-2-1" }

/-- **`[Falt1] Theorem 2.2`**。

内容 (Falt1 p.7、260dpi 目視): Suppose B=A+mB is an almost étale
covering of A, C an A-algebra, I⊂C a nilpotent ideal, and φ: B→C/I an
A-algebra morphism. Then φ lifts uniquely to B→C. -/
def theorem_2_2 (E : AlmostEtaleSetup) := E.thm22

def theorem_2_2.nonvacuous := theorem_2_2 AlmostEtaleSetup.example

def theorem_2_2.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 7, item := "Theorem 2.2", sectionId := "falt1-thm-2-2" }

/-- **`[Falt1] Theorem 2.3`**。

内容 (Falt1 p.7、260dpi 目視): Suppose I⊂A is a nilpotent ideal, B̄ an
almost étale covering of Ā=A/I, B̄=Ā+mB̄. Then there exists an almost
étale covering B of A such that B̄ ≅ B⊗_A Ā/(p-torsion). -/
@[reducible] def theorem_2_3 (E : AlmostEtaleSetup) := E.thm23

@[reducible] def theorem_2_3.nonvacuous := theorem_2_3 AlmostEtaleSetup.example

def theorem_2_3.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 7, item := "Theorem 2.3", sectionId := "falt1-thm-2-3" }

/-- **`[Falt1] Theorem 2.4`**。

内容 (Falt1 p.8、260dpi 目視): Suppose B is an almost étale covering of
A. (i) The map Ω_A⊗_A B → Ω_B is an almost isomorphism, that is, its
kernel and cokernel are annihilated by m. (ii) Suppose a finite group G
operates on B, such that B[1/p] is a Galois covering of A[1/p] with
group G. If M is any B-module with a semilinear G-action, then m
annihilates all higher cohomology H^i(G,M), i>0. The same holds for
M^G/tr_G(M). -/
def theorem_2_4 (E : AlmostEtaleSetup) :=
  (E.thm24i_ker, E.thm24i_coker, E.thm24ii, E.thm24ii_tr)

def theorem_2_4.nonvacuous := theorem_2_4 AlmostEtaleSetup.example

def theorem_2_4.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 8, item := "Theorem 2.4", sectionId := "falt1-thm-2-4" }

end ABC3.Skeleton.Falt1
