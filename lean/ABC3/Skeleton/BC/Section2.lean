import ABC3.Meta.Claim
import ABC3.Interface.BC.ComparisonIsomorphism

/-!
# [BC] Theorem 2.2.3 (Faltings)

原典: O. Brinon, B. Conrad, *CMI Summer School Notes on p-adic Hodge Theory*
(2009)。物理 p.13。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

`[LocProP] Proposition 2.4` が([Falt2] の代替として)引く箇所。
-/

namespace ABC3.Skeleton.BC

open ABC3.Interface.BC

/-- **`[BC] Theorem 2.2.3 (Faltings)`**。

原文 (BC p.13):
> Theorem 2.2.3 (Faltings). Let K be a p-adic field. For smooth proper
> K-schemes X, there is a canonical isomorphism -/
@[reducible] def theorem_2_2_3 (E : ComparisonIsoSetup) := Nonempty.intro E.thm223

@[reducible] def theorem_2_2_3.nonvacuous := theorem_2_2_3 ComparisonIsoSetup.example

def theorem_2_2_3.src : ABC3.Meta.Source :=
  { paper := "BC", pdfPage := 13, item := "Theorem 2.2.3", sectionId := "bc-thm-2-2-3" }

def theorem_2_2_3.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      "本Skeletonは式(2.2.1)の同型の存在だけを写した。スキーム論・
       エタールコホモロジー(H^n_et)・Hodgeコホモロジー(H^{n-q}(X,Ω^q))
       そのものはposit(抽象的な加法群の対)であり、構成しない。" 13 ]

end ABC3.Skeleton.BC
