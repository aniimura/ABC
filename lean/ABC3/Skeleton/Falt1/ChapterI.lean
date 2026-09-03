import ABC3.Meta.Claim
import ABC3.Interface.Falt1.PAdicHodgeTheory

/-!
# [Falt1] Chapter I, Theorem 4.4

原典: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299。
物理 p.17-18(印字 p.270-271)。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

`[LocProP] Lemma 2.1` が直接引用する箇所((i)・(iv))だけを Skeleton にする
((ii)(iii)(v)(vi)(vii) は `[LocProP]` が消費しないので省略、`.needs` に明記)。
-/

namespace ABC3.Skeleton.Falt1

open ABC3.Interface.Falt1

/-- **[Falt1] Theorem I.4.4, (i)(iv)**。

原文 (Falt1 p.270-271):
> 4.4. Theorem. (i) Hi(Δ,R^) ≅ Λi((R⊗V V̄)^(-1)d) ⊕ (rest), where "rest" is
> annihilated by p1/(p-1). ...
> (iv) The morphism ΩR/V → H1(Δ, ρ-1R^(1)) given by the extension Eρ is
> induced by a functorial isomorphism
> H1(Δ,R^)/(p-torsion) ≅ ΩR/V(dlog∞) ⊗R (R⊗V V̄)^(-1) -/
def theorem_I_4_4 (E : PAdicHodgeSetup) :=
  And.intro (fun j => Nonempty.intro (E.thm44i j)) (fun j => Nonempty.intro (E.thm44iv j))

def theorem_I_4_4.nonvacuous := theorem_I_4_4 PAdicHodgeSetup.example

def theorem_I_4_4.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 17, item := "Theorem 4.4", sectionId := "falt1-thm-i-4-4" }

def theorem_I_4_4.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      "本Skeletonは(i)(iv)だけを写した。原文は(ii)(iii)(v)(vi)(vii)も持つが
       [LocProP] Lemma 2.1は(i)(iv)しか引用しない(2026-09-04実測)。
       almost étale extension理論(Chapter I §2)がこの定理の証明の土台。" 17 ]

end ABC3.Skeleton.Falt1
