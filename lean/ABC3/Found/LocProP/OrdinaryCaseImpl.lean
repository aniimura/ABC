import ABC3.Meta.Claim
import ABC3.Interface.LocProP.OrdinaryCase

/-!
# [LocProP] §1 —— 実装(`OrdinaryCaseSetup` 上で)

`Skeleton/LocProP/Section1.lean` から委譲される実装本体。いずれも
`Interface/LocProP/OrdinaryCase.lean` が posit している主張そのものを
取り出すだけなので証明は空——本体は posit 側にある。
-/

namespace ABC3.Found.LocProP

open ABC3.Interface.LocProP

/-- **[LocProP] Lemma 1.2** の実装。 -/
def thetaUnique := @OrdinaryCaseSetup.theta_unique

/-- **[LocProP] Lemma 1.3** の実装。 -/
def kerEqIInter := @OrdinaryCaseSetup.ker_eq_iInter

/-- **[LocProP] Lemma 1.4** の実装。 -/
def alphaEtEqTheta := @OrdinaryCaseSetup.alphaEt_eq_theta

theorem thetaUnique.nonvacuous :
    thetaUnique OrdinaryCaseSetup.example OrdinaryCaseSetup.example.theta
      (fun _ => rfl) (fun _ _ => Subsingleton.elim _ _) = rfl :=
  rfl

theorem kerEqIInter.nonvacuous :
    kerEqIInter OrdinaryCaseSetup.example
      = OrdinaryCaseSetup.example.ker_eq_iInter :=
  rfl

theorem alphaEtEqTheta.nonvacuous (x : OrdinaryCaseSetup.example.PtX) :
    alphaEtEqTheta OrdinaryCaseSetup.example x
      = OrdinaryCaseSetup.example.alphaEt_eq_theta x :=
  rfl

def thetaUnique.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 18, item := "Lemma 1.2", sectionId := "locprop-lemma-1-2" }

def kerEqIInter.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 18, item := "Lemma 1.3", sectionId := "locprop-lemma-1-3" }

def alphaEtEqTheta.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 19, item := "Lemma 1.4", sectionId := "locprop-lemma-1-4" }

end ABC3.Found.LocProP
