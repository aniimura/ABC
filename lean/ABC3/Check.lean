import ABC3.Meta.Claim
import ABC3.Check.GenEll.HeightAxiomGap
import ABC3.Check.GenEll.RemarkAxiomGap
import ABC3.Check.GenEll.Def12Degenerate
import ABC3.Check.GenEll.Def12NonDegenerate
import ABC3.Check.GenEll.NorthcottProjModelNonvacuous
import ABC3.Check.GenEll.RemarkSigmaWitness
import ABC3.Check.GenEll.Prop17Witness
import ABC3.Check.GenEll.Prop17AxiomGap
import ABC3.Check.GenEll.Prop17Direction
import ABC3.Check.GenEll.BelyiGeneralVacuous
import ABC3.Check.GenEll.Thm21AxiomGap
import ABC3.Check.GenEll.Thm21Witness
import ABC3.Check.NCBelyi.Thm25AxiomGap
import ABC3.Check.NCBelyi.Thm25Witness
import ABC3.Check.GenEll.HeightInterfaceNondegenerate
import ABC3.Check.GenEll.EllModuliDegInfPos
import ABC3.Check.GenEll.ImageSL2NeedsL5
import ABC3.Check.GenEll.VeluQuotOKNeedsL5
import ABC3.Check.GenEll.AlphaNeedsOneBasis
import ABC3.Check.GenEll.VeluTateNeedsChange
import ABC3.Check.GenEll.QuotMuNeedsHypothesis
import ABC3.Check.GenEll.AdicCompleteMissing
import ABC3.Check.GenEll.CompletionBridgeWitness
import ABC3.Check.GenEll.LcyclicExcTooStrong
import ABC3.Check.GenEll.Section3NotProvable
import ABC3.Check.AbsTopIII.LogShellLanding
import ABC3.Check.PGC.Section1Discriminating
import ABC3.Check.PGC.InertiaDegeneracy
import ABC3.Check.PGC.InertiaDegeneracyMoved
import ABC3.Check.PGC.DischargeFiresNothing
import ABC3.Check.PGC.ResidueCardinalityNondegenerate
import ABC3.Check.IUTchIII.LogVolumeFillsInterface
import ABC3.Check.IUTchIII.Cor312Degenerate
import ABC3.Check.PGC.RefutationAttempts
import ABC3.Check.PGC.Theorem42Degenerate
import ABC3.Check.PGC.Def32Degenerate
import ABC3.Check.PGC.Cor33Degenerate
import ABC3.Check.PGC.Prop22Degenerate
import ABC3.Check.PGC.Cor13Degenerate
import ABC3.Check.FrdI.TwistedFrobenioid
import ABC3.Check.FrdI.AutAmpleGap
import ABC3.Check.FrdI.Prop21QuantifierGap
import ABC3.Check.FrdI.VLocFalse
import ABC3.Check.Arakelov.ArcSpaceNondegenerate
import ABC3.Check.Arakelov.PicNondegenerate
import ABC3.Check.Arakelov.PullbackNondegenerate
import ABC3.Check.Arakelov.MetricNondegenerate
import ABC3.Check.Arakelov.ProperFlatNondegenerate
import ABC3.Check.GaloisRep.OmegaNondegenerate
import ABC3.Check.GaloisRep.HtFaltPinned
import ABC3.Check.Arakelov.ProjectiveCaseWeak
import ABC3.Check.GaloisRep.TorsionEquivWeak
/-!
# Check — 我々のモデルについての検査

**主語の分離**。`Skeleton/` が「原典が何を主張しているか」を置く場所であるのに対し、
ここは「**我々の**型設計が空虚に潰れていないか」を検査する宣言を置く場所。

置くもの:
- 識別力の検査(述語が常に真でも常に偽でもないこと)
- 負の対照(G3)の witness — `NegControl.witness` が名指しする宣言
- 非空虚 witness の補助

## 規則

- ここの宣言は原典の主張ではないので **`.src` を要求されない**。
- **`Skeleton/` から `Check/` を import してはならない**(`tools/check.mjs` が検査する)。
  逆向き(Check → Skeleton)のみ許す。原典についての主張が、我々のモデルについての
  検査に依存する形を構造的に禁じるため。
- 「G1 を避けるための逃げ道」として使わないこと。原典の主張はここに置けない。
-/
