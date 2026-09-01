import ABC3.Meta.Claim
import ABC3.Skeleton.PGC.Setup
import ABC3.Skeleton.PGC.Section1
import ABC3.Skeleton.PGC.Section1Cor13
import ABC3.Skeleton.AbsTopIII.LogShell
import ABC3.Skeleton.IUTchI.InitialThetaData
import ABC3.Skeleton.IUTchIII.Cor312
import ABC3.Skeleton.GenEll.Heights
import ABC3.Skeleton.GenEll.EllModuliWitness
import ABC3.Skeleton.GenEll.GaloisLocal
import ABC3.Skeleton.GenEll.TateIsogeny
import ABC3.Skeleton.GenEll.TateLocalModel
import ABC3.Skeleton.GenEll.TateODE
import ABC3.Skeleton.GenEll.SigmaConvolution
import ABC3.Skeleton.GenEll.GaloisImage
import ABC3.Skeleton.GenEll.Section1
import ABC3.Skeleton.GenEll.DeligneHeight
import ABC3.Skeleton.GenEll.IsogenyHeight
import ABC3.Skeleton.GenEll.Section2
import ABC3.Skeleton.GenEll.Section2Converse
import ABC3.Skeleton.GenEll.Section3
import ABC3.Skeleton.GenEll.Uniformization
import ABC3.Skeleton.GenEll.AdditionTheorem
import ABC3.Skeleton.GenEll.LatticeFromInvariants
import ABC3.Skeleton.GenEll.Section4
import ABC3.Skeleton.NCBelyi.Theorem25
import ABC3.Skeleton.FrdI.Def28ProL
import ABC3.Skeleton.FrdI.Lemma65SixExp
import ABC3.Skeleton.FrdI.Thm64Deg
import ABC3.Skeleton.FrdI.Thm62Slim
import ABC3.Skeleton.GaloisRep.TateUniformization
import ABC3.Skeleton.GaloisRep.WeilFunctionField
import ABC3.Skeleton.GaloisRep.PointReduction
import ABC3.Skeleton.GaloisRep.WeilDivisor
import ABC3.Skeleton.GaloisRep.WeilRoot
import ABC3.Skeleton.GaloisRep.WeilPairingDef
import ABC3.Skeleton.GaloisRep.WeilPairing
import ABC3.Skeleton.Divisor.SchemeWeil
import ABC3.Skeleton.Divisor.Cartier
import ABC3.Skeleton.Divisor.ArithDivisor
import ABC3.Skeleton.Divisor.Normalization
import ABC3.Skeleton.Divisor.Hartogs
import ABC3.Skeleton.Divisor.NormalizationUniversal
import ABC3.Skeleton.NumberField.Chebotarev
/-!
# Skeleton — 論文の主張(証明しない)

Track A。原典の主張を statement として置く。証明は付けない(`sorry`)。

各宣言 `foo` は `foo.src : ABC3.Meta.Source` を伴う
(`tools/check.mjs` が G1 として検査し、`1_Structured` の該当 section と照合する)。
-/
