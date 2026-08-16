import ABC3.Meta.Claim
import ABC3.Skeleton.PGC.Setup
import ABC3.Skeleton.PGC.Section1
import ABC3.Skeleton.PGC.Section1Cor13
import ABC3.Skeleton.AbsTopIII.LogShell
import ABC3.Skeleton.IUTchI.InitialThetaData
import ABC3.Skeleton.IUTchIII.Cor312
import ABC3.Skeleton.GenEll.Heights
import ABC3.Skeleton.GenEll.GaloisImage
/-!
# Skeleton — 論文の主張(証明しない)

Track A。原典の主張を statement として置く。証明は付けない(`sorry`)。

各宣言 `foo` は `foo.src : ABC3.Meta.Source` を伴う
(`tools/check.mjs` が G1 として検査し、`1_Structured` の該当 section と照合する)。
-/
