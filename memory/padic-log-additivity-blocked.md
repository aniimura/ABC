---
name: padic-log-additivity-blocked
description: p進対数の準同型性を仕上げる2つの候補経路をどちらも試し、両方とも別の実装コストを要すると判明
metadata:
  type: project
---

pGC の p 進対数(`Found/PGC/PadicLog.lean`、級数の収束は sorry 無しで確立済み、
2026-09-04)について、残る「準同型性」`log((1+x)(1+y))=log(1+x)+log(1+y)` を
仕上げるための2つの候補経路をその場で検証した。**どちらも簡単な近道ではない**。

**経路①: `PowerSeries.log`(形式冪級数、`RingTheory/PowerSeries/Log.lean`)**
`coeff_log`/`deriv_log`/`map_log` はあるが、`logOf (f*g) = logOf f + logOf g` に
相当する積の加法性は**mathlib に無い**。2変数(`MvPowerSeries` か、係数環を
`PowerSeries ℚ` に取った入れ子の 1 変数冪級数)での構成が要り、形式微分による
「係数がすべて 0 → 冪級数として 0」という議論を自分で組む必要がある。

**経路②: `NormedSpace.exp`(`Analysis/Normed/Algebra/Exponential.lean`の一般 Banach 環指数関数)**
`NormedSpace.exp_add` は `[NormedAlgebra ℚ 𝔸]` を要求する。**この instance は
p進局所体 `K.carrier` には存在しない**(`infer_instance` で実測、2026-09-04)——
mathlib の `ℚ` の既定 `NormedField` インスタンスはアルキメデス的な絶対値であり、
`K.carrier` 上の ℚ の埋め込みが実際に持つ p 進絶対値と**両立しない**
(`‖1/p^m‖_ℚ→0` だが `‖1/p^m‖_{K,p進}=p^m→∞`)。使うには ℚ に p 進絶対値の
`NormedField` を **scoped で**別に立てる必要があり、それ自体が新しい instance
diamond のリスクを持つ小さくない作業。

**結論**: どちらの経路も「対数の収束」と同程度、あるいはそれ以上の実装コストを
要する——収束の確立(このセッションで sorry 無しに達成)より簡単な次の一手では
なかった。次にここへ戻るときは、経路②(scoped `NormedField ℚ` p進版を立てて
`NormedSpace.exp_add` を使う)の方が筋がよさそうだが、その instance が他の
無条件 `NormedAlgebra ℚ ?` 系の証明を汚染しないことを先に確認すること。

関連: [[measure-mathlib-before-skeleton]]
