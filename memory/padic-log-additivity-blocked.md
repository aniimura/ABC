---
name: padic-log-additivity-blocked
description: 対数側の準同型性は依然未解決だが、指数側(padicExp)を自前構築して加法性まで到達した(2026-09-04)
metadata:
  type: project
---

★★2026-09-04 追記: **経路③(自前で `padicExp` を組む)が指数の加法性まで到達した。**
`Found/PGC/PadicExp.lean` に `padicExp`・`padicExp_add`(`exp(x+y)=exp(x)·exp(y)`)を
`sorry` 無しで実装済み——`NormedAlgebra ℚ K.carrier` も形式冪級数も使わず、
`K.carrier` 上で直接 `expTerm`(`xⁿ/n!`)を定義し、Cauchy 積
(`Summable.tsum_mul_tsum_eq_tsum_sum_antidiagonal`)+ 二項定理(`add_pow`)
だけで閉じた。**対数の加法性より指数の加法性のほうが本質的に易しい**
(微分を経由しない、純粋に組み合わせ的な議論で済む)ことが判明——
下記①②の記録は「対数を直接攻める」経路の記録として残すが、**次にここへ
戻るときは③の続き(`exp`/`log` の互逆性)から入るのがよい**。

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

**結論(当時)**: ①②はどちらも「対数の収束」と同程度、あるいはそれ以上の実装
コストを要すると見えた。実際には③(対数を直接攻めず、指数を自前で組む)が
迂回路として機能した——**詰まったら「同じ対象を別の道具で直接攻める」のではなく
「双対のより易しい対象を作って橋渡しする」を検討すること**、が一般化できる教訓。

**残る作業**: `padicExp K (padicLog K x + padicLog K y) = padicExp K (padicLog K (x*y))`
のような互逆性(`exp∘log=id`・`log∘exp=id`)が示せれば、`padicExp_add` から
`log` の加法性が微分無しで出る。互逆性の証明は `exp`/`log` の係数どうしの合成公式
(スターリング数的な組み合わせ論)を要し、加法性(二項定理だけ)より重いが、
`Found/PGC/PadicExp.lean`/`PadicLog.lean` の道具立て(`expTerm`/`logTerm`・
Cauchy 積・`NonarchimedeanRing` instance)はすでに揃っている。

関連: [[measure-mathlib-before-skeleton]]
