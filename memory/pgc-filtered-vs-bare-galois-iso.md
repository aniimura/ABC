---
name: pgc-filtered-vs-bare-galois-iso
description: pGC Cor 3.1/3.3の形式化が「裸のΓ_K同型」を使っていた不整合を発見し訂正した(2026-09-04)
metadata:
  type: project
---

pGC の Corollary 3.1・Corollary 3.3(`lean/ABC3/Skeleton/PGC/Section3.lean`)は、
原文がどちらも明示的に「the **filtered** group Γ_K」から回復できると述べているのに、
形式化では `_α : ContinuousMulEquiv K.absGal K'.absGal`(**裸の**——フィルトレーション
と無関係な)同型を使っていた。`_RF : RamificationFiltration p` はパラメータとして
受け取りながら、`_α` の型を一切制約していない不整合(先頭の `_` が「未使用」の
意味そのままだった)——`Theorem 4.2`(`Section4.lean`)が正しく
`FilteredGroup.OuterIso` を使っているのと対照的。

**なぜこれが看過できないか(数論的根拠)**: `Γ_K` を裸の副有限群として見た同型類は
(奇素数 p では)`p` と `[K:Q_p]` だけで決まり、`K` 自身の他の情報(不分岐/分岐の
仕方など)には依らない(Iwasawa の型の古典的事実)。したがって裸の `α : Γ_K≅Γ_K'`
が存在するというだけでは、`K` に依存する任意の述語(`isHodgeTate` のような自由な
パラメータ)の不変性を導く根拠には**なりえない**——原文がわざわざ
「filtered group」と言っているのはこの理由による。この意味で旧来の形式化は
(証明できないという以前に)**主張として原文より強すぎ、おそらく偽**だった。

**訂正**: `_α` の型を `FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K')`
に直した(`filteredGroupOf` は `RamificationFiltration` から `FilteredGroup` を作る
小さな橋渡し——`ABC3.Interface.PGC.RamificationFiltration.filt` と同じ構成だが、
import の循環を避けるため `Section3.lean` に独立に複製した)。`RF` も `_RF` から
`RF` に変えて実際に型の中で使うようにした。まだ `sorry`(証明自体は Hodge-Tate
理論・uniformizing 判定条件がまだ無いのでできない)だが、**主張の形は原文に忠実に
なった**。

**一般化できる教訓**: 「`_` プレフィックスのパラメータが実は型の中で使われるべき
だった」という不整合は、原文が明示的に言及する構造(ここでは「filtered」という
一語)を見落とすと起きる。`Check/PGC/RefutationAttempts.lean` の監査は §1
(Prop 1.1・1.2・Cor 1.3)止まりで §3(Cor 3.1・3.3)は未監査だった——
セクションごとの反証可能性監査に抜けがないか、他のセクションも確認する価値がある。

関連: [[measure-mathlib-before-skeleton]]
