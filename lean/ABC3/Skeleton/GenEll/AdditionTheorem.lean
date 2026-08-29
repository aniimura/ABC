import ABC3.Found.GenEll.Uniformization

/-!
# スケルトン —— **`℘` の加法定理**（`Skeleton`）★★★★★★★★★★**閉じた**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★なぜこの節点を立てたのか——代数側と解析側を繋ぐ最後の一本

2026-08-29（第 586-605）に `Lemma 3.5` の周りが次まで進んだ:

| 側 | 状態 | どこ |
|---|---|---|
| アルキメデス | ★「次数 `l` の同種写像がある」だけから従う | `Found/GaloisRep/VeluNormalized.lean`（第 596） |
| 有限素点 | ★`E` が大域極小・`E′` が整なら自動 | 同上＋`neronExp_nonneg`（第 595） |
| Vélu 代数側 | ★定義・不変量・`l=2,3` の商・正規化（一般の `l`）・代表系不変 | `Found/GenEll/Velu.lean`（第 586-593） |
| Vélu 解析側 | ★`T` が `Λ′/Λ` の代表系であることだけから従う | `Found/GenEll/Uniformization.lean`（第 597-602） |
| 一様化 | ★**全射性が塞がった** | 同上（第 603-604） |
| 2-捩れ | ★代数側の場合分けが解析側で正当化された | 同上（第 605） |

☆★**残るのは 1 本**——**`℘` の加法定理**である。これが要るのは 2 か所:

1. 解析側の `℘_{Λ′}(z) = Σ_{w∈T} ℘_Λ(z+w) − c`（第 601-602）を
   代数側の `veluXGen`（第 591）へ翻訳するとき
2. 一様化の**単射性**と**群同型**（`ℂ/Λ ≅ E(ℂ)`）を出すとき
   ——これが `H ⊆ E(ℂ)` の位数 `l` の部分群を `Λ ⊆ Λ′` 指数 `l` に対応させる

## ★★★★★道具は揃っている

★`Found/GenEll/Uniformization.lean` の **`elliptic_liouville`**（第 598、
「整で二重周期的なら定数」）が本セッションで 4 度効いた（第 600・601・603・604）。
加法定理も同じ型の議論である:

* `F(z) ≔ ℘(z+w) + ℘(z) + ℘(w) − (1/4)·((℘′(z)−℘′(w))/(℘(z)−℘(w)))²` は `Λ`-周期的
* `z ∈ Λ` で `F` は解析的（`℘(z) ~ z⁻²` と `(1/4)(…)² ~ z⁻² + 2℘(w)` が打ち消し合う）
* `z ≡ −w` でも打ち消し合う（`℘(z+w)` の極と `(1/4)(…)²` の極）
* したがって `F` は整で二重周期的、Liouville で定数、`z → 0` で `F → 0`

☆★**残る手間は Laurent 展開**（`℘(t) = t⁻² + O(t²)`・`℘′(t) = −2t⁻³ + O(t)`）を
mathlib の `weierstrassPExcept`・`weierstrassPSeries` から取り出すところである。
★mathlib は `G n`（Eisenstein 級数）を持っているので、係数は書ける。

## ☆mathlib の在庫（2026-08-29 に測定）

★**ある**: `℘`・`℘′`・`g₂`・`g₃`・`G n`・周期性・偶奇・`deriv ℘ = ℘′`・
`meromorphic_weierstrassP`・`order_weierstrassP`・`derivWeierstrassP_sq`・
`weierstrassPExcept`（`l₀` 項を抜いた `℘`）。

☆**無い**: 加法定理、群同型 `ℂ/Λ ≅ E(ℂ)`、楕円関数の Liouville
（★本プロジェクトが第 598 で建てた）。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll

/-- ☆★★★★★★★★★★**`℘` の加法定理**

    `℘(z+w) = (1/4)·((℘′(z) − ℘′(w))/(℘(z) − ℘(w)))² − ℘(z) − ℘(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`latticeCurve P = ⟨0,0,0,−g₂/4,−g₃/4⟩` の群法則（`a₁ = a₂ = a₃ = 0`）は
`x₃ = λ² − x₁ − x₂`、`λ = (y₂−y₁)/(x₂−x₁)` である。
`x = ℘`・`y = ℘′/2` を入れると `λ = (℘′(w)−℘′(z))/(2(℘(w)−℘(z)))` なので、
上の等式は **`(℘, ℘′/2)` が群準同型である**ことと同値である。

☆**まだ証明していない**——これが `Lemma 3.5` に残る最後の一本である。 -/
theorem weierstrassP_add (P : PeriodPair) (z w : ℂ)
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice) (h2w : 2 * w ∉ P.lattice)
    (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z ≠ P.weierstrassP w) :
    P.weierstrassP (z + w)
      = ((P.derivWeierstrassP z - P.derivWeierstrassP w)
          / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4
        - P.weierstrassP z - P.weierstrassP w :=
  ABC3.Found.GenEll.weierstrassP_addition P w hw h2w hz hzw (sub_ne_zero.2 hne)

/-! ## ★出典の紐付け（`.src`）と、証明が要求するもの（`.needs`） -/

def weierstrassP_add.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘ の加法定理——代数側と解析側を繋ぐ最後の一本)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_add.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★★★2026-08-29 の到達点(第 586-605、16 ブロック): Lemma 3.5 の周りが" ++
       "アルキメデス側(第 596: 次数 l の同種写像があるだけから従う)、" ++
       "有限側(第 595: E が大域極小・E′ が整なら自動)、" ++
       "Vélu 代数側(第 586-593)、Vélu 解析側(第 597-602: 代表系だけから)、" ++
       "一様化の全射性(第 603-604)、2-捩れ(第 605)まで進んだ。" ++
       "☆残るのは本節点 1 本である") 17,
    .implicitStep
      ("★★★★★道具は揃っている: Found/GenEll/Uniformization.lean の elliptic_liouville" ++
       "(第 598、整で二重周期的なら定数)が本セッションで 4 度効いた" ++
       "(第 600 極の打ち消し・第 601 正規化・第 603 ℘ の全射性・第 604 一様化の全射性)。" ++
       "★加法定理も同じ型の議論で、F(z) ≔ ℘(z+w) + ℘(z) + ℘(w) − (1/4)((℘′(z)−℘′(w))/(℘(z)−℘(w)))²" ++
       "が整で二重周期的であることを示せばよい") 13,
    .implicitStep
      ("★★★★★★2026-08-29(第 607-608)で証明の核が取れた。" ++
       "Found/GenEll/Uniformization.lean の laurentB(= z²℘(z) の解析接続、原点で 1)・" ++
       "laurentA(= z³℘′(z) の解析接続、原点で −2)と、laurent_cancel:" ++
       "4·B·(B − z²x)² − (A − z³y)² = z²·M(z) —— ★ring 一発の恒等式。" ++
       "☆打ち消しの正体は B(0) = 1・A(0) = −2 から 4·1·1 − (−2)² = 0 である") 13,
    .implicitStep
      ("☆残る手順は 2 つだけである:" ++
       "(1) laurent_cancel から「F は原点を除いた近傍で解析関数と一致する」を出す" ++
       "(z² が約されるので有界、Riemann の除去可能特異点。" ++
       "★第 603 の wp_inv_differentiableAt_of_mem と同じ型の議論)。" ++
       "★周期性から格子の全点で従う。" ++
       "(2) z ≡ −w での極の打ち消し。★そちらは ℘ の 2 階 Taylor が要る" ++
       "(℘(t−w) − ℘(w) が t = 0 で 1 位の零点であることと、その 2 階係数)。" ++
       "☆mathlib の AnalyticAt.exists_eventuallyEq_pow_smul_nonzero_iff が使える") 13,
    .implicitStep
      ("★★★★★(1) は 2026-08-29(第 610)で塞がった: " ++
       "Found/GenEll/Uniformization.lean の addDefect・addDefectNear・" ++
       "addDefect_eq_near(z ≠ 0 で一致)・addDefectNear_zero(原点で 0)・" ++
       "analyticAt_addDefectNear。☆すなわち F_w の解析接続は原点で消える" ++
       "——加法定理が原点で成り立つことは Lean 上で確かめられた。" ++
       "★周期性から格子の全点で従う") 13,
    .implicitStep
      ("★(2) の最小形(2026-08-29 に紙の上で詰めた): t ≔ z + w と置き、" ++
       "u(t) ≔ ℘(t−w) − ℘(w)、v(t) ≔ ℘′(t−w) − ℘′(w)、û ≔ u/t とすると" ++
       "F_w(t−w) = (4û² − v²)/(4t²û²) + e(t) + ℘(t−w) + ℘(w) であり、" ++
       "★4û² − v² = (2û − v)(2û + v) なので 2û − v が 2 位の零点であればよい") 8,
    .implicitStep
      ("★★(2) をさらに落とすと: q(t) ≔ 2u(t) − t·u′(t) + t·℘′(w) と置けば " ++
       "2û − v = q/t であり、q は t = 0 で 3 位の零点である。" ++
       "★q(0) = 0(u(0) = 0)、q′(0) = u′(0) + ℘′(w) = ℘′(−w) + ℘′(w) = 0(℘′ は奇)、" ++
       "★★q″(t) = −t·u‴(t) なので q″(0) = 0 —— **自動的に消える**。" ++
       "☆Lean では mathlib の natCast_le_analyticOrderAt" ++
       "(n ≤ analyticOrderAt f z₀ ↔ ∃ g, f = (z−z₀)^n • g)で因数分解できる") 8,
    .implicitStep
      ("★なぜ要るか(2 か所): (1) 解析側の ℘_{Λ′}(z) = Σ_{w∈T} ℘_Λ(z+w) − c(第 601-602)を" ++
       "代数側の veluXGen(第 591)へ翻訳するとき、" ++
       "(2) 一様化の単射性と群同型 ℂ/Λ ≅ E(ℂ) を出すとき" ++
       "——これが H ⊆ E(ℂ) の位数 l の部分群を Λ ⊆ Λ′ 指数 l に対応させる") 15,
    .implicitStep
      ("☆mathlib の在庫(2026-08-29 に測定): ある——℘・℘′・g₂・g₃・G n・周期性・偶奇・" ++
       "deriv ℘ = ℘′・meromorphic_weierstrassP・order_weierstrassP・derivWeierstrassP_sq・" ++
       "weierstrassPExcept。無い——加法定理、群同型 ℂ/Λ ≅ E(ℂ)、楕円関数の Liouville" ++
       "(★本プロジェクトが第 598 で建てた)") 13,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, VI.3.6(℘ の加法定理)"
      (.absent "mathlib の Analysis/SpecialFunctions/Elliptic/Weierstrass.lean は ℘ の理論を 1080 行ぶん持つが加法定理は無い(2026-08-29 に全宣言名を確認)") 13,
    .implicitStep
      ("★★★★★★2026-08-29(第 612-614)で z ≡ −w の側の核も取れた: " ++
       "Found/GenEll/Uniformization.lean の addQ(= 2(℘(t−w)−℘(w)) − t(℘′(t−w)−℘′(w)))が " ++
       "t = 0 で 3 位の零点であること(three_le_analyticOrderAt_addQ)。" ++
       "★q(0) = 0(℘ は偶)、q′(0) = 0(℘′ は奇)、" ++
       "q″(t) = −12t·℘(t−w)·℘′(t−w) なので q″(0) = 0" ++
       "——既存の deriv_derivWeierstrassP(deriv ℘′ = 6℘² − g₂/2、" ++
       "Found/GenEll/WeierstrassODE.lean)がそのまま効いた") 13,
    .implicitStep
      ("☆★**訂正(2026-08-29、第 615)**: 第 606・第 611 で「残るのは Λ と −w+Λ の " ++
       "2 か所の極だけ」と書いたが、**それは足りない**。" ++
       "★F_w は {z : ℘(z) = ℘(w)} の点でも極を持ちうるので、" ++
       "Liouville を当てるには **℘(z) = ℘(w) ⟹ z ≡ ±w (mod Λ)** が要る。" ++
       "☆これは「℘ は各値をちょうど 2 回取る」(位数 2 の楕円関数)であり、" ++
       "★★古典的には**零点と極の個数の一致**(偏角の原理・留数定理)から出る。" ++
       "☆mathlib は Meromorphic/Divisor.lean(MeromorphicOn.divisor)を持つが、" ++
       "楕円関数の零点和 = 極和は無い(2026-08-29 に測定)") 21,
    .implicitStep
      ("★★これで加法定理と一様化の単射性が**同じ 1 つの事実**に帰着することが分かった: " ++
       "「℘ は各値をちょうど 2 回取る」。☆両者は絡み合っており、" ++
       "どちらか一方を先に取ることはできない。" ++
       "★**次の入口はそこ**——楕円関数の零点勘定である。" ++
       "☆道具の Liouville(第 598)・除去可能特異点(第 603)・" ++
       "解析的位数(mathlib の natCast_le_analyticOrderAt)はすべて揃っている") 21,
    .implicitStep
      ("★★2026-08-29(第 618)で零点勘定の第一の煉瓦が置けた: " ++
       "Found/GenEll/Uniformization.lean の elliptic_boundary_integral_zero" ++
       "(周期平行四辺形の境界積分は消える)。★機構は周期性だけで Cauchy は要らない" ++
       "——向かい合う辺が打ち消し合う") 13,
    .implicitStep
      ("★★★次の一手の候補は 2 つある。" ++
       "(A) 留数定理経路: ∮ f = 2πi·Σ res を平行四辺形で建て、第 618 と合わせて " ++
       "Σ res = 0、f = g′/g で偏角の原理。☆平行四辺形の輪郭変形(ホモトピー)が要り、" ++
       "mathlib の Cauchy は軸平行な長方形版なので直接は当たらない。" ++
       "(B) Jensen 経路: mathlib の MeromorphicOn.circleAverage_log_norm" ++
       "(Analysis/Complex/JensenFormula.lean)を大きな円に当てる。" ++
       "★℘ は格子の近傍を除いて有界なので circleAverage log‖℘−c‖ は O(1)、" ++
       "一方 Σ_{|z|<R} ord_z·log(R/|z|) は基本領域あたりの零点数を N とすると " ++
       "(N−2)·(πR²/covol) の程度で伸びる。★★したがって N = 2。" ++
       "☆こちらは mathlib の在庫にそのまま乗るが、評価の詰めが要る") 21,
    .implicitStep
      ("★★★★★★★★2026-08-29(第 622-627)で**極はすべて塞がった**。" ++
       "☆★訂正: 第 615 で「零点勘定(偏角の原理)が要る」と書いたが、" ++
       "**線型 2 階 ODE の一意性で迂回できた**:" ++
       "h ≔ ℘(·+a) − ℘ は h″ = 6(℘(·+a) + ℘)·h を満たし(第 622)、" ++
       "h(z₀) = h′(z₀) = 0 なら解析的位数の算術で近傍で 0(第 623)、" ++
       "一致の定理と ℘ の極から a ∈ Λ(第 624)——★一様化は単射") 17,
    .implicitStep
      ("★★★F_w の極の一覧(2026-08-29 時点、すべて Found/GenEll/Uniformization.lean):" ++
       "Λ —— 第 610(laurent_cancel で z² が約される、addDefect_eq_near)。" ++
       "z ≡ w —— ℘′(w) ≠ 0 なら比が有界(第 625 の derivWeierstrassP_eq_zero_iff)。" ++
       "z ≡ −w —— 第 627(addDefect_eq_nearNeg、u = tû と q = t³g で t² が約される)。" ++
       "その他 —— 第 624-625 の単射性で存在しない" ++
       "(℘(z₀)=℘(w) なら ℘′(z₀)=±℘′(w) で、いずれの場合も z₀ ≡ ±w)") 15,
    .implicitStep
      ("☆残るのは Lean 上の手当てだけである: mathlib の ℘ は格子点で" ++
       "**junk value**(正則部分の値)を取るので addDefect はそこで連続でない。" ++
       "★したがって Liouville(第 598、Differentiable を要求)に当てるには" ++
       "極を埋めた整関数 Ext を明示的に作る必要がある" ++
       "(Λ ∪ (±w+Λ) は可算閉離散なので、各点で Riemann の除去可能特異点" ++
       "——第 603 の wp_inv_differentiableAt_of_mem と同じ型——を当てる)。" ++
       "☆MeromorphicAt.analyticAt は ContinuousAt を要求するので直接は使えない") 8,
    .implicitStep
      ("★★★★★★2026-08-29(第 629-630)で**道具もすべて揃った**:" ++
       "analyticOrderAt_weierstrassP_sub_self(℘ − ℘(w) は w で 1 位の零点)、" ++
       "analyticAt_limUnder_of_eventuallyEq(除去可能特異点を埋めて整関数 Ext を作る)、" ++
       "analyticAt_limUnder_of_analyticAt。" ++
       "☆残るのは組み立てだけである: Ext ≔ fun z => limUnder (𝓝[≠] z) (addDefect P w) が" ++
       "各点で解析的であることを 4 通りに場合分けして示し" ++
       "(Λ は第 610、z ≡ w は第 629、z ≡ −w は第 627、その他は第 624-625 で解析的)、" ++
       "Λ-周期性(addDefect の定義から)と Ext 0 = 0(第 610 の addDefectNear_zero)を足して" ++
       "elliptic_liouville(第 598)に当てる") 8,
    .implicitStep
      ("★★★★★★★★★★2026-08-29(第 638-641)で**閉じた**。" ++
       "Found/GenEll/Uniformization.lean の weierstrassP_addition が実物である。" ++
       "★機構: Ext ≔ fun z => limUnder (𝓝[≠] z) (addDefect P w) は整で Λ-周期的" ++
       "(第 639)、Ext 0 = 0(第 640)、良い点で Ext = F_w(第 633)。" ++
       "★★極は 4 か所とも塞がっている: Λ(第 610・639)・z ≡ w(第 638)・" ++
       "z ≡ −w(第 637)・その他は単射性(第 624-625)で存在しない。" ++
       "☆★単射性は零点勘定(偏角の原理)を使わず、線型 2 階 ODE " ++
       "h″ = 6(℘(·+a)+℘)h の一意性(第 622-624)から出た" ++
       "——第 615 の見立て(「零点勘定が要る」)は迂回できた") 21,
    .otherPaper "GenEll" "Lemma 3.5(l-捩れの大域的な階数 1 の部分群)" 15 ]

end ABC3.Skeleton.GenEll
