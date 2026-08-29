import ABC3.Found.GenEll.Uniformization

/-!
# スケルトン —— **`℘` の加法定理**（`Skeleton`）☆★★★★★★★★**残り 1 本**

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
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice)
    (hne : P.weierstrassP z ≠ P.weierstrassP w) :
    P.weierstrassP (z + w)
      = ((P.derivWeierstrassP z - P.derivWeierstrassP w)
          / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4
        - P.weierstrassP z - P.weierstrassP w := by
  sorry

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
    .otherPaper "GenEll" "Lemma 3.5(l-捩れの大域的な階数 1 の部分群)" 15 ]

end ABC3.Skeleton.GenEll
