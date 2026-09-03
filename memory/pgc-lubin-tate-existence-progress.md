---
name: pgc-lubin-tate-existence-progress
description: pGC の局所類体論に要る Lubin-Tate 形式群の存在補題——帰納法1ステップの3部品(f側線形化・g側線形化・可解性)がsorry無しで揃った現状
metadata:
  type: project
---

pGC Proposition 1.2/Corollary 1.3/Proposition 2.1/Proposition 2.2/Theorem 4.2 の
`.needs` が要求する局所類体論の相互律は、Lubin-Tate 形式群の存在補題
(次数ごとの冪級数 `Φ : ℕ → MvPowerSeries (Fin 2) A` を再帰的に構成し、
関数等式 `Φ(g,g) = f(Φ)` を満たす極限を作る)を経由するのが標準的な道筋。
2026-09-04 時点でこの構成の**帰納法1ステップに要る3部品すべて**が
sorry無しで揃った:

1. **可除性**(`Found/PGC/LubinTateDivisibility.lean::residue_divides_R`):
   剰余項 `R_n := Φ(g,g)−f(Φ)` の剰余体への還元が0であること。任意の `Φ`
   について既存の `mvPowerSeries_pow_card_eq_expand` 1個から即座に出る
   ——「帰納的不変量が要る」という当初の見積りは誤りだった。
2. **f側の線形化**(`Found/PGC/LubinTateLinearization.lean::coeff_subst_linearize`):
   `f≡πX(mod deg2)` のとき `f(Φ+φ)−f(Φ)` が次数 `≤deg φ` の範囲で
   `π·φ` に一致する。`geom_sum₂_mul` の因数分解 + `MvPowerSeries.order` の
   次数勘定。
3. **g側の線形化**(`Found/PGC/LubinTateGLinearization.lean::coeff_subst_g_linearize`):
   `g≡πX(mod deg2)`・`φ` が次数 `n+1` の斉次式のとき `φ.subst(g,g)` が
   次数 `≤n+1` の範囲で `π^{n+1}•φ` に一致する。f側の道具
   (`order_pow_sub_pow_ge`)を `order(a-a')` を変数のまま持ち回れるよう
   一般化(`order_pow_sub_pow_ge'`)し、2変数の telescoping
   (`order_prod_pow_sub_prod_pow_ge`)に適用して得た。
4. さらに独立な基盤として `Found/PGC/ValuationRingDVR.lean::valuationRing_isDVR`
   ——任意の p 進局所体の付値環が `IsDiscreteValuationRing` であること
   (`Valued.integer.isDiscreteValuationRing_of_compactSpace` 経由、
   `IsCyclic (valueGroup)` を要求する重い経路を回避)。

**続報(同日、1ステップ全体が完成)**: 上の3部品を実際に組み合わせて、
「`Φ` の障害が次数 `n` まで消えている」→「`Φ+φ`(`φ` は次数 `n+1` の
斉次式)の障害が次数 `n+1` まで消えている」という**1ステップ全体**を
`Found/PGC/LubinTateStepAssembly.lean::exists_next_step` として sorry 無しで
完成させた。鍵になった簡略化: 可除性補題が返す解 `φ₀` は厳密に斉次とは
限らないが、`φ := homogeneousComponent(n+1) φ₀` と**取り出す**だけで
g側線形化が要求する厳密な斉次性と方程式 `(π−π^{n+1})•φ=R_n` の両方が
満たせる——`homogeneousComponent` が `R`-線型写像であることと、`R_n`
自身が既に次数 `n+1` の斉次式であること(不動点性)の2点だけから出た。

**続報(同日、近似列そのものが完成)**: `Found/PGC/LubinTateBaseCase.lean`
(出発点——加法的形式群 `X_0+X_1` が次数1まで解けていることを直接確認)と
`exists_next_step` を `Nat.rec` で繋ぎ、`Found/PGC/LubinTateSequence.lean::
ΦSeq` として「任意の `k` について次数 `k+1` まで関数等式を満たす近似
`Φ_k`」の列を sorry 無しで構成した。不変量に「障害が消えている」と
「定数項が0」の2条件を最初から組み込んだのが実装を数行に抑えた鍵。

**続報(同日、Lubin-Tate存在補題が完全に完成)**: `Found/PGC/
LubinTateLimit.lean` で近似列 `ΦSeq` の極限 `F` を実際に構成し、
**`F.subst(g,g)=f.subst(F)` が power series の等式として exact に成り立つ**
ことを sorry 無しで証明した(`LubinTateF_functional_equation`)。`F` は
`MvPowerSeries` が係数関数そのものであることを使い、各係数 `e` を
`ΦSeq(degree e)` から直接読み取る形で定義——極限操作や完備化を経由
しない。安定性(`F−ΦSeq m` の次数が `m+1` 以上)は、g側は
`MvPowerSeries.le_order_subst`、f側は既存の `coeff_subst_linearize` を
「`degree e<order φ` のとき線形化の誤差項が0になる」形で再利用して示し、
新しい一般補題を増やさずに済んだ。★★★これで**Lubin-Tate形式群の存在補題
(可除性・両側の線形化・1ステップの組み立て・近似列・極限)が全段に
わたって完全にsorry無しで完成**した——古典的な局所類体論の中核的構成の
一つが形式化されたことになる。

**残る作業**: この形式群の上に実際の局所類体論の相互律写像
(Frobenius元・Lubin-Tateの主定理)を構成する段——pGC の各項目
(Prop 1.2・Cor 1.3・Prop 2.1・Prop 2.2・Theorem 4.2)を閉じるには、
なお相互律写像そのものの構成・性質証明という大きな仕事が残っている。

詳細な発見の経緯は `ResearchPaper/blocked-leaves.json` の
`[pGC] Proposition 1.2 / ... —— 局所類体論の相互律` エントリの
`progress2026_09_04a`〜`j` に記録。[[padic-log-additivity-blocked]]・
[[pgc-ramification-naturality-gap]] は同じ blocked エントリの別方向の前進。
