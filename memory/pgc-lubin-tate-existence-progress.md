---
name: pgc-lubin-tate-existence-progress
description: pGC の局所類体論(Prop 1.2・Cor 1.3・Prop 2.1・Cor 3.3・Theorem 4.2 が共通に依存)に要る Lubin-Tate 相互律——★★★★★★★★★★★★★★★★★★★★ reciprocityMap の存在一意性・単数側全単射・群準同型は sorry 無しで完備。真に残る最後の柱は reciprocityMap の**全射性**(Galois が ψ_n の根に推移的に作用)——鍵となる irreducible_iteratedLubinTatePsi は既に確立済みで道筋は明確(AdjoinRoot 経由の共役構成 + IsAlgClosed.lift + Algebra.IsAlgebraic.bijective_of_isScalarTower')だが、K.closure に2つの Algebra 構造を与えようとして instance diamond に当たり未完成(型シノニムか AdjoinRoot を R に取る順序変更が次の一手)
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

**続報(同日、g=f特殊化で「fはF_fの自己準同型」を確立)**: `Found/PGC/
LubinTateFormalGroupLaw.lean`で`LubinTateF`を`g:=f`で特殊化した古典的な
形式群法則`formalGroupLaw`(F_f)を定義し、`formalGroupLaw_f_isEndomorphism`
(`F_f(f(X),f(Y))=f(F_f(X,Y))`)をsorry無しで示した——`LubinTateF_
functional_equation`の直接の特殊化。これはfが乗法元π∈𝒪_Kに対応する
自己準同型`[π]_{F_f}`を与える(𝒪_K→End(F_f)という環準同型の第一歩)、
という古典的な出発点。

**続報(同日、一意性補題が完成)**: `Found/PGC/LubinTateUniqueness.lean::
powerSeries_uniqueness`——古典的なLubin-Tateの一意性補題(Lubin-Tate
1965: 同じ次数1の係数を持ち同じ関数等式を満たす2つの冪級数は一致する)
をsorry無しで完成させた。1変数冪級数を「X_1に依存しない2変数冪級数」
として埋め込む`emb`写像を経由し、既存の2変数線形化(`coeff_subst_
linearize`)をそのまま再利用——新しい1変数専用の線形化を組み立てずに
済んだのが鍵。この補題はF_fが実際に形式群法則であること(単位元則
F_f(X,0)=X・結合律・可換律)を示す標準的な道具。

**続報(同日、単位元則F_f(X,0)=Xが完成)**: `Found/PGC/
LubinTateIdentityLaw.lean::formalGroupLaw_identity`——一意性補題を
実際に`formalGroupLaw`(F_f)に適用し、単位元則`F_f(X,0)=X`をsorry無しで
確立した。鍵は`PowerSeries.subst_def`(`PowerSeries.subst`が定数族による
`MvPowerSeries.subst`そのものであること)の発見——当初「finsumの
reindexingが要る」と見積もっていたが、これで不要になり、既存の
`MvPowerSeries.subst_comp_subst_apply`だけで閉じた。`ψ:=F_f(X,0)`が
`f`との関数等式を満たすこと(`psi_functional_equation`)を示し、`X`
(恒等射)も同じ等式を満たすことと合わせて一意性補題で`ψ=X`を結論する。

**続報(同日、対称な単位元則F_f(0,Y)=Yも完成)**: `Found/PGC/
LubinTateIdentityLawSymm.lean::formalGroupLaw_identity_left`——
`LubinTateIdentityLaw.lean`をX_0↔X_1で対称に反復し、`F_f(0,Y)=Y`も
sorry無しで確立した。`subst_family_comp_value`等の一般補題は
X_0/X_1によらないためそのまま再利用でき、2回目はほぼ一発で通った。
これで**両方の**単位元則(`F_f(X,0)=X`・`F_f(0,Y)=Y`)が完成。

**続報(2026-09-04、可換律F_f(X,Y)=F_f(Y,X)が完成)**: `Found/PGC/
LubinTateCommutativity.lean::formalGroupLaw_commutative`——1変数一意性
補題を2変数に一般化した`mvPowerSeries_uniqueness`と、X_0↔X_1の入れ替え
`swap`が関数等式・定数項0・次数1係数を保つことを示す一連の補題を組み、
`swap(F_f)`が`F_f`と同じ次数1係数(`F_f`の次数1部分がX_0+X_1で対称
だから)・同じ関数等式を満たすことから`F_f=swap(F_f)`を導いた。新規の
一般化`coeff_subst_g_linearize_order`(g側線形化を「次数≥n+1」へ拡張)
が2変数一意性補題の帰納法で必要になった。これで**単位元則2つ・可換律の
3性質が全てsorry無しで揃った**。

**続報(2026-09-04、結合律への布石——σ一般化の実測)**: 結合律
`F_f(F_f(X,Y),Z)=F_f(X,F_f(Y,Z))`は3変数の一意性補題(可換律で作った
`mvPowerSeries_uniqueness`のFin 3版)を要すると見て、その土台となる
`hasSubst_g_subst_X`(`Found/PGC/LubinTateGLinearization.lean`)と
`coeff_subst_linearize`(`Found/PGC/LubinTateLinearization.lean`)を
`Fin 2`から任意の添字型`{σ:Type*}`へ一般化した(前者は`[Finite σ]`も
追加)——両方とも証明が`Fin 2`に一切依存しておらず(`order_pow_sub_pow_ge`
は既に`{A σ}`一般だった)、シグネチャを広げるだけで既存の全呼び出し
箇所が無変更で通る後方互換な一般化だった。★一方`coeff_subst_g_linearize`
自身(`Finsupp.prod`を`Fin.prod_univ_two`で2項に展開する箇所)は
`Fin 2`固有の構造を使っており、3変数版(`Fin.prod_univ_three`)は
別途専用の補題が要ることも判明——次に結合律に戻るときの具体的な
出発点はここ(`coeff_subst_g_linearize`のFin 3版、`coeff_a_diff_order`・
`coeff_ad_eq_of_degree`も同様に3項版が要る)。

**続報(2026-09-04、一意性補題を任意有限変数へ一般化——結合律の土台完成)**:
`mvPowerSeries_uniqueness`(可換律のFin 2固定2変数版)の証明を精査し、
Fin 2への依存が「`Finsupp.prod`を`Fin.prod_univ_two`で2項展開する
箇所」だけだったと判明。その展開を「有限集合上の望遠鏡和」
(`order_prod_pow_sub_prod_pow_ge_finset`、`Finset.induction_on`で
1項ずつ剥がす)に置き換え、`Found/PGC/LubinTateGeneralUniqueness.lean::
mvPowerSeries_uniqueness_general`として**任意の有限添字型σ(Fin 3含む)**
への一般化をsorry無しで完成させた。結合律`F_f(F_f(X,Y),Z)=F_f(X,F_f(Y,Z))`
はσ:=Fin 3でこれを`G:=F_f(F_f(X,Y),Z)`・`H:=F_f(X,F_f(Y,Z))`に適用
すれば示せる見込み——残るのは次数1係数の一致とG・Hそれぞれの関数等式
(F_f自身の関数等式の2回合成)を示す段。
★診断: この一般化作業中、既存の`hasSubst_g_subst_X`呼び出し3箇所
(型注釈なしでσの遅延推論に依存)が「typeclass instance problem is
stuck」で再現性をもって落ちる現象を発見し、`(σ := Fin 2)`の明示注釈で
解消した——原因(前回のフルビルド成功時は通っていた)は未特定のまま
(並行セッションと`.lake/build`キャッシュを共有する環境固有の可能性)。

**続報(2026-09-04、結合律 F_f(F_f(X,Y),Z)=F_f(X,F_f(Y,Z)) が完成)**:
`Found/PGC/LubinTateAssociativity.lean::formalGroupLaw_associative`——
`G:=F_f(F_f(X,Y),Z)`・`H:=F_f(X,F_f(Y,Z))`を3変数の`MvPowerSeries`
として構成し、3変数一意性補題(`mvPowerSeries_uniqueness_general`、
`σ:=Fin 3`)を適用して`G=H`を確立した。鍵になった2つの新しい一般補題:
(1)`subst_preserves_functional_equation`——可換律の`swap_preserves_
functional_equation`を「`swap`という特定の代入」から「任意の代入族」
へ一般化したもの、(2)`coeff_single_subst_degree_one`——代入の次数1
係数についての一般公式(次数勘定で「次数≥2の項は効かない」ことから
出る)。この2つでG・Hそれぞれの関数等式・次数1係数(いずれも全て`1`)
がほぼ機械的に出た。★★★★★これで**古典的なLubin-Tate形式群法則の3性質
(単位元則2つ・可換律・結合律)が全てsorry無しで揃った**——存在補題・
自己準同型性・一意性補題と合わせ、この分野の形式化がひとつの節目に
到達した。

**続報(2026-09-04、𝒪_K作用への拡張・節目(1)の可除性が完成)**:
残り工程(1)([a]_f:f∘[a]_f=[a]_f∘fを満たす1変数冪級数の存在)の
可除性段を`Found/PGC/LubinTateEndoDivisibility.lean::residue_divides_
R_endo`としてsorry無しで確立した。2変数版(`residue_divides_R`)と
全く同じFrobenius恒等式の議論が、`PowerSeries R`が`MvPowerSeries
Unit R`のabbrevであることを利用してほぼそのまま移植できた——
`PowerSeries.expand`/`.map`/`.constantCoeff`はmathlibで文字通り
`MvPowerSeries`の対応物として定義されているため、既存の
`mvPowerSeries_pow_card_eq_expand`(σについて一般)を`σ:=Unit`で
そのまま適用できた。残るは両側の線形化・1ステップの組み立て・
近似列・極限という、2変数版と同型だが1変数へ作り直す工程。

**続報(2026-09-04、𝒪_K作用への拡張・両側の線形化が完成)**:
`Found/PGC/LubinTateEndoLinearization.lean`で両側の線形化をsorry無しで
確立した。f側は既存の一般補題`coeff_subst_linearize`(§9-1336でσ一般化
済み)を`σ:=Unit`で直接適用するだけ(1行)。g側は2変数版より**単純に
なる**ことを発見——`g≡πX(mod deg2)`のとき`g^{n+1}`の次数`n+1`の係数が
`π^{n+1}`という事実(`coeff_pow_self_eq_pow`)だけで足り、既にσ一般化
済みの`order_pow_sub_pow_ge'`を直接適用して出た。途中「`PowerSeries.
order`と`MvPowerSeries.order`はdefeq」という前回の見積りが誤りだったと
判明したが、mathlibの橋渡し補題`PowerSeries.order_eq_order`で解決した。
★これで節目(1)の1変数存在定理に要る3部品(可除性・f側線形化・g側
線形化)が全て揃った。

**続報(2026-09-04、次数ごとの1ステップの組み立てが完成)**:
`Found/PGC/LubinTateEndoStepAssembly.lean::exists_next_step_endo`——
3部品(可除性・f側線形化・g側線形化)を組み合わせ、「φの障害が次数≤n
で消えている」→「φ+c•X^{n+1}(cはスカラー)の障害が次数≤n+1で消える」
という1ステップ全体をsorry無しで完成させた。2変数版よりも簡略化された
——次数n+1の「斉次式」が単なるスカラーcなので、`homogeneousComponent`
を「取り出す」操作が一切不要だった。スカラー方程式`c(π-π^{n+1})=-r`は
`exists_step_solution`と全く同じ機構(単数の逆元)で解けた。

**続報(2026-09-04、出発点・近似列が完成——手続き部分すべて完了)**:
`Found/PGC/LubinTateEndoBaseCase.lean`(出発点`a•X`が次数≤1まで解けて
いることを確認)と`Found/PGC/LubinTateEndoSequence.lean::φSeq_endo`
(`Nat.rec`による近似列)をsorry無しで完成させた。2変数版`ΦSeq`と同じ
構造だが、各段の補正がスカラー`c`そのものなので`homogeneousComponent`
の抽出が不要という簡略化。★これで節目(1)の1変数存在定理は**次数ごとの
手続き部分すべて**(可除性・両側線形化・1ステップ・出発点・近似列)が
sorry無しで完成した。

**続報(2026-09-04、★★★★★★★★節目(1)完成——𝒪_K作用への拡張の1変数
存在定理そのものがsorry無しで完成)**: `Found/PGC/LubinTateEndoLimit.lean::
LubinTateEndo_functional_equation`——近似列`φSeq_endo`の極限
`LubinTateEndo`を構成し、関数等式`f.subst(LubinTateEndo)=LubinTateEndo.
subst(g)`(`f(φ(X))=φ(g(X))`)をexactに満たすことを証明した。2変数版
`LubinTateF`(存在補題本体)と完全に並行する議論(係数の安定性・g側/f側
の安定性・最終的な関数等式)。副産物として`constantCoeff_LubinTateEndo=0`・
`coeff_one_LubinTateEndo=a`(出発点のパラメータ`a`が実際に次数1係数として
実現されること)も確立した。★★★★★★★★これで**𝒪_K作用への拡張(節目1、
a↦[a]_fの存在)は1変数存在定理そのものが完全にsorry無しで完成した**
——可除性・両側線形化・1ステップ・出発点・近似列・極限という全工程が
揃った。前回(§9-1340)確立した通り、自己準同型性`F_f([a]X,[a]Y)=
[a]F_f(X,Y)`(節目2)はこの`LubinTateEndo`と`mvPowerSeries_uniqueness`
から既存の道具だけで形式的に出る見込み。

**続報(2026-09-04、★★★★★★★★節目(2)完成——`[a]_f`が`F_f`の自己準同型
であることがsorry無しで確立)**: `Found/PGC/LubinTateActionEndomorphism.
lean::formalGroupLaw_endomorphism_of_action`——`F_f([a]X,[a]Y)=[a]F_f(X,Y)`
を確立した。§9-1340で見立てた通り、`LubinTateEndo_functional_equation`
(g:=f特殊化=`LubinTateAction`)と可換律・結合律で確立済みの道具だけから
形式的に出た——新しい次数ごとの構成は一切不要だった。`Ψ_1:=F_f([a]X,
[a]Y)`・`Ψ_2:=[a]F_f(X,Y)`が同じ次数1係数(全て`a`)・同じ関数等式を
満たすことから、既存の2変数一意性補題で`Ψ_1=Ψ_2`を結論した。新たに得た
鍵の道具`subst_comp_subst_single_gen`(1変数どうしの代入の合成規則)は
今後も再利用できる見込み。★★★★★★★★これで**節目(1)(1変数存在定理)・
節目(2)(自己準同型性)の両方が完成**し、𝒪_K全体への作用拡張という
大きな節目が実質完了に近づいた。

**続報(2026-09-04、環準同型の乗法側`[ab]_f=[a]_f∘[b]_f`が完成)**:
`Found/PGC/LubinTateActionMul.lean::LubinTateAction_comp`——`α:=[ab]_f`・
`β:=[a]_f∘[b]_f`が同じ次数1係数(いずれも`ab`)・同じ関数等式(`f`と
可換)を満たすことから、**1変数の**一意性補題(`powerSeries_uniqueness`、
`[a]_f`自身が1変数の冪級数なので2変数版ではなくこちら)で`α=β`を
結論した。`β`の関数等式は`subst_comp_subst_single_gen`(節目2で確立
済み)を4回組み合わせるだけで出た。新規の補助補題`coeff_one_subst_
1var`(代入の次数1係数は係数の積になる、1変数版)も確立した。

**続報(2026-09-04、★★★★★★★環準同型の加法側`[a+b]_f=F_f([a]_f,[b]_f)`
が完成——節目2b全体が完成)**: `Found/PGC/LubinTateActionAdd.lean::
LubinTateAction_add`。当初「2変数の一意性補題+もう1段一般的な補題が
要る」と見立てていたが、実際には`δ:=F_f([a]_f,[b]_f)`(`family:=fun i=>
if i=0 then [a]_f else [b]_f`を`F_f`へ代入したもの)は`Unit`添字、
すなわち正真正銘**1変数**の冪級数だったので、乗法側と**同じ**1変数
一意性補題(`powerSeries_uniqueness`)で足りた。新規の一般補題
`subst_value_comp_family_general`(「1つの値の代入」と「代入の族」の
合成則のもう半分——`subst_family_comp_value_general`の`c(v(p))=p(c(v))`
に対し、こちらは`v(c(Φ))=Φ(fun i=>v(c_i))`の形、`PowerSeries.subst_def`
で両側から`MvPowerSeries.subst_comp_subst_apply`に橋渡しするだけ)が
鍵で、これと`subst_preserves_functional_equation`(結合律・節目2で
確立済み)・各成分の可換性(`LubinTateAction_functional_equation`)を
組み合わせて`δ`の関数等式を出した。新しい次数ごとの構成は一切不要
だった。★★★★★★★これで**`a↦[a]_f : 𝒪_K → End(F_f)`が環準同型である
こと(節目2b)が乗法側・加法側ともに完成**し、局所類体論の相互写像
へ向けた6節点ロードマップの節点(1)(2)(2b)がすべてsorry無しで確立
された。

**続報(2026-09-04、節目(3)への第一歩——`[π^n]_f` は `f` の `n` 回
自己合成)**: `Found/PGC/LubinTateActionPiPow.lean`。`LubinTateAction_
pi_eq_f`(`[π]_f=f`自身、`f`が「次数1係数π・自明な関数等式」を満たす
ことから一意性で直ちに従う——ほとんど同語反復)・`LubinTateAction_
one_eq_X`(`[1]_f=X`)・`LubinTateAction_pi_pow`(`[π^n]_f`は`f`の`n`回
自己合成`iteratedLubinTate`、`Nat.rec`)を確立。3つとも乗法側の環準
同型性(`LubinTateAction_comp`)・可換性(`LubinTateAction_functional_
equation`)という既存の道具だけから帰納法で出た。★これで節目(3)
(`Λ_n:=ker[π^n]_f`の構成)へ向けて、`[π^n]_f`が具体的に何であるか
(`f`のn回自己合成)が確立された。

**発見(mathlibにWeierstrass標準分解が存在)**: `node tools/decl-index.mjs
--mathlib`での探索で、`RingTheory/PowerSeries/WeierstrassPreparation.lean`
に`PowerSeries.weierstrassDistinguished`・`Polynomial.IsDistinguishedAt`・
`PowerSeries.IsWeierstrassFactorization`等、完備局所環上の冪級数の
Weierstrass標準分解(冪級数を「単数×distinguished多項式」に分解する
古典的定理)が**既に整備されている**ことを確認した(2026-08-14時点の
`blocked-leaves.json`項目8の調査ではこの発見はまだ無かった)。
`Polynomial.IsEisensteinAt.irreducible`(Eisenstein判定法)も存在。
★これは節目(3)で`[π^n]_f`(無限冪級数)から次数`q^n`の実際の多項式
(distinguished多項式、Eisenstein判定で既約性も出る)を取り出す際の
**大きな加速要因**になる見込み——Weierstrass標準分解そのものを自前で
構築する必要がなくなった。次にここへ戻るときは、まず`[π^n]_f`の
`IsWeierstrassFactorization`(または`IsDistinguishedAt`)条件の充足を
示すところから着手する。

**続報(2026-09-04、節目(3)への第二歩——`[π^n]_f` は mod `π` で
`X^(q^n)`)**: `Found/PGC/LubinTateActionPiPow.lean::iteratedLubinTate_
map_residue`——`n=1`の場合が仮定`hf`そのものである、古典的な基本事実。
`PowerSeries.map_subst`(写像は代入と可換)を軸に`f`のmod `π`での姿を
`n`回適用する帰納法で出た。途中、`PowerSeries.subst`の元の環が
`map(residue A)`経由で`ResidueField A`になっていることに気づかず`A`
のままと誤認し、`rw`が「パターンが見つからない」で数回落ちる罠に
はまった——`tools/lean-idioms.md`の「instances透明度」系ではなく、
単純な型の取り違えだったと後で判明(ゴールを注意深く読み直す以外に
近道は無かった)。プロジェクト内の既存補題`substXpow_eq_pow`
(`subst v (X^q)=v^q`)も再利用できた。

**続報(2026-09-04、★★★★★★★★節目(3)への第三歩——`[π^n]_f` の
Weierstrass標準分解、次数`q^n`のdistinguished多項式)**: `Found/PGC/
LubinTateActionWeierstrass.lean`。`A`が`π`-進完備(`[IsAdicComplete
(maximalIdeal A) A]`、古典的Lubin-Tate理論の標準設定`A=𝒪_K`、この節目
で新たに必要になった仮定)のもとで、mathlibのWeierstrass標準分解定理
をそのまま適用: `iteratedLubinTate_eq_distinguished_mul_unit`
(`[π^n]_f=D_n·U_n`)・`isDistinguishedAt_iteratedLubinTateDistinguished`
(`D_n`はdistinguished)・`isUnit_iteratedLubinTateUnit`(`U_n`は単元)・
★`natDegree_iteratedLubinTateDistinguished`(`D_n`の次数は`q^n`——古典的
理論で`|Λ_n|=q^n`の由来になる事実)。前提は`iteratedLubinTate_map_
residue`から従う非零性のみで、Weierstrass標準分解そのものはmathlibの
`exists_isWeierstrassFactorization`等をそのまま使い、自前で構築する
必要は一切無かった——前回の「大きな加速要因」という見立てが的中した。
次数の計算は`Polynomial.IsDistinguishedAt.coe_natDegree_eq_order_map`
をℕ∞値のまま比較する形で行い(`ENat.lift`の依存書き換え問題を回避)。
★これで節目(3)の核——`[π^n]_f`から次数`q^n`の具体的な多項式`D_n`を
取り出すこと——が完成し、次は`D_n`の根の集合を`Λ_n`として定義し、
`𝒪_K`加群構造・`K(Λ_n)/K`が完全分岐であることへ進む見込み。

**続報(2026-09-04、節目(3)への第四歩——`D_n`の定数項は正確に0、
`X∣D_n`)**: `[π^n]_f=D_n・U_n`の定数項を比較するだけの初等的な議論
(`[π^n]_f`自身の定数項が0であることと`U_n`の定数項が単元であること
から)で`D_n`の定数項が正確に0であることを示した
(`iteratedLubinTateDistinguished_coeff_zero`・`X_dvd_iteratedLubinTate
Distinguished`)。これは古典的な`D_n(X)=X・φ_n(X)`(`φ_n`が原始
`π^n`-捩れ点を統べる多項式)という構造の出発点。

**続報(2026-09-04、節目(3)への第五歩——`φ_n:=D_n/X`の次数は`q^n-1`)**:
`D_n=X・φ_n`(`X∣D_n`から選択関数で抽出、`iteratedLubinTateDistinguished_
eq_X_mul_primitive`)を確立し、`D_n`・`φ_n`ともモニックであること
(`monic_iteratedLubinTateDistinguished`・`monic_iteratedLubinTate
Primitive`)、★`φ_n`の次数は`q^n-1`(`natDegree_iteratedLubinTate
Primitive`、`natDegree_mul`+両辺非零)を示した。古典的理論で`φ_n`の根
が非零な`π^n`-捩れ点`Λ_n\{0}`を統べることの由来になる事実。

**続報(2026-09-04、節目(3)への第六歩——`[π^b]_f∣[π^{a+b}]_f`、
`Λ_n⊆Λ_m`の由来)**: `iteratedLubinTate_add`(`[π^{a+b}]_f=[π^a]_f∘
[π^b]_f`)・★`iteratedLubinTate_dvd_iteratedLubinTate_add`
(`[π^b]_f∣[π^{a+b}]_f`)・`iteratedLubinTate_dvd_of_le`(`n≤m`なら
`[π^n]_f∣[π^m]_f`)。Weierstrass標準分解と独立な、より基礎的な冪級数
レベルの整除関係で、`[π^a]_f`の定数項が0であることと`subst`が積を
保つこと(`PowerSeries.subst_mul`)を組み合わせるだけで出た——古典的
理論で「`n≤m`ならば`Λ_n⊆Λ_m`」の由来になる事実。

**続報(2026-09-04、★★★★★★★★★節目(3)への第七歩——`D_n`・`φ_n`の根は
代数閉体で重複度込み`q^n`・`q^n-1`個)**: `iteratedLubinTateAlgClosure
A:=AlgebraicClosure (FractionRing A)`(根を数える舞台、`A→`この
代数閉体への`algebraMap`はmathlibのインスタンス解決で自動的に得られた)
を固定し、`D_n`・`φ_n`をそこへ写した根の多重集合
(`iteratedLubinTateDistinguishedRoots`・`iteratedLubinTatePrimitive
Roots`)を定義。★`card_iteratedLubinTateDistinguishedRoots`
(`D_n`の根は重複度込みでちょうど`q^n`個)・`card_iteratedLubinTate
PrimitiveRoots`(`φ_n`の根は重複度込みで`q^n-1`個)——古典的理論の
`|Λ_n|=q^n`に対応する事実を、分離性・既約性を一切使わずに確立した。
鍵はmathlibの`IsAlgClosed.card_roots_eq_natDegree`(代数閉体上の
多項式の根の個数は次数に一致、0多項式も含め無条件)と`Polynomial.
Monic.natDegree_map`(モニック多項式の次数は任意の環準同型で保たれる)
——この2つと既存の次数公式を組み合わせるだけで、新しい構成は一切
不要だった。分離性(根が相異なること、`Λ_n`が真の`q^n`元集合になる
ために必要)は別途の課題として残る。

**続報(2026-09-04、★★★★★★★★★節目(3)への第八歩——`n≤m`ならば`D_n∣D_m`、
`Λ_n⊆Λ_m`の多項式版)**: `iteratedLubinTate_dvd_map_residue`(商`r`の
mod π像は`X^(q^m-q^n)`)を経て、★`iteratedLubinTateDistinguished_dvd_
of_le`(`n≤m`ならば`D_n∣D_m`)を確立した——商`r`自身をWeierstrass分解
し(`r=D'・U'`)、`(D_n・D')・(U_n・U')`が`[π^m]_f`のもう1つの
Weierstrass分解であることを`Polynomial.IsDistinguishedAt.mul`
(distinguished多項式の積はdistinguished)で確認、一意性
(`IsWeierstrassFactorization.unique`)から`D_m=D_n・D'`を結論する
という、トポロジー・分離性を一切使わない純粋代数的な議論。
`iteratedLubinTateDistinguishedRoots_le_of_le`で代数閉体上の根の
多重集合版(`Λ_n⊆Λ_m`)も得た。これで`Λ_n`が`n`について単調に増大
する族を成すことが確立された。

**続報(2026-09-04、`0∈Λ_n`の整合性確認+次の壁の見立て)**:
`zero_mem_iteratedLubinTateDistinguishedRoots`(`D_n(0)=0`の言い換え)。
★節目(3)の次の壁は**分離性**(`Λ_n`が真に`q^n`個の相異なる元を持つ
こと)——古典的な証明は導関数の付値評価(ramification-theoretic)を
要し、本セッションでは未着手の新規の数学的内容と見立てた。副産物として
mathlibに`RingTheory/PowerSeries/Evaluation.lean`(`PowerSeries.HasEval`・
`aeval`・`eval₂`、位相環での冪級数の「点での評価」)があることを発見
——将来`𝒪_K`加群作用を`Λ_n`の点に実際に作用させる際の加速要因になり
うるが、代数閉包に適切な位相・付値を載せる作業が別途必要。代替案
(`PowerSeries.weierstrassMod`経由の`A[X]/(D_n)`への純代数的な作用)も
検討したが、well-defined性の証明が新たな作業として残る。詳細は
`ResearchPaper/blocked-leaves.json`の`progress2026_09_04v`。

**続報(2026-09-04、★★★★★★★★★★節目(3)への核心的定理——`φ_1`は既約、
`K(α)/K`は次数`q-1`の完全分岐拡大)**: 直前の見立て(「分離性は導関数の
付値評価を要する新規の数学的内容」)を再検討したところ、`n=1`の場合
限定でEisenstein判定法が直接使えることに気づいた。鍵になったのは
`coeff_one_iteratedLubinTate`(`[π^n]_f`の次数1係数は`π^n`——`f'(0)=π`
の連鎖律を`coeff_one_subst_1var`で`n`回)——これと`D_n=X・φ_n`の
係数シフトを組み合わせ、`iteratedLubinTateDistinguished_coeff_one_mul`
(`D_n.coeff1・U_n定数項=π^n`)を経由して、★`iteratedLubinTatePrimitive_
coeff_zero_notMem_sq`(`φ_1`の定数項は`𝔪²`に属さない——属すと仮定すると
`π`自身が単元になる背理法)を確立した。`φ_n`が`D_n`(弱Eisenstein)の
係数シフトから自動的に弱Eisensteinであることと合わせ、
★★★★★★★★★★`irreducible_iteratedLubinTatePrimitive_one`
(`φ_1`は既約)を`Polynomial.IsEisensteinAt.irreducible`で結論した。
これは古典的なLubin-Tate理論の核心的な帰結——原始`π`-捩れ点を1つ
添加した拡大`K(α)/K`は次数`q-1`の**完全分岐拡大(Eisenstein拡大)**
である——の多項式版。`n≥2`の`φ_n`は同じ議論が通らない(定数項が
`π^n∈𝔪²`になってしまう、`n=1`固有の現象)ため、一般の`n`への拡張は
別の議論(原始捩れ点だけを分離する多項式`φ_n/φ_{n-1}`のEisenstein性
など)が必要——引き続きの課題として残る。

**続報(2026-09-04、★★★★★★★★★★★節目(3)の核心定理を任意の`n≥1`へ
一般化——`ψ_n:=D_n/D_{n-1}`は既約)**: 直前の`φ_1`(`n=1`限定)の議論を
再検討し、実は`ψ_n:=D_n/D_{n-1}`(原始`π^n`-捩れ点`Λ_n\Λ_{n-1}`を
統べる、次数`q^n-q^{n-1}`の多項式——`φ_1`とは異なる、より正確な
古典的対象)が**任意の`n≥1`**で既約になることに気づいた。鍵は
`[π^n]_f=[π^{n-1}]_f・r_n`という分解で、`r_n`の定数項が**ちょうど`π`**
になること(`n=1`の`φ_1`の議論での`π^n・単元`より単純——`f=X・h`の
`h(0)=f'(0)=π`が`n`に依らず直接効く)。`r_n`をWeierstrass分解し
(`r_n=ψ_n・U'_n`)、`(D_{n-1}・ψ_n)・(U_{n-1}・U'_n)`が`[π^n]_f`の
もう1つのWeierstrass分解であることから一意性で`D_n=D_{n-1}・ψ_n`を
得て(`iteratedLubinTateDistinguished_dvd_of_le`と同じ手法のa=1特化)、
`ψ_n`の定数項が`𝔪²`に属さないことを`n=1`と全く同じ背理法で示し、
★★★★★★★★★★★`irreducible_iteratedLubinTatePsi`(`ψ_n`は既約、任意の
`n≥1`)を`Found/PGC/LubinTateActionPsi.lean`(新規)で確立した。これは
`irreducible_iteratedLubinTatePrimitive_one`(`n=1`限定)の完全な
一般化——古典的なLubin-Tate理論の核心定理そのもの(原始`π^n`-捩れ点を
1つ添加した拡大`K(α)/K(Λ_{n-1})`は次数`q^n-q^{n-1}`の完全分岐拡大)
が任意の`n`で確立された。

**続報(2026-09-04、★★★★★★★★★ψ_nの根を添加した体は次数`q^n-q^{n-1}`)**:
`Found/PGC/LubinTateActionPsiField.lean`(新規)。`irreducible_
iteratedLubinTatePsi`(`A`上の既約性)をGaussの補題
(`Polynomial.IsPrimitive.irreducible_iff_irreducible_map_fraction_
map`)で`FractionRing A`上の既約性へ橋渡しし、`AdjoinRoot`の標準基底
(`PowerBasis`)の次元公式と組み合わせて、★`finrank_adjoinRoot_
iteratedLubinTatePsi`(`AdjoinRoot(ψ_n)`の`Frac(A)`上の次元は
`q^n-q^{n-1}`)を確立した——これで「原始`π^n`-捩れ点を1つ添加した
拡大の次数」という古典的な結論が**体の言葉**で得られた。★逸脱の記録:
Gaussの補題は`[UniqueFactorizationMonoid A]`を要求する——古典的
Lubin-Tate理論では`A=𝒪_K`はDVRなので自動的に成り立つ性質だが、
これまでの議論には含まれていなかった。既存の定理には一切触れず、
このファイルの新しい定理だけに使うので後続への影響は無いことを
確認済み(CLAUDE.md 逸脱の記録)。

**続報(2026-09-04、★★★★★★★★★★分離性——ψ_nの根は混標数の場合に
相異なる、`Λ_n`ギャップの解消)**: `[CharZero A]`を追加(古典的
Lubin-Tate理論の混標数の場合——`K`は`ℚ_p`の有限次拡大——に対応、
CLAUDE.md逸脱の記録)。標数0では既約多項式は自動的に分離的
(`Irreducible.separable`)なので、`separable_iteratedLubinTatePsi_
map_fractionRing`(`ψ_n`は`FractionRing A`上で分離的)を経て、
★`nodup_iteratedLubinTatePsiRoots`(`ψ_n`の根は互いに相異なる——
分離性を`iteratedLubinTateAlgClosure A`へ持ち上げ`Polynomial.nodup_
roots`を適用)を確立した。既出の`card_iteratedLubinTatePsiRoots`
(個数は`q^n-q^{n-1}`)と合わせて、混標数の場合には「原始`π^n`-捩れ点
が真に`q^n-q^{n-1}`個の**相異なる**元からなる」という古典的な
Lubin-Tate理論の帰結が完成した——数セッション前に「分離性は導関数の
付値評価を要する新規の数学的内容」と見立てていたギャップが、
混標数という自然な設定限定で完全に解消された。

**続報(2026-09-04、`ψ_n`の分解体は`Frac(A)`上Galois)**:
`isGalois_iteratedLubinTatePsi_splittingField`——分離多項式の分解体は
自動的にGalois(`IsGalois.of_separable_splitting_field`)という一般
事実を`separable_iteratedLubinTatePsi_map_fractionRing`に適用する
だけの軽量な追加。この分解体が`AdjoinRoot ψ_n`(1つの根を添加する
だけの体)と一致すること(𝒪_K作用による他の根への言及を要する、
古典的Lubin-Tate理論の実質的な内容)は別途の課題として残る。

**続報(2026-09-04、★★★★★★★★★大きな節目——実在のp進局所体でLubin-Tate
理論一式がすべて空虚でないと確定)**: 分離性の議論(付値拡大)で行き詰まった
ため、方針を転換し、これまで抽象的な`A`について構築してきたLubin-Tate
理論一式が**実在の対象に適用できるか**を確認した。3つの発見があった:

①`Found/PGC/ValuationRingDVR.lean::valuationRing_isDVR`(既出、任意の
p進局所体`K:PAdicLocalField p`の付値環`𝒪[K.carrier]`は離散付値環)を
起点に、★`Found/PGC/ValuationRingComplete.lean`(新規)で
`IsAdicComplete (maximalIdeal) (𝒪[K.carrier])`(コンパクト性→
`CompleteSpace`+距離位相が`maximalIdeal`-進位相と一致すること`IsAdic`
の確認、`Irreducible.maximalIdeal_pow_eq_closedBall_pow`・
`IsUltrametricDist.isOpen_closedBall`・`exists_pow_lt_of_lt_one`を
組み合わせる)・`UniqueFactorizationMonoid`(DVR→PID→UFD)・`CharZero`
(ℚ_pからの単射)を確立した。

②`Found/PGC/LocalFieldNorm.lean`(既出)が既に`Fintype (ResidueField)`・
`ExpChar ... p`・`hq`(`residueCard_isPrimePow`)を提供していたと判明。

③`Found/PGC/LubinTateSeriesExists.lean`(新規)で、`f(X):=πX+X^q`が
`f≡πX(mod deg2)`・`f≡X^q(mod π)`を満たすこと(`exists_lubinTateSeries`)
——Lubin-Tate級数`f`の存在そのもの——を確立した。

★これで**このセッションのLubin-Tate理論一式(節目1・2・2b・3の全定理)
が実在のp進局所体`K:PAdicLocalField p`の整数環`𝒪[K.carrier]`に対して
そのまま適用できる**ことが検証された——原典[pGC]が前提する対象で
理論全体が空虚でないと確定した、大きなマイルストーン。

**続報(2026-09-04、★★★★★★★★★★★★大きな節目——Newton polygon論法
完成、ψ_n・ψ_m(n≠m)は共通根を持たない)**: 前回のspectralNorm発見
(`Found/PGC/LubinTatePsiNorm.lean`)を押し進め、古典的なLubin-Tate理論
のNewton polygon論法を完全に確立した。`ψ_n.coeff0`のノルムが`‖π‖`
(一定)であることと、`K`の代数閉包`K.closure`での`ψ_n`の**任意の**根
`x`について(既約性から`minpoly K.carrier x=ψ_n`となるので)
`spectralNorm K.carrier K.closure x = ‖π‖^(1/(q^n-q^{n-1}))`が成り立つ
ことを示し(`spectralNorm_root_iteratedLubinTatePsi`)、`0<‖π‖<1`
(`norm_pi_pos_lt_one`)なので指数`1/(q^n-q^{n-1})`が`n`ごとに異なる値
を与えること(`torsionDegree_ne`・`rpow_ne_rpow_of_base_lt_one`)から、
★★★★★★★★★★★★`no_common_root_iteratedLubinTatePsi`(`ψ_n`・`ψ_m`、
`n≠m`は共通根を持たない)を結論した——共通根があればそのスペクトル
ノルムが2つの異なる値に等しくなり矛盾、というNewton polygon論法
そのもの。★これで「異なる段の捩れ点は互いに素な多項式の根である」
という古典的Lubin-Tate理論の核心的な補題が完全に確立された。
`D_n`全体の分離性(異なる`ψ_i`たちの積が互いに素な既約多項式の積で
あることからsquarefreeが従う見込み)・正規性(`K(Λ_n)`が`ψ_n`の分解体
と一致すること)への道が開いた。

(記録: このファイル群はコミット時に共有git索引の競合でピア(第4章
担当)のコミット`f0646fd3`に混入したが、`git diff HEAD`で内容が完全に
一致することを確認済み——既にorigin/masterへpush済み。)

**続報(2026-09-04、★★★★★★★★★★ψ_n・ψ_mは互いに素+各々分離的——
`D_n`分離性の最終部品)**: `no_common_root_iteratedLubinTatePsi`
(Newton polygon)とは別の、より初等的な経路も見つけた——
`isCoprime_iteratedLubinTatePsi`(`ψ_n・ψ_m`、`n≠m`は`IsCoprime`)は
既約性+次数の違い(`torsionDegree_ne`)だけから直接出る(`ψ_n∣ψ_m`と
仮定すると同伴→モニックなので等しい→次数矛盾)。さらに
`separable_iteratedLubinTatePsi_map_carrier`(`ψ_n`は`K.carrier`上
分離的、混標数なので`Irreducible.separable`)を確立した。これで
`D_n=X・ψ_1・…・ψ_n`の各因子が分離的かつ互いに素であることの材料が
出揃った——`D_n`全体の分離性(squarefree性)への組み立ては、
`iteratedLubinTateDistinguished_eq_mul_psi`の反復適用と
`Polynomial.separable_prod`(または`Separable.mul`の反復)で得られる
見込み(組み立て自体は今後の課題)。

**続報(2026-09-04、`D_0=X`——`D_n=D_{n-1}・ψ_n`の基底段)**:
`iteratedLubinTateDistinguished_zero`(`[π^0]_f=[1]_f=X`自身が既に
distinguishedであることから一意性で結論)。`D_n=X・ψ_1・…・ψ_n`という
明示的な積表示を組み立てる際の基底段になる。★試みたが完成しなかった
こと: `D_n`全体の分離性(squarefree性)を`D_n=D_{n-1}・ψ_n`の帰納法で
組み立てようとしたが、各段で`IsCoprime D_{n-1} ψ_n`を示すには
「`D_{n-1}`の根は`{0}`∪`ψ_1,…,ψ_{n-1}`の根」という積表示そのものを
経由する必要があり、単純な帰納法だけでは`D_{n-1}`の根の特徴づけを
先に確立しないと回らない(軽い循環)——正直に見送り、`D_0=X`という
綺麗な基底段だけを確立して記録した。次に戻るときは、`D_n`の積表示
そのものを先に(`D_n=X・∏_{k=1}^n ψ_k`という形で)確立し、そこから
分離性を導く順序で組み立てるのが筋が良さそう。

**続報(2026-09-04、★★★★★★★★★★★★★大きな節目——`D_n`は分離的、
`Λ_n`は真に`q^n`個の相異なる元、節目(3)の核心部分が完成)**: 前回
「循環で見送った」と記録した`D_n`全体の分離性を、実は循環ではなく
`D_n=D_{n-1}・ψ_n`・`D_0=X`という漸化式を「分離性」と「将来の段の
`ψ_j`(`j>n`)との互いの素性」を**同時に保つ**帰納法で押し進めるだけで
完全に組み立てられると気づいた——前回見立てた「先に積表示を確立する
必要がある」という判断は誤りで、より直接的な同時帰納法で十分だった。
`Found/PGC/LubinTateDistinguishedSeparable.lean`(新規)で:

- `X_isCoprime_iteratedLubinTatePsi`: `X`と`ψ_j`(`j≥1`)は互いに素
  (`ψ_j`の定数項が非零`0<‖ψ_j.coeff0‖=‖π‖`から)
- ★★★★★★★★★★★★★`separable_and_coprime_iteratedLubinTateDistinguished_
  map`/`separable_iteratedLubinTateDistinguished_map`: `D_n`は分離的
- `nodup_roots_iteratedLubinTateDistinguished_map`: `K.closure`での
  `D_n`の根は互いに相異なる
- `card_roots_iteratedLubinTateDistinguished_map`: 根の個数はちょうど
  `q^n`

★これで実在のp進局所体`K:PAdicLocalField p`について、古典的
Lubin-Tate理論の`|Λ_n|=q^n`(真の集合として)そのものが完全に確立
された——節目(3)(捩れ点構成)の核心部分がsorry無しで完成した。

**残る作業**: 節目(3)の残り(`Λ_n`への`𝒪_K`加群構造の付与——`[a]_f`を
実際に捩れ点へ作用させる、`PowerSeries.aeval`+`CompleteSpace`経由で
見通しが立っている)、節目(4)`K(Λ_n)`の正規性(`ψ_n`の分解体と一致
すること)・Galois群の具体的計算`Gal(L_n/K)≅(𝒪_K/π^n)^×`(Lubin-Tate
の主定理)、節目(5)`L_π:=∪L_n`・`Gal(L_π/K)≅𝒪_K^×`、節目(6)不分岐
部分と合わせた相互律写像`K^×→Gal(K^ab/K)`の構成——pGC の各項目
(Prop 1.2・Cor 1.3・Prop 2.1・Prop 2.2・Theorem 4.2)を閉じるには、
なお相互律写像そのものの構成・性質証明という大きな仕事が残っている。

詳細な発見の経緯は `ResearchPaper/blocked-leaves.json` の
`[pGC] Proposition 1.2 / ... —— 局所類体論の相互律` エントリの
`progress2026_09_04a`〜`j` に記録。[[padic-log-additivity-blocked]]・
[[pgc-ramification-naturality-gap]] は同じ blocked エントリの別方向の前進。

**続報(2026-09-04・続き、`Λ_n`を`Finset`として梱包)**: 上記の
`card_roots_iteratedLubinTateDistinguished_map`(個数)・
`nodup_roots_iteratedLubinTateDistinguished_map`(重複無し)を、実際に
使いやすい`Finset`のレベルへ梱包した(同じ
`LubinTateDistinguishedSeparable.lean`に追加、commit `1630ebdc`):

- `iteratedLubinTateTorsionPoints`: `D_n`(`K.closure`へ写したもの)の
  根の`Finset`(`Multiset.toFinset`、`Classical`で`DecidableEq`供給)
- `card_iteratedLubinTateTorsionPoints`: `|Λ_n|=q^n`
  (`Multiset.toFinset_card_of_nodup`で上記2定理から)
- `zero_mem_iteratedLubinTateTorsionPoints`: `0∈Λ_n`
  (`D_n`の定数項が`0`であることから`Polynomial.eval₂_at_zero`経由)
- `iteratedLubinTateTorsionPoints_subset_of_le`: `n≤m`ならば`Λ_n⊆Λ_m`
  (`D_n∣D_m`から`Polynomial.roots.le_of_dvd`で根の`Multiset`の包含が
  出て、`Multiset.mem_of_le`で個々の元の所属に落とす)

★命名の教訓: 当初`torsionPoints`という短い名前で書いたところ、
`check.mjs`のG1検査(出典`.src`の検証)がバケット横断で宣言名を
フラットに突き合わせるため、既存の無関係な
`ABC3.Interface.GaloisRep.torsionPoints`(楕円曲線のn-捩れ点、`.src`
付き)と名前が衝突し、「`torsionPoints.src`の中身を読めなかった」と
いう偽陽性のNGが出た——`iteratedLubinTateTorsionPoints`へ改名して
解消。この文脈で新しい短い名前(`torsionPoints`・`kernel`・`image`等
プロジェクト内で使われがちな一般名)を付ける前に、`grep`で既存の
`.src`付き宣言との衝突を確認するとよい。

**続報(2026-09-04・続き、`Λ_n`の元は位相的冪零)**: 同じファイルへ
`spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints`を追加
(commit `4b02a6b7`)——`Λ_n`の任意の元`x`について
`spectralNorm K.carrier K.closure x < 1`。`D_n=D_{n-1}・ψ_n`・`D_0=X`
の漸化式に沿った帰納法(基底: `D_0=X`の根は`0`のみ、
`spectralNorm 0=0`。帰納段: `D_{n+1}`の根は`D_n`の根(IH)か
`ψ_{n+1}`の根で、後者は既出の`spectralNorm_root_iteratedLubinTatePsi`
により`‖π‖^(1/(q^{n+1}-q^n))`に等しく`0<‖π‖<1`かつ指数が正なので
`Real.rpow_lt_one`で`1`未満)。これは`Λ_n`への`𝒪_K`加群構造の構成
(`PowerSeries.aeval`による`[a]_f`の評価が要求する
`HasEval`/位相的冪零性の橋渡し)へ向けた前段。

**続報(2026-09-04・続き、`Λ_n`の元は整)**: 同じファイルへ2つ追加
(commit `71a0d276`):

- `isIntegral_of_mem_iteratedLubinTateTorsionPoints`: `Λ_n`の元は
  `𝒪_K`上整(`D_n`がモニックなので`IsIntegral`の定義そのものを直接
  満たす)
- `isIntegral_carrier_of_mem_iteratedLubinTateTorsionPoints`: `Λ_n`の
  元は`K.carrier`(局所体本体)上でも整(`IsIntegral.tower_top`で
  中間の環へ持ち上げるだけ)

これは`K.carrier⟮x⟯`の有限次拡大性
(`IntermediateField.adjoin.finiteDimensional`)、ひいては将来の
`PowerSeries.aeval`による`[a]_f`の評価に必要な完備性(有限次拡大は
完備)へ向けた土台——`K.closure`自体は完備でない(`Q_p`の代数閉包が
完備でないのと同様)ので、`x`を含む有限次拡大`K.carrier⟮x⟯`の中で
評価する、というのが見通している経路。

**続報(2026-09-04・続き、`K.carrier⟮x⟯`は有限次)**: 同じファイルへ
`finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints`を
追加(commit `ade9db4d`)——`x∈Λ_n`ならば`K.carrier⟮x⟯`は`K.carrier`
上有限次。直前の`isIntegral_carrier_of_mem_iteratedLubinTateTorsionPoints`
から`IntermediateField.adjoin.finiteDimensional`で直接従う(1行の
帰結)。次の一歩は「有限次拡大は完備」という一般論(mathlibに
`FiniteDimensional.complete`的な定理があるはず)で`K.carrier⟮x⟯`へ
`CompleteSpace`を移し、そこで初めて`PowerSeries.aeval`による`[a]_f`
の評価が可能になる。

**続報(2026-09-04・続き、★★★★★★★★★★`K.closure`にもノルム延長・
`K.carrier⟮x⟯`は完備)**: 予告した「有限次拡大は完備」の一般論は
mathlibに`FiniteDimensional.complete`としてそのまま存在した
(`Topology/Algebra/Module/FiniteDimension.lean:511`)。これを使うには
`K.closure`自身にノルム体構造が要る——`LocalFieldNorm.lean`が
`ℚ_[p]→K.carrier`に施した手順(`spectralNorm.normedField`)を、今度は
`K.carrier→K.closure`に**そのまま繰り返せる**ことに気づいた。

★鍵となる発見: `spectralNorm.normedField (K)(L)[NontriviallyNormedField K]
[Algebra K L][Algebra.IsAlgebraic K L][IsUltrametricDist K][CompleteSpace K]`
は**有限次拡大を要求しない**——基点`K`が完備でありさえすれば任意の
代数拡大`L`に働く。`K.closure`は`K.carrier`上一般に無限次(代数閉包
なので)だが、それでも`NormedField K.closure`が直接作れた。

`LocalFieldNorm.lean`に追加(commit`652b4c61`):
- `closureNormedField`・`closureNormedAlgebra`・`closureIsUltrametric`

`LubinTateDistinguishedSeparable.lean`に追加(同commit):
- `completeSpace_adjoin_of_mem_iteratedLubinTateTorsionPoints`:
  `x∈Λ_n`ならば`K.carrier⟮x⟯`は**完備**。
  `finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints`・
  `K.carrier`自身の完備性(既存)・上記`closureNormedField`等を
  `FiniteDimensional.complete`に渡すだけ。

★これで`PowerSeries.aeval`による`[a]_f`の`x`での評価が可能になる
舞台が整った——`K.closure`自体は完備でないが(`ℚ_p`の代数閉包が
完備でないのと同様)、`Λ_n`の各点`x`を添加した有限次拡大
`K.carrier⟮x⟯`の中でなら評価できる。次の一歩は実際に`[a]_f`(または
より一般に`f`)を`x`で評価し、`𝒪_K`加群構造を`Λ_n`に与えること
(`PowerSeries.HasEval`=`IsTopologicallyNilpotent`の要件は
`spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints`で既に
満たされている——`K.carrier⟮x⟯`上の`spectralNorm`と`K.closure`上の
`spectralNorm`が両立することの確認は必要かもしれない)。

**続報(2026-09-04・続き、`hasEval`は確立・★重要な構造的発見)**:
`hasEval_of_mem_iteratedLubinTateTorsionPoints`(`x∈Λ_n`ならば
`PowerSeries.HasEval x`)を追加(commit`0b9f3b2c`)——
`spectralNorm_lt_one_...`とmathlibの
`tendsto_pow_atTop_nhds_zero_of_norm_lt_one`を組み合わせるだけ。

★★★重要な発見(REPLで`PowerSeries.aeval`を実際に組み立てようとして
判明、ファイル未反映): `PowerSeries.aeval`が要求するもう片方の条件
`IsLinearTopology S S`(評価先`S`自身の位相がイデアルで生成される
「線形位相」)は、**体**`K.closure`や`K.carrier⟮x⟯`それ自身では
成り立たない——体のイデアルは`{0}`と全体しかないので、非自明な位相
とは両立しない。mathlibでこの条件が出る自然な形は
`Ideal.isLinearTopology`(付値環の極大イデアルによるadic位相)であり、
これは**付値環**(体でなく)に対して成り立つ。実際、古典的な
Lubin-Tate理論でも`[a]_f`は「体」でなく「付値環」の間の写像として
評価される——`Λ_n`の元は`𝒪_K`上整(既出の`isIntegral_..._carrier_...`)
なので、そもそも付値環の中に住んでいるはずのものを、うっかり体の
レベルで評価しようとしていた。

★見通している次の一歩: `L:=K.carrier⟮x⟯`は`K.carrier`上有限次
(推移律で`ℚ_[p]`上も有限次)なので、**`L`自身を新たな
`PAdicLocalField p`として再構成**すれば、`Found/PGC/LocalFieldNorm.lean`・
`Found/PGC/ValuationRingDVR.lean`・`ValuationRingComplete.lean`の
機構(`𝒪[L]`・`Valued L`・その adic 完備性)が**そのまま流用**でき、
`IsLinearTopology 𝒪[L] 𝒪[L]`が`Ideal.isLinearTopology`から従うはず。
この再構成(`L`を`PAdicLocalField p`のインスタンスとして正しく
組み立て、`K.carrier`からの埋め込みと両立することを示す)が次の
節目——それができれば`PowerSeries.aeval`で`f`を`x`(`𝒪[L]`の元として)
に評価し、`Λ_n`への`𝒪_K`加群構造へ到達できる見通し。

**続報(2026-09-04・続き、★★★★★★★★★★見通していた再構成を実行)**:
新規ファイル`AdjoinPAdicLocalField.lean`に`adjoinPAdicLocalField`を
追加(commit`60ab7e35`)——`x:K.closure`を`K.carrier`に添加した単純
拡大`K.carrier⟮x⟯`自身を、新たな`PAdicLocalField p`として組み立てる
定義。予想通り`Field`・`Algebra ℚ_[p] _`は`AlgebraicClosure`/
`IntermediateField`の一般論から**自動的に**手に入り(`infer_instance`
で確認)、`FiniteDimensional ℚ_[p] _`だけ`FiniteDimensional.trans`で
明示的に組み立てれば十分だった——想定より簡潔に済んだ。
`LubinTateDistinguishedSeparable.lean`には
`torsionPointAdjoinPAdicLocalField`(torsion-point版、`x∈Λ_n`から
必要な`FiniteDimensional`を供給するだけ)を追加。

★★★未解決のまま次に確かめるべき課題として記録: この`L=
adjoinPAdicLocalField K x`を実際に使うと、`LocalFieldNorm.lean`の
`normedField` scoped instance が `spectralNorm.normedField ℚ_[p] L.carrier`
で**新しい**ノルム構造を`L.carrier`(=`K.carrier⟮x⟯`)に与える。一方
`K.carrier⟮x⟯`は`K.closure`の部分体でもあるので、mathlibの一般論に
より`K.closure`に載せた`NormedField`(`closureNormedField`)から**別の**
自動的なノルム構造も持つ。両者は(スペクトルノルムの延長の一意性
から)数学的には一致するはずだが、definitionally一致するとは限らず、
両方が同時にスコープに入るとinstance diamondになりうる——次の一歩は
①両者が一致することを示す橋渡し補題を用意するか、②片方だけを常に
使うよう規律を決めるかのいずれか。これが解決すれば
`Ideal.isLinearTopology`経由で`IsLinearTopology 𝒪[L] 𝒪[L]`が手に入り、
`PowerSeries.aeval`で実際に`[a]_f`を評価できる見通し。

**続報(2026-09-04・続き、★★★★★★★★★★訂正: 上記の懸念は杞憂だった)**:
実際にREPLで検証したところ、上で懸念した①(体レベルでの
instance diamond)は**存在しなかった**——
`IntermediateField.adjoin K.carrier {x}`が`K.closure`の部分体として
mathlibの一般論から自動的に持つ`NormedField`構造は、そのままの定義で

```
‖(⟨x,_⟩ : K.carrier⟮x⟯)‖ = spectralNorm K.carrier K.closure x  -- `rfl` で成立
```

を満たす(部分体の`NormedField`は単に周囲のノルムの制限として定義
されているため)。`hasEval_mem_adjoin_of_mem_iteratedLubinTateTorsionPoints`
(commit`0f687b3d`)として、`Λ_n`の元を`K.carrier⟮x⟯`の元として見た
ときも`PowerSeries.HasEval`(位相的冪零性)が成り立つことを記録した。

一方、`adjoinPAdicLocalField`(`ℚ_[p]`から`K.carrier⟮x⟯`を再構成する
経路)の方は、実際に`spectralNorm ℚ_[p] K.closure`と
`spectralNorm K.carrier K.closure`の一致(
`NormedAlgebra.norm_eq_spectralNorm`+`NormedAlgebra.restrictScalars`
で試みた)を証明しようとしたところ、`Algebra.IsAlgebraic ℚ_[p]
K.closure`の型クラス探索が**タイムアウトする**という、単独では速く
解決するのに文脈内では詰まるという厄介な現象に遭遇した(原因は
未特定——恐らく`ℚ_[p]→K.carrier→K.closure`の`Algebra`構造が2通りの
経路(既存の自動導出 vs `NormedAlgebra.restrictScalars`が内部で
作るもの)で非defeqになり、探索が迷走している)。しかし体レベルの
上記`rfl`一致で目的(`HasEval`の`K.carrier⟮x⟯`版)は既に達成できた
ので、この`adjoinPAdicLocalField`経由のℚ_p再構成ルートは**今のところ
不要**——依然として残る`IsLinearTopology`(付値環が要る)の解決には
有用かもしれないが、優先度を下げて良い。

★教訓: 「instance diamond になりうる」という懸念は、実際に`rfl`や
`decide`で試してみるまでは確定した事実として記録・行動指針にしない
方が良い——今回は杞憂だったが、逆に`adjoinPAdicLocalField`側では
本物の(未解決の)型クラス探索の詰まりに遭遇した。両方あり得るので、
「試してみて確認する」を省略しない。

**続報(2026-09-04・続き、★★★★★★★★★★★★★大きな節目——`PowerSeries.
aeval`が`𝒪[K.carrier]`上で実際に組み立てられることを確認)**:
残っていた`IsLinearTopology`の壁を突破した。`ValuationRingComplete.
lean`の既出`isAdic_maximalIdeal_valuationRing`(距離位相=
`maximalIdeal`-進位相)と mathlib の`Ideal.isLinearTopology`(イデアル
の adic 位相は線形位相)を組み合わせるだけで
`isLinearTopology_valuationRing : IsLinearTopology 𝒪[K.carrier]
𝒪[K.carrier]`が3行で出た(commit`e290ad91`)。さらに
`isAdicComplete_valuationRing`の証明に埋もれていた素の`CompleteSpace
𝒪[K.carrier]`を`completeSpace_valuationRing`として独立の定理に抽出
(commit`2835e7c2`)。

これで実在の任意の p進局所体`K`について

```
PowerSeries.aeval (R := 𝒪[K.carrier]) hHasEval
```

が実際に型検査を通ることをREPLで確認した——mathlibの
`PowerSeries.aeval`が要求する全条件(`UniformSpace`・
`IsUniformAddGroup`・`IsTopologicalSemiring`・`T2Space`・
`CompleteSpace`・`IsTopologicalRing`・`IsLinearTopology S S`・
`ContinuousSMul`)が`𝒪[K.carrier]`について**すべて**揃うことの直接的な
証拠——このセッションで数ヶ月?かけて積み上げてきた「体でなく付値環
で評価する」という見通しが、ついに実際に動くコードとして確認できた
瞬間。

★次の一歩: `Λ_n`の元`x`(`K.closure`の元、`K.carrier⟮x⟯`上に住む)を
実際に評価するには、`R:=𝒪[K.carrier]`・`S:=𝒪[K.carrier⟮x⟯]`として
同じ機構を`K.carrier⟮x⟯`(`adjoinPAdicLocalField`経由、または直接
構成)に適用する必要がある——`valuationRing_isDVR`・
`completeSpace_valuationRing`・`isLinearTopology_valuationRing`は
すべて`K : PAdicLocalField p`について証明した一般定理なので、
`adjoinPAdicLocalField K x`(すでに`PAdicLocalField p`のインスタンス
として構成済み!)にそのまま適用できるはず。残る課題は
`Algebra 𝒪[K.carrier] 𝒪[(adjoinPAdicLocalField K x).carrier]`
(`𝒪[K.carrier]`が`𝒪[L]`へ埋め込まれること)の確立——これができれば
`f`(または`[a]_f`)を実際に`x`で評価できる。

**続報(2026-09-04・続き、正直な記録: 今回は新規コミット無し——
2つの経路それぞれで型クラス探索/単一化の詰まりに遭遇)**:

上記の課題(`𝒪[K.carrier]`を`𝒪[L]`へ埋め込む)に取り組んだが、
**どちらの経路も実際に動かすところまで到達できなかった**——
正直に記録する(ファイルへは何もコミットしていない):

1. **`adjoinPAdicLocalField`経由**(`ℚ_[p]`再構成のノルム使用):
   `spectralNorm.eq_of_tower (K:=ℚ_[p]) (L:=(adjoinPAdicLocalField
   K x).carrier) y`を試みたところ、`Algebra K.carrier
   (adjoinPAdicLocalField K x).carrier`が自動的に見つからず(`def`が
   `@[reducible]`でないため`.carrier`が展開されない)、無理に進めると
   `whnf`でタイムアウトした。

2. **`K.closure`継承ノルム経由**(`adjoinPAdicLocalField`を経由しない、
   より直接的な経路):`ProperSpace (K.carrier⟮x⟯)`
   (`FiniteDimensional.proper K.carrier _`から)・
   `IsUltrametricDist (K.carrier⟮x⟯)`(自動継承)は**あっさり確認
   できた**(★これは確定した有用な事実として記録)。しかし、そこから
   `Valued (K.carrier⟮x⟯) NNReal := NormedField.toValued`を導入し
   (`letI`でも `scoped instance`でも同様)、`Valued.integer.mem_iff`
   と`Metric.closedBall`を結びつけて`isCompact_closedBall`を適用
   しようとしたところ、**120秒経っても終わらない**深刻な単一化の
   詰まりに遭遇した(`maxHeartbeats`を200万まで上げても解決せず、
   バックグラウンドタスクとして強制終了した)。原因は未特定——
   `IntermediateField extends Subfield extends Subring extends
   Submonoid ...`という何層にも重なった部分構造の射影を通した
   `PseudoMetricSpace`/`UniformSpace`の一致判定が高コストになって
   いる可能性が高い。

★教訓・次回への申し送り: この特定の「`IntermediateField.adjoin`型に
対して`Valued.integer`+`isCompact_closedBall`を直接適用する」という
経路は、**単純に試すと危険**(セッションを長時間ブロックしうる)。
次に試すべき代替案:
(a) `x`を先に「素の」独立した型(`AdjoinRoot`や新しい`structure`で
    包んだもの)へ移してから`Valued`/`Metric`の議論をする、
(b) `IntermediateField`のまま作業するのを諦め、既存の
    `PAdicLocalField`の抽象論(`valuationRing_isDVR`等)を先に
    **完全に汎用化**(`K:PAdicLocalField p`でなく`[NormedField L]
    [ProperSpace L]`などの直接の型クラス仮定を取る形に書き直す)
    した上で、`adjoinPAdicLocalField`(ℚ_[p]再構成版、defeq問題は
    `@[reducible]`を付けて再挑戦)に適用する、
(c) そもそも`Valued.integer`の一般論を経由せず、
    `spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints`の
    ような**具体的な数値評価**だけで`𝒪[L]`への所属を直接示す
    (`x`の`ノルム`が分かっているので、それが`≤1`であることを見せる
    だけなら`Valued.integer.mem_iff`一発で済むはず——コンパクト性の
    議論自体は今回の直接の目的には実は不要かもしれない、と気づいた:
    `PowerSeries.aeval`に必要なのは`𝒪[L]`の`CompleteSpace`・
    `IsLinearTopology`であって、`IsDiscreteValuationRing`(そこから
    来る一意化元の存在)は`isAdic_maximalIdeal_valuationRing`の証明
    でのみ使われていた——`CompleteSpace`は`Valued.integer`の
    コンパクト性を経由しない別証明があるかもしれない)。

**続報(2026-09-04・続き、★★★★★★★★★原因を特定し回避——`Valued`を
避ければ全部速い)**: 上記(c)の勘が当たった。`Valued`を一切使わず
`Subring.mk`で`{y | ‖y‖ ≤ 1}`を素朴に構成したところ
(`AdjoinIntegers.lean`(新規)、commit`5cc2332e`)、
`CompactSpace`・`IsClosed`・`CompleteSpace`のすべてが**1秒未満**で
確認できた——`Valued`特有の位相と既存の`NormedField`由来の位相の
一致検査のコストが詰まりの原因だったと確定できた。

- `adjoinIntegers`: `K.carrier⟮x⟯`の「整数環」(ノルム`≤1`の部分環、
  `Valued`不使用)
- `algebraMap_mem_adjoinIntegers`: `𝒪[K.carrier]`の元は`algebraMap`
  を通して`adjoinIntegers K x`に入る(`spectralNorm_extends`——
  基点の元のスペクトルノルムはもとのノルムそのもの、という完全に
  汎用的な補題——から)。★これが`𝒪[K.carrier]`を`𝒪[K.carrier⟮x⟯]`へ
  埋め込む、何段階も前から見通していた課題の核心部分。
- `mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints`: `Λ_n`の
  元`x`自身も`adjoinIntegers K x`に入る
- `isClosed_adjoinIntegers`・`completeSpace_adjoinIntegers`・
  `compactSpace_adjoinIntegers`: 閉・完備・コンパクト、すべて確立

★残る課題(narrower になった): `IsDiscreteValuationRing`(ひいては
`Ideal.isLinearTopology`経由の`IsLinearTopology`)を`Valued`を経由
せずに得る経路。`Valued.integer.isDiscreteValuationRing_of_
compactSpace`は`Valued`型クラスを要求するので、この素朴な`Subring`
には直接使えない——mathlibに`Valued`を経由しない「コンパクトな
付値環はDVR」に相当する定理があるか、無ければ`IsPrincipalIdealRing`
を直接(コンパクト性+非アルキメデス距離から)証明する必要がある。
これが解決すれば`Algebra 𝒪[K.carrier] (adjoinIntegers K x)`
(`algebraMap_mem_adjoinIntegers`から`RingHom.toAlgebra`的に構成可能)
と組み合わせて、`PowerSeries.aeval`で実際に`f`を`x`で評価できる。

**続報(2026-09-04・続き、★★★★★★★★★★★★★このセッションの最大の
節目——冪級数を`Λ_n`の元で実際に評価する`lubinTateEvalAtTorsionPoint`
が完成した)**: 上記の課題を、当初計画していた
「`IsDiscreteValuationRing`を`Valued`抜きで得る」というルートとは
**違う、もっと直接的な2つの経路**で解決した(commit`88879620`):

1. **`ValuationRing`(付値環)**: `adjoinIntegers K x`の任意の2元
   `a,b`について、ノルムの大小で場合分けし、小さい方を大きい方で
   割った商(ノルム`≤1`なので`adjoinIntegers K x`自身の元)を
   `PreValuationRing`の定義`∃c,a*c=b∨b*c=a`の証人として**直接**
   与えるだけ——`Valued`・コンパクト性のどちらも不要だった。
   `ValuationRing.iff_local_bezout_domain`で`IsLocalRing`も従う。
2. **`IsLinearTopology`(★これが本当の鍵だった)**: 半径`ε>0`の
   閉球`{y|‖y‖≤ε}`が、非アルキメデス三角不等式(`add_mem'`)と
   ノルム`≤1`の元によるスカラー倍(`smul_mem'`)から**そのまま
   イデアルの公理を満たす**ことに気づいた(`adjoinIntegersBall`)。
   `Metric.nhds_basis_closedBall`(距離空間の一般論)と組み合わせて
   `IsLinearTopology.mk_of_hasBasis`を直接適用——`Ideal.
   isLinearTopology`(adic位相経由、`IsDiscreteValuationRing`が要る)
   も`Valued`も一切経由しない、**当初想定より遥かに直接的な経路**
   だった。「一意化元1本で割り切れる」という DVR 特有の構造は
   実は不要で、ノルムがバナッハ空間的に閉球=イデアルを作るという
   より初等的な事実だけで十分だった。

さらに:
- `adjoinIntegersAlgebraMap`・`adjoinIntegersAlgebra`: `𝒪[K.carrier]`
  から`adjoinIntegers K x`への環準同型・代数構造
  (`algebraMap_mem_adjoinIntegers`から)——**`𝒪[K.carrier]`が
  `𝒪[K.carrier⟮x⟯]`へ埋め込まれる**という何段階も前から見通していた
  課題がついに解決した。
- `continuous_algebraMap_adjoinIntegers`・`continuousSMul_
  adjoinIntegers`: 上記の埋め込みの連続性(`continuousSMul_of_
  algebraMap`から)。

これらすべてを組み合わせ、
**`lubinTateEvalAtTorsionPoint : PowerSeries (𝒪[K.carrier]) →ₐ[𝒪[K.
carrier]] adjoinIntegers K x`**——`PowerSeries.aeval`で任意の冪級数
(特に`iteratedLubinTate`系列、`[a]_f`を表すもの)を`Λ_n`の実際の元
`x`で評価する代数準同型——を完成させた。★これは`sorry`無しで、
実在の任意のp進局所体`K`について成り立つ、完全に一般的な結果。

**次の一歩**: `lubinTateEvalAtTorsionPoint`を`iteratedLubinTate`系列
(具体的に`[a]_f`を表すもの)に適用し、評価結果が再び`Λ_∞=∪Λ_n`に
戻ることを示して、`Λ_n`への`𝒪_K`加群構造(`a·x := [a]_f(x)`)を
実際に定義すること。これができれば節目(3)が完全に完成し、節目(4)
(`K(Λ_n)`の正規性・Galois群計算)へ進める。

**続報(2026-09-04・続き、`a·x`の実装が完成)**: `AdjoinIntegers.
lean`に`lubinTateActionAtTorsionPoint`を追加(commit`eab1c178`)——
`Λ_n`の元`x`について`a·x:=[a]_f(x)`を実際の元として計算する。
`LubinTateAction`(`[a]_f`を表す形式冪級数、`Found/PGC/
LubinTateActionEndomorphism.lean`で既に`sorry`無しで確立済み)を
`lubinTateEvalAtTorsionPoint`で`x`へ評価するだけ——新しい数学的
内容は要らない、純粋な「代入」。

**残る仕事(正直な記録、まだ手を付けていない)**: この作用`a·x`が
①`Λ_∞`(または少なくとも`Λ_n`自身)に**戻ること**(まだ示していない
——`a·x`は現状ただ「`adjoinIntegers K x`のどこかの元」というだけで、
それが実は既知の捩れ点の1つに一致することは未証明)、②加法性
`(a+b)·x = F_f(a·x,b·x)`(形式群の加法)、③乗法性`a·(b·x)=(ab)·x`
——これらはいずれも`LubinTateAction`自身の既に確立済みの関数
等式・準同型性を、`PowerSeries.subst`(形式的代入、`[a]_f`の構成に
使われた)ベースの等式から`PowerSeries.aeval`(位相的評価、
`lubinTateEvalAtTorsionPoint`が使う)を経由した等式へ**橋渡し**する
必要がある——`PowerSeries.substAlgHom_eq_aeval`(`DiscreteUniformity`
の元でsubstとaevalが一致するという趣旨の補題、mathlibに存在を確認
済み)がこの橋渡しの鍵になりそうだが、まだ実際には試していない。

**続報(2026-09-04・続き、単位律`1·x=x`を確認)**: 深いsubst/aeval
橋渡し理論を経由せずに得られる、もっと手前の具体的な事実を先に
押さえた(commit`bc5cbe76`):

- `aeval_X_eq_self`: `PowerSeries.aeval`は`X`を評価点そのものへ送る
  (`X`を`Polynomial.X`の像として書き直し、`PowerSeries.aeval_coe`・
  `Polynomial.aeval_X`と組み合わせるだけ。汎用的な事実)。
- `lubinTateActionAtTorsionPoint_one`: **`1·x=x`**。既に確立済みの
  `LubinTateAction_one_eq_X`(`[1]_f=X`)と`aeval_X_eq_self`を組み
  合わせるだけで従った。

★これは加群作用の公理のうち最初の1つ(単位律)が実際に成り立つ
ことを確認した最初の具体例——`subst`/`aeval`の深い橋渡し理論は
「加法性」「乗法性」(一般の`a,b`について)には依然として必要だが、
`a=1`という特殊ケースは`[1]_f=X`という既存の具体的な等式だけで
迂回できた。次に狙うべき比較的手前の目標: `a=0`の場合(`0·x=0`、
`LubinTateAction`に`0`での挙動を示す既存補題があれば同様に迂回
できるかもしれない)。

**続報(2026-09-04・続き、`a=0`の場合も解決——`0·x=0`)**: 予告通り
既存の補題だけで迂回できた(commit`8bf3e26b`)。
`LubinTateActionPiPow.lean`に**`LubinTateAction_zero_eq_zero`**
(`[0]_f=0`)を追加——`LubinTateAction_one_eq_X`と全く同じ
`powerSeries_uniqueness`の型で、鍵は`coeff_subst_zero_eq_zero_1var`
(`Found/PGC/LubinTateEndoBaseCase.lean`に**既に存在していた**、
`subst 0 f = 0`を`f`の定数項が0であることから示す1変数版補題)。
`AdjoinIntegers.lean`に**`lubinTateActionAtTorsionPoint_zero`**
(`0·x=0`)も追加——上記と`PowerSeries.aeval`が`AlgHom`(`0↦0`)で
あることを組み合わせるだけ。

★これで`lubinTateActionAtTorsionPoint`が満たすべき加群作用の公理の
うち、**単位律(`1·x=x`)・零元の吸収律(`0·x=0`)の2つ**が実際に
成り立つことを確認できた。教訓: 一般の`a,b`についての加法性・乗法性
には依然として`subst`/`aeval`の深い橋渡しが必要という見通しは
変わらないが、**特殊値(`a=0,1`)は既存の具体的な等式だけで完全に
迂回できる**——このプロジェクトの膨大な既存資産(`LubinTateAction`
関連ファイル群)には、探せばこの種の「特殊値での挙動」を示す補題が
既に眠っていることが多い、と再確認した。

**続報(2026-09-04・続き、★★★★★★★★★★定量的な事実`‖a·x‖≤‖x‖`を
subst/aeval橋渡し無しで確立)**: 一般の`a,b`についての加法性・乗法性
には依然として`subst`/`aeval`の深い橋渡しが必要だが、**「作用が
ノルムを増やさない」という定量的な事実は全く別の、遥かに単純な
経路で得られる**と気づいた(commit`b0cfe199`)。鍵:

1. `[a]_f`の定数項が`0`(`constantCoeff_LubinTateAction`、既出)
2. `PowerSeries.X_dvd_iff`(定数項`0`の冪級数は`X`で割り切れる、
   mathlib既存)で`[a]_f = X * h`と分解
3. `PowerSeries.aeval`が(積を保つ)**環準同型**であることから
   `a·x = aeval x(X*h) = x * (aeval x h)`——ここは単なる`map_mul`、
   `subst`は一切登場しない
4. `aeval x h`はそもそも`aeval`の**値域**である`adjoinIntegers K x`
   の元なので自動的にノルム`≤1`

これで`‖a·x‖=‖x‖*‖aeval x h‖≤‖x‖*1=‖x‖`(`norm_
lubinTateActionAtTorsionPoint_le`)、さらに`x`自身が位相的冪零
(`‖x‖<1`)と合わせて`‖a·x‖<1`(`norm_
lubinTateActionAtTorsionPoint_lt_one`)も得た。

★これで`a·x`自身も位相的冪零(`PowerSeries.HasEval`の条件を満たす
見込み)であることが分かり、**作用の反復適用**(`a·(b·x)`のような
合成)への土台が整った——次の一歩は実際に`a·(b·x)`を意味づけて
`a·(b·x)=(ab)·x`(乗法性)を、この反復可能性を使って狙うこと。
加法性(`F_f`経由)は依然として`subst`/`aeval`橋渡しが必要な見通し。

**続報(2026-09-04・続き、正直な記録: 乗法性への「近道」を探したが
見つからず)**: `a·(b·x)=(ab)·x`を`subst`/`aeval`の一般論を経由せず
得る近道として、`LubinTateAction_comp`(`[ab]_f=subst([b]_f)([a]_f)`、
`Found/PGC/LubinTateActionMul.lean`に既存)と`PowerSeries.aeval_
unique`(連続な`AlgHom`は`aeval`とその評価点で一意に特徴づけられる、
mathlib既存)を組み合わせる経路を検討した:
`ε:=(aeval x)∘(substAlgHom hg)`(`g:=[b]_f`、`hg:HasSubst g`は
`constantCoeff g=0`から自動)が連続な`AlgHom`であれば、
`aeval_unique`で`ε=aeval(ε X)`となり、`p:=[a]_f`で評価すれば
`aeval x([ab]_f)=aeval x(subst g [a]_f)=ε([a]_f)=aeval(ε X)([a]_f)
=aeval(aeval x g)([a]_f)=aeval(b·x)([a]_f)=a·(b·x)`(狙い通り)が
出るはずだった。

★しかし`ε`の**連続性**を示す段階で行き詰まった——mathlibで`subst`の
連続性を保証する唯一の補題`MvPowerSeries.continuous_subst`は
`[DiscreteUniformity R][DiscreteUniformity S]`(両側の環を**離散**
位相として扱う)を要求する。これは`aeval`/`HasEval`の枠組みが前提
とする「値付き体の位相」(非離散、`IsLinearTopology`)とは**別物**
——`subst`は元来「離散(=formal/代数的)」な設定でしか連続性が
保証されていない。したがって`substAlgHom hg`が`aeval`の使う
標準位相(次数ごとの収束による位相)について連続であることを示す
既存の道具はmathlibに見当たらなかった。

★これは「探せば近道が見つかる」という前回までの楽観的なパターンが
**今回は当てはまらなかった**という正直な記録——`subst`(形式的・
代数的)と`aeval`(位相的・収束による評価)を繋ぐ一般的な「連鎖律」
は、この特定のケース(値付き体への評価)について、mathlibにはまだ
無い**本物の理論的ギャップ**であり、埋めるには(a)`substAlgHom`が
次数ごとの収束位相について連続であることを独自に証明する(`subst`
の係数ごとの有限和公式から直接、`DiscreteUniformity`を経由せず)、
または(b)加法性・乗法性の証明を`subst`/`aeval`の一般論を経由しない
別の(まだ見えていない)経路で組み立て直す、のいずれかが必要——
どちらも本セッションの残り時間で片付く規模の作業ではなく、次回以降
への持ち越しとして正直に記録する。

**続報(2026-09-04・続き、★★★★★★★★★★★★★★★★理論的な突破——
上記のギャップを実際に埋めた)**: 「次回以降への持ち越し」と記録した
直後、実際に(a)の方向(独自の連続性議論)を、当初想定していた
「`substAlgHom`自体の連続性」ではなく**トランケーション経由の
別経路**で構成できた(commit`7f09476c`)。鍵となった発見の連鎖:

1. `PowerSeries.coeff_subst'`の`finsum`公式`∑_d coeff d p •
   coeff e(g^d)`は、`g`の定数項が`0`なら`g=X*h`と分解できることから
   `coeff e(g^d)=0`(`d>e`のとき)——つまり**`d≤e`の有限個の項だけ**
   で決まる(`coeff_pow_eq_zero_of_lt`)。
2. これで`subst g p`の次数`e`の係数は`p`を次数`e+1`でトランケート
   しても変わらない(`coeff_subst_trunc_eq`)——`subst`自身の連続性
   ではなく、その**係数公式の有限性**を直接使っただけ。
3. mathlibの`PowerSeries.WithPiTopology.tendsto_iff_coeff_tendsto`
   (各次数の収束の言い換え)と組み合わせ、`subst g(trunc N p)`が
   `N→∞`で`subst g p`へ収束することを示した(`tendsto_subst_trunc`)。
4. 既存の`PowerSeries.continuous_aeval`をこの収束に適用して極限を
   `aeval`の中へ運び、`Polynomial.aeval_algHom_apply`(代数準同型は
   多項式評価と可換、mathlib既存)で**有限段(多項式)での等式**を
   確立し、もう一度`PowerSeries.WithPiTopology.tendsto_trunc_atTop`
   で極限に戻すことで結論した:

```
★★★★★★★★★★★★★aeval_subst_eq_aeval_aeval:
  aeval x(subst g p) = aeval(aeval x g)p
```

★これで`subst`(形式的代入、`LubinTateAction`の構成に使われる)と
`aeval`(位相的評価、`lubinTateEvalAtTorsionPoint`が使う)を繋ぐ
**本物の連鎖律**が手に入った——`subst`自身の連続性は一切経由せず、
その係数ごとの有限性だけを使う、当初想定より遥かに直接的な経路
だった。教訓: 「ある操作Xの連続性」という入口で詰まったら、
「Xの出力の各成分は有限個の入力成分だけで決まる」という**離散的な
安定性**の議論で位相的な収束を回避できないか探す価値がある——
このセッション序盤で得た「入口を1つに決め打ちしない」という教訓
(`ValuationRingDVR.lean`の docstring)が、今回も別の形で当てはまった。

**次の一歩**: `aeval_subst_eq_aeval_aeval`と既存の
`LubinTateAction_comp`([ab]_f=subst([b]_f)([a]_f))を組み合わせて
乗法性`a·(b·x)=(ab)·x`を、`LubinTateAction_add`と組み合わせて
加法性`(a+b)·x=F_f(a·x,b·x)`を、それぞれ実際に証明すること。
`a·(b·x)`の意味づけ(`b·x`を新たな評価点として、同じ`adjoinIntegers
K x`の中で`aeval`を再度組み立てる——`norm_lubinTateActionAtTorsionPoint_
lt_one`から`b·x`のHasEvalは既に得られている)には、コアーション
(`adjoinIntegers K x`↔周囲の体)の慎重な追跡が必要で、まだ完成
していない。

**続報(2026-09-04・続き、★★★★★★★★★★★★★★★★節目達成——乗法性
`a·(b·x)=(ab)·x`を確立)**: 上で予告した「コアーションの慎重な追跡」
を実際にやり遂げた(commit`cb53f8ab`)。鍵となった追加補題:

- `hasEval_iff_coe`: `adjoinIntegers K x`の元`z`の位相的冪零性は、
  周囲の体での位相的冪零性と同値——`adjoinIntegers K x`の位相は
  周囲の体からの**誘導位相**(部分環)なので、mathlibの
  `tendsto_subtype_rng`(部分型での収束は包含写像を通した収束と
  同値、という一般論)から**直ちに**従った——想定より簡単だった。
- `hasEval_lubinTateActionAtTorsionPoint`: `b·x`は位相的冪零
  (`norm_lubinTateActionAtTorsionPoint_lt_one`+上記から)。
- `lubinTateEvalAtPoint`: `adjoinIntegers K x`の**任意**の位相的
  冪零な元`z`を評価点として冪級数を評価する——`x`に特化していた
  `lubinTateEvalAtTorsionPoint`の一般化。`CompleteSpace`・
  `IsLinearTopology`・`ContinuousSMul`は`x`(座標系)のみに依存し
  `z`には依存しないので、**座標変換無しでそのまま流用できた**
  (これが「慎重な追跡」の実態——恐れていたほど複雑ではなかった)。

そして本体:

```
★★★★★★★★★★★★★★★★lubinTateAction_mul: a·(b·x) = (ab)·x
```

証明は`LubinTateAction_comp`(`[ab]_f=subst([b]_f)([a]_f)`、既存)の
両辺を`x`で評価し、前回構成した連鎖律`aeval_subst_eq_aeval_aeval`で
右辺を`aeval(aeval x[b]_f)([a]_f)=a·(b·x)`へ変形するだけ——短い
証明だが、これまで積み上げてきた全ての層(`adjoinIntegers`の構成・
`IsLinearTopology`・連鎖律)が実際に噛み合って初めて可能になった。

★これで前回の「理論的な突破」(連鎖律の構成)が単なる技術的達成に
終わらず、**Lubin-Tate加群構造の核心的な性質(乗法性)の確立**に
直結したことが実証された。次の一歩は加法性
`(a+b)·x=F_f(a·x,b·x)`を、`LubinTateAction_add`と同じ連鎖律を
組み合わせて確立すること——`formalGroupLaw`(2変数)を経由する分
だけ乗法性より複雑になる見込みだが、基本戦略(既存の`subst`ベースの
恒等式の両辺を`aeval`で評価し、連鎖律で変形する)は同じはず。

**続報(2026-09-04・続き、★★★★★★★★★★★★★★★★★★★★最大の節目——
加法性を確立、𝒪_K加群公理が全て揃った)**: 予告通り、1変数連鎖律
(`aeval_subst_eq_aeval_aeval`)を2変数(`MvPowerSeries`)へ拡張する
だけで加法性も確立できた(commit`106726d7`)。想定していた複雑さは
**実際には現れなかった**——鍵となった簡略化:

★`family i = X * h_i`(各成分が1変数の冪級数)という事実のおかげで、
`d.prod(family s)^(d s)`という2変数の積が、実は**1変数の積**
`X^|d| * (h_0^(d 0) * h_1^(d 1))`に潰れる。これで「多変数`order`
理論」を一切使わず、1変数版と**全く同じ`X`分解トリック**
(`coeff_prod_pow_eq_zero_of_lt`)で係数の有限性が示せた——2変数化は
見た目ほど本質的な難所ではなかった。

追加した部品(1変数版の完全な2変数アナログ):
- `coeff_prod_pow_eq_zero_of_lt`・`coeff_subst_family_trunc_eq`・
  `tendsto_subst_family_trunc`: 上記の簡略化を使い、1変数版と同じ
  「`finsum`の有限性→係数の安定性→収束」という論法をそのまま踏襲。
- ★★★★★★★★★★★★★★★★`aeval_subst_family_eq_aeval_aeval`:
  連鎖律の2変数版。`MvPolynomial.comp_aeval_apply`(1変数版
  `Polynomial.aeval_algHom_apply`のmathlib既存の多変数アナログ、
  探したらそのまま見つかった)・`MvPowerSeries.continuous_aeval`・
  `MvPowerSeries.WithPiTopology.tendsto_trunc'_atTop`を組み合わせる。
- `hasEval_actionFam2`: `family:=fun i=>if i=0 then a·x else b·x`が
  `MvPowerSeries.HasEval`を満たすこと——`MvPowerSeries.HasEval`は
  「各成分の位相的冪零性」+「`cofinite`フィルターでの`0`への収束」
  の2条件からなるが、**`Fin 2`は有限型なので`cofinite`フィルターは
  自動的に`⊥`になり**、`Tendsto _ ⊥ _`は常に真——第2条件は完全に
  無償で手に入った。
- `lubinTateEvalFormalGroupAt`: `formalGroupLaw`を`adjoinIntegers
  K x`の点の族で評価する、`lubinTateEvalAtPoint`の2変数版。

そして本体:

```
★★★★★★★★★★★★★★★★★★★★lubinTateAction_add: (a+b)·x = F_f(a·x,b·x)
```

★★★これで`lubinTateActionAtTorsionPoint`が満たすべき`𝒪_K`加群の
構造公理——**単位律・零元の吸収律・乗法性・加法性のすべて**——が
確立された。実在の任意のp進局所体`K`について、実在の捩れ点
`x∈Λ_n`に対する古典的なLubin-Tate加群構造が、`sorry`無しで完全に
構築されたことになる——このセッションを通じて積み上げてきた
Lubin-Tate理論の実装が、ここで1つの大きな到達点に達した。

**次の一歩**: この作用が実際に`Λ_∞(=∪Λ_n)`に戻ることを示す(現状、
`a·x`は「`adjoinIntegers K x`のどこかの元」というだけで、それが
既知の捩れ点の1つに一致することは未証明)。これができれば
`K(Λ_n)`の正規性・Galois群計算`Gal(L_n/K)≅(𝒪_K/π^n)^×`(古典的な
Lubin-Tate理論の主定理)へ進める見通し。

**続報(2026-09-04・続き、★★★★★★★★★★★★★★★★★★★★節目達成——
予告した「作用が`Λ_n`に戻る」ことを実際に確立した)**: commit
`4744a07a`。鍵となった証明の筋(古典的な議論をそのまま形式化した):

1. **`x∈Λ_n⟹π^n·x=0`**(`pi_pow_action_eq_zero`): `D_n(x)=0`
   (`Λ_n`の定義そのもの、`adjoinIntegers K x`から`K.closure`への
   単射環準同型`g`を通して`Polynomial.aeval`のレベルへ引き戻す
   ——`aeval_iteratedLubinTateDistinguished_eq_zero`)と
   `[π^n]_f=D_n*U_n`(既存)を掛け合わせるだけ(`D_n(x)*U_n(x)=
   0*U_n(x)=0`、`U_n`が単位かどうかは不要——**この向きは簡単**)。
2. **`π^n·(a·x)=a·(π^n·x)=a·0=0`**(`pi_pow_action_action_eq_zero`):
   乗法性(`lubinTateAction_mul`)を2回・可換律を使うだけ。
3. **`[π^n]_f(z)=0⟹D_n(z)=0`**(`eq_zero_of_pi_pow_action_eq_zero`、
   **逆向きはここが本質的に難しい**): `[π^n]_f=D_n*U_n`の`U_n(z)`
   が**単位**であること(単位を環準同型`aeval`で送った像は自動的に
   単位)を使い、`D_n(z)*U_n(z)=0`から単位`U_n(z)`を約分
   (`IsUnit.mul_left_eq_zero`)して`D_n(z)=0`を得る。
4. 3.の結果(`D_n(a·x)=0`、`adjoinIntegers K x`の中での話)を
   `Polynomial.hom_eval₂`で`K.closure`のレベルへ押し出し
   (同じ単射環準同型`g`を使う、1.の逆向きの計算)、
   `iteratedLubinTateTorsionPoints`の定義そのものに一致させる
   ——**`lubinTateActionAtTorsionPoint_mem`**: `a·x∈Λ_n`。

★★★これで、実在の任意のp進局所体`K`について、古典的なLubin-Tate
理論の**`Λ_n`が真に`𝒪_K`-加群である**という中核事実が`sorry`無しで
完全に確立された——単位律・零元律・乗法性・加法性(前コミット)に
続き、作用が実際に`Λ_n`自身に留まることまで示せたことになる。
`Λ_n`の`𝒪_K`加群としての構造は、これで数学的に完全な形になった。

**次の一歩**: `K(Λ_n)`(`Λ_n`の元をすべて`K`に添加した体)の
`K`上の正規性、そしてGalois群の具体的計算
`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`(古典的なLubin-Tate理論の主定理)。
これは今回確立した「`Λ_n`は`𝒪_K`-加群」という事実と、既に確立済みの
`|Λ_n|=q^n`(`card_iteratedLubinTateTorsionPoints`)を組み合わせて
狙う——`𝒪_K/π^n`が`Λ_n`に(単位元`ζ`を1つ選べば)**単純推移的に**
作用することを示せば、Galois群の計算に繋がる見通し。

**続報(2026-09-04、`Λ_n=Λ_{n-1}∪(ψ_nの根)`、`Found/PGC/
LubinTateDistinguishedSeparable.lean`、commit`12d4b722`)**: 上の
「次の一歩」への最初の具体的な布石として、`D_n`の`K.closure`上の
分解`D_n=D_{n-1}・ψ_n`(既出`iteratedLubinTateDistinguished_eq_mul_psi`)
を`Λ_n`自身のFinsetレベルへ運んだ:

1. `iteratedLubinTatePsiTorsionPoints`: `ψ_n`(`K.closure`へ写したもの)
   の根のなす`Finset`——「原始的な」π^n-捩れ点全体。`Λ_n`
   (`iteratedLubinTateTorsionPoints`)と同じ`Multiset.toFinset`パターン。
2. `nodup_roots_iteratedLubinTatePsi_map`・`card_roots_iteratedLubinTatePsi_map`:
   `D_n`について確立済みの2定理(`nodup_roots_iteratedLubinTateDistinguished_map`・
   `card_roots_iteratedLubinTateDistinguished_map`)と**全く同じ手筋**
   (`algebraMap`の分解+`Polynomial.map_map`、`separable_iteratedLubinTatePsi_map_carrier`
   は`LubinTatePsiNorm.lean`に既出)を`ψ_n`へ複製するだけ——
   `|ψ_nの根|=q^n-q^{n-1}`。
3. ★**`iteratedLubinTateTorsionPoints_eq_union`**: `Λ_n=Λ_{n-1}∪
   (ψ_nの根)`。`D_n=D_{n-1}・ψ_n`を`K.closure`へ写して
   `Polynomial.roots_mul`(積の根はmultisetの加法)を適用し、
   `Multiset.toFinset`の言葉に戻すだけ——新しい数学的道具は不要、
   既出の事実の組み合わせのみ。
4. `card_iteratedLubinTatePsiTorsionPoints`: `|ψ_nの根|=q^n-q^{n-1}`
   (Finsetの濃度版、`nodup`+`card_roots`から)。

これで「原始的な」π^n-捩れ点(`Λ_n`にあって`Λ_{n-1}`に無いもの)が
ちょうど`ψ_n`の根であることが、Finsetの言葉で明示された——
`Gal(K(Λ_n)/K)`の作用が原始的捩れ点上どう振る舞うか(単純推移性)
を調べるための土台。

**続報(同日、非交和として確定、commit`11d5e24b`)**: 上で残した課題
(交わりが無いことの明示)をすぐに解消した:

1. ★**`iteratedLubinTateTorsionPoints_disjoint_iteratedLubinTatePsiTorsionPoints`**:
   `Λ_{n-1}`と`ψ_nの根`は交わらない。`D_n=D_{n-1}・ψ_n`が重複無しの根を
   持つ(`nodup_roots_iteratedLubinTateDistinguished_map`)ことから、
   両因子に共通根`x`があれば`Multiset.count x`が両側から`≥1`ずつ足し
   合わさって`(D_n).roots`での重複度`≥2`になり`Nodup`に矛盾——という
   `Multiset.count_add`+`omega`だけの短い議論。
2. **`iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints`**:
   合併(前コミット)+非交和を`Finset.union_sdiff_cancel_left`で
   組み合わせ、`Λ_n\Λ_{n-1}=ψ_nの根`を確定。

これで「原始的なπ^n-捩れ点=`ψ_n`の根」という特徴づけがFinsetレベルで
完全に閉じた——次は`Λ_n`(あるいは`ψ_nの根`)への`𝒪_K/π^n`(あるいは
`(𝒪_K/π^n)^×`)の作用の単純推移性、そしてそこから`Gal(K(Λ_n)/K)`の
構造への橋渡しに進む。

**続報(同日、非空性・非自明性、commit`9476ba69`)**: 「生成元を1つ選ぶ」
議論の前提となる基本事実を3つ追加した:

1. `iteratedLubinTatePsiTorsionPoints_nonempty`: `|ψ_nの根|=
   q^{n-1}(q-1)>0`(`q>1`、既出`Fintype.one_lt_card`)から、
   「原始的な」π^n-捩れ点が実際に存在することを確認。
2. `coeff_zero_iteratedLubinTatePsi_ne_zero`: `ψ_n`の定数項
   (`𝒪[K.carrier]`の元)が非零——`‖ψ_n.coeff0‖=‖π‖≠0`
   (既出`norm_iteratedLubinTatePsi_coeff_zero`)から。
3. ★`zero_not_mem_iteratedLubinTatePsiTorsionPoints`: `0`は`ψ_n`の
   根ではない——2.と`𝒪[K.carrier]→K.closure`の単射性(`K.carrier`
   への包含+体拡大の単射性、`Subtype.coe_injective`と
   `(algebraMap K.carrier K.closure).injective`の合成、
   `FaithfulSMul`系の重いインスタンス探索は使わず直接構成)から。
   `Λ_0={0}`なので、「原始的な」捩れ点が`Λ_0`の元と重ならない
   非自明な元であることの確認になる。

次は`(𝒪_K/π^n)`(あるいは`(𝒪_K/π^n)^×`)の`Λ_n`(または`ψ_nの根`)
への作用の単純推移性——ただし**注意**: Lubin-Tate加群の加法は
`AdjoinIntegers.lean::lubinTateAction_add`にある通り形式群則`F_f`に
よるものであり、環の通常の加法とは一般に一致しない(`F_f(X,Y)≡X+Y
mod deg2`だが高次では異なる)。したがって「`a≡b mod π^n`ならば
`a·x=b·x`」(well-definedness)を示すには通常の減法ではなく`F_f`
加法群としての逆元・結合律などの形式群公理を経由する必要があり、
これは次の大きめの一歩になる見通し。

**続報(同日、★核の確定、`AdjoinIntegers.lean`、commit`4e4499f4`)**:
上の懸念(`F_f`加法の逆元・結合律が要る)を実際には**回避できる**
ことが分かった——「`a≡b mod π^n`なら`a·x=b·x`」という**加法群
としての**well-definednessではなく、写像`a↦a·x`(`𝒪_K→Λ_n`)の
**核**を直接特徴づける、というより直接的な道筋で古典的な同型
`𝒪_K/π^n≅Λ_n`に迫れる:

1. `lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_
   iteratedLubinTatePsiTorsionPoints`: 「原始的な」π^n-捩れ点
   (`ψ_n`の根)は`π^{n-1}`で消えない——`x∉Λ_{n-1}`(前々回の非交和)
   の対偶を、`aeval_iteratedLubinTateDistinguished_eq_zero`と同じ
   `Polynomial.hom_eval₂`橋渡しで使う。
2. `irreducible_of_maximalIdeal_eq_span`: `hπmax`から`π`が既約元
   であること(`IsDiscreteValuationRing.exists_irreducible`+同伴)。
3. ★★★**`lubinTateActionAtTorsionPoint_eq_zero_imp_dvd_of_mem_
   iteratedLubinTatePsiTorsionPoints`**(核心、⟹方向): `a·x=0⟹
   π^n∣a`。付値環の一意分解`a=u*π^k`(`u`単位、mathlib
   `IsDiscreteValuationRing.eq_unit_mul_pow_irreducible`)を取り、
   `k<n`と仮定して矛盾を導く——`lubinTateAction_mul`(乗法性)**だけ**
   で`u⁻¹`を作用させて単位を約分し`π^k·x=0`を復元、さらに
   `π^{n-1-k}`を掛けて`π^{n-1}·x=0`を導き1.に矛盾する。**`F_f`加法の
   逆元・結合律は一切使わなかった**——乗法性のみで核の非自明な
   半分が閉じた、という見立て違いの良い意味での訂正。
4. ★**`lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_mem_
   iteratedLubinTatePsiTorsionPoints`**(iff): ⟸は既出
   `pi_pow_action_action_eq_zero`そのもの、⟹は3.。

**続報(同日、単位元則の評価版+主単数がxを固定、commit`736f086c`)**:
「単射性一般には`F_f`の引き算が要る」という壁は残るが、**特別な
場合**(主単数`1+π^n𝒪_K`が`x`を固定すること)は`F_f`の**単位元則
だけ**で示せることを発見し、実際に確立した:

1. `Found/PGC/LubinTateIdentityLaw.lean`に**単位元則を評価レベルへ
   持ち上げる**新しい章を追加:
   - `coeff_single0_formalGroupLaw`: `F_f`の「`X_0`のみ」の係数
     (`X_1`次数`0`)は`X`自身の係数と一致——`psi`の定義(`coeff_
     restrictR_eq_of_e1_zero`が`X_1↦0`の代入で`X_1`次数`0`の係数を
     保つことを保証)と`formalGroupLaw_identity`(`psi=X`)を組み合わ
     せるだけ、`restrictR`/`psi`自体を経由した「形式的な」代入の
     世界を出ずに得られる係数レベルの帰結。
   - ★★★**`aeval_formalGroupLaw_eq_of_snd_eq_zero`**: `F_f(y₀,0)=y₀`
     (`y₁=0`のとき)。`y₁=0`だと`d₁≥1`の項は`y₁^{d₁}=0`で消え、
     `d₁=0`の項は`coeff_single0_formalGroupLaw`により`d₀=1`(係数`1`)
     以外すべて`0`——`Finset.sum_eq_single`で単項だけ残る。
     truncationの評価値が(`n₀≥1`である限り)**常に**`y₀`という
     「一定値」の列になることと、`aeval`の定義そのもの(truncationの
     極限)への収束を`tendsto_nhds_unique`で結びつける——このセッションで
     `aeval_subst_eq_aeval_aeval`を作った時と全く同じtruncation-limit
     手法の再利用。
2. `Found/PGC/AdjoinIntegers.lean`に★★★**`one_add_mul_pi_pow_action_eq_self`**:
   **`(1+c*π^n)·x = x`**(`x∈Λ_n`、任意の`n`・`c∈𝒪_K`)。
   `lubinTateAction_add`(`(1+cπ^n)·x=F_f(1·x,(cπ^n)·x)`)・
   `lubinTateActionAtTorsionPoint_one`(`1·x=x`)・`pi_pow_action_
   action_eq_zero`(`(cπ^n)·x=0`)を代入し、上の評価レベル単位元則
   `F_f(x,0)=x`で結論する。**`F_f`加法の逆元・結合律は一切不要**
   ——単位元則だけで閉じた。

**教訓(証明の作り方)**: 「`hfam_eq : family1=family2`を`have`で
別立てし`congr 1`で橋渡しする」という当初の方針は、`Fin 2`の
`if i=0`の`Decidable`インスタンスが箇所ごとに一致しない(`rw`/
`simp only [if_pos rfl]`では閉じない、`exact`の defeq 判定でしか
閉じない場合がある)という技術的な罠にはまり続けた。最終的に
`simp only [lubinTateActionAtTorsionPoint_one, hcpi0] at hadd`
——**既存の等式を直接書き換える**——という遥かに単純な経路に
気づいてから一気に片付いた。次に似た「`if`の入った関数の等式」に
当たったら、まずこちらを試す。

**続報(同日、乗法的不変性、commit`0bd61b30`)**: `one_add_mul_pi_pow_
action_eq_self`(`(1+cπ^n)·x=x`)の直接の系として、待望していた
「well-definedness」の**乗法版**をすぐに確立できた:

1. `lubinTateEvalAtPoint_congr`: `lubinTateEvalAtPoint`は評価点が
   (等式として)一致すれば同じ値を返す——既存の`lubinTateEvalAtPoint_
   eq_zero_of_eq_zero`の`w=0`限定を一般の`w`へ広げた、再利用可能な
   補題(依存型の`rw`が失敗する箇所をすべて`subst`で回避できる)。
2. ★★★**`mul_one_add_pi_pow_action_eq`**: **`(a*(1+cπ^n))·x = a·x`**
   (任意の`a,c∈𝒪_K`)。`lubinTateAction_mul`+`one_add_mul_pi_pow_
   action_eq_self`+`lubinTateEvalAtPoint_congr`を組み合わせるだけ。

これで「`a↦a·x`は`a`を`(1+π^n𝒪_K)`を法として見たとき不変」——
`(𝒪_K)^×`から`ψ_nの根`への作用が`(𝒪_K/π^n)^×`を経由してwell-defined
であることの**乗法版のwell-definedness**が確立された。残る壁は
「加法版」(`a≡b mod π^n`ならば`a·x=b·x`、こちらは依然`F_f`の
引き算が要る)だけであり、`(𝒪_K)^×`(乗法群)に限った話であれば
この乗法版だけで十分な場合が多い——次の一歩は`(𝒪_K)^×→ψ_nの根`
の写像が本当に`(𝒪_K/π^n)^×→ψ_nの根`の写像として well-defined に
factor through することの正式な組み立てと、単射性(こちらは依然
残された課題)。

**続報(同日、`principalUnits`部分群、commit`ed578fe1`)**: 上の
「正式な組み立て」の第一歩として、`1+π^n𝒪_K`の形の単数のなす
`(𝒪_K)^×`の**部分群**を明示的に構成した:

1. `principalUnits K π n : Subgroup (𝒪_K)^×`: `(v-1)`の式変形だけで
   単位元・積・逆元の閉性が出る、純粋に環論的な構成
   (`vw-1=v(w-1)+(v-1)`・`v⁻¹-1=-v⁻¹(v-1)`)——`F_f`・Lubin-Tate
   固有の議論は一切不要。
2. `mem_principalUnits_iff`: 定義の直接の言い換え(`v=1+cπ^n`の形)。
3. ★★★**`mul_principalUnits_action_eq`**: `v∈principalUnits K π n`
   ならば`(u*v)·x=u·x`——前回の`mul_one_add_pi_pow_action_eq`を
   部分群の言葉で言い換えただけ。

これで「`(𝒪_K)^×`の`Λ_n`への作用は`principalUnits K π n`を法として
well-defined」ことが**群論的に正しい形**(`Subgroup`+その法での不変性)
で定式化できた。次の一歩: `(𝒪_K)^×⧸principalUnits K π n`
(標準的には`(𝒪_K/π^n)^×`に同型のはずだが、その同型自体はまだ
構築していない)から`ψ_nの根`への**誘導された写像**
(`QuotientGroup.lift`)を実際に構成すること——`mul_principalUnits_
action_eq`はまさにこの`lift`が要求する「well-definedness」の
証明そのものになっている。

**続報(同日、`QuotientGroup.lift`の実現、commit`978198a6`)**: 予告
通り、この誘導写像を実際に構成した:

1. ★★★**`unitActionQuotientLift`**: `(𝒪_K)^×⧸principalUnits K π n
   → adjoinIntegers K x`。`Quotient.lift`に`mul_principalUnits_
   action_eq`をそのままwell-definedness証明として渡すだけ——
   `QuotientGroup.leftRel_apply`(`a≈b ↔ a⁻¹*b∈N`)で`Quotient`の
   同値関係`≈`を部分群の言葉へ変換する橋渡しが唯一の技術的な段
   (`rw`では`≈`のパターンにマッチせず、`.mp`を項レベルで適用する
   ことで解決——`change`によるdefeq判定でも通らなかった、また
   1つの`if`の`Decidable`インスタンス不一致に似た軽い罠)。
2. `unitActionQuotientLift_mk`: `QuotientGroup.mk u`での値が文字通り
   `u·x`であること(`rfl`)——定義の直接の確認。

これで「単数の作用が有限群`(𝒪_K)^×⧸principalUnits K π n`上の写像
として矛盾なく定義できる」ことがmathlibの`QuotientGroup`インフラの
上で明示的に構成された。次の一歩(いずれも本質的な数学的内容が
要る、未解決): (a) この写像の**像**が本当に`ψ_nの根`全体
(`iteratedLubinTatePsiTorsionPoints`)に一致すること(全射性)、
(b) この写像が**単射**であること(`F_f`の引き算/キャンセレーション
が要る、これまでの最大の壁)、(c) `(𝒪_K)^×⧸principalUnits K π n`が
実際に`(𝒪_K/π^n)^×`に同型であること(標準的だが未構築)。

**続報(同日、(c)を解決、commit`563ccebb`)**: 課題(c)——
`(𝒪_K)^×⧸principalUnits K π n`が古典的な`(𝒪_K/π^n)^×`に同型である
こと——を、`F_f`・`Λ_n`・`ψ_n`固有の議論を一切経由せず、**純粋に
環論的な事実**として解決した:

1. `principalUnits_eq_ker`: `principalUnits K π n`は、剰余写像
   `𝒪_K→𝒪_K/π^n𝒪_K`が誘導する単数群の写像`Units.map`の**核**に
   一致する——`v↦v-1`の言葉での定義を`Ideal.Quotient.eq_zero_iff_
   mem`で言い換えるだけ。
2. `nontrivial_quotient_span_pi_pow`: `𝒪_K/π^n𝒪_K`は非自明
   (`n≥1`のとき)——`π^n∈maximalIdeal`から`span{π^n}⊊⊤`。
3. ★★★★★★★★★★**`principalUnitsQuotientEquiv`**: **`(𝒪_K)^×⧸
   principalUnits K π n ≃* (𝒪_K/π^n𝒪_K)^×`**。鍵はmathlib
   `IsLocalRing.surjective_units_map_of_local_ringHom`——**局所環
   からの全射環準同型が誘導する単数群の写像は全射**、という
   ドンピシャの既存定理が見つかったこと(`𝒪_K`が局所環であることと
   `IsLocalHom.of_surjective`だけで適用できる)。あとは1.と組み合わせ
   `QuotientGroup.quotientKerEquivOfSurjective`(第一同型定理)で
   即座に結論した。

これで`unitActionQuotientLift`の**定義域**`(𝒪_K)^×⧸principalUnits
K π n`が、古典的な`(𝒪_K/π^n)^×`と正式に同一視できるようになった
——`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`という主定理の**左辺**(定義域)が
完全に確立された形。残る課題は(a)全射性・(b)単射性(いずれも
Lubin-Tate固有の深い議論——`F_f`のキャンセレーションまたは
successive approximationが要る)のみ。

**続報(同日、★★★濃度の一致を達成、新規ファイル`Found/PGC/
QuotientCardinality.lean`(333行)、commit`2c94e367`)**: 全射性・
単射性の議論を直接攻めるのではなく、**両者を一発で結ぶ「濃度の
一致」**を、`F_f`・`Λ_n`・`ψ_n`固有の議論を一切経由せず、**純粋に
環論的な事実として**確立した:

1. `gradedPieceMap`/`maxIdealCardMap`(`x↦π^n*x`・`x↦π*x`が誘導する
   次数付き商の写像、`Submodule.mapQ`)とその核(`⊥`)・像(次数付き
   ピース、`π^n𝒪_K/π^{n+1}𝒪_K`または`π𝒪_K/π^{n+1}𝒪_K`)を計算し、
   `Submodule.card_quotient_mul_card_quotient`(第三同型定理の濃度版、
   mathlib既存)で帰納法を回す。
2. ★★★★★★★★★★★★★★★★★★★★★★★★**`card_quotient_span_pi_pow`**:
   **`|𝒪_K/π^n𝒪_K| = q^n`**。
3. ★★★★★★★★★★★★★★★★★★★★★★★★**`card_units_quotient_span_pi_pow`**:
   **`|(𝒪_K/π^n𝒪_K)^×| = q^n-q^{n-1}`**——局所環`R=R^×⊔maximalIdeal R`
   (集合として、`Equiv.sumCompl IsUnit`)の分解と`maxIdealCardMap`
   (`maximalIdeal R`の濃度が`𝒪_K/π^{n-1}𝒪_K`に一致すること、
   `IsLocalRing.map_maximalIdeal_of_surjective`で最大イデアルの像を
   同定)を組み合わせる。
4. ★★★★★★★★★★★★★★★★★★★★★★★★**`card_principalUnitsQuotient`**:
   **`|(𝒪_K)^×⧸principalUnits K π n| = q^n-q^{n-1}`**——
   `card_iteratedLubinTatePsiTorsionPoints`(`|ψ_nの根|=q^n-q^{n-1}`、
   既出)と**ぴったり一致する濃度**。

これで`unitActionQuotientLift`(定義域`(𝒪_K)^×⧸principalUnits K π n`、
値域`ψ_nの根`を含む`adjoinIntegers K x`)が**単射・全射のどちらか
一方さえ示せば**、有限集合・同じ濃度の写像は「単射⟺全射⟺全単射」
という一般論だけでもう一方が自動的に従う——数え上げの土台が完成した。
残る課題は「単射または全射のどちらか一方」に絞られた(以前は両方が
独立した課題だった)。

**教訓(技術的な罠)**: `n-1+1=n`が変数`n`に対して定義的に成り立たない
(自然数減法は`omega`等の証明を要する、`rfl`では通らない)ことから、
`LinearMap.range (何かK π (n-1))`のような式の周辺で`isDefEq`が
**極端に遅くなる**(数十秒〜タイムアウト)罠に複数回はまった——
今セッションで繰り返し遭遇してきた「`Valued`インスタンスの重複」
と同じ「深くネストした型の再帰的unify」という根本原因のファミリー。
回避策: 該当の補題を`n`ではなく`m+1`(`m:=n-1`)の形で最初から
述べ、`n,hn:1≤n`版は最後に`Nat.exists_eq_add_of_le`で言い換える
だけにする——型が最初から構文的に一致するので`isDefEq`のコストが
かからない。次に似た「`n-1`が型に登場する」状況に当たったら、
まずこの`m+1`書き直しを試す。

## ★★★戦略的確認(サブエージェント調査、同日): Proposition 2.1 の
残タスクも「同じ相互律」だったと確定

`.needs`欄が「p進対数(解消済み)+Verlagerung(mathlib既存`MonoidHom.
transfer`)+最後の1段(`.implicitStep`)」としか書いておらず、この
「1段」がbookkeeping程度で済むなら、このLubin-Tate trackより先に
`/goal`の項目を1つ閉じられるかもしれない、と考えて別のサブエージェント
に調査させた。結論: **同じ相互律に依存していた**——`ResearchPaper/
blocked-leaves.json`の既存監査(`progress2026_09_04n`)が既にこの
結論に到達しており、原文([pGC] p.4)自身も「`U_K`を`Γ_K^ab`の
部分群とみなす」の一文が実質的に相互律そのものであることを裏付けた。
`padicLog_bijOn`(半径1/4の球上の全単射)は`.needs`の「対数」部分
(解析的な半分)を解消しただけで、`U_K↪Γ_K^ab`という群論的な
識別(相互律の核心)には一切触れていない——**近道は無かった**。

**結論**: この Lubin-Tate track こそが、8項目中5〜6項目
(Prop 1.2・Cor 1.3・Prop 2.1・Cor 3.1(部分)・Cor 3.3・Theorem 4.2)
に共通する、唯一の実質的なボトルネックであることが、2回の独立した
調査で再確認された。回り道は無い——この track を最後まで押し切る
以外に近道は無いと判断し、続行する。

## 次の技術的な大きな一歩(まだ未着手): 単射性/全射性を解決する
3つの選択肢

濃度の一致(`card_principalUnitsQuotient`)により、単射性か全射性の
**どちらか一方**さえ示せば数え上げでもう一方が自動的に従うところ
まで来た。今日1日で検討した3つの候補経路とその評価:

1. **`F_f`の形式逆元冪級数**(陰関数定理型構成)——`F_f(z,w1)=F_f(z,w2)
   ⟹w1=w2`という完全なキャンセレーションが要る。既存の`Φ`存在補題
   (このセッション以前の大工事)と同規模の新規構築。
2. **`F_f(z,w)=z⟹w=0`の軽量版**(Y-線形係数が単位になることを使う)
   ——`coeff_single01_formalGroupLaw`は既にあるが、これを実際の
   評価に使うには「2変数冪級数をY次数でテイラー展開する」という
   一般論が要り、これも`aeval_subst_eq_aeval_aeval`と同規模。
3. **Lubin-Tate対数**`λ_f`(`λ_f(f(X))=π·λ_f(X)`で特徴づけられる、
   `K[[X]]`(`𝒪_K`ではなく`K`——`π`で割る係数を含む)の冪級数、
   `λ_f(X)=lim_n f^{(n)}(X)/π^n`)——`λ_f(F_f(X,Y))=λ_f(X)+λ_f(Y)`
   (加法性)・`λ_f([a]_f(X))=a·λ_f(X)`(線形性)を関数等式から導出。
   古典的に最も標準的な経路だが、収束性・関数等式の証明は
   「今セッション以前のLubin-Tate存在補題の大工事」と同規模。

いずれも今日1日で片付けられる規模ではなく、独立した大きめの
セッション(このセッション以前のLubin-Tate存在補題の構築に匹敵する
規模)が要ると判断した。次にこのtrackへ戻るときは、上記3つの
どれか1つを選んで正式に着手する——おそらく3(Lubin-Tate対数)が
最も古典的で見通しが良い。

これで「`a·x=0 ↔ π^n∣a`」(`x`が原始的なπ^n-捩れ点のとき)が
sorry無しで確立された——`𝒪_K/π^n≅Λ_n`の**核**が確定。
`|𝒪_K/π^n|=q^n=|Λ_n|`(既出`card_iteratedLubinTateTorsionPoints`)
と組み合わせれば、有限集合間の単射写像は同じ濃度なら全単射、という
一般論だけで**全単射性**(したがって同型そのもの)に到達できる
見通し——次の一歩はこの「単射→全単射」の橋渡しと、`𝒪_K`の乗法群
`(𝒪_K)^×`への制限、そして`Gal(K(Λ_n)/K)`との関係付け。

## ★★★戦略地図(2026-09-04、サブエージェント調査): この track が
どの `/goal` 項目を実際に支えるか

`/goal` の残り8項目それぞれの`.needs`(下界)を読み合わせた結果:

- **Proposition 1.1**: 局所 Tate 双対性(絶対不在)——**別系統**、この
  Lubin-Tate track とは無関係。
- **Proposition 1.2**: `Γ_K^ab≅(K^×)^`(局所類体論の相互律)そのもの
  ——**この track が直接埋める対象**。
- **Corollary 1.3**: `SubgroupCorrespondence`(開部分群の対応)経由で
  同じ相互律に依存——**この track が直接埋める対象**。
- **Proposition 2.1**: 表向きの`.needs`は「p進対数(解消済み)+
  Verlagerung(mathlib既存)+最後の1段」だけに見えるが、原文
  (`section-2.html`実測)を読むと「`U_K`を`Γ_K^ab`の部分群とみなす」
  という一文が実質的に相互律(の少なくとも単数群部分)を前提にしている
  ——**この track に実質的に依存する**(`.needs`の★下界には未だ
  明示されていない隠れた依存、CLAUDE.md「省略」の実例)。
- **Proposition 2.2**: 独立の絶対不在(上付き番号付け分岐群・
  Herbrand の定理、mathlib皆無)——**別系統**。ただし基礎対象
  `O_K̄`・`K̄^`の具体構成にはおそらくProp 2.1の構成を再利用する見通し。
- **Corollary 3.1**: Hodge-Tate 表現論(絶対不在)が主だが、docstring
  自身が「異なる`K`間の`Γ_K`同型の構成/反証」という**Prop 1.1/1.2/
  Cor1.3と同じ根本原因**にも触れている——**部分的にこの track**。
- **Corollary 3.3**: 「相互律に相当する非自明な数学的内容が要る」と
  docstringが明言——**この track に直接依存**(+独立のHodge-Tate
  `d_V(i)`判定基準も要る)。
- **Theorem 4.2**: 単射性のために`Γ_K^ab≅(K^×)^`を明示的に引用
  (§1 Prop 1.2と同じ箇所)——**この track に直接依存**。加えて
  `RamificationFiltration`に自然性公理が無いという**別の・新発見の**
  Interface設計上の穴もある(`memory/pgc-ramification-naturality-gap.md`)。

**結論**: 8項目中 **5項目**(Prop 1.2・Cor 1.3・Prop 2.1・Cor 3.1(部分)・
Cor 3.3・Theorem 4.2 ——数え方により5〜6)が、程度の差はあれこの
Lubin-Tate 相互律 track に懸かっている。残る真に独立した障害は
「局所 Tate 双対性」(Prop 1.1)・「上付き番号付け分岐群」(Prop 2.2)・
「Hodge-Tate 表現論」(Cor 3.1・Cor 3.3 の一部)——いずれも mathlib
絶対不在で、この track とは別に一から構築が要る。この track を
最後まで押し切る価値は高いと判断し、続行する。

## 次の技術的な壁: `a·x=b·x ⟹ π^n∣(a-b)`(単射性)に要る `F_f` の
「引き算」——回避できないと判明

前回「`F_f`加法の逆元・結合律は不要だった」と書いたが、それは
「**核**(`a·x=0`の場合)」に限った話——**一般の単射性**
(`a·x=b·x⟹π^n∣(a-b)`)には、`(a-b)·x`を`a·x`と`b·x`から復元する
何らかの手段が要り、これは`F_f`の**引き算**(の代わりになる何らかの
キャンセレーション)を避けられないと確認した。試した/検討した経路:

1. 一般の`F_f(z,w1)=F_f(z,w2)⟹w1=w2`(完全なキャンセレーション)
   ——`F_f`の形式逆元冪級数(陰関数定理型構成、既存の`Φ`存在補題と
   同規模の新規構築)か、Lubin-Tate対数(`λ(F_f(X,Y))=λ(X)+λ(Y)`
   となる冪級数、これも新規構築)のどちらかが要る——大掛かり。
2. 特殊ケース`F_f(z,w)=z⟹w=0`だけならもっと軽く済む見通しを発見:
   `H(Y):=F_f(z,Y)-z`を`Y`の1変数冪級数とみて(`z`はパラメータ)、
   その**Y-線形係数**`c(z):=1+Σ_{i≥1}(F_fのX^iY^1係数)z^i`が
   `‖z‖<1`のとき`‖c(z)‖=1`(単位、非零)になる(`Σ_{i≥1}`の項は
   `‖z‖<1`から`<1`に収まるので`1+(<1)`の項が支配項——非アルキメデス
   距離の性質)ことを示せれば、`H(w)=0∧w≠0`から`‖w‖=‖c(z)‖^{-1}‖
   (次数≥2の項)‖≤‖w‖^2`となり`‖w‖<1`と矛盾——`w=0`が出る。
   この経路は「2変数のGaussNorm一般論」(mathlib`MvPowerSeries.
   GaussNorm`——`HasGaussNorm`・`gaussNorm`・`le_gaussNorm`等、
   実在するが`aeval`との橋渡し定理は無い)より、既存の`aeval`の
   truncation-limit手法(このセッションで`aeval_subst_eq_aeval_aeval`
   の構築に使った手法と同型)を再利用する方が筋が良さそうと判断した
   ——ただし「1変数への制限(`X:=z`を代入してYだけ残す)」を厳密に
   冪級数として構成する部分がまだ具体化できていない。
   **まだ着手できていない、次に狙う技術的な要——`F_f`の逆元/対数を
   フルに作る大工事は避けられる見通しだが、この「Y-線形係数が単位」
   の構成自体は新しい仕事(このセッションでは完成させられなかった)。**

**続報(同日、`unit_action_mem_iteratedLubinTatePsiTorsionPoints`、
`AdjoinIntegers.lean`、commit`f9739556`)**: 上の技術的な壁(単射性
一般には`F_f`の引き算が要る)は残るが、**単射性を経由しない別の
入口**を見つけて前進した——`(𝒪_K)^×`が「原始的な捩れ点」の集合
(`ψ_nの根`)に**作用する**こと自体は、単射性抜きでも直接示せる:

★★★**`unit_action_mem_iteratedLubinTatePsiTorsionPoints`**:
`u∈(𝒪_K)^×`、`x`が`ψ_n`の根ならば`u·x`も`ψ_n`の根。証明は
`u·x∈Λ_{n-1}`と仮定して矛盾を導く背理法——`aeval_iteratedLubinTate
Distinguished_eq_zero`と同じ`Polynomial.hom_eval₂`橋渡しで
`D_{n-1}(u·x)=0`(`x`自身の`adjoinIntegers K x`のレベルで、`u·x`の
ための**新しい**adjoinIntegersは一切構築しない)を得て、
`[π^{n-1}]_f=D_{n-1}*U_{n-1}`から`(π^{n-1}*u)·x=0`、単位`u`を
`lubinTateAction_mul`で約分して`π^{n-1}·x=0`——前々回の primitivity
補題に矛盾。**`F_f`加法もキャンセレーションも一切不要**——乗法性
だけで閉じた、単射性の議論とは独立に成立する軌道の話。

これで「`(𝒪_K)^×`は`ψ_nの根`の集合に(閉じた)群作用を持つ」ことが
確立された。単純推移性(全単射性)にはまだ届かないが——`(𝒪_K)^×`
の軌道が`ψ_nの根`全体(単一の軌道)に一致すること、あるいは少なくとも
軌道が真部分集合でないこと——は依然として単射性・全射性いずれかの
議論を要する。次の一歩の選択肢: (a) 上述のY-線形係数論法を完成させる、
(b) 軌道の濃度(`|(𝒪_K)^×の軌道|=|(𝒪_K)^×|/|安定化群|`)を安定化群の
評価だけで`q^{n-1}(q-1)`に持ち上げられないか検討する(単射性を全射性
無しで済ませられる可能性)、のいずれか。

## 続報(2026-09-04、`norm_algEquiv_eq`、`AdjoinIntegers.lean`、
commit `277510eb`): Galois同変性の前提を独立に確立

単射性/全射性の3経路(`F_f`の形式逆元、Y-線形係数、Lubin-Tate対数)
のどれを選ぶにせよ、いずれ`Gal(K.closure/K.carrier)`が`Λ_n`に
どう作用するかという話に合流する。その最初の一歩として、**単射性/
全射性の議論そのものとは独立に**、`σ:K.closure≃ₐ[K.carrier]
K.closure`が`spectralNorm`の等長写像であること
(`spectralNorm K.carrier K.closure (σ x) = spectralNorm K.carrier
K.closure x`)を証明した——`norm_algEquiv_eq`。

証明の骨子: `σ x`は`x`と同じ`minpoly K.carrier x`の根になる
(`Polynomial.aeval_algHom_apply`で`aeval`とσの可換性を得てから
`minpoly.aeval`で`0`に潰す)。mathlibの
`spectralNorm.spectralMulAlgNorm_eq_of_mem_roots`(同じ最小多項式の
根はすべて等しい`spectralMulAlgNorm`を持つ、一般に`K⊆L⊆E`の3層で
書かれている)を`E=L=K.closure`という退化した場合に適用すると、
`algebraMap K.closure K.closure`が恒等写像になるので
`spectralMulAlgNorm_def`で展開したあと`simp`だけで閉じる——
`AlgEquiv.commutes'`のような一般の`E≠L`用の補題は不要だった
(最初の試みではこれを使おうとして失敗した——退化ケース特有の罠)。

このステップ自体は`a·x`(Lubin-Tate作用)に一切触れていない、純粋
なノルムの話であり、単射性/全射性のどの経路を最終的に選んでも
共通して要る土台になる。

## 続報(同日、`algEquiv_mem_iteratedLubinTatePsiTorsionPoints_of_mem`、
`AdjoinIntegers.lean`、commit `fc845a9b`): Galois群が`ψ_n`の根を保つ
——`spectralNorm`を経由しない、もっと直接的な経路

上記で「`spectralNorm`が保たれる⟹次数の一致から根の集合内に留まる」
という筋を検討していたが、実際にはもっと初等的な経路で足りると
気づいた:**`ψ_n`自身が`𝒪[K.carrier]`(したがって`K.carrier`の像)
に係数を持つ多項式である**ので、`σ:K.closure≃ₐ[K.carrier]
K.closure`が`K.carrier`を固定するという定義そのものから、
`x`が`ψ_n`の根なら`σ x`も同じ多項式の根になる——
`Polynomial.aeval_algHom_apply`(σと`aeval`の可換性)だけで閉じる。
`algebraMap 𝒪[K.carrier] K.closure`を`K.carrier`経由に分解する
`IsScalarTower.algebraMap_eq`は、既存の
`nodup_roots_iteratedLubinTatePsi_map`等と全く同じ手筋の再利用。

これで**`(𝒪_K)^×`と`Gal(K.closure/K.carrier)`の双方が、独立に、
「原始的な捩れ点全体」(`ψ_n`の根)という同じ有限集合上に閉じた
作用を持つ**ことが確立された(単数側は`unit_action_mem_
iteratedLubinTatePsiTorsionPoints`、既出)。`norm_algEquiv_eq`
(spectralNormのσ不変性)は依然として、`σ(a·x)=a·σ(x)`という
**Galois同変性そのもの**(2つの作用が可換であること)を示す際に
truncation-limit手法の中で使う見込みで、無駄にはならない——ただし
「σがψ_nの根の集合を保つ」という一歩自体は、それを経由せずに
直接示せた、という発見。

次の一歩:(a)`σ(a·x)=a·σ(x)`(Galois同変性そのもの)を
truncation-limit手法で示す、(b)あるいは単射性/全射性(3経路の
いずれか)を先に完成させる——どちらから手をつけても、最終的な
`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`には両方必要になる見通し。

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★ 続報(2026-09-04、
最大の技術的突破): 単射性が、`F_f`の減法・逆元・対数を一切経由せず
に完成した——当初の3候補はすべて不要だった

これまで「単射性(`a·x=b·x⟹π^n∣(a-b)`)には`F_f`の形式逆元・
Y-線形係数・Lubin-Tate対数のいずれかのフルな構成が要り、既存の
存在補題と同規模の新規構築になる」と繰り返し判断されてきたが、
実際には**もっと軽い第4の経路**で完全に閉じることが判明した。

**鍵となる洞察**: 一般の`F_f(z,w1)=F_f(z,w2)⟹w1=w2`(フルな
キャンセレーション)は不要——`a·x=b·x`という仮定と既存の加法公式
`lubinTateAction_add`(`(b+c)·x=F_f(b·x,c·x)`、`c:=a-b`)を
組み合わせると、`F_f(z,w)=z`という**特殊な**形の式(`z:=b·x`、
`w:=c·x`)が自動的に出てくる。この特殊形から`w=0`を出すだけなら、
`F_f`の単位元則(両側、`F_f(X,0)=X`・`F_f(0,Y)=Y`、既出)から
「次数`≤1`部分がちょうど`X+Y`」であることを使った**評価レベルの
不等式**`‖F_f(z,w)-z-w‖≤‖z‖*‖w‖`だけで十分——`z,w`の役割を
入れ替えて一般化する必要も、真の逆元を構成する必要も無い。

**新規に確立した内容(すべて`sorry`無し、ゲート通過済み)**:

1. **`coeff_single1_formalGroupLaw`**(`LubinTateIdentityLawSymm.lean`、
   commit`05619d4e`): `coeff_single0_formalGroupLaw`のX_0↔X_1対称版
   ——`F_f`の「`X_1`のみ」の係数の一般`n`版。既存の`formalGroupLaw_
   identity_left`(`F_f(0,Y)=Y`)から同じ手筋で出す。

2. **★★★★★★★★★★★★★★★★★★★★★★★★`norm_aeval_formalGroupLaw_sub_le`**
   (新ファイル`LubinTateFormalGroupLawEstimate.lean`、commit
   `05619d4e`): **`‖F_f(z,w)-z-w‖≤‖z‖*‖w‖`**(評価レベル、
   `‖z‖,‖w‖≤1`の任意の点)。`aeval_formalGroupLaw_eq_of_snd_eq_zero`
   と全く同じtruncation-limit手法の**不等式版**——各truncationの
   support から`X_0^1`・`X_1^1`の2つの特別な単項式(係数ちょうど1、
   `coeff_single0/1_formalGroupLaw`)を`Finset.add_sum_erase`で2回
   取り除き、残りの各単項式(すべて`i≥1∧j≥1`、さもなくば係数`0`)を
   `‖coeff‖≤1`(係数の埋め込み先での評価)・`‖z‖^i≤‖z‖`
   (`pow_le_of_le_one`、`i≥1`)から`‖z‖*‖w‖`で抑え
   (`IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg`)、
   最後に`le_of_tendsto`+`Filter.eventually_ge_atTop`で
   `n→∞`の極限へ持ち上げる。抽象的な環`S`(`NormedCommRing`+
   `IsUltrametricDist`+`IsLinearTopology`等)について述べており、
   `adjoinIntegers K x`に限らず再利用可能な一般補題。
   **重要な技術的注意**: `S:=adjoinIntegers K x`のまま`aeval`を
   含む文を直接ステートメントに書くと`CompleteSpace ↥(adjoinIntegers
   K x)`のインスタンス解決が原因不明に失敗する(`haveI`で明示的に
   注入しても直らない、しかし`haveI`を消しても後続が壊れない
   ——ステートメント自体の型検査時点で自動探索が先に失敗する
   ことが原因と判明)。**回避策**: 抽象的な`S`で述べ、呼び出し側
   (`lubinTateActionAtTorsionPoint_injective_of_eq`)で`S:=
   adjoinIntegers K x`に**具体化して適用する**(この具体化された
   呼び出しでは`haveI`で注入した具体的なインスタンスがちゃんと
   使われるので問題が起きない)——`aeval_formalGroupLaw_eq_of_
   snd_eq_zero`が最初からこの設計だった理由がここで腑に落ちた。

3. **★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
   `lubinTateActionAtTorsionPoint_injective_of_eq`**(新ファイル
   `LubinTateActionInjective.lean`、commit`05619d4e`):
   **`a·x=b·x⟹π^n∣(a-b)`**(`x`が原始的な`π^n`-捩れ点のとき)。
   `c:=a-b`、`lubinTateAction_add`で`a·x=F_f(b·x,c·x)`、仮定と
   合わせて`b·x=F_f(b·x,c·x)`。上の不等式を`(z,w):=(b·x,c·x)`に
   適用すると`‖c·x‖=‖F_f(b·x,c·x)-b·x-c·x‖≤‖b·x‖*‖c·x‖`。`b·x`は
   捩れ点で`‖b·x‖<1`(`spectralNorm_lt_one_of_mem_
   iteratedLubinTateTorsionPoints`)、`c·x≠0`と仮定すると
   `‖c·x‖≤‖b·x‖*‖c·x‖<‖c·x‖`で矛盾——`c·x=0`が出て、既存の核の
   特徴づけ(`lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_
   mem_iteratedLubinTatePsiTorsionPoints`)から結論。

4. **`unitActionQuotientLift_injective`**(`LubinTateActionInjective.lean`、
   commit`c87960d0`): 上を`u,v:(𝒪_K)^×`(単数)に適用し、
   `π^n∣(u-v)`から`u⁻¹*v∈principalUnits K π n`
   (`QuotientGroup.eq`+`mem_principalUnits_iff`+`u⁻¹*u=1`の
   純代数計算、`linear_combination`で閉じる)を導くだけ。

**現状**: `unitActionQuotientLift:(𝒪_K)^×⧸principalUnits K π n→
adjoinIntegers K x`は**単射**であることが確立された。既存の
`card_principalUnitsQuotient`=`card_iteratedLubinTatePsiTorsionPoints`
(`QuotientCardinality.lean`、`Nat.card`ベース)と合わせれば、同じ
有限濃度の集合間の単射——**全単射のはず**(有限集合の鳩の巣原理)。

**残る最後の一歩(次のセッションへの最有力候補)**: 「単射+濃度
一致⟹全単射」を`Function.Bijective`の言葉できちんとパッケージする
だけの**純粋な有限組合せ論の仕上げ**——新しい数学的困難は無い
見込みだが、`Nat.card`(quotient側)と`Finset.card`
(`iteratedLubinTatePsiTorsionPoints`側、`K.closure`に住む)を
橋渡しする型の付け替え(`adjoinIntegers K x`の元を`K.closure`へ
2段階coerceする)が必要で、今回は時間の都合で着手を見送った。
`unit_action_mem_iteratedLubinTatePsiTorsionPoints`(単数の作用は
`ψ_nの根`に留まる、既出)を使えば「像が`ψ_nの根`に含まれる」ことは
既に分かっているので、`Set.BijOn`/`Function.Injective.bijective`
系のmathlib補題(`Finite.injOn_iff_bijOn_of_mapsTo`等、要調査)を
探すところから始めるとよい。

**当初の3つの候補(`F_f`の形式逆元のフル構成・Y-線形係数論法・
Lubin-Tate対数)はすべて不要になった**——このセッションで実際に
機能した経路は、既存の加法公式`lubinTateAction_add`と、単位元則
から出る「次数`≤1`部分=`X+Y`」の評価レベル不等式を組み合わせる
だけの、遥かに軽い第4の道だった。今後Galois同変性の完成や他の
局面で同様の「大掛かりな構成を避ける軽い道」がないか、まず疑って
みる価値がある教訓。

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★ 続報
(2026-09-04、同日、セッション最終到達点): 「残る最後の一歩」も
即日で完成——`(𝒪_K)^×⧸principalUnits K π n ≃ ψ_nの根`(全単射)

上で「次のセッションへ」と記録した「単射+濃度一致⟹全単射」の
純粋な有限組合せ論の仕上げは、**同じセッション内で**完成した。
懸念していた「`Nat.card`と`Finset.card`の橋渡し」も、実際には
`Fintype.card_coe`+`Nat.card_eq_fintype_card`の2つの標準補題を
並べるだけで、想定より遥かに軽かった。

**新規に確立した内容(すべて`sorry`無し、ゲート通過済み)**:

1. **`finite_principalUnitsQuotient`**(`QuotientCardinality.lean`、
   commit`c1ff72c8`): `(𝒪_K)^×⧸principalUnits K π n`が有限。
   `principalUnitsQuotientEquiv`経由で`(𝒪_K/π^n)^×`の有限性
   (`Finite Rˣ`は`Finite R`から自動、`infer_instance`一発)を
   移すだけ(`Finite.of_equiv`)。

2. **`unitActionQuotientLift_mem_iteratedLubinTatePsiTorsionPoints`・
   `unitActionQuotientBijOn`**(新ファイル`LubinTateActionBijective.lean`、
   commit`c1ff72c8`): `unitActionQuotientLift`の値域を`ψ_nの根`
   (`Finset`の`Sort`強制)へ制限した写像。値域の制限は既出の
   `unit_action_mem_iteratedLubinTatePsiTorsionPoints`を
   `QuotientGroup.induction_on`で一般の`U`へ持ち上げるだけ。

3. **`unitActionQuotientBijOn_injective`**: 上の制限写像の単射性
   ——`unitActionQuotientLift_injective`+「`adjoinIntegers K x→
   K.carrier⟮x⟯→K.closure`の二重coeが単射」(`Subtype.coe_
   injective`を2回)。

4. **★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
   `unitActionQuotientBijOn_bijective`**: **`(𝒪_K)^×⧸
   principalUnits K π n ≃ ψ_nの根`という全単射**。
   `Fintype.card_eq.mp`(濃度が等しければ`Nonempty(α≃β)`)+
   `Function.Injective.surjective_of_finite`(そのequivの証拠が
   あれば単射⟹全射)を組み合わせるだけ——`card_principalUnits
   Quotient`=`card_iteratedLubinTatePsiTorsionPoints`
   (ともに`q^n-q^{n-1}`、既出)を`Fintype.card_coe`+
   `Nat.card_eq_fintype_card`で橋渡しする。

**到達点のまとめ**: `principalUnitsQuotientEquiv`
(`(𝒪_K)^×⧸principalUnits K π n≃*(𝒪_K/π^n)^×`、既出)と組み合わせ
れば、**`(𝒪_K/π^n)^×≃ψ_nの根`**——古典的なLubin-Tate理論の核心的な
事実(原始的な`π^n`-捩れ点の全体が`(𝒪_K/π^n)^×`と自然に一対一
対応すること)が完全に確立された。`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`
という主定理に向けて、**単数側の道はここで完結**した。

**残る唯一の大きな柱**: Galois側——`Gal(K.closure/K.carrier)`が
`ψ_nの根`を置換として保つこと(`image_algEquiv_iteratedLubinTatePsi
TorsionPoints`、既出)は分かっているが、`Gal(K.closure/K.carrier)`
の各`σ`に対して「`σ(x)=u_σ·x`となる一意な単数`u_σ`」を対応させる
写像(これが目的の同型`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`の実体)を構成
するには、いま確立した全単射(単数の作用が`ψ_nの根`上で単純推移的)
を使えばよい——`u_σ`の一意存在性は今回の全単射から直ちに従う
(`ψ_nの根`上の`(𝒪_K)^×`の作用が全単射的なので、任意の2点間に
一意な単数が存在する)。真に残っているのは、この対応`σ↦u_σ`が
**準同型**であること(`u_{στ}=u_σ*u_τ`)を示す部分——これには
`σ(a·x)=a·σ(x)`(Galois同変性そのもの、cross-point instance
bridgingという既知の難所)が必要になる見通し。次のセッションでは
まずこの同変性の構成、あるいはそれを回避する別の道が無いか
(このセッションで単射性の壁を回避できたのと同様に)を最初に
検討する価値がある。

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★ 続報
(同日、さらに続き): 相互写像`σ↦u_σ`が集合論的に完成——
`reciprocityMap`として具体的な名前付き対象になった

上の「`u_σ`の一意存在性」を実際に定理として仕上げ、さらに
`Classical.choose`で具体的な写像として定義した(すべて`sorry`無し、
ゲート通過済み、commit`a988ffda`・`82c3f6bd`):

1. **`existsUnique_unitActionQuotient_eq_algEquiv`**
   (`LubinTateActionBijective.lean`): `∀σ, ∃!U, U·x=σ(x)`。
   `σ(x)∈ψ_nの根`(`algEquiv_mem_iteratedLubinTatePsiTorsionPoints_
   of_mem`、既出)+`unitActionQuotientBijOn`の全単射性(直前)を
   組み合わせるだけ——新しい数学的困難は無い、既存結果の直接の系。

2. **`reciprocityMap`**: 上の`ExistsUnique`の`Classical.choose`に
   よる定義——`Gal(K.closure/K.carrier)→(𝒪_K)^×⧸principalUnits
   K π n`という、目的の同型の実体となる具体的な写像。

3. **`reciprocityMap_spec`**: 定義性質`reciprocityMap σ · x = σ(x)`。

4. **`reciprocityMap_one`**: 整合性チェック
   `reciprocityMap (AlgEquiv.refl) = 1`——`1·x=x`(既出)と一意性
   から、Galois同変性を一切経由せず独立に示せる。

**現状のまとめ**: `Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`という主定理の
**写像そのもの**(`reciprocityMap`)が、`sorry`無しで具体的に
構成された。残るのは**準同型性**
(`reciprocityMap(σ*τ)=reciprocityMap σ * reciprocityMap τ`)の
証明だけ——これにはGalois同変性`σ(a·x)=a·σ(x)`(cross-point
instance bridging)が必要になる見通し(全単射から`u_σ`の存在は
出るが、`u_{στ}=u_σ*u_τ`を出すには`σ(u_τ·x)=u_τ·σ(x)`という
同変性の1インスタンスが要る——`lubinTateAction_mul`(乗法性、既出)
と組み合わせれば`(στ)(x)=σ(u_τ·x)=u_τ·σ(x)=u_τ·(u_σ·x)=(u_τ*u_σ)·x`
となり、単数群の可換性で`u_{στ}=u_σ*u_τ`が出る、という計算の
筋道自体は明確)。次のセッションでは、この「`a=u_τ`という**単数**
限定のGalois同変性」だけで十分なので、一般の`a∈𝒪_K`に対する
フルな同変性より的を絞った、より軽い経路がないか(単射性で
そうだったように)を最初に検討する価値がある。

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
続報(同日、セッション最大の到達点): **Galois同変性 `σ(a·x)=a·σ(x)`
が cross-point instance bridging を一切経由せず完成した**

上で「一般のaに対するフルな同変性より的を絞った軽い経路」を検討
すると書いたが、実際には**一般の`a∈𝒪_K`に対する完全な同変性**が、
cross-point bridgingという壁ごと回避する形で、丸ごと確立できた。

**鍵となる洞察**: `adjoinIntegers K x`と`adjoinIntegers K (σx)`
という異なる2つの環を橋渡しする必要があると警戒されてきたが、実は
**σ(x)自身がK.carrier⟮x⟯の中に留まる**(単数の作用の全射性、
`unitActionQuotientBijOn`、既出)という事実だけで、この「2つの環」
問題自体が消える——`σ`は`K.carrier⟮x⟯`を自分自身へ写す自己同型に
制限でき(単射性+有限次元性から)、その制限された自己同型を
`adjoinIntegers K x`自身への自己写像として扱えば、`σ(x)`を
「`x`自身の座標系の中の1点」として直接扱える。`adjoinIntegers
K (σx)`という別の環は一度も構築しない。

**新規に確立した内容(すべて`sorry`無し、ゲート通過済み、commit
`1341fcfd`・`b6c1c88e`)**:

1. **`algHom_aeval_powerSeries_comm`**(新ファイル`PowerSeriesAeval
   Comm.lean`): 連続な代数準同型`σ:S→ₐ[A]S`は`PowerSeries.aeval`と
   可換(`σ(aeval hz g)=aeval hz' g`)。Lubin-Tate理論に一切依存
   しない純粋な位相環論の一般補題——`Polynomial.aeval_algHom_apply`
   (多項式レベルの可換性)+`PowerSeries.WithPiTopology.tendsto_
   trunc_atTop`(truncationの収束)+`tendsto_nhds_unique`。
   `aeval_formalGroupLaw_eq_of_snd_eq_zero`と同じ手法の、より単純な
   (等式のみ、不等式を経由しない)版。

2. **`iteratedLubinTatePsiTorsionPoints_subset_adjoin`**(新ファイル
   `LubinTateActionEquivariance.lean`): `ψ_nの根`は`K.carrier⟮x⟯`
   に含まれる——単数の作用の全射性から、任意の`ψ_nの根`の元は
   `u·x`の形に書け、常に`adjoinIntegers K x⊆K.carrier⟮x⟯`に入る。

3. **`algEquiv_map_adjoin_eq`**: `IntermediateField.map σ.toAlgHom
   K.carrier⟮x⟯=K.carrier⟮x⟯`。`σ(x)∈K.carrier⟮x⟯`
   (上記+Galois保存性)から`≤`が`IntermediateField.adjoin_map`+
   `IntermediateField.adjoin_simple_le_iff`で出て、単射性+有限
   次元の次元の等式(`IntermediateField.eq_of_le_of_finrank_eq`)
   で`=`まで持ち上がる。

4. **`algEquivRestrictSelf`**: `σ`を`K.carrier⟮x⟯`自身への自己
   同型に制限したもの(`IntermediateField.equivMap`+
   `IntermediateField.equivOfEq`の合成)。`coe_algEquivRestrictSelf`
   (座標は`σ`そのもの、`rfl`)・`norm_algEquivRestrictSelf`
   (`norm_algEquiv_eq`からノルムを保つ)も確立。

5. **`adjoinIntegersRestrictSelf`・`adjoinIntegersRestrictSelfAlgHom`**:
   ノルムを保つことから`adjoinIntegers K x`(ノルム`≤1`の部分環)
   への制限が構成でき、加法・乗法・単位元・`algebraMap`との可換性
   (いずれも`K.carrier⟮x⟯`の中でのσ自身の準同型性に帰着)から
   `𝒪_K`-代数準同型として束ねられる。

6. **`continuous_algEquivRestrictSelf`・`continuous_adjoinIntegers
   RestrictSelfAlgHom`**: `K.carrier⟮x⟯`が`K.carrier`上有限次元
   なので、その上の`K.carrier`-線形自己写像は自動的に連続
   (`LinearMap.continuous_of_finiteDimensional`)——`K.closure`
   自身の無限次元性を経由する必要が無い。

7. **`hasEval_adjoinIntegersRestrictSelfAlgHom_mk`**: `σ(x)`も
   位相的冪零(Galoisが捩れ点全体を保つことから)。

8. **`lubinTateActionAtAlgEquivPoint`**: `a`を`σ(x)`で評価した
   もの——`adjoinIntegers K (σx)`を一切経由せず、`x`自身の座標系の
   中だけで計算する。

9. **★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
   `algEquiv_lubinTateActionAtTorsionPoint_comm`**: **`σ(a·x)=
   a·σ(x)`**(Galois同変性そのもの)。`algHom_aeval_powerSeries_comm`
   を`adjoinIntegersRestrictSelfAlgHom`に適用するだけで閉じる。

**到達点のまとめ**: `Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`という主定理を
証明するための数学的道具が**すべて**揃った——単射性・全単射性・
`reciprocityMap`の存在・そしてGalois同変性。残るのは
`reciprocityMap`の**準同型性**という最後の1つの計算だけ:
`(στ)(x)=σ(τ(x))=σ(u_τ·x)`(τ(x)=u_τ·x、`reciprocityMap_spec`)
`=u_τ·σ(x)`(**新しい同変性定理**、`a:=u_τ`)`=u_τ·(u_σ·x)`
(`σ(x)=u_σ·x`、`reciprocityMap_spec`)`=(u_τ*u_σ)·x`
(`lubinTateAction_mul`)`=(u_σ*u_τ)·x`(`𝒪_K`の可換性)。
`existsUnique_unitActionQuotient_eq_algEquiv`の一意性から
`reciprocityMap(στ)=u_σ*u_τ=reciprocityMap σ*reciprocityMap τ`。

**この最後の計算に要る技術的な注意点**(次のセッションで着手する際
の具体的な足がかり): (a) `lubinTateActionAtAlgEquivPoint K...σ u_τ`
(σ(x)での`u_τ`の評価)を`lubinTateActionAtTorsionPoint K...u_σ`
(すなわち`u_σ·x`、`x`自身での評価)へ書き換える必要があるが、
これは両者の**coeが等しい**(ともに`σ(x)`)ことから
`Subtype.coe_injective`で「`adjoinIntegers K x`の元として等しい」
へ持ち上げてから、`lubinTateEvalAtPoint`側の`HasEval`証明項の
食い違いを`lubinTateEvalAtPoint_congr`(既出、まさにこの用途の
ために用意されていた補題)で吸収する、という手筋が必要になる
見込み——「異なる証明項を持つ同じ値」問題は本セッションで何度も
遭遇した`Fin 2 ite`・依存引数`rw`の罠と同系統なので、
`tools/lean-idioms.md`のパターンをまず参照すること。

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
続報(同日、セッションの最終到達点・集大成): **`reciprocityMap`の
準同型性が完成し、Lubin-Tate相互律の核心部分が完結した**

上記「この最後の計算」を、まさにその足がかり通りの手筋
(`lubinTateEvalAtPoint_congr`+`lubinTateAction_mul`)で、
**同じセッション内で**完成させた(新ファイル
`LubinTateReciprocityHomomorphism.lean`、commit`4776a0d2`)。

**`reciprocityMap_mul`**: **`reciprocityMap(σ*τ)=reciprocityMap σ*
reciprocityMap τ`**。証明の骨子は事前に見立てた通り: `u_σ,u_τ`
(単数代表元)を取り、`(σ*τ)(x)=σ(τ(x))=σ(u_τ·x)`に**Galois同変性**
(`algEquiv_lubinTateActionAtTorsionPoint_comm`)を`a:=u_τ`で適用
すると`σ(u_τ·x)="u_τ·σ(x)"`。`σ(x)=u_σ·x`に対応する
`adjoinIntegersRestrictSelfAlgHom σ`の値と`u_σ·x`が(座標が一致
するので)同じ`adjoinIntegers K x`の元であることを
`lubinTateEvalAtPoint_congr`で処理し(`Prop`の証明無関係性のおかげ
で、`▸`で書き換えた`HasEval`証明項が`lubinTateAction_mul`の
LHSの証明項と`exact`で直接マッチする——`rw`ではなく`exact`を使う
のが鍵)、`lubinTateAction_mul`で`"u_τ·σ(x)"=(u_τ*u_σ)·x`まで計算。
一意性(`existsUnique_unitActionQuotient_eq_algEquiv`)と`𝒪_K`の
可換性(`mul_comm`、`ring`ではなく明示的に——`Units`の乗法は
`CommRing`の`ring`タクティクの対象外)から結論する。

**AlgEquivの積の順序に関する注意**: `(σ*τ) x = σ (τ x)`(`rfl`で
確認済み)——`Function.comp`と同じ「右から先に適用」の順序。

**到達点の総括(このセッション全体)**: `Gal(K.closure/K.carrier)→
(𝒪_K)^×⧸principalUnits K π n`という写像(`reciprocityMap`)が、
①well-defined(存在・一意性)、②(`ψ_nの根`への制限を経由して)
全単射、③**群準同型**——という3つの性質すべてを`sorry`無しで
兼ね備えることが証明された。`principalUnitsQuotientEquiv`
(`(𝒪_K)^×⧸principalUnits K π n≃*(𝒪_K/π^n)^×`、既出)と組み合わ
せれば、古典的なLubin-Tate理論の主定理
**`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`**へ向けた数学的な道具が完全に
揃った状態にある(`Gal(K.closure/K.carrier)`から`Gal(K(Λ_n)/K)`
への制限・核の同定という、標準的だが未着手の「体拡大のGalois理論
の一般論」への還元だけが残る——新しい数学的困難ではなく、既存の
mathlibのGalois理論API`IsGalois`/`AlgEquiv.restrictNormal`系を
探索するだけの仕上げ作業になる見込み)。

このセッションで、当初「`F_f`の形式逆元・Y-線形係数・Lubin-Tate
対数のいずれかのフルな構成が要る」「異なる2つの環を橋渡しする
cross-point instance bridgingが避けられない」と繰り返し警戒されて
きた**2つの最大の技術的難所**を、どちらも**当初の想定より遥かに
軽い経路**(前者は既存の加法公式+評価不等式、後者はσ(x)自身が
同じ座標系に留まるという事実の活用)で突破し、Lubin-Tate相互律の
核心部分を完成させた——本プロジェクトのこのトラックにおける
最大の到達点。

## 続報(同日、さらに続き): `reciprocityMap`の**全射性**——真に残る
最後の柱を特定、道筋は明確化したが instance diamond で今回は
完成に至らず(honest記録)

上で「残るはGalois理論の標準的な仕上げ」と書いたが、より正確に
掘り下げると、真に欠けている事実は**`reciprocityMap`の全射性**
(同値に:`Gal(K.closure/K.carrier)`が`ψ_nの根`に**推移的**に作用
すること)だと判明した——これは`unitActionQuotientBijOn`の全単射性
(単数側)とは**独立**の、Galois側だけの追加の事実で、まだ確立
していない。

**朗報**: この全射性を出すのに必要な核心的な事実——**`ψ_n`が
`𝒪_K`上既約**——は、実は**既に確立済み**だった
(`irreducible_iteratedLubinTatePsi`、`LubinTateActionPsi.lean`、
以前のセッションで完成)。`LubinTatePsiNorm.lean`(`spectralNorm_
root_iteratedLubinTatePsi`の証明)には、これをGaussの補題で
`K.carrier`上へ持ち上げ(`Polynomial.IsPrimitive.irreducible_iff_
irreducible_map_fraction_map`)、`minpoly.eq_of_irreducible_of_
monic`で「`ψ_n`の根`x`の`minpoly K.carrier x`は`ψ_n`自身に一致
する」ところまで**既に橋渡し済み**のコードがある。

**残る論証**: `x,y`が`ψ_n`の(異なる)根なら`minpoly K.carrier x =
ψ_n = minpoly K.carrier y`(同じ既約多項式)なので、標準的なGalois
理論により`σ(x)=y`となる`σ∈Gal(K.closure/K.carrier)`が存在する
(「同じ最小多項式を持つ2元は共役」という古典的事実)。

**この標準的事実に要るmathlib補題を特定した(すべて実在確認済み)**:
- `AdjoinRoot.liftAlgHom p (i:R→ₐ[S]T) (x:T) (heval) : AdjoinRoot p→ₐ[S]T`
  (根を`y`に送る準同型)
- `IntermediateField.adjoinRootEquivAdjoin F hx : AdjoinRoot(minpoly F x)
  ≃ₐ[F] F⟮x⟯`(`AdjoinRoot`と`F⟮x⟯`の同一視)
- 上記2つを合成して`φ0:K.carrier⟮x⟯→ₐ[K.carrier]K.closure`
  (`x↦y`)を得る
- `IsAlgClosed.lift`(algebraic拡大からalg閉体への埋め込みの存在)
  で`φ0`を`K.closure`全体へ延長する
- **`Algebra.IsAlgebraic.bijective_of_isScalarTower' φ`**:
  `K.carrier`-代数自己準同型`φ:K.closure→ₐ[K.carrier]K.closure`は
  **自動的に全単射**(algebraic拡大の自己埋め込みは常に全単射、
  という一般論)——これで`AlgEquiv.ofBijective`により目的の`σ`を
  作れる。

**今回ぶつかった技術的な壁(instance diamond、未解決)**: `φ0`を
使って`K.closure`に`Algebra K.carrier⟮x⟯ K.closure`という
(φ0経由の)**別の**代数構造を`letI`で局所的に与えようとしたが、
`K.closure`は`IntermediateField K.carrier K.closure`の元として
`K.carrier⟮x⟯`から**既に自然な`Algebra`インスタンスを持つ**ため、
`letI`で上書きしようとすると型クラス解決が衝突する(本セッション
序盤で対処した`CompleteSpace`インスタンス絡みの罠や、過去セッション
の`Valued`インスタンスダイアモンドと同系統の問題)。`IsAlgClosed.
lift`を`R:=K.carrier⟮x⟯`で直接使おうとする代わりに、**`K.carrier⟮x⟯`
と可換な、独立な型シノニム**(あるいは`AdjoinRoot(minpoly K.carrier
x)`自身をRとして使い、`IntermediateField.adjoinRootEquivAdjoin`の
合成を`IsAlgClosed.lift`の**後**に行う、という順序の入れ替え)を
検討すべき——`AdjoinRoot p`は`K.closure`の部分型ではない独立した
型なので、それを`R`として`IsAlgClosed.lift`(`S:=AdjoinRoot p`、
`M:=K.closure`、`Algebra R M`は`φ0'`(`AdjoinRoot p→ₐ[K.carrier]
K.closure`)経由で作る)を適用すれば、既存の`Algebra L K.closure`
インスタンスと衝突しない可能性が高い——次回はこちらから試すこと。

**現状**: ファイルへの書き込み・コミットは無し(REPL内の探索のみ、
正直に記録)。次のセッションでの最有力候補は、この「`AdjoinRoot`を
`R`として`IsAlgClosed.lift`を適用する」経路の実装。これが通れば
`reciprocityMap`の全射性(ひいては`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`の
完全な証明)に到達する見込み。

## 続報(同日、`algEquiv_mem_iteratedLubinTateTorsionPoints_of_mem`、
`AdjoinIntegers.lean`、commit `93108293`): `Λ_n`全体もσで保たれる

上の`ψ_n`の根版と全く同じ証明パターン(`D_n`も`𝒪[K.carrier]`に
係数を持つので、σがK.carrierを固定することから直接出る)を
`iteratedLubinTateDistinguished`(`D_n`)に適用するだけで、
`Gal(K.closure/K.carrier)`が`Λ_n`(π^n-捩れ点全体、`ψ_nの根の
帰納的合併`を経由する必要すら無い)そのものを保つことが出た。

**(a)の`σ(a·x)=a·σ(x)`(Galois同変性そのもの)についての判断**:
`a·x`は`x`固有の`adjoinIntegers K x`という座標系で定義されており、
`a·σ(x)`は`adjoinIntegers K (σ x)`という別の座標系になる——
これらをK.closure経由で橋渡しする「cross-point instance bridging」
は、本セッション中に`unit_action_mem_iteratedLubinTatePsiTorsion
Points`の証明で意図的に**避けてきた**技術的な難所(以前の`Valued`
インスタンス衝突と同系統のリスク)であり、着手すると長時間の
不生産的な格闘になる可能性が高いと判断し、今回は見送った。
代わりに(σがψ_n根・Λ_n全体を保つという)より安全に完成できる
2つの補題を確立した——これらは「cross-point bridging」を一切
要さない、`x`自身の中だけで閉じる議論だったため。

現状のまとめ: `(𝒪_K)^×`と`Gal(K.closure/K.carrier)`の双方が、
独立に、`Λ_n`・`ψ_nの根`という同じ有限集合上に閉じた作用を持つ
ことが確立された。次の一歩は、単射性/全射性(3経路のいずれか)を
先に完成させるほうが、Galois同変性(cross-point bridgingを要する)
より安全な優先順位と判断する。

## 続報(同日、`image_algEquiv_iteratedLubinTatePsiTorsionPoints`・
`image_algEquiv_iteratedLubinTateTorsionPoints`、`AdjoinIntegers.lean`、
commit `5a85063e`): 「保つ」を「並べ替える」に強化

上の2つの「σは写像として保つ」(`⊆`)を、σとσ.symmの**両方**に
適用して両側から挟むだけで、`Finset.image σ S = S`という強い形
(単射性・濃度の議論は一切不要——`σ.symm`も同じ`AlgEquiv over
K.carrier`なので同じ補題がそのまま適用できる)が出た。これで
`Gal(K.closure/K.carrier)`が有限集合`ψ_nの根`・`Λ_n`上に**置換**
として作用することが確立された——将来`Gal(K.closure/K.carrier)
→ Equiv.Perm (ψ_nの根)`のような準同型を構成する際、直接使える形。

これで今回のセッションでの「Galois側の作用の基礎固め」は
ひとまず一区切りとし、次はcross-point bridgingを要するGalois
同変性そのもの、または単射性/全射性(3経路)のどちらかに進む
判断になる。後者(特にY-線形係数の経路)は、`F_f`の次数1部分が
`X+Y`であること自体は既に`coeff_single01_formalGroupLaw`・
`coeff_single10_formalGroupLaw`(`LubinTateIdentityLawSymm.lean`・
`LubinTateCommutativity.lean`、既出)で確立済みと判明——形式的な
土台の半分は既にある。残るのは「評価レベルでΣ_{i≥2}の項の
ノルムが`<1`に収まる」という非アルキメデス的な極限の見積もり
(既存の`aeval_formalGroupLaw_eq_of_snd_eq_zero`と同系統の
truncation-limit手法の再利用が見込めるが、まだ着手していない)。
