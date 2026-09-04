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
