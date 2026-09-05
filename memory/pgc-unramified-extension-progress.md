---
name: pgc-unramified-extension-progress
description: 不分岐拡大理論(Gal(K^ur/K)≅Ẑ)の構築進捗——節目(5)をProposition 1.2/2.1へ接続するための第二の柱。2026-09-05にスケルトン(定義のみ)着手
metadata:
  type: project
---

`memory/pgc-prop12-reciprocity-gap.md`で確認したとおり、節目(5)
(`pgc-lubin-tate-existence-progress.md`——古典的Lubin-Tate理論の
射影極限、完全分岐アーベル拡大の部分は完成済み)だけでは
`Γ_K^ab≅(K^×)^∧`(ひいてはProposition 1.2・2.1)に届かない——
**不分岐拡大理論**(`K^ur`・`Gal(K^ur/K)≅Ẑ`)という、mathlib・本
プロジェクトともに未着手の第二の柱が要る。本ファイルはこの構築の
進捗を記録する。ユーザーの承認(2026-09-05、AskUserQuestionで確認)
を得て着手。

## 数学的な見取り図

1. **不分岐性の定義**: `x∈K.closure`(`K(x)/K`有限)が不分岐 ⟺
   剰余次数`f=[残余体(adjoinIntegers K x):残余体(𝒪_K)]`が
   `[K(x):K]`に一致(`e·f=[K(x):K]`のうち`e=1`)。
2. **存在性**(各次数`n`の不分岐拡大が存在): 分離的多項式の分解体は
   不分岐、という一般論(`X^{q^n}-X`または`X^{q^n-1}-1`の分解体)+
   `HenselianLocalRing`(完備局所環はHenselian)による`𝒪_K`への
   持ち上げ。
3. **一意性**: 次数`n`の不分岐拡大は`K.closure`内でただ1つ——残余体
   の拡大`F_{q^n}/F_q`の一意性(`FiniteField.algEquivOfCardEq`)から。
4. **`K^ur`とGalois群**: `K^ur:=⋃n(次数nの不分岐拡大)`、
   `Gal(K^ur/K)≅Gal(F̄_q/F_q)`(剰余体還元)、`Gal(F̄_q/F_q)≅Ẑ`
   (Frobenius`x↦x^q`が位相的生成元)。

## 調査済みのmathlib道具(2026-09-05実測)

- `HenselianLocalRing`(`Mathlib/RingTheory/Henselian.lean:108`)——
  Hensel's lemmaの標準的特徴づけ(`HenselianLocalRing.TFAE`)。
  `𝒪[K.carrier]`が`IsAdicComplete`(既出)からこのinstanceを持つ
  はず——**本プロジェクトでの確認はまだ**。
- `Ideal.ramificationIdx`・`Ideal.inertiaDeg`(`Mathlib/NumberTheory/
  RamificationInertia/`)——Dedekind整域一般の分岐理論。
  `Algebra.isUnramifiedAt_iff_map_eq`(`RingTheory/Unramified/
  LocalRing.lean:134`)が`IsUnramifiedAt ↔ 分離性∧極大イデアルが写る`
  という特徴づけを与える——`𝒪[K.carrier]→adjoinIntegers K x`への
  適用はまだ試していない。
- `FiniteField.isSplittingField_sub`(`X^{Card K}-X`の分解体が
  `K`自身)・`GaloisField`(`F_{p^n}`の標準構成)・`FiniteField.
  algEquivOfCardEq`(同じ濃度の有限体は同型)——残余体側の道具は
  比較的揃っている。
- `InfiniteGalois`(`FieldTheory/Galois/Profinite.lean`、既出、
  節目(5)のコンパクト性論法で使用)——`Gal(F̄_q/F_q)`の副有限群
  構造自体はこの一般論の特殊化として得られる見込みだが、
  「≅Ẑ」という**具体的な**同型はまだ探していない
  (`ZMod`・`ProfiniteGrp`絡みの検索が必要)。
- **mathlibに不在と確認済み**: 「有限体の代数閉包のGalois群≅Ẑ」
  という直接の定理・「局所体の不分岐拡大」という文脈での`Frobenius`
  ・`Gal(K^ur/K)`(`NumberField.InfinitePlace`関連の`IsUnramified`は
  archimedean placeの話で無関係)。

## 現時点の進捗(2026-09-05、スケルトンの第一歩)

`Found/PGC/UnramifiedExtension.lean`(新規、sorry無し・現時点では
**定義のみ**):
- `residueDegree K x`: `adjoinIntegers K x`の剰余体の濃度
  (`Nat.card`、有限とは限らないので`0`になる余地を許容)。
- `finiteDimensional_adjoin_zero`: 退化した基点`x=0`での
  `FiniteDimensional`instanceの健全性チェック(`K(0)=K`、`⊥`)。

まだ**不分岐性のProp自体を切り出していない**——次の一歩は:
1. `IsUnramifiedPoint K x : Prop := residueDegree K x = (𝒪_Kの剰余体
   濃度) ^ (finrank K.carrier (K(x)))`を定義する。
2. `x=0`(または`x∈range(algebraMap K.carrier K.closure)`一般)が
   この条件を満たすこと(`finrank=1`・`residueDegree`が`𝒪_K`自身の
   剰余体濃度に一致すること、後者には`adjoinIntegers K 0≃+*𝒪[K.
   carrier]`という環同型が要る——まだ構築していない)を確認する。
3. `HenselianLocalRing (𝒪[K.carrier])`のinstanceを実際に確認・
   構築する(`IsAdicComplete`から従うはずの一般論を探すか、
   `IsDiscreteValuationRing`+完備性から直接示す)。
4. 存在性(次数`n`の不分岐拡大の構成)へ進む。

## 次に戻るときの最優先タスク

上記2(退化した基点の不分岐性確認、`adjoinIntegers K 0≃+*𝒪[K.
carrier]`の構築)が最も具体的で低リスクな次の一歩。3(Henselian
instance)と並行して進められる見込み。

★2026-09-05続報(完成): 退化した基点`x=0`の不分岐性チェックが
**sorry無しで完成した**——`Found/PGC/UnramifiedExtension.lean`に
`adjoinIntegersZeroRingHom`・`adjoinIntegersZeroRingHom_bijective`・
`adjoinIntegersZeroEquiv`(`𝒪[K.carrier]≃+*adjoinIntegers K 0`)・
`residueDegree_zero`(`residueDegree K 0=Nat.card(ResidueField 𝒪_K)`)
をすべて確立しcommit済み。`[K(0):K]=1`(`finiteDimensional_adjoin_
zero`、既出)と合わせ、退化した基点が「剰余次数`=`剰余体濃度の
`[K(x):K]`乗」という不分岐性の条件を自明に満たすことが確認できた。
新しい罠は無かった——唯一、`residueDegree K x`の定義が`[FiniteDimensional
...]`を暗黙のinstance引数に取るため、`residueDegree_zero`の**主張
自体**にも同じinstance引数を明示的に添える必要があった(本セッション
で何度も遭遇した「STATEMENT自身の型検査にinstanceが要る」パターン、
`haveI`では手遅れ)。

以下は当時の発見の記録(上記の完成に使われた): 「`K.closure`上の
spectralNormが`K.carrier`上の元では元のノルムに一致する」という
橋渡しは、**既に`norm_algebraMap' K.closure k : ‖algebraMap K.carrier
K.closure k‖ = ‖k‖`としてmathlibに存在する**(`exact?`で即発見、
`spectralNorm_extends`を経由する必要すら無かった)ことを確認した。
これで`adjoinIntegers K 0`(`{y:K.carrier⟮0⟯|‖y‖≤1}`)と`𝒪[K.carrier]`
(`{y:K.carrier|‖y‖≤1}`)が、`IntermediateField.botEquiv K.carrier
K.closure`(`(⊥:IntermediateField K.carrier K.closure)≃ₐ[K.carrier]
K.carrier`、既存)とノルムの両立性を組み合わせれば対応することの
**数学的な裏付けは揃った**——残るのは`IntermediateField.adjoin
K.carrier {0}=⊥`という**型レベルの書き換え**を経由して`adjoinIntegers
K 0`(`IntermediateField.adjoin K.carrier {0}`上のSubringとして定義
されている)を`⊥`上のSubringとして正しく`show`/`cast`する技術的な
組み立てのみ(本セッションのLubin-Tate作業で何度も遭遇した「型注釈と
既存の項の食い違い」と同種の罠が予想される——`tools/lean-idioms.md`
の該当エントリを先に見直すこと)。

★上記は完成した(冒頭の「続報(完成)」参照)——予想された罠は
**実際には起きなかった**(直接`RingHom`を構成する方式を採ったため、
`Subring`同士の`cast`という難しい経路を最初から避けられた)。

## 次に戻るときの最優先タスク(更新、2026-09-05)

退化した基点(`x=0`、次数1)の健全性チェックが完成したので、次は
**次数`n≥2`の不分岐拡大の存在性**——本格的な数学的内容が要る最初の
関門:
1. `HenselianLocalRing (𝒪[K.carrier])`のinstanceを確認・構築する
   (`IsAdicComplete`から従う一般論を探すか、`IsDiscreteValuationRing`
   +完備性から直接示す)——これがまだ手つかず。

   ★2026-09-05実測(重要): **mathlibには「完備局所環はHenselian」
   という一般定理もinstanceも無い**——`Mathlib/RingTheory/Henselian.
   lean`にあるのは`HenselianLocalRing`の**定義**と特徴づけ
   (`HenselianLocalRing.TFAE`)のみで、instanceを生成する補題はゼロ件
   (`grep "instance.*Henselian"`で確認)。一方`Mathlib/NumberTheory/
   Padics/Hensel.lean`には**`ℤ_[p]`専用の**Hensel's lemma(Newton法の
   収束、`newton_cau_seq`・`ih_gen`等の補助定義を含む約450行)がある
   が、一般の`𝒪_K`には使えない。
   **⟹ `𝒪[K.carrier]`のHensel's lemmaは、それ自体がmathlibの1ファイル
   分に匹敵する独立した仕事**(`IsPrecomplete.prec`——本プロジェクトで
   `unitReductionHom_surjective`に使った道具——を使ってNewton逐次近似を
   自前で組む筋が有力)。不分岐拡大理論の中でも、これが最初の重い
   関門になる。

   ★★2026-09-05さらに実測(Hensel回避ルートの検討結果): Hensel's
   lemmaを**回避する**ルートとして、mathlibの**分岐理論の基本等式**
   `Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing`
   (`e·f=[L:K]`、`NumberTheory/RamificationInertia/Basic.lean:650`)を
   使う筋を検討した(「剰余次数`f≥n`さえ言えれば`e=1,f=n`が従う」)。
   しかしこの補題が要求するinstanceを実測したところ:
   - `IsDedekindDomain (𝒪[K.carrier])`: **OK**(`valuationRing_isDVR`
     から`infer_instance`で自動)。
   - `IsFractionRing (𝒪[K.carrier]) K.carrier`: **OK**だが自動では
     見つからない——`ValuationRing.instIsFractionRingInteger Valued.v`
     を明示的に与える必要がある(`exact?`で発見)。
   - `IsFractionRing (adjoinIntegers K x) K.carrier⟮x⟯`: **`exact?`が
     `isDefEq`でタイムアウト**(200000 heartbeats)——`adjoinIntegers`が
     `{y|‖y‖≤1}`という**ノルムでの定義**であり、`Valued.v`由来の
     `v.integer`と構文的に一致しないためと思われる。橋渡しが要る。
   - `Module.Finite (𝒪[K.carrier]) (adjoinIntegers K x)`: **`exact?`
     で見つからない**——別途証明が要る。
   - 他に`IsDedekindDomain (adjoinIntegers K x)`・各種`IsScalarTower`
     も要る(未確認)。

   ⟹ **Hensel回避ルートも「instanceの地ならし」という別種の重い作業**
   を伴う。どちらの筋を選ぶにせよ、不分岐拡大理論は節目(5)と同等以上の
   規模になることが実測から確認できた。

   ★★★2026-09-05さらに実測(instance地ならしの詳細な現状——次に
   戻る人はここから読むこと):
   - **`adjoinIntegers K x = 𝒪[K.carrier⟮x⟯]`は`rfl`で成り立つ**
     (実測済み、重要な発見)——つまり`adjoinIntegers`は拡大体に同じ
     `𝒪[...]`記法を適用したものと**定義的に同一**であり、`Valued`系の
     mathlib機構がそのまま適用できる可能性が開けた。
   - `IsFractionRing (adjoinIntegers K x) K.carrier⟮x⟯`: **OK**——
     ただし`ValuationRing.instIsFractionRingInteger (K := K.carrier⟮x⟯)
     Valued.v`のように**体を明示的に指定**する必要がある。指定しないと
     `exact?`も明示適用も`whnf`/`isDefEq`でタイムアウトする(200000
     heartbeats)が、指定すれば**0.3秒**で通る。★この「体を明示すれば
     速い」という違いは大きい——`𝒪[...]`絡みのinstanceを扱うときは
     まず体を明示すること。
   - `IsDiscreteValuationRing (adjoinIntegers K x)`: **未達**。
     `Valued.integer.isDiscreteValuationRing_of_compactSpace`
     (`valuationRing_isDVR`が基礎体で使っているのと同じ補題)に
     `compactSpace_adjoinIntegers`(既出)を与えるところまでは進むが、
     最後に**`Valued.v.IsNontrivial`(拡大体の付値の非自明性)**の
     instanceが無くて止まる。基礎体`K.carrier`では自動的に見つかって
     いる(`valuationRing_isDVR`が通っている)ので、その導出経路を
     追って拡大体版を作るのが次の一歩。
   - `Module.Finite (𝒪[K.carrier]) (adjoinIntegers K x)`: **未達**
     (`exact?`で見つからない)。
   - `adjoinPAdicLocalField K x`(既存)の`.carrier`は`K.carrier⟮x⟯`と
     `rfl`で一致するが、**`𝒪[(adjoinPAdicLocalField K x).carrier]`と
     `adjoinIntegers K x`は`rfl`では一致しない**(instance diamond、
     `tools/lean-idioms.md` #31/#33の類型)——この経路でDVR性を借りて
     くるのは避けたほうがよい。

   ⟹ 残る具体的な障害は**2つだけ**(`Valued.v.IsNontrivial`の拡大体版
   と`Module.Finite`)。これらが埋まれば`Ideal.ramificationIdx_mul_
   inertiaDeg_of_isLocalRing`(基本等式`e·f=[L:K]`)が使えるようになり、
   Hensel's lemmaを自作せずに不分岐拡大の理論を進められる見込み。

   ★★★★2026-09-05(2つの障害のうち**1つを解消**、commit済み):
   - **`isNontrivial_valued_adjoin`**(新規、sorry無し): 拡大体
     `K.carrier⟮x⟯`の付値の非自明性。`p`の像のノルムが`1`未満
     (`norm_natCast_p_lt_one`+`norm_algebraMap'`)かつ`0`でない
     (`Valuation.ne_zero_iff`+標数`0`)ことから構成した。★鍵になった
     発見: **`Valued.v y = ‖y‖₊`が`NNReal.eq rfl`で成り立つ**
     (基礎体でも拡大体でも)——`NormedField.v_eq_valuation`経由の
     `rw`は(instance経路の違いで)効かないが、この`rfl`一発で済む。
   - **`isDiscreteValuationRing_adjoinIntegers`**(新規、sorry無し):
     `adjoinIntegers K x`が離散付値環——基礎体の`valuationRing_isDVR`と
     同じ`Valued.integer.isDiscreteValuationRing_of_compactSpace`に
     `compactSpace_adjoinIntegers`(既出)と上を与えるだけ。これで
     `IsDedekindDomain (adjoinIntegers K x)`も従う。
     ★配管の罠(記録): `compactSpace_adjoinIntegers K x`は
     `CompactSpace (adjoinIntegers K x)`という形だが補題側は
     `CompactSpace 𝒪[K.carrier⟮x⟯]`を要求する——`rfl`で一致するのに
     **instance探索は`def`の壁を越えられない**(lean-idioms #23/#31)。
     `haveI : CompactSpace 𝒪[...] := compactSpace_adjoinIntegers K x`と
     **補題側の形で書いて`:=`で渡す**ことで解決(defeqは`exact`が受ける)。

   ⟹ **残る障害は`Module.Finite (𝒪[K.carrier]) (adjoinIntegers K x)`
   ただ1つ**(`exact?`では見つからない——`adjoinIntegers K x`が
   `𝒪[K.carrier]`の`K.carrier⟮x⟯`における整閉包であることを経由して
   `IsIntegralClosure.finite`系の補題を使う筋が有力だが未着手)。
   これが埋まれば分岐理論の基本等式が使える。

   ★★★★★2026-09-05(残る1つの障害の**正体を特定**): mathlibに
   **`IsIntegralClosure.finite [IsIntegrallyClosed A][IsNoetherianRing A]
   : Module.Finite A C`**(`RingTheory/DedekindDomain/IntegralClosure.
   lean:174`)がある——`A:=𝒪[K.carrier]`(DVRなので両hypothesisとも
   自動)に適用できれば`Module.Finite`が即座に出る。必要なのは
   **`IsIntegralClosure (adjoinIntegers K x) (𝒪[K.carrier])
   K.carrier⟮x⟯`**、すなわち「`adjoinIntegers K x`(ノルム≤1の元)が
   ちょうど`𝒪_K`上整な元全体である」という特徴づけ。
   - 易しい向き(整⟹ノルム≤1): 付値の乗法性+モニック方程式から。
   - **難しい向き(ノルム≤1⟹整)**: これが本丸。完備非アルキメデス体の
     標準的な事実だが、mathlibに直接の補題は見当たらない
     (`Valuation.Integers.isIntegrallyClosed_integers`は「`v.integer`
     が自分の商体で整閉」であって別物)。使えそうな材料:
     `spectralNorm.spectralNorm_eq_norm_coeff_zero_rpow`
     (`spectralNorm K L x=‖(minpoly K x).coeff 0‖^(1/deg)`、mathlib)
     ・`norm_root_le_spectralValue`・Gauss補題/Newton多角形の議論。
   ⟹ 次に戻るときは**この「ノルム≤1⟹整」1本に集中する**のがよい
   ——これさえ示せれば`IsIntegralClosure`→`Module.Finite`→分岐理論の
   基本等式、と一気に繋がる。
2. 剰余体`F_q`の次数`n`拡大`F_{q^n}`を`GaloisField p (f*n)`
   (`hq:Fintype.card(ResidueField 𝒪_K)=p^f`から`p^{f*n}=q^n`)として
   具体的に構成し、その原始元の最小多項式(次数`n`、分離的)を
   `Henselian`で`𝒪_K`へ持ち上げる、という筋を検討する。
3. 持ち上げた多項式の`K.closure`内での根`x`が実際に`residueDegree
   K x=q^n`(不分岐性の条件)を満たすことを確認する——ここで
   `Ideal.ramificationIdx`・`Algebra.isUnramifiedAt_iff_map_eq`
   (mathlib、既出)の適用を試す。

この3ステップが次数`n`の不分岐拡大1つを構成する核心。一意性・
`K^ur`全体・Galois群の構造はその後。

## ★★★ 2026-09-05 突破: 「ノルム≤1 ⟹ 整」はスペクトルノルムで数行、
`e·f=[L:K]` まで一気に到達した(すべて`sorry`無し・全ゲート通過)

上で「次に戻るときはこの1本に集中するのがよい」と書いた**難しい向き
(ノルム≤1⟹整)**が、mathlib の **spectral norm** の道具立てで
**数行**で書けた。Hensel の補題も共役の対称式も自前で組む必要が無い。

### 使った鍵(すべて mathlib、`Analysis/Normed/Unbundled/SpectralNorm.lean`)

* `NormedAlgebra.norm_eq_spectralNorm K x : ‖x‖ = spectralNorm K L x`
  ——`[NontriviallyNormedField K][IsUltrametricDist K][NormedField L]
  [NormedAlgebra K L][Algebra.IsAlgebraic K L][CompleteSpace K]`。
  本プロジェクトの`K.carrier`・`K.carrier⟮x⟯`で**インスタンスが
  すべて自動で揃う**(0.9秒で通る)。
* `spectralNorm K L y = spectralValue (minpoly K y)` は**`rfl`**。
* `spectralValue_le_one_iff (hP : P.Monic) :
  spectralValue P ≤ 1 ↔ ∀ n, ‖P.coeff n‖ ≤ 1`。
* 逆向きは `norm_root_le_spectralValue` に `f := spectralAlgNorm K L`
  (`spectralAlgNorm_isPowMul`・`isNonarchimedean_spectralNorm`が既存)。
* 係数を部分環に落とすのは `Polynomial.toSubring`・
  `Polynomial.monic_toSubring`・`Polynomial.map_toSubring`。

### この日に`Found/PGC/UnramifiedExtension.lean`へ入った定理(sorry無し)

1. `isIntegral_of_norm_le_one` (難しい向き)
2. `norm_le_one_of_isIntegral` (易しい向き)
3. `isIntegral_iff_norm_le_one` : `IsIntegral 𝒪[K.carrier] y ↔ ‖y‖ ≤ 1`
4. `isNontrivial_valued_carrier` / `isDiscreteValuationRing_carrierIntegers`
   (基礎体版。`IsNoetherianRing 𝒪[K.carrier]`がここから出る)
5. `isIntegralClosure_adjoinIntegers` :
   `IsIntegralClosure (adjoinIntegers K x) 𝒪[K.carrier] K.carrier⟮x⟯`
6. `module_finite_adjoinIntegers` :
   **`Module.Finite 𝒪[K.carrier] (adjoinIntegers K x)`**
   ——長らく最後の1点だった前提。`IsIntegralClosure.finite`に(5)を
   与えるだけ(分離性は標数0から自動)。
7. `ramificationIndex` / `inertiaDegree` (薄い`def`ラッパー)
8. **`ramificationIndex_mul_inertiaDegree : e·f = [K(x):K]`**
   ——`Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing`が
   そのまま使えるようになった。局所体の分岐理論の基本等式が
   本プロジェクトの`PAdicLocalField`設定で**利用可能**になった。
9. `IsUnramifiedAdjoin K x := (ramificationIndex K x = 1)` と
   `inertiaDegree_eq_finrank_of_isUnramified`。

### ★配管(記録、`tools/lean-idioms.md`の類型)

* `rw [← NormedAlgebra.norm_eq_spectralNorm ...]`は**発火しない**
  (`Valued`由来の`NormedField`インスタンス経路と補題側が syntactic
  に一致しない)。`exact`/`le_of_eq_of_le`にすると defeq で通る(#37)。
* `Ideal.ramificationIdx`/`inertiaDeg`は`IsLocalRing (adjoinIntegers
  K x)`を**主張の型の段階**で要求する(`haveI`を証明の中に置いても
  遅い)。`residueDegree`と同じく`haveI := isLocalRing_adjoinIntegers
  K x`を**`def`の本体**に置いた薄いラッパーを挟むと、利用側に
  インスタンス束縛を波及させずに済む。
* `IsScalarTower 𝒪[K.carrier] K.carrier K.carrier⟮x⟯` と
  `IsScalarTower 𝒪[K.carrier] (adjoinIntegers K x) K.carrier⟮x⟯` は
  **インスタンス探索では見つからないが**
  `IsScalarTower.of_algebraMap_eq (fun _ => rfl)` で即座に出る。
  一方 `Algebra 𝒪[K.carrier] K.carrier⟮x⟯` と
  `Algebra 𝒪[K.carrier] (adjoinIntegers K x)` は**自動で見つかる**。

### 次の一歩(不分岐拡大の存在)

`e·f=[L:K]`が手に入ったので、当初の見取り図の 2.→3. に進める:
剰余体`F_q`の次数`n`拡大を`GaloisField`で構成し、その原始元の
最小多項式を`𝒪_K`へ持ち上げ、根`x`に対して
`IsUnramifiedAdjoin K x`(= `ramificationIndex K x = 1`)と
`Module.finrank K.carrier K.carrier⟮x⟯ = n`を示す。
持ち上げに Hensel が要る点は変わらない(mathlibに完備局所環の
Henselian性の一般インスタンスは無いことを確認済み)が、
**`e·f=[L:K]`があるので`e=1`を`f=n`から出す**という逆向きの
経路も取れる——`f = inertiaDegree` の計算(剰余体の拡大次数)に
持ち込めれば Hensel を迂回できる可能性がある。

## 2026-09-05(続き): 剰余体側からの不分岐判定まで到達

`e·f=[L:K]` に続けて、**剰余体の側から `e=1` を判定する**道具を
`Found/PGC/UnramifiedExtension.lean` に入れた(すべて`sorry`無し、
`lake build ABC3` 6590 jobs 成功・`check.mjs` NG 13 件・文字化けなし)。

* `instIsLocalRingAdjoinIntegers`(インスタンス化)
* `isLocalHom_adjoinIntegersAlgebraMap` :
  `IsLocalHom (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))`
  ——`norm_algebraMap'` でノルムが保たれるので非単元は非単元へ。
  これで `Ideal.LiesOver` と剰余体間の `Algebra` が自動で出る。
* `nontriviallyNormedField_adjoin`(`@[implicit_reducible]` 必須)
* `finite_residueField_adjoinIntegers` : 拡大体の剰余体も有限
  ——`Found/ResidueFieldFinite.lean` の一般結果を拡大体に適用。
* `inertiaDegree_eq_finrank_residueField` : `f` は剰余体の拡大次数
* `residueDegree_eq_residueCard_pow` : `q_L = q^f`
* `isUnramifiedAdjoin_of_inertiaDegree_eq_finrank` : `f=[L:K] ⟹ e=1`
* `one_lt_residueCard` : `q ≥ 2`
* `isUnramifiedAdjoin_of_residueDegree` :
  **`q_L = q^{[L:K]}` ⟹ 不分岐**

### ★配管(記録)

* `nontriviallyNormedField_adjoin` に `@[implicit_reducible]` を
  付けないと、`letI` で入れた瞬間に
  `IsUltrametricDist ↥K.carrier⟮x⟯` の探索が失敗する
  (ノルム構造の経路が変わるため)。
* `Valuation.Integer.not_isUnit_iff_valuation_lt_one` は
  `adjoinIntegers K x` の元に対して `apply` では unify しない
  (`↥(adjoinIntegers K x)` と `↥v.integer` の `def` の壁)。
  `(v := ...)` と `(x := ...)` を**明示**して `refine` すると通る。
* `‖·‖₊`(nnnorm)版の `norm_algebraMap'` は無いので
  `NNReal.eq (norm_algebraMap' _ _)` で作る。

### 次の一歩と、そこで見えている選択

不分岐拡大の**存在**(次数`n`ごとに1つ)には、
`isUnramifiedAdjoin_of_residueDegree` の仮定
`q_L = q^{[L:K]}` を満たす `x` を作ればよい。古典的な作り方は
`x = ζ_{q^n-1}`(1の`q^n-1`乗根)だが、
`[K(ζ):K] ≤ n` を出す段で結局 Hensel の補題(`X^m-1` の
`𝒪_K` 上での分解が剰余体上の分解を持ち上げること)が要る。

★2026-09-05 実測(再確認): mathlib の `RingTheory/Henselian.lean`
には `HenselianLocalRing` の**定義と TFAE** はあるが、
`.cache/mathlib-index.txt` を "henselian" で引いた結果は
`HenselianLocalRing`(class)・`HenselianRing`(class)・
`HenselianLocalRing.TFAE`・
`IsLocalRing.eq_of_eval_eq_zero_of_not_isUnit_sub`・
`isLocalHom_of_le_jacobson_bot` の**5 件のみ**——
「完備局所環は Henselian」に相当するインスタンスは**無い**。
`ℤ_[p]` 専用の Hensel(`NumberTheory/Padics/Hensel.lean`、約450行)
は一般化しない。

したがって不分岐拡大の存在を出すには
(a) `𝒪[K.carrier]` が `HenselianLocalRing` であることを
    自分で示す(完備性+極大イデアルでの逐次近似、Newton法)、
(b) それを避けて別経路を探す、
のどちらか。(a) は `HenselianLocalRing.TFAE` の 2 番目
(「剰余体での単根は持ち上がる」)を目標にすれば、
`𝒪[K.carrier]` の完備性(既存)と極大イデアル進の収束で書けるはず
——`ℤ_[p]` 版 450 行より短くなる見込みは薄いが、
`e·f=[L:K]` と剰余体判定がすでにあるぶん、
「持ち上げたあと何を示すか」は明確になっている。

## ★★★訂正(2026-09-05、同日中): Hensel の補題は mathlib に**ある**

上で二度「完備局所環が Henselian であることに相当するものは mathlib
に無い」と記録したが、**これは誤測定だった**。`HenselianLocalRing` を
**結論とする**宣言だけを引いていた(`Field.henselian` しか無い)。
正しい入口は **`HenselianRing R I`** の方で、

```
IsAdicComplete.henselianRing (R) (I) [IsAdicComplete I R] : HenselianRing R I
```

が `Mathlib/RingTheory/Henselian.lean` にある。そして本プロジェクトは
すでに `ABC3.Found.PGC.isAdicComplete_valuationRing`
(`Found/PGC/ValuationRingComplete.lean`)を持っていた。

⟹ `HenselianLocalRing 𝒪[K.carrier]` は**6 行**で出た
(`Found/PGC/UnramifiedExtension.lean::henselianLocalRing_carrierIntegers`、
`IsUnit.map` を一つ挟むだけ)。`ℤ_[p]` 専用の 450 行は不要だった。

★測定の教訓(`Found/ResidueFieldFinite.lean` の docstring と同じ轍):
**「無い」という測定は探索範囲を書かないと再現できない**。今回は
`grep -i henselian .cache/mathlib-index.txt` の 5 件だけを見て
「インスタンスが無い」と結論したが、必要だったのは
「その 5 件のうち `HenselianRing` という**別のクラス**を結論とする
`IsAdicComplete.henselianRing` を実際に `#check` すること」だった。
索引の grep は**宣言名**にしか当たらないので、
`run_cmd` で「結論の head constant が C である宣言」を列挙する方が
確実(今回はこれで一発だった)。

## 不分岐拡大の存在——Hensel すら要らない筋が見えた

`isUnramifiedAdjoin_of_inertiaDegree_eq_finrank`(`f=[L:K] ⟹ e=1`)が
あるので、次数 `n` の不分岐拡大は次の手順で作れるはず(Hensel 不要):

1. 剰余体 `𝓀` 上の**モニック既約** `n` 次多項式 `g` を取る。
2. `Polynomial.lifts_and_degree_eq_and_monic` で `g` を
   `𝒪[K.carrier]` 上のモニック `n` 次 `f` に持ち上げる。
3. `Polynomial.Monic.irreducible_of_irreducible_map` で
   `Irreducible f`(`𝒪_K` 上)。
4. `Polynomial.Monic.irreducible_iff_irreducible_map_fraction_map`
   (`IsIntegrallyClosed 𝒪_K` は既出)で `K.carrier` 上でも既約。
5. `K.closure` の中の根 `x` を取る ⟹ `[K(x):K] = n`。
6. `norm_le_one_of_isIntegral` で `‖x‖ ≤ 1`、剰余 `x̄` は `g` の根。
7. `g` 既約モニックだから `minpoly 𝓀 x̄ = g`、よって
   `n ≤ f = inertiaDegree K x`。`e·f=n` と `e≥1` から `f=n`・`e=1`。

残る材料の所在: 1. の「有限体上の任意次数のモニック既約多項式の存在」
——mathlib での名前は未確認(`GaloisField` の原始元の `minpoly` を
使う経路が確実)。

## ★★★★2026-09-05: 不分岐拡大の**構成**が完成(Hensel を使わずに)

`Found/PGC/UnramifiedExtension.lean::exists_isUnramifiedAdjoin_of_irreducible`
(sorry 無し、全ゲート通過):

```
剰余体 𝓀[K.carrier] 上のモニック既約 g  ⟹  ∃ x : K.closure,
  [K(x):K] = deg g  ∧  IsUnramifiedAdjoin K x
```

証明の骨格(Hensel の補題は**使っていない**):
1. `Polynomial.lifts_and_degree_eq_and_monic` で `g` を `𝒪_K` 上の
   モニック `f`(同次数)へ持ち上げる。
2. `Polynomial.Monic.irreducible_of_irreducible_map` で `f` も既約。
3. Gauss(`Monic.irreducible_iff_irreducible_map_fraction_map`、
   `IsIntegrallyClosed 𝒪_K` は付値環だから自動)で `K` 上でも既約。
4. `IsAlgClosed.exists_root` で根 `x`、`minpoly.eq_of_irreducible_
   of_monic` から `[K(x):K] = deg g`。
5. `f` モニック ⟹ `x` は整 ⟹ `‖x‖≤1`
   (`norm_le_one_of_isIntegral`)⟹ `x ∈ adjoinIntegers K x`。
6. `Polynomial.hom_eval₂` を二度(`Subring.subtype` と
   `IsLocalRing.residue`)使って剰余体へ落とし、`x̄` が `g` の根。
7. `minpoly 𝓀 x̄ = g` ⟹ `deg g ≤ f = inertiaDegree K x`
   (`minpoly.natDegree_le`)。
8. `e·f = [K(x):K] = deg g` と `e ≥ 1` から `f ≤ deg g`。
   両方合わせて `f = deg g = [K(x):K]`、すなわち **`e = 1`**。

また `finiteDimensional_adjoin_closure` をインスタンス化した
(`K.closure` は代数閉包なので単項生成中間体はつねに有限次元)——
以降 `[FiniteDimensional ...]` 束縛は自動で埋まる。

### 残り 1 点: 有限体上の任意次数のモニック既約多項式の存在

これだけが「次数 `n` の不分岐拡大が各 `n` に存在する」を言うために
足りない。純粋に有限体の話で、この節とは独立。mathlib での所在は
未確認(`GaloisField p n` の原始元の `minpoly` を使う経路が有力)。

## ★★★★★2026-09-05: 各次数の不分岐拡大の**存在**が完成

最後に残っていた「有限体上の任意次数のモニック既約多項式の存在」を
新規ファイル `Found/FiniteFieldIrreducible.lean` に書いた
(`ABC3.Found.exists_monic_irreducible_natDegree_eq`、一般の結果として
mathlib へ出せる形。`Found/ResidueFieldFinite.lean` と同じ位置づけ)。

筋: `F` の標数 `p`、`m := [F : ZMod p]` として `GaloisField p (m*n)` を
取り、`m ∣ m*n` から `FiniteField.nonempty_algHom_of_finrank_dvd` で
`F`-代数とみなす。塔の公式で `[GaloisField p (m*n) : F] = n`、
有限体は完全体だから `Field.exists_primitive_element` で原始元 `α`、
`minpoly F α` が求めるもの。

これを流し込んで
`Found/PGC/UnramifiedExtension.lean::exists_isUnramifiedAdjoin`:

```
n ≠ 0  ⟹  ∃ x : K.closure, [K(x):K] = n ∧ IsUnramifiedAdjoin K x
```

**不分岐拡大理論の第一の柱(存在)が立った**。

### 次の柱

1. **一意性**: 同じ次数の不分岐拡大は(`K.closure` の中で)一致する。
   剰余体が `𝔽_{q^n}` で決まることと `e=1` から、`K(x)` の元は
   `𝒪_{K(x)}` の Teichmüller 持ち上げで決まる——ここで
   `HenselianLocalRing`(本日入った)が効くはず。
2. `K^ur` の構成(すべての不分岐拡大の合併=`IntermediateField` の
   `⨆`)と `Gal(K^ur/K) ≅ Gal(𝔽̄_q/𝔽_q) ≅ Ẑ`。
   `Gal(𝔽̄_q/𝔽_q) ≅ Ẑ` そのものは mathlib に不在(2026-09-05 実測)
   なので、Frobenius の位相的生成性から自前で組む必要がある。

## 2026-09-05: 拡大体側も Henselian(一意性への準備)

`Found/PGC/UnramifiedExtension.lean`(すべて sorry 無し):
- `isAdic_maximalIdeal_adjoinIntegers` : 拡大体の整数環でも
  「`maximalIdeal`-進位相 = 距離位相」
- `isAdicComplete_adjoinIntegers` : 同じく `maximalIdeal`-進完備
- `henselianLocalRing_adjoinIntegers`(instance):
  **`HenselianLocalRing (adjoinIntegers K x)`**

基礎体版(`ValuationRingComplete.lean`)と同じ筋がそのまま通ったが、
2 つの配管の違いがあった(`tools/lean-idioms.md` #55 に記録):
* `↥(adjoinIntegers K x)` と `↥𝒪[K.carrier⟮x⟯]` は `CommRing` の
  インスタンス**経路**が違う(`CommRing.toCommSemiring` vs
  `SubsemiringClass.toCommSemiring`)ので `rw` が型検査で落ちる。
  `▸` の項レベルキャストにすると defeq で通る。
* `Valued.integer.norm_irreducible_pos` 等は
  `[NontriviallyNormedField K]` を要求する。拡大体にはこの instance が
  無いので `letI := nontriviallyNormedField_adjoin K x` を先に置く
  (置き忘れると `whnf` で 200000 heartbeats タイムアウト)。

### 一意性の残り

`HenselianLocalRing (adjoinIntegers K y)` が揃ったので、classical な
一意性の議論(「`𝒪_{K(y)}` の剰余体 `𝔽_{q^n}` にある `ḡ` の根を
Hensel で持ち上げ、`K(x) ⊆ K(y)`、次数が等しいので一致」)の最後の
道具が入った。残る材料:
1. 「次数 `n` の既約多項式は `𝔽_{q^n}` に根を持つ」
   ——`FiniteField.nonempty_algHom_of_finrank_dvd` から出るはず。
2. `K⟮x'⟯ = K⟮x⟯`(`g` の根はどれも同じ体を生成する)——不分岐拡大が
   Galois であること(剰余体側が Galois で Hensel の持ち上げが一意)
   から出る。ここが一番手数が要りそう。

## ★★★★★★2026-09-05: 不分岐拡大は **Galois**——Hensel で分裂を持ち上げた

新規 `Found/HenselianSplits.lean`(一般の結果、mathlib へ出せる形):
- `exists_root_of_residue_root` : Henselian 局所環で、剰余体の**単根**
  (`F̄(b)=0` かつ `F̄'(b)≠0`)は `A` の根に持ち上がる。
- `splits_map_of_residue_splits` : `A` Henselian 局所整域、`F` モニックで
  `F̄` が**分離的かつ剰余体で分裂**するなら、`A` を埋め込む任意の体で
  `F` も分裂する。骨格は「`residue` が `F.roots → F̄.roots` の全射
  (Hensel)⟹ `A` の中に `deg F` 個以上の相異なる根」+
  `roots.card ≤ natDegree` の挟み撃ち。

`Found/PGC/UnramifiedExtension.lean`:
- `splits_adjoin_of_lift` : 剰余体の既約多項式の持ち上げ `f` は `K(x)` で
  分裂する(`𝓀_{K(x)}/𝓀` は有限体の有限次拡大なので normal、分離性は
  完全体から)。
- `normal_of_splits_in_adjoin` : `K(x)` で `F` が分裂し `x` が根なら
  `K(x)` は `F` の分裂体——`Normal K.carrier K(x)`。
- `exists_isUnramifiedAdjoin_of_irreducible` / `exists_isUnramifiedAdjoin`
  の結論に **`∧ Normal K.carrier K(x)`** を追加。
- `exists_isGalois_isUnramifiedAdjoin` : 標数 0 なので **Galois**。

★ファイル内の順序を入れ替えた: 拡大体側の Henselian 節
(`isAdic_...`/`isAdicComplete_...`/`henselianLocalRing_adjoinIntegers`)を
構成節の**前**へ移した(`splits_adjoin_of_lift` が Henselian を要求する
ため)。`finiteDimensional_adjoin_closure`(instance)はさらにその前。

### 次: 一意性 → `K^ur` → `Gal(K^ur/K) ≅ Ẑ`

Galois 性が付いたので、同じ次数の不分岐拡大が一致することは
「`𝒪_{K(y)}` の剰余体に `ḡ` の根がある(次数が割り切るから)⟹ Hensel で
持ち上げ ⟹ `K(x)` の共役が `K(y)` に入る ⟹ normal だから `K(x) ⊆ K(y)`」
で言えるはず。材料はすべて揃った。

## ★★★★★★★2026-09-05: `Gal(K(x)/K) ≃* Gal(𝓀_{K(x)}/𝓀)`——不分岐理論の核心

`Found/PGC/UnramifiedExtension.lean`(すべて sorry 無し):
- `norm_algEquiv` : `K(x)` の `K`-自己同型はノルムを保つ
  (`spectralNorm_eq_of_equiv` + `NormedAlgebra.norm_eq_spectralNorm`)
- `algEquivIntegers` : したがって整数環の環同型を誘導
- `residueAlgEquiv` / `residueAlgEquiv_apply` : 剰余体の `𝓀`-代数同型
- `residueGalHom` : **還元射 `Gal(K(x)/K) →* Gal(𝓀_{K(x)}/𝓀)`**
- `algEquiv_eq_one_of_fixes_gen` : 生成元を固定する自己同型は恒等
  (`IntermediateField.algHom_ext_of_eq_adjoin`)
- `isUnit_eval_derivative` : `f̄` の根で `f'` は単元(分離性)
- `residueGalHom_injective` : **単射**——Hensel の一意性
  (`IsLocalRing.eq_of_eval_eq_zero_of_not_isUnit_sub`)。`σ` が剰余体で
  恒等なら `σ(x)` と `x` は同じ剰余を持つ `f` の 2 根、`f'(x)` が単元
  なので `σ(x)=x`、`x` は生成元だから `σ=1`。
- `residueGalHom_bijective` : **全単射**——位数の比較。
  `|Gal(K(x)/K)| = [K(x):K]`(Galois)、`|Gal(𝓀_{K(x)}/𝓀)| = f`
  (有限体は Galois)、不分岐なので `[K(x):K] = f`。
- `exists_isUnramifiedAdjoin_of_irreducible` / `exists_isUnramifiedAdjoin`
  の結論に `∧ Function.Bijective (residueGalHom K x)` を追加。
- **`exists_mulEquiv_residueGal`** : 各次数 `n` に対し
  `Gal(K(x)/K) ≃* Gal(𝓀_{K(x)}/𝓀)`。

右辺は有限体の Galois 群なので Frobenius が生成する巡回群 `ℤ/n`
——**これが `Gal(K^ur/K) ≅ Ẑ` の出どころ**。

★ファイル内の順序をもう一度入れ替えた: 還元射の節を構成節の前へ。

### 次

1. `Gal(𝓀_{K(x)}/𝓀) ≃ ℤ/n`(Frobenius が生成)——mathlib に有限体の
   Galois 群の巡回性はあるはず(`FiniteField`/`GaloisField` 周辺)。
2. 一意性(同じ次数の不分岐拡大は `K.closure` の中で一致)。
3. `K^ur` の構成と `Gal(K^ur/K) ≅ Ẑ`(射影極限)。

## 2026-09-05: 不分岐拡大の Galois 群は `ℤ/n`

`Found/PGC/UnramifiedExtension.lean`:
- `exists_isCyclic_gal` : 次数 `n` の不分岐拡大の Galois 群は
  **位数 `n` の巡回群**(mathlib のインスタンス
  `IsCyclic (E ≃ₐ[F] E)`(有限体)を `residueGalHom` の同型で移す)
- `exists_gal_mulEquiv_zmod` : したがって
  **`Gal(K(x)/K) ≃* Multiplicative (ZMod n)`**
  (`zmodCyclicMulEquiv`)

`Gal(K^ur/K) ≅ Ẑ` は、これらを `n` について射影極限に組み上げたもの
——残る仕事は (a) 一意性(次数 `n` の不分岐拡大が `K.closure` の中で
一意)と (b) 極限の構成。

## ★★★★★★★★2026-09-05: 不分岐拡大の**一意性**——第二の柱が立った

### 構造の整理(リファクタ)

中心を `isUnramifiedAdjoin_of_lift`(**不分岐拡大の判定子**)に一本化した:

```
f : 𝒪[K.carrier][X] モニック、f̄ 既約、x は f_K の根
 ⟹ [K(x):K] = deg f ∧ IsUnramifiedAdjoin K x ∧ Normal ∧ Bijective (residueGalHom K x)
```

存在(`exists_isUnramifiedAdjoin_of_irreducible`)はこれに流し込むだけの
30 行に縮んだ。

### 逆向き——不分岐性から定義多項式を取り出す

`exists_integral_generator` : `IsUnramifiedAdjoin K x` なら、剰余体の
原始元(`Field.exists_primitive_element`)を Hensel で持ち上げた `θ` が
`K(x)` の生成元になり、その最小多項式は剰余体で既約な `f` の持ち上げ。
——`[K(θ):K] = deg f`(既約性から直接)と `K(θ) ⊆ K(x)` の次数一致で
`K(θ) = K(x)`。これで**任意の**不分岐拡大に判定子を適用できる。

### 結果

- `normal_of_isUnramifiedAdjoin` : **不分岐 ⟹ normal**(標数 0 で Galois)
- `mem_adjoin_of_root_of_splits` : `F` が `K(x)` で分裂するなら
  `K.closure` の `F` の根はすべて `K(x)` に入る
- `exists_root_lift_of_residue_root` : 剰余体の `f̄` の根は Hensel で
  `𝒪_{K(x)}` の `f` の根に持ち上がる
- `ABC3.Found.exists_root_of_finrank_eq`(`Found/FiniteFieldIrreducible.lean`)
  : 次数が一致する有限体拡大には既約多項式の根がある(`AdjoinRoot` +
  `FiniteField.nonempty_algHom_of_finrank_dvd`)
- **`adjoin_eq_of_isUnramified`** : **同じ次数の不分岐拡大は
  `K.closure` の中で一致する**(一意性)
- **`exists_unique_adjoin_isUnramified`** : 各次数 `n ≥ 1` に不分岐拡大が
  **存在して一意**

一意性の筋: `K(y)` の定義多項式 `f` を取る。`𝓀_{K(x)}` も `n` 次拡大
だから `f̄` の根 `b` を含む。`b` を Hensel で持ち上げて `f` の根
`z ∈ K(x)` を得ると `[K(z):K] = n` で `K(z) = K(x)`。一方 `K(y)` は
normal で `f_K` が分裂するから `z ∈ K(y)`、同様に `K(z) = K(y)`。

### 残り: `K^ur` と `Gal(K^ur/K) ≅ Ẑ`

存在・一意性・`Gal ≃ ℤ/n` が揃ったので、あとは `n` について
`K_n`(次数 `n` の不分岐拡大)を選び、`K^ur := ⨆ n, K_n` を作って
`Gal(K^ur/K) ≅ lim ℤ/n! = Ẑ` を出す段。

## ★★★★★★★★★2026-09-05: 最大不分岐拡大 `K^ur` を構成

`Found/PGC/UnramifiedExtension.lean`(sorry 無し):
- `adjoin_le_of_dvd` : **次数が割り切れば包含する**
  (`[K(x):K] ∣ [K(y):K] ⟹ K(x) ⊆ K(y)`)。`f̄` の根が `𝓀_{K(y)}` にある
  こと(`exists_root_of_natDegree_dvd`)を Hensel で持ち上げ、一意性で
  `K(x)` と同定する。
- `exists_isUnramified_ge` : 不分岐拡大は**有向系**をなす
  (次数 `m*n` の不分岐拡大が両方を含む)
- `unramifiedClosure`(**`K^ur`**): 不分岐な単項拡大すべての上限
- `nonempty_isUnramifiedAdjoin` / `directed_isUnramifiedAdjoin`
- **`mem_unramifiedClosure_iff`** : `z ∈ K^ur ⟺ ∃ 不分岐 x, z ∈ K(x)`
  (`IntermediateField.coe_iSup_of_directed`)
- `adjoin_le_unramifiedClosure`

`Found/FiniteFieldIrreducible.lean`:
- `exists_root_of_natDegree_dvd`(`exists_root_of_finrank_eq` の可除版)

### 残り: `Gal(K^ur/K) ≅ Ẑ`

`Gal(K_n/K) ≅ ℤ/n`(`exists_gal_mulEquiv_zmod`)を `n` について射影極限に
組み上げる段。mathlib に `Ẑ` があるか(`ZHat`?)は未確認——無ければ
`lim_n ℤ/n!` を自前で組むか、`Gal(K^ur/K)` の位相群としての構造
(profinite、`InfiniteGalois` 周辺)経由で述べる。

## ★★★★★★★★★★2026-09-05: `K^ur/K` は Galois、そして **`Γ_K ↠ ℤ/n`**

`Found/PGC/UnramifiedExtension.lean`(sorry 無し):
- `normal_unramifiedClosure` : `K^ur/K` は normal
  (`IntermediateField.normal_iSup`——normal な単項拡大の上限)
- `isGalois_unramifiedClosure` : **`K^ur/K` は Galois**
- **`exists_surjective_absGal_to_zmod`** : 任意の `n ≥ 1` に対し
  `Γ_K = Gal(K̄/K)` から `Multiplicative (ZMod n)` への**全射**準同型が
  存在する。制限射 `Γ_K →* Gal(K_n/K)` が全射
  (`AlgEquiv.restrictNormalHom_surjective`、`K_n` は normal)で、
  `Gal(K_n/K) ≃ ℤ/n` だから。

★2026-09-05 実測: mathlib に `Ẑ`(`ZHat`)は**無い**
(`.cache/mathlib-index.txt` を "ZHat"/"Zhat" で引いて 0 件)。
そこで `Gal(K^ur/K) ≅ Ẑ` を「各 `n` で `ℤ/n` へ全射」という形に
言い換えて述べた。射影極限としての `Ẑ` を作るなら
`lim_n ℤ/n!` を自前で組む必要がある。

### この節の到達点(まとめ)

古典的局所類体論の**不分岐半分**は、次の形で本プロジェクトに入った:
1. `e·f=[L:K]`(分岐理論の基本等式)
2. 不分岐拡大の**存在**(各次数)・**一意性**
3. `Gal(K_n/K) ≃ Gal(𝓀_n/𝓀) ≃ ℤ/n`
4. `K^ur` の構成と Galois 性
5. `Γ_K ↠ ℤ/n`(各 `n`)

残るのは「完全分岐部分(節目(5)、既に完成)と不分岐部分を**組み合わせて**
`Γ_K^ab ≅ (K^×)^∧` を出す」段——`pgc-prop12-reciprocity-gap.md` の
(2)(3)(4) にあたる。
