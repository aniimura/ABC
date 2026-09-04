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
