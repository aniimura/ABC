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

★2026-09-05続報: 上記2で要る「`K.closure`上のspectralNormが
`K.carrier`上の元では元のノルムに一致する」という橋渡しは、
**既に`norm_algebraMap' K.closure k : ‖algebraMap K.carrier K.closure
k‖ = ‖k‖`としてmathlibに存在する**(`exact?`で即発見、
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
