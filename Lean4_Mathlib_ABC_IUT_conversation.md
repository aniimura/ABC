# Lean 4 / Mathlib拡張とABC・IUT形式化に関する会話履歴

## 1. ユーザー質問

Lean4のMathlibをさらに拡張したものは公開されていますか？

## 2. アシスタント回答

はい。しかも「Mathlibをさらに拡張したライブラリ」という意味では、かなり多くのプロジェクトが公開されています。

ただし、**Mathlibそのものを包含して「Mathlibの上位互換」として一つに統合された巨大ライブラリ**がある、という状況ではありません。現在は、Mathlibを基盤として分野ごとのライブラリを積み上げる構造が主流です。

### 代表的なもの

| プロジェクト | 主な対象 | Mathlibとの関係 |
|---|---|---|
| **Mathlib** | 数学全般 | 基盤 |
| **Batteries** | データ構造・プログラミング・基礎機能 | Lean 4自体を拡張 |
| **AddCombi** | 加法的組合せ論 | Mathlibのサブライブラリ |
| **FLT** | フェルマーの最終定理 | Mathlib上で大規模形式化 |
| **Formal Conjectures** | 未解決問題・予想 | Mathlib上で命題を形式化 |
| **PhysLib** | 物理学 | Mathlibを利用 |
| **CSLib** | コンピュータサイエンス | Mathlibを利用 |

Lean Community自身が、Mathlibを利用するこうしたプロジェクトを一覧化しています。現在もかなり多数のプロジェクトが登録されています。

また、Mathlibに依存する公開パッケージはReservoir上で500以上確認できます。

### 「Mathlibより巨大な数学ライブラリ」という観点

現在のLeanエコシステムは、

```text
                 ┌─ AddCombi
                 ├─ FLT
Lean ── Mathlib ─┼─ PhysLib
                 ├─ Formal Conjectures
                 ├─ CSLib
                 └─ その他多数
```

という**「Mathlibを共通基盤として、その上に専門分野を積み上げる」**形に近いです。

Mathlib自体も非常に巨大で、数学だけでなく、証明自動化やプログラミング基盤まで含んでいます。

---

そして、質問の意図がABC予想やIUTの形式化との関係であれば、

**Mathlib → 数論 → 代数幾何 → 算術幾何 → IUTに必要な理論**

という巨大な依存関係を構築することになります。

Mathlibには既に代数幾何、数論、圏論、ホモロジーなどのかなり深い部分が入っています。

---

# 3. ユーザー質問

提案の階層全て

## 4. アシスタント回答

はい。では、前回提案した **「Mathlib → 代数幾何 → 数論 → 楕円曲線 → モジュラー形式 → ABC → IUT」**を、Lean 4の状況に合わせて、**「何が既に形式化されているか」「どこが未整備か」「ABC/IUTまで行くには何が足りないか」**という観点で整理します。

結論を先に言うと、

> **Mathlibは既にかなり深いところまで来ています。**
> 特に代数幾何・圏論・エタール理論などは相当進んでいます。
> 一方、**ABC → IUT**に近づくにつれて、Mathlibの既存理論だけでは足りない巨大な専門領域が急激に増えます。

そして現在は、その「足りない部分」を埋めようとするプロジェクトが複数動いています。

# 全体像

```text
Lean 4
 │
 └── Mathlib
      │
      ├── 基礎数学
      │    ├── 集合・論理
      │    ├── 群・環・体
      │    ├── 線形代数
      │    └── 位相・解析
      │
      ├── ① 代数幾何
      │    ├── Scheme
      │    ├── Spec
      │    ├── 射
      │    ├── 局所化
      │    ├── étale
      │    └── étale site
      │
      ├── ② 数論
      │    ├── 整数論
      │    ├── 素数
      │    ├── p進的構造
      │    ├── 数体
      │    ├── Dedekind domain
      │    └── Galois theory
      │
      ├── ③ 楕円曲線
      │    ├── Weierstrass equation
      │    ├── Mordell curves
      │    ├── rational points
      │    └── arithmetic geometry
      │
      ├── ④ モジュラー形式
      │    ├── modular forms
      │    ├── Hecke algebra
      │    ├── Galois representations
      │    └── modularity
      │
      ├── ⑤ ABC
      │    ├── Mason–Stothers
      │    ├── function-field ABC
      │    ├── number-field ABC
      │    └── arithmetic ABC
      │
      └── ⑥ IUT
           ├── Frobenioid
           ├── Hodge theaters
           ├── Θ-link
           ├── log-shell
           ├── theta data
           └── ABC implication
```

重要なのは、**上に行くほど「Mathlibに既にあるもの」から「専門プロジェクトで開発中のもの」へ変わっていく**ことです。

---

# ① 代数幾何

Mathlibには既に、

- Scheme
- Affine Scheme
- Spec
- morphism
- sheaf
- stalk
- residue field
- fiber
- subscheme
- smooth morphism
- flat morphism
- étale morphism

などが形式化されています。

例えば現在のMathlibには、étale morphismについて、

\[
\text{étale}
\iff
\text{flat}
\land
\text{formally unramified}
\land
\text{locally of finite presentation}
\]

といった数学的結果まで定理として存在します。

さらに、**étale site**まで来ています。

Mathlibにはsmall affine étale siteが実装され、

- Grothendieck topology
- sheaf
- dense subsite
- hypercover

なども扱われています。

つまり、

```text
環
 ↓
Spec
 ↓
Scheme
 ↓
morphism
 ↓
étale morphism
 ↓
étale site
 ↓
sheaf
```

という、かなり本格的な代数幾何の言語をLean上で扱えるところまで来ています。

### ただし重要な空白

IUTや現代的算術幾何で重要になる、

- étale cohomology
- crystalline cohomology
- p-adic Hodge theory
- arithmetic fundamental group
- deep intersection theory

などは、**「Schemeがある」＝「代数幾何が完成している」ではありません。**

ここにはまだ巨大な距離があります。

---

# ② 数論

数論についてはMathlibはかなり成熟しています。

例えば、

- 素数
- divisibility
- congruence
- quadratic reciprocity
- Dirichlet characters
- algebraic number theory
- number fields
- ideals
- Dedekind domains
- valuations
- Galois theory

などが整備されています。

したがって、

```text
整数論
  ↓
代数的整数論
  ↓
数体
  ↓
イデアル
  ↓
Dedekind domain
  ↓
Galois theory
```

というところまではかなり強いです。

---

# ③ 楕円曲線

ここから一気に重要になります。

ABC、Fermat、IUTを考える場合、

\[
E:y^2=x^3+ax+b
\]

のような楕円曲線が中心的な役割を持ちます。

Leanでは楕円曲線そのものに関する形式化も進んでいます。

さらに2026年には、**LeanBridge**というプロジェクトが非常に興味深いです。

これは、

> LMFDB（L-functions and Modular Forms Database）の数学的概念をMathlibに形式化する

ことを目指しています。

対象として、

- Number fields
- Elliptic curves
- Modular forms

を掲げています。

つまり、

```text
Mathlib
   ↓
Number fields
   ↓
Elliptic curves
   ↓
Modular forms
```

という、今回の話に直結する拡張が進んでいます。

---

# ④ モジュラー形式

ここが**Wilesの証明との接続点**です。

Fermatの最終定理は、概略すると、

\[
\text{elliptic curve}
\longrightarrow
\text{modular form}
\]

という巨大な橋を使います。

より具体的には、

```text
Fermat
   ↓
Frey curve
   ↓
Galois representation
   ↓
Modularity
   ↓
Modular form
   ↓
矛盾
```

という構造です。

ここについて重要なのが、Imperial College Londonの**FLT formalization project**です。

これは単にFLTの命題を書くプロジェクトではありません。

Wiles/Taylor-Wiles系の証明をLeanで形式化しようとしており、特に、

\[
R \cong T
\]

という**変形環とHecke代数の同型**を中心とした現代的な証明ルートを採用しています。

プロジェクト自身が、

> Mathlibから下から構築し、FLT側から上から構築して途中で合流する

という方針を明示しています。

これは、

```text
                 FLT
                  ↑
          modularity lifting
                  ↑
             R ≅ T
                  ↑
          Hecke algebra
                  ↑
          modular forms
                  ↑
        elliptic curves
                  ↑
             Mathlib
```

という巨大な塔をLeanで作ろうとしていることを意味します。

---

# ⑤ ABC予想

ここからさらに重要になります。

ABCは、

\[
a+b=c,\qquad \gcd(a,b)=1
\]

に対して、

\[
c \ll_\epsilon \operatorname{rad}(abc)^{1+\epsilon}
\]

という形の不等式を主張します。

この形式化には、

- 整数論
- 素因数
- radical
- valuations
- 高さ
- Diophantine equations
- function fields
- elliptic curves

などが必要になります。

## 既に形式化されている重要なABC関連結果

特に面白いのが**Mason–Stothers theorem**です。

これは「多項式版ABC」と考えられる定理です。

研究では、Mason–Stothers theoremとその系をLean 4で形式化しています。

さらに、

- Fermat–Cartan type equations
- elliptic curveの非パラメトリゼーション
- Davenport theorem

なども形式化されています。

つまり、

```text
ABC
 │
 ├── function field analogue
 │        ↓
 │   Mason–Stothers
 │        ↓
 │   polynomial FLT
 │
 └── number field / integer version
          ↓
        本来のABC
```

という「ABCに向かう下位理論」の形式化が既に始まっています。

---

# Formal Conjectures

ここで別の重要プロジェクトがあります。

Google DeepMindの**Formal Conjectures**です。

これはMathlibを基盤にして、数学の未解決予想をLeanで「命題として正確に記述」するプロジェクトです。

重要なのは、

> **「証明を形式化する」のではなく、「予想そのものを機械が理解できる形にする」**

という段階が大規模に整備されていることです。

これは将来的に、

```text
数学論文
   ↓
自然言語
   ↓
Lean形式化
   ↓
自動証明
```

という流れを作るための基盤にもなります。

---

# ⑥ IUT

そして最上段がIUTです。

ここは明確に区別する必要があります。

**IUT全体がLeanで証明されたわけではありません。**

しかし、2026年にはかなり具体的な形式化プロジェクトが存在します。

特に **LANA AgentsのIUT project** が重要です。

リポジトリには、

- IUT IV Section 1
- Corollary 3.12 variant
- log-shell
- theta data
- processions
- tensor-packets
- large volume container
- log-volume
- holomorphic hull

などが入っています。

さらに重要なのが、このプロジェクト自身が、

> **It does not verify IUT.**

としていることです。

つまり、

```text
IUTの数学的主張
        ↓
Leanに形式化
        ↓
定義・構造・依存関係を明示
        ↓
どこまで既存理論で証明できるか
        ↓
未形式化部分を特定
```

という段階です。

---

# IUTプロジェクトの「正直な境界」

特に重要なのが、

> Mathlibに存在しない数学的インフラを、勝手にaxiomとして導入して証明したことにしない

という方針です。

IUT projectでは、IUT I–IIIから持ち込まれる主張やMathlibに存在しないインフラについて、**explicit interfaces / certificates**として境界を作っています。

これは形式化では非常に重要です。

例えば、

```lean
axiom Mochizuki_theorem : ...
```

と書けば「証明」は一瞬で完成します。

しかし、それでは何も検証していません。

そこで、

```text
┌────────────────────────────┐
│ Lean kernelで検証できる部分 │
├────────────────────────────┤
│ interface / certificate     │
├────────────────────────────┤
│ 外部から仮定した数学        │
└────────────────────────────┘
```

という境界を明示しているわけです。

---

# 現在の全体像を「成熟度」で見る

| 階層 | 分野 | 状況 |
|---|---|---|
| 0 | 基礎数学 | **非常に成熟** |
| 1 | 代数・解析・位相 | **非常に成熟** |
| 2 | 代数幾何 | **かなり進展** |
| 3 | 代数的整数論 | **かなり進展** |
| 4 | 楕円曲線 | **発展中** |
| 5 | モジュラー形式 | **大規模形式化中** |
| 6 | Wiles型FLT | **大規模プロジェクト進行中** |
| 7 | ABC関連 | **部分的に形式化** |
| 8 | IUT | **研究段階** |
| 9 | IUT → ABCの完全証明 | **未達成** |

特に5～8の間で、大きな進歩が起きています。

---

# 「ABC」と「FLT」は少し違う

FLTは、

\[
\boxed{\text{巨大だが、既に確立された証明を形式化する問題}}
\]

です。

一方ABCは、

\[
\boxed{\text{未解決予想そのものを証明する問題}}
\]

です。

したがって、

```text
FLT
 ↓
既存のWiles/Taylor-Wiles証明
 ↓
Lean formalization
```

なのに対して、

```text
ABC
 ↓
証明そのものが存在しない
 ↓
IUT等の証明戦略
 ↓
その数学をLeanに形式化
 ↓
さらに証明を完成
```

となります。

だからABCのほうが圧倒的に難しいです。

---

# 横方向の発展：LMFDB → LeanBridge

縦の階層とは別に、

**LMFDB → LeanBridge**

という方向があります。

```text
                    LMFDB
                      │
                      ↓
                  LeanBridge
                      │
       ┌──────────────┼──────────────┐
       ↓              ↓              ↓
  Number fields   Elliptic curves  Modular forms
       │              │              │
       └──────────────┼──────────────┘
                      ↓
                   Mathlib
```

LMFDBは現代数論の巨大なデータベースなので、これをLeanの数学的構造と接続できれば、

**「論文に書かれた個別の定理」だけでなく、「現代数論のデータベースそのものを形式数学の対象にする」**

方向へ進めます。

---

# ABC/IUTについての2026年時点の状況

現在の状況を一本の線にすると、

```text
                       ┌───────────────┐
                       │    IUT        │
                       │   研究段階    │
                       └───────┬───────┘
                               │
                         ABC target
                               │
                 ┌─────────────┴─────────────┐
                 │                           │
           IUT formalization           Classical ABC
                 │                           │
                 │                    Mason–Stothers
                 │                           │
                 └─────────────┬─────────────┘
                               │
                         Arithmetic
                         Geometry
                               │
              ┌────────────────┼────────────────┐
              │                │                │
          Elliptic        Number Fields    Modular Forms
              │                │                │
              └────────────────┼────────────────┘
                               │
                         Algebraic Geometry
                               │
                            Mathlib
                               │
                            Lean 4
```

という状態です。

特に2026年時点のLANAのIUTリポジトリには、**ABC target statementそのもの**が含まれています。

したがって、

> 「LeanでABC予想を形式化する人がいる」

という段階を超えて、

> **「IUTからABCへ至る部分を、Lean上で具体的な数学的構造として切り分け始めている」**

という段階に入っています。

ただし、これは**ABCの証明が完成したという意味ではありません**。

---

# 5. 全体を俯瞰した結論

現在のLean数学では、Mathlibそのものを巨大化させることだけが重要なのではありません。

むしろ、

> **Mathlibを「共通数学基盤」として、その上にFLT、LeanBridge、Formal Conjectures、IUTなどが専門領域ごとに積み上がっている**

という構造になっています。

そしてABC/IUTは、ちょうどその**最前線**に位置しています。

2026年時点のLANAのIUTリポジトリはLean 4.32.0 + Mathlib 4.32.0に固定されており、具体的な開発段階にあります。

今後さらに重要なのは、

**「この階層のうち、ABCの完全形式化に本当に不足している数学は何なのか」**

を、

1. Mathlibに既にある
2. 少し拡張すれば使える
3. 専門プロジェクトで開発中
4. まだ巨大な未形式化領域

に分けて棚卸しすることです。

それによって、**IUT → ABCの形式化が実際にどれくらい遠いのか**をかなり具体的に測ることができます。
