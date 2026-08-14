はい、**可能です**。Lean 4では、ある `.lean` ファイルで定義した `def`、`theorem`、`lemma`、`structure`、`inductive` などを、別の `.lean` ファイルから `import` して利用できます。

Lean 4のモジュールシステムでは、基本的に**1つのLeanファイルが1つのモジュール**になり、`import`するとそのモジュールの宣言が利用可能になります。さらに、そのファイルがimportしているモジュールも推移的に利用できます。([Lean Language][1])

例えば、

```text
MyProject/
├── Basic.lean
└── Main.lean
```

### `Basic.lean`

```lean
namespace MyMath

def double (n : Nat) : Nat :=
  n * 2

theorem double_zero : double 0 = 0 := by
  rfl

end MyMath
```

### `Main.lean`

```lean
import MyProject.Basic

#check MyMath.double
#check MyMath.double_zero

example : MyMath.double 5 = 10 := by
  rfl
```

これで `Main.lean`から`Basic.lean`の定義・定理を利用できます。

### Namespaceも重要です

例えば、

```lean
namespace MyMath

def double (n : Nat) := n * 2

end MyMath
```

なら、別ファイルでは

```lean
import MyProject.Basic

#check MyMath.double
```

となります。

あるいは、

```lean
open MyMath

#check double
```

のように`open`して短い名前で使うこともできます。

---

## 数学の形式化では、この仕組みがかなり重要です

例えば大きな定理を形式化するとき、

```text
Definitions.lean
    ↓
BasicLemmas.lean
    ↓
MainTheorem.lean
```

のように分割できます。

```lean
-- Definitions.lean
def MyObject := ...
def myProperty := ...
```

↓

```lean
-- BasicLemmas.lean
import MyProject.Definitions

theorem lemma1 : ... := by
  ...
```

↓

```lean
-- MainTheorem.lean
import MyProject.BasicLemmas

theorem main_theorem : ... := by
  ...
```

この構造にすると、**定義・基本補題・大きな定理を分離**できます。実際、Lean/mathlib系のプロジェクトでも、`Defs`のような定義用モジュールと、それを利用する定理用モジュールに分ける構成が一般的です。([Reservoir][2])

また、Leanはimportされたモジュールをビルド時に管理し、コンパイル済みのモジュール情報（`.olean`等）を利用する仕組みになっています。([Lean Language][3])

**ただし注意点が一つあります。**

`A.lean`から`B.lean`をimportし、`B.lean`から`A.lean`をimportするような**循環importはできません**。

したがって、大規模な形式化では、

```text
基礎定義
   ↓
補題
   ↓
理論
   ↓
最終定理
```

という**依存関係が一方向になるようにファイルを設計する**のが重要です。

もし現在、**ABC予想などの大きな数学の形式化をLean 4で分割して作ること**を考えているのであれば、「どの単位で `.lean` ファイルを分割するとよいか」まで含めて、具体的なプロジェクト構成例を示せます。

[1]: https://lean-lang.org/papers/thesis-sebastian.pdf?utm_source=chatgpt.com "An Extensible Theorem Proving"
[2]: https://reservoir.lean-lang.org/%40SamuelSchlesinger/complexitylib?utm_source=chatgpt.com "complexitylib | Reservoir"
[3]: https://lean-lang.org/doc/api/Lake/Config/Module.html?utm_source=chatgpt.com "Lake.Config.Module"


単純にはそうなりません。むしろ、**「論文を上から順番にLeanへ翻訳する」より、「論文の数学的依存関係を分解して、Leanが扱える形に再構成する」**という作業に近いです。

特に数学論文の形式化では、次の違いが大きいです。

### 1. 論文の順番 ≠ Leanの依存関係

論文では、

> Definition → Lemma → Theorem

と書かれていても、あるLemmaの証明に、論文の後半で証明される結果と同値な既知定理が使えることがあります。

Leanでは**依存関係が循環しないように**しなければならないので、

```text
論文
  Definition A
      ↓
  Lemma B
      ↓
  Lemma C
      ↓
  Main Theorem
```

ではなく、

```text
基本定義
 ├── 基本補題1
 ├── 基本補題2
 │      ↓
 └── 補題3
        ↓
     中間定理
        ↓
     Main Theorem
```

という**依存グラフ**を作ることになります。

---

### 2. 論文には「暗黙の数学」が大量にあります

人間の数学論文では、

> 明らかに〜である

> よく知られた事実から

> 適切な変形を行えば

> この議論を繰り返せば

といった記述が非常に多いです。

Leanでは、これらを具体的な定理・補題として用意する必要があります。

例えば論文に

> (x>0) なので (x^2>0) は明らか。

と書かれていても、Leanでは

```lean
have h : 0 < x ^ 2 := by
  positivity
```

のように、**機械が検証できる証明に落とす**必要があります。

つまり形式化では、

> 「著者が省略した部分を発見する」

という作業が非常に大きな割合を占めます。

---

### 3. さらに「定義そのもの」をLean向けに作り直すことがあります

例えば論文で

> (f) は連続関数である

と書かれていても、Leanでは「連続」をどういう型・定義として扱うかが必要です。

mathlibに既存の概念があれば、

```lean
Continuous f
```

を利用できます。

一方、論文独自の概念なら、

```lean
def MyConcept (f : α → β) : Prop :=
  ...
```

のように定義します。

したがって形式化の最初には、

**「この論文で使われている概念はmathlibに既に存在するか？」**

という調査が非常に重要です。

---

## 4. 実際には「論文 → Lean」ではなく「数学 → Lean → 論文」の往復になる

実務的には、こんな流れになります。

```text
        数学論文
           │
           ▼
    数学的依存関係を解析
           │
           ▼
    既存mathlibを調査
           │
           ▼
      定義を形式化
           │
           ▼
      基本補題を証明
           │
           ▼
      中間定理を証明
           │
           ▼
      最終定理を証明
           │
           └────→ 問題発見
                     │
                     ▼
                定義を修正
```

つまり**一方向の翻訳ではなく、かなり反復的な作業**になります。

---

## 5. そして興味深いのが「論文では見えなかった穴」が見つかること

これは形式化の非常に重要な価値です。

例えば論文では、

> By Lemma 3, we immediately obtain Theorem 4.

となっていたとします。

Leanにすると、

```lean
exact lemma3
```

では型が合わない。

すると、

> Lemma 3からTheorem 4を導くには、実はもう一つ仮定が必要なのでは？

ということが分かります。

つまりLeanは単なる**証明の翻訳機**ではなく、

> **数学の論理構造を露出させる装置**

として機能します。

---

## 6. 大規模な数学ほど、この傾向が強くなります

例えば以前話していた**ABC予想のような巨大な形式化**を考えると、論文のページを1ページずつLeanにするという方法ではかなり厳しいです。

むしろ、

```text
数論
│
├── 整数論の基本定理
│
├── 素因数分解
│
├── valuation
│
├── 高さ
│
├── logarithmic height
│
├── Diophantine approximation
│
├── elliptic curves
│
├── ...
│
└── ABC theorem
```

という**数学的依存関係のDAG（有向非巡回グラフ）**を構築することになります。

そして各ノードをLeanのmodule / theorem群として実装します。

これは以前ご質問されていた**AIを使った数学研究・形式化**ともかなり相性がいい部分です。

AIに、

> 「この論文を上からLeanに変換して」

と依頼するより、

> 「この定理の証明に必要な数学的依存関係を抽出して、Leanで利用可能な既存mathlibの定理との対応を作って」

と依頼する方が、はるかに実用的です。

### まとめると

**論文の順番に形式化することはできますが、それが最適な方法ではありません。**

むしろ、

> **論文を読む → 数学的依存関係を抽出 → 既存ライブラリを探す → 定義・補題を構築 → Leanで証明 → 不足していた仮定や補題を発見 → 構造を修正**

というプロセスになります。

そしてこの観点から見ると、**「数学論文をLean化するAIエージェント」では、コード生成能力よりも「数学的依存関係の抽出」と「既存定理とのマッチング」の能力が非常に重要**になります。これは、以前お話しされていたAIを数学研究に使うというテーマにも直結するところです。


はい。**大規模な数学論文をLean形式化するなら、「大きな1ファイル」よりも、主張・概念・依存関係ごとに細かくモジュール化する方が基本的に有利**です。

ただし、**「1主張 = 1ファイル」まで細かくする必要はありません**。おすすめは「数学的な依存関係」を単位にすることです。

### 例えば

ABC予想のような形式化なら、

```text
ABC/
├── Basic/
│   ├── Definitions.lean
│   ├── Integers.lean
│   └── Valuation.lean
│
├── Height/
│   ├── Definitions.lean
│   ├── Basic.lean
│   └── Properties.lean
│
├── EllipticCurve/
│   ├── Definitions.lean
│   ├── Basic.lean
│   └── Properties.lean
│
├── Lemmas/
│   ├── Lemma01.lean
│   ├── Lemma02.lean
│   └── Lemma03.lean
│
└── ABC.lean
```

のような構造が考えられます。

重要なのは、例えば

```lean
-- Height/Definitions.lean
def height ...
```

を定義したら、

```lean
-- Height/Properties.lean
import ABC.Height.Definitions

theorem height_nonneg : ...
```

のように**依存関係を明示する**ことです。

---

## 「主張単位」はかなり良い発想です

特にAIを使って形式化する場合は有効です。

例えば論文の一節を、

```text
Claim 1
Claim 2
Claim 3
```

と分解し、それぞれについて

```text
Claim
 ↓
必要な定義
 ↓
必要な補題
 ↓
Lean proof
```

を独立して管理します。

そして最終的に、

```lean
import ABC.Lemma01
import ABC.Lemma02
import ABC.Lemma03

theorem main_theorem : ... := by
  ...
```

と統合する。

こうすると、**途中で1つの証明方法を変更しても、他の部分への影響を限定できます。**

---

## ただし「細かければ細かいほど良い」わけではありません

例えば、

```text
Lemma01.lean
Lemma02.lean
Lemma03.lean
Lemma04.lean
...
```

と100個の小さなファイルを作ると、逆に依存関係の管理が難しくなります。

私なら目安として、

> **1ファイル = 1つの数学的テーマ・概念・補題群**

とします。

例えば、

```text
PrimeFactorization.lean
```

の中に、

```lean
theorem factorization_exists
theorem factorization_unique
theorem factorization_finite
```

など、強く関連する定理をまとめます。

---

## AI形式化では、さらに「中間成果物」を保存するのが重要

ここがかなり重要です。

単にLeanコードだけを保存するのではなく、

```text
paper/
   original.pdf

analysis/
   dependency_graph.md
   definitions.md
   claims.md

lean/
   Basic/
   Lemmas/
   Main/
```

のように、

**① 論文の記述
② 数学的解釈
③ 依存関係
④ Lean実装**

を分離して保存するとよいです。

例えば、

```text
Claim 17
────────────
論文:
  ...

数学的意味:
  ...

必要な仮定:
  A, B, C

依存:
  Lemma 4
  Lemma 9
  Definition 12

Lean:
  theorem claim17 : ... := by
    ...
```

という情報を持たせます。

これをやっておくと、AIに後から

> 「Claim 17の証明が失敗した。どの仮定が不足している？」

と聞いたときにも追跡できます。

---

## さらに重要なのはGitとの組み合わせです

この種の作業では、

```text
Paper
  ↓
Claim decomposition
  ↓
Lean formalization
  ↓
Proof verified
```

という**小さな単位ごとにcommit**するのがおすすめです。

例えば、

```text
commit 1
Add definition of logarithmic height

commit 2
Prove basic properties of logarithmic height

commit 3
Formalize Lemma 4

commit 4
Formalize Lemma 5

commit 5
Connect Lemma 4 and Lemma 5

commit 6
Formalize Main Theorem
```

とします。

すると、AIが後から証明を壊しても、

> 「どこから壊れたか」

が明確になります。

---

### そして、これは数学AIにとってかなり重要な構造になります

最終的には論文を

```text
                    Main Theorem
                    /          \
                Lemma A       Lemma B
                /    \          /   \
             Def 1   Lemma C  Def 2 Lemma D
```

という**依存DAG**として表現できます。

Leanの`import`構造は、この数学的DAGをかなり直接的にコードへ反映できます。

そのため、以前話されていた**「AIエージェントによる数学研究」**を考えるなら、

> **論文 → 主張分解 → 依存DAG → Lean module → proof verification**

というパイプラインにするのが非常に面白いと思います。

特に、AIに「論文全体を一気にLean化させる」のではなく、**1つのClaimを独立した小さな検証可能単位として処理させ、Lean compilerをフィードバックループにする**方が、現在のLLMの能力にもかなり適した設計です。


依存関係を把握する目的で不完全なLeanファイル作成から開始することは考えられますか？骨組みのようなファイル体系や、容易に形式化できる部分は記述しておき、厳密化や省略されている部分の補完を後続の作業とするような方法です。


はい。むしろ**大規模な数学論文のLean形式化では、かなり合理的な方法**だと思います。

おっしゃっている方法は、単なる「未完成コード」ではなく、**論文の数学的依存関係を先にLean上へ投影する「スケルトン形式化」**と考えるとよいです。

### 典型的には3段階に分けられます

```text
Phase 1                  Phase 2                 Phase 3
構造化                   厳密化                   完全証明
────────────────────────────────────────────────────
論文の主張を抽出         定義を厳密化             証明を完成
依存関係を整理           仮定を明示               sorryを除去
ファイルを分割           型を確定                 Leanで検証
容易な定義を実装         省略部分を補完
```

例えば最初は、

```lean
-- ABC/Height/Basic.lean

def Height (x : ℚ) : ℝ := sorry

theorem height_nonneg (x : ℚ) :
    0 ≤ Height x := by
  sorry

theorem height_add (x y : ℚ) :
    Height (x + y) ≤ Height x + Height y + ... := by
  sorry
```

のような**「型と主張は確定しているが、証明は未完成」**の状態で構いません。

---

## これは依存関係を発見するために非常に有効です

例えば論文から、

> Theorem A follows from Lemma B and Lemma C.

を発見したら、

```lean
theorem theoremA :
    P := by
  sorry
```

としておきます。

さらに、

```text
Theorem A
 ├── Lemma B
 │    ├── Definition X
 │    └── Lemma D
 │
 └── Lemma C
      └── Definition Y
```

という依存関係を作っていく。

この段階では**証明できるかどうかは問題にしない**。

これが重要です。

---

# Leanを「依存関係データベース」として使う

これはかなり面白い考え方です。

例えば、

```lean
import ABC.Height.Basic
import ABC.Radical.Basic
import ABC.Valuation.Basic

theorem abc_main
    (a b c : ℕ)
    (h : a + b = c) :
    ... := by
  sorry
```

まで作ってしまう。

すると、

```text
ABC.Main
 │
 ├── Height.Basic
 │
 ├── Radical.Basic
 │
 └── Valuation.Basic
```

という依存関係が**Leanのmodule systemそのものに記録されます**。

つまり最初から完成した証明を作る必要はなく、

> **「この定理は、この定義・この補題に依存する」**

という情報だけを先に形式化できます。

---

# さらに「論文の曖昧さ」を発見できます

これはスケルトン方式の非常に大きなメリットです。

例えば論文では、

> By applying the previous lemma, we obtain (X).

と書いてある。

最初の段階で、

```lean
theorem previous_lemma : A → B := by
  sorry

theorem next_claim : C := by
  sorry
```

までは書けます。

ところが後から実際に

```lean
exact previous_lemma ...
```

を完成させようとすると、

```text
Aが必要
↓
しかし現在の仮定にはAがない
↓
論文の別のLemmaが必要？
↓
論文では暗黙に仮定されている？
```

という問題が発見されます。

つまり、

**スケルトン → 厳密化**

という段階を分けることで、論文の「暗黙の数学」を体系的に発見できます。

---

# `sorry` を使う方法はLeanでは特に適しています

Leanでは、

```lean
theorem foo : P := by
  sorry
```

のように未証明の定理を暫定的に置けます。

したがって、

```text
論文の主張
     ↓
Leanのstatementにする
     ↓
sorry
     ↓
依存関係を構築
     ↓
後でproofを埋める
```

という開発ができます。

ただし、**最終成果物では`sorry`をなくす**ことが目標になります。

また、`sorry`を残した状態でもLeanはファイルの型チェックを行えるので、**「statementの型が整合しているか」「依存関係が成立しているか」**を先に検証できます。

---

# さらに一歩進めるなら「3種類」に分けるとよいです

私は大規模形式化なら、単純な`完了/未完了`ではなく、

```text
[FORMALIZED]
完全にLeanで証明済み

[SKELETON]
statementと依存関係は確定、証明未完成

[TODO]
論文から数学的意味をまだ十分に解釈できていない
```

くらいに分けるのをおすすめします。

例えばコメントで、

```lean
/-
Status: SKELETON

Paper: Proposition 4.2
Depends on:
  - Lemma 3.7
  - Definition 2.4

TODO:
  - clarify implicit compactness assumption
  - formalize measurability condition
-/
theorem proposition_4_2 : ... := by
  sorry
```

とする。

これならAIエージェントにも扱いやすいです。

---

# AIによる形式化なら、さらに有効です

ここが今回の発想の特に面白いところです。

AIに論文全体を一気に形式化させるのではなく、

### Agent 1：論文解析

```text
論文
 ↓
Definition
Lemma
Proposition
Theorem
 ↓
依存関係DAG
```

### Agent 2：Skeleton作成

```text
DAG
 ↓
Lean modules
 ↓
theorem statements
 ↓
sorry
```

### Agent 3：既存mathlib探索

```text
各statement
 ↓
mathlib検索
 ↓
既存lemmaとの対応
```

### Agent 4：証明

```text
sorry
 ↓
proof generation
 ↓
Lean compiler
 ↓
error
 ↓
修正
```

という分業ができます。

---

## そして「ファイル体系そのもの」が研究成果になり得ます

最初は、

```text
Paper
│
├── Definition A
├── Lemma B
├── Lemma C
└── Theorem D
```

だったものが、

```text
Definitions
     │
     ├──────┐
     ▼      ▼
  Lemma A  Lemma B
     │      │
     └──┬───┘
        ▼
     Lemma C
        │
        ▼
    Theorem D
```

という構造になっていく。

そして形式化が進むにつれて、

> **「論文著者が提示した証明の構造」**

と

> **「Leanで実際に必要だった依存関係」**

の差が見えてきます。

これは非常に重要です。

---

### 私ならこの方式を採用します

特に大きな数学の場合は、

**① 論文を主張単位に分解
② Leanのstatementを作る
③ `sorry`で暫定的に接続
④ import/module構造で依存DAGを構築
⑤ mathlibとの対応を調査
⑥ 暗黙の仮定を発見
⑦ 定義を厳密化
⑧ 下位lemmaから順に証明
⑨ 最後に`sorry`をゼロにする**

という順番にします。

これは「論文をLeanに翻訳する」というより、**論文から形式的な数学の依存グラフを抽出し、それを徐々に証明可能な数学へ精錬していく**方法です。

そして以前話されていた**ABC予想のような巨大な形式化**では、こちらの方が「最初から完全な証明を書こうとする」より現実的だと思います。特にAIを大量に使う場合、**Skeletonを共通の中間表現（IR）のように扱う**設計はかなり相性が良いです。

