# 01. 論文は何種類あるか —— FrdI を証拠に

## 問い

> FrdI は大きな主張に向けて証明していく型か、いずれ必要になる道具を作っている型か。

## 答え —— **どちらでもない。第 3 の型である**

FrdI は

> **「望む操作が既存の枠組みでは*不可能である*ことを確認したうえで、その操作が可能になる新しい枠組みを作り、
> 新枠組みが旧枠組みを*忘れていない*ことを証明する」型**

の論文である。以下は原典の自己記述（物理 p.7、§I3）。

> the morphism p → pⁿ [where p is a prime number] clearly **does not extend to a ring homomorphism ℤ → ℤ**!
> That is to say, it is **difficult to see how to accommodate** such a "Frobenius morphism for number fields"
> **within the framework of scheme theory**.
> On the other hand, **if one works with monoids** as in the theory of log schemes,
> then such a morphism "p → pⁿ" **does indeed make sense**.

そして結論の 1 文（同 p.7）:

> One important motivation for the author in developing the theory of Frobenioids was the goal of developing
> **a geometric framework — i.e., roughly speaking, a geometry built up solely from "Galois theory" and "monoids" —
> in which a "Frobenius morphism on number fields" may be constructed.**

### この型の論理構造は 5 段である

| 段 | 内容 | FrdI での実体 |
|---|---|---|
| 1 | **遠い目標**を述べる | 数体版の Teichmüller 理論（§I3） |
| 2 | それに要る**操作**を特定する | 数体上の Frobenius 射 |
| 3 | ★**その操作が現行の枠組みで不可能なことを示す** | `p ↦ pⁿ` は環準同型 `ℤ → ℤ` にならない |
| 4 | 操作が生き残る**部分構造だけ**で枠組みを組み直す | Galois 圏 ＋ 単系 = Frobenioid |
| 5 | ★**新枠組みが旧枠組みを忘れていないことを証明する** | §3・§4・§6 の「category-theoreticity」定理群 |

★★**段 3 と段 5 がこの型の指紋である。**
「道具型」の論文には段 3 が無い（不可能性を経由しない）。
「主張型」の論文には段 5 が無い（枠組みを作り直さない）。

### 段 5 が論文の分量の大半を占める

これは形式化して初めて肌でわかった。FrdI の定理は、ほとんどが次のどちらかの形をしている。

* **「X は圏論的である」** —— 圏同値は X を保ち、かつ反射する
  （`Theorem 3.4` = Base と Frobenius 次数、`Theorem 4.2` / `Corollary 4.10` / `Corollary 4.11`）
* **「余計な自己同型が無い」** —— slim / Frobenius-slim / Div-slim
  （`Theorem 6.4, (i)`、`Theorem 6.2, (iv)`）

そして §6 の存在理由がここでわかる。§6「Some Motivating Examples」は
**枠組みが空虚でないことと、何も失っていないことの 2 つを同時に示す章**である:

| §6 の主張 | 段 5 のどの役割か |
|---|---|
| `Example 6.1` / `Example 6.3` —— 幾何・算術の Frobenioid が実在する | **非空虚性**（枠組みに中身がある） |
| `Theorem 6.4, (i)` —— それが rationally standard・slim・Div-slim | **剛性**（余計な自己同型が無い） |
| `Theorem 6.4, (iv)` —— `deg(Ψ) = 1` かつ **`L₁ ≅ L₂`** | ★★**復元**（数体そのものが取り戻せる） |

★`Theorem 6.4, (iv)` が段 5 の頂点である。**「Frobenioid が同型なら数体が同型」**
＝ 抽象化で数体を失っていない、という宣言そのものだからである。

★★そして `δ_A : Pic_Φ(A) ≅ ℝ`（Dirichlet 単数定理）も同じ役割を持つ ——
古典的な不変量（大域次数）が新枠組みの中でそのまま読めることを示している。
原文が §6 の冒頭（物理 p.109）で

> an interesting "category-theoretic interpretation" of some results of classical number theory,
> such as the **Dirichlet unit theorem** and **Tchebotarev's density theorem**

と書くのは、この意味である。

---

## 論文の 3 型（本プロジェクトで実際に見たもの）

| 型 | 成果物 | 「打点する」とは何をすることか | 例 |
|---|---|---|---|
| **A. 枠組み型** | 定義 ＋ 剛性定理 ＋ 非空虚な実例 | **定義が空虚でないこと・何も失っていないことを確かめる** | FrdI |
| **B. 橋型** | 2 つの主張の**同値** | 翻訳を**実効的に**する（どちらの向きに何が要るかを測る） | GenEll `Theorem 2.1`（abc ⟺ ℙ¹∖{0,1,∞} 版 abc） |
| **C. 帰結型** | 目標の主張そのもの | 前提を集めて**組み上げる** | GenEll `Corollary 4.4`、IUTchIV `Corollary 2.2` |

★★**型が違えば打点の作業内容が違う。** これは方法論の設計に直接効く:

* **A 型**では「実例を 1 つ作る」ことが定理の証明と同じくらい重い。
  実際、`Theorem 6.4` の 4 サブタスクのうち 3 つ（`𝒞^rlf` の実例・Div-slim・δ_A）は
  **すべて「算術の実例で確かめる」作業**だった。
* **A 型**では**仮定の充足性**が命取りになる。
  本セッションで `Skeletal 𝒟` を仮定に置きかけたが、算術の `𝒟 = (FinSub F K̄)ᵒᵖ` では
  **共役な中間体が同型だが相異なる**ので偽であり、置いていたら定理が空虚になっていた。
* **C 型**では逆に、実例は 1 つあれば足りる。組み上げの順序が主戦場になる。

---

## Wiles の比喩との突き合わせ

> 暗闇の中を手探りで家具の位置を探り、ある日電気がパッとつく。

本プロジェクトの状況は**これとは違う暗さ**である。区別しておく価値がある。

| | Wiles の暗闇 | 本プロジェクトの暗闇 |
|---|---|---|
| 家具はあるか | **わからない**（定理が成り立つかも未知） | **ある**（論文に書いてある） |
| 何が見えないか | 対象の存在と構造 | ★**著者が畳んだもの**と★**我々の読みが正しいか** |
| 明かりがつく瞬間 | 全体像の理解 | ★**抽象的な定義を具体形に置き換えられた瞬間**（後述） |

★★したがって「打点して道を炙り出す」は**半分正しい**。
正確には次である。

> **打点の帰結は 5 値であり、そのうち 4 つは地図を書き換える。**

| 帰結 | 意味 | 本セッションでの実例 |
|---|---|---|
| ① 閉じた | 節点が済になる | `pairing-vanishes` 以外の大半 |
| ② **割れた** | 見積りが甘く、節点が新しい鎖に分解された | `𝒞^rlf` は 3 回分解し直した |
| ③ **読み違えていた** | 原典の読みが違った（が主張は真） | `𝒞^rlf` はテンソルではなく `effR ⊤` だった |
| ④ **偽だった** | 我々の読みが**反証された** | `v-loc`（「`𝒟` は `𝒞` の局所化」）は `Example 4.3` で反例 |
| ⑤ 逸脱で閉じる | 仮定を足せば閉じる（債務として記録） | `Proposition 1.6` の (v)(vi) |

★①だけが「道が伸びた」。②③④⑤は**地図の修正**であり、
数としては②〜⑤のほうがはるかに多かった。

---

## 「明かりがつく瞬間」の正体（本セッションで 2 回起きた）

2 回とも**同じ形**だった。書いておく価値がある。

### 例 1 —— `𝒞^rlf`（実現化）

* 詰まっていた形: `Φ^rlf = ℝ≥0 ⊗_ℕ Φ`（テンソル）。
  これを divisorial にするには **`ℝ≥0` が ℕ 上平坦**が要り、凸幾何の議論になる。
* 明かり: 実現化した対象は**具体的に書ける** ——
  `Φ^rlf(L) = ℝ≥0[V(L)]` ＝「全素点で実係数の有効因子」＝ `grp := ⊤` の `ArithDatum`。
  ★テンソルも錐も平坦性も要らなくなった。
* ★★**そして原典自身が物理 p.115 でそう書いていた**:
  > the set of elements of (Φ^rlf)^gp(L) **with finite support** whose image under deg_L^arith is 0

  我々は**抽象的な定義**（テンソル積）を読んでいて、
  **著者が別の場所で書いた具体的な記述**を読んでいなかった。

### 例 2 —— `exponent θ ≠ 0`

* 詰まっていた形: 「`ℤ[θ]` が `𝓞 N` の中で有限指数」という仮定を運び続けていた。
* 明かり: **判別式**。`disc(θ)·𝓞 N ⊆ ℤ[θ]` は古典的事実で、
  mathlib に `Algebra.discr_mul_isIntegral_mem_adjoin` として**在庫があった**。

### 共通する形

> **抽象的な定義で詰まったら、その対象が*意図された模型で何であるか*を書き下す。
> 障害は数学ではなく、我々が選んだ符号化にあることが多い。**

★★これは「電気がつく」というより「**部屋を間違えていたと気づく**」に近い。
そして間違いに気づかせたのは、どちらの場合も**外部の審判**だった ——
例 1 では引用照合ゲートが p.115 を指し、例 2 では在庫索引が mathlib の補題を出した。
