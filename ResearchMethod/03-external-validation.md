# 03. 外部検証 —— 文献・既存プロジェクトとの突き合わせ

`02-method.md` は Math_ABC3 の**内部経験だけ**から書いた。
本書はそれを外部の仕事と突き合わせ、**(1) 確認できたもの / (2) 本プロジェクト固有と思われるもの /
(3) 外部から学んで直すべきもの**に分ける。

調査日: 2026-08-26。

---

# 1. 外部で確認できたもの（＝方法論として一般性がある）

## 1-1. 依存グラフを中心に据えること —— **確認**

Lean Blueprint（Patrick Massot）が確立した方式で、
[PFR（Polynomial Freiman–Ruzsa）予想の形式化](https://terrytao.wordpress.com/2023/11/18/formalizing-the-proof-of-pfr-in-lean4-using-blueprint-a-short-tour/)、
Buzzard の Fermat 最終定理 5 年計画などが採用している。
blueprint は「人間可読の証明を Lean の形式化に紐づけた高水準の依存グラフ」であり、
**計画文書と分散チームの調整ツールを兼ねる**。

★本プロジェクトの「主張グラフ ＋ 層 ＋ 葉」はこれと同じ思想である。
★★ただし本プロジェクトは **`.src` locator を Lean 側に置いて原典ページと逐語を機械照合する**点が違う
（blueprint は LaTeX 側に置く）。
[LeanArchitect](https://arxiv.org/abs/2601.22554) が同方向（Lean のアノテーションから
blueprint を生成し依存を自動推論）に動いているので、この差は収束しつつある。

## 1-2. ★★「打点は地図を書き換える」—— **強く確認**

[LeanMarathon](https://arxiv.org/abs/2606.05400)（2026）が
**"evolving proof DAG"** を中核に据えている。同論文の記述:

> Long-horizon autoformalization of research mathematics fails not only at hard lemmas, but **at scale**:
> **statements drift, dependencies tangle, context decays, and local repairs corrupt distant work.**

> the system lets it **evolve, splitting an over-large node or repairing a misformalized one**,
> until each piece is small enough for one agent to discharge.

★★これは `01-paper-as-object.md` の 5 値のうち **② 割れた**（over-large node の分割）と
**③ 読み違えた**（misformalized node の repair）そのものである。
**独立に同じ結論に達している**ので、この部分は一般性があると考えてよい。

★同論文はさらに **"agent durability"** という枠で問題を立て、
「1 本の脆い多日運転を、多数の短く復旧可能で並列な運転に変える」と述べる。
★★これは本プロジェクトの多セッション運用の**より良い言語化**である（第 3 節で取り込む）。

## 1-3. 決定的な CI ゲート ＋ 役割を絞ったエージェント —— **確認**

LeanMarathon は **"four contract-scoped agents and a deterministic CI gate"**
（construct / audit / prove / repair）を採る。
★本プロジェクトの `check.mjs` ゲート＋区画分離（`Skeleton` / `Interface` / `Found` / `Check` / `Gap`）は
同じ役割を**役割分担ではなく型と規則**で果たしている。

## 1-4. ★★★検証の厳しさが使える自動化の量を決める —— **確認（Tao）**

Terence Tao の言明が `02-method.md` 第 IV 部の主張そのものである:

> **The level of automation and AI power that you can profitably use before it becomes slop is
> roughly proportionate to how stringent your verification is.**
> （[the-decoder](https://the-decoder.com/terence-tao-argues-ai-could-bring-division-of-labor-to-math-for-the-first-time-in-history/)）

★★★**この 1 文が本方法論の存在理由である。** 第 IV 部（計測器の 3 原則）を先に置く根拠になる。

## 1-5. 形式化は「抑圧された仮定」を暴く —— **確認**

> Natural-language proofs routinely **suppress measurability and topological assumptions**,
> conflate almost-sure and pointwise statements, and **compress multi-step arguments into informal phrases**,
> and formalization forces each of these gaps to be made explicit and resolved.
> （[AI4SLT](https://arxiv.org/pdf/2602.02285)）

★「informal phrases への圧縮」は本プロジェクトの**折り畳み**そのものである。
★★ただし外部は**それを暴く**とは言うが、**着手前に数える**とは言っていない（第 2 節参照）。

## 1-6. 空虚な形式化は名前のついた常習犯 —— **確認、かつ本プロジェクトの経験を補正**

> **Vacuous formalization** can occur where the Lean statement diverges from the textbook,
> with a common failure being **signature drift toward an overly loose statement** that admits
> unintended witnesses and makes downstream lemmas trivially provable but unfaithful.

★★★これは `02-method.md` の**手 3（仮定の充足性を先に確かめる）**と
ゲート規則 **G2（非空虚 witness）**を正当化する。

★★**ただし向きが 2 つある**。本プロジェクトが踏んだのは外部の記述と**逆向き**だった:

| 向き | 症状 | 実例 |
|---|---|---|
| **仮定が強すぎる** | 意図された模型で仮定が偽 → 定理が空虚 | ★本プロジェクト（`Skeletal 𝒟`、`IsLocalization`） |
| **主張が緩すぎる** | 意図しない witness を許し、下流が自明に通る | 外部が報告する "signature drift" |

★**両方を検査する必要がある。** `02-method.md` にこの 2 向きを明記する（第 3 節で修正）。

## 1-7. 読みの誤りが支配的なリスクである —— **確認、数値つき**

* 公開された Lean のベンチマークで **31.8% のエントリに形式化の誤り**があったという分析
  （[miniF2F-Lean Revisited](https://arxiv.org/pdf/2511.03108)）
* [FormalRx](https://arxiv.org/pdf/2607.04655) は形式化の意味的失敗に
  **28 カテゴリの誤り分類**を与えている

★★これは**読みの台帳**（`02-method.md` 第 I 部 3）を別台帳として持つことの強い根拠である。
「主張が正しいか」と「読みが正しいか」は別の失敗であり、別に管理すべきである。

## 1-8. 分業がはじめて可能になった —— **確認、かつ摩擦も確認**

Tao は「数学者は従来すべてを自分でやるほかなく、専門分化は選択肢ではなかった。
AI がはじめて分業をもたらしうる」と述べる。

★本プロジェクトの並行 2 セッション運用はまさにその実践であり、
**そして摩擦も実際に起きた**（一方が `git add -A` して他方のコミット鎖 37 本を飲み込んだ）。
★★分業の利得は本物だが、**作業領域の分割と生存検証を運用規則にしないと履歴が壊れる**。

---

# 2. 外部に見当たらなかったもの（＝本プロジェクト固有の候補）

★以下は「見つからなかった」であって「存在しない」ではない。ただし主要なサーベイと
2026 年の主要論文には現れなかった。

## 2-1. ★★★折り畳み台帳（省略の合図を**数える**）

LeanMarathon に対して明示的に確認した結果:

> **(5) Hedge Words & Omitted Steps** — No information provided about measuring or handling
> imprecise language ("routine," "immediate," "well-known") from source mathematics.

★外部は「形式化が抑圧された仮定を暴く」とは言うが、
**着手前に `routine` / `immediately` / `well-known` を数えて工数を見積もる**という運用は見当たらない。

★★★そして本プロジェクトが観測した**逆相関**（語が希少なほど重い。`routine` は 6 件しかないのに
最後まで残った唯一の未着手の葉がそこ）も、外部に対応する報告を見つけられなかった。

**→ これは方法論として持ち出す価値がある。**

## 2-2. ★★負の対照（negative control）をゲート規則にする

> load-bearing な宣言には、**性質を 1 つだけ落とした対照**を書き、それが破れることを確認する。
> **破れないならその性質は効いていない。**

★実験科学の negative control をそのまま形式化に持ち込んだもの。
外部の形式化プロジェクトで**ゲート規則として強制している**例は見つからなかった。

★★エージェントは仮定を足しがちなので、この規則は AI 主導の形式化で特に効くはずである。

## 2-3. ★穴の 3 分類 ＋ falsifier 必須

「閉じた / 原典の穴の疑い / 反証済み」の 3 分類と、
**②を名乗るには falsifier（何が示せれば③になるか）を書かねばならない**という規則。

★外部には「gap を見つけた」「誤りを見つけた」の報告
（例: [two Higgs doublet model の文献の誤りを Lean で特定](https://arxiv.org/pdf/2603.08139)）はあるが、
**未解決の穴を分類して falsifier つきで台帳に載せる**運用は見つからなかった。

## 2-4. ★逸脱台帳に「消費者の実測」を必須にする

仮定を足したとき、**その仮定に依存する下流を機械で数えて記録する**。
本プロジェクトの `Proposition 1.6` では「我々のツリー内に消費者はいない。原典側の消費者は [IUTchI]」
まで確定させた。★これにより逸脱の**危険範囲が確定**する。

## 2-5. ★計測器の自己検査を毎回走らせる

本プロジェクトのゲートは毎回まず**19 種類の壊れた入力**を食わせて落とせることを確かめる。
★ソフトウェア工学では常識だが、形式化プロジェクトの**規範として**明記された例は見つからなかった。

---

# 3. 外部から学んで直すべきもの（`02-method.md` への修正）

## 3-1. ★★形式化には「品質の段」がある —— 本方法論は最上段を前提にしている

Tao の整理:

> Formalization is not one activity but **a spectrum of quality tiers** ...
> from reusable, human-interpretable libraries built to publication standards at the top,
> down to **bare formalizations that merely certify a statement is true**, are not meant to be read or reused,
> and can be produced with heavy automation.

★★**本方法論は最上段（再利用可能・人間可読・原典への逐語照合つき）を前提にしている。**
`.src` locator、引用照合、`.needs`、Gap 台帳はすべてその段のためのコストである。

★**下の段を選ぶなら、6 つの台帳のうち 1・6 だけでよい。** 明記すべきだった。

| 段 | 要る台帳 | 使う場面 |
|---|---|---|
| 上（再利用可能・原典追跡可能） | 1〜6 すべて | ★本プロジェクト（最終目標が遠く、多論文、長期） |
| 中（正しいことは示すが読ませない） | 1（グラフ）＋ 6（在庫）＋ ゲート | 単一定理の検証 |
| 下（真であることの証明書だけ） | ゲートのみ | 探索・ふるい分け |

## 3-2. ★blueprint 先行 vs 押し切り —— **どちらも実在する**

[Liquid Tensor Experiment](https://leanprover-community.github.io/blog/posts/lte-final/) では
参加者は **「補題を次々に形式化して頂上まで押し切り」、blueprint は完成後に書いた**。
★本方法論（グラフ先行・葉から）とは**逆順**である。

**使い分けの見立て**（本プロジェクトの経験と LTE の記録から）:

| | blueprint 先行（本方法論） | 押し切り（LTE） |
|---|---|---|
| 目標 | 遠い・複数論文 | **単一の頂上が見えている** |
| 人員 | 多数 / エージェント / 長期 | 全体像を頭に持つ少人数 |
| 効く理由 | 葉が並列化でき、進捗が測れる | 計画のオーバーヘッドが無い |

★★**エージェントには blueprint 先行が要る。** 全体像を頭に保持できないからである
（LeanMarathon の "context decays" がこれを裏づける）。

## 3-3. ★折り畳みの分類に「暗黙の定数と量化子の順序」を足す

LTE の記録:

> "a very delicate argument fighting with estimates against homological algebra" with
> **"a hierarchy of implicit constants that have to be chosen"** in exactly the right order

> the main theorem ... **∀∃∀∃∀∃** ... "the most logically involved statement" the author had ever proved

★★**これは合図の語では捕まらない種類の折り畳みである。**
著者は `routine` とも `immediately` とも書かない —— ただ定数の選び順を書かないだけである。

**→ 折り畳み台帳に第 2 の軸を足す:**

| 軸 | 検出法 | 例 |
|---|---|---|
| **語による折り畳み** | 合図を数える | `routine` / `immediately` |
| ★**構造による折り畳み** | 主張の**量化子の交替数**と**暗黙の定数の個数**を数える | LTE の `∀∃∀∃∀∃`、定数の選び順 |

★本プロジェクトでは後者が問題にならなかった（FrdI は圏論的で不等式が少ない）が、
解析寄りの論文では**こちらが主戦場**になるはずである。

## 3-4. ★空虚性の検査を**両向き**にする（1-6 参照）

`02-method.md` の手 3 を「仮定が強すぎないか」だけでなく
**「主張が緩すぎないか（意図しない witness を許していないか）」**も含む形に読み替える。

## 3-5. ★「agent durability」という言語を借りる

LeanMarathon の

> turns **one brittle multi-day run into many short, recoverable, parallel ones**

は、本プロジェクトの多セッション運用の狙いを本方法論より正確に言っている。
★**「長時間走る」ではなく「短く区切って復旧可能にする」**が正しい定式化である。

---

# 4. 審判が無い場合について（射程外の見立て）

`02-method.md` 第 VI 部で「審判（形式化系）が無いと読みの誤りは検出できない」と書いた。
外部を見た結果、この見立ては**強まった**:

* 公開ベンチマークの **31.8%** に形式化の誤りがあった
* 意味的失敗に **28 カテゴリ**の分類が要る
* Tao —— 検証の厳しさが使える自動化の量を決める

★★★したがって **AI エージェントに数学研究作業をさせるなら、審判の構築が先である。**
審判の無い探索は Tao の言う "slop" になる。

★ただし審判は形式化系でなくてもよい。要件は
**(a) エージェントの信念から独立、(b) 自分が壊れていないことを示せる、(c) 主張の外に置かれている**
の 3 つである（`02-method.md` 第 IV 部）。

---

# 5. 総括 —— 本方法論のどこが持ち出せるか

| 項目 | 判定 |
|---|---|
| 依存グラフ中心・葉から着手 | 外部で確立。**そのまま使える** |
| 打点が地図を書き換える（5 値） | 外部（LeanMarathon）と独立に一致。**強い** |
| 決定的ゲート ＋ 区画分離 | 外部で確立。**そのまま使える** |
| 検証の厳しさ ∝ 使える自動化 | Tao が定式化済み。**方法論の第一原理に置く** |
| 読みを別台帳にする | 外部の誤り率データが裏づけ。**強い** |
| ★折り畳み台帳（合図を数える） | 外部に見当たらず。**持ち出す価値がある** |
| ★負の対照をゲート規則に | 外部に見当たらず。**持ち出す価値がある** |
| ★穴の 3 分類 ＋ falsifier 必須 | 外部に見当たらず。**持ち出す価値がある** |
| ★逸脱台帳に消費者の実測 | 外部に見当たらず。**持ち出す価値がある** |
| ★計測器の自己検査 | 工学では常識だが形式化の規範としては見当たらず。**明文化する価値がある** |
| 品質の段の明示 | ★**欠けていた**（3-1 で補う） |
| blueprint 先行 vs 押し切りの使い分け | ★**欠けていた**（3-2 で補う） |
| 構造による折り畳み（量化子・定数） | ★**欠けていた**（3-3 で補う） |

---

## 出典

- [Formalizing the proof of PFR in Lean4 using Blueprint: a short tour — Terence Tao](https://terrytao.wordpress.com/2023/11/18/formalizing-the-proof-of-pfr-in-lean4-using-blueprint-a-short-tour/)
- [LeanMarathon: Toward Reliable AI Co-Mathematicians through Long-Horizon Lean Autoformalization](https://arxiv.org/abs/2606.05400)
- [LeanArchitect: Automating Blueprint Generation for Humans and AI](https://arxiv.org/abs/2601.22554)
- [Terence Tao argues AI could bring division of labor to math for the first time in history](https://the-decoder.com/terence-tao-argues-ai-could-bring-division-of-labor-to-math-for-the-first-time-in-history/)
- [‘The job description is changing’: mathematician Terence Tao on the rise of AI — Nature](https://www.nature.com/articles/d41586-026-01246-9)
- [miniF2F-Lean Revisited: Reviewing Limitations and Charting a Path Forward](https://arxiv.org/pdf/2511.03108)
- [FormalRx: Rectify and eXamine Semantic Failures in Autoformalization](https://arxiv.org/pdf/2607.04655)
- [AI4SLT: Empirical Processes in Lean 4 for Formal Statistical Learning Theory](https://arxiv.org/pdf/2602.02285)
- [Formalizing the stability of the two Higgs doublet model potential into Lean: identifying an error in the literature](https://arxiv.org/pdf/2603.08139)
- [Completion of the Liquid Tensor Experiment — Lean community blog](https://leanprover-community.github.io/blog/posts/lte-final/)
