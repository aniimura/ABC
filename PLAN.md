# Math_ABC3 — 数学論文のLean形式化: 進め方の検討

2026-08-14 / 入力: `idea.md`(形式化の方式)・`idea2.md`(研究の姿勢)・`ResearchPaper/0_Source`(116論文・310,070行・PDF 114本)・`ResearchPaper/1_Structured`(21文書)

根拠は、本検討中に**実際に実行して確かめた**ものに限る(再現手順を各節に付す)。

---

## 0. 北極星

**望月新一氏のIUT全体——abc予想の証明を含む——のLean形式化。**

原典に飛躍が含まれる可能性を前提に置く。したがって本計画は純粋な転写作業ではなく、**未解決問題の研究を含みうる**。その局面での姿勢は `idea2.md`(Wilesの研究方法)を土台とする——要点は「**目的には頑固、方法には柔軟**」「失敗は問題の否定ではなく方法の否定」。運用への落とし込みは §5。

**「正しい」と結論することは目的ではない。** 原典の主張を忠実にLeanへ移し、証明できればそれが成果、移せない箇所が見つかればそれも等価な成果。ただし後者の場合、①こちらのモデル化の誤り ②必要な数学が未構築 ③原典側の飛躍、を厳格に区別して記録する(§5)。

---

## 1. 実測した3つの事実

設計はこの3点からほぼ一意に決まる。

### 事実1 — Leanは statement の中身を検査しない(実演済み)

`idea.md`の中心的な提案は「statement を先に置いて `sorry`、証明は後で」である。この方式が何を検査し何を検査しないかを、素のLean 4(mathlib不要)で確かめた。

論文の定義を写す際に条件を1つ取り違えた、という設定:

```lean
structure Tower where
  f          : Nat → Nat
  strictMono : ∀ n, f n < f (n + 1)
  bounded    : ∀ n, f n < 10        -- 誤って写した条件

theorem prop_paper_skeleton (T : Tower) : T.f 0 = 7 := by
  sorry
```

`Tower` は「狭義単調増加かつ有界な自然数列」を要求しており、実例が存在しない(`strictMono` から `n ≤ f n`、よって `f 10 ≥ 10` が `bounded 10` と矛盾。`no_tower : ¬ Nonempty Tower` として証明済み)。しかし:

- **skeleton はビルドが通る**(警告は `declaration uses 'sorry'` のみ)。
- 後日 `sorry` を埋めようとすると、**埋まってしまう**。矛盾からは任意の結論が出るため:

```lean
theorem prop_paper_proved   (T : Tower) : T.f 0 = 7   := absurd (⟨T⟩ : Nonempty Tower) no_tower
theorem prop_paper_nonsense (T : Tower) : (0:Nat) = 1 := absurd (⟨T⟩ : Nonempty Tower) no_tower
```

実際のビルド出力:

```
warning: Vacdemo/Basic.lean:20:8: declaration uses `sorry`
info: 'Demo.prop_paper_proved'   does not depend on any axioms
info: 'Demo.prop_paper_nonsense' does not depend on any axioms
Build completed successfully (4 jobs).
```

`sorry`無し・**公理依存ゼロ**。想定しうる最も強い受理ゲートを全部通り、それでいて同じ資格で `0 = 1` が証明できる。

**最も重要な点**: 原典の対象は存在する(望月氏は実際に構成している)。空虚なのは**我々の転写**である。そしてLeanはその差を報告しない。

> **結論**: skeleton方式は「Leanが検査する側(証明)」を後回しにして「Leanが検査しない側(statement)」に内容の100%を置く。**statement専用のゲートが要る。**

検出する検査は具体的に書ける——「仮説を満たす具体例を1つ構成せよ」:

```lean
example : Nonempty Tower' := ⟨{ f := id, strictMono := fun n => Nat.lt_succ_self n }⟩
```

**派生する設計規則**: 未構築の基礎を `axiom` で埋めてはならない。

```lean
axiom bounded_ax : ∀ n : Nat, n < 10
theorem anything_from_axiom : (0:Nat) = 1 := absurd (bounded_ax 10) (by decide)
-- info: depends on axioms: [bounded_ax]
```

`#print axioms` は依存を表示するが、その公理が**偽である**ことは教えてくれない。対して `structure` 版なら「`Nonempty` が示せるか」という機械にかけられる問いに変わる。

> **`axiom` は検出不能な破滅、`structure` は検出可能な空虚。** 未構築の基礎は必ず `structure`(条件付き形式化)で受ける。

なぜこれが人工的な例でないか: IUTの定義は同時条件の束である(mono-theta環境 = 位相群 + `Out` の部分群 + 共役類 + 両立条件)。5条件のうち1つを強く写せば型は付くのに実例は消える。しかも**本物の実例を1つ作るには手前の理論をほぼ全部作る必要がある**ため、「定義と定理だけ書いて実例は作らない」のが自然な進み方になる——それが空虚化の温床である。

再現: `scratchpad/vacdemo/`(`lake build` で上記出力)。Phase 0 で `lean/ABC3/Meta/Calibration.lean` としてリポジトリに入れる。

### 事実2 — `.txt` は対象の同定に使えない(PDF目視で確認済み)

| | `.txt`(`pdftotext`) | PDF目視 |
|---|---|---|
| [EtTh] **p.11** | `Z d=ef Gal(Y /X) (∼= Z).` | `Gal(Y/X)` の Y・X は**下線なし** |
| [EtTh] **p.43** Def 2.13 (i) | `..., Gal(Y /X) (∼= l·Z).` | `Gal(Y̲̲/X̲̲)` の Y・X は**下線2本** |

`.txt`上は同一文字列でありながら、片方は `≅ Z`、もう片方は `≅ l·Z`。テキストだけを見れば**論文が自己矛盾しているように見える**が、PDFを見れば別の対象についての別の主張であって矛盾ではない。p.43には `Ÿ̲̲`(二点＋下線2本)・`Y̲̲`・`Y̲` が同一ページに共存している。下線はベクター図形として描画されており、抽出精度を上げても原理的に取れない。

> - locator は **PDFのページ番号**を正とする。`.txt`の行番号は検索補助。
> - 記法が対象同定に関わる箇所は **PDF画像の目視確認**(`pdftoppm -png -r 150 -f N -l N`)。**自動化できない固定費**として見積りに乗せる。
> - 確認結果は失われない形で保存する。既存`1_Structured`のHTML(`<span class="ul1">`/`<span class="ul2">`)はこの目的に合っており、`grep 'class="ul2"'` で機械的に洗い出せる。

### 事実3補 — 「mathlib に無い」は Lean エコシステム全体では偽でありうる(2026-08-14 実測)

mathlib だけを測るのは不十分だった。公開されている形式化プロジェクトを clone して実測した(記録: `ResearchPaper/lean-ecosystem.json`、測定日つき——**測定は古くなる**)。

| | Lean | .lean | 定理 | sorry | 最終更新 |
|---|---|---|---|---|---|
| [ClassFieldTheory](https://github.com/kbuzzard/ClassFieldTheory)(Clay 夏の学校 2025) | v4.33.0-rc1 | 108 | 483 | **25** | 2026-07-31 |
| [LocalClassFieldTheory](https://github.com/mariainesdff/LocalClassFieldTheory)(CPP 2024) | v4.22.0-rc2 | 66 | 234 | **114** | 2025-07-02 |

**結果**: 我々の義務に対して **drop-in で使えるものは 0 件**。ただし「ゼロから作る」でもない——

- 局所類体論: ClassFieldTheory に群コホモロジー・Tate コホモロジー・Herbrand 商・抽象的な `reciprocityIso` がある。ただし具体形 Γ_K^ab ≅ (K^×)^∧ は未完成で、`LocalCFT/` は2ファイルのみ。
- ℚ_p の有限次拡大 → 局所体のインスタンス: `Instances.lean` にあるが **sorry 11 件**——我々が要る箇所がちょうど未完成。
- 局所 Tate 双対性・**高次分岐群**: 両プロジェクトとも **0 件**。mathlib にも無い。

**★予測は外れた**: 調査前の見込みは「2件は既存で埋まる」だったが、実測では 0 件だった。類推でなく実測を優先する規律(`idea2.md` ⑤「文献を答えでなく道具箱として読む」)の実例として残す。

**設計への反映**: `Meta/Claim.lean` の `MathlibStatus` を `LeanStatus` に置き換え、`inMathlib` / `inProject`(公開プロジェクトに sorry 無しで存在)/ `inProgress`(**作業中**——独立に作ると重複投資)/ `absent` / `unmeasured` の5値にした。`check.mjs` は `inProgress` の件数を警告として印字する。

**戦略への含意**: 局所類体論はコミュニティの進行中作業がある領域なので、Track B で独立に作るのは重複投資になりうる。一方 **高次分岐群は誰も作っていない**(mathlib 側に TODO コメントがあるだけ)——ここは我々が作る意味がある。

### 事実3 — mathlibに望月語彙は無いが、古典的な土台はある(実測済み)

mathlib v4.31.0-rc2(8,176ファイル)を実測:

| 概念 | 状態 |
|---|---|
| 混標数局所体(MLF) | **有**: `NumberTheory/LocalField/Basic.lean` の `IsNonarchimedeanLocalField` |
| 円分指標 | **有**: `NumberTheory/Cyclotomic/CyclotomicCharacter.lean` |
| 付値論一般 | **有**: 378ファイル |
| 分解群・惰性群 | **有**: `RingTheory/Valuation/RamificationGroup.lean` |
| 高次分岐群(下付き番号付け) | **無**。同ファイルに `TODO: Define higher ramification groups in lower numbering` |
| 上付き番号付け・Herbrandの定理 | **無**(0件) |
| 局所類体論・local Tate双対性 | 実質**無** |
| anabelian / Frobenioid / log-shell / semi-graph / tempered基本群 | **すべて0件** |

> 作業は質の違う2種類に分かれる。**(a) 古典的だがmathlibに無い数学の構築**(難しいが既知、確実に完成する)と、**(b) 望月氏独自の機構の転写**((a)が揃うまで型が付かない)。
> **(b)から始められない。しかし(a)を全部揃えてから(b)に入ると、何が本当に必要かが何年も分からない。** これを解くのが §3。

---

## 2. `idea.md` の評価

| 提案 | 判定 | 根拠 |
|---|---|---|
| 論文の順序でなく依存DAG順 | **採用** | ただしDAGは論文から読むのでなくTrack Aで**機械生成**する |
| `import`構造で依存DAGを表現 | **採用(限定付き)** | importは「どのファイルを見たか」であって「何を消費したか」ではない。実消費は依存宣言の走査で別途取る |
| `sorry`で骨組みを先に作る | **採用、ゲート必須** | 事実1。方式自体は正しいが、ゲート無しでは無意味な骨組みを大量生産する |
| 3ステータス | **修正: 2軸にする** | 「証明済み」と「意味がある」は独立(事実1)。1軸だと `prop_paper_proved` が最上位に見える |
| コメントで Status/Paper/Depends | **採用＋機械化必須** | locator は PDF に対して機械照合できる(事実2) |
| Agent 1〜4 の分業 | **部分採用** | 既存mathlib探索は望月語彙に対し全件0ヒットで無情報(事実3)。古典層(a)には有効、独自機構(b)には無効 |
| 小さい単位でcommit | **採用** | ★現状 `D:\Math_ABC3` はgitリポジトリでない。Phase 0 で `git init` |
| 中間成果物の分離保存 | **採用** | `0_Source`・`1_Structured` は良い。§8-2 の要判断あり |

`idea.md`に**書かれていないが必須**のもの: ①非空虚性の検査(事実1)、②PDFページ単位のlocatorと記法の明示(事実2)、③`axiom`の禁止、④飛躍の扱い(§5)。

---

## 3. 構成: 2トラックで上下から挟む

```
        [IUTchIII] Cor 3.12  ← 北極星
                │
   Track A ▼ 骨格(上から)          statement のみ。証明しない。
        Skeleton/                  未構築の基礎は Interface/ の structure で受ける
                │
        ┌───────┴───────┐
        │  Interface/   │  ← 「まだ無い基礎」の型だけの定義(axiom 禁止)
        └───────┬───────┘          ★両トラックの接続点
   Track B ▲ 基礎(下から)          Interface の実装を1つずつ本物にする
        Found/                     sorry無し。mathlib品質。
                │
             mathlib
```

**Track A(骨格・上から)**: Cor 3.12 から後ろ向きに statement を置く。証明しない。未構築の基礎は `Interface/` の `structure` で受ける。成果物は **import グラフ = 依存DAG**、`Interface/` 一覧 = **「下に何を作らねばならないか」の機械生成された目録**。

**Track B(基礎・下から)**: Track A が「実際に載っている」と示した `Interface` から本物を作る。sorry無し、mathlibへPRできる形。

**この形を採る理由**(★2026-08-14 訂正、下記):

`Interface` で受けることの効果は、**退化の可能性を statement から `Interface` へ局在させる**ことである(実証: `Check/PGC/InertiaDegeneracy.lean` → `InertiaDegeneracyMoved.lean`)。局在した退化が実際に消えるのは、その statement が依存する `Interface` が**すべて**本物になったときであって、1つ discharge した時点ではない。

また `Nonempty` は退化 witness で満たされるので、**G2 の通過は内容の保証ではない**——内容を保証するのは「本物の実例を入れて下流が自明化しないこと」であり、それは G3(負の対照)の側の問いである。

`axiom` との対比は生き残る: 偽の公理はそもそも「実例はあるか」という問いにできない。

> ### ★訂正の記録(消さない)
>
> ここには当初こう書いてあった——
>
> > Track B が `Interface` の実例を1つ供給した瞬間、それに依存する Track A の statement 群が一斉に非空虚性の検査を受ける。**空虚さが後から自動的に暴かれる構造**になっている
>
> `ResidueCardinality` を実際に discharge した後で検証したところ、**3点で不正確だった**(検証: `Check/PGC/DischargeFiresNothing.lean`、すべて `sorry` 無し・標準3公理のみ):
>
> 1. **非空虚性の検査は discharge 前から通っていた。** `Nonempty (ResidueCardinality p)` は退化した実例(`card := fun _ => p`)で示せる——Track B の実装は要らない。「実例を供給した瞬間に検査が走る」という出来事は**起こりようがなかった**。G2 が測るのは「我々が書き下した条件が充足可能か」であって「意図した対象と一致するか」ではない。
> 2. **discharge しても退化は消えなかった。** Corollary 1.3 は `RD` と `SC` の**2つ**の `Interface` に依存する。`degenerateSC_inertia_eq_top (RD) (K) : inertia RD degenerateSC K = ⊤` が示すとおり、**`RD` が何であっても**惰性群は `⊤` に潰れる。しかも潰れる理由は `Interface` の条件そのもの(`isPrimePow` から `card K > 1` が出て、`q = q^[Γ_K:H]` が `[Γ_K:H] = 1` を強制する)にあり、**`RD` の側では原理的に塞げない**。
> 3. **したがって「一斉に」は成り立たない。** 検査が意味を持つのは、その statement が依存する `Interface` が**すべて**本物になったときだけである。
>
> この訂正は、自分で書いた設計の中心的な主張が、自分で作った器具によって反証された例として残す。

**副次的だが重要な性質(`idea2.md` ①)**: Wilesは「解けなくても数学的な成果が残る問題を選ぶ」姿勢を取った。この2トラック構成はその性質を持つ——Track B の成果(高次分岐群・局所類体論等)はIUTの当否に関わらずmathlibへの貢献として残り、Track A の成果(依存DAGと`Interface`目録)は「IUT形式化に何がどれだけ要るか」の検証可能な見積りとして残る。**最終的な接続だけが、IUTが正しいことに依存する。**

---

## 4. skeleton の受理条件(ゲート)

| | 内容 | 対応 |
|---|---|---|
| **G1 出典** | `paper` / **PDFページ番号** / 項目番号 / **原文逐語** / **記法**(下線・点、PDF目視確認済みか)。逐語はPDF該当ページのテキストに機械照合 | 事実2 |
| **G2 非空虚witness** | 仮説を満たす実例を構成する。まだ作れない場合は「どの `Interface` の実装を待っているか」を明示(空欄不可)。★**退化 witness でも通る**(実測)——G2 は条件の充足可能性しか測らない。実例が意図した対象かは G3 の側で測る | 事実1 |
| **G3 負の対照** | 主張している性質を**1つだけ**落とした対照が**破れる**ことを確認。破れないなら、その性質は statement に効いていない | 事実1 |
| **G4 `axiom` ゼロ** | プロジェクト内 `axiom` 宣言を禁止。`#print axioms` は mathlib標準の3つ以外を許さない | 事実1 |
| **G5 弱化禁止** | 証明できないとき、原典の主張を**導ける形に弱めない**。不足は §5 の機構で明示的な仮説として型に出す | §5 |
| **G6 規模の記録** | `Skeleton` の `theorem`/`lemma` は `X.needs : List ProofObligation` を持つ。**原文の証明文から**、証明が要求するもの(引用・典拠なしの事実・暗黙の段・別論文への枝・原文内の導出)を抽出する。**空リストは省略ではなく主張** | 下記 |

**荷重配分**: G1〜G5 全部を全skeletonに課すのは非現実的で、有害でもある(薄く均一に儀式を撒くと、効かせるべき箇所の検査が相対的に軽くなる)。

> **主定理へ実際に消費される系列(load-bearing set)にのみ G1〜G5 を課す。それ以外は G1・G4 のみ。**

G4 は全件に課してよい(コストゼロで破滅的失敗を構造的に締め出せる)。

**ステータスは2軸**:

```
証明軸:  TODO ──▶ SKELETON(型が付いた, sorry) ──▶ PROVED(sorry無し)
較正軸:  UNGATED ──▶ G1 ──▶ GATED(G1-G5)
```

`PROVED × UNGATED` が事実1の `prop_paper_proved` の座標。1軸のステータスではこの状態が最上位に見えてしまう。

---

## 5. ★飛躍(gap)の扱い

北極星が「飛躍がありうる」を前提に置く以上、これは例外処理ではなく本工程の一部である。

### 5-1. 飛躍は「証明できない」ではなく「追加仮説」として型に出す

原典が「Lemma 3.5 より Proposition 3.6 は直ちに従う」と述べ、実際には従わない場合。

**やってはいけないこと**: Proposition 3.6 を導ける形に弱める(G5違反)。これは事実1と同型の失敗——型は付くが原典について何も言わない。

**やること**:

```lean
/-- [IUTchIII] p.NNN では Lemma 3.5 から直ちに従うとされるが、
    実際には次の追加仮説を要すると判明した箇所。 -/
structure Gap_3_6 (d : Lemma35Data) : Prop where
  extra : (足りない性質を、原典の語彙で正確に書いたもの)

theorem prop_3_6 (d : Lemma35Data) (H : Gap_3_6 d) : Conclusion d := ...
```

得られるもの:

- 何が足りないかが **Lean の項として曖昧さなく**述べられる(散文の「〜が不明瞭」ではない)。
- 下流の全定理が `H` に依存することが機械的に追える——**依存の伝播が可視化される**。
- `H` が後に証明されれば、下流は自動的に無条件になる。`H` が偽と判明すれば、下流が一斉に落ちる。
- **`H` 自体が独立した well-posed な数学の問題になる**——`idea2.md` の「問題を別の問題へ写す」の実装。Wilesが `Fermat` を `modularity` へ写したのと同じ操作を、機械可読な形で行う。

### 5-2. トリアージと、その非対称性

`sorry` が埋まらないとき、原因は3つに1つ:

- **① モデル化の誤り** — 我々の `Interface`・型設計が原典を反映していない
- **② 未構築の数学** — 必要な基礎がまだ無い(mathlib不在を含む)
- **③ 原典側の飛躍** — 原典が実際にその段を確立していない

**③ を主張するには ①② を排除しなければならない。この非対称性を運用に埋める**:

- **既定の分類は ①**。次に ②。③ は最後。これは `idea2.md` ④「失敗を問題の否定でなく方法の否定と考える」の直接の実装。
- ③ と分類するには最低限: (i) 該当箇所の全 locator を PDF で目視確認済み、(ii) **複数の独立な型設計**で試して同じ壁に当たる、(iii) `H` が原典のどの主張からも導けないことを具体的に述べられる。
- ③ の記録には必ず **falsifier**(何が起きればこれが①②に落ちるか)を書く。書けないなら、それはまだ③ではない。

### 5-3. 公開されている争点の所在と、ゲートとの一致

`0_Source` には二次文献が2本ある:

- `Scholze-Stix - Why abc is still a conjecture.txt`(693行)
- `Joshi - Final Report on the Mochizuki-Scholze-Stix Controversy.txt`(551行)

争点は **[IUTchIII] Theorem 3.11 → Corollary 3.12** に集中している。注目すべきは指摘の**形**である:

- Scholze–Stix(l.631): 「the critical [IUTT-3, Theorem 3.11] does not become false, but **trivial**」——偽ではなく**自明**になる。
- Joshi(l.28): Arithmetic Holomorphic Structure が「2つ以上あることを定量的に主張するには不十分」——**存在のwitness**が足りない。

この2つは、事実1から**独立に**導いたゲートと同じものである:

| 公開されている指摘 | 対応するゲート |
|---|---|
| 「偽ではなく自明になる」 | **G3**(負の対照: 性質を落としても結論が変わらないなら、その性質は効いていない) |
| 「witnessが足りない」 | **G2**(非空虚witness: 仮説を満たす実例を構成せよ) |

偶然ではない。**形式化が寄与しうるのは元々この種の問い**だからである。逆に言えばここが本プロジェクトの最大のレバレッジ点であり、Phase 2 の骨格は Cor 3.12 / Thm 3.11 まで降りることを優先すべき。

**正直な限定(この節の価値を守るために必須)**: Thm 3.11 を忠実に述べるには、その下の塔(mono-theta環境・log-shell・Frobenioid・tempered群)が要る。`Interface` で条件付きに述べることはできるが、**その `Interface` を我々が設計する以上、自明化が「原典のせい」か「我々の `Interface` が弱すぎるせい」かは、各フィールドが原典の逐語に裏付けられている(G1)ときにのみ切り分けられる**。ここを緩めると本節の価値は消える。

### 5-4. Wiles流の運用(`idea2.md` の落とし込み)

- **目的には頑固、方法には柔軟**。北極星(§0)は変えない。工程・ゲート設計・`Interface` の型設計・2トラック構成そのものは、測定に基づいて変えてよい/変えるべき。
- **一つの型設計に固執しない**(idea2.md ③)。同じ対象に複数の `Interface` 設計がありうるなら並べて比較し、捨てるのを躊躇しない。§5-2 の (ii) はこれを検査手順としても要求している。
- **文献を「答え」ではなく「道具箱」として読む**(idea2.md ⑤)。行き詰まったら、望月氏以外の文献・他の形式化プロジェクト・mathlibの周辺分野を「この構造を自分の問題のどこに移植できるか」の目で読む。
  **ただし二次文献(Scholze–Stix・Joshi等)の結論を、自分の結論として引用しない。** 彼らの指摘は「どこを見るべきか」の道具として使い、判定は自分の機械検査で行う。
- **「ほとんどの場合、アイデアはない」**(idea2.md 7)。仮説→計算→失敗→別の理論、という泥臭い往復が実態である前提で工程を組む。stop-loss は「進んでいないから諦める」ためでなく、**同じ道具で無限に粘るのを防ぐ**ために置く。

---

## 6. ディレクトリと道具

```
D:\Math_ABC3\
├─ CLAUDE.md                ← 北極星と規律。短く保つ(長い規律は守られない)
├─ PLAN.md                  ← 本書
├─ memory\                  ← 確定した前提の記録。MEMORY.md が索引。★作業開始前に読む
├─ ResearchPaper\
│   ├─ 0_Source\            ← .pdf = 権威 / .txt = 検索補助(事実2)
│   └─ 1_Structured\        ← 構造化HTML。記法を明示。§8-2 の要判断あり
└─ lean\
    ├─ lakefile.toml        ← mathlib v4.31.0-rc2
    └─ ABC3\
        ├─ Meta\            ← 台帳の型(Source/LoadBearing/NegControl/…)・★較正デモ
        ├─ Interface\       ← 「まだ無い基礎」の型(structure のみ、axiom 禁止)
        ├─ Skeleton\        ← Track A: 論文の主張(sorry)。**原典の主張だけ**
        ├─ Check\           ← ★我々のモデルについての検査(識別力・負の対照の witness)
        ├─ Found\           ← Track B: 実装済みの基礎(sorry無し)
        └─ Gap\             ← §5 の追加仮説(飛躍)。falsifier 必須
```

**`Check/` を分けた理由(2026-08-14、最初の skeleton を書いて判明)**: 識別力の検査や負の対照は「**我々のモデル**についての事実」であり、「**原典**が何を主張しているか」ではない。同じ場所に置くと主語が混ざる——本計画が最も警戒している失敗の形。`Skeleton/` の全宣言に出典(G1)を要求するゲートが、この混在を実際に検出した。

`Skeleton/` から `Check/` を import することは禁止(`check.mjs` が検査)。原典についての主張が、我々のモデルの検査に依存する形を構造的に塞ぐ。

**唯一のゲート** `tools/check.mjs`:

1. `lake build`
2. `sorry` 台帳(件数を必ず印字。隠さない)
3. **`axiom` 検査**(G4。プロジェクト内 `axiom` 宣言 = 即NG)
4. **locator照合**(G1。PDFページ番号 + 逐語引用を、そのページの `pdftotext` 出力に照合。記法で落ちる文字は照合対象から外す)
5. `Interface` 実装待ち一覧(= Track B の作業キュー)を生成
6. `Gap` 一覧(現在の分類①②③・falsifier・下流依存の件数)を生成

散文の進捗ログは書かない。現在地は `check.mjs` の出力(生成物)とする。

---

## 7. フェーズ

### Phase 0 — 器具を先に作り、器具自体を較正する ★2026-08-14 完了(selftest 13/13)

> **経緯(残す)**: 一度「Phase 0 完了」と記録したが**過大な申告だった**——完了していたのは 1_Structured 側のゲートだけで、Lean 側(G2・G3・Claim構造・Gap)は未実装だった。器具を作り切る前に Phase 1 の内容作業へ進んだのが原因で、**計画自身が警告していた順序の誤り**。ユーザーの確認(「Phase 0 は完了ですか？」)で発覚し、残りを実装して達成した。この経緯を消さずに残す——受理条件を自分で緩める失敗の実例として。

| | 内容 | 状態 |
|---|---|---|
| 1 | `git init` | 完了 |
| 2 | `lean/` 初期化、mathlib v4.31.0-rc2 | 完了・`lake build` 成功 |
| 3 | `Meta/Calibration.lean` — 事実1の実演を常設の証拠として配置 | 完了 |
| 4 | `Meta/Claim.lean` — `Source`/`LoadBearing`/`NegControl`/`WaitingFor`/`GapRecord` | 完了 |
| 5 | `Interface/` `Skeleton/` `Found/` `Gap/` — 各 bucket の規則を module docstring に | 完了 |
| 6 | `ResearchPaper/papers.json` — 論文の登記簿 | 完了(pGC・EtTh) |
| 7 | `1_Structured/README.md` — 構造化の規約 | 完了(初版) |
| 8 | `tools/check.mjs` — 唯一のゲート(下記の全項目) | 完了 |
| 9 | **器具の較正**(下記 D1–D13) | 完了・13/13 PASS |

**`check.mjs` が実装している検査**

| 対象 | 検査 |
|---|---|
| 1_Structured | S1 必須属性 / S2 論文の登記 / S3 ページ範囲 / S4 逐語の PDF 照合 / S5 目視確認日 / S6 id の一意性 |
| Lean | `lake build` / `sorry` 台帳(件数を必ず印字) |
| Lean | **G1** `Skeleton` の各宣言に `X.src : Source`。その paper・pdfPage・sectionId が `papers.json` と `1_Structured` に実在するか照合 |
| Lean | **G2** `Interface` の各 `structure` に `X.nonvacuous : Nonempty X` か `X.waiting : WaitingFor` |
| Lean | **G3** `X.loadBearing` を付けた宣言に `X.negControl : NegControl` |
| Lean | **G4** プロジェクト内の `axiom` 宣言を禁止(較正デモのみ免除) |
| Lean | `Gap` の各 `structure` に `X.record : GapRecord`(`falsifier` は型が必須にしている) |
| 出力 | `Interface` 実装待ち一覧(= Track B の作業キュー)/ `Gap` 一覧(分類・falsifier) |

**分業の原則**: check.mjs が見るのは**宣言の存在**だけ。`X.nonvacuous` が本当に `Nonempty X` を証明しているかは `lake build` が保証する。この境界を出力に明記させている。

**器具の較正(`node tools/check.mjs --selftest`)** — 意図的に壊した入力を落とし、正しい入力は通すか。

| | ダミー | 期待 | 結果 |
|---|---|---|---|
| D1 | 存在しないPDFページを指す locator | 落ちる | ok |
| D2 | 逐語が該当ページに無い | 落ちる | ok |
| D3 | 必須属性の欠落 | 落ちる | ok |
| D4 | 未登記の論文タグ | 落ちる | ok |
| D5 | 目視確認日が日付でない | 落ちる | ok |
| **D6** | **正しい入力(pGC §1)** | **通る** | ok |
| D7 | `Interface` に witness も waiting も無い | 落ちる | ok |
| D8 | プロジェクト内の `axiom` 宣言 | 落ちる | ok |
| D9 | `Gap` に `GapRecord` が無い | 落ちる | ok |
| **D10** | **非空虚 witness を持つ `Interface`** | **通る** | ok |
| **D11** | **waiting を書いた `Interface`** | **通る** | ok |
| D12 | load-bearing なのに負の対照が無い | 落ちる | ok |
| **D13** | **load-bearing + 負の対照あり** | **通る** | ok |

D6・D10・D11・D13(通るべきもの)を入れているのは、**「全部落とす」器具は「全部通す」器具と同じくらい無情報**だから。両方向を検査して初めて識別力があると言える。

**較正が実際に欠陥を捕まえた実例(3件)**

1. 初回実行で pGC §1 の逐語 **8件中8件が落ちた**。原因は原文でなく器具側——(a) `pdftotext` は下付き文字の直後に空白を挿入する(`Γ_K.` → `ΓK .`)、(b) 既定モードは行の順序を入れ替えることがある(pGC p.3 の Corollary 1.3)。空白を無視した照合と `-layout` 併用で解消。**器具を先に作らず Lean を書き始めていたら、この2つは「原文と合わない」という誤った結論として現れていた。**
2. `axiom` 免除パスの相対パス基準が食い違っており、較正デモ自身が G4 で落ちた。
3. **D13(通るべき fixture)が落ちた**——台帳の付随宣言(`X.loadBearing` 等)自身に出典を要求してしまう偽陽性。**落とす側の検査だけを並べていたら発見できなかった。**

### Phase 1 — 1本を完走させる

**1つの結果を skeleton → G1〜G5 → 証明完成(`sorry` 0)まで、工程を飛ばさずに通す。**

第一の成果物は定理そのものではなく**実測値**である: G2/G3 の設計にかかる時間、PDF目視の1命題あたりの固定費、`Interface` に切り出す粒度、locator照合の実用性。**この実測なしに「IUT形式化に何年」と言うのは根拠のない数字になる。**

対象の選定は §8-1。

### Phase 2 — 骨格を張る(証明しない)

Cor 3.12 から後ろ向きに `Skeleton/` を展開。**この段階では G1・G4 のみ**。

- 成果物: import グラフ = 依存DAG、`Interface/` 一覧 = 基礎の目録。
- ここで初めて、IUT形式化の規模が**推測でなく計数で**言える。
- §5-3 より、**Thm 3.11 / Cor 3.12 への到達を優先**する。

> **★2026-08-14 訂正(Phase 1 の最初の skeleton から判明)**: 「Track A の成果物 = `Interface` 一覧 = 基礎の目録」は**言い過ぎだった**。
>
> pGC Proposition 1.1 の skeleton を書いたところ、`Interface` が **1件も立たなかった**。statement が言及する対象(絶対 Galois 群・Krull 位相・円分指標・`ℚ_[p]`)は全て mathlib にあり、不在なのは**証明に要る**もの(局所 Tate 双対性)だったため。
>
> **statement は、それが言及する対象しか表面化させない。** 証明に要るだけの数学は skeleton 段階では見えない。したがって Phase 2 が生む目録は「基礎の全目録」ではなく「**statement を書くために要る対象の目録**」であり、規模の下界にしかならない。IUT 方向では statement 自体が Frobenioid・log-shell を名指しするので状況は変わるが、それでも**証明側の不足は Phase 4 まで見えない**。規模の見積りを Phase 2 の出力だけで出さないこと。

### Phase 3 — 下から埋める

`Interface` を Track B で本物にする。load-bearing なものから。実装1件ごとに上の skeleton 群が非空虚性の検査を受ける(§3)。ここで `Gap` が立ち始める(§5)。

### Phase 4 — 証明を埋める

`sorry` を下位から消す。最も長い。Phase 1 の実測値で見積りを更新しながら進む。

---

## 8. 決定事項と、残る検討課題

> **2026-08-14 ユーザー決定**
> - **8-1 Phase 1 の対象 = [pGC]**(決定)
> - **8-2 `1_Structured` = 再作成で対応**(決定)。`0_Source` は今後増えるため、`1_Structured` も本プロジェクトの中で作成・拡充を続ける。切れた参照は「除去」ではなく「再作成」で解消する。
> - **8-3 Phase 1 の stop-loss = 継続検討課題**(未決)

### 8-1. Phase 1 の対象 — 決定: [pGC]

| 候補 | 見通し |
|---|---|
| **[pGC] §1**(p進局所体のGrothendieck予想、1997、8ページ) | 対象(Galois群・円分指標・惰性群)はmathlibに**有**。不在は local Tate双対性・局所類体論・高次分岐群という**古典的で確実に作れる**もの。論文が短く完走が現実的 |
| **[AbsTopIII] §3 Def 3.1/Prop 3.2**(MLF-Galois T-pair) | mathlibに `IsNonarchimedeanLocalField` があり入口はある。周辺(monoid cyclotome等)の不在が深い |
| **[SemiGraph] Def 3.1-3.4**(anabelioidの半グラフ) | mathlibに関連物ゼロ。土台から全部作ることになる |

**推奨: [pGC] §1**。IUTに近いからではなく、**Phase 1 の目的が手法の較正だから**——完走できる見込みが最も高く、途中で必要になるもの(高次分岐群)はmathlib側でも明示的に求められている(事実3の TODO)ので副産物が無駄にならない。`idea2.md` ①(解けなくても成果が残る問題を選ぶ)にも合う。臨界パス上かどうかは Phase 2 で機械的に判定できるので、Phase 1 の選定基準にする必要はない。

### 8-2. `1_Structured` — 決定: 再作成

`0_Source` は今後増える。したがって `1_Structured` は「引き継いだ成果物」ではなく、**本プロジェクトが所有し、繰り返し実行する工程**である。旧構成のファイル(全21件が存在しないフォルダを参照していた)は、参照を除去するのではなく**本プロジェクトの規約で作り直す**。

規約は [`ResearchPaper/1_Structured/README.md`](ResearchPaper/1_Structured/README.md)(2026-08-14 初版)。要点:

- **ゲート G1 の入力になる形**で作る——読みやすいHTMLを作ることが目的ではない。`data-paper` / `data-pdf-page` / `data-item` / `data-notation-checked` を必須属性とし、`.verbatim` に原文逐語、`.reading` に我々の読みを**分けて**入れる。
- **`data-txt`** で「その部分の `pdftotext` 出力」を明示する(`<span class="bar" data-txt="K">K</span>` 等)。これにより、装飾を復元した逐語でも PDF に対して機械照合できる。
- 論文の登記簿は [`ResearchPaper/papers.json`](ResearchPaper/papers.json)(タグ・ファイル名・総ページ数・`pageOffset`・`notationRisk`)。`0_Source` に論文を足したら、まずここに1件追加する。
- 旧構成のファイルは **`*.legacy.*`** にリネームし、`check.mjs` の対象外にする。内容は再作成の下書きとして残すが、**そこに書かれた事実を検証せずに引き継がない**。

### 8-3. Phase 1 の stop-loss — 継続検討課題(未決)

`idea2.md` の姿勢(方法には柔軟)を実装するには、「同じ道具で粘り続けない」ための閾値が要る。候補: 経過時間 / `Interface` の増殖数 / 同じ `sorry` への再挑戦回数。

Phase 1 を進めながら実測値を集め、**根拠を持って**決める(いま数字を決めても根拠が無い)。決めるまでの間は、閾値の代わりに「`Interface` が増えたら必ず `check.mjs` の作業キューに現れる」という可視化で代用する。

---

## 9. 正直な制約

- **PDF目視は自動化できない**(事実2)。対象を増やすたびに人手の固定費が乗る。「後で効率化できる」と見込んではいけない。
- **Track B の規模は Phase 2 が終わるまで分からない**。現在言えるのは「望月語彙はmathlibに0件」(事実3)までで、person-year 換算は `Interface` 目録が出るまで推測にすぎない。**見積りを Phase 2 の前に出さない**のが正直な扱い。
- **`Interface` の `structure` は実装が来るまで空虚でありうる**。これは設計通り(空虚は検出可能・`axiom`の偽は検出不能)。ただし「`Interface` を置いた」ことを「形式化した」と呼ばない。
- **本計画はIUTの当否を判定する装置ではない**。§5-3 の一致は、この方式が争点と同じ形の問いを機械にかけられることを意味するだけで、判定できることを意味しない。判定には忠実性(G1)の担保が前提で、それは目視確認という人手の工程に依存している。
- **③(原典側の飛躍)は、これまで一度も確認されていない**。既定は①、次に②。§5-2 の要件を満たさない③の主張はしない。
- **`Interface` を1つ discharge しても、同じ statement が依存する他の `Interface` が退化を再導入できる**(実証済み、§3 の訂正記録)。discharge の効果は「依存する `Interface` が全部揃うまで測れない」。
- **G2 は名前で回避できる。** `X.nonvacuous` という名前を付けなければ実装待ちに数えられ、付ければ(退化 witness でも)キューから消える。実際 `SubgroupCorrespondence` の `Nonempty` は退化 witness で証明済みだが、**意図的に `.nonvacuous` と名付けていない**——名付けた瞬間、実装ゼロのまま作業キューが空になるため。この判断は人手であり、機械では担保されていない。

---

## 10. 一行まとめ

> `idea.md` の skeleton 方式を採る。ただし **Leanは statement の中身を検査しない**(実演済み)ので、statement 専用のゲートを先に作り、その較正を通してから1本を完走させる。基礎は `axiom` でなく `structure` で受け、上からの骨格と下からの実装で挟む。飛躍は「証明できない」ではなく**追加仮説として型に出し**、①②③を厳格に区別する——`idea2.md` の「目的には頑固、方法には柔軟」を、機械が強制できる形にしたものが本計画である。
