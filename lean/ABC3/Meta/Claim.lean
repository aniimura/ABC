/-!
# 台帳の型 — 出典・実装待ち・飛躍

本プロジェクトは台帳を JSON ではなく **Lean の宣言そのもの**として持つ。
ノードとファイルが別々に生まれないので、被覆の乖離が構造的に起こりえない。

`tools/check.mjs` は、ここで定義した型の**存在**を機械的に検査する。
**型が付くかどうかは Lean が検査する**——この分業を守ること。
check.mjs は「`Foo.nonvacuous` という宣言があるか」しか見ておらず、
それが本当に `Nonempty Foo` を証明しているかは `lake build` が保証する。
-/

namespace ABC3.Meta

/-- 原典の一つの主張への参照。`1_Structured` の構造単位と 1:1 に対応する。

`sectionId` は `1_Structured/<論文>/section-N.html` の `<section id="...">`。
これにより Lean の宣言から、PDF 目視確認済みの逐語まで一本で辿れる。 -/
structure Source where
  /-- `ResearchPaper/papers.json` のタグ(例: `"pGC"`) -/
  paper : String
  /-- **物理**ページ番号(`pdftoppm -f N` が指すページ) -/
  pdfPage : Nat
  /-- 項目名(例: `"Proposition 1.1"`) -/
  item : String
  /-- `1_Structured` の `<section id="...">` -/
  sectionId : String
  deriving Repr, DecidableEq

/-- ある結果が Lean のエコシステムに存在するか。**実測値のみ**を書く(推測を書かない)。

mathlib だけを見るのは不十分だった(2026-08-14 実測)——mathlib に無いものが、
進行中の公開プロジェクトにはあることがある。逆に「公開プロジェクトにある」と
聞いて安心するのも危険で、`sorry` が残っていれば使えない。両者を型で分ける。

**測定は時とともに古くなる**。何をいつ測ったかは
`ResearchPaper/lean-ecosystem.json` に記録する。

## ★`absent` は探索範囲を伴う(2026-08-14 の実失敗から)

`Found/PGC/LocalFieldNorm.lean` に「mathlib に `charP_of_prime_eq_zero` 相当の補題は
**無い**(実測)」と書いたが、**誤りだった**——`CharP.charP_iff_prime_eq_zero` が
`Mathlib/Algebra/CharP/Basic.lean:103` に存在する。実際に見ていたのは
`Mathlib/Algebra/CharP/Defs.lean` の `ringChar` 名前空間だけで、`Basic.lean` は
見ていなかった。

**探索範囲を書かない「無い」は、再現も反証もできない。** それは測定ではなく印象である。
ゆえに `absent` は `searched` を要求する形にした——「どこを・どのパターンで探して
0 件だったか」を書けないなら、それはまだ `unmeasured` である。 -/
inductive LeanStatus
  /-- mathlib にある。宣言名を書く -/
  | inMathlib (decl : String)
  /-- mathlib 外の公開プロジェクトに **sorry 無しで**ある。移植/依存の判断が要る -/
  | inProject (repo item : String)
  /-- 公開プロジェクトで**作業中**(sorry が残る)。★独立に作ると重複投資になる -/
  | inProgress (repo note : String)
  /-- 実測して見つからなかった。

  `searched` には**探索範囲**を書く——どのリポジトリ/ディレクトリを、
  どのパターンで探して 0 件だったか。書けないなら `unmeasured` を使うこと。
  この引数が無かったために誤った「無い」が1件通っている(上記 docstring)。 -/
  | absent (searched : String)
  /-- まだ測っていない。**暫定であり放置しない**——`check.mjs` が件数を印字する -/
  | unmeasured
  deriving Repr

/-- **原典の証明が要求するもの**。原文の証明文から抽出する(推測で足さない)。

これを skeleton の時点で書くのは、規模を statement だけから見積もると
**下界にしかならない**ため——statement は、それが言及する対象しか表面化させない。
原文の証明文は「何に依拠するか」を既に書いているので、そこを拾う。

`page` は物理ページ番号(`pdftoppm -f N` が指すページ)。 -/
inductive ProofObligation
  /-- 原文が明示的に引用している外部文献の結果 -/
  | citation (ref item : String) (status : LeanStatus) (page : Nat)
  /-- 原文が典拠なしに「well-known」等として使っている事実。**大きさは未知** -/
  | folklore (what : String) (page : Nat)
  /-- 原文が段を飛ばしている箇所。`Gap` の候補 -/
  | implicitStep (what : String) (page : Nat)
  /-- **番号付き項目への依拠**(依存グラフの辺)。

  ## ★★辺の意味(2026-08-15 制定)

  **辺は「原文の証明・主張が実際に依拠しているもの」を意味する。**

  **辺である**:
  - 使う対象・記号・演算の**定義**
  - その主張を**確立している**結果
  - 証明の段が明示的に引く結果

  **辺ではない**:
  - **前方参照**(`cf. …, below`)——後で述べるものに、先に置かれた主張の証明は依拠できない
  - **解説への案内**(`the discussion of …`, `see also`)
  - **導入部から本体への指し**——導入部は証明を持たないので「依拠」という関係が無い

  ★★**種別では判定しない。役割で判定する。** `Remark` だから辺ではない、は**誤り**である
  ——`[IUTchIII] Remark 3.9.5, (i)` は**正則包を定義**しており、我々は実際に依拠している。

  ★この判定は**機械では検査できない**(`tools/check.mjs` 冒頭 A 群と同じ性質)。
  迷ったら `.implicitStep` に迷いを書くこと。**辺を残したいために「依拠している」と
  判定してはならない。**

  ## なぜ定めたか(実測)

  辺の意味が未定義だったため、原文側のグラフに **262 節点の強連結成分**が現れていた
  (`ResearchPaper/cycle-analysis.md`)。前方参照・解説案内・注釈からの辺を落とすと
  **262 → 8** になり、**残る循環はすべて1論文の中に閉じる**。
  すなわち**循環は理論の循環ではなく、辺の定義の副作用**だった。

  ★★**名前は `otherPaper` だが、`paper` は「同じ論文」でよい。**
  むしろ規模の本体はそちら側にある——2026-08-15 の実測で、
  現在の 11 本の辺のうち **6 本が同一論文内**(`[IUTchIII]` → `[IUTchIII]`)である。
  `tools/check.mjs` はこの2種を分けて数え、印字する。

  ★この名前は**誤解を招いた**: 名前から「別論文への枝」だけを書くものと読めるため、
  同一論文内の依存を `.derivation`(指す先を持たない自由文字列)で書いてしまい、
  結果として辺にならない、という書き方に誘導しうる。実際その誤読が起きた。
  改名は既存の使用箇所を壊すので行わず、ここに明記して機械側で分けて数える形にした。

  ★★`page` は **辺の先の論文における、その item の物理ページ**である
  (引用している側のページではない)。2026-08-15 の監査で、この2つが
  混在していたことが判明した——`cor_3_12.needs` の 6 本のうち
  `Theorem 3.11` だけが辺の先のページ(153)で、他の 5 本は
  引用している側のページ(173/174)を書いていた。
  三つ組 `(paper, item, page)` が自己完結して初めて
  `tools/check.mjs` が **辺の先を検査**できる(登記の有無・ページ範囲)。
  混在していると、検査は所有論文の範囲と比べるしかなく、事実上何も見ていない。 -/
  | otherPaper (paper item : String) (page : Nat)
  /-- 原文の中で導出されている段(外部依存ではないが作業は要る) -/
  | derivation (what : String) (page : Nat)
  deriving Repr

/-- **load-bearing** の印——主定理へ実際に消費される系列であること。

`Skeleton` の宣言 `foo` に `foo.loadBearing` を付けると、`tools/check.mjs` は
その宣言に **G3(負の対照)** `foo.negControl` を要求する。

荷重配分の規律: G1〜G5 全部を全 skeleton に課すのは有害(薄く均一に儀式を撒くと、
効かせるべき箇所の検査が相対的に軽くなる)。**この印を付けたものにだけ全ゲートを課す。** -/
structure LoadBearing where
  /-- どの結果へ向けて消費されるか(消費者を名指しできないなら load-bearing ではない) -/
  consumer : String
  deriving Repr

/-- **負の対照**(G3)の記録。

主張している性質を **1つだけ** 落とした対照が、実際に破れることを確認した記録。
破れないなら、その性質は statement に効いていない——飾りである。

`witness` には、対照が破れることを示す Lean の宣言名を書く。 -/
structure NegControl where
  /-- 落とした性質(1つだけ) -/
  dropped : String
  /-- 対照が破れることを示す宣言名 -/
  witness : String
  deriving Repr

/-- `Interface` の構造体が、まだ非空虚 witness を持てない理由。

**空欄で済ませないための型**——「まだ作れない」は許されるが、
「何を待っているか言えない」は許されない。 -/
structure WaitingFor where
  /-- 何の実装を待っているか(例: `"局所類体論の相互法則"`) -/
  what : String
  /-- Track B のどの作業に対応するか -/
  trackB : String
  deriving Repr

/-- **論文に属さないノードの所属**——既存理論の登記(2026-09-03、第 1452)。

## ★なぜ要るか(実測)

`lean/ABC3/<bucket>/<X>/` の `<X>` のうち 6 つは論文名ではなく**理論名**であり、
そこに `Found/` の 1,676 本のうち 951 本(57%)が入っていた
(`GaloisRep` 459・`Arakelov` 375・`Divisor` 65・`NumberField` 20・`SixExp` 20・`ProL` 9)。
`.src` を見れば所属**論文**は分かるが、ディレクトリからは**何の理論か**が分からず、
Explorer 上でも依存グラフ上でも所属が見えなかった。

## ★★依存は包含ではない

理論を消費する論文の下へ移してはならない。
`GenEll` が mathlib を使っていても mathlib は `GenEll/` の下に無いのと同じで、
**使う側と使われる側は兄弟として並ぶ**。したがって理論は理論のまま置き、
「何の理論で、どの論文が消費し、mathlib に在庫があるか」をここに書く。

## ☆置き方

理論ディレクトリごとに `Theory.lean` を 1 本置き、`<X>.theory : Theory` を宣言する。
`tools/check.mjs` が G10 として「登記の無いディレクトリ」を落とす。 -/
structure Theory where
  /-- 何の理論か(例: `"楕円曲線の Galois 表現・Tate 一意化・Vélu の同種"`) -/
  what : String
  /-- この理論を消費している原典の論文タグ(例: `["GenEll"]`)。★複数あってよい。 -/
  consumers : List String
  /-- mathlib に在庫があるか。**実測値のみ**を書く(推測を書かない)。 -/
  mathlibStatus : String
  deriving Repr

/-- 飛躍の分類。**既定は `modelError`**——`sourceGap` は最後の手段。 -/
inductive GapClass
  /-- ① こちらのモデル化の誤り -/
  | modelError
  /-- ② 必要な数学が未構築 -/
  | missingMath
  /-- ③ 原典側の飛躍。`falsifier` を書けないなら、これを名乗ってはならない -/
  | sourceGap
  deriving Repr, DecidableEq

/-- 原典が確立していないと判明した段の記録。

`Gap/` の各構造体はこの記録を伴う。**`falsifier` は必須**——
「何が起きればこれが①②に落ちるか」を書けないうちは、まだ③ではない。 -/
structure GapRecord where
  source : Source
  classification : GapClass
  /-- 何が起きればこの分類が覆るか -/
  falsifier : String
  deriving Repr

end ABC3.Meta
