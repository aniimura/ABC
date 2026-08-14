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
  /-- 望月氏の別論文への依拠。後ろ向き追跡の枝。

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
