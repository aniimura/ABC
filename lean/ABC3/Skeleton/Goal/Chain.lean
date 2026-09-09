import ABC3.Meta.Stmt
import ABC3.Skeleton.PGC.Section2
import ABC3.Skeleton.PGC.Section3
import ABC3.Skeleton.PGC.Section4
import ABC3.Skeleton.GenEll.Section3
import ABC3.Skeleton.GenEll.Section4
import ABC3.Skeleton.GenEll.GaloisImage

/-!
# 還元の登記 —— 「証明の前に依存を型で出す」

## なぜこのファイルが要るか

Skeleton の主張どうしは statement のうえでほぼ繋がっていない
(2026-09-09 実測: 247 件中 **2 件**しか他の主張に触れない。`tools/stmt-edges.py`)。
★原典に忠実な結果である——原典でも「Lemma 3.1」は Proposition 3.4 の
**証明**に出て**主張**には出ない。⇒ 依存は証明を書くまで機械に見えない。

★本ファイルは、原典が**明示している**依存を、証明を書く前に**含意として型にする**。

```
X_of_Y (h : stmt% Y) : stmt% X := sorry     -- ★還元。未証明
chain_X : stmt% X := X_of_Y @Y              -- ★ここが型検査され、定数 Y を消費する
```

`chain_X` が `Y` を消費するので、依存グラフに辺 `Y → X` が**証明 0 のまま**出る。
測定は `node tools/goal-chain.mjs`。

## ★このファイルが主張していないこと

- 各 `_of_` の**数学的内容は空欄**である(`sorry`)。本体は原典の該当証明を読んでいない。
- 辺は「原典がその依存を明示している」という**書誌的事実**であり、
  形式化された証明が実際にその経路を通る保証ではない
  (例: 木の `prop_2_2` は原典と違って有限次降下を使っていない)。

## 辺の根拠(すべて測ってから書いた)

| 辺 | 根拠 |
|---|---|
| `prop_2_1 → prop_2_2` | 木 `Skeleton/PGC/Section2.lean` の docstring が引く原文「有限次拡大 `L/K` へ降りて Prop 2.1 を使い」 |
| `prop_2_2 → cor_3_1` | 原典 §3「Next, let us observe the following **formal consequence of Proposition 2.2**」(`1_Structured/…/section-3.html`) |
| `cor_3_1 → cor_3_3` | 同 §3「Corollary 3.1 と同型の論法: uniformizing 性は `d_V(1)`・`d_V(0)` という(Corollary 3.1 により回復可能な)不変量だけで判定できる」 |
| `cor_3_3 → theorem_4_2` | 原典 §4 本文「**by Corollary 3.3**, it follows that `V` is also a uniformizing …」 |
| `lemma_4_1 → theorem_4_2` | 原典 §4 本文「generates `K` … (**by Lemma 4.1 below**)」 |
| `theorem_3_8 / lemma_4_1 / lemma_4_2 / prop_3_4 → cor_4_3` | `cor_4_3.needs` の `.otherPaper` 4 件 |
| `theorem_3_8 / cor_4_3 / lemma_4_1 → cor_4_4` | `cor_4_4.needs` の `.otherPaper` 3 件 |
-/

namespace ABC3.Skeleton.Goal.Chain

open ABC3.Meta

/-- ★同じ主張への 2 つの経路を 1 つにまとめる。返すのは第 1 引数(木の定理)。

第 2 引数(原典が明示する依存経路)は使わないが、**適用項に残る**ので
依存グラフには**両方の辺**が出る。★これが目的である——
「木の頂点」と「原典が言う経路」のどちらが欠けても、繋がっていることにならない。 -/
def viaBoth {P : Prop} (tree _paper : P) : P := tree

/-! ## ① pGC —— §2 → §3 → §4 -/

/-- **[pGC] Prop 2.1 ⟹ Prop 2.2**。★未証明。 -/
theorem prop_2_2_of_prop_2_1 (_h : stmt% ABC3.Skeleton.PGC.prop_2_1) :
    stmt% ABC3.Skeleton.PGC.prop_2_2 := sorry

theorem chain_prop_2_2 : stmt% ABC3.Skeleton.PGC.prop_2_2 :=
  prop_2_2_of_prop_2_1 @ABC3.Skeleton.PGC.prop_2_1

/-- **[pGC] Prop 2.2 ⟹ Cor 3.1**。原文「formal consequence of Proposition 2.2」。★未証明。 -/
theorem cor_3_1_of_prop_2_2 (_h : stmt% ABC3.Skeleton.PGC.prop_2_2) :
    stmt% ABC3.Skeleton.PGC.cor_3_1 := sorry

theorem chain_cor_3_1 : stmt% ABC3.Skeleton.PGC.cor_3_1 :=
  cor_3_1_of_prop_2_2 chain_prop_2_2

/-- **[pGC] Cor 3.1 ⟹ Cor 3.3**。★未証明。 -/
theorem cor_3_3_of_cor_3_1 (_h : stmt% ABC3.Skeleton.PGC.cor_3_1) :
    stmt% ABC3.Skeleton.PGC.cor_3_3 := sorry

theorem chain_cor_3_3 : stmt% ABC3.Skeleton.PGC.cor_3_3 :=
  cor_3_3_of_cor_3_1 chain_cor_3_1

/-- **[pGC] Cor 3.3 と Lemma 4.1 ⟹ Theorem 4.2**。原文 §4 本文が両方を名指ししている。★未証明。 -/
theorem theorem_4_2_of_cor_3_3 (_h33 : stmt% ABC3.Skeleton.PGC.cor_3_3)
    (_h41 : stmt% ABC3.Skeleton.PGC.lemma_4_1) :
    stmt% ABC3.Skeleton.PGC.theorem_4_2 := sorry

/-- ★**pGC の頂点**。`chain_cor_3_3` 経由なので、§2 まで辺が繋がっている。 -/
theorem chain_theorem_4_2 : stmt% ABC3.Skeleton.PGC.theorem_4_2 :=
  theorem_4_2_of_cor_3_3 chain_cor_3_3 @ABC3.Skeleton.PGC.lemma_4_1

/-! ## ② GenEll —— §3 → §4 -/

/-- **[GenEll] Theorem 3.8 / Lemma 4.1 / Lemma 4.2 / Prop 3.4 ⟹ Corollary 4.3**。
`cor_4_3.needs` の `.otherPaper` 4 件をそのまま仮説にした。★未証明。 -/
theorem cor_4_3_of_parts
    (_h38 : stmt% ABC3.Skeleton.GenEll.theorem_3_8)
    (_h41 : stmt% ABC3.Skeleton.GenEll.lemma_4_1)
    (_h42 : stmt% ABC3.Skeleton.GenEll.lemma_4_2)
    (_h34 : stmt% ABC3.Skeleton.GenEll.prop_3_4) :
    stmt% ABC3.Skeleton.GenEll.cor_4_3 := sorry

theorem chain_cor_4_3 : stmt% ABC3.Skeleton.GenEll.cor_4_3 :=
  cor_4_3_of_parts @ABC3.Skeleton.GenEll.theorem_3_8 @ABC3.Skeleton.GenEll.lemma_4_1
    @ABC3.Skeleton.GenEll.lemma_4_2 @ABC3.Skeleton.GenEll.prop_3_4

/-- **[GenEll] Theorem 3.8 / Corollary 4.3 / Lemma 4.1 ⟹ Corollary 4.4**。
原文の証明は 3 行で「Corollary 4.3 とまったく同様」。★未証明。 -/
theorem cor_4_4_of_cor_4_3
    (_h38 : stmt% ABC3.Skeleton.GenEll.theorem_3_8)
    (_h43 : stmt% ABC3.Skeleton.GenEll.cor_4_3)
    (_h41 : stmt% ABC3.Skeleton.GenEll.lemma_4_1) :
    stmt% ABC3.Skeleton.GenEll.cor_4_4 := sorry

/-- ★**GenEll の頂点**。`chain_cor_4_3` 経由なので §3 まで辺が繋がっている。 -/
theorem chain_cor_4_4 : stmt% ABC3.Skeleton.GenEll.cor_4_4 :=
  cor_4_4_of_cor_4_3 @ABC3.Skeleton.GenEll.theorem_3_8 chain_cor_4_3
    @ABC3.Skeleton.GenEll.lemma_4_1

/-- ★**pGC の頂点を 2 経路から取る**。`theorem_4_2` は先頭に暗黙束縛 `{p}` を持つので、
`viaBoth` を頭で当てると `p` が先に具体化されて型が合わない(2026-09-09 実測、#356)。
★結論まで eta 展開してから当てる。 -/
theorem chain_theorem_4_2_both : stmt% ABC3.Skeleton.PGC.theorem_4_2 :=
  fun {p : ℕ} [inst : Fact (Nat.Prime p)] (K K' : ABC3.Skeleton.PGC.PAdicLocalField p) =>
    viaBoth (ABC3.Skeleton.PGC.theorem_4_2 (p := p) K K') (chain_theorem_4_2 (p := p) K K')


/-! ## 出典(G1)と、証明が要求するもの(G6)

★`.src` の item は「原典の項目そのもの」ではなく**還元**である。
辺の根拠はファイル冒頭の表に書いた。 -/

def prop_2_2_of_prop_2_1.src : Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2 (Prop 2.1 からの還元)",
    sectionId := "prop-2-2" }

def prop_2_2_of_prop_2_1.needs : List ProofObligation :=
  [ .implicitStep
      ("原典 §2 が「有限次拡大 L/K へ降りて Prop 2.1 を使い、上付き→下付き→上付きと番号付けを往復する」と述べる段。★本体は読んでいない。") 5 ]

def chain_prop_2_2.src : Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2 (還元の合成)",
    sectionId := "prop-2-2" }

def chain_prop_2_2.needs : List ProofObligation := []

def cor_3_1_of_prop_2_2.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1 (Prop 2.2 からの還元)",
    sectionId := "cor-3-1" }

def cor_3_1_of_prop_2_2.needs : List ProofObligation :=
  [ .implicitStep
      ("原典 §3「Next, let us observe the following formal consequence of Proposition 2.2」。★本体は読んでいない。") 6 ]

def chain_cor_3_1.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1 (還元の合成)",
    sectionId := "cor-3-1" }

def chain_cor_3_1.needs : List ProofObligation := []

def cor_3_3_of_cor_3_1.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.3 (Cor 3.1 からの還元)",
    sectionId := "cor-3-3" }

def cor_3_3_of_cor_3_1.needs : List ProofObligation :=
  [ .implicitStep
      ("原典 §3「uniformizing 性は d_V(1)・d_V(0) という(Cor 3.1 により回復可能な)不変量だけで判定できる」。★本体は読んでいない。") 6 ]

def chain_cor_3_3.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.3 (還元の合成)",
    sectionId := "cor-3-3" }

def chain_cor_3_3.needs : List ProofObligation := []

def theorem_4_2_of_cor_3_3.src : Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2 (Cor 3.3 と Lemma 4.1 からの還元)",
    sectionId := "theorem-4-2" }

def theorem_4_2_of_cor_3_3.needs : List ProofObligation :=
  [ .implicitStep
      ("原典 §4 本文が両方を名指ししている——「by Corollary 3.3, it follows that V is also a uniformizing」および「generates K … (by Lemma 4.1 below)」。★本体は読んでいない。") 7 ]

def chain_theorem_4_2.src : Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2 (還元の合成)",
    sectionId := "theorem-4-2" }

def chain_theorem_4_2.needs : List ProofObligation := []

def chain_theorem_4_2_both.src : Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2 (2 経路の記録)",
    sectionId := "theorem-4-2" }

def chain_theorem_4_2_both.needs : List ProofObligation := []

def viaBoth.src : Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2 (2 経路の記録)",
    sectionId := "theorem-4-2" }

def viaBoth.needs : List ProofObligation := []

def cor_4_3_of_parts.src : Source :=
  { paper := "GenEll", pdfPage := 22, item := "Corollary 4.3 (Thm 3.8 / Lem 4.1 / Lem 4.2 / Prop 3.4 からの還元)",
    sectionId := "genell-cor-4-3" }

def cor_4_3_of_parts.needs : List ProofObligation :=
  [ .implicitStep
      ("cor_4_3.needs の .otherPaper 4 件をそのまま仮説にした。★本体は原典の証明を読んでいない。") 22 ]

def chain_cor_4_3.src : Source :=
  { paper := "GenEll", pdfPage := 22, item := "Corollary 4.3 (還元の合成)",
    sectionId := "genell-cor-4-3" }

def chain_cor_4_3.needs : List ProofObligation := []

def cor_4_4_of_cor_4_3.src : Source :=
  { paper := "GenEll", pdfPage := 23, item := "Corollary 4.4 (Thm 3.8 / Cor 4.3 / Lem 4.1 からの還元)",
    sectionId := "genell-cor-4-4" }

def cor_4_4_of_cor_4_3.needs : List ProofObligation :=
  [ .implicitStep
      ("原文の証明は 3 行で「Corollary 4.3 とまったく同様(むしろ易しい)」。cor_4_4.needs の .otherPaper 3 件をそのまま仮説にした。★本体は読んでいない。") 23 ]

def chain_cor_4_4.src : Source :=
  { paper := "GenEll", pdfPage := 23, item := "Corollary 4.4 (還元の合成)",
    sectionId := "genell-cor-4-4" }

def chain_cor_4_4.needs : List ProofObligation := []

end ABC3.Skeleton.Goal.Chain
