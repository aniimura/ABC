---
name: genell-track-b
description: [GenEll] トラック(B)の長期ゴールと、実測で分かった律速要因・臨界パスの到達点。
metadata:
  type: project
---

**北極星(B トラック)**: `[GenEll] Corollary 4.4` まで、IUTchIV が要る範囲の [GenEll] を
`Found/` に `sorry` 無しで載せる。**同じ土台に立つ第二の頂は `Theorem 2.1`**
(「abc ⟺ ℙ¹∖{0,1,∞} 版 abc」。★**IUT を一切使わない**)。
計画本体は `D:\Math_ABC3\ResearchPaper\genell-goal.md`、進捗計数は
`node tools/genell-progress.mjs`(分母 = IUT 側からの需要の推移閉包 **24 件**)。

## ★★★2026-08-17: §1 の律速判定を撤回した

**2026-08-16 に「§1 の律速は Arakelov 理論(算術直線束)の不在」と測ったのは誤りだった。**
★★実際に要るのは `X^arc` の**位相とコンパクト性だけ**で、
正則構造・連接層・GAGA はどこにも現れない。

同日 33 ファイルを実装し(すべて sorry 0・標準 3 公理のみ)、**高さそのものを構成した**:

| 到達点 | 場所 |
|---|---|
| `X^arc` の位相とコンパクト性(GAGA 不要) | `ProjTopology` / `ProjClosed` / `ArcModel` |
| ℂ-点と Green 関数(解析空間不要) | `ArchPoint` |
| 高さ `htU : U_X(ℚ̄) → ℝ` | `HeightConstruction` … `UPoint` |
| `Prop 1.4, (i)(ii)(iii)` が**無条件** | `HeightAdditive` / `HeightNonneg` / `HeightClass` |
| `Prop 1.6` が組み上がった | `Prop16` |

★★`Prop 1.4` は posit された `HeightTheoryData` の上では**偽**である
(反例は `Check/GenEll/HeightAxiomGap.lean`。`Remark 1.4.1` / `1.5.1` も同様)。
**構成に置き換えると真になる**——しかも (ii)(iii) は原文より強い形
(`≳` でなく `≥`、BD-class でなく等式)で出る。

## ★★§1 に残るのは mathlib の在庫 3 本だけ

| 残り | 効く先 |
|---|---|
| `ℙⁿ` の**点の関手**(`ℙⁿ(ℂ) ≅ ℙ ℂ V`) | `Prop 1.4 (ii)` + `Prop 1.6` の 2 件 |
| **Auslander–Buchsbaum**(正則局所 ⇒ UFD) | `Def 1.5 (ii)` |
| 可逆層のテンソル積 or 移動補題 | `Def 1.2 (i)` の全域化 |

★★★**どれも「原文の数学」ではなく「mathlib の在庫」の問題である。**

**Why:** 分母 24 のうち完成は **5 件**。カウンタが動かないのは怠慢ではなく、
`.src` の 2 値規則(条つきは数えない)による。★本日 §1 の実質は大きく進んだが、
1 件に数えるには命題**全体**が原文どおりの仮定で揃う必要がある。

**How to apply:**
- ★最短の再開点は **`Proposition 1.6`**。残るのは `ℙⁿ` の点の関手 1 本で、
  それが入れば `Prop 1.4 (ii)` と合わせて **5/24 → 7/24** も見える。
- ★★§3(Galois 作用)は未着手のまま。`Lemma 3.1 (iv)` は 9 段中 6 段が済。
- ★★★**「mathlib に無い」と書く前に別の語で検索する**——本日 2 度外した
  ([[frdi-split-nonisotropic-not-derivable]] に記録)。

関連: [[abc3-plan-two-track]] / [[challenger-audit-without-context]] /
[[genell-bd-class-direction]] / [[frdi-split-nonisotropic-not-derivable]]
