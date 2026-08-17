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

## ★★★2026-08-17 深夜: 条ごとに測り直したら、見積りが 2 つ誤っていた

**分子の規則**: `Found/` にある**条なし**の `.src`(`item := "Kind N.M"` 完全一致)。
★★`"Definition 1.5, (i)"` のような**サブ項目も「条つき」**で分子に入らない。

| 命題 | `Found/` に条なしで在る条 | 残る条 | 要るもの |
|---|---|---|---|
| Definition 1.1 | — | 全体 | 射影埋め込み(`ArcModel` の構成) |
| Definition 1.2 | (ii) | (i) | 可逆層のテンソル積 or 移動補題 |
| Example 1.3 | (i) | (ii) | コンパクト領域(非アルキメデス側) |
| Proposition 1.4 | (iv) | (i)(ii)(iii) | (i)(iii) 可逆層、(ii) 射影埋め込み |
| Definition 1.5 | (i)(iii)(iv) | (ii) | **Auslander–Buchsbaum** |
| Proposition 1.6 | — | 全体 | 射影埋め込み(`ArcModel` の構成) |

★★★**誤り 1**: 「`ℙⁿ` の点の関手で 5/24 → 7/24」は**過大**。
`Proposition 1.4` は (i)(ii)(iii)(iv) 全部揃って初めて数に入るので、
`ℙⁿ` が解くのは `Proposition 1.6` の **1 件だけ**(5/24 → **6/24**)。

★★★**誤り 2**: 「Auslander–Buchsbaum 1 本」は**規模の過小評価**。
mathlib の在庫は `RingTheory/RegularLocalRing/Defs.lean` **1 ファイルだけ**で、
定理本体は Serre の特徴づけ + Nagata の補題 + 次元の帰納——**mathlib 級の企画**。
★ただし **`Definition 1.5` はこの 1 本だけで 1 件動く**(位置は良い)。

## ★§3・§4 の残りも測った

- §3 済 = `Lemma 3.1` / `Lemma 3.6`。残り 7 件は **Faltings 高さ・半安定還元・`M_ell`**。
  ★`Lemma 3.5` / `Lemma 3.7` を「補題だから自己完結」と当てたが**外れ**——`htFalt` が要る。
- §4 済 = `Lemma 4.1` / `Remark 4.1.1` / `Lemma 4.2`。残り 2 件は
  `Cor 4.3` / `Cor 4.4` = **論文の最終結果**(IUT 依存)。
  ★★3/5 と見えるのは残りが最終結果だからで、**近いのではない**。

★★★**24 件のうち安い項目は残っていない。**単一障害で 1 件動くのは
`Definition 1.5`(Auslander–Buchsbaum)と `Proposition 1.6`(射影埋め込み)の 2 つだけ。

**Why:** 分母 24 のうち完成は **5 件**。カウンタが動かないのは怠慢ではなく、
`.src` の 2 値規則(条つきは数えない)による。★§1 の実質は大きく進んだが、
1 件に数えるには命題**全体**が原文どおりの仮定で揃う必要がある。

**How to apply:**
- ★最短の再開点は **`Proposition 1.6`**(`ℙⁿ` の点の関手 → `ArcModel` の構成)。
  ★★ただし**動くのは 1 件**である。7/24 と書いた古い見積りを使わないこと。
- ★★`Example 1.3, (i)`(`Found/GenEll/DegSubset.lean`)は済——
  `X(ℚ̄)^{≤d}` / `X(ℚ̄)^{=d}` / `E^{≤d}` / Galois-finite。`Lemma 3.7` が要求するもの。
  ★商 `UPoint` の上で「次数」は代表元ごとに違うので、`≤ d` は
  **「そういう代表元が存在する」**と読む(存在量化なので well-defined が自動)。
- ★★★**「mathlib に無い」と書く前に別の語で検索する**——2026-08-17 に 2 度外した
  ([[frdi-split-nonisotropic-not-derivable]] に記録)。
  ★同日、`ProjectiveSpectrum/Functor.lean` は**ディレクトリを見て**当てた——
  概念名(`functor of points`)では引けなかった。

関連: [[abc3-plan-two-track]] / [[challenger-audit-without-context]] /
[[genell-bd-class-direction]] / [[frdi-split-nonisotropic-not-derivable]]
