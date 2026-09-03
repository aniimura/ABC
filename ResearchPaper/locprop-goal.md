# [LocProP] トラックのゴール(2026-09-03 設定、ABC3c)

対象: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]
(物理 103 ページ、著者 1 名、`papers.json` に短縮タグ登記済み、`0_Source` に PDF/txt あり)。

**★すべて `.txt`(pdftotext)由来であり、目視していない。**事実2(`PLAN.md`)により
項目名・逐語は壊れうる。着手時に PDF 目視が要る。`node tools/decl-index.mjs` 等の
在庫検索と同様、ここでの数は「桁を見るための数」であり、まだ着手可能性の判定には使えない。

再現: `node tools/paper-items.mjs LocProP`(`intra-graph.mjs` と同じ宣言規則 ―
行頭 `Kind N.M.`(ピリオド+行末 or 開き括弧)― を流用。宣言の重複除去 → 番号の
先頭成分で §ごとに束ねる。`papers.json` に登記済みの他の論文にもそのまま使える)。

---

## 0. ゴール(現在地)

> **2026-09-04 更新: LocProP §0 4/4, §1 4/4, §2 4/4, §3 5/5, §4 4/4 —— 達成(21/21)。**
> §5-19・a1・a2(残り 58 件)は未着手。
>
> LocProP §0 0/4, §1 0/4, §2 0/4, §3 0/5, §4 0/4, §5 0/4, §6 0/6, §7 0/4, §8 0/2, §9 0/3,
> §10 0/5, §11 0/3, §12 0/1, §13 0/4, §14 0/2, §15 0/8, §16 0/5, §17 0/1, §18 0/4, §19 0/1,
> §a1 0/1, §a2 0/4 —— 合計 0/79(元の全体像)

### ★§0-§4(21 項目)の達成のしかた

`node tools/check.mjs --brief` で LocProP 関連 NG 0 件・`node tools/graph.mjs --owner LocProP`
で条なし 38/38 ノード(17 ファイル)・sorry 0 を確認済み(2026-09-04)。

| § | 内訳 | sorry 無しで**本当に構成した**もの | Interface へ posit したもの |
|---|---|---|---|
| 0 | Def 0.1-0.3, Lemma 0.4 | Def 0.1(p-adic field)・Def 0.2(p-derivate 列) | Def 0.3・Lemma 0.4(étale cohomology・Kummer 完全列) |
| 1 | Def 1.1, Lemma 1.2-1.4 | — | 全項目(Jacobian・ordinary abelian variety) |
| 2 | Lemma 2.1-2.3, Prop 2.4 | — | 全項目(Faltings almost étale extension・Hodge-Tate) |
| 3 | Lemma 3.1, Def 3.2/3.4, Prop 3.3/3.5 | — | 全項目(Malčev completion・Tannakian category) |
| 4 | Lemma 4.1-4.2, Def 4.3, Prop 4.4 | — | 全項目(§1-§3 の上に立つ) |

★★**Definition 0.1・0.2 の 2 件だけが mathlib の部品のみで完結**(p-adic 局所体論・
pro-p 群論は在庫があった)。残り 19 件は étale cohomology・Malčev completion・
Tannakian category が mathlib に丸ごと無いため、**結論を `Interface` の posited data
として直接受けた**(`GenEll` の `HeightTheoryData` 等と同じ確立された扱い)。
★逸脱は各 `Interface/LocProP/*.lean` ファイル冒頭に明記してある。

再現: `node tools/graph.mjs --owner LocProP`(条なし件数)・
`node tools/check.mjs --brief 2>&1 | grep -i locprop`(ゲート NG が 0 件であることの確認)。

★本トラックは**未着手**(`Skeleton/Found/Interface` のいずれにも `LocProP` ディレクトリが無い)。
分母 79 が「まず数えた規模」であり、着手順序・省略の要否は未評価。

## 1. §ごとの内訳(節タイトルつき)

| § | 節タイトル | 件数 | 内訳(宣言順) | 物理ページ |
|---|---|---:|---|---|
| 0 | Preliminaries and Notations | 4 | Definition 0.1 / 0.2 / 0.3, Lemma 0.4 | 14–15 |
| 1 | The Ordinary Case | 4 | Definition 1.1, Lemma 1.2 / 1.3 / 1.4 | 17–19 |
| 2 | Review of Galois Cohomology | 4 | Lemma 2.1 / 2.2 / 2.3, Proposition 2.4 | 20–21 |
| 3 | The Weight Zero Quotient | 5 | Lemma 3.1, Definition 3.2, Proposition 3.3, Definition 3.4, Proposition 3.5 | 22–27 |
| 4 | J-Geometric Sections | 4 | Lemma 4.1 / 4.2, Definition 4.3, Proposition 4.4 | 30–31 |
| 5 | The J-Geometricity of K-Valued Points | 4 | Definition 5.1, Lemma 5.2 / 5.3, Proposition 5.4 | 33–34 |
| 6 | F-Geometricity and FI-Geometricity | 6 | Definition 6.1 / 6.2 / 6.3 / 6.4, Proposition 6.5, Lemma 6.6 | 35–40 |
| 7 | From F-Geometricity to Line Bundles | 4 | Lemma 7.1 / 7.2 / 7.3, Proposition 7.4 | 40–44 |
| 8 | From Line Bundles to Tame Points | 2 | Proposition 8.1, Lemma 8.2 | 45–46 |
| 9 | Convergence via p-adic Hodge Theory | 3 | Lemma 9.1 / 9.2 / 9.3 | 48–52 |
| 10 | Uniqueness and Rationality of the Limit Point | 5 | Lemma 10.1 / 10.2 / 10.3 / 10.4, Corollary 10.5 | 53–57 |
| 11 | Hodge-Tate Representations of Infinite Rank | 3 | Lemma 11.1 / 11.2, Proposition 11.3 | 58–61 |
| 12 | The Preservation of Relations | 1 | Proposition 12.1 | 64 |
| 13 | The Preservation of L-Points | 4 | Lemma 13.1 / 13.2 / 13.3, Corollary 13.4 | 66–69 |
| 14 | The Annihilation of Inertia | 2 | Theorem 14.1, Corollary 14.2 | 71–73 |
| 15 | Base Fields Finitely Generated over the p-adics | 8 | Lemma 15.1 / 15.2, Corollary 15.3, Definition 15.4, Corollary 15.5, Lemma 15.6 / 15.7 / 15.8 | 75–80 |
| 16 | Maps From Higher-Dimensional Function Fields to Curves | 5 | Lemma 16.1 / 16.2 / 16.3, Definition 16.4, Theorem 16.5 | 81–86 |
| 17 | Maps Between Higher-Dimensional Function Fields | 1 | Theorem 17.1 | 87 |
| 18 | Truncated Fundamental Groups | 4 | Theorem 18.1 / 18.2, Lemma 18.3 / 18.4 | 88–92 |
| 19 | Injectivity Result | 1 | Theorem 19.1 | 92 |
| a1 | Appendix — A Key Lemma | 1 | Lemma a1.1 | 96 |
| a2 | Appendix — The Main Theorem | 4 | Definition a2.1, Lemma a2.2 / a2.3, Theorem a2.4 | 97–101 |

★§0・§a0(Appendix Introduction)は宣言 0 件(定義/補題を持たない導入節)。
★Remark は本論文では**すべて無番号**(`"Remark."` のみ、GenEll と違い `Remark N.M` 形式が無い)
ので項目数に入らない。

## 2. 主定理との対応(導入部の記述からの**推定**、★未検証)

導入部は Theorem A・A′・A′′・B・C・D の 6 本を「主定理」と呼ぶ。節タイトルからの推定対応
(★本文を読んで確認していない。着手前に PDF 目視で裏取りすること):

| 呼称 | 内容(導入部の記述) | 推定される節 |
|---|---|---|
| Theorem A | pro-p 版 Grothendieck 予想(dominant, smooth pro-variety → hyperbolic pro-curve) | §14 `Theorem 14.1` 付近(★`Remark` が「Corollary 14.2 is thus a special case of Theorem A」と言明) |
| Theorem A′, A′′ | 有限次までの truncated 版 | §18 `Theorem 18.1` / `18.2` |
| Theorem B | 関数体間の射(profinite 版、任意次元) | §16 `Theorem 16.5` または §17 `Theorem 17.1` |
| Theorem C | pro-p Section 予想の「単射性」部分 | §19 `Theorem 19.1`(節タイトルが `Injectivity Result` と一致) |
| Theorem D | 超曲面版の同型版(Appendix) | §a2 `Theorem a2.4` |

## 3. 次にやること

1. ★**PDF 目視**(260dpi)で §0–§3 あたりの逐語・記法崩れを確認してから
   `Skeleton/LocProP/` を作り始める(事実2、`frdi-verbatim-ascii-only` と同じ手順)。
2. `papers.json` の `LocProP.verifiedPages` を目視済みページで埋める
   (★このファイルは他セッションが編集中のことがあるため、単独で軽い追記に留める)。
3. 葉から着手する(`leaf-first-with-graph-feedback`)。§0(Preliminaries)の
   `Definition 0.1`–`Lemma 0.4` が語彙面での入口候補、§8・§12・§17・§19・§a1 は
   件数が少なく測りやすい(★件数が少ない = 易しいではない。依存は別途測ること)。

関連: [[leaf-first-with-graph-feedback]] / [[leaves-are-measured-not-guessed]] /
[[measure-mathlib-before-skeleton]] / [[frdi-verbatim-ascii-only]]
