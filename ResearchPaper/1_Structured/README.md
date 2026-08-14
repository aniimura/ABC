# 1_Structured — 構造化の規約

`0_Source` の論文を、Lean 形式化のゲート **G1(出典)** がそのまま機械照合できる形へ構造化する。

`0_Source` は今後増えていく。本フォルダも、それに応じて**このプロジェクトの中で作成・拡充を続ける**。以下は固定仕様ではなく、実作業を通じて改良する生きた規約とする(変更したら理由をこのファイルに追記する)。

---

## 0. この工程が答えるべき問い

**「Lean のある宣言は、原典のどの主張を写したものか」を、機械が検証できる形で記録すること。**

読みやすいHTMLを作ることが目的ではない。`1_Structured` は `tools/extract_claims.mjs` の入力であり、そこから出た JSON を `tools/check.mjs` が PDF と Lean の両方に対して照合する:

```
0_Source/*.pdf  (権威)
      │ 目視確認 + 逐語抽出
      ▼
1_Structured/<論文>/section-N.html      ← data-* 属性 + .verbatim
      │ tools/extract_claims.mjs
      ▼
   claims.json
      │ tools/check.mjs が両方向に照合
      ├──▶ 該当 PDF ページの pdftotext 出力に逐語が含まれるか
      └──▶ Lean の宣言が指す locator が実在するか
```

したがって**属性の欠落は「読みにくい」ではなく「検証不能」を意味する**。

---

## 1. なぜ PDF が権威で `.txt` が補助なのか

`.txt`(`pdftotext`)は、望月論文が対象の区別に使っている記法を保持できない。実測した2例:

| | `.txt` | PDF目視 |
|---|---|---|
| [EtTh] p.11 | `Gal(Y /X) (∼= Z)` | `Gal(Y/X)`——Y・X は**下線なし** |
| [EtTh] p.43 Def 2.13(i) | `Gal(Y /X) (∼= l·Z)` | `Gal(Y̲̲/X̲̲)`——Y・X は**下線2本** |

`.txt` 上は同一文字列でありながら、片方は `≅ Z`、もう片方は `≅ l·Z`。テキストだけを見れば論文が自己矛盾しているように見えるが、PDF では別の対象についての別の主張であって矛盾ではない。下線はPDF内でベクター図形として描画されており、**抽出精度を上げても原理的に取れない**。

装飾の脱落は [EtTh] 特有ではない。[pGC] でも `Ẑ`(ハット)→`Z`、`K′`(プライム)→`K `、`𝒪_K`→`OK` の脱落を確認済み。

> **規約**: locator は **物理ページ番号**(`pdftoppm -f N` / `pdftotext -f N` が指すページ)を正とする。`.txt` の行番号は検索補助にのみ使い、locator にしない。

確認コマンド:

```bash
pdftoppm -png -r 150 -f <ページ> -l <ページ> "0_Source/<論文>.pdf" out
```

---

## 2. 手順

1. **登記**: `ResearchPaper/papers.json` にその論文の項目があるか確認する。無ければ先に追加する(タグ・ファイル名・総ページ数・`pageOffset`・`notationRisk`)。
2. **範囲を決める**: 1論文1フォルダ。フォルダ内は原則 `index.html`(全体の目次)+ 節ごとに `section-N.html`。節が長い場合や、特定の定義群だけを扱う場合は `<主題>.html` としてよい。
3. **PDF目視**: 対象範囲の各ページを画像化して見る。`notationRisk: high` の論文では**全ページ必須**。`medium`/`low` でも、数式が現れるページは見る。
4. **逐語抽出**: `.verbatim` に原文をそのまま入れる。記法は §3 のクラスで表現する。
5. **読みを分ける**: 我々の解釈・注釈・Lean 化に向けたメモは `.reading` に入れる。**`.verbatim` に混ぜない**(混ぜると照合が壊れ、かつ「原典が言っていること」と「我々が読み取ったこと」の区別が消える)。
6. **確認日を記録**: `data-notation-checked` に目視確認した日付を入れる。装飾記号が一切現れない項目のみ `none` を許す。
7. **未解決を印す**: 判断がついていない箇所は `<p class="open">` で明示する。空欄で済ませない。

---

## 3. HTML の形

### 構造単位

```html
<section class="statement proposition"
         id="prop-1-1"
         data-paper="pGC"
         data-pdf-page="3"
         data-item="Proposition 1.1"
         data-notation-checked="2026-08-14">
  <div class="verbatim">The cyclotomic character χ : Γ<sub>K</sub> → <span class="bb">Z</span><sub>p</sub><sup>×</sup> can be recovered entirely group-theoretically from Γ<sub>K</sub>.</div>
  <div class="reading">…我々の読み…</div>
</section>
```

`class` の第2語: `definition` / `proposition` / `lemma` / `theorem` / `corollary` / `remark` / `setup`(地の文)。

**必須属性**: `id` / `data-paper` / `data-pdf-page` / `data-item` / `data-notation-checked`。1つでも欠けると `check.mjs` が落とす。

### 証明

```html
<div class="proof" data-proves="prop-1-1" data-pdf-page="3">…</div>
```

### 記法クラス

Unicode の結合文字はフォント依存で不安定なため、**装飾はクラスで表す**(合成済み文字が確実に存在する場合のみ直接使ってよい)。

| クラス | 意味 | CSS |
|---|---|---|
| `ul1` | 下線1本 | `text-decoration: underline` |
| `ul2` | 下線2本 | `text-decoration: underline; text-decoration-style: double` |
| `dot1` / `dot2` | 上の点1個 / 2個 | `::before` で ˙ / ¨ |
| `bar` | 上線 | `text-decoration: overline` |
| `hat` / `tilde` | ハット / チルダ | 疑似要素 |
| `bb` | 黒板太字・`\mathbf`(**Z**, **Q** 等) | `font-weight: bold` |
| `scr` | script体(𝒪 等) | `font-style: italic` |

利点: `grep 'class="ul2"'` でその対象への言及箇所を機械的に洗い出せる。`.txt` では不可能。

### `data-txt` — 照合用の投影を明示する

装飾を復元すると、`.verbatim` の表示テキストは `pdftotext` の出力と一致しなくなる(`K̄` に対して `.txt` は `K`、`K′` に対しては何も出力しない)。そこで**任意の `<span>` に `data-txt` を付けて、その部分の `pdftotext` 出力を明示できる**ものとする。

```html
Fix an algebraic closure <span class="bar" data-txt="K">K</span> of K.
another local p-adic field K<span class="prime" data-txt="">&prime;</span>,
&rarr; <span class="hat" data-txt="Z"><span class="bb">Z</span></span> &rarr; 0
&Gamma;<sub>K</sub> <span class="glyph" data-txt="d=ef">&#8797;</span> Gal(...)
```

`data-txt` が無い場合はその要素のテキスト内容がそのまま使われる。`data-txt=""` は「`.txt` には何も出力されない」を意味する。

**作業手順としては**: (1) `pdftotext -f N -l N` の出力をそのまま `.verbatim` に貼る → (2) PDF 画像を見て装飾を復元し、復元した箇所に `data-txt` で元の出力を残す。この順序で作れば、照合は定義上必ず通る——**通らなければ、転写の途中でどこかを書き換えた証拠**になる。

### 地の文に埋め込まれた暗黙の定義

見出し語を持たないが実質的に定義である箇所は、`<section class="statement definition" data-implicit="true">` として明示的に切り出す。**印を付けるだけで済ませない**——Lean 側は定義を要求するので、ここで切り出さなければ後工程で必ず作り直しになる。

---

## 4. 照合の規則(`check.mjs` が実際に行うこと)

`.verbatim` の**照合射影**——タグを剥がし、`data-txt` を持つ要素はその値に置換し、HTML実体参照を解決し、連続空白を1つに畳んだもの——が、`data-pdf-page` のページの `pdftotext -enc UTF-8 -f N -l N` 出力に(同じ正規化の後で)**部分文字列として含まれること**。

正規化: 合字(ﬁ ﬂ)の展開、改行の空白化、連続空白の畳み込み、両端の空白除去。

検査項目(1構造単位あたり):

| | 検査 | 落ちたときの意味 |
|---|---|---|
| S1 | 必須属性が揃っている | 記録が不完全。G1 を満たさない |
| S2 | `data-paper` が `papers.json` に存在する | 登記されていない論文 |
| S3 | `data-pdf-page` が総ページ数の範囲内 | locator が実在しない |
| S4 | 照合射影が該当ページに含まれる | 転写の誤り、またはページ番号の誤り |
| S5 | `data-notation-checked` が日付、または `none` | 目視確認が未実施 |
| S6 | `id` が文書内で一意 | 参照が曖昧になる |

**S4 が通っても、それは「そのページにその文字列がある」ことしか意味しない。** 装飾の正しさは目視確認(S5 の日付)だけが担保する——この2つは別の検査であり、片方が他方を代替しない。

---

## 5. やらないこと

- **原文の言い換えを `.verbatim` に入れない。** 読みやすさのための整形も不可。
- **他の工程・他のプロジェクトのフォルダを参照しない。** 参照は `0_Source` の PDF と `papers.json` に閉じる。切れたリンクは誤読を生む(実際に生んだ)。
- **数式の完全な再現を目指さない。** 目的は対象の同定であり、組版ではない。同定に関わらない装飾(イタリック等)は落としてよい。

---

## 6. 旧構成ファイルの扱い(`*.legacy.*`)

本フォルダには、この規約が定まる前に**別の作業体系**で作られたファイルが残っている。それらは存在しないフォルダを参照しており(`2_LocatorMap`・`4_CalibrationPlan`・`IUT_DependencyMap`・`MAP-N` 等)、本規約の属性も持たない。

> **規約**: 旧構成のファイルは `<名前>.legacy.<拡張子>` にリネームする。`check.mjs` はこれらを検査対象から**除外**する。

- ファイル名そのものが状態を表すので、開いた人が即座に判別できる。
- 内容は破棄しない——再作成時の下書きとして参照してよい。ただし**そこに書かれた事実(ページ番号・記法・引用)を検証せずに引き継がない**。必ず PDF に当たり直す。
- 再作成が済んだら `.legacy.` 版は削除する(git に残る)。

## 7. 現在の状態

| 論文 | 範囲 | 状態 |
|---|---|---|
| pGC | §1(物理 p.2–4) | **本規約で作成済み**(2026-08-14、目視確認 p.2・p.3・p.4冒頭) |
| pGC | §0・§2・§3・§4 | 未着手 |
| その他8論文 | — | `*.legacy.*` として残置。必要になった時点で本規約で再作成 |

## 8. 改訂履歴

- **2026-08-14 初版**。`PLAN.md` の事実2(`.txt` は対象の同定に使えない)とゲート G1 から逆算して設計した。`data-txt` による照合射影、`papers.json` 登記簿、`*.legacy.*` 規約を導入。最初の適用は pGC §1。
