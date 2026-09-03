# [CorrHyp] トラックのゴール(2026-09-04 設定、ABC3b)

対象: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]
(物理 18 ページ、著者 1 名、`papers.json` に短縮タグ登記済み、`0_Source` に PDF/txt あり)。

**★2026-09-04、p.3-13・p.15・p.17 を 200dpi 目視確認済み**(`papers.json`
`notationRisk: "medium"`)。実害は overline(X̄ → X)と script 文字(𝒳・𝒴 → X・Y)の
脱落の 2 点(GenEll のような行列の順序入れ替えは無い、詳細は `papers.json` を見よ)。

再現: `node tools/paper-items.mjs CorrHyp`(ABC3c が LocProP のために新設した汎用ツール。
行頭 `Kind N.M.`(ピリオド+行末 or 開き括弧)の宣言規則で番号付き項目を §ごとに数える)。
24 件・重複なしを確認済み(手作業での逐次読解と突き合わせて一致)。

---

## 0. ゴール(現在地)

> **CorrHyp §1 5/5, §2 6/6, §3 3/3, §4 2/2, §5 7/7, §6 1/1 —— 合計 24/24(Skeleton)**

★★**2026-09-04、S0(Skeleton = statement を型で固定する段)を完了した。**
`Interface/CorrHyp/HyperbolicCurve.lean`(`HyperbolicCurveData`・`StackType` を posit)
と `Skeleton/CorrHyp/Section1.lean`〜`Section6.lean` に 24 項目すべてを置き、
`lake build` 0 エラー・`tools/check.mjs` で G1(出典・逐語照合)を全項目パスさせた。
残る `sorry` は 19 件(Definition 5 件は sorry 無しで構成済み)。G9(非空虚性の対照)は
14 件が未着手のまま——プロジェクト全体でトラッキング対象の既知 debt であり、
`check.mjs` 自身が「新規は落とす(ブロックしない)」と明記している。
次の段(S1: Track B で `Found/CorrHyp/` を積む)はまだ手つかず。

## 1. §ごとの内訳(節タイトルつき)

| § | 節タイトル | 件数 | 内訳(宣言順) | 物理ページ |
|---|---|---:|---|---|
| 1 | Basic Definitions | 5 | Definition 1.1 / 1.2 / 1.3 / 1.4 / 1.5 | 3–4 |
| 2 | Review of Results of Margulis and Takeuchi | 6 | Definition 2.1 / 2.2 / 2.3, Proposition 2.4, Theorem 2.5 / 2.6 | 5–7 |
| 3 | The Non-arithmetic Case | 3 | Definition 3.1, Proposition 3.2, Theorem 3.3 | 8 |
| 4 | The Main Theorem | 2 | Lemma 4.1, Theorem 4.2 | 9–10 |
| 5 | Isogenies of General Curves | 7 | Lemma 5.1, Definition 5.2, Theorem 5.3, Lemma 5.4 / 5.5 / 5.6, Theorem 5.7 | 11–15 |
| 6 | Interpretation of a Theorem of Royden | 1 | Theorem 6.1 | 17 |

★`Corollary` は本論文に 0 件。`Remark` は本論文では**すべて無番号**(`"Remark."` のみ、
`Remark N.M` 形式が無い)ので LocProP と同様に項目数に入らない。

★§0(Introduction)は分母 24 から**除外した**(LocProP の §0 とは扱いが異なる——
LocProP の §0 は新規の Definition/Lemma を持つが、本論文の §0 は `Theorem A`/`B`/`C` の
3 件しかなく、そのいずれも本文中で明示的に後続定理の再掲と言明されている。二重計上を
避けるための除外であり、GenEll の §1 起点と同じ扱い):

| §0 の呼称 | 本文中の再掲先の明言(逐語) | ページ |
|---|---|---|
| Theorem A | 「cf. Lemma 4.1 and Theorem 4.2 in the text」 | p.1 |
| Theorem B | 「Theorem B follows from Theorem 5.3 in the text」 | p.2 |
| Theorem C | 「(given as Theorem 6.1 in the text)」 | p.2 |

## 2. 主定理との対応(★本文の言明そのもの、LocProP と違い推定ではない)

| 呼称 | 内容(導入部の記述) | 対応する節 |
|---|---|---|
| Theorem A | 双曲曲線に isogenous な曲線の有限性(第一主定理) | `Lemma 4.1` + `Theorem 4.2`(§4) |
| Theorem B | 一般の曲線は自身の hyperbolic core に一致——非自明な相関を持たない | `Theorem 5.3`(§5) |
| Theorem C | `M_{g,r}` は非自明な自己同型・相関を持たない(Royden の定理の帰結) | `Theorem 6.1`(§6) |

§2(Margulis・Takeuchi の結果の review、`Theorem 2.5`/`2.6`)は他論文 `[Marg]`/`[Take]` の
引用結果であり、本プロジェクトでは証明せず posit する対象になる可能性が高い(要検討——
`Interface` 行き。`genell-track-b` で `Theorem 3.8` 等の外部結果を扱った前例を参照)。

## 3. Track B(Found/CorrHyp/)——2026-09-04 に着手、実測に基づく規模感

★★**Found(sorry 無しの証明)まで24件全部を1セッションで終える見込みは無い**
——GenEll(数か月・現在も§3-§4に3本のsorryが残る)と同型の規模になることを、
実際に mathlib を調べたうえで確認した。「壁」ではなく道として、葉から実装する。

### 使える資産(mathlib、2026-09-04 実測)

| 材料 | 場所 |
|---|---|
| `ℍ`(上半平面)・`SL(2,ℝ)` の Möbius 作用 | `Analysis.Complex.UpperHalfPlane.*` |
| `Subgroup 𝒢 ≤ SL(2,ℝ)`(離散)→ `ℍ` に固有不連続作用 | `UpperHalfPlane.instProperlyDiscontinuousSL2RSubgroup` |
| Schreier の補題(有限指数部分群は有限生成) | `Subgroup.fg_of_index_ne_zero` |
| 被覆空間・モノドロミー | `Topology.Covering.*` |
| `PSL(n,K)` の構成 | `Matrix.ProjectiveSpecialLinearGroup` |

### 無い(§2・§6 の律速)

Margulis/Shimura-arithmeticity(代数群・非可換 Galois コホモロジー・Brauer 群)、
四元数環の分類理論、モジュライスタック `M_{g,r}`、Teichmüller 空間、
Royden の定理——いずれも GenEll のモジュラー多項式ギャップと同型の
「未構築の数学」で、§2(`Definition 2.2`/`2.3`・`Proposition 2.4`・`Theorem 2.5`-`2.6`)
と §6(`Theorem 6.1`)はこれ単体で人年規模になりうる。

### 実装した分(`Found/CorrHyp/FuchsianGroup.lean`、すべて sorry 無し・標準3公理のみ)

★2026-09-04、`Space := FuchsianGroup`(`SL(2,ℝ)` の離散部分群、原文の
`PSL₂(ℝ)⁰` からの逸脱を記録)というモデルで §1・§3 の骨格を実装した:

| 節点 | 実装した宣言 | 状態 |
|---|---|---|
| §1 `Definition 1.1`-`1.4` の土台 | `FEt := IsFiniteIndexIn`・`comp`・`pullback`(`inter`)・`pbFst`・`pbSnd`・`idFEt` | ✅ 圏構造が揃った |
| §1 `Definition 1.5` 直後(isogeny は同値関係) | `isIsogenous_refl`/`_symm`/`_trans` | ✅ **完成** |
| §2 `Definition 2.1` 直前(`Γ ⊆ Comm(Γ)`) | `self_le_commensurator`(mathlib の `Subgroup.Commensurable.commensurator` を発見・流用) | ✅ **完成** |
| §3 `Proposition 3.2`(`Γ_Z ⊆ Comm(Γ_X)`) | `prop_3_2`(`Ncore` = `Z.Γ` の中での `C.Γ` の normal core、を経由) | ✅ **完成**——numbered item の本体として初 |
| §3 `Theorem 3.3` | `fg_of_finiteIndexIn`(Schreier の補題、`Theorem 3.3` が要る片方向) | 部分——下記参照 |

### `Theorem 3.3` に残る 2 つの壁(実測、まだ未着手)

1. **群論側**: 原文は「`Comm(Γ_X)` は有限生成部分群 `Γ_X` を有限指数に持つので、
   それ自身も有限生成」と1行で済ませているが、これは
   `fg_of_finiteIndexIn`(有限指数**部分群**が有限生成、Schreier)の**逆向き**
   ——「有限生成な有限指数部分群を持つ群は有限生成」——であり、`exact?` で
   mathlib に見当たらなかった(2026-09-04 実測)。Reidemeister–Schreier 型の
   議論(`H` の生成元 + `G/H` の代表系)を自分で組む必要がある。
2. **幾何側(§3 単体では閉じない)**: 「次数が `g',r'` と `e_Y` だけで抑えられる」
   という原文の主張は、`FuchsianGroup` を実際の種数・穴の数(§1 冒頭の `type`)
   に結び付ける Riemann–Hurwitz 型の計算を要求する——これは §5 の `StackType`/
   `e_Y` の装置を先に作ることに等しく、`Theorem 3.3` 単体としては閉じない
   (`Theorem 4.2` も同型の壁を持つ)。

### 次の一手(依存の少なさ順)

1. 上記1(f.g. の逆向き伝播)を単独の補題として仕上げる——`Theorem 3.3`
   本体はまだ届かないが、Reidemeister–Schreier の実装は他の節(`Lemma 5.5`)
   にも要る資産になる。
2. §5 の `StackType`/`e_Y`(Euler 標数)を `FuchsianGroup` の商 `ℍ/Γ` の
   幾何と結び付ける——これができて初めて `Theorem 3.3`/`Theorem 5.3`/
   `Lemma 5.4`-`5.7` に届く。
3. §4(`Lemma 4.1`・`Theorem 4.2`)は係数拡大 `Ext`(k → K)が要り、
   スキーム論(有限型・降下理論)が要るためさらに重い。
4. §2 の `MargulisArithmetic`/`ShimuraArithmetic` 本体と §6 の Royden の定理は
   上記のとおり大きな未構築の数学(代数群・Teichmüller 空間)を要する。
5. G9(非空虚性の対照)は Track B が本物のデータを供給して初めて
   意味のある形で作れる(`.waiting` に記録済み)。

関連: [[leaf-first-with-graph-feedback]] / [[leaves-are-measured-not-guessed]] /
[[measure-mathlib-before-skeleton]] / [[genell-track-b]] / [[corrhyp-track-goal]] /
[[stale-status-read-lean-first]] / [[no-wall-decompose-instead]]
