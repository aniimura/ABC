import ABC3.Found.FrdI.Lemma65

/-!
# Gap — [FrdI] §6 で mathlib に在庫が無い 1 段

★**`Found/FrdI/Lemma65.lean` で `Lemma 6.5, (i)` を証明した。**
★**(ii) は超越数論の定理を要し、mathlib に在庫が無い。**

## ★原文

原文 (FrdI p.117):
> Assertion (i) is a formal consequence of the fact that Z is a unique factor-

★(i) は「`ℤ` が一意分解環であること」の形式的帰結 —— **実装済み**
(`log_primes_linearIndependent`)。

★★(ii) は原文が **Lang の定理**([Lang1]、Baker の本 p.119)に送っている。
★これは **six exponentials theorem**(6 指数定理)である:

> `x₁, x₂` が ℚ 上一次独立、`y₁, y₂, y₃` が ℚ 上一次独立な複素数なら、
> `exp(xᵢyⱼ)`(6 個)の少なくとも 1 つは超越数である。

★原文の (ii) はこれを `{log p₂, log p₄, log p₆}` と `{1, log p₃/log p₄}` に当てる。

## ★★測定(2026-08-17)

★mathlib を検索した:

| 検索語 | 件数 |
|---|---|
| `sixExponentials` / `six_exponentials` | **0** |
| `Baker`(超越論の文脈) | **0** |
| `Real.log` の素数についての一次独立性 | **0**(我々が最初) |

★★**したがって (ii) は「原典の飛躍」ではなく「在庫の不足」**である。
★分類は **② `missingMath`** —— `Definition 2.8, (ii)` の pro-`l` 分解と同じ形である。

## ★★★2026-08-18 の訂正 —— **上の規模の見積もりは誤りだった**

★以下は 2026-08-17 に書いた見積もりであり、**取り消す**:

> ★mathlib には `Mathlib/Analysis/Transcendental/` に **Liouville 数**と
> **`e` の超越性**しか無く、Gelfond–Schneider すら無い。
> ★★**mathlib の PR 何本ぶんかでは済まない** —— 解析的整数論の一分野を
> 入れる作業になる。

★★**mathlib に `Analysis/Transcendental/` というディレクトリは無い。**
実体は **`NumberTheory/Transcendental/`** である。**存在しないパスを探して 0 件を得ていた。**
★`Meta/Claim.lean` の `LeanStatus.absent` が `searched` を要求するのは、まさに
この型の誤りを捕まえるためだった —— ★**`searched` に書いたパスが実在するかを
確かめる段が抜けていた。**

★★改めて測った在庫(2026-08-18):

| 要るもの | mathlib | 判定 |
|---|---|---|
| **数体上の Siegel の補題**(house 版) | `NumberField.house.exists_ne_zero_int_vec_house_le`(`NumberTheory/NumberField/House.lean`) | ★★**ある**(算術側の心臓部) |
| Siegel の補題(整数版) | `Int.Matrix.exists_ne_zero_int_vec_norm_le`(`NumberTheory/SiegelsLemma.lean:181`) | ★**ある** |
| 代数的数の house とその評価 | `NumberField.house` / `house_mul_le` / `norm_embedding_le_house` | ★**ある** |
| 最大値原理・Schwarz(1 点・高位)・Cauchy 評価・Hadamard 三線・Jensen | `Analysis/Complex/` | ★**ある** |
| **多零点版の小ささ評価** | `Analysis/Complex/` を目視して該当なし | ★**無い** |
| Schneider–Lang / six exponentials | `NumberTheory/Transcendental/` は `Liouville/` と `Lindemann/AnalyticalPart.lean` のみ | ★**無い** |

★★**「解析的整数論の一分野を入れる作業」ではない。**
statement は `Skeleton/FrdI/Lemma65SixExp.lean` に型で固定し、
道は `ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン
(8 節点・6 層・葉 3 個)に割った。`node tools/frdi-newleaves.mjs` が層と葉を印字する。

★★★**ただし [FrdI] の本筋への影響は限定的**である —— `Lemma 6.5, (ii)` は
`§6` の**例**(`Example 6.3` の非同型性)を言うために使われるのであって、
Frobenioid の理論そのものには入らない。
-/

namespace ABC3.Gap.FrdI

/-- ★★**`Lemma 6.5, (ii)` に不足しているもの**。

★分類は **② `missingMath`** —— 原文の主張は標準的な(しかし深い)数学であり真である。
**足りないのは mathlib の在庫**であって、原典の飛躍ではない。 -/
structure Gap_6_5_ii : Prop where
  /-- 不足: **six exponentials theorem**(Lang / Ramachandra)。 -/
  sixExponentials : True

def Gap_6_5_ii.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 116, item := "Lemma 6.5, (ii)",
        sectionId := "frdi-lemma-6-5" },
    classification := ABC3.Meta.GapClass.missingMath,
    falsifier :=
      "★★★★★2026-08-20: **解消した** —— " ++
      "`Found/SixExp/Theorem.lean` で six exponentials theorem を**証明した**(sorry なし)。" ++
      "`Found/FrdI/Lemma65ii.lean` の `lemma_6_5_ii` が `Lemma 6.5, (ii)` そのものである。" ++
      "★この記録は「何が足りなかったか」の履歴として残す。 " ++
      "★**これが ① に落ちる条件**: mathlib(または我々)が six exponentials " ++
      "theorem を実装すること。★前提として Gelfond-Schneider の方法" ++
      "(Siegel の補題・補間行列式・Schwarz の補題)が要る。" ++
      "★★**2026-08-18 の訂正**: 前日ここに書いた「mathlib の " ++
      "`Analysis/Transcendental/` には Liouville 数と e の超越性しか無い」は" ++
      "**探索先のパスが存在しなかった**ための誤りである(実体は " ++
      "`NumberTheory/Transcendental/`)。★改めて測ると、" ++
      "**数体上の Siegel の補題**(`NumberField.house.exists_ne_zero_int_vec_house_le`)も" ++
      "整数版(`Int.Matrix.exists_ne_zero_int_vec_norm_le`)も**在庫にある**。" ++
      "残るのは多零点版の小ささ評価と外挿の帰納であり、" ++
      "道は `ResearchPaper/frdi-decomposition.json` の `sixexp` チェーンに割った。" ++
      "★**③ に上がることはない** —— " ++
      "six exponentials theorem は 1960 年代に確立した定理である。" ++
      "★影響範囲は限定的で、Frobenioid の理論本体ではなく " ++
      "`Example 6.3` の非同型性を言う所にだけ使われる。" }

end ABC3.Gap.FrdI
