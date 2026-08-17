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

## ★規模の見積もり

★six exponentials theorem は **Gelfond–Schneider の方法**(補間行列式・
Siegel の補題・Schwarz の補題)で証明される。★mathlib には
`Mathlib/Analysis/Transcendental/` に **Liouville 数**と
**`e` の超越性**しか無く、Gelfond–Schneider すら無い。
★★**mathlib の PR 何本ぶんかでは済まない** —— 解析的整数論の一分野を
入れる作業になる。

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
      "★**これが ① に落ちる条件**: mathlib(または我々)が six exponentials " ++
      "theorem を実装すること。★前提として Gelfond-Schneider の方法" ++
      "(Siegel の補題・補間行列式・Schwarz の補題)が要り、" ++
      "mathlib の `Analysis/Transcendental/` には Liouville 数と e の超越性しか無い" ++
      "(2026-08-17 に検索)。★**③ に上がることはない** —— " ++
      "six exponentials theorem は 1960 年代に確立した定理である。" ++
      "★影響範囲は限定的で、Frobenioid の理論本体ではなく " ++
      "`Example 6.3` の非同型性を言う所にだけ使われる。" }

end ABC3.Gap.FrdI
