import ABC3.Found.PGC.Prop22FilteredHypothesis
import ABC3.Found.PGC.QpResidueField

/-!
# [pGC] Theorem 4.2 の現行形（濾過を実物に固定した形）で、退化の逃げ道が塞がったこと

`Check/PGC/Theorem42NaiveGC.lean` は、`RF` を全称量化した旧形について

* `topFiltration p`（`Gv ≡ ⊤`）が `IsNaturalFiltration` を満たす、
* ゆえに `RF := topFiltration p` を代入でき、
* 代入すると `OutFilt` は「ただの連続外部同型の全体」になり、
* 全射性がそのまま「素朴 Grothendieck 予想」になる（原典が偽と述べている命題）

ことを `sorry` 無しで示した。

2026-09-08 に `Skeleton/PGC/Section4.lean::theorem_4_2` は `RF` を実物
`ABC3.Found.PGC.ramificationFiltration p` に固定する形へ直された。本ファイルは

  その代入が現行形では実行できないこと

を機械検査する。すなわち `ramificationFiltration p` は `topFiltration p` ではない。

## 何を示していて、何を示していないか（測ってあることだけ書く）

* 示したこと: `ramificationFiltration p ≠ topFiltration p`。
  根拠は `Found/PGC/Prop22FilteredHypothesis.lean::ramificationFiltration_Gv_zero_ne_top`
  （仮定ゼロ）——実物は `v = 0` で惰性群であり、`⊤` ではない。
  ゆえに `Check/PGC/Theorem42NaiveGC.lean` の代入 `RF := topFiltration p` は
  現行形の statement には当てはまらない。
* 示していないこと: 現行形が真であること。
  とくに `v > 0` の範囲で `Γ_K^v` が `⊤` でないことは測っていない
  （本ファイルが使うのは `v = 0` の情報だけである）。
  上付き番号付けの高次分岐群が `v > 0` でどう減るかは別の節点の仕事。

## 原文（該当箇所、逐語）

原文 (pGC 物理 p.1, Introduction):

> On the one hand, one knows (cf. the Remark in [4] following Theorem 4.2) that the
> Grothendieck Conjecture cannot hold in the naive sense (i.e., if one removes the condition
> of "compatibility with the filtrations" from the outer isomorphisms considered - see, e.g.,
> [8]), so one must put some sort of condition on the outer isomorphisms of Galois groups that
> one considers.

現行形は濾過を実物に固定しているので、この「条件を外した版」を特殊化として含まない。

これは原典の主張ではなく我々のモデルについての事実なので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Found.PGC

/-- 実物の高次分岐濾過は退化濾過ではない。

`Check/PGC/Theorem42NaiveGC.lean` が使った代入 `RF := topFiltration p` は、
`RF` を実物に固定した現行形の `theorem_4_2` には当てはまらない。 -/
theorem ramificationFiltration_ne_topFiltration (p : ℕ) [Fact p.Prime] :
    ramificationFiltration p ≠ topFiltration p := by
  intro h
  refine ramificationFiltration_Gv_zero_ne_top (selfField p) ?_
  rw [h]
  rfl

#print axioms ramificationFiltration_ne_topFiltration

end ABC3.Check.PGC
