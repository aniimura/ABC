import ABC3.Meta.Claim
import Mathlib.Data.Set.Finite.Basic

/-!
# [NCBelyi] Theorem 2.5 の語彙(`Interface`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5–p.6。
**260 dpi 目視確認 2026-08-17**。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★★なぜこの論文が要るのか

`[GenEll] Theorem 2.1`(**IUT を経由しない第 2 の頂上**)の証明が
**唯一引く外部文献**がこれである(`[GenEll]` 物理 p.11)。

★★**そして原典は `0_Source` にある**(9 ページ)。
mathlib に `Belyi` は 0 件だが、**原典は手元にある**——
ゆえにこの依存は「不在」ではなく「未転写」である。

## ★posit するのは語彙だけ

★`Theorem 2.5` の結論(射 `φ` の存在)は**入れない**。
本構造体は**公理を 1 つも持たない**(型と述語だけ)ので、
`Skeleton/NCBelyi/Theorem25.lean` の `sorry` は本当に残る。

## ★mathlib 実測(2026-08-17)

| statement の語 | mathlib |
|---|---|
| `Belyi` | ★**0 件**(2026-08-16、全体を case-insensitive で検索) |
| `ℚ̄` 上の滑らかな固有連結曲線 | `AlgebraicGeometry` はあるが「曲線」の型は無い |
| **Serre 双対性**(証明が使う) | ★**未測定**——「まだ測っていない」であって「無い」ではない |
-/

namespace ABC3.Interface.NCBelyi

open ABC3.Meta

/-- **[NCBelyi] Theorem 2.5** の statement を書くのに要る語彙。

★フィールドはすべて原文 p.5–p.6 に**実際に現れる語**と 1 対 1 に対応する。 -/
structure BelyiSetup where
  /-- 原文「a smooth, proper, connected curve over `ℚ̄`」。 -/
  Curve : Type
  /-- `X(ℚ̄)` —— `ℚ̄`-有理点。 -/
  Point : Curve → Type
  /-- `ℙ¹_ℚ̄` の `ℚ̄`-有理点。 -/
  ProjPoint : Type
  /-- 原文の 3 点 `{0, 1, ∞}`。 -/
  three : Set ProjPoint
  /-- 射 `φ : X → ℙ¹_ℚ̄`。 -/
  Mor : Curve → Type
  /-- 射が点に及ぼす作用 `φ(−)`。 -/
  app : {X : Curve} → Mor X → Point X → ProjPoint
  /-- 原文 (a)「`φ` is **unramified** over `ℙ¹_ℚ̄ \ {0,1,∞}`」。 -/
  UnramifiedOutsideThree : {X : Curve} → Mor X → Prop
  /-- 原文「a number field `F`」。 -/
  NumField : Type
  /-- 原文「`X` … defined over a number field `F`」。 -/
  CurveOver : Curve → NumField → Prop
  /-- 原文「`S`, and `T` are defined over a number field `F`」。 -/
  PointsOver : {X : Curve} → Set (Point X) → NumField → Prop
  /-- 原文「`φ` may be taken to be **defined over** `F`」。 -/
  MorOver : {X : Curve} → Mor X → NumField → Prop

/-- ★Track B は何を作らねばならないか。 -/
def BelyiSetup.waiting : WaitingFor :=
  { what := "ℚ̄ 上の滑らかな固有連結曲線と ℙ¹_ℚ̄ への射、分岐の概念、および数体上の定義体。★原文の証明はさらに **Serre 双対性**・曲線上の直線束の次数・線型系を使う(物理 p.6 目視)"
    trackB := "Found/NCBelyi — ★mathlib に `Belyi` は 0 件(2026-08-16 実測)。ただし**原典 9 ページは `0_Source` にある**ので、これは『不在』ではなく『未転写』である。★Serre 双対性が mathlib にあるかは**未測定**——測ってから判断すること" }

/-! ## ★出典の紐付け(`.src`) -/

def BelyiSetup.src : Source :=
  { paper := "NCBelyi", pdfPage := 5, item := "Theorem 2.5",
    sectionId := "ncbelyi-thm-2-5" }

end ABC3.Interface.NCBelyi
