import ABC3.Interface.GenEll.HeightTheory

/-!
# [GenEll] §2 —— Theorem 2.1 の設定を受ける `Interface`

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.11。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★なぜ `Interface` が要るのか

`Theorem 2.1` の statement には、`HeightTheoryData` の語彙に無いものが 3 つ現れる:

| statement の語 | mathlib(2026-08-16 実測) |
|---|---|
| 数体上の**滑らかな固有幾何的連結曲線** | `AlgebraicGeometry` はあるが「数体上の曲線」の型は無い |
| **標準層** `ω_X` と `ω_X(D)` | ★`LineBundle` が **0 件**なので `ω_X(D)` も無い |
| **双曲的**(`deg ω_X(D) > 0`) | 直線束の次数が無い |

★ここでも **posit するのは語彙だけ**である。
`Theorem 2.1` の結論(2 つの予想の同値性)は一切入れていない。

## ★★この定理は「IUT を使わない」が「容易」ではない

★`Theorem 2.1` は `[IUTchIV]` の結果を**一切使わない**——§1 の枠組みと
**noncritical Belyi maps**([Mzk1])だけで証明される。
★しかし mathlib に `Belyi` は **0 件**であり、`Arakelov` も **0 件**である。
**IUT 非依存であることは、形式化が容易であることを意味しない。**
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta

/-- **`Theorem 2.1` の設定**を受ける `Interface`。

★フィールドはすべて原文 p.11 に**実際に現れる語**と 1 対 1 に対応する。 -/
structure AbcSetup where
  /-- 原文「Let `X` be a smooth, proper, geometrically connected curve over a number field」。 -/
  Curve : Type
  /-- 各曲線が定める §1 の高さ理論。 -/
  data : Curve → HeightTheoryData
  /-- 原文「`D ⊆ X` a **reduced** divisor」。 -/
  Reduced : (X : Curve) → (data X).Divisor → Prop
  /-- 原文「Write `ω_X` for the canonical sheaf on `X`」から作られる `ω_X(D)`。 -/
  omegaOf : (X : Curve) → (data X).Divisor → (data X).ABundle
  /-- 原文「`U_X` is a **hyperbolic curve** — i.e., that the degree of the line bundle
  `ω_X(D)` is **positive**」。 -/
  Hyperbolic : (X : Curve) → (data X).Divisor → Prop
  /-- 原文「`P ≝ ℙ¹_ℚ` be the **projective line** over `ℚ`」。 -/
  projLine : Curve
  /-- 原文「`C ⊆ P` the divisor consisting of the three points "0", "1", and "∞"」。 -/
  threePoints : (data projLine).Divisor
  /-- 原文「a **compactly bounded subset** … whose **support contains** `Σ`」。
  ★`Example 1.3, (ii)` の support は「compactly bounded subset 自身から完全に決まる」と
  原文 p.6 が述べているので、集合と素数集合の 2 引数述語として持てば足りる。 -/
  SupportContains : Set (data projLine).Point → Finset ℕ → Prop

/-- ★Track B は何を作らねばならないか。 -/
def AbcSetup.waiting : WaitingFor :=
  { what := "数体上の滑らかな固有幾何的連結曲線、標準層 ω_X と ω_X(D)、直線束の次数と双曲性、および ℙ¹ 上の 3 点因子 {0,1,∞}。★加えて証明には noncritical Belyi maps([Mzk1])が要る"
    trackB := "Found/GenEll — ★mathlib に `Belyi` は 0 件、`LineBundle` は 0 件(2026-08-16 実測)。★ただし本定理は [IUTchIV] を一切使わないので、**IUT 側の進捗と独立に進められる**——GenEll の 2 つ目の頂上である" }

/-! ## ★出典の紐付け(`.src`) -/

def AbcSetup.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1",
    sectionId := "genell-thm-2-1" }

end ABC3.Interface.GenEll
