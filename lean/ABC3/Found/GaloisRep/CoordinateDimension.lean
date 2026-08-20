import ABC3.Found.GaloisRep.TorsionIdealIntegral
import Mathlib.RingTheory.DedekindDomain.Basic

/-!
# Galois (G5) 第 127 ブロック —— **★★★座標環は Krull 次元 ≤ 1**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★Dedekind への 2 条件のうち 1 つ

§9-434 の再実測で、mathlib は **Dedekind 環の因子機構を完備している**ことが分かった。
★欠けているのは `IsDedekindDomain W.CoordinateRing` の instance であり、
その内訳は:

| 条件 | 状態 |
|---|---|
| `IsDomain W.CoordinateRing` | ✅ mathlib |
| `IsNoetherianRing W.CoordinateRing` | ✅ mathlib |
| **`Ring.DimensionLEOne W.CoordinateRing`** | ✅ **本ブロック** |
| `IsIntegrallyClosed W.CoordinateRing` | ❌ 残り |

★★次元の側は**第 116 ブロックの整拡大**がそのまま効く:

    F[X] は PID(次元 ≤ 1)、F[W] は F[X] 上整
      ⟹ F[W] も次元 ≤ 1        (mathlib の `DimensionLEOne.of_isIntegral`)

## ★★残るのは整閉性だけ

`IsIntegrallyClosed` は**アフィン曲線の正則性**である。
★mathlib には局所判定 `IsIntegrallyClosed.of_localization_maximal` があるので、
各極大イデアルでの局所化が DVR であることを示す形になる
(非特異点なら `x − x₀` か `y − y₀` が極大イデアルを生成する)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `coordinateRing_dimensionLEOne` | ★★★**座標環は Krull 次元 ≤ 1** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★**座標環は Krull 次元 ≤ 1 である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`F[X]` は PID(次元 ≤ 1)で、`F[W]` はその上整(第 116 ブロック)なので、
mathlib の `DimensionLEOne.of_isIntegral` がそのまま効く。 -/
instance coordinateRing_dimensionLEOne (W : WeierstrassCurve.Affine F) :
    Ring.DimensionLEOne W.CoordinateRing :=
  Ring.DimensionLEOne.of_isIntegral (Polynomial F) W.CoordinateRing

/-! ## ★出典の紐付け(`.src`) -/

def coordinateRing_dimensionLEOne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——座標環の Krull 次元が 1 以下であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
