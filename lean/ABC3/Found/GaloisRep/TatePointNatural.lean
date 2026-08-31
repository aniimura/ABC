/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateCoordNatural
import ABC3.Found.GaloisRep.TatePt

/-!
# Tate 一意化の**点は自然**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★段 2c が閉じる

`galois-equivariant-tate-uniformization` の段 2c は
「`tatePtPair` の自然性（級数と `σ` の交換）」であり、同変性の**本体**であった。

| 段 | 内容 | ファイル |
|---|---|---|
| 2c-1 | `evalAdic` の自然性 | `AdicEvalNatural.lean` |
| 2c-2 | `adicSum` 〜 `tateXK` の自然性 | `TateSeriesNatural` / `TateCoordNatural` |
| **2c-3** | **`tatePtPair` の自然性** | ★本ファイル |

★★到達点:

> `tatePtPair (σa) (σw) q … = Point.some (σ_K X) (σ_K Y) …`

すなわち **`σ` を施したデータの点は、元の点の座標に `σ` を施したもの**である。

## ★★★★★★母数を固定する場合に絞る

`Point.some` の合同を取るには**同じ曲線の上**でなければならない。
★`tateCurveAt q` は `q` で決まるので、`σ q = q` に絞る
（`G_K` は `q ∈ K` を固定するので、これが実際の場合である）。

★★`σ q = q` から `tateXK … (σ q) … = tateXK … q …` へ移るのは
**証明の非関与**（`congr 1`）だけである。

## ★残っている段（明示）

★段 3（有限拡大 `K ⊆ L` に対する自然性）と段 4（`K̄` への帰納極限）が残る。
★★どちらも**新しい数学ではなく組み立て**である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★母数を固定した形の座標の自然性 -/

/-- ★★**`σ q = q` の場合の `X` 座標の自然性**（同じ曲線の上）。 -/
theorem tateXK_map_fixed {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [Algebra R K]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σ r))
    (a w q : R) (hq : q ∈ I) (hσq : σ q = q) (hw : w ∈ I) :
    σK (tateXK a w q hq : K) = tateXK (σ a) (σ w) q hq := by
  have hq' : σ q ∈ I := by rw [hσq]; exact hq
  rw [tateXK_map σ hσI σK hcompat a w q hq hq' hw]
  congr 1

/-- ★★**`σ q = q` の場合の `Y` 座標の自然性**。 -/
theorem tateYK_map_fixed {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [Algebra R K]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σ r))
    (a w q : R) (hq : q ∈ I) (hσq : σ q = q) (hw : w ∈ I) :
    σK (tateYK a w q hq : K) = tateYK (σ a) (σ w) q hq := by
  have hq' : σ q ∈ I := by rw [hσq]; exact hq
  rw [tateYK_map σ hσI σK hcompat a w q hq hq' hw]
  congr 1

/-! ## ★★★非特異性も移る -/

/-- ★★★**`σ` を施した座標も非特異**。

★元の非特異性（`nonsingular_tateK`）を座標の自然性で移すだけ。 -/
theorem nonsingular_tateK_map {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [Algebra R K]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σ r))
    (a w q : R) (hq : q ∈ I) (hσq : σ q = q) (hw : w ∈ I)
    (haw : σ a * σ w = q) (hwu : IsUnit (1 - σ w))
    (hne : algebraMap R K (1 - σ a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (σK (tateXK a w q hq)) (σK (tateYK a w q hq)) := by
  rw [tateXK_map_fixed σ hσI σK hcompat a w q hq hσq hw,
    tateYK_map_fixed σ hσI σK hcompat a w q hq hσq hw]
  exact nonsingular_tateK (σ a) (σ w) q hq haw hwu hne hΔ

/-! ## ★★★★★★★★★★★段 2c の到達点 -/

/-- ★★★★★★★★★★★**Tate 一意化の点は `σ` で自然**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> `tatePtPair (σa) (σw) q … = Point.some (σ_K X) (σ_K Y) …`

★★すなわち **`σ` を施したデータの点は、元の点の座標に `σ` を施したもの**である。

★★★これが `galois-equivariant-tate-uniformization` の段 2c
（同変性の**本体**）の到達点である。
★機構は座標の自然性（`TateCoordNatural.lean`）＋ `Point.some` の合同だけ。 -/
theorem tatePtPair_map {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σ r))
    (a w q : R) (hq : q ∈ I) (hσq : σ q = q) (hw : w ∈ I)
    (haw : σ a * σ w = q) (hwu : IsUnit (1 - σ w))
    (hne : algebraMap R K (1 - σ a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (hns : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (σK (tateXK a w q hq)) (σK (tateYK a w q hq))) :
    tatePtPair (σ a) (σ w) q hq haw hwu hne hΔ
      = Point.some (σK (tateXK a w q hq)) (σK (tateYK a w q hq)) hns := by
  unfold tatePtPair
  congr 1
  · exact (tateXK_map_fixed σ hσI σK hcompat a w q hq hσq hw).symm
  · exact (tateYK_map_fixed σ hσI σK hcompat a w q hq hσq hw).symm

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**段 2c の到達点**である。
段 3（有限拡大）・段 4（帰納極限）は含まない。 -/

def tateXK_map_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——母数を固定した場合の X 座標の自然性)",
    sectionId := "genell-def-3-3" }

def nonsingular_tateK_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——σ を施した座標も非特異)",
    sectionId := "genell-def-3-3" }

def tatePtPair_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——点の自然性。有限拡大と帰納極限は含まない)",
    sectionId := "genell-def-3-3" }

def tatePtPair_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tateXK_map / tateYK_map(座標の自然性——段 2c-2)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateXK_map") 15,
    .citation "[ABC3]" "nonsingular_tateK(座標が曲線の上にあること)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.nonsingular_tateK") 15,
    .implicitStep
      ("★★★★**同変性の本体(段 2c)がこれで閉じた。**" ++
       "★積み上げは evalAdic → adicSum → tateXtail → tateXpairE → tateXK → tatePtPair " ++
       "の 5 段で、どれも「一意 ⟹ 自然」か環演算だけであった") 15,
    .implicitStep
      ("★母数を固定する場合(σ q = q)に絞った——Point.some の合同を取るには" ++
       "同じ曲線の上でなければならないからである。" ++
       "★★G_K は q ∈ K を固定するので、これが実際の場合である") 15,
    .implicitStep
      ("★★残る段 3(有限拡大 K ⊆ L)・段 4(K̄ への帰納極限)は" ++
       "**新しい数学ではなく組み立て**である") 15 ]

end ABC3.Found.GaloisRep
