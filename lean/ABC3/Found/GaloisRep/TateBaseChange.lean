/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TatePointNatural

/-!
# Tate 一意化は**体の拡大と両立する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★段 3b —— 点の水準での両立

`galois-equivariant-tate-uniformization` の段 3 は
「有限拡大 `K ⊆ L` に対する自然性」であり、
★座標側（段 3a）は `TateCoordNatural.lean` で閉じた
（補題が `σ : R →+* R′`・`σ_K : K →+* L` の形で一般化されている）。

★★本ファイルはそれを**点の水準**へ上げる:

> `K`-代数準同型 `φ : K →ₐ[R] L` に対し、
> **`L` の上の点は `K` の上の点の座標に `φ` を施したもの**である。

## ★★★★★★`Point.map` を経由しない理由（逸脱の記録）

mathlib には `WeierstrassCurve.Affine.Point.map`（`F →ₐ[S] K` に沿った点の写像）が
**ある**。★しかし `tatePtPair` は `(tateCurveAt q).map (algebraMap R K)` の上にあり、
`Point.map` は `W.baseChange F` の上にある。
★★両者は**定義上等しい**が `instances` 透明度では合わず、
`rw [Point.map_some]` も `simp only` も当たらない
（2026-08-27 実測。`Application type mismatch: … (tateCurveAt q hq)⁄K`）。

★★★そこで**`Point.map` を経由せず、結論を直接述べる**:

    `tatePtPair a w q … (L の上) = Point.some (φ X_K) (φ Y_K) …`

★これは「`L`-点は `K`-点の像である」ことそのものであり、数学的内容は同じである。
★★`Point.map` での言い換えは mathlib 側の配管であって、消費側は要らない。

## ★残っている段（明示）

★段 4（`K̄ = ∪_L L` への帰納極限）が残る。
★★段 3b（本ファイル）があれば形式的な組み立てである。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★★拡大した体でも非特異 -/

/-- ★★★**`φ` を施した座標は `L` の上でも非特異**。

★座標の自然性（段 3a）で `L` の座標に直してから
`nonsingular_tateK` を当てるだけ。 -/
theorem nonsingular_tateK_baseChange {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K L : Type} [Field K] [Field L] [Algebra R K] [Algebra R L] (φ : K →ₐ[R] L)
    (a w q : R) (hq : q ∈ I) (hw : w ∈ I) (haw : a * w = q) (hwu : IsUnit (1 - w))
    (hne' : algebraMap R L (1 - a) ≠ 0)
    (hΔ' : ((tateCurveAt q hq).map (algebraMap R L)).toAffine.Δ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R L)).toAffine.Nonsingular
      (φ (tateXK a w q hq : K)) (φ (tateYK a w q hq : K)) := by
  have hx : φ (tateXK a w q hq : K) = (tateXK a w q hq : L) :=
    tateXK_map (RingHom.id R) (fun x hx => hx) (φ : K →+* L)
      (fun r => φ.commutes r) a w q hq hq hw
  have hy : φ (tateYK a w q hq : K) = (tateYK a w q hq : L) :=
    tateYK_map (RingHom.id R) (fun x hx => hx) (φ : K →+* L)
      (fun r => φ.commutes r) a w q hq hq hw
  rw [hx, hy]
  exact nonsingular_tateK a w q hq haw hwu hne' hΔ'

/-! ## ★★★★★★★★★★段 3b の到達点 -/

/-- ★★★★★★★★★★**Tate 一意化は体の拡大と両立する**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> `K`-代数準同型 `φ : K →ₐ[R] L` に対し
> **`L` の上の点 = `Point.some (φ X_K) (φ Y_K)`**

★★すなわち `L` の上の Tate 点は、`K` の上の点の座標を `φ` で送ったものである。
★★★これが `G_K` の作用を `K̄ = ∪_L L` へ持ち上げるための両立条件である。

★機構は座標の自然性（段 3a、`TateCoordNatural.lean`）＋ `Point.some` の合同だけ。 -/
theorem tatePtPair_baseChange {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K L : Type} [Field K] [Field L] [DecidableEq L] [Algebra R K] [Algebra R L]
    (φ : K →ₐ[R] L)
    (a w q : R) (hq : q ∈ I) (hw : w ∈ I) (haw : a * w = q) (hwu : IsUnit (1 - w))
    (hne' : algebraMap R L (1 - a) ≠ 0)
    (hΔ' : ((tateCurveAt q hq).map (algebraMap R L)).toAffine.Δ ≠ 0)
    (hns : ((tateCurveAt q hq).map (algebraMap R L)).toAffine.Nonsingular
      (φ (tateXK a w q hq : K)) (φ (tateYK a w q hq : K))) :
    tatePtPair a w q hq haw hwu hne' hΔ'
      = Point.some (φ (tateXK a w q hq : K)) (φ (tateYK a w q hq : K)) hns := by
  have hx : φ (tateXK a w q hq : K) = (tateXK a w q hq : L) :=
    tateXK_map (RingHom.id R) (fun x hx => hx) (φ : K →+* L)
      (fun r => φ.commutes r) a w q hq hq hw
  have hy : φ (tateYK a w q hq : K) = (tateYK a w q hq : L) :=
    tateYK_map (RingHom.id R) (fun x hx => hx) (φ : K →+* L)
      (fun r => φ.commutes r) a w q hq hq hw
  unfold tatePtPair
  congr 1
  · exact hx.symm
  · exact hy.symm

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**段 3b** である。
段 4（`K̄` への帰納極限）は含まない。 -/

def nonsingular_tateK_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——拡大した体でも非特異)",
    sectionId := "genell-def-3-3" }

def tatePtPair_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——体の拡大との両立。帰納極限は含まない)",
    sectionId := "genell-def-3-3" }

def tatePtPair_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tateXK_map / tateYK_map(座標の自然性——段 3a)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateXK_map") 15,
    .implicitStep
      ("★★★逸脱: mathlib の Point.map(F →ₐ[S] K に沿った点の写像)は**あるが**、" ++
       "tatePtPair は (tateCurveAt q).map (algebraMap R K) の上にあり " ++
       "Point.map は W.baseChange F の上にある。" ++
       "★両者は定義上等しいが instances 透明度では合わず rw も simp も当たらない" ++
       "(2026-08-27 実測)。" ++
       "★★そこで Point.map を経由せず結論を直接述べた——数学的内容は同じである") 15,
    .implicitStep
      ("★★残る段 4(K̄ = ∪_L L への帰納極限)は、本ファイルがあれば" ++
       "**形式的な組み立て**である") 15 ]

end ABC3.Found.GaloisRep
