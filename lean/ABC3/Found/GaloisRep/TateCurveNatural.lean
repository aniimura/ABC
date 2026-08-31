/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateUniformizationEquivariant

/-!
# Tate 曲線そのものが**環について自然**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★消費側へ繋ぐための 2 本

`tatePhi_map`（`Φ` の同変性）を `Lemma 3.2` へ繋ぐには、
**「どの環の上の Tate 曲線か」を揃える**必要がある:

* 曲線 `tateCurveAt q` は `q` を含む環の上にある。
* `G_K` は `L` の上で作用するが `R_L` を点ごとには固定しない
  （固定するのは `R_K` である）。

★本ファイルは **`tateCurveAt` が環について自然**であることを取る:

> **`(tateCurveAt q hq).map σ = tateCurveAt (σ q) hq′`**

★★これで `R_K` の上の曲線と `R_L` の上の曲線が `L` の上で一致し、
mathlib の `Point.map`（`F →ₐ[S] K` に沿った点の写像）が `S = R_K` で使える。

## ★★★★★★`Point.map` は `rw` では当たらないが `exact` なら通る

前ブロック（`TateBaseChange.lean`）では
`rw [Point.map_some]` が `instances` 透明度で落ちたため `Point.map` を避けた。
★**しかし `refine (Point.map_some …).trans ?_` なら通る**（2026-08-27 実測）——
`rw` はパターン照合、`exact`/`refine` は既定透明度での `isDefEq` だからである。

★★したがって本ファイルでは `Point.map` を**経由できる**:

> `Point.map φ (tatePtPair a w q …) = tatePtPair a w q …`（`L` の上）

## ★機構

`tateCurveAt q hq = tateCurve.map (evalAdicHom q hq)` なので、
`WeierstrassCurve.map_map` と **`evalAdic_map`（段 2c-1）** だけで出る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★評価準同型の自然性 -/

/-- ★★**`evalAdicHom` は環準同型で自然**（`σ I ⊆ I′` なら）。

★`evalAdic_map`（段 2c-1）を環準同型の水準で述べ直しただけ。 -/
theorem evalAdicHom_map {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {I' : Ideal R'}
    [IsAdicComplete I R] [IsAdicComplete I' R']
    (σ : R →+* R') (hσI : ∀ x ∈ I, σ x ∈ I') (q : R) (hq : q ∈ I) (hq' : σ q ∈ I') :
    σ.comp (evalAdicHom q hq) = evalAdicHom (σ q) hq' := by
  ext f
  show σ (evalAdic f q hq) = evalAdic f (σ q) hq'
  exact (evalAdic_map σ hσI f q hq hq').symm

/-! ## ★★★★★★★★Tate 曲線の自然性 -/

/-- ★★★★★★★★**Tate 曲線は環について自然**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> **`(tateCurveAt q hq).map σ = tateCurveAt (σ q) hq′`**

★★これで「どの環の上の Tate 曲線か」を揃えられる——
`R_K` の上の曲線と `R_L` の上の曲線が `L` の上で一致する。

★機構は `tateCurveAt q hq = tateCurve.map (evalAdicHom q hq)` と
`WeierstrassCurve.map_map`、それに `evalAdic_map`（段 2c-1）だけ。 -/
theorem tateCurveAt_map {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {I' : Ideal R'}
    [IsAdicComplete I R] [IsAdicComplete I' R']
    (σ : R →+* R') (hσI : ∀ x ∈ I, σ x ∈ I') (q : R) (hq : q ∈ I) (hq' : σ q ∈ I') :
    (tateCurveAt q hq).map σ = tateCurveAt (σ q) hq' := by
  unfold tateCurveAt
  rw [WeierstrassCurve.map_map, evalAdicHom_map σ hσI q hq hq']

/-! ## ★★★★★★★★★★`Point.map` を経由した形 -/

/-- ★★★★★★★★★★**mathlib の `Point.map` で述べた体拡大との両立**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> **`Point.map φ (tatePtPair a w q …) = tatePtPair a w q …`**（`L` の上）

★★`TateBaseChange.lean` では `rw [Point.map_some]` が
`instances` 透明度で落ちたため `Point.map` を避けていたが、
**`refine (Point.map_some …).trans ?_` なら通る**（2026-08-27 実測）。

★★★これで `G_K` の作用を mathlib の `Point.map` として扱える
——`Lemma 3.2` が要求する `M_l(E)` の `G_K`-加群構造そのものである。 -/
theorem tatePtPair_pointMap {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K L : Type} [Field K] [Field L] [DecidableEq K] [DecidableEq L]
    [Algebra R K] [Algebra R L] (φ : K →ₐ[R] L)
    (a w q : R) (hq : q ∈ I) (hw : w ∈ I) (haw : a * w = q) (hwu : IsUnit (1 - w))
    (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (hne' : algebraMap R L (1 - a) ≠ 0)
    (hΔ' : ((tateCurveAt q hq).map (algebraMap R L)).toAffine.Δ ≠ 0) :
    (Point.map (S := R) φ) (tatePtPair a w q hq haw hwu hne hΔ)
      = tatePtPair a w q hq haw hwu hne' hΔ' := by
  have hx : φ (tateXK a w q hq : K) = (tateXK a w q hq : L) :=
    tateXK_map (RingHom.id R) (fun x hx => hx) (φ : K →+* L)
      (fun r => φ.commutes r) a w q hq hq hw
  have hy : φ (tateYK a w q hq : K) = (tateYK a w q hq : L) :=
    tateYK_map (RingHom.id R) (fun x hx => hx) (φ : K →+* L)
      (fun r => φ.commutes r) a w q hq hq hw
  refine (Point.map_some (S := R) φ (nonsingular_tateK a w q hq haw hwu hne hΔ)).trans ?_
  unfold tatePtPair
  congr 1

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）を消費側へ繋ぐための配管である。 -/

def tateCurveAt_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線は環について自然)",
    sectionId := "genell-def-3-3" }

def tatePtPair_pointMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Point.map で述べた体拡大との両立)",
    sectionId := "genell-def-3-3" }

def tateCurveAt_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "evalAdic_map(冪級数の adic 値の自然性——段 2c-1)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdic_map") 15,
    .citation "[mathlib]" "WeierstrassCurve.map_map"
      (.inMathlib "WeierstrassCurve.map_map") 15,
    .implicitStep
      ("★★★**`Point.map` は `rw` では当たらないが `exact` なら通る**" ++
       "(2026-08-27 実測)。tatePtPair は (tateCurveAt q).map (algebraMap R K) の上、" ++
       "Point.map は W.baseChange F の上にあり、定義上等しいが " ++
       "instances 透明度では合わない。" ++
       "★rw はパターン照合、refine/exact は既定透明度の isDefEq だからである") 15,
    .implicitStep
      ("★★残る段: M_l(E) を (Kˣ/q^ℤ)[l] と同定すること" ++
       "(有限拡大 L = K(E[l]) の選択と、Φ の全単射性からの捩れの対応)。" ++
       "★消費側(Lemma32StableLine.lean)は Kˣ/q^ℤ の上で既に閉じている") 15 ]

end ABC3.Found.GaloisRep
