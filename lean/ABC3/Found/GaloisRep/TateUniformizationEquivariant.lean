/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateBaseChange
import ABC3.Found.GaloisRep.TatePhiNatural
import ABC3.Found.GaloisRep.TatePhi

/-!
# Tate 一意化は **`G_K`-同変**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★到達点 —— `Φ` が同変である

`galois-equivariant-tate-uniformization` の積み上げの**帰結**:

> **`Φ(σ_* c) = σ_K(Φ(c))`**（座標ごと）

すなわち Tate 一意化 `Φ : Kˣ/q^ℤ → E_q(K)` は
`q` を固定し付値を保つ `σ` に対して**同変**である。

## ★★★★★★★★★★帰納極限は**要らない**

台帳は段 4 に「`K̄ = ∪_L L` への帰納極限」を置いていたが、
★**消費側（`Lemma 3.2`）には要らない**:

* `M_l(E) = E(K̄)[l]` の `l`-捩れは**有限拡大 `L = K(E[l])` の中に載る**。
* `L` は完備離散付値体なので段 1（`tatePhiAddEquivAll`）がそのまま使える。
* `G_K` は `E[l]` に **`Gal(L/K)` を通して**作用し、
  その同変性が本ファイルの `tatePhi_map` である（`σ q = q` は `q ∈ K` だから）。

★★したがって `K̄` を直接扱う必要はない。

## ★★★★★★積み上げの全体（本日 6 ブロック）

| 段 | 内容 | ファイル |
|---|---|---|
| 1 | 固定した完備体での**選ばれた**同型 | `TateDoubling.lean`（既存） |
| 2a | `normRep` の自然性 | `NormRepNatural.lean` |
| 2b | `tateAOf` / `tateWOf` の自然性 | `TatePhiNatural.lean` |
| 2c-1 | `evalAdic` の自然性 | `AdicEvalNatural.lean` |
| 2c-2 | `adicSum` 〜 `tateXK` の自然性 | `TateSeriesNatural` / `TateCoordNatural` |
| 2c-3 | `tatePtPair` の自然性 | `TatePointNatural.lean` |
| 3a | 行き先が別の環・別の体でもよい | 上の 3 ファイルを一般化 |
| 3b | 体の拡大との両立 | `TateBaseChange.lean` |
| **帰結** | **`Φ` の同変性** | ★本ファイル |

★★★どの段も「**一意 ⟹ 自然**」か環演算だけであった
（`eq_normRep`・`tateAOf_spec`・`evalAdic_unique`・`adicSum_unique`）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★単位類は原点に行く（同変性の自明な側） -/

/-- ★**単位類は `σ` で動かず、原点に行く**。 -/
theorem tatePhi_map_one {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K] (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (σU : Kˣ →* Kˣ) (hσq : σU S.Q = S.Q) :
    tatePhi S hΔ (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) 1) = 0 := by
  rw [map_one, tatePhi_one]

/-! ## ★★★★★★★★★★★`Φ` の同変性 -/

/-- ★★★★★★★★★★★**Tate 一意化 `Φ` は `σ` で同変**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> **`Φ(σ_* c) = Point.some (σ_K X) (σ_K Y)`**、ただし `Point.some X Y _ = Φ(c)`

すなわち `Φ(σ_* c)` は `Φ(c)` の座標に `σ_K` を施したものである。

★★仮説は 3 つだけ:
* `σ_U S.Q = S.Q`（`q` を固定する——`q ∈ K` なので `G_K` では自動）
* `vAdd (σ_U u) = vAdd u`（付値を保つ——付値の拡大の一意性から）
* `σ I ⊆ I`・両立（`σ_K ∘ algebraMap = algebraMap ∘ σ`）

★★★これが `Lemma 3.2` が要求する **`G_K`-同変な Tate 一意化**である
——`l`-捩れは有限拡大に載るので、`K̄` への帰納極限は要らない。 -/
theorem tatePhi_map {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K] (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σ r))
    (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σK (u : K))
    (hσq : σU S.Q = S.Q) (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    {c : Kˣ ⧸ Subgroup.zpowers S.Q}
    (hc' : QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) c ≠ 1)
    (hns : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Nonsingular
      (σK (tateXK (tateAOf S c) (tateWOf S c) S.q S.hq))
      (σK (tateYK (tateAOf S c) (tateWOf S c) S.q S.hq))) :
    tatePhi S hΔ (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) c)
      = Point.some (σK (tateXK (tateAOf S c) (tateWOf S c) S.q S.hq))
          (σK (tateYK (tateAOf S c) (tateWOf S c) S.q S.hq)) hns := by
  have hqq : σ S.q = S.q := tateSetup_q_map S σ σK hcompat σU hσU hσq
  have hA := tateAOf_map S σ σK hcompat σU hσU hσq hσv c
  have hW := tateWOf_map S σ σK hcompat σU hσU hσq hσv c
  rw [tatePhi_eq S hΔ hc']
  unfold tatePtPair
  congr 1
  · rw [← hA, ← hW]
    exact (tateXK_map_fixed σ hσI σK hcompat (tateAOf S c) (tateWOf S c) S.q S.hq hqq
      (tateWOf_mem S c)).symm
  · rw [← hA, ← hW]
    exact (tateYK_map_fixed σ hσI σK hcompat (tateAOf S c) (tateWOf S c) S.q S.hq hqq
      (tateWOf_mem S c)).symm

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**同変性**である。
★★★`Lemma 3.2` 全体ではない——`M_l(E)` を `(Kˣ/q^ℤ)[l]` と同定する段
（有限拡大 `L = K(E[l])` の選択）は含まない。 -/

def tatePhi_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Φ の G_K-同変性)",
    sectionId := "genell-def-3-3" }

def tatePhi_map_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単位類は σ で動かない)",
    sectionId := "genell-def-3-3" }

def tatePhi_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tateAOf_map / tateWOf_map(係数の自然性——段 2b)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateAOf_map") 15,
    .citation "[ABC3]" "tateXK_map_fixed / tateYK_map_fixed(座標の自然性——段 2c)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateXK_map_fixed") 15,
    .citation "[ABC3]" "tatePhiAddEquivAll(選ばれた同型——段 1)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhiAddEquivAll") 15,
    .otherPaper "[FC]" "Degenerations of Abelian Varieties, Chapter III, Corollary 7.3" 15,
    .implicitStep
      ("★★★★★**帰納極限は要らない**——M_l(E) の l-捩れは" ++
       "有限拡大 L = K(E[l]) の中に載るので、L(完備)で段 1 が使え、" ++
       "G_K は Gal(L/K) を通して作用する。★q ∈ K なので σ q = q は自動である") 15,
    .implicitStep
      ("★★積み上げは 8 段: normRep → tateAOf/tateWOf → evalAdic → adicSum → " ++
       "tateXtail → tateXK → tatePtPair → Φ。" ++
       "★**どの段も「一意 ⟹ 自然」か環演算だけ**であった" ++
       "(eq_normRep・tateAOf_spec・evalAdic_unique・adicSum_unique)") 15,
    .implicitStep
      ("★残る段: M_l(E) を (Kˣ/q^ℤ)[l] と同定すること" ++
       "(有限拡大 L = K(E[l]) の選択)。" ++
       "★★消費側(Lemma32StableLine.lean)は Kˣ/q^ℤ の上で既に閉じている") 15 ]

end ABC3.Found.GaloisRep
