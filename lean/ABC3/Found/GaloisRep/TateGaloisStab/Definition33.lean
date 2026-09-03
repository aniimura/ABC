/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateUniformizationEquivariant

/-!
# TateGaloisStab —— `[GenEll] Definition 3.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★★★母数を固定することは自動 -/

/-- ★★★★**`R`-代数射は Tate 母数 `Q` を固定する**。

★`S.Q` の台は `algebraMap R L S.q` であり、`R`-代数射はそれを動かさない。
★★したがって `tatePhi_map` の仮説 `hσq` は**仮説ではない**。 -/
theorem tateSetup_Q_fixed {R : Type} [CommRing R] {I : Ideal R}
    {K : Type} [Field K] [Algebra R K] (S : TateSetup R I K)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K)) :
    σU S.Q = S.Q := by
  apply Units.ext
  rw [hσU, ← S.hQq, σA.commutes]

/-- ★★★**単数群への持ち上げは単射**——体の準同型が単射だからである。 -/
theorem tateSetup_σU_injective {R : Type} [CommRing R] {I : Ideal R}
    {K : Type} [Field K] [Algebra R K] (_S : TateSetup R I K)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K)) :
    Function.Injective σU := by
  intro u v h
  apply Units.ext
  have : σA (u : K) = σA (v : K) := by rw [← hσU, ← hσU, h]
  exact (σA : K →+* K).injective this

/-- ★★★**単位類でない類は `σ` で単位類に落ちない**。

★`tatePhi_map` は `σ_* c ≠ 1` を要求するので、それを `c ≠ 1` から出す段が要る。 -/
theorem tateSetup_quotMap_ne_one {R : Type} [CommRing R] {I : Ideal R}
    {K : Type} [Field K] [Algebra R K] (S : TateSetup R I K)
    (σU : Kˣ →* Kˣ) (hσq : σU S.Q = S.Q) (hinj : Function.Injective σU)
    {c : Kˣ ⧸ Subgroup.zpowers S.Q} (hc : c ≠ 1) :
    QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) c ≠ 1 := by
  intro h
  refine hc ?_
  obtain ⟨x, rfl⟩ := QuotientGroup.mk_surjective c
  have hx : (QuotientGroup.mk (σU x) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := h
  obtain ⟨n, hn⟩ := (QuotientGroup.eq_one_iff _).1 hx
  have hn' : σU (S.Q ^ n) = σU x := by
    rw [map_zpow, hσq]; exact hn
  have hxn : x = S.Q ^ n := (hinj hn').symm
  rw [hxn]
  exact (QuotientGroup.eq_one_iff _).2 ⟨n, rfl⟩

/-! ## ★★★★★★★★★★★同変性を `Point.map` の形で -/

/-- ★★★★★★★★★★★**Tate 一意化 `Φ` は `Point.map` と可換**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

    `Point.map σ_A (Φ c) = Φ (σ_* c)`

★★左辺の `Point.map` は **mathlib の加法群準同型**である
——すなわち本定理は「`Φ` は `G_K`-加群の同型」を意味する。

★★★機構は 2 段だけ:
* 単位類は原点に行き、`Point.map` は原点を原点に送る（`map_zero`）
* それ以外では `tatePhi_map`（座標の式）に `Point.map_some` を当てる

★`rw [Point.map_some]` は `instances` 透明度で落ちる
（`(tateCurveAt q hq).map (algebraMap R K)).toAffine` と `(tateCurveAt q hq)⁄K` の食い違い）。
★★`refine (Point.map_some (S := R) σA _).trans ?_` なら通る
（`tools\lean-idioms.md`)。 -/
theorem tatePhi_pointMap {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K] (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K))
    (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    Point.map (S := R) σA (tatePhi S hΔ c)
      = tatePhi S hΔ (QuotientGroup.map _ _ σU
          (zpowers_le_comap_self S.Q σU (tateSetup_Q_fixed S σA σU hσU)) c) := by
  by_cases hc : c = 1
  · subst hc
    rw [map_one, tatePhi_one]
    exact map_zero _
  · have hc' := tateSetup_quotMap_ne_one S σU (tateSetup_Q_fixed S σA σU hσU)
      (tateSetup_σU_injective S σA σU hσU) hc
    rw [tatePhi_eq S hΔ hc, tatePtPair]
    refine (Point.map_some (S := R) σA _).trans ?_
    exact (tatePhi_map S hΔ (RingHom.id R) (fun x hx => hx) (σA : K →+* K)
      (fun r => σA.commutes r) σU hσU (tateSetup_Q_fixed S σA σU hσU) hσv hc' _).symm

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`(Tate 一意化)の同変性を、消費側(`Lemma 3.2, (i)`)が使う形へ
翻訳する段である。 -/

def tateSetup_Q_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(R-代数射は Tate 母数を固定すること)",
    sectionId := "genell-def-3-3" }

def tatePhi_pointMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化は Point.map と可換——同変性の群の形)",
    sectionId := "genell-def-3-3" }

def tatePhi_pointMap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tatePhi_map(同変性——座標の式、2026-08-27)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhi_map") 15,
    .citation "[mathlib]" "WeierstrassCurve.Affine.Point.map(点の加法群準同型)"
      (.inMathlib "WeierstrassCurve.Affine.Point.map") 15,
    .implicitStep
      ("★仮説として残るのは hσv(付値が Galois で不変)だけである" ++
       "——これは付値の拡大の一意性であり、本ファイルの守備範囲ではない") 15 ]

end ABC3.Found.GaloisRep
