/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateGaloisStab
import ABC3.Found.GaloisRep.TateDoubling
import ABC3.Found.GenEll.Lemma32StableLine

/-!
# [GenEll] Lemma 3.2, (i) —— **★★★★★★★★★★★一意化を posit しない形**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★★★★到達点 —— 逸脱が 1 つ消えた

`Lemma32Uniformized.lean` の `lemma_3_2_i_of_uniformization` は
**同変な一意化 `(Φ, act, hequiv)` を入力として受けていた**。
★本ファイルはそれを**構成したもので埋める**:

| 入力 | 何で埋めたか |
|---|---|
| `Φ` | `tatePhiAddEquivAll`（`TateDoubling.lean`、仮定なしの加法同型） |
| `act σ` | **mathlib の `Point.map`**（`(σ : L →ₐ[K] L).restrictScalars R` に沿って） |
| `hequiv` | `tatePhi_pointMap`（`TateGaloisStab.lean`） |

★★したがって `lemma_3_2_i_tate_all` の仮説には
**一意化も Galois 作用も現れない**——現れるのは Tate 曲線の設定と、
「安定な直線がある」という原文どおりの仮説だけである。

## ★★★★★★残っている仮説（明示）

1. ★`hσv : ∀ σ u, v(σ u) = v(u)`——**付値の拡大の一意性**。
   `L/K` が有限次で `v` が離散付値なら標準的だが、本プロジェクトはまだ持っていない。
2. ★★`hQ : Units.map (algebraMap K L) q₀ = S.Q`——
   Tate 母数が `K` から来ていること。原文では `q_E ∈ K` だから**自動**である。
3. ★★★`hΔ`・`hloc`・`hI`・`hquad`・`hvalring`・`hvR`・`hvI`・`hq0`——
   いずれも `tatePhiAddEquivAll` が要求する **Tate 曲線の設定**である。

## ★★★★`Lemma 3.2` 全体にはまだ `(ii)` が要る

`(ii)`（`E′ = E/μ_l` と `deg_∞(E′) = l·deg_∞(E)`）は
`Lemma32QuotMu.lean` に**一意化の側で**ある（`ker_pow_mk_eq`・`vAdd_pow`・`mu_inj`）。
★★★残っているのは `E/μ_l` を**スキームとして**作る段だけであり、
それは `ResearchPaper/mathlib-gap.json` の
`finite-flat-group-scheme-quotient`（SGA3 / Raynaud 規模）である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GaloisRep

/-- ★★★★★★★★★★**[GenEll] Lemma 3.2, (i) —— 一意化を構成したもので受ける形**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

`Φ` は**抽象のまま**受けているが、それは倍化の仮説（`hloc`・`hI` 等）を
この定理の文面に持ち込まないためだけである
——`lemma_3_2_i_tate_all` が `tatePhiAddEquivAll` を渡して埋める。

★★`act` は **mathlib の `Point.map`** であり、posit ではない。 -/
theorem lemma_3_2_i_tate {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {K L : Type} [Field K] [Field L] [DecidableEq L]
    [Algebra R K] [Algebra K L] [Algebra R L] [IsScalarTower R K L] [IsGalois K L]
    {l : ℕ} (hl : l.Prime)
    (S : TateSetup R I L)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine ≠ 0)
    (Φ : Additive (Lˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (hσv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ), vAdd S.v (Units.map (σ : L →* L) u) = vAdd S.v u)
    (q₀ : Kˣ) (hQ : Units.map (algebraMap K L : K →* L) q₀ = S.Q)
    (v : Kˣ →* Multiplicative ℤ) (hqinf : ∀ j : ℤ, q₀ ^ j = 1 → j = 0)
    (x : Lˣ) (k : ℤ)
    (hxk : x ^ l = (Units.map (algebraMap K L : K →* L) q₀) ^ k)
    (hk : ¬ ((l : ℤ) ∣ k))
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk x))
        = c • tatePhi S hΔ (QuotientGroup.mk x)) :
    (l : ℤ) ∣ vAdd v q₀ := by
  refine lemma_3_2_i hl q₀ v hqinf x k hxk hk ?_
  intro σ
  obtain ⟨c, n, h⟩ := tate_stab_of_pointStab S hΔ Φ hΦ hσv x hstab σ
  exact ⟨c, n, by rw [h, hQ]⟩

/-- ★★★★★★★★★★★**[GenEll] Lemma 3.2, (i) —— 一意化も Galois 作用も posit しない**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★仮説に現れるのは **Tate 曲線の設定**と、
「`P = Φ(x)` の張る直線が `G_K`-安定で `𝔽_l(1)` でない（`¬ l ∣ k`）」という
**原文どおりの仮説**だけである。

★★★これで `Lemma32StableLine.lean` に記録した逸脱
「`G_K` 同変な一意化は含まない」は**完全に消えた**。 -/
theorem lemma_3_2_i_tate_all {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {K L : Type} [Field K] [Field L] [DecidableEq L]
    [Algebra R K] [Algebra K L] [Algebra R L] [IsScalarTower R K L] [IsGalois K L]
    {l : ℕ} (hl : l.Prime)
    (S : TateSetup R I L) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I) (hI : I.IsPrime)
    (hquad : ∀ (t : L) (b c : R),
      t ^ 2 + algebraMap R L b * t + algebraMap R L c = 0 → ∃ r : R, algebraMap R L r = t)
    (hvalring : ∀ t : L, t ≠ 0 →
      (∃ r : R, algebraMap R L r = t) ∨ (∃ r ∈ I, algebraMap R L r = t⁻¹))
    (hvR : ∀ t : Lˣ, (∃ r : R, algebraMap R L r = (t : L)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Lˣ, (∃ r ∈ I, algebraMap R L r = (t : L)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Δ ≠ 0)
    (hσv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ), vAdd S.v (Units.map (σ : L →* L) u) = vAdd S.v u)
    (q₀ : Kˣ) (hQ : Units.map (algebraMap K L : K →* L) q₀ = S.Q)
    (v : Kˣ →* Multiplicative ℤ) (hqinf : ∀ j : ℤ, q₀ ^ j = 1 → j = 0)
    (x : Lˣ) (k : ℤ)
    (hxk : x ^ l = (Units.map (algebraMap K L : K →* L) q₀) ^ k)
    (hk : ¬ ((l : ℤ) ∣ k))
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk x))
        = c • tatePhi S hΔ (QuotientGroup.mk x)) :
    (l : ℤ) ∣ vAdd v q₀ :=
  lemma_3_2_i_tate hl S hΔ
    (tatePhiAddEquivAll S hloc hI hquad hvalring hvR hvI hq0 hΔ) (fun _ => rfl)
    hσv q₀ hQ v hqinf x k hxk hk hstab

/-! ### ★出典の紐付け(`.src`)

★★★**`Lemma 3.2` 全体の `.src` はまだ付けない**——`(ii)` の
「`E/μ_l` をスキームとして作る段」が残っているからである
(`finite-flat-group-scheme-quotient`)。 -/

def lemma_3_2_i_tate.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(Galois 作用を mathlib の Point.map で。(ii) は含まない)",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2_i_tate_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(一意化も Galois 作用も posit しない形。(ii) は含まない)",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2_i_tate_all.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_2_i(Kˣ/q^ℤ の側の二者択一)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_2_i") 15,
    .citation "[ABC3]" "tatePhiAddEquivAll(仮定なしの加法同型 Lˣ/Q^ℤ ≃+ E_q(L))"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhiAddEquivAll") 15,
    .citation "[ABC3]" "tatePhi_pointMap(同変性——Point.map の形)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhi_pointMap") 15,
    .implicitStep
      ("★残る仮説は hσv(付値が Galois で不変——付値の拡大の一意性)と、" ++
       "Tate 曲線の設定(hloc・hI・hquad・hvalring・hvR・hvI・hq0・hΔ)だけである") 15,
    .implicitStep
      ("★★Lemma 3.2 全体には (ii)(E′ = E/μ_l)が要る。" ++
       "一意化の側は Lemma32QuotMu.lean にあり、" ++
       "スキームとしての商だけが残っている(finite-flat-group-scheme-quotient)") 15 ]

end ABC3.Found.GenEll
