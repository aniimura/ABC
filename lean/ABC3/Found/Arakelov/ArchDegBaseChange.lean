/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.EmbeddingCount
import ABC3.Found.Arakelov.PullbackNorm
import ABC3.Found.Arakelov.SpanPullSec
import ABC3.Found.Arakelov.ArchDeg
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★底変換の**アルキメデス側が閉じた** `archDeg_K(f^*L̄) = archDeg_F(L̄)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> Note that if we set degF = (1/[F:Q])·degF, then for any finite extension K of F, it follows that degK(L|Spec(OK)) = degF (L)

## ★★★★★★★★★★これは何か

原文が「`deg_K(L|_{Spec 𝓞_K}) = deg_F(L)`」と書く**正規化次数の底変換不変性**の、
アルキメデス部分である。

## ★★★★★機構（3 段、すべて在庫）

    embSpecPoint_comp   : `embSpecPoint K τ ≫ f = embSpecPoint F (τ|_F)`（★本ファイル、`congr 1` で出る）
    norm_pullSec        : `|f^*s|(q) = |s|(q ≫ f)`（`§9-784`）
    sum_over_extensions : `Σ_τ h(τ|_F) = [K:F]·Σ_σ h(σ)`（`§9-795`）

★あとは `[K:ℚ] = [F:ℚ]·[K:F]`（`Module.finrank_mul_finrank`）で
`[K:F]` 倍と正規化が打ち消し合う。

## ★測定の記録

`embSpecPoint K τ ≫ Spec.map (algebraMap 𝓞_F 𝓞_K) = embSpecPoint F (τ|_F)` は
`Spec.map_comp` のあと **`congr 1` だけで閉じる**——
`𝓞_F → 𝓞_K → K` と `𝓞_F → F → K` の合成が**定義的に等しい**からである（2026-08-28 実測）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

set_option maxHeartbeats 1000000 in
/-- ★★★★**複素点は底変換と両立する**。

    `embSpecPoint K τ ≫ Spec.map (algebraMap 𝓞_F 𝓞_K) = embSpecPoint F (τ|_F)`

★`Spec.map_comp` のあと `congr 1` だけで閉じる
——`𝓞_F → 𝓞_K → K` と `𝓞_F → F → K` の合成が定義的に等しいからである。 -/
theorem embSpecPoint_comp (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (τ : K →+* ℂ) :
    embSpecPoint K τ ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))
      = embSpecPoint F (τ.comp (algebraMap F K)) := by
  show Spec.map (CommRingCat.ofHom (τ.comp (algebraMap (𝓞 K) K)))
      ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) = _
  rw [← Spec.map_comp]
  congr 1

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**底変換のアルキメデス側**——`archDeg` は底変換で不変である。

原文 (GenEll p.4):
> Note that if we set degF = (1/[F:Q])·degF, then for any finite extension K of F, it follows that degK(L|Spec(OK)) = degF (L)

★機構: `norm_pullSec`（`§9-784`）で各点のノルムを元のノルムに直し、
`sum_over_extensions`（`§9-795`）で `Σ_τ = [K:F]·Σ_σ` にまとめ、
`[K:ℚ] = [F:ℚ]·[K:F]` で正規化と打ち消す。 -/
theorem archDeg_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (L : AMetric (Spec (CommRingCat.of (𝓞 F)))) (s : (L.sheaf.obj (op ⊤) : Type)) :
    archDeg K (AMetricPullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L)
        (pullSecT (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L.sheaf s)
      = archDeg F L s := by
  have hnorm : ∀ τ : K →+* ℂ,
      (AMetricPullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L).norm
          (pullSecT (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L.sheaf s)
          (embSpecPoint K τ)
        = L.norm s (embSpecPoint F (τ.comp (algebraMap F K))) := by
    intro τ
    have h := norm_pullSec (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L s
      (embSpecPoint K τ)
    rw [embSpecPoint_comp F K τ] at h
    exact h
  have hfin : (Module.finrank ℚ F) * (Module.finrank F K) = Module.finrank ℚ K :=
    Module.finrank_mul_finrank ℚ F K
  have hFpos : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have hKFpos : (0 : ℝ) < (Module.finrank F K : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := F) (M := K)
  show -(∑ τ : K →+* ℂ, Real.log
      ((AMetricPullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L).norm
        (pullSecT (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L.sheaf s)
        (embSpecPoint K τ))) / (Module.finrank ℚ K : ℝ) = _
  simp only [hnorm]
  rw [sum_over_extensions F K (fun σ => Real.log (L.norm s (embSpecPoint F σ)))]
  show _ = -(∑ σ : F →+* ℂ, Real.log (L.norm s (embSpecPoint F σ))) / (Module.finrank ℚ F : ℝ)
  rw [← hfin]
  push_cast
  field_simp

/-! ### ★出典の紐付け(`.src`) -/

def embSpecPoint_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(複素点は底変換と両立する)",
    sectionId := "genell-def-1-1-ii" }

def archDeg_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(底変換のアルキメデス側——archDeg は底変換で不変)",
    sectionId := "genell-def-1-1-ii" }

def archDeg_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "norm_pullSec(|f^*s|(q) = |s|(q ≫ f)、§9-784)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.norm_pullSec") 3,
    .citation "[ABC3]" "sum_over_extensions(Σ_τ h(τ|_F) = [K:F]·Σ_σ h(σ)、§9-795)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.sum_over_extensions") 4,
    .citation "[mathlib]" "Module.finrank_mul_finrank([K:ℚ] = [F:ℚ]·[K:F])"
      (.inMathlib "Module.finrank_mul_finrank") 4,
    .implicitStep
      ("★残るのは有限側の組み立てである: §9-791 の同型 Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L) と " ++
       "§9-794 の数え上げ #(𝓞_K ⊗ Q) = #(Q)^{[K:F]} を組んで " ++
       "degFin_K = [K:F]·degFin_F を出す") 4 ]

end ABC3.Found.Arakelov
