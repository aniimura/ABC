/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma32Tate
import ABC3.Found.GaloisRep.TateMuPoint
import ABC3.Found.GaloisRep.TateVeluPoints

/-!
# 第 906 ブロック —— **★★★★★★★★★★★★★★★★`l` が局所高さと互いに素なら
安定な直線は `μ_l`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★★★これは何か

`lemma_3_2_i_tate_all`（証明済み）は

    安定な直線が `𝔽_l(1)` でない（`¬ l ∣ k`） ⇒ `l ∣ v_K(q_E)`

を与える。★本ブロックは**対偶を取って結論側を出す**:

    `l` が局所高さと互いに素 ⇒ `l ∣ k` ⇒ `[x] = [ζ]`（`ζ^l = 1`）

☆後半は `exists_rootOfUnity_mk_eq`（第 905）である。

★これが原文の「the assumption on l implies, by Lemma 3.2, (i), that at all the
primes of multiplicative reduction, HF corresponds to the subspace Fl(1)」の内容であり、
`Lemma 3.5` に残る 3 つのうち (3) である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GaloisRep

/-- ★★★★★★★★★★★★★★★★**`l` が局所高さと互いに素なら、
`G_K`-安定な直線は `μ_l` の元の類である**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★仮説は `lemma_3_2_i_tate_all` と同じで、`hk : ¬ l ∣ k` の代わりに
**`hcop : ¬ l ∣ v_K(q_E)`**（原文の「`l` is prime to the local heights」）を受ける。 -/
theorem stableLine_is_mu_of_coprime {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
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
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk x))
        = c • tatePhi S hΔ (QuotientGroup.mk x))
    (hcop : ¬ ((l : ℤ) ∣ vAdd v q₀)) :
    ∃ ζ : Lˣ, ζ ^ l = 1 ∧
      (QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers S.Q) = QuotientGroup.mk ζ := by
  -- ★対偶: `¬ l ∣ v_K(q₀)` なら `l ∣ k`
  have hdvd : (l : ℤ) ∣ k := by
    by_contra hk
    exact hcop (lemma_3_2_i_tate_all hl S hloc hI hquad hvalring hvR hvI hq0 hΔ hσv
      q₀ hQ v hqinf x k hxk hk hstab)
  -- ★第 905
  rw [← hQ]
  exact exists_rootOfUnity_mk_eq _ hl.pos x k hxk hdvd

/-! ## ★★★★★★★★★★★★点の側で述べた形 -/

/-- ★★★★★★★★★★★★**位数 `l` の `G_K`-安定な点は `μ_l` の点である**
（`l` が局所高さと互いに素なら）。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★これが `Lemma 3.5` の残り (b)「`Q ⊗ Lv` を `Φ(ζ)` と同定する段」である。

☆道は 3 段:

1. `Φ` は加法同型（全射）なので `P = Φ([x])` と書ける
2. `l • P = 0` から `x^l = q₀^k`（第 916）
3. `stableLine_is_mu_of_coprime`（第 906）で `[x] = [ζ]` -/
theorem exists_mu_point_of_stable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
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
    (hcop : ¬ ((l : ℤ) ∣ vAdd v q₀))
    (Φ : Additive (Lˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (P : ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Point)
    (hP : l • P = 0)
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R) P = c • P) :
    ∃ ζ : Lˣ, ζ ^ l = 1 ∧ P = tatePhi S hΔ (QuotientGroup.mk ζ) := by
  -- ★段 1: `P = Φ([x])`
  obtain ⟨c, hc⟩ := Φ.surjective P
  obtain ⟨x, hx⟩ := QuotientGroup.mk_surjective (s := Subgroup.zpowers S.Q)
    (Additive.toMul c)
  have hPx : tatePhi S hΔ (QuotientGroup.mk x) = P := by
    rw [← hΦ, hx]
    simpa using hc
  -- ★段 2: `x^l = q₀^k`
  obtain ⟨k, hk⟩ := exists_zpow_of_nsmul_tatePhi_eq_zero S hΔ Φ hΦ x l
    (by rw [hPx]; exact hP)
  have hxk : x ^ l = (Units.map (algebraMap K L : K →* L) q₀) ^ k := by
    rw [hQ, hk]
  -- ★段 3: 安定性を `[x]` の側へ
  have hstab' : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk x))
        = c • tatePhi S hΔ (QuotientGroup.mk x) := by
    intro σ
    rw [hPx]
    exact hstab σ
  obtain ⟨ζ, hζl, hζc⟩ := stableLine_is_mu_of_coprime hl S hloc hI hquad hvalring hvR hvI
    hq0 hΔ hσv q₀ hQ v hqinf x k hxk hstab' hcop
  exact ⟨ζ, hζl, by rw [← hPx, hζc]⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_mu_point_of_stable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(位数 l の G_K-安定な点は μ_l の点である。★無条件)",
    sectionId := "genell-lemma-3-2" }

def stableLine_is_mu_of_coprime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(l が局所高さと互いに素なら安定な直線は μ_l。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GenEll
