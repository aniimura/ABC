/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateUniformizationEquivariant
import ABC3.Found.GaloisRep.TateGaloisStab.Definition33

/-!
# TateGaloisStab —— `[GenEll] Lemma 3.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★★★★★★安定な直線の翻訳 -/

/-- ★★★★★★★**点の側の安定性を単数の側へ移す**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

`P = Φ(x)` の張る直線が `G_K`-安定、すなわち各 `σ` について
`Point.map σ P = c·P` なら、単数の側では

    `σ(x) = x^c · Q^n`

が成り立つ。★これが `lemma_3_2_i` が要求する形そのものである。

★★入力の `Φ` は**構成したもの**でよい——`tatePhiAddEquivAll` と `fun _ => rfl` を渡せる。
★★★`Φ` を抽象のまま受けているのは、倍化の仮説（`hloc`・`hI` 等）を
本ファイルに持ち込まないためだけである。 -/
theorem tate_stab_of_pointStab {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {K L : Type} [Field K] [Field L] [DecidableEq L]
    [Algebra R K] [Algebra K L] [Algebra R L] [IsScalarTower R K L]
    (S : TateSetup R I L)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine ≠ 0)
    (Φ : Additive (Lˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (hσv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ), vAdd S.v (Units.map (σ : L →* L) u) = vAdd S.v u)
    (x : Lˣ)
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk x))
        = c • tatePhi S hΔ (QuotientGroup.mk x)) :
    ∀ σ : L ≃ₐ[K] L, ∃ c n : ℤ, Units.map (σ : L →* L) x = x ^ c * S.Q ^ n := by
  intro σ
  obtain ⟨c, hc⟩ := hstab σ
  have h1 : ∀ u : Lˣ, ((Units.map (σ : L →* L) u : Lˣ) : L)
      = ((σ : L →ₐ[K] L).restrictScalars R) (u : L) := fun u => rfl
  have hmap := tatePhi_pointMap S hΔ ((σ : L →ₐ[K] L).restrictScalars R)
    (Units.map (σ : L →* L)) h1 (hσv σ) (QuotientGroup.mk x)
  have key : tatePhi S hΔ (QuotientGroup.mk (Units.map (σ : L →* L) x))
      = c • tatePhi S hΔ (QuotientGroup.mk x) := by
    rw [← hc, hmap]; rfl
  have e1 : Φ (Additive.ofMul ((QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers S.Q) ^ c))
      = c • Φ (Additive.ofMul (QuotientGroup.mk x)) := by
    rw [← map_zsmul]
    rfl
  have h2 : Φ (Additive.ofMul (QuotientGroup.mk (Units.map (σ : L →* L) x)))
      = Φ (Additive.ofMul ((QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers S.Q) ^ c)) := by
    rw [hΦ, key, e1, hΦ]
  have h3 : (QuotientGroup.mk (Units.map (σ : L →* L) x) : Lˣ ⧸ Subgroup.zpowers S.Q)
      = QuotientGroup.mk (x ^ c) := Additive.ofMul.injective (Φ.injective h2)
  have h4 : (Units.map (σ : L →* L) x) * (x ^ c)⁻¹ ∈ Subgroup.zpowers S.Q := by
    rw [← QuotientGroup.eq_one_iff, QuotientGroup.mk_mul, h3, QuotientGroup.mk_inv,
      mul_inv_cancel]
  obtain ⟨n, hn⟩ := h4
  refine ⟨c, n, ?_⟩
  have hn' : S.Q ^ n = (Units.map (σ : L →* L) x) * (x ^ c)⁻¹ := hn
  rw [hn', ← mul_assoc, mul_comm (x ^ c) _, mul_assoc, mul_inv_cancel, mul_one]

def tate_stab_of_pointStab.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(点の側の安定性を単数の側へ移す段)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
