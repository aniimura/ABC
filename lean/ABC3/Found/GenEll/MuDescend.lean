/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma32MuLine
import ABC3.Found.GaloisRep.TateGaloisStab
import Mathlib.FieldTheory.Galois.Infinite

/-!
# 第 945 ブロック —— **★★★★★★★★★★★★★★★★有理点なら `ζ` は下の体にある**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★★★★★★★★★これは何か——**拡大の層を消す**

第 906 の `exists_mu_point_of_stable` は、位数 `l` の `G_K`-安定な点 `P` に対し
`ζ : Lˣ`（`ζ^l = 1`、`P = tatePhi([ζ])`）を与える。
☆しかし `ζ` は**大きい体 `L` の元**であり、Vélu の計算は `K` の上でやりたい。

★★**だが `P` が `K`-有理なら、`ζ` は自動的に `K` に入る**。

☆理由は短い:

1. `P` が `σ` で不変 ⇒ `tatePhi([σζ]) = tatePhi([ζ])`（`tatePhi_pointMap`、第 907）
2. `tatePhi` は単射（第 863）なので `[σζ] = [ζ]`、すなわち `σζ = ζ·Q^k`
3. 付値を取ると `k·v(Q) = 0`、`v(Q) > 0` なので `k = 0`、つまり `σζ = ζ`
4. `IsGalois K L` なので `ζ ∈ K`

★★★★**これで第 943 の「残る層 (A)（`μ_l` を含む有限 Galois 拡大を立てる）」が
不要になる**——`Lemma 3.5` の `Q` ははじめから `L`-有理な点だからである。

| 定理 | 内容 |
|---|---|
| `mu_fixed_of_pointMap_eq` | ★★★★★★★★★★★★**有理点の `ζ` は `σ` で動かない** |
| `exists_mu_base_of_rational` | ★★★★★★★★**`ζ` は `K` から来る** |
| `exists_mu_point_rational` | ★★★★★★★★★★★★★★★★**第 906 と繋いだ形** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K L : Type} [Field K] [Field L] [DecidableEq L]
  [Algebra R K] [Algebra K L] [Algebra R L] [IsScalarTower R K L]

/-! ## ★★★★★★★★★★★★有理点の `ζ` は `σ` で動かない -/

/-- ★★★★★★★★★★★★**`tatePhi([ζ])` が `σ` で不変なら `σζ = ζ`**。

☆`[σζ] = [ζ]` から `σζ = ζ·Q^k` を取り、付値で `k = 0` を出す。
`v(σζ) = v(ζ)`（`hσv`）と `v(Q) > 0`（`TateSetup`）だけである——
★**`ζ` が 1 の冪根であることさえ使わない**。 -/
theorem mu_fixed_of_pointMap_eq (S : TateSetup R I L) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Δ ≠ 0)
    (hσv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ), vAdd S.v (Units.map (σ : L →* L) u) = vAdd S.v u)
    (ζ : Lˣ)
    (hrat : ∀ σ : L ≃ₐ[K] L,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk ζ))
        = tatePhi S hΔ (QuotientGroup.mk ζ)) :
    ∀ σ : L ≃ₐ[K] L, σ (ζ : L) = (ζ : L) := by
  intro σ
  -- ★段 1: `tatePhi` の σ-同変性（第 907）
  have hσU : ∀ u : Lˣ, ((Units.map (σ : L →* L) u : Lˣ) : L)
      = ((σ : L →ₐ[K] L).restrictScalars R) (u : L) := fun _ => rfl
  have hmapeq := tatePhi_pointMap S hΔ ((σ : L →ₐ[K] L).restrictScalars R)
    (Units.map (σ : L →* L)) hσU (hσv σ) (QuotientGroup.mk ζ)
  rw [hrat σ] at hmapeq
  -- ★段 2: `tatePhi` の単射性（第 863）で類の等式にする
  have hcls : (QuotientGroup.mk ζ : Lˣ ⧸ Subgroup.zpowers S.Q)
      = QuotientGroup.mk (Units.map (σ : L →* L) ζ) :=
    tatePhi_injective S hloc hΔ hmapeq
  have hmem := QuotientGroup.eq.1 hcls
  rw [Subgroup.mem_zpowers_iff] at hmem
  obtain ⟨k, hk⟩ := hmem
  -- ★段 3: 付値で `k = 0`
  have hmul : ζ * S.Q ^ k = Units.map (σ : L →* L) ζ := by
    rw [hk, mul_inv_cancel_left]
  have hval := congrArg (vAdd S.v) hmul
  rw [vAdd_mul, vAdd_zpow, hσv σ ζ] at hval
  have hk0 : k * vAdd S.v S.Q = 0 := by linarith
  have hkz : k = 0 := by
    rcases mul_eq_zero.1 hk0 with h | h
    · exact h
    · exact absurd h S.hQ.ne'
  rw [hkz, zpow_zero, mul_one] at hmul
  -- ★段 4: `Units` の等式を体の等式に戻す
  simpa using congrArg (Units.val) hmul.symm

/-! ## ★★★★★★★★`ζ` は `K` から来る -/

/-- ★★★★★★★★**`ζ` は `K` から来る**——Galois 降下。

☆mathlib の `InfiniteGalois.mem_range_algebraMap_iff_fixed` を当てるだけである。 -/
theorem exists_mu_base_of_rational [IsGalois K L]
    (S : TateSetup R I L) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Δ ≠ 0)
    (hσv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ), vAdd S.v (Units.map (σ : L →* L) u) = vAdd S.v u)
    (ζ : Lˣ)
    (hrat : ∀ σ : L ≃ₐ[K] L,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk ζ))
        = tatePhi S hΔ (QuotientGroup.mk ζ)) :
    ∃ z : K, algebraMap K L z = (ζ : L) :=
  (InfiniteGalois.mem_range_algebraMap_iff_fixed (ζ : L)).2
    (mu_fixed_of_pointMap_eq S hloc hΔ hσv ζ hrat)

/-! ## ★★★★★★★★★★★★★★★★第 906 と繋いだ形 -/

/-- ★★★★★★★★★★★★★★★★**`K`-有理な位数 `l` の点は `K` の中の `μ_l` の点**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★★★★**2026-09-01（第 945）**——これが第 943 の「残る層 (A)」を消す。
`Lemma 3.5` の `Q` ははじめから `L`-有理な点だから、
`μ_l` を含む有限 Galois 拡大を別に立てる必要はなく、
☆**`ζ` はもとの体に自動的に入る**。 -/
theorem exists_mu_point_rational [IsGalois K L] {l : ℕ} (hl : l.Prime)
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
    (hrat : ∀ σ : L ≃ₐ[K] L,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R) P = P) :
    ∃ ζ : Lˣ, ζ ^ l = 1 ∧ (∃ z : K, algebraMap K L z = (ζ : L))
      ∧ P = tatePhi S hΔ (QuotientGroup.mk ζ) := by
  -- ★段 1: 有理性は安定性（`c = 1`）を含む
  have hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R) P = c • P :=
    fun σ => ⟨1, by rw [hrat σ, one_zsmul]⟩
  -- ★段 2: 第 906
  obtain ⟨ζ, hζl, hζP⟩ := exists_mu_point_of_stable hl S hloc hI hquad hvalring hvR hvI
    hq0 hΔ hσv q₀ hQ v hqinf hcop Φ hΦ P hP hstab
  -- ★段 3: 第 945
  refine ⟨ζ, hζl, ?_, hζP⟩
  refine exists_mu_base_of_rational S hloc hΔ hσv ζ ?_
  intro σ
  rw [← hζP]
  exact hrat σ

/-! ## ★出典の紐付け(`.src`) -/

def mu_fixed_of_pointMap_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(有理点の ζ は σ で動かない)",
    sectionId := "genell-lemma-3-2" }

def exists_mu_base_of_rational.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(ζ は K から来る——Galois 降下)",
    sectionId := "genell-lemma-3-2" }

def exists_mu_point_rational.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(K-有理な位数 l の点は K の中の μ_l の点)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GenEll
