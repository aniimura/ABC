/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.MuRational
import ABC3.Meta.Claim

/-!
# 第 1410 ブロック —— **★★★★★★★★★★★★★★★★二者択一：`μ_l` 型か、深い代表か**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★★★★★★★★★これは何か——`hcop` を**仮定から外す**ための骨

第 946 の `exists_mu_of_rational` は `l ∤ vAdd(Q)`（＝`hcop`）を仮定して
「位数 `l` の類は `μ_l` の類」を出していた。
★第 1408 で分かったとおり、`hcop` は**底変換で保たれない**
（`jExp P (E⊗M) = e(P∣p)·jExp p E` なので `l ∣ e(P∣p)` で壊れる）。
☆したがって `l ∣ vAdd(Q)` の場合を正面から扱わねばならない。

★★本ブロックは **`hcop` を落として、代わりに二者択一を出す**:

`x^l = Q^k` のとき

| 場合 | 結論 |
|---|---|
| `l ∣ k` | ☆`[x] = [ζ]`（`ζ^l = 1`）——**`μ_l` 型**（第 1404-1406 の道） |
| `l ∤ k` | ★**`0 < vAdd(y) < vAdd(Q)` な代表 `y` が取れる**——**深い代表** |

☆後者の道: `l·vAdd(x) = k·vAdd(Q)` と `l` 素・`l ∤ k` から `l ∣ vAdd(Q)`。
`vAdd(Q) = l·d` と置くと `vAdd(x) = k·d` で、`y := x·Q^{-⌊k/l⌋}` は

    `vAdd(y) = (k - l⌊k/l⌋)·d = (k mod l)·d`

を満たす。`0 < k mod l < l` なので `0 < vAdd(y) < l·d = vAdd(Q)`。

★★★これが効く理由: `vAdd(y) > 0` は **`y ∈ 𝔪`** を意味し、そのとき核の
`x` 座標がすべて `𝔪` に入るので `veluV2 ≡ b₄ ≡ 0 (mod 𝔪)`、
したがって `c₄(veluCurve) = c₄(E_q) + 240v` は単元——**半安定性が出る**。
☆`μ_l` 型のときは逆に `c₄′ ≡ l⁴` が単元（第 1388 の `h4`）。
★どちらの場合も `c₄′` が単元なので、悪い素点は `hcop` なしで閉じる見込みである。

| 定理 | 内容 |
|---|---|
| `mu_or_deep_of_pow_eq_zpow` | ★★★★★★★★類の側の二者択一 |
| `mu_or_deep_point` | ★★★★★★★★★★★★点の側の二者択一 |
| `mu_or_deep_point_dvr` | ★★★★★★★★★★★★★★★★**完備 DVR だけで述べた形** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GaloisRep

/-! ## ★★★★★★★★類の側の二者択一 -/

section Rational

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

omit [IsDomain R] [IsAdicComplete I R] [DecidableEq K] in
/-- ★★★★★★★★**`x^l = Q^k` なら、`[x]` は `μ_l` の類か、深い代表を持つ**（第 1410）。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

☆`l ∤ vAdd(Q)`（第 946 の `hcop`）を**仮定しない**のが要点である。

★`l ∤ k` の側では `l ∣ vAdd(Q)` が自動的に従い（`l·vAdd(x) = k·vAdd(Q)` と `l` の素性）、
`vAdd(Q) = l·d` として `y := x·Q^{-⌊k/l⌋}` の付値が `(k mod l)·d` になる。 -/
theorem mu_or_deep_of_pow_eq_zpow (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    (x : Kˣ) (k : ℤ) (hxk : x ^ l = S.Q ^ k) :
    (∃ ζ : Kˣ, ζ ^ l = 1 ∧
        (QuotientGroup.mk x : Kˣ ⧸ Subgroup.zpowers S.Q) = QuotientGroup.mk ζ)
      ∨ (∃ y : Kˣ, (QuotientGroup.mk x : Kˣ ⧸ Subgroup.zpowers S.Q) = QuotientGroup.mk y
          ∧ 0 < vAdd S.v y ∧ vAdd S.v y < vAdd S.v S.Q) := by
  by_cases hdvd : (l : ℤ) ∣ k
  · -- ★`l ∣ k` なら第 905 がそのまま `μ_l` 型を与える
    exact Or.inl (exists_rootOfUnity_mk_eq S.Q hl.pos x k hxk hdvd)
  · right
    -- ★付値の関係式 `l·vAdd(x) = k·vAdd(Q)`
    have hv : (l : ℤ) * vAdd S.v x = k * vAdd S.v S.Q := by
      have h1 : vAdd S.v (x ^ l) = (l : ℤ) * vAdd S.v x := by
        rw [← zpow_natCast x l, vAdd_zpow]
      have h2 : vAdd S.v (S.Q ^ k) = k * vAdd S.v S.Q := vAdd_zpow S.v S.Q k
      rw [← h1, ← h2, hxk]
    -- ★`l` は素で `l ∤ k` なので `l ∣ vAdd(Q)`
    have hlQ : (l : ℤ) ∣ vAdd S.v S.Q := by
      have hdvd2 : (l : ℤ) ∣ k * vAdd S.v S.Q := ⟨vAdd S.v x, hv.symm⟩
      rcases (Nat.prime_iff_prime_int.1 hl).dvd_mul.1 hdvd2 with h | h
      · exact absurd h hdvd
      · exact h
    obtain ⟨d, hd⟩ := hlQ
    have hlpos : (0 : ℤ) < (l : ℤ) := by exact_mod_cast hl.pos
    have hQpos : 0 < vAdd S.v S.Q := S.hQ
    have hd0 : 0 < d := by nlinarith [hQpos, hd, hlpos]
    have hvx : vAdd S.v x = k * d := by
      rw [hd] at hv
      have hcan : (l : ℤ) * vAdd S.v x = (l : ℤ) * (k * d) := by linarith [hv]
      exact mul_left_cancel₀ (by omega) hcan
    -- ★`k = l·⌊k/l⌋ + (k mod l)` と `0 < k mod l < l`
    have hkr : k = (l : ℤ) * (k / (l : ℤ)) + k % (l : ℤ) :=
      (Int.mul_ediv_add_emod k (l : ℤ)).symm
    have hr0 : 0 < k % (l : ℤ) := by
      rcases (Int.emod_nonneg k (by omega : (l : ℤ) ≠ 0)).lt_or_eq with h | h
      · exact h
      · exact absurd (Int.dvd_of_emod_eq_zero h.symm) hdvd
    have hrl : k % (l : ℤ) < (l : ℤ) := Int.emod_lt_of_pos k hlpos
    -- ★代表 `y := x·Q^{-⌊k/l⌋}` の付値は `(k mod l)·d`
    have hvy : vAdd S.v (x * (S.Q ^ (k / (l : ℤ)))⁻¹) = d * (k % (l : ℤ)) := by
      rw [vAdd_mul, vAdd_inv, vAdd_zpow, hvx, hd]
      nlinarith [hkr]
    refine ⟨x * (S.Q ^ (k / (l : ℤ)))⁻¹, ?_, ?_, ?_⟩
    · refine (QuotientGroup.eq (s := Subgroup.zpowers S.Q)).2 ?_
      have hval : x⁻¹ * (x * (S.Q ^ (k / (l : ℤ)))⁻¹) = (S.Q ^ (k / (l : ℤ)))⁻¹ := by
        rw [← mul_assoc, inv_mul_cancel, one_mul]
      rw [hval]
      exact Subgroup.inv_mem _ (Subgroup.zpow_mem _ (Subgroup.mem_zpowers S.Q) _)
    · rw [hvy]; positivity
    · rw [hvy, hd]; nlinarith [hrl, hd0]

/-! ## ★★★★★★★★★★★★点の側の二者択一 -/

/-- ★★★★★★★★★★★★**位数 `l` の `K`-有理点は、`μ_l` の点か、深い代表を持つ**（第 1410）。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

☆道は第 946 の `exists_mu_point_of_rational` と同じ 3 段だが、
最後の段で `hcop` を使わず `mu_or_deep_of_pow_eq_zpow` を当てる:

1. `Φ` は全射なので `P = tatePhi([x])` と書ける
2. `l • P = 0` から `x^l = Q^k`（第 916）
3. ★二者択一（`μ_l` か深い代表か） -/
theorem mu_or_deep_point (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (P : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) :
    (∃ ζ : Kˣ, ζ ^ l = 1 ∧ P = tatePhi S hΔ (QuotientGroup.mk ζ))
      ∨ (∃ y : Kˣ, P = tatePhi S hΔ (QuotientGroup.mk y)
          ∧ 0 < vAdd S.v y ∧ vAdd S.v y < vAdd S.v S.Q) := by
  -- ★段 1
  obtain ⟨c, hc⟩ := Φ.surjective P
  obtain ⟨x, hx⟩ := QuotientGroup.mk_surjective (s := Subgroup.zpowers S.Q)
    (Additive.toMul c)
  have hPx : tatePhi S hΔ (QuotientGroup.mk x) = P := by
    rw [← hΦ, hx]
    simpa using hc
  -- ★段 2
  obtain ⟨k, hk⟩ := exists_zpow_of_nsmul_tatePhi_eq_zero S hΔ Φ hΦ x l
    (by rw [hPx]; exact hP)
  -- ★段 3
  rcases mu_or_deep_of_pow_eq_zpow S hl x k hk.symm with ⟨ζ, hζl, hζc⟩ | ⟨y, hy, hy0, hy1⟩
  · exact Or.inl ⟨ζ, hζl, by rw [← hPx, hζc]⟩
  · exact Or.inr ⟨y, by rw [← hPx, hy], hy0, hy1⟩

end Rational

/-! ## ★★★★★★★★★★★★★★★★完備 DVR だけで述べた形 -/

section Dvr

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★★★★★★★★★★★★★**完備な DVR だけで、位数 `l` の `K`-有理点は
`μ_l` の点か、深い代表を持つ**（第 1410）。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★★★★**2026-09-02（第 1410）**——引数は `q`・`Δ ≠ 0`・`l` が素・`l • P = 0` の
4 つだけである。☆第 946 の `exists_mu_point_dvr` から **`hcop` が落ちた**
（その代わり結論が二者択一になった）。 -/
theorem mu_or_deep_point_dvr (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime)
    (P : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) :
    (∃ ζ : Kˣ, ζ ^ l = 1 ∧
        P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk ζ))
      ∨ (∃ y : Kˣ, P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk y)
          ∧ 0 < vAdd (mkTateSetup (K := K) q hq hq0).v y
          ∧ vAdd (mkTateSetup (K := K) q hq hq0).v y
            < vAdd (mkTateSetup (K := K) q hq hq0).v (mkTateSetup (K := K) q hq hq0).Q) :=
  mu_or_deep_point (mkTateSetup q hq hq0) hl hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) P hP

end Dvr

/-! ## ★出典の紐付け(`.src`) -/

def mu_or_deep_of_pow_eq_zpow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(x^l = Q^k の二者択一——μ_l 型か深い代表か。★l ∤ v(Q) 不要)",
    sectionId := "genell-lemma-3-2" }

def mu_or_deep_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(位数 l の有理点の二者択一。★l ∤ v(Q) 不要)",
    sectionId := "genell-lemma-3-2" }

def mu_or_deep_point_dvr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(完備 DVR だけで述べた二者択一。★l ∤ v(Q) 不要)",
    sectionId := "genell-lemma-3-2" }

def mu_or_deep_point_dvr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_rootOfUnity_mk_eq(第 905、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_rootOfUnity_mk_eq") 1,
    .citation "[ABC3]" "exists_zpow_of_nsmul_tatePhi_eq_zero(第 916、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_zpow_of_nsmul_tatePhi_eq_zero") 1,
    .citation "[ABC3]" "dvrTatePhiAddEquiv(第 899、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.dvrTatePhiAddEquiv") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1410）**——第 946 の `hcop`（`l ∤ vAdd(Q)`）を落とし、" ++
       "代わりに**二者択一**を出した。" ++
       "☆`l ∤ k` の側では `l ∣ vAdd(Q)` が自動で従い、`vAdd(Q) = l·d` として " ++
       "`y := x·Q^{-⌊k/l⌋}` が `vAdd(y) = (k mod l)·d ∈ (0, vAdd(Q))` を満たす。" ++
       "★★★これは `hcopDoesNotDescend2026_09_02`（底変換で `l ∤ jExp` が壊れる）を" ++
       "**根本から回避する**ための骨である" ++
       "——`vAdd(y) > 0` は `y ∈ 𝔪` を意味し、核の `x` 座標がすべて `𝔪` に入るので " ++
       "`veluV2 ≡ b₄ ≡ 0`、したがって `c₄(veluCurve) = c₄(E_q) + 240v` は単元になる。") 15 ]

end ABC3.Found.GenEll
