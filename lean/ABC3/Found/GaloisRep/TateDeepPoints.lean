/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluPoints
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Meta.Claim

/-!
# 第 1411 ブロック —— **深い代表の道具立て**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★これは何か——第 1410 の二者択一の「深い側」を使う道具

第 1410 の `mu_or_deep_point` は、位数 `l` の点が **`μ_l` 型**か
**深い代表を持つ**かの二者択一を与えた。
☆本ブロックは後者の場合に必要な道具を並べる。

| 定理 | 内容 |
|---|---|
| `normRep_vAdd_pos_of_not_dvd` | ★★`v(Q) ∤ v(z)` なら正規化代表の付値は正 |
| `mk_pow_ne_one_of_not_dvd` | ★★`v(Q) ∤ i·v(y)` なら `[yⁱ] ≠ 1` |
| `mk_pow_injOn_of_ne_one` | ★★★`[yᵏ] ≠ 1`（`0 < k < l`）から冪の単射性 |
| `pointCoords_tatePhi_of_mem` | ★★★★`a(c) ∈ 𝔪` なら座標は `R` の元の像 |
| `tateCurveAt_a₄_mem` | ★`a₄(E_q) = −5s₃(q) ∈ 𝔪` |
| `veluV2_tateCurveAt_mem` | ★★★★★★**`x, y ∈ 𝔪` なら `v_Q ∈ 𝔪`** |

★★★最後の 1 本が要点である——深い側では核の座標がすべて `𝔪` に入るので
`v = Σ v_Q ∈ 𝔪`、したがって `c₄(veluCurve) = c₄(E_q) + 240v ≡ c₄(E_q)` は単元になる。

☆在庫にあるものはそのまま使う（重複して作らない）:
`tateAOf_mem_of_pos`（第 292）・`tateXpair_mem` / `tateYpair_mem`（第 300 系）・
`isUnit_one_sub`（第 268）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GenEll Finset

/-! ## ★★正規化代表の付値が正であること -/

section Rep

variable {K : Type} [Field K]

/-- ★★**`v(Q) ∤ v(z)` なら `[z]` の正規化代表の付値は正**。

☆`normRep [z] = z·Qⁿ` なので付値は `v(z) + n·v(Q)`。
これが `0` なら `v(Q) ∣ v(z)` である。 -/
theorem normRep_vAdd_pos_of_not_dvd (v : Kˣ →* Multiplicative ℤ)
    (Q : Kˣ) (hQ : 0 < vAdd v Q) (z : Kˣ) (hz : ¬ (vAdd v Q ∣ vAdd v z)) :
    0 < vAdd v (normRep v Q hQ (QuotientGroup.mk z)) := by
  rcases (normRep_nonneg v Q hQ (QuotientGroup.mk z)).lt_or_eq with h | h
  · exact h
  · exfalso
    obtain ⟨n, hn⟩ := (QuotientGroup.eq (s := Subgroup.zpowers Q)).1
      (normRep_mk v Q hQ (QuotientGroup.mk z))
    have hv : n * vAdd v Q = - vAdd v (normRep v Q hQ (QuotientGroup.mk z)) + vAdd v z := by
      have hc := congrArg (vAdd v) hn
      rwa [vAdd_zpow, vAdd_mul, vAdd_inv] at hc
    exact hz ⟨n, by rw [mul_comm]; omega⟩

/-- ★★**`v(Q) ∤ i·v(y)` なら `[yⁱ]` は単位類でない**。 -/
theorem mk_pow_ne_one_of_not_dvd (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ)
    (y : Kˣ) (i : ℕ) (hi : ¬ (vAdd v Q ∣ (i : ℤ) * vAdd v y)) :
    (QuotientGroup.mk (y ^ i) : Kˣ ⧸ Subgroup.zpowers Q) ≠ 1 := by
  intro h
  obtain ⟨n, hn⟩ := (QuotientGroup.eq_one_iff _).1 h
  have hval := congrArg (vAdd v) hn
  rw [vAdd_zpow, show vAdd v (y ^ i) = (i : ℤ) * vAdd v y by
    rw [← zpow_natCast y i, vAdd_zpow]] at hval
  exact hi ⟨n, by rw [mul_comm (vAdd v Q) n]; exact hval.symm⟩

/-- ★★★**`[yᵏ] ≠ 1`（`0 < k < l`）なら `i ↦ [yⁱ]` は `range l` 上で単射**。

☆`[yᵃ] = [yᵇ]`（`b ≤ a`）なら `[y^{a−b}] = 1` である。 -/
theorem mk_pow_injOn_of_ne_one (Q : Kˣ) {l : ℕ} (y : Kˣ)
    (hne : ∀ k : ℕ, 0 < k → k < l →
      (QuotientGroup.mk (y ^ k) : Kˣ ⧸ Subgroup.zpowers Q) ≠ 1)
    {i j : ℕ} (hi : i < l) (hj : j < l)
    (h : (QuotientGroup.mk (y ^ i) : Kˣ ⧸ Subgroup.zpowers Q) = QuotientGroup.mk (y ^ j)) :
    i = j := by
  have key : ∀ a b : ℕ, b ≤ a → a < l →
      (QuotientGroup.mk (y ^ a) : Kˣ ⧸ Subgroup.zpowers Q) = QuotientGroup.mk (y ^ b) →
      a = b := by
    intro a b hba hal hab
    by_contra hne'
    refine hne (a - b) (by omega) (by omega) ?_
    have hsplit : y ^ (a - b) * y ^ b = y ^ a := by
      rw [← pow_add]
      congr 1
      omega
    have hmem := (QuotientGroup.eq (s := Subgroup.zpowers Q)).1 hab
    have heq : y ^ (a - b) = ((y ^ a)⁻¹ * y ^ b)⁻¹ := by
      rw [mul_inv_rev, inv_inv, ← hsplit, mul_comm (y ^ (a - b)) (y ^ b),
        ← mul_assoc, inv_mul_cancel, one_mul]
    refine (QuotientGroup.eq_one_iff (y ^ (a - b))).2 ?_
    rw [heq]
    exact Subgroup.inv_mem _ hmem
  rcases le_total j i with hle | hle
  · exact key i j hle hi h
  · exact (key j i hle hj h.symm).symm

end Rep

/-! ## ★★★★深い類の点の座標 -/

section Coords

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

/-- ★★★★**代表が `I` に入っていれば `Φ(c)` の座標は `R` の元の像**。

☆`1 − a(c)` は `a(c) ∈ I` なら単元（`isUnit_one_sub`）なので、
`pointCoords_tatePtPair` がそのまま効く。 -/
theorem pointCoords_tatePhi_of_mem (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {c : Kˣ ⧸ Subgroup.zpowers S.Q} (hc : c ≠ 1) (ha : tateAOf S c ∈ I) :
    pointCoords (tatePhi S hΔ c)
      = ((algebraMap R K (tateXpair (tateAOf S c) (tateWOf S c) S.q S.hq),
          algebraMap R K (tateYpair (tateAOf S c) (tateWOf S c) S.q S.hq)) : K × K) := by
  rw [tatePhi_eq S hΔ hc]
  exact pointCoords_tatePtPair _ _ _ S.hq _ _ _ hΔ (isUnit_one_sub ha)

end Coords

/-! ## ★★★★★★深い座標では Vélu の `v_Q` も `I` に入る -/

section Velu

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★**`a₄(E_q) = −5·s₃(q) ∈ I`**。 -/
theorem tateCurveAt_a₄_mem (q : R) (hq : q ∈ I) : (tateCurveAt q hq).a₄ ∈ I := by
  rw [tateCurveAt_a₄]
  exact evalAdic_mem _ coeff_zero_tateA4 q hq

/-- ★★★★★★**`x, y ∈ I` なら `v_Q = 3x² + a₄ − y ∈ I`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが「深い核なら `c₄` が単元」の核心である
——`v = Σ v_Q ∈ I` なので `c₄(veluCurve) = c₄(E_q) + 240v ≡ c₄(E_q)`。 -/
theorem veluV2_tateCurveAt_mem {q : R} (hq : q ∈ I) {x y : R} (hx : x ∈ I) (hy : y ∈ I) :
    veluV2 (tateCurveAt q hq) x y ∈ I := by
  rw [veluV2, veluGx_tateCurveAt]
  refine Ideal.sub_mem _ (Ideal.add_mem _ (Ideal.mul_mem_left _ _ ?_)
    (tateCurveAt_a₄_mem q hq)) hy
  rw [sq]
  exact Ideal.mul_mem_left _ _ hx

end Velu

/-! ## ★出典の紐付け(`.src`) -/

def normRep_vAdd_pos_of_not_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(v(Q) ∤ v(z) なら正規化代表の付値は正)",
    sectionId := "genell-def-3-3" }

def mk_pow_ne_one_of_not_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(v(Q) ∤ i·v(y) なら [y^i] は単位類でない)",
    sectionId := "genell-def-3-3" }

def mk_pow_injOn_of_ne_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(単位類を避ける冪は range l 上で単射)",
    sectionId := "genell-def-3-3" }

def pointCoords_tatePhi_of_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(代表が I に入る類の点の座標は R の元の像)",
    sectionId := "genell-def-3-3" }

def veluV2_tateCurveAt_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い座標では Vélu の v_Q も I に入る)",
    sectionId := "genell-lemma-3-5" }

def veluV2_tateCurveAt_mem.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tateCurveAt_a₄(第 850 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateCurveAt_a₄") 1,
    .citation "[ABC3]" "veluGx_tateCurveAt(第 850 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluGx_tateCurveAt") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1411）**——第 1410 の二者択一の「深い側」を使う道具を並べた。" ++
       "☆要点は `veluV2_tateCurveAt_mem`——Tate 曲線では `v_Q = 3x² + a₄ − y` で " ++
       "`a₄ = −5s₃(q) ∈ 𝔪` なので、`x, y ∈ 𝔪` なら `v_Q ∈ 𝔪` である。" ++
       "★したがって核の座標がすべて深ければ `v = Σ v_Q ∈ 𝔪` で " ++
       "`c₄(veluCurve) = c₄(E_q) + 240v` は単元になる。") 17 ]

end ABC3.Found.GaloisRep
