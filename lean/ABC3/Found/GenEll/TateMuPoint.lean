/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38PiFromPhi
import ABC3.Meta.Claim

/-!
# 第 1297 ブロック —— **`Φ[ζ]` は位数ちょうど `l` の点**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`μ_l ⊂ E(K₀)` の本体

第 1296 は「`σ` が位数 `l` の点を固定する」ことを要求する。
★Tate 曲線では `μ_l ⊂ E[l]` がそれを与える——`ζ ∈ K₀` なら
`Φ[ζ]` は**基礎体の上の**位数 `l` の点だからである。

☆本ブロックはその位数がちょうど `l` であることを示す:

* `l · Φ[ζ] = Φ[ζˡ] = Φ[1] = 0`
* `Φ[ζ] ≠ 0`——さもなくば `ζ ∈ ⟨Q⟩`、しかし `v(ζ) = 0 < v(Q)` なので `ζ = 1` になる

★★★これで**基礎局所体の上の Tate 一意化だけ**から `μ_l` の点が取れる。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★**`l` 乗根の付値は `0`**——★**無条件**（第 1297）。 -/
theorem vAdd_eq_zero_of_pow_eq_one (S : TateSetup R I K) {l : ℕ} (hl : 0 < l)
    (ζ : Kˣ) (hζl : ζ ^ l = 1) : vAdd S.v ζ = 0 := by
  have h1 : vAdd S.v (ζ ^ l) = 0 := by rw [hζl]; simp [vAdd]
  have h2 : vAdd S.v (ζ ^ l) = (l : ℤ) * vAdd S.v ζ := by
    rw [← zpow_natCast ζ l, vAdd_zpow]
  rw [h2] at h1
  have : (l : ℤ) ≠ 0 := by exact_mod_cast hl.ne'
  exact by
    rcases mul_eq_zero.mp h1 with h | h
    · exact absurd h this
    · exact h

/-- ★★★★★★★★★★★★★★★★
**`Φ[ζ]` は位数ちょうど `l`**——★**無条件**（第 1297）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`v(ζ) = 0 < v(Q)` なので `ζ ∈ ⟨Q⟩` なら `ζ = 1` になってしまう。 -/
theorem addOrderOf_tatePhi_zeta (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    {P : Type*} [AddCommGroup P] (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (ζ : Kˣ) (hζ : IsPrimitiveRoot ((ζ : K)) l) :
    addOrderOf (Φ (Additive.ofMul (QuotientGroup.mk ζ))) = l := by
  have hζl : ζ ^ l = 1 := units_pow_eq_one_of_isPrimitiveRoot hζ
  have hvz : vAdd S.v ζ = 0 := vAdd_eq_zero_of_pow_eq_one S hl.pos ζ hζl
  have hkill : l • (Φ (Additive.ofMul (QuotientGroup.mk ζ))) = 0 := by
    have h1 : l • (Additive.ofMul (QuotientGroup.mk ζ : Kˣ ⧸ Subgroup.zpowers S.Q))
        = Additive.ofMul (QuotientGroup.mk (ζ ^ l)) := by
      rw [← ofMul_pow]
      rfl
    rw [← map_nsmul, h1, hζl]
    have h2 : (QuotientGroup.mk (1 : Kˣ) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
      simp
    rw [h2]
    exact map_zero Φ
  have hne : Φ (Additive.ofMul (QuotientGroup.mk ζ)) ≠ 0 := by
    intro hcon
    have hcls : (QuotientGroup.mk ζ : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
      have h3 : Φ (Additive.ofMul (QuotientGroup.mk ζ)) = Φ 0 := by rw [hcon, map_zero]
      have h4 := Φ.injective h3
      exact h4
    obtain ⟨n, hn⟩ := (QuotientGroup.eq_one_iff _).1 hcls
    -- ζ = Q ^ n, valuations force n = 0
    have hv : vAdd S.v ζ = n * vAdd S.v S.Q := by
      rw [← hn, vAdd_zpow]
    rw [hvz] at hv
    have hQ : 0 < vAdd S.v S.Q := S.hQ
    have hn0 : n = 0 := by
      rcases mul_eq_zero.mp hv.symm with h | h
      · exact h
      · omega
    have hQ0 : S.Q ^ (0 : ℤ) = ζ := by rw [← hn0]; exact hn
    have hζu : ζ = 1 := by simpa using hQ0.symm
    have hζ1 : ((ζ : K)) = 1 := by rw [hζu]; simp
    have h1lt : 1 < l := hl.one_lt
    exact (hζ.ne_one h1lt) hζ1
  have hdvd : addOrderOf (Φ (Additive.ofMul (QuotientGroup.mk ζ))) ∣ l :=
    addOrderOf_dvd_of_nsmul_eq_zero hkill
  rcases (Nat.Prime.eq_one_or_self_of_dvd hl _ hdvd) with h | h
  · exact absurd (AddMonoid.addOrderOf_eq_one_iff.mp h) hne
  · exact h

/-! ## ★出典の紐付け(`.src`) -/

def vAdd_eq_zero_of_pow_eq_one.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 乗根の付値は 0。★無条件)",
    sectionId := "genell-thm-3-8" }

def addOrderOf_tatePhi_zeta.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Φ[ζ] は位数ちょうど l。★無条件)",
    sectionId := "genell-thm-3-8" }

def addOrderOf_tatePhi_zeta.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1297）**——`μ_l ⊂ E(K₀)` の本体である。" ++
       "☆`v(ζ) = 0 < v(Q)` なので `ζ ∈ ⟨Q⟩` なら `ζ = 1` になってしまう。" ++
       "★★★これで**基礎局所体の上の Tate 一意化だけ**から `μ_l` の点が取れる。") 2 ]

end ABC3.Found.GenEll
