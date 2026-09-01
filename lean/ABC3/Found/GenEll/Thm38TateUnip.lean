/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38ZetaPi
import ABC3.Meta.Claim

/-!
# 第 1272 ブロック —— **Tate 一意化の `σ` は `E[l]` に幂単かつ非自明に作用する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——II 側の**心臓**

第 1270（`T_l E` → `E[l]`）と第 1271（埋め込みで運ぶ）で配管は揃った。
★残っていたのは**「`σ` が `E[l]` に `α` として作用する」の中身**である。

☆本ブロックはそれを**完全に抽象な形**で取る:

    Φ : Additive (Kˣ ⧸ ⟨Q⟩) ≃+ P     （Tate 一意化）
    τ : P →+ P                        （`σ` の点の側の作用）
    hτ : τ (Φ [x]) = Φ [σ x]          （同変性——`tatePhi_pointMap`、在庫）

このとき `σ(ζ) = ζ`・`σ(π) = ζπ` なら

| # | 結論 | 中身 |
|---|---|---|
| 1 | `τ²p + p = 2τp`（`l·p = 0` のとき） | **幂単** |
| 2 | `∃ p, l·p = 0 ∧ τp ≠ p` | **非自明** |

★★★これが `alpha_mem_map_of_galTate`（第 1237）が要る 2 条件の**局所での中身**である。
☆機構は第 1174（`l`-捩れは `[ζᵃπᵇ]`）と `α = (1 1 / 0 1)` の指数の動きだけである。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★**単数群の準同型としての `α` の作用**——★**無条件**（第 1272）。

☆`sigma_acts_as_alpha`（第 993）の `Kˣ →* Kˣ` 版である。 -/
theorem sigmaU_zeta_pi (σU : Kˣ →* Kˣ) {ζ π : Kˣ} (hσζ : σU ζ = ζ) (hσπ : σU π = ζ * π)
    (a b : ℤ) : σU (ζ ^ a * π ^ b) = ζ ^ (a + b) * π ^ b := by
  rw [map_mul, map_zpow, map_zpow, hσζ, hσπ, mul_zpow, zpow_add, mul_assoc]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**Tate 一意化の `σ` は `E[l]` で幂単**——★**無条件**（第 1272）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`l`-捩れ点は `Φ[ζᵃπᵇ]`（第 1174）。`σ` は指数を `(a,b) ↦ (a+b,b)` で動かすので、
`σ²` は `(a,b) ↦ (a+2b,b)`。したがって `σ²p + p = 2σp` である。

★★★これが `alpha_mem_map_of_galTate`（第 1237）の `h2` の中身である。 -/
theorem tate_unipotent_of_sigma (S : TateSetup R I K) {l : ℕ} [NeZero l]
    {ζ π : Kˣ} (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    {P : Type*} [AddCommGroup P] (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (τ : P →+ P) (σU : Kˣ →* Kˣ)
    (hτ : ∀ x : Kˣ, τ (Φ (Additive.ofMul (QuotientGroup.mk x)))
      = Φ (Additive.ofMul (QuotientGroup.mk (σU x))))
    (hσζ : σU ζ = ζ) (hσπ : σU π = ζ * π)
    (p : P) (hp : l • p = 0) :
    τ (τ p) + p = τ p + τ p := by
  obtain ⟨c, hc⟩ := Φ.surjective p
  obtain ⟨x, hx⟩ := QuotientGroup.mk_surjective (Additive.toMul c)
  have hcx : Φ (Additive.ofMul (QuotientGroup.mk x)) = p := by
    rw [hx]; exact hc
  have hp' : l • Φ (Additive.ofMul (QuotientGroup.mk x)) = 0 := by rw [hcx]; exact hp
  obtain ⟨a, b, hab⟩ := tate_torsion_eq_phi_zeta_pi S hζ hπl Φ x hp'
  have hpz : p = Φ (Additive.ofMul (QuotientGroup.mk (ζ ^ a * π ^ b))) := by
    rw [← hcx, hab]
  have h1 : τ p = Φ (Additive.ofMul (QuotientGroup.mk (ζ ^ (a + b) * π ^ b))) := by
    rw [hpz, hτ, sigmaU_zeta_pi σU hσζ hσπ]
  have h2 : τ (τ p) = Φ (Additive.ofMul (QuotientGroup.mk (ζ ^ (a + b + b) * π ^ b))) := by
    rw [h1, hτ, sigmaU_zeta_pi σU hσζ hσπ]
  have hz : ζ ^ (a + b + b) * ζ ^ a = ζ ^ (a + b) * ζ ^ (a + b) := by
    rw [← zpow_add, ← zpow_add]
    have he : a + b + b + a = a + b + (a + b) := by ring
    rw [he]
  have key : (ζ ^ (a + b + b) * π ^ b) * (ζ ^ a * π ^ b)
      = (ζ ^ (a + b) * π ^ b) * (ζ ^ (a + b) * π ^ b) := by
    rw [mul_mul_mul_comm (ζ ^ (a + b + b)) (π ^ b) (ζ ^ a) (π ^ b),
      mul_mul_mul_comm (ζ ^ (a + b)) (π ^ b) (ζ ^ (a + b)) (π ^ b), hz]
  rw [h2, h1, hpz, ← Φ.map_add, ← Φ.map_add]
  refine congrArg Φ ?_
  rw [← ofMul_mul, ← ofMul_mul, ← QuotientGroup.mk_mul, ← QuotientGroup.mk_mul, key]

/-- ★★★★★★★★★★★★★★★★★★★★
**Tate 一意化の `σ` は `E[l]` で非自明**——★**無条件**（第 1272）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆証人は `p = Φ[π]` である——`l·p = Φ[πˡ] = Φ[Q] = 0` で、
`τp = Φ[ζπ] ≠ Φ[π]`（さもなくば `ζ ∈ ⟨Q⟩` となり第 1161 に反する）。

★★★これが `alpha_mem_map_of_galTate`（第 1237）の `h1` の中身である。 -/
theorem tate_exists_ne_of_sigma (S : TateSetup R I K) {l : ℕ} [NeZero l] (hl : 1 < l)
    {ζ π : Kˣ} (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    {P : Type*} [AddCommGroup P] (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (τ : P →+ P) (σU : Kˣ →* Kˣ)
    (hτ : ∀ x : Kˣ, τ (Φ (Additive.ofMul (QuotientGroup.mk x)))
      = Φ (Additive.ofMul (QuotientGroup.mk (σU x))))
    (hσπ : σU π = ζ * π) :
    ∃ p : P, l • p = 0 ∧ τ p ≠ p := by
  refine ⟨Φ (Additive.ofMul (QuotientGroup.mk π)), ?_, ?_⟩
  · have h1 : l • (Additive.ofMul (QuotientGroup.mk π : Kˣ ⧸ Subgroup.zpowers S.Q))
        = Additive.ofMul (QuotientGroup.mk (π ^ l)) := by
      rw [← ofMul_pow]
      rfl
    have hQ1 : (QuotientGroup.mk S.Q : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 :=
      (QuotientGroup.eq_one_iff _).2 (Subgroup.mem_zpowers S.Q)
    rw [← map_nsmul, h1, hπl, hQ1]
    exact map_zero Φ
  · intro hcon
    rw [hτ π, hσπ] at hcon
    have heq : (QuotientGroup.mk (ζ * π) : Kˣ ⧸ Subgroup.zpowers S.Q)
        = QuotientGroup.mk π := Φ.injective hcon
    have hz : (QuotientGroup.mk ζ : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
      have := congrArg (fun z : Kˣ ⧸ Subgroup.zpowers S.Q =>
        z * (QuotientGroup.mk π)⁻¹) heq
      simpa [← QuotientGroup.mk_mul, ← QuotientGroup.mk_inv, mul_inv_cancel_right] using this
    obtain ⟨n, hn⟩ := (QuotientGroup.eq_one_iff _).1 hz
    have hmem : ∃ n : ℤ, ζ ^ (1 : ℤ) * π ^ (0 : ℤ) = S.Q ^ n := ⟨n, by simpa using hn.symm⟩
    have hdvd := (zeta_pi_mem_zpowers_iff (by omega : 0 < l)
      (units_pow_eq_one_of_isPrimitiveRoot hζ)
      (units_pow_ne_one_of_isPrimitiveRoot hζ) hπl
      (tateSetup_Q_zpow_eq_one S) 1 0).1 hmem
    have : (l : ℤ) ∣ 1 := hdvd.1
    have hle : (l : ℤ) ≤ 1 := Int.le_of_dvd one_pos this
    omega

/-! ## ★出典の紐付け(`.src`) -/

def sigmaU_zeta_pi.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(単数群の準同型としての α の作用。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_unipotent_of_sigma.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 一意化の σ は E[l] で幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_exists_ne_of_sigma.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 一意化の σ は E[l] で非自明。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_unipotent_of_sigma.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_torsion_eq_phi_zeta_pi(第 1174、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.tate_torsion_eq_phi_zeta_pi") 1,
    .citation "[ABC3]" "zeta_pi_mem_zpowers_iff(第 1162、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.zeta_pi_mem_zpowers_iff") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1272）**——`alpha_mem_map_of_galTate`（第 1237）が要る" ++
       "2 条件の**局所での中身**である。☆同変性 `hτ` は `tatePhi_pointMap`（在庫）が与える。" ++
       "★★★これで II 側に残るのは、Tate 曲線を実際の悪い素点の `E` に" ++
       "当てはめる段（変数変換と完備化）と、大域へ運ぶ段（第 1271）だけになった。") 3 ]

end ABC3.Found.GenEll
