/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.ProL.Decomposition
import ABC3.Found.ProL.ZetaPow
import ABC3.Found.ProL.DivByP.Definition28

/-!
# DivByP —— `[FrdI] Proposition 5.6` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.ProL

universe u
variable {N : Type u} [CommGroup N] [TopologicalSpace N] [IsTopologicalGroup N]
  [CompactSpace N] [TotallyDisconnectedSpace N]

/-! ## ★★`Proposition 5.6` の `u` の構成で使う 2 本

原文 (FrdI p.107):
> — which [by the total epimorphicity of C] implies that wl

★原文は `w_l ≈ w_l^p` から `w_l ≈ 1` を出し、
`u · φ_p · u⁻¹ ≈ u^{1-p} · φ_p` を解いて `u_p` を作る。
★どちらも「`p-1` は `p` と素」だけで出る。 -/

theorem coprime_pred_self {p : ℕ} (hp : p.Prime) : Nat.Coprime (p - 1) p := by
  have h2 := hp.two_le
  have hs : p = (p - 1) + 1 := by omega
  rw [hs]
  simp [Nat.Coprime]

/-- ★★pro-`p` 群では `w = w^p` なら `w = 1`。 -/
theorem eq_one_of_eq_pow_self {p : ℕ} (hp : p.Prime)
    (hpro : ∀ U : OpenNormalSubgroup N, IsPGroup p (N ⧸ U.toSubgroup))
    {w : N} (hw : w = w ^ p) : w = 1 := by
  have hinj := pow_injective_of_proL hp hpro (coprime_pred_self hp)
  refine hinj ?_
  show w ^ (p - 1) = (1 : N) ^ (p - 1)
  rw [one_pow]
  have h2 : w ^ (p - 1) * w = w ^ p := by
    rw [← pow_succ]
    congr 1
    have := hp.two_le
    omega
  rw [← hw] at h2
  calc w ^ (p - 1) = w ^ (p - 1) * w * w⁻¹ := by group
    _ = w * w⁻¹ := by rw [h2]
    _ = 1 := mul_inv_cancel w

/-- ★★pro-`p` 群では `u ↦ u^{p-1}` は全単射(`p-1` は `p` と素)。 -/
theorem pow_pred_bijective {p : ℕ} (hp : p.Prime)
    (hpro : ∀ U : OpenNormalSubgroup N, IsPGroup p (N ⧸ U.toSubgroup)) :
    Function.Bijective (fun x : N => x ^ (p - 1)) :=
  pow_bijective_of_proL hp hpro (coprime_pred_self hp)

/-- ★★★pro-`p` 群では `u · (u^l)⁻¹` を任意の値に取れる（`l-1` が `p` と素なとき）。

★`u · (u^l)⁻¹ = v` ⇔ `u^{l-1} = v⁻¹`。原文 p.107 の `u^{1-p} = v_p` を解く段。 -/
theorem exists_pow_one_sub {p l : ℕ} (hp : p.Prime)
    (hpro : ∀ U : OpenNormalSubgroup N, IsPGroup p (N ⧸ U.toSubgroup))
    (hcop : Nat.Coprime (l - 1) p) (hl : 1 ≤ l) (v : N) :
    ∃ u : N, u * (u ^ l)⁻¹ = v := by
  obtain ⟨u, hu⟩ := (pow_bijective_of_proL hp hpro hcop).2 v⁻¹
  refine ⟨u, ?_⟩
  have hu' : u ^ (l - 1) = v⁻¹ := hu
  have hsplit : u ^ l = u ^ (l - 1) * u := by
    rw [← pow_succ]
    congr 1
    omega
  rw [hsplit, hu', mul_inv_rev, inv_inv]
  calc u * (u⁻¹ * v) = (u * u⁻¹) * v := by rw [mul_assoc]
    _ = v := by rw [mul_inv_cancel, one_mul]

variable {M : Type u} [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
  [CompactSpace M] [TotallyDisconnectedSpace M]

/-- ★`lPart M l` は pro-`l` である(部分群の言葉で)。 -/
theorem lPart_isProL (l : ℕ) :
    ∀ U : OpenNormalSubgroup ↥(lPart M l), IsPGroup l (↥(lPart M l) ⧸ U.toSubgroup) := by
  have h := isProL_lPartGrp (M := M) l
  rw [isProL_iff] at h
  exact h

/-- ★★★★★**`p` のすべての冪で割り切れる ⟺ `p` 成分が自明**。

★★これが原文の `M_p`(`l ≠ p` の pro-`l` 部分が生成する閉部分群)と、
**抽象的な定義 `⋂_k M^{p^k}` との一致**である。
★左辺は冪しか使わないので、**任意の抽象自己同型で保たれる** ——
すなわち `Proposition 5.6` の `M_p` の共役安定性に
**Nikolov–Segal は要らない**。 -/
theorem forall_pow_pow_iff_pComp_eq_one (p : Nat.Primes) (y : M) :
    (∀ k : ℕ, ∃ x : M, x ^ ((p : ℕ) ^ k) = y) ↔ (decompEquiv M y) p = 1 := by
  constructor
  · intro h
    refine eq_one_of_forall_pow_pow p.2 (lPart_isProL (M := M) (p : ℕ)) ?_
    intro k
    obtain ⟨x, hx⟩ := h k
    refine ⟨(decompEquiv M x) p, ?_⟩
    have h2 : (decompEquiv M (x ^ ((p : ℕ) ^ k))) p = (decompEquiv M y) p := by rw [hx]
    rw [map_pow] at h2
    exact h2
  · intro h k
    have hroot : ∀ l : Nat.Primes,
        ∃ z : ↥(lPart M l.1), z ^ ((p : ℕ) ^ k) = (decompEquiv M y) l := by
      intro l
      by_cases hlp : l = p
      · subst hlp
        exact ⟨1, by rw [one_pow, h]⟩
      · have hne : (l : ℕ) ≠ (p : ℕ) := fun hc => hlp (Subtype.ext hc)
        have hcop : Nat.Coprime ((p : ℕ) ^ k) l.1 :=
          Nat.Coprime.pow_left k (Nat.Coprime.symm ((Nat.coprime_primes l.2 p.2).mpr hne))
        obtain ⟨z, hz⟩ :=
          (pow_bijective_of_proL l.2 (lPart_isProL (M := M) l.1) hcop).2 ((decompEquiv M y) l)
        exact ⟨z, hz⟩
    choose z hz using hroot
    refine ⟨(decompEquiv M).symm z, ?_⟩
    refine (decompEquiv M).injective ?_
    rw [map_pow, ContinuousMulEquiv.apply_symm_apply]
    funext l
    exact hz l

/-- ★★locator —— `Proposition 5.6` の `M_p` の抽象群としての特徴づけ。 -/
def forall_pow_pow_iff_pComp_eq_one.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 106,
    item := "Proposition 5.6 — M_p の抽象群としての特徴づけ",
    sectionId := "frdi-prop-5-6" }

end ABC3.Found.ProL
