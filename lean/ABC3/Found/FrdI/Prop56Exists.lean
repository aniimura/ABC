/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop56Sec
import ABC3.Found.FrdI.Prop56Limit
import ABC3.Found.ProL.OneSubPow

/-!
# `Proposition 5.6` の `u` の存在(鎖 `prop56` の `p56-limit` / `p56`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.106–107。

原文 (FrdI p.107):
> p, as well as a wl ∈O×(A)[p] such that wl · up · φl · u−1

## ★★残っていた 1 段

`Prop56Sec.lean` の `prop_5_6_claim` は
「素数のところで共役が合う `u` が **1 つあれば**、`(σ, φ)` の両方が合う」
を与えていた。残るのは**その `u` の存在**である。

## ★構成(原文との対応)

1. `φ′_l = v_l · φ_l`(`exists_otimes_of_frobType_same_deg`、既済)
2. `φ, φ′` がともに単系準同型であることから `v` の**両立条件**
   `v_{l₁} · v_{l₂}^{l₁} = v_{l₂} · v_{l₁}^{l₂}` が出る(`otimes_compat`、本ファイル)。
   ★消去に使うのは `𝒞` の **totally epimorphic** 性だけである。
3. `u^{1-l} = v_l` を全素数で同時に解く(`exists_forall_pow_one_sub`、`ProL/OneSubPow.lean`)。
4. Frobenius-normalized で `u · φ_l · u⁻¹ = u^{1-l} · φ_l = v_l · φ_l = φ′_l`。

★★原文の「`E_{M_p}` の中で `p≈` を使い、最後に `u_p` の無限積を取る」は、
`ProL/OneSubPow.lean` の冒頭に記した読み替え(積分解の逆像)で置き換えてある。
-/

namespace ABC3.Found.FrdI

open CategoryTheory ABC3.Found.ProL

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★素数を `ℕ≥1` の元として見る。 -/
def pnatOfPrime (l : Nat.Primes) : ℕ+ := ⟨(l : ℕ), l.2.pos⟩

@[simp] theorem pnatOfPrime_coe (l : Nat.Primes) : ((pnatOfPrime l : ℕ+) : ℕ) = (l : ℕ) := rfl

/-! ## ★1. 両立条件 -/

/-- ★★★★**`v` の両立条件** —— `φ`, `φ′` がともに可換(単系準同型の像)なら

  `v₁ · v₂^{d₁} = v₂ · v₁^{d₂}`

(`dᵢ = deg φᵢ`)。★消去に使うのは **totally epimorphic** 性だけである。

原文 (FrdI p.107):
> [cf. Definition 2.8, (ii)], it thus follows that to prove the existence of a u ∈O×(A) -/
theorem otimes_compat (htot : IsTotallyEpimorphic C) {A : C}
    (hfn : IsFrobeniusNormalized P A)
    {phi1 phi2 v1 v2 : End A}
    (hphi1 : IsBaseIdentity P ((phi1 : A ⟶ A))) (hphi2 : IsBaseIdentity P ((phi2 : A ⟶ A)))
    (hv1 : v1 ∈ OTri P A) (hv2 : v2 ∈ OTri P A)
    (hphicomm : phi1 * phi2 = phi2 * phi1)
    (hpsicomm : (v1 * phi1) * (v2 * phi2) = (v2 * phi2) * (v1 * phi1)) :
    v1 * v2 ^ ((P.degFr ((phi1 : A ⟶ A)) : ℕ+) : ℕ)
      = v2 * v1 ^ ((P.degFr ((phi2 : A ⟶ A)) : ℕ+) : ℕ) := by
  set d1 := ((P.degFr ((phi1 : A ⟶ A)) : ℕ+) : ℕ) with hd1
  set d2 := ((P.degFr ((phi2 : A ⟶ A)) : ℕ+) : ℕ) with hd2
  have hL : (v1 * phi1) * (v2 * phi2) = (v1 * v2 ^ d1) * (phi1 * phi2) := by
    calc (v1 * phi1) * (v2 * phi2) = v1 * (phi1 * v2) * phi2 := by simp only [mul_assoc]
      _ = v1 * (v2 ^ d1 * phi1) * phi2 := by rw [frobNorm_end hfn hphi1 hv2]
      _ = (v1 * v2 ^ d1) * (phi1 * phi2) := by simp only [mul_assoc]
  have hR : (v2 * phi2) * (v1 * phi1) = (v2 * v1 ^ d2) * (phi1 * phi2) := by
    calc (v2 * phi2) * (v1 * phi1) = v2 * (phi2 * v1) * phi1 := by simp only [mul_assoc]
      _ = v2 * (v1 ^ d2 * phi2) * phi1 := by rw [frobNorm_end hfn hphi2 hv1]
      _ = (v2 * v1 ^ d2) * (phi2 * phi1) := by simp only [mul_assoc]
      _ = (v2 * v1 ^ d2) * (phi1 * phi2) := by rw [hphicomm]
  rw [hL, hR] at hpsicomm
  haveI : Epi (((phi1 * phi2 : End A)) : A ⟶ A) := htot A A _
  exact (cancel_epi (((phi1 * phi2 : End A)) : A ⟶ A)).mp hpsicomm

/-! ## ★2. `u` の存在 -/

/-- ★★★★★★**[FrdI] Proposition 5.6 —— `u` の存在**。

★★2 つの Frobenius-section `F`, `F′` に対し、`𝒪^×(A)` の元 `u` があって
**すべての素数 `l`** で `u · φ_l · u⁻¹ = φ′_l` となる。

★原文 p.107 の議論(`M_p` を潰して `p` 成分ごとに解き、`u_p` の無限積を取る)を、
`M ≅ ∏_l M[l]`(`decompEquiv`)＋各成分で `x ↦ x^{q-1}` が全単射
(`pow_pred_bijective`)に読み替えたもの。

原文 (FrdI p.107):
> as well as a wl ∈O×(A)[p] such that wl · up · φl · u−1 -/
theorem exists_conj_frobenius [IsConnected D]
    (G : Frobenioid P) (htot : IsTotallyEpimorphic C) (hup : IsOfUnitProfiniteType P)
    {A : C} (hfn : IsFrobeniusNormalized P A)
    {S S2 : BaseSection P} {Fs : ℕ+ →* SectionEnd S} {Fs2 : ℕ+ →* SectionEnd S2}
    (hFs : IsFrobeniusSection S Fs) (hFs2 : IsFrobeniusSection S2 Fs2)
    (hA : S.objP A) (hA2 : S2.objP A) :
    ∃ u ui : End A, u ∈ OTimes P A ∧ ui ∈ OTimes P A ∧ u * ui = 1 ∧ ui * u = 1 ∧
      ∀ l : ℕ+, Nat.Prime (l : ℕ) →
        u * ((Fs l).app ⟨A, hA⟩) * ui = (Fs2 l).app ⟨A, hA2⟩ := by
  classical
  obtain ⟨t, hcomm, htg, hcpt, htd, -⟩ := hup A
  letI := t
  letI : CommGroup ↥(OTimes P A) :=
    { (inferInstance : Group ↥(OTimes P A)) with mul_comm := hcomm }
  haveI : IsTopologicalGroup ↥(OTimes P A) := htg
  haveI : CompactSpace ↥(OTimes P A) := hcpt
  haveI : TotallyDisconnectedSpace ↥(OTimes P A) := htd
  have hdegphi : ∀ n : ℕ+, P.degFr (((Fs n).app ⟨A, hA⟩ : A ⟶ A)) = n := fun n =>
    ((Fs n).deg_eq ⟨A, hA⟩).symm.trans (hFs.degSection n)
  have hdegpsi : ∀ n : ℕ+, P.degFr (((Fs2 n).app ⟨A, hA2⟩ : A ⟶ A)) = n := fun n =>
    ((Fs2 n).deg_eq ⟨A, hA2⟩).symm.trans (hFs2.degSection n)
  -- ★1. `psi_l = v_l · phi_l`
  have hex : ∀ l : Nat.Primes, ∃ v : ↥(OTimes P A),
      (Fs2 (pnatOfPrime l)).app ⟨A, hA2⟩
        = (v : End A) * (Fs (pnatOfPrime l)).app ⟨A, hA⟩ := by
    intro l
    obtain ⟨v, hvmem, hv⟩ := exists_otimes_of_frobType_same_deg G
      ((Fs (pnatOfPrime l)).app ⟨A, hA⟩) ((Fs2 (pnatOfPrime l)).app ⟨A, hA2⟩)
      (hFs.baseIdentity _ ⟨A, hA⟩) (hFs2.baseIdentity _ ⟨A, hA2⟩)
      (hFs.frobType _ ⟨A, hA⟩) (hFs2.frobType _ ⟨A, hA2⟩)
      (by rw [hdegphi, hdegpsi])
    exact ⟨⟨v, hvmem⟩, hv⟩
  choose v hv using hex
  -- ★2. 両立条件
  have hcompat : ∀ l1 l2 : Nat.Primes,
      v l1 * (v l1 ^ (l2 : ℕ))⁻¹ = v l2 * (v l2 ^ (l1 : ℕ))⁻¹ := by
    intro l1 l2
    have hkey : (v l1 : End A) * (v l2 : End A) ^ (l1 : ℕ)
        = (v l2 : End A) * (v l1 : End A) ^ (l2 : ℕ) := by
      have h := otimes_compat (P := P) htot hfn
        (hFs.baseIdentity (pnatOfPrime l1) ⟨A, hA⟩)
        (hFs.baseIdentity (pnatOfPrime l2) ⟨A, hA⟩)
        (OTimes_le_OTri P A (v l1).2) (OTimes_le_OTri P A (v l2).2)
        (by
          have e1 : (Fs (pnatOfPrime l1)).app ⟨A, hA⟩ * (Fs (pnatOfPrime l2)).app ⟨A, hA⟩
              = (Fs (pnatOfPrime l1 * pnatOfPrime l2)).app ⟨A, hA⟩ := by
            rw [map_mul]; rfl
          have e2 : (Fs (pnatOfPrime l2)).app ⟨A, hA⟩ * (Fs (pnatOfPrime l1)).app ⟨A, hA⟩
              = (Fs (pnatOfPrime l2 * pnatOfPrime l1)).app ⟨A, hA⟩ := by
            rw [map_mul]; rfl
          rw [e1, e2, mul_comm])
        (by
          rw [← hv l1, ← hv l2]
          have e1 : (Fs2 (pnatOfPrime l1)).app ⟨A, hA2⟩ * (Fs2 (pnatOfPrime l2)).app ⟨A, hA2⟩
              = (Fs2 (pnatOfPrime l1 * pnatOfPrime l2)).app ⟨A, hA2⟩ := by
            rw [map_mul]; rfl
          have e2 : (Fs2 (pnatOfPrime l2)).app ⟨A, hA2⟩ * (Fs2 (pnatOfPrime l1)).app ⟨A, hA2⟩
              = (Fs2 (pnatOfPrime l2 * pnatOfPrime l1)).app ⟨A, hA2⟩ := by
            rw [map_mul]; rfl
          rw [e1, e2, mul_comm])
      rw [hdegphi, hdegphi] at h
      exact h
    -- ★部分群の中に降ろす
    have hsub : v l1 * v l2 ^ (l1 : ℕ) = v l2 * v l1 ^ (l2 : ℕ) := by
      refine Subtype.ext ?_
      rw [Submonoid.coe_mul, Submonoid.coe_mul, SubmonoidClass.coe_pow, SubmonoidClass.coe_pow]
      exact hkey
    -- ★`a·b^m = b·a^n` から `a·(a^n)⁻¹ = b·(b^m)⁻¹`(可換群)
    refine mul_right_cancel (b := v l1 ^ (l2 : ℕ) * v l2 ^ (l1 : ℕ)) ?_
    calc v l1 * (v l1 ^ (l2 : ℕ))⁻¹ * (v l1 ^ (l2 : ℕ) * v l2 ^ (l1 : ℕ))
        = v l1 * v l2 ^ (l1 : ℕ) := by group
      _ = v l2 * v l1 ^ (l2 : ℕ) := hsub
      _ = v l2 * (v l2 ^ (l1 : ℕ))⁻¹ * (v l1 ^ (l2 : ℕ) * v l2 ^ (l1 : ℕ)) := by
          rw [mul_comm (v l1 ^ (l2 : ℕ)) (v l2 ^ (l1 : ℕ))]; group
  -- ★3. 連立を解く
  obtain ⟨U, hU⟩ := exists_forall_pow_one_sub v hcompat
  refine ⟨(U : End A), ((U⁻¹ : ↥(OTimes P A)) : End A), U.2, (U⁻¹ : ↥(OTimes P A)).2,
    congrArg Subtype.val (mul_inv_cancel U), congrArg Subtype.val (inv_mul_cancel U), ?_⟩
  intro l hl
  have hui : ((U⁻¹ : ↥(OTimes P A)) : End A) ∈ OTri P A :=
    OTimes_le_OTri P A (U⁻¹ : ↥(OTimes P A)).2
  have hconj := conj_frobNorm (P := P) hfn (φ := (Fs l).app ⟨A, hA⟩) (w := (U : End A))
    (wi := ((U⁻¹ : ↥(OTimes P A)) : End A)) (hFs.baseIdentity l ⟨A, hA⟩) hui
  rw [hdegphi] at hconj
  have hlp : pnatOfPrime ⟨(l : ℕ), hl⟩ = l := rfl
  have hstep : (U : End A) * ((U⁻¹ : ↥(OTimes P A)) : End A) ^ ((l : ℕ+) : ℕ)
      = (v ⟨(l : ℕ), hl⟩ : End A) := by
    calc (U : End A) * ((U⁻¹ : ↥(OTimes P A)) : End A) ^ ((l : ℕ+) : ℕ)
        = ((U * (U⁻¹ : ↥(OTimes P A)) ^ ((l : ℕ+) : ℕ) : ↥(OTimes P A)) : End A) := by
          rw [Submonoid.coe_mul, SubmonoidClass.coe_pow]
      _ = (v ⟨(l : ℕ), hl⟩ : End A) := by
          congr 1
          rw [inv_pow]
          exact hU ⟨(l : ℕ), hl⟩
  rw [hconj, hstep]
  have := hv ⟨(l : ℕ), hl⟩
  rw [hlp] at this
  exact this.symm

/-! ## ★3. `Proposition 5.6` 本体 -/

/-- ★★★★★★**[FrdI] Proposition 5.6**(Base-Sections of Frobenius-Trivial Objects)。

★★2 つの base-Frobenius pair `(P, F)`、`(P′, F′)` が `A` の上に定める対 `(σ, φ)` は、
**`𝒪^×(A)` の 1 つの元による共役**の違いを除いて一致する ——
すなわち `(C, A)` だけで定まり `(P, F)` に依らない。

| 段 | 実装 |
|---|---|
| `u` の存在 | `exists_conj_frobenius`(本ファイル) |
| `u` から `(σ, φ)` の両方へ | `prop_5_6_claim`(`Prop56Sec.lean`) |

原文 (FrdI p.105):
> — where σ is a group homomorphism whose composite with the natural surjec- -/
theorem prop_5_6 [IsConnected D]
    (G : Frobenioid P) (htot : IsTotallyEpimorphic C) (hup : IsOfUnitProfiniteType P)
    {A : C} (hfn : IsFrobeniusNormalized P A)
    {S S2 : BaseSection P} {Fs : ℕ+ →* SectionEnd S} {Fs2 : ℕ+ →* SectionEnd S2}
    (hFs : IsFrobeniusSection S Fs) (hFs2 : IsFrobeniusSection S2 Fs2)
    (hA : S.objP A) (hA2 : S2.objP A) :
    ∃ u ui : End A, u ∈ OTimes P A ∧ ui ∈ OTimes P A ∧ u * ui = 1 ∧ ui * u = 1 ∧
      (∀ n : ℕ+, u * ((Fs n).app ⟨A, hA⟩) * ui = (Fs2 n).app ⟨A, hA2⟩) ∧
      (∀ α : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base, IsIso α →
        u * S.sigma hA α * ui = S2.sigma hA2 α) := by
  obtain ⟨u, ui, hu, hui, huui, huiu, hconj⟩ :=
    exists_conj_frobenius G htot hup hfn hFs hFs2 hA hA2
  refine ⟨u, ui, hu, hui, huui, huiu, ?_⟩
  exact prop_5_6_claim htot hfn hFs hA hA2 (OTimes_le_OTri P A hu) (OTimes_le_OTri P A hui)
    huui huiu (fun p _ => hconj p ‹_›)

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Proposition 5.6`(条なし)。

| 条 | 実装 |
|---|---|
| `u` の存在 | `exists_conj_frobenius` |
| 両立条件 | `otimes_compat` |
| `φ` 側 | `prop_5_6_phi_conj`(`Prop56Sec.lean`) |
| `σ` 側 | `prop_5_6_sigma_conj`(`Prop56Sec.lean`) |
| `M_p` | `MpDiv`(`Prop56MpAbs.lean`、Nikolov–Segal 抜き) | -/
def prop_5_6.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.6",
    sectionId := "frdi-prop-5-6" }

/-- ★★locator —— `Proposition 5.6` の両立条件。 -/
def otimes_compat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 107,
    item := "Proposition 5.6 — v の両立条件(φ と φ' がともに準同型)",
    sectionId := "frdi-prop-5-6" }

/-- ★★★locator —— `Proposition 5.6` の `u` の存在。 -/
def exists_conj_frobenius.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 107,
    item := "Proposition 5.6 — u の存在(全素数で同時に共役)",
    sectionId := "frdi-prop-5-6" }

end ABC3.Found.FrdI
