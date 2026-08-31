/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35Concrete

/-!
# 第 901 ブロック —— **★★★★★★★★★★★★★★★★★★★★`Lemma 3.5` に等号は要らない
——不等号で足り、良い素点の義務が消える**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か

これまで `Lemma 3.5` の入力は

    `deg∞(E′) = l·deg∞(E)`（**等号**）

だと思っていた。等号を出すには**良い還元の素点**でも
`v_p(Δ_min(E′)) = 0` を言わねばならず、それは
Néron–Ogg–Shafarevich（同種な曲線は同じ還元型をもつ）を要求する——mathlib に無い。

★★★★**しかし証明を見ると、使っているのは片側だけである**:

    `l·deg∞(E) ≤ deg∞(E′) ≤ 12(1+ϵ)·ht^Falt(E′) + C`

☆したがって局所で要るのは `l·v_p(Δ_min(E)) ≤ v_p(Δ_min(E′))` だけであり、
**良い素点では左辺が `0`** なので `minDeltaExp_nonneg` だけで自動的に成り立つ。

★★★★★★★★**これで `jExp_velu_good`（Néron–Ogg–Shafarevich）の義務が消える**。

| 定理 | 内容 |
|---|---|
| `degInfOf_ge_of_local` | ★局所の不等式を `deg∞` に足し上げる |
| `lemma_3_5_of_isogeny_estimate_le` | ★不等号で受ける組み立て |
| `lemma_3_5_velu_le` | ★`hfalt` を外した形（不等号） |
| `lemma_3_5_velu_local_le` | ★局所の不等式だけを受ける形 |
| `minDeltaExp_le_of_jExp_bad` | ★★**悪い素点だけで足りる** |
| `lemma_3_5_velu_bad_only` | ★★★★★★★★★★**悪い素点の `v_p(j′) = l·v_p(j)` だけ** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★★局所の不等式を足し上げる -/

/-- ★★★★★★★**局所の不等式を大域の `deg∞` へ足し上げる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`log N(p) ≥ log 2 > 0` なので各項で向きが保たれる。 -/
theorem degInfOf_ge_of_local (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ)
    (hloc : ∀ p : HeightOneSpectrum (𝓞 L),
      (l : ℤ) * minDeltaExp p E ≤ minDeltaExp p E') :
    (l : ℝ) * degInfOf L E ≤ degInfOf L E' := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by
    exact_mod_cast Module.finrank_pos
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  have hΔ' : E'.Δ ≠ 0 := E'.isUnit_Δ.ne_zero
  -- ★各項の不等式
  have hterm : ∀ p : HeightOneSpectrum (𝓞 L),
      (l : ℝ) * ((minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (minDeltaExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
    intro p
    have hlog : (0:ℝ) ≤ Real.log (Ideal.absNorm p.asIdeal) :=
      le_trans (Real.log_nonneg (by norm_num)) (log_two_le_log_absNorm p)
    have hcast : ((l : ℝ)) * (minDeltaExp p E : ℝ) ≤ (minDeltaExp p E' : ℝ) := by
      exact_mod_cast hloc p
    calc (l : ℝ) * ((minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        = ((l : ℝ) * (minDeltaExp p E : ℝ)) * Real.log (Ideal.absNorm p.asIdeal) := by ring
      _ ≤ (minDeltaExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal) :=
        mul_le_mul_of_nonneg_right hcast hlog
  -- ★有限台
  have hfin : (Function.support
      (fun p : HeightOneSpectrum (𝓞 L) => (l : ℝ)
        * ((minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)))).Finite := by
    refine Set.Finite.subset (hasFiniteSupport_degInf E hΔ) ?_
    intro p hp
    simp only [Function.mem_support, ne_eq] at hp ⊢
    intro h0
    exact hp (by rw [h0, mul_zero])
  have hfin' := hasFiniteSupport_degInf E' hΔ'
  -- ★和の不等式
  have hsum : ∑ᶠ p : HeightOneSpectrum (𝓞 L),
      (l : ℝ) * ((minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
      ≤ ∑ᶠ p : HeightOneSpectrum (𝓞 L),
        (minDeltaExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal) :=
    finsum_le_finsum' hfin hfin' hterm
  rw [← mul_finsum, ← finrank_mul_degInfOf, ← finrank_mul_degInfOf] at hsum
  nlinarith [hsum, hd]

/-! ## ★★★★★★不等号で受ける組み立て -/

/-- ★★★★★★**[GenEll] Lemma 3.5 —— `deg∞` を不等号で受ける形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★証明が使っているのは `l·deg∞(E) ≤ deg∞(E′)` だけである。 -/
theorem lemma_3_5_of_isogeny_estimate_le (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ),
      (∀ p, SemistableAt p E') →
      (l : ℝ) * degInfOf L E ≤ degInfOf L E' →
      htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := prop_3_4_chain_semistable eps heps
  have hpos : (0:ℝ) < 12 * (1 + eps) := by linarith
  refine ⟨C / (12 * (1 + eps)), fun L _ _ E E' _ _ l hss hdeg hfalt => ?_⟩
  obtain ⟨h1, h2, _⟩ := hC L E' hss
  have hkey : degInfOf L E' ≤ 12 * (1 + eps) * htFaltOf L E' + C := le_trans h1 h2
  have hstep : (l : ℝ) * degInfOf L E
      ≤ 12 * (1 + eps) * (htFaltOf L E + 2 * Real.log l) + C := by
    refine le_trans (le_trans hdeg hkey) ?_
    have := mul_le_mul_of_nonneg_left hfalt hpos.le
    linarith
  rw [show (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
      = ((l : ℝ) * degInfOf L E) / (12 * (1 + eps)) by ring,
    div_le_iff₀ hpos]
  have hC' : C / (12 * (1 + eps)) * (12 * (1 + eps)) = C := by field_simp
  nlinarith [hstep, hC']

/-! ## ★★★★★★★★★★悪い素点だけで足りる -/

/-- ★★**良い素点では不等式は自動的に成り立つ**。

☆半安定で `v_p(j) ≥ 0` なら `v_p(Δ_min(E)) = 0` なので左辺は `0`、
右辺は `minDeltaExp_nonneg` で非負である。
★ここに Néron–Ogg–Shafarevich は要らない。 -/
theorem minDeltaExp_le_of_jExp_bad (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : SemistableAt p E) (hss' : SemistableAt p E') (l : ℕ)
    (hbad : jExp p E < 0 → jExp p E' = (l : ℤ) * jExp p E) :
    (l : ℤ) * minDeltaExp p E ≤ minDeltaExp p E' := by
  rcases lt_or_ge (jExp p E) 0 with h | h
  · rw [minDeltaExp_eq_maxJ_of_semistable p E' hss', minDeltaExp_eq_maxJ_of_semistable p E hss,
      hbad h, ← mul_neg, max_zero_mul _ (by positivity)]
  · have h2 : minDeltaExp p E = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E hss]
      exact max_eq_left (by omega)
    rw [h2, mul_zero]
    exact minDeltaExp_nonneg p E'

/-! ## ★★★★★★★★★★★★★★★★★★★★連鎖——悪い素点の `v_p(j′) = l·v_p(j)` だけ -/

/-- ★★★★★★**[GenEll] Lemma 3.5 —— `hfalt` を外した不等号の形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let -/
theorem lemma_3_5_velu_le (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E') →
      (l : ℝ) * degInfOf L E ≤ degInfOf L E' →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_of_isogeny_estimate_le eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hss hdeg => ?_⟩
  exact hC L E E' l hss hdeg
    (htFalt_veluQuotientFull_le E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint)

/-- ★★★★★★★★**[GenEll] Lemma 3.5 —— 局所の不等式だけを受ける形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let -/
theorem lemma_3_5_velu_local_le (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L),
        (l : ℤ) * minDeltaExp p E ≤ minDeltaExp p E') →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_le eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hss hloc => ?_⟩
  exact hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hss
    (degInfOf_ge_of_local E E' l hloc)

/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 残る入力は**悪い素点だけ**の `v_p(j′) = l·v_p(j)`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-08-31（第 902）**——以前の最終形 `lemma_3_5_velu_j` は
`hgood : ∀ p, 0 ≤ v_p(j) → 0 ≤ v_p(j′)` も受けていた。
それは**同種な曲線は同じ還元型をもつ**（Néron–Ogg–Shafarevich、mathlib に無い）
を要求していた。

☆しかし `Lemma 3.5` の証明が使うのは `l·deg∞(E) ≤ deg∞(E′)` の**片側だけ**であり、
良い素点では左辺が `0` なので自動的に成り立つ。
★★★★**したがって Néron–Ogg–Shafarevich の義務は消えた**。 -/
theorem lemma_3_5_velu_bad_only (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L),
        jExp p E < 0 → jExp p E' = (l : ℤ) * jExp p E) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_local_le eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hssE hssE' hbad => ?_⟩
  exact hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hssE'
    (fun p => minDeltaExp_le_of_jExp_bad p E E' (hssE p) (hssE' p) l (hbad p))

/-! ## ★★★★★★★★★★★★`Δ_min` の言葉で直接受ける形 -/

/-- ★★**悪い素点だけの `Δ_min` の関係から不等式へ**。

☆良い素点（`v_p(j) ≥ 0`）では半安定性から `v_p(Δ_min(E)) = 0` なので自動的。 -/
theorem minDeltaExp_le_of_bad_delta (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : SemistableAt p E) (l : ℕ)
    (hbad : jExp p E < 0 → minDeltaExp p E' = l * minDeltaExp p E) :
    (l : ℤ) * minDeltaExp p E ≤ minDeltaExp p E' := by
  rcases lt_or_ge (jExp p E) 0 with h | h
  · rw [hbad h]
  · have h2 : minDeltaExp p E = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E hss]
      exact max_eq_left (by omega)
    rw [h2, mul_zero]
    exact minDeltaExp_nonneg p E'

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 残る入力は**悪い素点だけ**の
`v_p(Δ_min(E′)) = l·v_p(Δ_min(E))`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`minDeltaExp_eq_mul_of_tateParamR`（第 892）が**そのまま与える形**である。
☆`jExp` を往復する必要がないので、局所の連鎖と直接繋がる。 -/
theorem lemma_3_5_velu_bad_delta (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L),
        jExp p E < 0 → minDeltaExp p E' = l * minDeltaExp p E) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_local_le eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hssE hssE' hbad => ?_⟩
  exact hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hssE'
    (fun p => minDeltaExp_le_of_bad_delta p E E' (hssE p) l (hbad p))

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_5_velu_bad_delta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力は悪い素点だけの Δ_min の関係——局所の連鎖と直接繋がる)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(hfalt を外した不等号の形)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_local_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(局所の不等式だけを受ける形)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_bad_only.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力は悪い素点だけの v_p(j′) = l·v_p(j)——NOS は要らない)",
    sectionId := "genell-lemma-3-5" }

def degInfOf_ge_of_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(局所の不等式を deg∞ に足し上げる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_of_isogeny_estimate_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(deg∞ を不等号で受ける形——等号は要らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_le_of_jExp_bad.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点では不等式は自動的——NOS は要らない。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
