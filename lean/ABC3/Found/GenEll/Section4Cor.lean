/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.Section3
import ABC3.Skeleton.GenEll.GaloisImage
import ABC3.Found.GenEll.PrimeConstants
import ABC3.Found.GenEll.PrimesOfSize
import ABC3.Meta.Claim

/-!
# 第 1444 ブロック —— **★★★★★★★★★★★★★★★★
`Corollary 4.3` と `Corollary 4.4` の証明本体を `Found/` へ**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.22–p.23。

原文 (GenEll p.23):
> Corollary 4.4. (Full Galois Actions for Compactly Bounded Subsets)

## ★★これは何か——場所を規約どおりに直しただけである

`Corollary 4.3` / `Corollary 4.4` の証明は 2026-08-26（第 367）に `sorry` 0 になっていたが、
本体が `Skeleton/GenEll/Section4.lean` に置かれたままであった。
☆本プロジェクトの規約は

* `Skeleton/` = **原典の主張を写す場所**（`sorry` があってよい）
* `Found/` = **`sorry` の無い実装**

であり、`Lemma 4.1` / `Lemma 4.2` は既にその形（`Skeleton` は statement を写して
`Found` の定理へ委譲する）になっている。
★本ファイルは `Corollary 4.3` / `Corollary 4.4` を**同じ形に揃える**。

## ★★★★中身は 1 文字も変えていない

`#print axioms` は移動の前後どちらでも
`[propext, Classical.choice, Quot.sound]` である（`sorryAx` は無い）。

## ☆何に依っているか

* `theorem_3_8`（`Skeleton/GenEll/GaloisImage.lean`、`sorry` 0）
* `prop_3_4`（`Skeleton/GenEll/Section3.lean`）
* `EllModuliData`（`Interface/GenEll/EllModuli.lean` の欄）
* 数値の鎖 `cor43_numeric` / `cor44_numeric`（`Found/GenEll/PrimeConstants.lean`）

★★**`EllModuliData` の欄そのものは未構成である**——
`Interface` で受けている以上、これは「`Corollary 4.4` が絶対的に証明された」ことを意味しない。
☆欄を実際に埋めるのが `Skeleton/GenEll/EllModuliWitness.lean` の仕事であり、
そこには `VeluSemistable.lean` の `j(E′)` の整性と
`GaloisLocal.lean` の Tate 加群の Galois 像が残っている。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Skeleton.GenEll

/-! ## Corollary 4.3 —— 退化する楕円曲線の完全 Galois 作用 -/

/-- **[GenEll] Corollary 4.3**(Full Galois Actions for Degenerating Elliptic Curves)——
★**証明本体**（第 1444 で `Skeleton` から移した）。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

☆statement の逐語の読みは `Skeleton/GenEll/Section4.lean` の `cor_4_3` にある。 -/
theorem cor_4_3 (D : EllModuliData) (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (S : Finset ℕ), (∀ p ∈ S, p.Prime) →
        D.MinimalField E → D.cls E ∉ Exc → D.HasPotMultRed E →
        ∃ lo lb : ℕ, Nat.Prime lo ∧ Nat.Prime lb ∧ lo ∉ S ∧ lb ∉ S
          ∧ (D.PrimeToMultPrimes E lo ∧ D.PrimeToLocalHeights E lo
              ∧ D.PrimeToMultPrimes E lb ∧ D.PrimeToLocalHeights E lb
              ∧ D.PrimeToRamification E lb)
          ∧ (D.ImageContainsSL2 E lo ∧ D.ImageSurjective E lb)
          ∧ (lo : ℝ) ≤ 23040 * 900 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p)
              + C * (D.degOfDefinition E : ℝ) ^ (1 + eps)
          ∧ (lb : ℝ) ≤ 23040 * 900 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 6 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p)
              + C * (D.degOfDefinition E : ℝ) ^ (1 + eps) := by
  classical
  obtain ⟨C₈, hC₈pos, Exc, hExc, h38⟩ := theorem_3_8 D ∅ D.compactlyBounded_empty eps heps
  obtain ⟨hdegbd, hhtbd, -, -⟩ := prop_3_4 D (1/5 : ℝ) (by norm_num)
  obtain ⟨A₁, hA₁⟩ := hdegbd
  obtain ⟨A₂, hA₂⟩ := hhtbd
  obtain ⟨Bl, hBl⟩ := D.faltingsHeight_bddBelow
  obtain ⟨xeps, Ceps, hxepos, hCepos, hCx, hi1, hi2, hii⟩ :=
    ABC3.Found.GenEll.exists_cond_i_ii 1
  obtain ⟨T, hTprime, hTgt⟩ := ABC3.Found.GenEll.exists_finset_primes_sum_log_gt xeps
  have hxTnn : (0:ℝ) ≤ ∑ p ∈ T, Real.log p :=
    Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
  have hunion : ∀ U V : Finset ℕ,
      (∑ p ∈ U ∪ V, Real.log p) ≤ (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := by
    intro U V
    have hui : (∑ p ∈ U ∪ V, Real.log p) + ∑ p ∈ U ∩ V, Real.log p
        = (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := Finset.sum_union_inter
    have hnn : (0:ℝ) ≤ ∑ p ∈ U ∩ V, Real.log p :=
      Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
    linarith
  refine ⟨23040 * (5 * |A₁ + A₂| + 800 * C₈ + 828 * |Bl|) + 2 * (∑ p ∈ T, Real.log p) + 1,
    by nlinarith [abs_nonneg (A₁ + A₂), abs_nonneg Bl, hC₈pos, hxTnn], Exc, hExc, ?_⟩
  intro E S hSprime _hmin hExcE hpot
  have hd1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) := by exact_mod_cast D.degOfDefinition_pos E
  have hdpos : (0:ℝ) < (D.degOfDefinition E : ℝ) := by linarith
  have hdeps1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) ^ eps := Real.one_le_rpow hd1 heps.le
  have hdd : (D.degOfDefinition E : ℝ) * (D.degOfDefinition E : ℝ) ^ eps
      = (D.degOfDefinition E : ℝ) ^ (1 + eps) := by
    rw [Real.rpow_add hdpos, Real.rpow_one]
  have hdP : (D.degOfDefinition E : ℝ) ≤ (D.degOfDefinition E : ℝ) ^ (1 + eps) := by
    rw [← hdd]; nlinarith
  have hP1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) ^ (1 + eps) := le_trans hd1 hdP
  have hBnn : (0:ℝ) ≤ |Bl| := abs_nonneg Bl
  have hFB : -|Bl| ≤ D.faltingsHeight (D.cls E) := le_trans (neg_abs_le Bl) (hBl (D.cls E))
  obtain ⟨H, hHdef⟩ : ∃ H : ℝ, H = max 0 (23040 * 100 * (D.degOfDefinition E : ℝ)
      * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps)) := ⟨_, rfl⟩
  have hHnn : (0:ℝ) ≤ H := by rw [hHdef]; exact le_max_left _ _
  have hXH : 23040 * 100 * (D.degOfDefinition E : ℝ)
      * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps) ≤ H := by
    rw [hHdef]; exact le_max_right _ _
  have hFBd : (0:ℝ) ≤ (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
      + (D.degOfDefinition E : ℝ) * |Bl| := by nlinarith [hFB, hdpos]
  have hHle : H ≤ 23040 * 100 * ((D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E))
      + 23040 * 100 * (C₈ * (D.degOfDefinition E : ℝ) ^ (1 + eps))
      + 23040 * 100 * ((D.degOfDefinition E : ℝ) * |Bl|) := by
    have hXeq : 23040 * 100 * (D.degOfDefinition E : ℝ)
          * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps)
        = 23040 * 100 * ((D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E))
          + 23040 * 100 * (C₈ * (D.degOfDefinition E : ℝ) ^ (1 + eps)) := by
      rw [← hdd]; ring
    have hdB : (0:ℝ) ≤ (D.degOfDefinition E : ℝ) * |Bl| := mul_nonneg hdpos.le hBnn
    have hCP : (0:ℝ) < C₈ * (D.degOfDefinition E : ℝ) ^ (1 + eps) :=
      mul_pos hC₈pos (Real.rpow_pos_of_pos hdpos _)
    rw [hHdef]
    rcases le_or_gt 0 (23040 * 100 * (D.degOfDefinition E : ℝ)
        * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps)) with hX | hX
    · rw [max_eq_right hX, hXeq]; linarith
    · rw [max_eq_left hX.le]; linarith
  -- ★★`Lemma 4.2` を `h = 23040·d·deg∞` に適用する(原文 p.23)
  have hL42 := lemma_4_2 (D.multCard E) (D.multCard_pos E) (D.multPrime E) (D.multPrime_prime E)
      (D.localHt E) (D.localHt_pos E)
  have hsum := D.sum_localHt_eq E
  have hxbadle : (∑ p ∈ D.badPrimes E, Real.log p)
      ≤ 5/2 * (23040 * ((D.degOfDefinition E : ℝ) * D.degInf (D.cls E))) := by
    have h1 := hL42.1
    have h3 := hL42.2.2
    rw [hsum] at h1 h3
    have hb := D.sum_log_badPrimes_le E
    linarith
  have hAineq : D.degInf (D.cls E)
      ≤ 12 * (1 + 1/5) * D.faltingsHeight (D.cls E) + (A₁ + A₂) := by
    have h1 := hA₁ (D.cls E)
    have h2 := hA₂ (D.cls E)
    simp only at h1 h2
    linarith
  -- ★★★`Lemma 4.1` を `M = 1`･`1 + 6ϵ = 2` で適用する
  have hget : ∀ Bad : Finset ℕ, (∀ p ∈ Bad, Nat.Prime p) →
      ∃ l : ℕ, Nat.Prime l ∧ l ∉ S ∧ l ∉ Bad ∧ H ≤ (l:ℝ) ∧
        (l:ℝ) ≤ 2 * ((∑ p ∈ S, Real.log p) + (∑ p ∈ Bad, Real.log p)
          + (∑ p ∈ T, Real.log p)) + 8 * H := by
    intro Bad hBadp
    have hAllp : ∀ p ∈ S ∪ Bad ∪ T, Nat.Prime p := by
      intro p hp
      rcases Finset.mem_union.1 hp with hp' | hp'
      · rcases Finset.mem_union.1 hp' with hp'' | hp''
        · exact hSprime p hp''
        · exact hBadp p hp''
      · exact hTprime p hp'
    have hge : (∑ p ∈ T, Real.log p) ≤ ∑ p ∈ S ∪ Bad ∪ T, Real.log p :=
      Finset.sum_le_sum_of_subset_of_nonneg Finset.subset_union_right
        (fun p _ _ => Real.log_natCast_nonneg p)
    have hAx : xeps < ∑ p ∈ S ∪ Bad ∪ T, Real.log p := lt_of_lt_of_le hTgt hge
    obtain ⟨P, hcard, hPp, hPnot, hPb⟩ :=
      lemma_4_1 1 Nat.one_pos (1/6 : ℝ) xeps Ceps (by norm_num) hxepos hCepos (by norm_num) hCx
        (fun x hx _ => hi1 x hx) (fun x hx _ hxx => hi2 x hx hxx) hii H hHnn
        (S ∪ Bad ∪ T) hAllp hAx
    obtain ⟨a, ha⟩ := Finset.card_eq_one.1 hcard
    have hmem : a ∈ P := by rw [ha]; exact Finset.mem_singleton_self a
    have hnot := hPnot a hmem
    refine ⟨a, hPp a hmem, ?_, ?_, (hPb a hmem).1, ?_⟩
    · exact fun hS => hnot (Finset.mem_union_left _ (Finset.mem_union_left _ hS))
    · exact fun hBad => hnot (Finset.mem_union_left _ (Finset.mem_union_right _ hBad))
    · have h2 := (hPb a hmem).2
      have hnum : (1 + 6 * (1/6 : ℝ)) = 2 := by norm_num
      rw [hnum] at h2
      have hu1 := hunion (S ∪ Bad) T
      have hu2 := hunion S Bad
      linarith
  obtain ⟨lo, hlop, hloS, hlobad, hloH, hloub⟩ := hget (D.badPrimes E) (D.badPrimes_prime E)
  obtain ⟨lb, hlbp, hlbS, hlbram, hlbH, hlbub⟩ := hget (D.ramPrimes E) (D.ramPrimes_prime E)
  have hlbbad : lb ∉ D.badPrimes E := fun hc => hlbram (D.badPrimes_subset_ramPrimes E hc)
  obtain ⟨hlo1, hlo2⟩ := D.primeTo_badPrimes E lo hlop hlobad
  obtain ⟨hlb1, hlb2⟩ := D.primeTo_badPrimes E lb hlbp hlbbad
  have hlb3 : D.PrimeToRamification E lb := D.primeTo_ramPrimes E lb hlbp hlbram
  have hSL2lo : D.ImageContainsSL2 E lo :=
    h38 E lo hlop hExcE (Or.inl ⟨le_trans hXH hloH, hpot⟩)
  have hSL2lb : D.ImageContainsSL2 E lb :=
    h38 E lb hlbp hExcE (Or.inl ⟨le_trans hXH hlbH, hpot⟩)
  have hsurj : D.ImageSurjective E lb := D.imageSurjective_of_containsSL2 E lb hlbp hlb3 hSL2lb
  have hramle := D.sum_log_ramPrimes_le E
  have hclo := ABC3.Found.GenEll.cor4_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p) 0 H (lo : ℝ)
      (A₁ + A₂) |Bl| C₈ ((D.degOfDefinition E : ℝ) ^ (1 + eps))
      hd1 hP1 hdP hAineq hxbadle hFB hBnn hC₈pos hxTnn hHle (by linarith [hloub])
  have hclb := ABC3.Found.GenEll.cor4_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p)
      (3 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)) H (lb : ℝ)
      (A₁ + A₂) |Bl| C₈ ((D.degOfDefinition E : ℝ) ^ (1 + eps))
      hd1 hP1 hdP hAineq hxbadle hFB hBnn hC₈pos hxTnn hHle (by linarith [hlbub, hramle])
  exact ⟨lo, lb, hlop, hlbp, hloS, hlbS, ⟨hlo1, hlo2, hlb1, hlb2, hlb3⟩,
    ⟨hSL2lo, hsurj⟩, by linarith [hclo], by linarith [hclb]⟩

/-! ## Corollary 4.4 —— compactly bounded 部分集合の完全 Galois 作用(★北極星) -/

/-- **[GenEll] Corollary 4.4**(Full Galois Actions for Compactly Bounded Subsets)——
★**証明本体**（第 1444 で `Skeleton` から移した）。

原文 (GenEll p.23):
> Corollary 4.4. (Full Galois Actions for Compactly Bounded Subsets)

☆statement の逐語の読みと `Corollary 4.3` との差（`900d` → `100d`、`C·d^{1+ϵ}` → `C·d`）は
`Skeleton/GenEll/Section4.lean` の `cor_4_4` にある。 -/
theorem cor_4_4 (D : EllModuliData) (KV : Set D.EllClass) (hKV : D.CompactlyBounded KV) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (S : Finset ℕ), (∀ p ∈ S, p.Prime) →
        D.MinimalField E → D.cls E ∈ KV → D.cls E ∉ Exc →
        ∃ lo lb : ℕ, Nat.Prime lo ∧ Nat.Prime lb ∧ lo ∉ S ∧ lb ∉ S
          ∧ (D.PrimeToMultPrimes E lo ∧ D.PrimeToLocalHeights E lo
              ∧ D.PrimeToMultPrimes E lb ∧ D.PrimeToLocalHeights E lb
              ∧ D.PrimeToRamification E lb)
          ∧ (D.ImageContainsSL2 E lo ∧ D.ImageSurjective E lb)
          ∧ (lo : ℝ) ≤ 23040 * 100 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p) + C * (D.degOfDefinition E : ℝ)
          ∧ (lb : ℝ) ≤ 23040 * 100 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 6 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p) + C * (D.degOfDefinition E : ℝ) := by
  classical
  -- ★原文 p.23:『`Theorem 3.8` の条件 (a) の代わりに**条件 (b)** を使う』
  -- ——条件 (b) は `ϵ` を使わないので `ϵ := 1` でよい
  obtain ⟨C₈, hC₈pos, Exc, hExc, h38⟩ := theorem_3_8 D KV hKV 1 one_pos
  obtain ⟨hdegbd, hhtbd, -, -⟩ := prop_3_4 D (1/5 : ℝ) (by norm_num)
  obtain ⟨A₁, hA₁⟩ := hdegbd
  obtain ⟨A₂, hA₂⟩ := hhtbd
  obtain ⟨Bl, hBl⟩ := D.faltingsHeight_bddBelow
  obtain ⟨xeps, Ceps, hxepos, hCepos, hCx, hi1, hi2, hii⟩ :=
    ABC3.Found.GenEll.exists_cond_i_ii 1
  obtain ⟨T₀, hT₀prime, hT₀gt⟩ := ABC3.Found.GenEll.exists_finset_primes_sum_log_gt xeps
  -- ★★条件 (b) の『`30` と素』のために `2, 3, 5` を除外集合に入れておく
  obtain ⟨T, hTdef⟩ : ∃ T : Finset ℕ, T = T₀ ∪ ({2, 3, 5} : Finset ℕ) := ⟨_, rfl⟩
  have hTprime : ∀ p ∈ T, Nat.Prime p := by
    intro p hp
    rw [hTdef] at hp
    rcases Finset.mem_union.1 hp with hp' | hp'
    · exact hT₀prime p hp'
    · fin_cases hp' <;> norm_num
  have hT₀sub : T₀ ⊆ T := by rw [hTdef]; exact Finset.subset_union_left
  have hTgt : xeps < ∑ p ∈ T, Real.log p :=
    lt_of_lt_of_le hT₀gt (Finset.sum_le_sum_of_subset_of_nonneg hT₀sub
      (fun p _ _ => Real.log_natCast_nonneg p))
  have hxTnn : (0:ℝ) ≤ ∑ p ∈ T, Real.log p :=
    Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
  have hm2 : (2:ℕ) ∈ T := by rw [hTdef]; exact Finset.mem_union_right _ (by decide)
  have hm3 : (3:ℕ) ∈ T := by rw [hTdef]; exact Finset.mem_union_right _ (by decide)
  have hm5 : (5:ℕ) ∈ T := by rw [hTdef]; exact Finset.mem_union_right _ (by decide)
  have hunion : ∀ U V : Finset ℕ,
      (∑ p ∈ U ∪ V, Real.log p) ≤ (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := by
    intro U V
    have hui : (∑ p ∈ U ∪ V, Real.log p) + ∑ p ∈ U ∩ V, Real.log p
        = (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := Finset.sum_union_inter
    have hnn : (0:ℝ) ≤ ∑ p ∈ U ∩ V, Real.log p :=
      Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
    linarith
  refine ⟨23040 * (5 * |A₁ + A₂| + 28 * |Bl|) + 2 * (∑ p ∈ T, Real.log p) + 1,
    by nlinarith [abs_nonneg (A₁ + A₂), abs_nonneg Bl, hxTnn], Exc, hExc, ?_⟩
  intro E S hSprime _hmin hKVE hExcE
  have hd1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) := by exact_mod_cast D.degOfDefinition_pos E
  have hBnn : (0:ℝ) ≤ |Bl| := abs_nonneg Bl
  have hFB : -|Bl| ≤ D.faltingsHeight (D.cls E) := le_trans (neg_abs_le Bl) (hBl (D.cls E))
  have hL42 := lemma_4_2 (D.multCard E) (D.multCard_pos E) (D.multPrime E) (D.multPrime_prime E)
      (D.localHt E) (D.localHt_pos E)
  have hsum := D.sum_localHt_eq E
  have hxbadle : (∑ p ∈ D.badPrimes E, Real.log p)
      ≤ 5/2 * (23040 * ((D.degOfDefinition E : ℝ) * D.degInf (D.cls E))) := by
    have h1 := hL42.1
    have h3 := hL42.2.2
    rw [hsum] at h1 h3
    have hb := D.sum_log_badPrimes_le E
    linarith
  have hAineq : D.degInf (D.cls E)
      ≤ 12 * (1 + 1/5) * D.faltingsHeight (D.cls E) + (A₁ + A₂) := by
    have h1 := hA₁ (D.cls E)
    have h2 := hA₂ (D.cls E)
    simp only at h1 h2
    linarith
  -- ★★★`Lemma 4.1` を `M = 1`･`h = 0` で適用する(原文 p.23)
  have hget : ∀ Bad : Finset ℕ, (∀ p ∈ Bad, Nat.Prime p) →
      ∃ l : ℕ, Nat.Prime l ∧ l ∉ S ∧ l ∉ Bad ∧ l ∉ T ∧
        (l:ℝ) ≤ 2 * ((∑ p ∈ S, Real.log p) + (∑ p ∈ Bad, Real.log p)
          + (∑ p ∈ T, Real.log p)) := by
    intro Bad hBadp
    have hAllp : ∀ p ∈ S ∪ Bad ∪ T, Nat.Prime p := by
      intro p hp
      rcases Finset.mem_union.1 hp with hp' | hp'
      · rcases Finset.mem_union.1 hp' with hp'' | hp''
        · exact hSprime p hp''
        · exact hBadp p hp''
      · exact hTprime p hp'
    have hge : (∑ p ∈ T, Real.log p) ≤ ∑ p ∈ S ∪ Bad ∪ T, Real.log p :=
      Finset.sum_le_sum_of_subset_of_nonneg Finset.subset_union_right
        (fun p _ _ => Real.log_natCast_nonneg p)
    have hAx : xeps < ∑ p ∈ S ∪ Bad ∪ T, Real.log p := lt_of_lt_of_le hTgt hge
    obtain ⟨P, hcard, hPp, hPnot, hPb⟩ :=
      lemma_4_1 1 Nat.one_pos (1/6 : ℝ) xeps Ceps (by norm_num) hxepos hCepos (by norm_num) hCx
        (fun x hx _ => hi1 x hx) (fun x hx _ hxx => hi2 x hx hxx) hii 0 le_rfl
        (S ∪ Bad ∪ T) hAllp hAx
    obtain ⟨a, ha⟩ := Finset.card_eq_one.1 hcard
    have hmem : a ∈ P := by rw [ha]; exact Finset.mem_singleton_self a
    have hnot := hPnot a hmem
    refine ⟨a, hPp a hmem, ?_, ?_, ?_, ?_⟩
    · exact fun hS => hnot (Finset.mem_union_left _ (Finset.mem_union_left _ hS))
    · exact fun hBad => hnot (Finset.mem_union_left _ (Finset.mem_union_right _ hBad))
    · exact fun hT => hnot (Finset.mem_union_right _ hT)
    · have h2 := (hPb a hmem).2
      have hnum : (1 + 6 * (1/6 : ℝ)) = 2 := by norm_num
      rw [hnum] at h2
      have hu1 := hunion (S ∪ Bad) T
      have hu2 := hunion S Bad
      linarith
  obtain ⟨lo, hlop, hloS, hlobad, hloT, hloub⟩ := hget (D.badPrimes E) (D.badPrimes_prime E)
  obtain ⟨lb, hlbp, hlbS, hlbram, hlbT, hlbub⟩ := hget (D.ramPrimes E) (D.ramPrimes_prime E)
  have hlbbad : lb ∉ D.badPrimes E := fun hc => hlbram (D.badPrimes_subset_ramPrimes E hc)
  obtain ⟨hlo1, hlo2⟩ := D.primeTo_badPrimes E lo hlop hlobad
  obtain ⟨hlb1, hlb2⟩ := D.primeTo_badPrimes E lb hlbp hlbbad
  have hlb3 : D.PrimeToRamification E lb := D.primeTo_ramPrimes E lb hlbp hlbram
  have hcoplo : Nat.Coprime lo 30 :=
    ABC3.Found.GenEll.coprime_thirty_of_prime hlop
      (fun h => hloT (h ▸ hm2)) (fun h => hloT (h ▸ hm3)) (fun h => hloT (h ▸ hm5))
  have hcoplb : Nat.Coprime lb 30 :=
    ABC3.Found.GenEll.coprime_thirty_of_prime hlbp
      (fun h => hlbT (h ▸ hm2)) (fun h => hlbT (h ▸ hm3)) (fun h => hlbT (h ▸ hm5))
  have hSL2lo : D.ImageContainsSL2 E lo := h38 E lo hlop hExcE (Or.inr ⟨hKVE, hlo2, hcoplo⟩)
  have hSL2lb : D.ImageContainsSL2 E lb := h38 E lb hlbp hExcE (Or.inr ⟨hKVE, hlb2, hcoplb⟩)
  have hsurj : D.ImageSurjective E lb := D.imageSurjective_of_containsSL2 E lb hlbp hlb3 hSL2lb
  have hramle := D.sum_log_ramPrimes_le E
  have hclo := ABC3.Found.GenEll.cor44_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p) 0 (lo : ℝ)
      (A₁ + A₂) |Bl| (D.degOfDefinition E : ℝ)
      hd1 hd1 le_rfl hAineq hxbadle hFB hBnn hxTnn (by linarith [hloub])
  have hclb := ABC3.Found.GenEll.cor44_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p)
      (3 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)) (lb : ℝ)
      (A₁ + A₂) |Bl| (D.degOfDefinition E : ℝ)
      hd1 hd1 le_rfl hAineq hxbadle hFB hBnn hxTnn (by linarith [hlbub, hramle])
  exact ⟨lo, lb, hlop, hlbp, hloS, hlbS, ⟨hlo1, hlo2, hlb1, hlb2, hlb3⟩,
    ⟨hSL2lo, hsurj⟩, by linarith [hclo], by linarith [hclb]⟩

/-! ## ★出典の紐付け(`.src`) -/

def cor_4_3.src : Source :=
  { paper := "GenEll", pdfPage := 22, item := "Corollary 4.3",
    sectionId := "genell-cor-4-3" }

def cor_4_3.needs : List ProofObligation :=
  [ .citation "[ABC3]" "theorem_3_8(Skeleton/GenEll/GaloisImage.lean、sorry 0)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.theorem_3_8") 1,
    .citation "[ABC3]" "prop_3_4(Skeleton/GenEll/Section3.lean)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.prop_3_4") 1,
    .citation "[ABC3]" "cor43_numeric(数値の鎖、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.cor43_numeric") 1 ]

def cor_4_4.src : Source :=
  { paper := "GenEll", pdfPage := 23, item := "Corollary 4.4",
    sectionId := "genell-cor-4-4" }

def cor_4_4.needs : List ProofObligation :=
  [ .citation "[ABC3]" "theorem_3_8(条件 (b) の側を使う。sorry 0)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.theorem_3_8") 1,
    .citation "[ABC3]" "cor44_numeric(数値の鎖、h = 0 で 8h の項が消える)"
      (.inProject "ABC3" "ABC3.Found.GenEll.cor44_numeric") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1444）——場所を規約どおりに直した**。" ++
       "☆証明は 2026-08-26（第 367）に `sorry` 0 になっていたが本体が `Skeleton/` にあった。" ++
       "`Lemma 4.1` / `Lemma 4.2` と同じ形（`Skeleton` は写して `Found` へ委譲）に揃えた。" ++
       "★中身は 1 文字も変えていない。" ++
       "★★ただし `EllModuliData` の欄は `Interface` で受けたままであり、" ++
       "「絶対的に証明された」ことを意味しない——欄を埋めるのは `EllModuliWitness.lean` の仕事で、" ++
       "そこには `j(E′)` の整性と Tate 加群の Galois 像が残っている。") 1 ]

end ABC3.Found.GenEll
