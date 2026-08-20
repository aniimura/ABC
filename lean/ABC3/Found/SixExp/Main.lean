/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.ParamsGap

/-!
# six exponentials の本体(背理法の中身)

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★組み立て

これまでに用意した部品をすべてつなぐ。

| 段 | 部品 |
|---|---|
| 定数の設定 | `X = max‖xⱼ‖`、`Y = max‖yᵢ‖`、`A₀`(house の上界)、`Bβ = max 1 ‖β‖` |
| パラメータ | `ParamsGap.lean` の `exists_params_gap` |
| 係数 | `SiegelPoly.lean` の `exists_siegel_coeffsG_poly` |
| 底ベクトル | `DenomClear.lean` の `denomBase`(分母を払う) |
| 零点の鎖 | `LatticeBox.lean` の `latBox`(`|T k| = (N+k)³`) |
| 解析側 | `LatticeBox.lean` の `norm_auxFun_le_of_box` |
| 外挿 | `Assembly.lean` の `extrapolation_induction_scaled` |
| 締め | `Extrapolation.lean` の `coeff_eq_zero_of_auxFun_latticePt` |

★**数え上げがちょうど合う**のが要点である:

  `Δ k · M k · (H k)^e = P^{1+e} · Z₀^{Lb(n+1)}`   (`gap_algebra` / `pow_e_split`)

で、右辺が `exists_params_gap` の左辺そのものになる。

## ★`Gt` の切り詰めについて

`house_auxMatrixG_le` は**すべての** `(p,q)` で底ベクトルの house の上界を要求するが、
`denomBase` の評価は `p+q ≤ Lb` のときしか効かない。そこで

  `Gt pq := denomBase b Eh Lb (min pq.1 L', min pq.2 L')`

と切り詰める。`s` の上では `min` が外れて元に戻るので、`latValG` の値は変わらない。
-/

namespace ABC3.Found.SixExp

open Metric Complex NumberField Finset

variable {K : Type*} [Field K] [NumberField K]

/-- ★★数え上げの代数 —— 段ごとの 3 つの因子がちょうど `P^{1+e}·Z₀^m` にまとまる。 -/
theorem gap_algebra (P u v w : ℝ) (m e : ℕ) :
    u ^ m * (P * v ^ m) * (P ^ e * w ^ m) = P ^ (1 + e) * (u * v * w) ^ m := by
  rw [mul_pow, mul_pow, pow_add, pow_one]
  ring

/-- ★`(P · (A³)^m)^e = P^e · (A^{3e})^m`。 -/
theorem pow_e_split (P A : ℝ) (m e : ℕ) :
    (P * (A ^ 3) ^ m) ^ e = P ^ e * (A ^ (3 * e)) ^ m := by
  rw [mul_pow, ← pow_mul, ← pow_mul, ← pow_mul]
  congr 2
  ring

/-- ★★★★★★**six exponentials の本体** —— `exp(x_j y_i)` が 6 個とも
(分母を払って)代数的整数なら矛盾する。 -/
theorem sixExp_contradiction [DecidableEq (K →+* ℂ)]
    (σ : K →+* ℂ) (x : Fin 2 → ℂ) (y : Fin 3 → ℂ)
    (hx : LinearIndependent ℚ x) (hy : LinearIndependent ℚ y)
    (b : 𝓞 K) (hb : b ≠ 0) (Eh : Fin 2 → Fin 3 → 𝓞 K)
    (hE : ∀ j i, σ ((Eh j i : 𝓞 K) : K) = σ ((b : 𝓞 K) : K) * sixExpVals x y j i) :
    False := by
  classical
  -- ★0. 定数
  set β : ℂ := σ ((b : 𝓞 K) : K) with hβdef
  have hbK : ((b : 𝓞 K) : K) ≠ 0 := by
    simpa using (RingOfIntegers.coe_ne_zero_iff (K := K)).mpr hb
  have hβ0 : β ≠ 0 := by
    rw [hβdef]; exact fun h => hbK (σ.injective (by simpa using h))
  set X : ℝ := max ‖x 0‖ ‖x 1‖ with hXdef
  set Y : ℝ := max ‖y 0‖ (max ‖y 1‖ ‖y 2‖) with hYdef
  have hXb : ∀ j, ‖x j‖ ≤ X := by
    intro j; fin_cases j
    · exact le_max_left _ _
    · exact le_max_right _ _
  have hYb : ∀ i, ‖y i‖ ≤ Y := by
    intro i; fin_cases i
    · exact le_max_left _ _
    · exact le_trans (le_max_left _ _) (le_max_right _ _)
    · exact le_trans (le_max_right _ _) (le_max_right _ _)
  have hX0 : 0 ≤ X := le_trans (norm_nonneg _) (hXb 0)
  have hY0 : 0 < Y := lt_of_lt_of_le (norm_pos_iff.mpr (hy.ne_zero 0)) (hYb 0)
  set A₀ : ℝ := max 1 (max (house ((b : 𝓞 K) : K))
    ((Finset.univ : Finset (Fin 2 × Fin 3)).sup' Finset.univ_nonempty
      (fun p => house ((Eh p.1 p.2 : 𝓞 K) : K)))) with hA₀def
  have hA₀1 : (1:ℝ) ≤ A₀ := le_max_left _ _
  have hA₀b : house ((b : 𝓞 K) : K) ≤ A₀ := le_trans (le_max_left _ _) (le_max_right _ _)
  have hA₀E : ∀ j i, house ((Eh j i : 𝓞 K) : K) ≤ A₀ := by
    intro j i
    refine le_trans ?_ (le_trans (le_max_right _ _) (le_max_right _ _))
    exact Finset.le_sup' (fun p : Fin 2 × Fin 3 => house ((Eh p.1 p.2 : 𝓞 K) : K))
      (Finset.mem_univ (j, i))
  set Bβ : ℝ := max 1 ‖β‖ with hBβdef
  have hBβ1 : (1:ℝ) ≤ Bβ := le_max_left _ _
  have hBβn : ‖β‖ ≤ Bβ := le_max_right _ _
  set cS : ℝ := siegelCpos K with hcSdef
  have hcS1 : (1:ℝ) ≤ cS := one_le_siegelCpos
  set e : ℕ := Module.finrank ℚ K - 1 with hedef
  set Z₀ : ℝ := Bβ ^ 3 * Real.exp (30 * X * Y) * A₀ ^ (3 * e) with hZ₀def
  have hZ₀1 : (1:ℝ) ≤ Z₀ := by
    rw [hZ₀def]
    have h1 : (1:ℝ) ≤ Bβ ^ 3 := one_le_pow₀ hBβ1
    have h2 : (1:ℝ) ≤ Real.exp (30 * X * Y) := Real.one_le_exp (by positivity)
    have h3 : (1:ℝ) ≤ A₀ ^ (3 * e) := one_le_pow₀ hA₀1
    calc (1:ℝ) = 1 * 1 * 1 := by ring
      _ ≤ Bβ ^ 3 * Real.exp (30 * X * Y) * A₀ ^ (3 * e) := by gcongr
  obtain ⟨N, L', hN, hcard, hgapAll⟩ := exists_params_gap A₀ Z₀ cS hA₀1 hZ₀1 hcS1 e
  -- ★1. `s`・`box`・`Gt`・Siegel
  set Lb : ℕ := 2 * L' with hLbdef
  set s : Finset (ℕ × ℕ) := (Finset.range L') ×ˢ (Finset.range L') with hsdef
  have hscard : s.card = L' ^ 2 := by
    rw [hsdef, Finset.card_product, Finset.card_range]; ring
  have hsL : ∀ pq ∈ s, pq.1 + pq.2 ≤ Lb := by
    intro pq hpq
    rw [hsdef, Finset.mem_product, Finset.mem_range, Finset.mem_range] at hpq
    omega
  have hsLe : ∀ pq ∈ s, pq.1 ≤ Lb ∧ pq.2 ≤ Lb := by
    intro pq hpq; have := hsL pq hpq; omega
  have hL'1 : 1 ≤ L' := by
    rcases Nat.eq_zero_or_pos L' with h | h
    · rw [h] at hcard; norm_num at hcard
    · exact h
  have hL1R : (1:ℝ) ≤ (L':ℝ) := by exact_mod_cast hL'1
  set box : Finset (Fin 3 → ℕ) := natBox N with hboxdef
  have hboxcard : box.card = N ^ 3 := card_natBox N
  have hzero : (fun _ => 0) ∈ box := zero_mem_natBox hN
  have hboxN : ∀ m ∈ box, ∑ i : Fin 3, m i ≤ 3 * N := fun m hm => sum_le_of_mem_natBox hm
  have hcard2 : 2 * box.card ≤ s.card := by rw [hboxcard, hscard]; omega
  set Gt : ℕ × ℕ → Fin 3 → 𝓞 K :=
    fun pq => denomBase b Eh Lb (min pq.1 L', min pq.2 L') with hGtdef
  have hGtB : ∀ pq i, house ((Gt pq i : 𝓞 K) : K) ≤ A₀ ^ Lb := by
    intro pq i
    refine house_denomBase_le b Eh Lb hA₀1 hA₀b hA₀E ?_ i
    show min pq.1 L' + min pq.2 L' ≤ 2 * L'
    omega
  have hGts : ∀ pq ∈ s, Gt pq = denomBase b Eh Lb pq := by
    intro pq hpq
    rw [hsdef, Finset.mem_product, Finset.mem_range, Finset.mem_range] at hpq
    have h1 : min pq.1 L' = pq.1 := min_eq_left (by omega)
    have h2 : min pq.2 L' = pq.2 := min_eq_left (by omega)
    show denomBase b Eh Lb (min pq.1 L', min pq.2 L') = denomBase b Eh Lb pq
    rw [h1, h2]
  obtain ⟨c, hcne, hczero, hcbound⟩ :=
    exists_siegel_coeffsG_poly Gt s box hzero hcard2 (one_le_pow₀ hA₀1) hGtB hboxN
  set Cb : ℝ := cS ^ 2 * ((L':ℝ) ^ 2) * A₀ ^ (6 * L' * N) with hCbdef
  have hcb : ∀ pq, house ((c pq : 𝓞 K) : K) ≤ Cb := by
    intro pq
    refine le_trans (hcbound pq) ?_
    rw [hCbdef, hcSdef, hscard]
    have h1 : (((L' ^ 2 : ℕ)) : ℝ) = ((L':ℝ) ^ 2) := by push_cast; ring
    rw [h1, ← pow_mul]
    have h2 : Lb * (3 * N) = 6 * L' * N := by rw [hLbdef]; ring
    rw [h2]
  have hCb1 : (1:ℝ) ≤ Cb := by
    rw [hCbdef]
    have h1 : (1:ℝ) ≤ cS ^ 2 := one_le_pow₀ hcS1
    have h2 : (1:ℝ) ≤ ((L':ℝ) ^ 2) := one_le_pow₀ hL1R
    have h3 : (1:ℝ) ≤ A₀ ^ (6 * L' * N) := one_le_pow₀ hA₀1
    calc (1:ℝ) = 1 * 1 * 1 := by ring
      _ ≤ cS ^ 2 * ((L':ℝ) ^ 2) * A₀ ^ (6 * L' * N) := by gcongr
  -- ★2. 補助関数と鎖
  set c' : ℕ × ℕ → ℂ := fun pq => σ ((c pq : 𝓞 K) : K) with hc'def
  set F : ℂ → ℂ := auxFun x s c' with hFdef
  set P : ℝ := ((L':ℝ) ^ 2) * Cb with hPdef
  have hP1 : (1:ℝ) ≤ P := by
    rw [hPdef]
    have h2 : (1:ℝ) ≤ ((L':ℝ) ^ 2) := one_le_pow₀ hL1R
    calc (1:ℝ) = 1 * 1 := by ring
      _ ≤ ((L':ℝ) ^ 2) * Cb := by gcongr
  set T : ℕ → Finset ℂ := fun k => latBox y (N + k) with hTdef
  set rr : ℕ → ℝ := fun k => 3 * ((N + k + 1 : ℕ) : ℝ) * Y with hrrdef
  set RR : ℕ → ℝ := fun k => 5 * rr k with hRRdef
  set mm : ℕ → ℕ := fun k => Lb * ((N + k) + 1) with hmmdef
  set MM : ℕ → ℝ := fun k => P * (Real.exp (30 * X * Y)) ^ (mm k) with hMMdef
  set HH : ℕ → ℝ := fun k => P * (A₀ ^ 3) ^ (mm k) with hHHdef
  set DD : ℕ → ℝ := fun k => (Bβ ^ 3) ^ (mm k) with hDDdef
  have hrpos : ∀ k, 0 < rr k := by
    intro k
    rw [hrrdef]
    have h : (0:ℝ) < ((N + k + 1 : ℕ) : ℝ) := by positivity
    positivity
  have hr0 : ∀ k, 0 ≤ rr k := fun k => (hrpos k).le
  have h5 : ∀ k, 5 * rr k ≤ RR k := fun k => le_refl _
  have hrR : ∀ k, rr k < RR k := by
    intro k; rw [hRRdef]; have := hrpos k; linarith
  have hFdiff : ∀ k, DifferentiableOn ℂ F (closedBall 0 (RR k)) := fun k =>
    (auxFun_differentiable x s c').differentiableOn
  have hH1 : ∀ k, 1 ≤ HH k := by
    intro k
    rw [hHHdef]
    have h1 : (1:ℝ) ≤ (A₀ ^ 3) ^ (mm k) := one_le_pow₀ (one_le_pow₀ hA₀1)
    calc (1:ℝ) = 1 * 1 := by ring
      _ ≤ P * (A₀ ^ 3) ^ (mm k) := by gcongr
  have hD0 : ∀ k, 0 ≤ DD k := by intro k; rw [hDDdef]; positivity
  have hTcard : ∀ k, (T k).card = (N + k) ^ 3 := fun k => card_latBox hy _
  have hmono : ∀ k, T k ⊆ T (k + 1) := fun k => latBox_mono y (by omega)
  have hTr : ∀ k, ∀ w ∈ T (k + 1), ‖w‖ ≤ rr k := by
    intro k w hw
    have hw' : w ∈ latBox y (N + k + 1) := by
      rw [hTdef] at hw; simpa [Nat.add_assoc] using hw
    have h := norm_le_of_mem_latBox hYb hw'
    rw [hrrdef]; exact h
  -- ★3. σ 像
  have hlatG : ∀ m : Fin 3 → ℕ, latValG Gt m s c = latValG (denomBase b Eh Lb) m s c := by
    intro m
    rw [latValG, latValG]
    refine Finset.sum_congr rfl (fun pq hpq => ?_)
    simp only [auxMatrixG]
    rw [hGts pq hpq]
  have hσlat : ∀ m : Fin 3 → ℕ,
      σ ((latValG Gt m s c : 𝓞 K) : K) = β ^ (Lb * ∑ i : Fin 3, m i) * F (latticePt y m) := by
    intro m
    rw [hlatG m]
    exact map_latValG_denomBase σ b Eh x y Lb hE s hsL c m
  -- ★4. 解析側の上界
  have hsum : (∑ pq ∈ s, ‖c' pq‖) ≤ P := by
    calc (∑ pq ∈ s, ‖c' pq‖) ≤ ∑ _pq ∈ s, Cb :=
          Finset.sum_le_sum (fun pq _ => le_trans (norm_embedding_le_house _ σ) (hcb pq))
      _ = (s.card : ℝ) * Cb := by rw [Finset.sum_const, nsmul_eq_mul]
      _ = P := by rw [hscard, hPdef]; push_cast; ring
  have hMk : ∀ k, ∀ ζ : ℂ, ‖ζ‖ = RR k → ‖F ζ‖ ≤ MM k := by
    intro k ζ hζ
    have h1 := norm_auxFun_le_of_box hXb s hsLe c' hζ
    refine le_trans h1 ?_
    have hexp : Real.exp (2 * (Lb:ℝ) * X * RR k) = (Real.exp (30 * X * Y)) ^ (mm k) := by
      rw [← Real.exp_nat_mul]
      congr 1
      rw [hRRdef, hrrdef, hmmdef]
      push_cast
      ring
    rw [hexp, hMMdef]
    exact mul_le_mul_of_nonneg_right hsum (by positivity)
  -- ★5. 初段
  have h0 : ∀ w ∈ T 0, F w = 0 := by
    intro w hw
    rw [hTdef] at hw
    simp only [Nat.add_zero, latBox, Finset.mem_image] at hw
    obtain ⟨m, hm, rfl⟩ := hw
    have hz := hczero m hm
    have hσ := hσlat m
    rw [hz] at hσ
    simp only [map_zero, ZeroMemClass.coe_zero] at hσ
    have hβp : β ^ (Lb * ∑ i : Fin 3, m i) ≠ 0 := pow_ne_zero _ hβ0
    exact (mul_eq_zero.mp hσ.symm).resolve_left hβp
  -- ★6. 各段の代数的整数
  have halg : ∀ k, ∀ w ∈ T (k + 1), ∃ (α : 𝓞 K) (δ : ℂ), δ ≠ 0 ∧ ‖δ‖ ≤ DD k ∧
      σ ((α : 𝓞 K) : K) = δ * F w ∧ house ((α : 𝓞 K) : K) ≤ HH k := by
    intro k w hw
    have hw' : w ∈ latBox y (N + (k + 1)) := by rw [hTdef] at hw; exact hw
    rw [latBox, Finset.mem_image] at hw'
    obtain ⟨m, hm, rfl⟩ := hw'
    have hm3 : ∑ i : Fin 3, m i ≤ 3 * (N + (k + 1)) := sum_le_of_mem_natBox hm
    refine ⟨latValG Gt m s c, β ^ (Lb * ∑ i : Fin 3, m i), pow_ne_zero _ hβ0, ?_, hσlat m, ?_⟩
    · rw [norm_pow, hDDdef]
      calc ‖β‖ ^ (Lb * ∑ i : Fin 3, m i) ≤ Bβ ^ (Lb * ∑ i : Fin 3, m i) :=
            pow_le_pow_left₀ (norm_nonneg _) hBβn _
        _ ≤ Bβ ^ (Lb * (3 * (N + (k + 1)))) :=
            pow_le_pow_right₀ hBβ1 (Nat.mul_le_mul_left _ hm3)
        _ = (Bβ ^ 3) ^ (mm k) := by
            rw [← pow_mul, hmmdef]
            congr 1
            ring
    · have h := house_latValG_le Gt m s c (one_le_pow₀ hA₀1) hGtB (fun pq _ => hcb pq) hm3
      refine le_trans h ?_
      have hpoweq : ((A₀ ^ Lb) ^ (3 * (N + (k + 1)))) = (A₀ ^ 3) ^ (mm k) := by
        rw [← pow_mul, ← pow_mul, hmmdef]
        congr 1
        ring
      rw [hpoweq, hHHdef, hPdef, hscard]
      push_cast
      ring_nf
      rfl
  -- ★7. 数え上げ
  have hgap : ∀ k, DD k * (MM k * (1/2 : ℝ) ^ (T k).card)
      < (HH k ^ (Module.finrank ℚ K - 1))⁻¹ := by
    intro k
    have hHpos : (0:ℝ) < HH k := lt_of_lt_of_le zero_lt_one (hH1 k)
    rw [← hedef, ← mul_assoc, gap_iff hHpos e ((T k).card), hTcard k]
    have hprod : (DD k * MM k) * HH k ^ e = P ^ (1 + e) * Z₀ ^ (mm k) := by
      rw [hDDdef, hMMdef, hHHdef, pow_e_split, hZ₀def]
      exact gap_algebra P (Bβ ^ 3) (Real.exp (30 * X * Y)) (A₀ ^ (3 * e)) (mm k) e
    rw [hprod, hmmdef]
    exact hgapAll (N + k) (Nat.le_add_right N k)
  -- ★8. 外挿
  have hall := extrapolation_induction_scaled σ rr RR MM HH DD hr0 h5 hrR hFdiff hMk T
    hTr hmono h0 hH1 hD0 halg hgap
  -- ★9. 締め
  have hFall : ∀ m : Fin 3 → ℕ, F (latticePt y m) = 0 := by
    intro m
    have hmem : m ∈ natBox (N + ((∑ j : Fin 3, m j) + 1)) := by
      refine mem_natBox.mpr (fun i => ?_)
      have h : m i ≤ ∑ j : Fin 3, m j :=
        Finset.single_le_sum (f := fun j : Fin 3 => m j) (fun _ _ => Nat.zero_le _)
          (Finset.mem_univ i)
      omega
    have hin : latticePt y m ∈ T ((∑ j : Fin 3, m j) + 1) := by
      simp only [hTdef, latBox, Finset.mem_image]
      exact ⟨m, hmem, rfl⟩
    exact hall _ _ hin
  have hzero' := coeff_eq_zero_of_auxFun_latticePt x y s c'
    (latVec_injective hx hy).injOn hFall
  obtain ⟨pq, hpq, hne⟩ := hcne
  have hc0 : c' pq = 0 := hzero' pq hpq
  rw [hc'def] at hc0
  simp only at hc0
  have h1 : ((c pq : 𝓞 K) : K) = 0 := σ.injective (by simpa using hc0)
  exact hne ((RingOfIntegers.coe_eq_zero_iff (K := K)).mp h1)

end ABC3.Found.SixExp
