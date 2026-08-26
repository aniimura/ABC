import ABC3.Found.GaloisRep.CollDenomFree

/-!
# Galois (G6) 第 278 ブロック —— **★★★★★★★★★★原点近傍を含む共線性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点——単元条件が消えた

> `u·v·w = q` なら **`collDefectE u v w q hq = 0`**(`collDefectE_eq_zero`)

★★★★仮定は `q ∈ I` と `u·v·w = q` **だけ**である。第 256 の `tate_collinear` が
要求していた 6 つの単元条件 `hcp` は**もう要らない**。
したがって **`u ≡ 1` の点(原点近傍)も 3 つ組に入れてよい**。

## ★★★★★★★★分子の整除性は分母を増やしても変わらない(2 度目)

第 274 と同じ手である:

| | 第 274(方程式) | 本ブロック(共線性) |
|---|---|---|
| 抜く分母 | `1 − A` | `1 − (6 つの底)` |
| 掛ける量 | `(1−A)⁶` | `M₁M₂M₃`(行ごと) |
| 既存の分子 | `univA_pow_dvd_numerator` | `cU/cV/cW_pow_dvd_numerator` |

★★★`CollUnivE → CollUniv` の像で `collDefectTruncE = M₁M₂M₃·collDefectTrunc` なので、
分母を `collMbase·d` に取り替えれば既存の分子の整除性がそのまま使える。
★★★★**万有な環を作り直したのではなく、同じ分子の主張を別の分母で読んだ**。

## ★★`1 − q^{m+1}·(底)` は落とせない

★抜けるのは `1 − (底)` だけである。`1 − q^{m+1}·(底)` の側は
**どの領域でも自動的に単元**(`q^{m+1}·x ∈ I`)なので、抜く必要が無い。
★★★尾の部分が最初から素直だったことが、ここで効いている。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `map_collDefectTruncE` | ★★★環準同型で移る |
| `CollUnivE`・`collSpecializeE` | ★★★★★★**`1 − (底)` を反転しない万有な環** |
| `collMbase` | ★★6 点の分母の積 |
| `collDefectTruncE_univ_dvd` | ★★★★★★★★万有な環での整除性 |
| `collDefectTruncE_mem` | ★★★★★★★★切り詰めは `I^n` に入る |
| `collDefectE_eq_zero` | ★★★★★★★★★★**原点近傍を含む共線性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★環準同型で移る -/

section Map

variable {R R' : Type} [CommRing R] [CommRing R']

theorem map_tateM (φ : R →+* R') (p r : R) : φ (tateM p r) = tateM (φ p) (φ r) := by
  rw [tateM, tateM]
  simp only [map_mul, map_pow, map_sub, map_one]

theorem map_tateXtruncEE (φ : R →+* R') (n : ℕ) (p r q : R)
    (hqp : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * p))
    (hqr : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * r)) :
    φ (tateXtruncEE n p r q) = tateXtruncEE n (φ p) (φ r) (φ q) := by
  have hsp : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * p)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ p)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqp m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hsr : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * r)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ r)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqr m (Finset.mem_range.1 hm)), map_mul, map_pow]
  rw [tateXtruncEE, tateXtruncEE]
  simp only [map_add, map_sub, map_mul, map_pow, map_one, map_ofNat]
  rw [map_tateM, hsp, hsr, map_partialEval]

theorem map_tateYtruncEE (φ : R →+* R') (n : ℕ) (p r q : R)
    (hqp : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * p))
    (hqr : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * r)) :
    φ (tateYtruncEE n p r q) = tateYtruncEE n (φ p) (φ r) (φ q) := by
  have hsp : φ (partialSum (fun m => tateYterm (q ^ (m + 1) * p)) n)
      = partialSum (fun m => tateYterm (φ q ^ (m + 1) * φ p)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateYterm φ (hqp m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hsrX : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * r)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ r)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqr m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hsrY : φ (partialSum (fun m => tateYterm (q ^ (m + 1) * r)) n)
      = partialSum (fun m => tateYterm (φ q ^ (m + 1) * φ r)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateYterm φ (hqr m (Finset.mem_range.1 hm)), map_mul, map_pow]
  rw [tateYtruncEE, tateYtruncEE]
  simp only [map_add, map_sub, map_mul, map_pow, map_one]
  rw [map_tateM, hsp, hsrX, hsrY, map_partialEval]

theorem map_collDefectTruncE (φ : R →+* R') (n : ℕ) (u v w q : R)
    (hcq : ∀ (m : ℕ), m < n → ∀ i, IsUnit (1 - q ^ (m + 1) * collPts u v w i)) :
    φ (collDefectTruncE n u v w q) = collDefectTruncE n (φ u) (φ v) (φ w) (φ q) := by
  have hu : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * u) := fun m hm => hcq m hm 0
  have hv : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * v) := fun m hm => hcq m hm 1
  have hw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w) := fun m hm => hcq m hm 2
  have hvw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * (v * w)) := fun m hm => hcq m hm 3
  have huw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * (u * w)) := fun m hm => hcq m hm 4
  have huv : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * (u * v)) := fun m hm => hcq m hm 5
  rw [collDefectTruncE, collDefectTruncE]
  simp only [map_add, map_sub, map_mul]
  rw [map_tateXtruncEE φ n u (v * w) q hu hvw, map_tateYtruncEE φ n u (v * w) q hu hvw,
    map_tateXtruncEE φ n v (u * w) q hv huw, map_tateYtruncEE φ n v (u * w) q hv huw,
    map_tateXtruncEE φ n w (u * v) q hw huv, map_tateYtruncEE φ n w (u * v) q hw huv,
    map_tateM, map_tateM, map_tateM]
  simp only [map_mul]

end Map

/-! ## ★★★★★★`1 − (底)` を反転しない万有な環 -/

/-- ★分母の集合から `1 − (底)` を抜いたもの。 -/
noncomputable def collDenomSetE : Set CollBase :=
  Set.range fun p : ℕ × Fin 6 => 1 - cQ ^ (p.1 + 1) * collDenomBases p.2

noncomputable def collDenomsE : Submonoid CollBase := Submonoid.closure collDenomSetE

/-- ★★★★★★**`1 − (底)` を反転しない共線性の万有な環**。 -/
abbrev CollUnivE : Type := Localization collDenomsE

section Specialize

variable {R : Type} [CommRing R] {I : Ideal R}

theorem collEvalE_isUnit_denomSet (u v w : R)
    (hcq : ∀ (m : ℕ) (i : Fin 6), IsUnit (1 - (u * v * w) ^ (m + 1) * collPts u v w i))
    {x : CollBase} (hx : x ∈ collDenomSetE) : IsUnit (collEval u v w x) := by
  obtain ⟨⟨m, i⟩, rfl⟩ := hx
  rw [map_sub, map_one, map_mul, map_pow, collEval_Q, collEval_denomBases]
  exact hcq m i

theorem collEvalE_isUnit_denoms (u v w : R)
    (hcq : ∀ (m : ℕ) (i : Fin 6), IsUnit (1 - (u * v * w) ^ (m + 1) * collPts u v w i))
    (y : collDenomsE) : IsUnit (collEval u v w (y : CollBase)) := by
  obtain ⟨x, hx⟩ := y
  refine Submonoid.closure_induction (s := collDenomSetE)
    (motive := fun z _ => IsUnit (collEval u v w z)) ?_ ?_ ?_ hx
  · intro z hz
    exact collEvalE_isUnit_denomSet u v w hcq hz
  · simp
  · intro z₁ z₂ _ _ h₁ h₂
    rw [map_mul]
    exact h₁.mul h₂

/-- ★★★★★★**`1 − (底)` を反転しない万有な環からの特殊化**。 -/
noncomputable def collSpecializeE (u v w : R)
    (hcq : ∀ (m : ℕ) (i : Fin 6), IsUnit (1 - (u * v * w) ^ (m + 1) * collPts u v w i)) :
    CollUnivE →+* R :=
  IsLocalization.lift (collEvalE_isUnit_denoms u v w hcq)

theorem collSpecializeE_base (u v w : R)
    (hcq : ∀ (m : ℕ) (i : Fin 6), IsUnit (1 - (u * v * w) ^ (m + 1) * collPts u v w i))
    (x : CollBase) :
    collSpecializeE u v w hcq (algebraMap CollBase CollUnivE x) = collEval u v w x :=
  IsLocalization.lift_eq _ x

end Specialize

/-! ## ★万有な環の側の元と単元性 -/

/-- ★万有な環 `CollUnivE` の中の `U`。 -/
noncomputable def fU : CollUnivE := algebraMap CollBase CollUnivE cU
/-- ★万有な環 `CollUnivE` の中の `V`。 -/
noncomputable def fV : CollUnivE := algebraMap CollBase CollUnivE cV
/-- ★万有な環 `CollUnivE` の中の `W`。 -/
noncomputable def fW : CollUnivE := algebraMap CollBase CollUnivE cW

theorem fQ_eq : fU * fV * fW = algebraMap CollBase CollUnivE cQ := by
  rw [cQ, map_mul, map_mul, fU, fV, fW]

theorem mem_collDenomsE_of_mem (x : CollBase) (hx : x ∈ collDenomSetE) : x ∈ collDenomsE :=
  Submonoid.subset_closure hx

theorem collPts_algebraMap {S : Type} [CommRing S] (f : CollBase →+* S) (i : Fin 6) :
    collPts (f cU) (f cV) (f cW) i = f (collDenomBases i) := by
  fin_cases i <;> simp [collPts, collDenomBases, map_mul]

theorem isUnit_one_sub_fq (m : ℕ) (i : Fin 6) :
    IsUnit (1 - (algebraMap CollBase CollUnivE cQ) ^ (m + 1)
      * algebraMap CollBase CollUnivE (collDenomBases i)) := by
  have h : (1 : CollBase) - cQ ^ (m + 1) * collDenomBases i ∈ collDenomsE :=
    mem_collDenomsE_of_mem _ ⟨(m, i), rfl⟩
  have h2 := IsLocalization.map_units (M := collDenomsE) CollUnivE (⟨_, h⟩ : collDenomsE)
  simpa using h2

theorem hcq_f : ∀ (m : ℕ) (i : Fin 6),
    IsUnit (1 - (fU * fV * fW) ^ (m + 1) * collPts fU fV fW i) := by
  intro m i
  rw [fQ_eq]
  simp only [fU, fV, fW]
  rw [collPts_algebraMap (algebraMap CollBase CollUnivE) i]
  exact isUnit_one_sub_fq m i

theorem hcq_k : ∀ (m : ℕ) (i : Fin 6),
    IsUnit (1 - (kU * kV * kW) ^ (m + 1) * collPts kU kV kW i) := by
  intro m i
  rw [kQ_eq]
  simp only [kU, kV, kW]
  rw [collPts_algebraMap (algebraMap CollBase CollUniv) i]
  exact isUnit_one_sub_kq m i

theorem hcp_k : ∀ i, IsUnit (1 - collPts kU kV kW i) := by
  intro i
  simp only [kU, kV, kW]
  rw [collPts_algebraMap (algebraMap CollBase CollUniv) i]
  exact isUnit_one_sub_kbase i

/-- ★`CollUniv` への評価は構造射に一致する。 -/
theorem collEval_kU_kV_kW (x : CollBase) :
    collEval kU kV kW x = algebraMap CollBase CollUniv x := by
  have h : (collEval kU kV kW : CollBase →+* CollUniv)
      = (algebraMap CollBase CollUniv : CollBase →+* CollUniv) := by
    apply MvPolynomial.ringHom_ext
    · intro r
      rw [collEval, MvPolynomial.eval₂Hom_C]
      simp
    · intro i
      fin_cases i
      · rw [show (⟨0, by norm_num⟩ : Fin 3) = 0 from rfl, collEval, MvPolynomial.eval₂Hom_X']
        rw [show (![kU, kV, kW] 0) = kU from rfl, kU, cU]
      · rw [show (⟨1, by norm_num⟩ : Fin 3) = 1 from rfl, collEval, MvPolynomial.eval₂Hom_X']
        rw [show (![kU, kV, kW] 1) = kV from rfl, kV, cV]
      · rw [show (⟨2, by norm_num⟩ : Fin 3) = 2 from rfl, collEval, MvPolynomial.eval₂Hom_X']
        rw [show (![kU, kV, kW] 2) = kW from rfl, kW, cW]
  exact DFunLike.congr_fun h x

theorem collDenomsE_le : collDenomsE ≤ collDenoms := by
  refine Submonoid.closure_le.2 fun x hx => ?_
  obtain ⟨⟨m, i⟩, rfl⟩ := hx
  exact Submonoid.subset_closure (Or.inr ⟨(m, i), rfl⟩)

/-- ★★★`CollUnivE → CollUniv`。 -/
noncomputable def iotaCE : CollUnivE →+* CollUniv := collSpecializeE kU kV kW hcq_k

theorem iotaCE_base (x : CollBase) :
    iotaCE (algebraMap CollBase CollUnivE x) = algebraMap CollBase CollUniv x := by
  rw [iotaCE, collSpecializeE_base, collEval_kU_kV_kW]

theorem iotaCE_fU : iotaCE fU = kU := by rw [fU, iotaCE_base, kU]
theorem iotaCE_fV : iotaCE fV = kV := by rw [fV, iotaCE_base, kV]
theorem iotaCE_fW : iotaCE fW = kW := by rw [fW, iotaCE_base, kW]

theorem dvd_pow_of_baseCE (n : ℕ) (x : CollUnivE)
    (H : ∀ (P : CollBase) (d : collDenomsE),
      x * algebraMap CollBase CollUnivE (d : CollBase) = algebraMap CollBase CollUnivE P →
      cQ ^ n ∣ P) :
    (fU * fV * fW) ^ n ∣ x := by
  obtain ⟨⟨P, d⟩, hPd⟩ := IsLocalization.surj collDenomsE x
  obtain ⟨c, hc⟩ := H P d hPd
  obtain ⟨y, hy⟩ := IsLocalization.map_units CollUnivE d
  refine ⟨algebraMap CollBase CollUnivE c * ↑y⁻¹, ?_⟩
  have hkey : x * (y : CollUnivE)
      = (fU * fV * fW) ^ n * algebraMap CollBase CollUnivE c := by
    rw [hy, hPd, hc, map_mul, map_pow, fQ_eq]
  calc x = x * (y : CollUnivE) * ↑y⁻¹ := by rw [mul_assoc, y.mul_inv, mul_one]
    _ = (fU * fV * fW) ^ n * algebraMap CollBase CollUnivE c * ↑y⁻¹ := by rw [hkey]
    _ = (fU * fV * fW) ^ n * (algebraMap CollBase CollUnivE c * ↑y⁻¹) := by ring

/-! ## ★★★★★★★★万有な環での整除性 -/

/-- ★6 点の分母の積。 -/
noncomputable def collMbase : CollBase :=
  (1 - collDenomBases 0) ^ 3 * (1 - collDenomBases 3) ^ 3
    * ((1 - collDenomBases 1) ^ 3 * (1 - collDenomBases 4) ^ 3)
    * ((1 - collDenomBases 2) ^ 3 * (1 - collDenomBases 5) ^ 3)

theorem collMbase_mem : collMbase ∈ collDenoms := by
  have hb : ∀ i : Fin 6, (1 : CollBase) - collDenomBases i ∈ collDenoms :=
    fun i => Submonoid.subset_closure (Or.inl ⟨i, rfl⟩)
  refine collDenoms.mul_mem (collDenoms.mul_mem
    (collDenoms.mul_mem (collDenoms.pow_mem (hb 0) 3) (collDenoms.pow_mem (hb 3) 3))
    (collDenoms.mul_mem (collDenoms.pow_mem (hb 1) 3) (collDenoms.pow_mem (hb 4) 3)))
    (collDenoms.mul_mem (collDenoms.pow_mem (hb 2) 3) (collDenoms.pow_mem (hb 5) 3))

theorem collMbase_map : algebraMap CollBase CollUniv collMbase
    = tateM kU (kV * kW) * tateM kV (kU * kW) * tateM kW (kU * kV) := by
  have h : ∀ i : Fin 6, algebraMap CollBase CollUniv (collDenomBases i) = collPts kU kV kW i := by
    intro i
    simp only [kU, kV, kW]
    exact (collPts_algebraMap (algebraMap CollBase CollUniv) i).symm
  have e0 : collPts kU kV kW 0 = kU := rfl
  have e1 : collPts kU kV kW 1 = kV := rfl
  have e2 : collPts kU kV kW 2 = kW := rfl
  have e3 : collPts kU kV kW 3 = kV * kW := rfl
  have e4 : collPts kU kV kW 4 = kU * kW := rfl
  have e5 : collPts kU kV kW 5 = kU * kV := rfl
  rw [collMbase, tateM, tateM, tateM]
  simp only [map_mul, map_pow, map_sub, map_one, h, e0, e1, e2, e3, e4, e5]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**`1 − (底)` を反転しない環でも `(UVW)^n` で割れる**。

★★分子の整除性(第 256)は分母を `collMbase·d` に取り替えてもそのまま使える。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collDefectTruncE_univ_dvd (n : ℕ) :
    (fU * fV * fW) ^ n ∣ collDefectTruncE n fU fV fW (fU * fV * fW) := by
  refine dvd_pow_of_baseCE n _ fun P d hPd => ?_
  have h2 : iotaCE (collDefectTruncE n fU fV fW (fU * fV * fW))
      = collDefectTruncE n kU kV kW (kU * kV * kW) := by
    rw [map_collDefectTruncE iotaCE n fU fV fW (fU * fV * fW) (fun m _ i => hcq_f m i),
      iotaCE_fU, iotaCE_fV, iotaCE_fW, map_mul, map_mul, iotaCE_fU, iotaCE_fV, iotaCE_fW]
  have h1 := congrArg iotaCE hPd
  rw [map_mul, iotaCE_base, iotaCE_base, h2,
    collDefectTruncE_eq n kU kV kW (kU * kV * kW) hcp_k, ← collMbase_map] at h1
  have hmem : (collMbase * (d : CollBase)) ∈ collDenoms :=
    collDenoms.mul_mem collMbase_mem (collDenomsE_le d.2)
  have h3 : collDefectTrunc n kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv ((⟨collMbase * (d : CollBase), hmem⟩ : collDenoms)
        : CollBase) = algebraMap CollBase CollUniv P := by
    rw [show (((⟨collMbase * (d : CollBase), hmem⟩ : collDenoms) : CollBase))
      = collMbase * (d : CollBase) from rfl, map_mul]
    rw [← h1]
    ring
  cases n with
  | zero => simp
  | succ m =>
    exact cQ_pow_dvd (m + 1) (cU_pow_dvd_numerator m P _ h3) (cV_pow_dvd_numerator m P _ h3)
      (cW_pow_dvd_numerator m P _ h3)

/-! ## ★★★★★★★★★★原点近傍を含む共線性 -/

section Final

variable {R : Type} [CommRing R] {I : Ideal R}

set_option maxHeartbeats 1600000 in
theorem collDefectTruncE_mem [IsAdicComplete I R] (n : ℕ) (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) : collDefectTruncE n u v w q ∈ I ^ n := by
  have hcq : ∀ (m : ℕ) (i : Fin 6), IsUnit (1 - (u * v * w) ^ (m + 1) * collPts u v w i) := by
    intro m i
    rw [huvw]
    exact isUnit_one_sub (I := I) (pow_succ_mul_mem hq m)
  set φ := collSpecializeE u v w hcq with hφ
  have hU : φ fU = u := by rw [hφ, fU, collSpecializeE_base, collEval_U]
  have hV : φ fV = v := by rw [hφ, fV, collSpecializeE_base, collEval_V]
  have hW : φ fW = w := by rw [hφ, fW, collSpecializeE_base, collEval_W]
  obtain ⟨c, hc⟩ := collDefectTruncE_univ_dvd n
  have hmap : φ (collDefectTruncE n fU fV fW (fU * fV * fW)) = collDefectTruncE n u v w q := by
    rw [map_collDefectTruncE φ n fU fV fW (fU * fV * fW) (fun m _ i => hcq_f m i),
      hU, hV, hW, map_mul, map_mul, hU, hV, hW, huvw]
  rw [← hmap, hc, map_mul, map_pow, map_mul, map_mul, hU, hV, hW, huvw]
  exact Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)

/-- ★★★★★★★★★★**原点近傍を含む共線性**——単元条件 `hcp` が要らない。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collDefectE_eq_zero [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) : collDefectE u v w q hq = 0 := by
  refine eq_zero_of_mem_pow (I := I) fun n => ?_
  have h1 := collDefectE_sub_trunc (I := I) n u v w q hq
  have h2 := collDefectTruncE_mem (I := I) n u v w q hq huvw
  have h3 := Ideal.add_mem (I ^ n) h1 h2
  simpa using h3

end Final

/-! ## ★出典の紐付け(`.src`) -/

def CollUnivE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——1-(底) を反転しない共線性の万有な環)",
    sectionId := "genell-def-3-3" }

def collDefectE_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍を含む共線性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
