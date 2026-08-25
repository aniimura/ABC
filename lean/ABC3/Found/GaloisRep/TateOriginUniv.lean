import ABC3.Found.GaloisRep.TateUnitGroup

/-!
# Galois (G6) 第 274 ブロック —— **★★★★★★★★★`1 − a` を反転しない万有な環**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★原点近傍のために

葉 (e) の残りは **`u ≡ 1 mod 𝔪`** の領域である。そこでは `1 − u` が単元でないので
`X(u) = u/(1−u)² + …` は `R` に入らない——**`K` に極を持つ**。
`1 + 𝔪` は開かつ閉なので、外から近づけることもできない(近似では届かない)。

★★★したがって `K` の水準で方程式を示すほかない。本ブロックはその**土台**を作る。

## ★★★★★★分母を払った切り詰め

    XE_n(a,w,q) := a + (1−a)²·(尾の部分)          (`= (1−a)²·X_n`)
    YE_n(a,w,q) := a² + (1−a)³·(尾の部分)          (`= (1−a)³·Y_n`)
    DE_n := YE² + XE·YE·(1−a) − XE³ − a₄·XE·(1−a)⁴ − a₆·(1−a)⁶   (`= (1−a)⁶·D_n`)

★★**`Ring.inverse (1−a)` を一切含まない**。`(1−a)²·f(a) = a`(第 262)を使って
極を先に払っておくのが要点である。

## ★★★★★★★★`1 − A` を反転しない万有な環

`tateDenomSet` から `1 − A` を抜いた `tateDenomSetE` で局所化した `TateUnivE` を作る。

    (AW)^n ∣ DE_n(A, W, AW)        (`TateUnivE` の中で)

★★★これは既存の整除性(第 224・第 240)から**ただで出る**:
`TateUnivE → TateUniv` の像で `DE_n = (1−A)⁶·D_n` であり、
分母を `(1−A)⁶·d` に取り替えれば `univA_pow_dvd_numerator` がそのまま使える。
★★★★**分子の整除性は分母を増やしても変わらない**——ここが効いた。

## ★★★★★★★★★到達点

> `a·w = q`、`1 − w` が単元なら **`DE_n(a,w,q) ∈ I^n`**——`1 − a` の単元性は要らない。

★これで `u ≡ 1` でも「`(1−u)⁶ ×(方程式の差)` が `I^n` に入る」が言える。
`K` の付値で `v((1−u)⁶) < ∞` なので、`n → ∞` で差は 0 になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXtruncE`・`tateYtruncE`・`tateDefectTruncE` | ★★★★★★分母を払った切り詰め |
| `tateDefectTruncE_eq` | ★★★★`DE_n = (1−a)⁶·D_n` |
| `map_tateDefectTruncE` | ★★★環準同型で移る |
| `TateUnivE`・`tateSpecializeE` | ★★★★★★**`1 − A` を反転しない万有な環** |
| `tateDefectTruncE_univ_dvd` | ★★★★★★★★万有な環での整除性 |
| `tateDefectTruncE_mem` | ★★★★★★★★★**`1 − a` の単元性が要らない `I^n` 所属** |
-/

namespace ABC3.Found.GaloisRep

open MvPolynomial

/-! ## ★★★★★★分母を払った切り詰め -/

section Trunc

variable {R : Type} [CommRing R]

/-- ★★★★★★**`(1−a)²` を掛けた `X` の切り詰め**——`Ring.inverse (1−a)` を含まない。 -/
noncomputable def tateXtruncE (n : ℕ) (a w q : R) : R :=
  a + (1 - a) ^ 2 * (partialSum (fun m => tateXterm (q ^ (m + 1) * a)) n
    + (tateXterm w + partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
    - 2 * partialEval (sigmaSeries 1) q n)

/-- ★★★★★★**`(1−a)³` を掛けた `Y` の切り詰め**。 -/
noncomputable def tateYtruncE (n : ℕ) (a w q : R) : R :=
  a ^ 2 + (1 - a) ^ 3 * (partialSum (fun m => tateYterm (q ^ (m + 1) * a)) n
    - (tateXterm w + partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
    - (tateYterm w + partialSum (fun m => tateYterm (q ^ (m + 1) * w)) n)
    + partialEval (sigmaSeries 1) q n)

theorem tateXtruncE_eq (n : ℕ) (a w q : R) (ha : IsUnit (1 - a)) :
    tateXtruncE n a w q = (1 - a) ^ 2 * tateXtrunc n a w q := by
  rw [tateXtruncE, tateXtrunc]
  linear_combination -mul_tateXterm' ha

theorem tateYtruncE_eq (n : ℕ) (a w q : R) (ha : IsUnit (1 - a)) :
    tateYtruncE n a w q = (1 - a) ^ 3 * tateYtrunc n a w q := by
  rw [tateYtruncE, tateYtrunc]
  linear_combination -mul_tateYterm' ha

/-- ★★★★★★**`(1−a)⁶` を掛けた方程式の差**——`Ring.inverse (1−a)` を含まない。 -/
noncomputable def tateDefectTruncE (n : ℕ) (a w q : R) : R :=
  tateYtruncE n a w q ^ 2 + tateXtruncE n a w q * tateYtruncE n a w q * (1 - a)
    - (tateXtruncE n a w q ^ 3 + partialEval tateA4 q n * tateXtruncE n a w q * (1 - a) ^ 4
      + partialEval tateA6 q n * (1 - a) ^ 6)

/-- ★★★★**`DE_n = (1−a)⁶·D_n`**。 -/
theorem tateDefectTruncE_eq (n : ℕ) (a w q : R) (ha : IsUnit (1 - a)) :
    tateDefectTruncE n a w q = (1 - a) ^ 6 * tateDefectTrunc n a w q := by
  rw [tateDefectTruncE, tateDefectTrunc, tateXtruncE_eq n a w q ha, tateYtruncE_eq n a w q ha]
  ring

end Trunc

/-! ## ★★★環準同型で移る -/

section Map

variable {R R' : Type} [CommRing R] [CommRing R']

theorem map_tateXtruncE (φ : R →+* R') (n : ℕ) (a w q : R) (hw : IsUnit (1 - w))
    (hqa : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * a))
    (hqw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w)) :
    φ (tateXtruncE n a w q) = tateXtruncE n (φ a) (φ w) (φ q) := by
  have hsa : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * a)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ a)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqa m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hsw : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ w)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqw m (Finset.mem_range.1 hm)), map_mul, map_pow]
  rw [tateXtruncE, tateXtruncE]
  simp only [map_add, map_sub, map_mul, map_pow, map_one, map_ofNat]
  rw [map_tateXterm φ hw, hsa, hsw, map_partialEval]

theorem map_tateYtruncE (φ : R →+* R') (n : ℕ) (a w q : R) (hw : IsUnit (1 - w))
    (hqa : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * a))
    (hqw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w)) :
    φ (tateYtruncE n a w q) = tateYtruncE n (φ a) (φ w) (φ q) := by
  have hsa : φ (partialSum (fun m => tateYterm (q ^ (m + 1) * a)) n)
      = partialSum (fun m => tateYterm (φ q ^ (m + 1) * φ a)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateYterm φ (hqa m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hswX : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ w)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqw m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hswY : φ (partialSum (fun m => tateYterm (q ^ (m + 1) * w)) n)
      = partialSum (fun m => tateYterm (φ q ^ (m + 1) * φ w)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateYterm φ (hqw m (Finset.mem_range.1 hm)), map_mul, map_pow]
  rw [tateYtruncE, tateYtruncE]
  simp only [map_add, map_sub, map_mul, map_pow, map_one]
  rw [map_tateXterm φ hw, map_tateYterm φ hw, hsa, hswX, hswY, map_partialEval]

theorem map_tateDefectTruncE (φ : R →+* R') (n : ℕ) (a w q : R) (hw : IsUnit (1 - w))
    (hqa : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * a))
    (hqw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w)) :
    φ (tateDefectTruncE n a w q) = tateDefectTruncE n (φ a) (φ w) (φ q) := by
  rw [tateDefectTruncE, tateDefectTruncE]
  simp only [map_add, map_sub, map_mul, map_pow, map_one]
  rw [map_tateXtruncE φ n a w q hw hqa hqw, map_tateYtruncE φ n a w q hw hqa hqw,
    map_partialEval, map_partialEval]

end Map

/-! ## ★★★★★★`1 − A` を反転しない万有な環 -/

/-- ★分母の集合から `1 − A` を抜いたもの。 -/
noncomputable def tateDenomSetE : Set TateBase :=
  {1 - univW} ∪ (Set.range fun m : ℕ => 1 - univQ ^ (m + 1) * univA)
    ∪ (Set.range fun m : ℕ => 1 - univQ ^ (m + 1) * univW)

noncomputable def tateDenomsE : Submonoid TateBase := Submonoid.closure tateDenomSetE

/-- ★★★★★★**`1 − A` を反転しない万有な環**。 -/
abbrev TateUnivE : Type := Localization tateDenomsE

section Specialize

variable {R : Type} [CommRing R] {I : Ideal R}

theorem tateEvalE_isUnit_denomSet (a w : R) (hw : IsUnit (1 - w))
    (hqa : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * a))
    (hqw : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * w))
    {x : TateBase} (hx : x ∈ tateDenomSetE) : IsUnit (tateEval a w x) := by
  rcases hx with (hx | ⟨m, rfl⟩) | ⟨m, rfl⟩
  · rw [Set.mem_singleton_iff] at hx
    subst hx
    simpa using hw
  · have hval : tateEval a w (1 - univQ ^ (m + 1) * univA) = 1 - (a * w) ^ (m + 1) * a := by
      simp [univQ]
    rw [hval]
    exact hqa m
  · have hval : tateEval a w (1 - univQ ^ (m + 1) * univW) = 1 - (a * w) ^ (m + 1) * w := by
      simp [univQ]
    rw [hval]
    exact hqw m

theorem tateEvalE_isUnit_denoms (a w : R) (hw : IsUnit (1 - w))
    (hqa : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * a))
    (hqw : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * w))
    (y : tateDenomsE) : IsUnit (tateEval a w (y : TateBase)) := by
  obtain ⟨x, hx⟩ := y
  refine Submonoid.closure_induction (s := tateDenomSetE)
    (motive := fun z _ => IsUnit (tateEval a w z)) ?_ ?_ ?_ hx
  · intro z hz
    exact tateEvalE_isUnit_denomSet a w hw hqa hqw hz
  · simp
  · intro z₁ z₂ _ _ h₁ h₂
    rw [map_mul]
    exact h₁.mul h₂

/-- ★★★★★★**`1 − A` を反転しない万有な環からの特殊化**——`1 − a` の単元性が要らない。 -/
noncomputable def tateSpecializeE (a w : R) (hw : IsUnit (1 - w))
    (hqa : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * a))
    (hqw : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * w)) : TateUnivE →+* R :=
  IsLocalization.lift (tateEvalE_isUnit_denoms a w hw hqa hqw)

theorem tateSpecializeE_base (a w : R) (hw : IsUnit (1 - w))
    (hqa : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * a))
    (hqw : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * w)) (x : TateBase) :
    tateSpecializeE a w hw hqa hqw (algebraMap TateBase TateUnivE x) = tateEval a w x :=
  IsLocalization.lift_eq _ x

end Specialize

/-! ## ★万有な環の側の元と単元性 -/

/-- ★万有な環 `TateUnivE` の中の `A`。 -/
noncomputable def eA : TateUnivE := algebraMap TateBase TateUnivE univA
/-- ★万有な環 `TateUnivE` の中の `W`。 -/
noncomputable def eW : TateUnivE := algebraMap TateBase TateUnivE univW

theorem eA_mul_eW : eA * eW = algebraMap TateBase TateUnivE univQ := by
  rw [eA, eW, ← map_mul, univQ]

theorem mem_tateDenomsE_of_mem (x : TateBase) (hx : x ∈ tateDenomSetE) : x ∈ tateDenomsE :=
  Submonoid.subset_closure hx

theorem isUnit_one_sub_eW : IsUnit (1 - eW) := by
  have h : (1 : TateUnivE) - eW = algebraMap TateBase TateUnivE (1 - univW) := by
    rw [map_sub, map_one, eW]
  rw [h]
  exact IsLocalization.map_units TateUnivE
    ⟨1 - univW, mem_tateDenomsE_of_mem _ (by left; left; rfl)⟩

theorem isUnit_one_sub_eQA (m : ℕ) : IsUnit (1 - (eA * eW) ^ (m + 1) * eA) := by
  have h : (1 : TateUnivE) - (eA * eW) ^ (m + 1) * eA
      = algebraMap TateBase TateUnivE (1 - univQ ^ (m + 1) * univA) := by
    rw [map_sub, map_one, map_mul, map_pow, ← eA_mul_eW, eA]
  rw [h]
  exact IsLocalization.map_units TateUnivE
    ⟨_, mem_tateDenomsE_of_mem _ (by left; right; exact ⟨m, rfl⟩)⟩

theorem isUnit_one_sub_eQW (m : ℕ) : IsUnit (1 - (eA * eW) ^ (m + 1) * eW) := by
  have h : (1 : TateUnivE) - (eA * eW) ^ (m + 1) * eW
      = algebraMap TateBase TateUnivE (1 - univQ ^ (m + 1) * univW) := by
    rw [map_sub, map_one, map_mul, map_pow, ← eA_mul_eW, eW]
  rw [h]
  exact IsLocalization.map_units TateUnivE
    ⟨_, mem_tateDenomsE_of_mem _ (by right; exact ⟨m, rfl⟩)⟩

/-! ## ★★★★局所化から分子へ -/

theorem dvd_pow_of_baseE (n : ℕ) (x : TateUnivE)
    (H : ∀ (P : TateBase) (d : tateDenomsE),
      x * algebraMap TateBase TateUnivE (d : TateBase) = algebraMap TateBase TateUnivE P →
      univQ ^ n ∣ P) :
    (eA * eW) ^ n ∣ x := by
  obtain ⟨⟨P, d⟩, hPd⟩ := IsLocalization.surj tateDenomsE x
  obtain ⟨c, hc⟩ := H P d hPd
  obtain ⟨v, hv⟩ := IsLocalization.map_units TateUnivE d
  refine ⟨algebraMap TateBase TateUnivE c * ↑v⁻¹, ?_⟩
  have hkey : x * (v : TateUnivE) = (eA * eW) ^ n * algebraMap TateBase TateUnivE c := by
    rw [hv, hPd, hc, map_mul, map_pow, eA_mul_eW]
  calc x = x * (v : TateUnivE) * ↑v⁻¹ := by rw [mul_assoc, v.mul_inv, mul_one]
    _ = (eA * eW) ^ n * algebraMap TateBase TateUnivE c * ↑v⁻¹ := by rw [hkey]
    _ = (eA * eW) ^ n * (algebraMap TateBase TateUnivE c * ↑v⁻¹) := by ring

/-! ## ★★★★★★★★万有な環での整除性 -/

/-- ★`TateUniv` への評価は構造射に一致する。 -/
theorem tateEval_uA_uW (x : TateBase) :
    tateEval uA uW x = algebraMap TateBase TateUniv x := by
  have h : (tateEval uA uW : TateBase →+* TateUniv)
      = (algebraMap TateBase TateUniv : TateBase →+* TateUniv) := by
    apply MvPolynomial.ringHom_ext
    · intro r
      rw [tateEval, MvPolynomial.eval₂Hom_C]
      simp
    · intro i
      fin_cases i
      · rw [show (⟨0, by norm_num⟩ : Fin 2) = 0 from rfl, tateEval, MvPolynomial.eval₂Hom_X']
        rw [show (![uA, uW] 0) = uA from rfl, uA, univA]
      · rw [show (⟨1, by norm_num⟩ : Fin 2) = 1 from rfl, tateEval, MvPolynomial.eval₂Hom_X']
        rw [show (![uA, uW] 1) = uW from rfl, uW, univW]
  exact DFunLike.congr_fun h x

theorem tateDenomsE_le : tateDenomsE ≤ tateDenoms := by
  refine Submonoid.closure_le.2 fun x hx => ?_
  refine Submonoid.subset_closure ?_
  rcases hx with (hx | ⟨m, rfl⟩) | ⟨m, rfl⟩
  · rw [Set.mem_singleton_iff] at hx
    subst hx
    exact Or.inl (Or.inl (Or.inr rfl))
  · exact Or.inl (Or.inr ⟨m, rfl⟩)
  · exact Or.inr ⟨m, rfl⟩

/-- ★★★`TateUnivE → TateUniv`。 -/
noncomputable def iotaE : TateUnivE →+* TateUniv :=
  tateSpecializeE uA uW isUnit_one_sub_uW isUnit_one_sub_uQA isUnit_one_sub_uQW

theorem iotaE_base (x : TateBase) :
    iotaE (algebraMap TateBase TateUnivE x) = algebraMap TateBase TateUniv x := by
  rw [iotaE, tateSpecializeE_base, tateEval_uA_uW]

theorem iotaE_eA : iotaE eA = uA := by rw [eA, iotaE_base, uA]

theorem iotaE_eW : iotaE eW = uW := by rw [eW, iotaE_base, uW]

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★**`1 − A` を反転しない環でも `(AW)^n` で割れる**。

★★分子の整除性(第 240)は分母を `(1−A)⁶·d` に取り替えてもそのまま使える。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateDefectTruncE_univ_dvd (n : ℕ) :
    (eA * eW) ^ n ∣ tateDefectTruncE n eA eW (eA * eW) := by
  refine dvd_pow_of_baseE n _ fun P d hPd => ?_
  have h2 : iotaE (tateDefectTruncE n eA eW (eA * eW))
      = tateDefectTruncE n uA uW (uA * uW) := by
    rw [map_tateDefectTruncE iotaE n eA eW (eA * eW) isUnit_one_sub_eW
      (fun m _ => isUnit_one_sub_eQA m) (fun m _ => isUnit_one_sub_eQW m),
      iotaE_eA, iotaE_eW, map_mul, iotaE_eA, iotaE_eW]
  have h1 := congrArg iotaE hPd
  rw [map_mul, iotaE_base, iotaE_base, h2,
    tateDefectTruncE_eq n uA uW (uA * uW) isUnit_one_sub_uA] at h1
  have hd6 : ((1 : TateUniv) - uA) ^ 6 = algebraMap TateBase TateUniv ((1 - univA) ^ 6) := by
    rw [map_pow, map_sub, map_one, uA]
  rw [hd6] at h1
  have hA1 : (1 - univA) ∈ tateDenoms :=
    Submonoid.subset_closure (Or.inl (Or.inl (Or.inl rfl)))
  have hmem : ((1 - univA) ^ 6 * (d : TateBase)) ∈ tateDenoms :=
    tateDenoms.mul_mem (tateDenoms.pow_mem hA1 6) (tateDenomsE_le d.2)
  have h3 : tateDefectTrunc n uA uW (uA * uW)
      * algebraMap TateBase TateUniv ((⟨(1 - univA) ^ 6 * (d : TateBase), hmem⟩ : tateDenoms)
        : TateBase) = algebraMap TateBase TateUniv P := by
    rw [show (((⟨(1 - univA) ^ 6 * (d : TateBase), hmem⟩ : tateDenoms) : TateBase))
      = (1 - univA) ^ 6 * (d : TateBase) from rfl, map_mul]
    rw [← h1]
    ring
  cases n with
  | zero => simp
  | succ m =>
    exact univQ_pow_dvd (m + 1) (univA_pow_dvd_numerator m P _ h3)
      (univW_pow_dvd_numerator m P _ h3)

/-! ## ★★★★★★★★★`1 − a` の単元性が要らない `I^n` 所属 -/

section Final

variable {R : Type} [CommRing R] {I : Ideal R}

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★★**`(1−a)⁶` 倍した差は `I^n` に入る**——`1 − a` が単元でなくてよい。

★これが原点近傍(`u ≡ 1`)へ渡る唯一の橋である。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateDefectTruncE_mem [IsAdicComplete I R] (n : ℕ) (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (hw : IsUnit (1 - w)) : tateDefectTruncE n a w q ∈ I ^ n := by
  have hqa : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * a) := fun m => by
    rw [haw]; exact isUnit_one_sub (I := I) (pow_succ_mul_mem hq m)
  have hqw : ∀ m : ℕ, IsUnit (1 - (a * w) ^ (m + 1) * w) := fun m => by
    rw [haw]; exact isUnit_one_sub (I := I) (pow_succ_mul_mem hq m)
  set φ := tateSpecializeE a w hw hqa hqw with hφ
  have hA : φ eA = a := by rw [hφ, eA, tateSpecializeE_base, tateEval_A]
  have hW : φ eW = w := by rw [hφ, eW, tateSpecializeE_base, tateEval_W]
  obtain ⟨c, hc⟩ := tateDefectTruncE_univ_dvd n
  have hmap : φ (tateDefectTruncE n eA eW (eA * eW)) = tateDefectTruncE n a w q := by
    rw [map_tateDefectTruncE φ n eA eW (eA * eW) isUnit_one_sub_eW
      (fun m _ => isUnit_one_sub_eQA m) (fun m _ => isUnit_one_sub_eQW m),
      hA, hW, map_mul, hA, hW, haw]
  rw [← hmap, hc, map_mul, map_pow, map_mul, hA, hW, haw]
  exact Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)

end Final

/-! ## ★出典の紐付け(`.src`) -/

def TateUnivE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——1-A を反転しない万有な環)",
    sectionId := "genell-def-3-3" }

def tateDefectTruncE_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——1-a の単元性が要らない I^n 所属)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
