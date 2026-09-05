import ABC3.Found.PGC.QpRootsOfUnity
import ABC3.Found.PGC.UnramifiedSubextension

/-!
# `Γ_{ℚ_p}` は非可換——`ℚ_p(p^{1/ℓ})` は normal でない

`Check/PGC/RefutationAttempts.lean` は「共役による反証には**非正規な開部分群**が
要る」と記録していた。本ファイルはそれを作る。

## 筋

`ℓ > p` なる素数を取る(`exists_prime_gt`)。すると `ℓ ≠ p` かつ `ℓ ∤ p-1`。

1. `X^ℓ - p` は `ℚ_p` 上既約(`QpRootsOfUnity.lean::irreducible_X_pow_sub_C_p`
   ——値群が `p^ℤ` だから `b^ℓ = p` は不可能)。その根 `x` は
   `[ℚ_p(x):ℚ_p] = ℓ` を与える(`exists_root_pow_eq_p`)。
2. `ℚ_p(x)` の中の `ℓ` 乗根は `1` だけ(`eq_one_of_pow_eq_one_mem_adjoin`):
   `ζ ≠ 1` なら `ζ` は `∑_{i<ℓ} X^i`(次数 `ℓ-1`)の根なので
   `[ℚ_p(ζ):ℚ_p] ≤ ℓ-1`。一方それは `ℓ`(素数)を割る
   (`finrank_dvd_of_adjoin_le`)ので `1`、つまり `ζ ∈ ℚ_p`。
   すると `QpRootsOfUnity.lean::eq_one_of_pow_eq_one_qp`(`|μ^{(p')}| = p-1`)で
   `ζ = 1`——矛盾。
3. ゆえに `Gal(ℚ_p(x)/ℚ_p)` は自明(`gal_adjoin_trivial`)。
   `σ x / x` は `ℓ` 乗根なので `1`。
4. `ℚ_p(x)` が normal なら(標数 0 なので)Galois で
   `|Gal| = [ℚ_p(x):ℚ_p] = ℓ ≥ 2`——矛盾。**normal でない**。
5. `Γ_{ℚ_p}` が可換なら任意の部分群が正規、特に `ℚ_p(x)` の固定部分群が正規で、
   `InfiniteGalois.normal_iff_isGalois` により `ℚ_p(x)/ℚ_p` は Galois、
   したがって normal——矛盾。**`Γ_{ℚ_p}` は非可換**。

## 何に使うか

`Skeleton/PGC/Section2.lean::prop_2_2` は `IntKbar`・`CompKbar` を**自由な
型族**(しかも `SMul` は公理なしの自由なインスタンス)で受け取っている。
非中心元 `g₀` と `c g₀ c⁻¹ ≠ g₀` なる `c` があれば、
`IntKbar K := ℤ`・`g • n := if g = g₀ then 0 else n`・`α := c による共役`
で反例が作れる。本ファイルはその `g₀`・`c` の存在を与える。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC Polynomial
open scoped Valued

/-! ## 補助 -/

theorem exists_prime_gt (p : ℕ) : ∃ ℓ : ℕ, ℓ.Prime ∧ p < ℓ := by
  obtain ⟨ℓ, hle, hℓ⟩ := Nat.exists_infinite_primes (p + 1)
  exact ⟨ℓ, hℓ, by omega⟩

theorem not_dvd_sub_one_of_gt {p ℓ : ℕ} (hp : 2 ≤ p) (hℓ : p < ℓ) : ¬ (ℓ ∣ (p - 1)) := by
  intro h
  have h1 : 0 < p - 1 := by omega
  have h2 := Nat.le_of_dvd h1 h
  omega

/-! ## `X^ℓ - p` の根 -/

/-- `X^ℓ - p` の根 `x` は `[ℚ_p(x):ℚ_p] = ℓ` を与える。 -/
theorem exists_root_pow_eq_p (p : ℕ) [Fact p.Prime] {ℓ : ℕ} (hℓ : ℓ.Prime) :
    ∃ x : (selfField p).closure,
      x ^ ℓ = algebraMap (selfField p).carrier (selfField p).closure
          (((p : ℕ) : (selfField p).carrier))
      ∧ Module.finrank (selfField p).carrier
          (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) = ℓ := by
  set K := selfField p with hK
  set f : Polynomial K.carrier := X ^ ℓ - C (((p : ℕ) : K.carrier)) with hf
  have hmonic : f.Monic := monic_X_pow_sub_C _ hℓ.ne_zero
  have hirr : Irreducible f := irreducible_X_pow_sub_C_p p hℓ
  have hdeg : f.natDegree = ℓ := natDegree_X_pow_sub_C
  have hdeg0 : f.degree ≠ 0 := by
    rw [degree_eq_natDegree hmonic.ne_zero, hdeg]
    exact_mod_cast Nat.cast_ne_zero.mpr hℓ.ne_zero
  obtain ⟨x, hx⟩ := IsAlgClosed.exists_root (k := K.closure)
    (f.map (algebraMap K.carrier K.closure)) (by rw [degree_map]; exact hdeg0)
  have haeval : (aeval x) f = 0 := by
    rw [aeval_def, ← eval_map]
    exact hx
  have hmin : minpoly K.carrier x = f := (minpoly.eq_of_irreducible_of_monic hirr haeval hmonic).symm
  have hint : IsIntegral K.carrier x := ⟨f, hmonic, haeval⟩
  refine ⟨x, ?_, ?_⟩
  · have h := haeval
    rw [hf] at h
    simp only [map_sub, map_pow, aeval_X, aeval_C] at h
    linear_combination h
  · rw [IntermediateField.adjoin.finrank hint, hmin, hdeg]

/-! ## `ℚ_p(x)` の中の `ℓ` 乗根は自明 -/

/-- **★★★★★`[ℚ_p(x):ℚ_p] = ℓ` なら `ℚ_p(x)` の中の `ℓ` 乗根は `1` だけ**。 -/
theorem eq_one_of_pow_eq_one_mem_adjoin (p : ℕ) [Fact p.Prime] {ℓ : ℕ} (hℓ : ℓ.Prime)
    (hℓgt : p < ℓ) {x : (selfField p).closure}
    (hrank : Module.finrank (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) = ℓ)
    {ζ : (selfField p).closure}
    (hζE : ζ ∈ IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))
    (hζ : ζ ^ ℓ = 1) : ζ = 1 := by
  by_contra hne
  have hgeom : ∑ i ∈ Finset.range ℓ, ζ ^ i = 0 := by
    have h := geom_sum_mul ζ ℓ
    rw [hζ, sub_self] at h
    rcases mul_eq_zero.mp h with h1 | h1
    · exact h1
    · exact absurd h1 (sub_ne_zero.mpr hne)
  have hg0 : (∑ i ∈ Finset.range ℓ, (X : Polynomial (selfField p).carrier) ^ i) ≠ 0 := by
    intro h0
    have hev : ((∑ i ∈ Finset.range ℓ, (X : Polynomial (selfField p).carrier) ^ i).eval 0) = 1 := by
      rw [eval_finsetSum, Finset.sum_eq_single 0]
      · simp
      · intro b _ hb; simp [zero_pow hb]
      · intro hb; exact absurd (Finset.mem_range.mpr hℓ.pos) hb
    rw [h0] at hev
    simp at hev
  have hgdeg : (∑ i ∈ Finset.range ℓ, (X : Polynomial (selfField p).carrier) ^ i).degree
      ≤ ((ℓ - 1 : ℕ) : WithBot ℕ) := by
    refine (degree_sum_le _ _).trans (Finset.sup_le ?_)
    intro i hi
    rw [Finset.mem_range] at hi
    refine (degree_X_pow_le i).trans ?_
    exact_mod_cast Nat.cast_le.mpr (Nat.le_sub_one_of_lt hi)
  have haeval :
      (aeval ζ) (∑ i ∈ Finset.range ℓ, (X : Polynomial (selfField p).carrier) ^ i) = 0 := by
    rw [map_sum]
    simpa using hgeom
  have hint : IsIntegral (selfField p).carrier ζ :=
    IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic _)
  have hle : Module.finrank (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({ζ} : Set (selfField p).closure))
      ≤ ℓ - 1 := by
    rw [IntermediateField.adjoin.finrank hint]
    have hd := minpoly.degree_le_of_ne_zero (selfField p).carrier ζ hg0 haeval
    have hmne : minpoly (selfField p).carrier ζ ≠ 0 := minpoly.ne_zero hint
    rw [degree_eq_natDegree hmne] at hd
    exact_mod_cast hd.trans hgdeg
  have hsub : IntermediateField.adjoin (selfField p).carrier ({ζ} : Set (selfField p).closure)
      ≤ IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure) :=
    IntermediateField.adjoin_simple_le_iff.mpr hζE
  have hdvd := finrank_dvd_of_adjoin_le (selfField p) hsub
  rw [hrank] at hdvd
  have hpos : 0 < Module.finrank (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({ζ} : Set (selfField p).closure)) :=
    Module.finrank_pos
  have hd1 : Module.finrank (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({ζ} : Set (selfField p).closure)) = 1 := by
    rcases (Nat.Prime.eq_one_or_self_of_dvd hℓ _ hdvd) with h | h
    · exact h
    · omega
  have hbot := IntermediateField.finrank_eq_one_iff.mp hd1
  have hmem : ζ ∈ (⊥ : IntermediateField (selfField p).carrier (selfField p).closure) := by
    rw [← hbot]
    exact IntermediateField.mem_adjoin_simple_self _ _
  obtain ⟨c, hc⟩ := IntermediateField.mem_bot.mp hmem
  have hcpow : c ^ ℓ = 1 := by
    have h1 : (algebraMap (selfField p).carrier (selfField p).closure) (c ^ ℓ)
        = (algebraMap (selfField p).carrier (selfField p).closure) 1 := by
      rw [map_pow, hc, hζ, map_one]
    exact (algebraMap (selfField p).carrier (selfField p).closure).injective h1
  have hc1 : c = 1 :=
    eq_one_of_pow_eq_one_qp p hℓ (by omega)
      (not_dvd_sub_one_of_gt (Fact.out : p.Prime).two_le hℓgt) hcpow
  rw [hc1, map_one] at hc
  exact hne hc.symm

/-! ## `Gal(ℚ_p(x)/ℚ_p)` は自明 -/

/-- **★★★★★`Gal(ℚ_p(x)/ℚ_p) = 1`**——`σ x / x` は `ℓ` 乗根なので `1`。 -/
theorem gal_adjoin_trivial (p : ℕ) [Fact p.Prime] {ℓ : ℕ} (hℓ : ℓ.Prime) (hℓgt : p < ℓ)
    {x : (selfField p).closure}
    (hxpow : x ^ ℓ = algebraMap (selfField p).carrier (selfField p).closure
      (((p : ℕ) : (selfField p).carrier)))
    (hrank : Module.finrank (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) = ℓ)
    (σ : (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))
      ≃ₐ[(selfField p).carrier]
        (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))) :
    σ = 1 := by
  have hpne : ((p : ℕ) : (selfField p).carrier) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Fact.out : p.Prime).ne_zero
  have hx0 : x ≠ 0 := by
    intro h0
    rw [h0, zero_pow hℓ.ne_zero] at hxpow
    exact hpne ((map_eq_zero_iff _ (algebraMap _ _).injective).mp hxpow.symm)
  have hxEpow : (⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩ :
      IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) ^ ℓ
      = algebraMap (selfField p).carrier _ (((p : ℕ) : (selfField p).carrier)) := by
    apply Subtype.ext
    show x ^ ℓ = _
    exact hxpow
  have hσpow : ((σ ⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩ :
      IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))
      : (selfField p).closure) ^ ℓ = x ^ ℓ := by
    have h1 : (σ ⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩) ^ ℓ
        = (⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩ :
          IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) ^ ℓ := by
      rw [← map_pow, hxEpow, AlgEquiv.commutes]
    have h2 := congrArg (fun z :
      IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure) =>
        ((z : (selfField p).closure))) h1
    simpa using h2
  have hζpow : (((σ ⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩ :
      IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))
      : (selfField p).closure) / x) ^ ℓ = 1 := by
    rw [div_pow, hσpow, div_self (pow_ne_zero ℓ hx0)]
  have hζE : (((σ ⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩ :
      IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))
      : (selfField p).closure) / x)
      ∈ IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure) :=
    div_mem (σ ⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩).2
      (IntermediateField.mem_adjoin_simple_self _ _)
  have hζ1 := eq_one_of_pow_eq_one_mem_adjoin p hℓ hℓgt hrank hζE hζpow
  refine algEquiv_eq_one_of_fixes_gen (selfField p) x σ (Subtype.ext ?_)
  show ((σ ⟨x, IntermediateField.mem_adjoin_simple_self (selfField p).carrier x⟩ :
    IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))
    : (selfField p).closure) = x
  rwa [div_eq_one_iff_eq hx0] at hζ1

/-! ## normal でない拡大、そして非可換性 -/

/-- **★★★★★★`ℚ_p(x)` は normal でない**。 -/
theorem not_normal_adjoin (p : ℕ) [Fact p.Prime] {ℓ : ℕ} (hℓ : ℓ.Prime) (hℓgt : p < ℓ)
    {x : (selfField p).closure}
    (hxpow : x ^ ℓ = algebraMap (selfField p).carrier (selfField p).closure
      (((p : ℕ) : (selfField p).carrier)))
    (hrank : Module.finrank (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) = ℓ) :
    ¬ Normal (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) := by
  intro hnorm
  haveI := hnorm
  haveI : CharZero (selfField p).carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] (selfField p).carrier).injective
  haveI : Algebra.IsSeparable (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) :=
    IntermediateField.isSeparable_tower_bot (selfField p).carrier _
  haveI : IsGalois (selfField p).carrier
      (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure)) := ⟨⟩
  have hcard := IsGalois.card_aut_eq_finrank (selfField p).carrier
    (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))
  have hsub : Subsingleton ((IntermediateField.adjoin (selfField p).carrier
      ({x} : Set (selfField p).closure)) ≃ₐ[(selfField p).carrier]
      (IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure))) :=
    ⟨fun a b => by
      rw [gal_adjoin_trivial p hℓ hℓgt hxpow hrank a, gal_adjoin_trivial p hℓ hℓgt hxpow hrank b]⟩
  rw [Nat.card_eq_one_iff_unique.mpr ⟨hsub, ⟨1⟩⟩, hrank] at hcard
  exact hℓ.ne_one hcard.symm

/-- **★★★★★★★`Γ_{ℚ_p}` は非可換**——normal でない有限拡大があるから。 -/
theorem not_commutative_absGal (p : ℕ) [Fact p.Prime] :
    ¬ (∀ a b : (selfField p).absGal, a * b = b * a) := by
  intro hcomm
  obtain ⟨ℓ, hℓ, hℓgt⟩ := exists_prime_gt p
  obtain ⟨x, hxpow, hrank⟩ := exists_root_pow_eq_p p hℓ
  haveI := isGalois_closure (selfField p)
  set E := IntermediateField.adjoin (selfField p).carrier ({x} : Set (selfField p).closure) with hE
  haveI hnormal : (E.fixingSubgroup).Normal := by
    refine ⟨fun n hn g => ?_⟩
    have : g * n * g⁻¹ = n := by
      rw [hcomm g n, mul_assoc, mul_inv_cancel, mul_one]
    rw [this]
    exact hn
  have hgal : IsGalois (selfField p).carrier E :=
    (InfiniteGalois.normal_iff_isGalois E).mp hnormal
  haveI := hgal
  exact not_normal_adjoin p hℓ hℓgt hxpow hrank inferInstance

/-- **非中心元が存在する**——`prop_2_2` の反証に必要な `g₀`・`c`。 -/
theorem exists_not_commute_absGal (p : ℕ) [Fact p.Prime] :
    ∃ c g : (selfField p).absGal, c * g * c⁻¹ ≠ g := by
  by_contra h
  refine not_commutative_absGal p (fun a b => ?_)
  have hab : a * b * a⁻¹ = b := not_not.mp (fun hc => h ⟨a, b, hc⟩)
  calc a * b = (a * b * a⁻¹) * a := by group
    _ = b * a := by rw [hab]

end ABC3.Found.PGC
