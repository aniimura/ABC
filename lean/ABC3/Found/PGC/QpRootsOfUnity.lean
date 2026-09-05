import ABC3.Found.PGC.PrimeToPTorsion
import ABC3.Found.PGC.QpResidueField

/-!
# `ℚ_p` の中の `ℓ` 乗根は自明——そして `X^ℓ - p` は既約

本ファイルは「`Γ_{ℚ_p}` は非可換」を作るための2本の柱を用意する。

* **`eq_one_of_pow_eq_one_qp`**: `ℓ` が素数で `ℓ ≠ p`・`ℓ ∤ p-1` なら、
  `ℚ_p` の中の `ℓ` 乗根は `1` だけ。
  ——`PrimeToPTorsion.lean::card_primeToPTorsion`(`p` と素な捩れ部分群の
  位数がちょうど `q-1`)の**最初の応用**。`ℚ_p` では `q = p` なので `p-1`。
* **`irreducible_X_pow_sub_C_p`**: `X^ℓ - p` は `ℚ_p` 上既約。
  ——`ℚ_p` の値群が `p^ℤ` なので `b^ℓ = p` は `ℓ·v(b) = 1` を要求し、
  `ℓ ≥ 2` では不可能(`not_isPow_p`)。あとは mathlib の
  `X_pow_sub_C_irreducible_iff_of_prime`。

## なぜ要るか(次の段)

`ℓ > p` なる素数を取れば `ℓ ≠ p` かつ `ℓ ∤ p-1`。すると
`L := ℚ_p(p^{1/ℓ})` は次数 `ℓ` の拡大で、**normal でない**:
normal なら `X^ℓ - p` が `L` で分解し、2根の比として原始 `ℓ` 乗根 `ζ` が
`L` に入る。`[ℚ_p(ζ):ℚ_p]` は `ℓ`(素数)を割り、かつ `ζ` は
`∑_{i<ℓ} X^i`(次数 `ℓ-1`)の根なので `ℓ` 未満——ゆえに `1`、
すなわち `ζ ∈ ℚ_p`。上の第一の柱で `ζ = 1`、矛盾。

normal でない有限拡大があれば、その固定部分群は**正規でない開部分群**で
あり、`Γ_{ℚ_p}` は非可換である。これは `Skeleton/PGC/Section2.lean::prop_2_2`
(`IntKbar`・`CompKbar` が自由な型族)の反証に必要な材料
——`Check/PGC/RefutationAttempts.lean` が「共役による反証には非正規な部分群が
要る」と記録していた穴。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found Polynomial
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## `p` と素な捩れから: `ℓ` 乗根の自明性 -/

theorem eq_one_of_pow_eq_one_real {x : ℝ} {n : ℕ} (hn : n ≠ 0) (hx : 0 ≤ x) (h : x ^ n = 1) :
    x = 1 := by
  rcases lt_trichotomy x 1 with hlt | heq | hgt
  · exact absurd (h ▸ pow_lt_one₀ hx hlt hn) (lt_irrefl 1)
  · exact heq
  · exact absurd (h ▸ one_lt_pow₀ hgt hn) (lt_irrefl 1)

/-- `‖z‖ = 1` の元を `𝒪_K` の単数として見る。 -/
noncomputable def unitOfNormOne' (K : PAdicLocalField p) {z : K.carrier} (hz : ‖z‖ = 1) :
    (𝒪[K.carrier])ˣ := by
  have hmem : z ∈ 𝒪[K.carrier] := by
    rw [Valuation.mem_integer_iff]
    have hv : Valued.v z = (‖z‖₊ : NNReal) := NNReal.eq rfl
    rw [hv]
    have h1 : ‖z‖ ≤ 1 := le_of_eq hz
    exact_mod_cast h1
  refine (?_ : IsUnit (⟨z, hmem⟩ : 𝒪[K.carrier])).unit
  rw [Valued.integer.isUnit_iff_norm_eq_one]
  exact hz

@[simp] theorem unitOfNormOne'_val (K : PAdicLocalField p) {z : K.carrier} (hz : ‖z‖ = 1) :
    (((unitOfNormOne' K hz : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) = z := rfl

/-- **`ℓ` 乗根が自明になる十分条件**——`ℓ` 素数、`ℓ ≠ p`、`ℓ ∤ q-1`。
`p` と素な捩れ部分群の位数が `q-1`(`card_primeToPTorsion`)なので、
位数が `ℓ` を割りかつ `q-1` を割る元は自明。 -/
theorem eq_one_of_pow_eq_one_of_norm_one (K : PAdicLocalField p) {ℓ : ℕ} (hℓ : ℓ.Prime)
    (hℓp : ℓ ≠ p) (hnd : ¬ (ℓ ∣ (Nat.card 𝓀[K.carrier] - 1)))
    {z : K.carrier} (hz : ‖z‖ = 1) (hpow : z ^ ℓ = 1) : z = 1 := by
  set u := unitOfNormOne' K hz with hu
  have hupow : u ^ ℓ = 1 := by
    apply Units.ext
    apply Subtype.ext
    show ((((u ^ ℓ : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier)) = 1
    rw [Units.val_pow_eq_pow_val]
    show (((u : 𝒪[K.carrier]) : K.carrier)) ^ ℓ = 1
    rw [hu, unitOfNormOne'_val K hz, hpow]
  have hnpd : ¬ (p ∣ ℓ) := by
    intro h
    exact hℓp (((Nat.prime_dvd_prime_iff_eq Fact.out hℓ).mp h).symm)
  have hmem : u ∈ primeToPTorsion K := ⟨ℓ, hnpd, hupow⟩
  set v : primeToPTorsion K := ⟨u, hmem⟩ with hv
  have hvpow : v ^ ℓ = 1 := Subtype.ext hupow
  have hd1 : orderOf v ∣ ℓ := orderOf_dvd_of_pow_eq_one hvpow
  have hd2 : orderOf v ∣ Nat.card (primeToPTorsion K) := orderOf_dvd_natCard v
  rw [card_primeToPTorsion K] at hd2
  have hone : orderOf v = 1 := by
    rcases (Nat.Prime.eq_one_or_self_of_dvd hℓ _ hd1) with h | h
    · exact h
    · exact absurd (h ▸ hd2) hnd
  have hv1 : v = 1 := orderOf_eq_one_iff.mp hone
  have hu1 : u = 1 := congrArg Subtype.val hv1
  have h2 : z = (((1 : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) := by
    rw [← unitOfNormOne'_val K hz, ← hu, hu1]
  simpa using h2

theorem card_primeToPTorsion_selfField (p : ℕ) [Fact p.Prime] :
    Nat.card (primeToPTorsion (selfField p)) = p - 1 := by
  rw [card_primeToPTorsion]
  rw [show Nat.card 𝓀[(selfField p).carrier] = residueCard (selfField p) from rfl,
    residueCard_selfField]

/-- **★★★★★`ℚ_p` の中の `ℓ` 乗根は自明**(`ℓ` 素数、`ℓ ≠ p`、`ℓ ∤ p-1`)。 -/
theorem eq_one_of_pow_eq_one_qp (p : ℕ) [Fact p.Prime] {ℓ : ℕ} (hℓ : ℓ.Prime)
    (hℓp : ℓ ≠ p) (hnd : ¬ (ℓ ∣ (p - 1))) {ζ : ℚ_[p]} (hζ : ζ ^ ℓ = 1) : ζ = 1 := by
  have hnorm : @norm ℚ_[p] _ ζ = 1 :=
    eq_one_of_pow_eq_one_real hℓ.ne_zero (norm_nonneg _) (by rw [← norm_pow, hζ, norm_one])
  have hz := (norm_selfField p ζ).trans hnorm
  have hcard : Nat.card 𝓀[(selfField p).carrier] - 1 = p - 1 := by
    rw [show Nat.card 𝓀[(selfField p).carrier] = residueCard (selfField p) from rfl,
      residueCard_selfField]
  exact eq_one_of_pow_eq_one_of_norm_one (selfField p) hℓ hℓp (by rw [hcard]; exact hnd) hz hζ

/-! ## 値群からの障害: `X^ℓ - p` は既約 -/

/-- **`ℚ_p` の値群は `p^ℤ`**——`n ≥ 2` なら `p` は `n` 乗にならない
(`n · v(b) = 1` は整数解を持たない)。 -/
theorem not_isPow_p (p : ℕ) [Fact p.Prime] {n : ℕ} (hn : 2 ≤ n) (b : ℚ_[p]) :
    b ^ n ≠ (p : ℚ_[p]) := by
  intro h
  have hp1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hb : b ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega)] at h
    have hz : ‖(p : ℚ_[p])‖ = 0 := by rw [← h, norm_zero]
    rw [Padic.norm_p] at hz
    have hpos : (0:ℝ) < (p:ℝ)⁻¹ := by positivity
    linarith
  have hn1 : ‖b ^ n‖ = ‖(p : ℚ_[p])‖ := by rw [h]
  rw [norm_pow, Padic.norm_eq_zpow_neg_valuation hb, Padic.norm_p] at hn1
  rw [← zpow_natCast ((p:ℝ) ^ (-b.valuation)) n, ← zpow_mul] at hn1
  have hinv : ((p : ℝ))⁻¹ = (p : ℝ) ^ (-1 : ℤ) := by rw [zpow_neg, zpow_one]
  rw [hinv] at hn1
  have hv := zpow_right_injective₀ (by positivity) (ne_of_gt hp1) hn1
  have h2 : b.valuation * (n : ℤ) = 1 := by linear_combination -hv
  have hdvd : (n : ℤ) ∣ 1 := Dvd.intro_left _ h2
  have hle : (n : ℤ) ≤ 1 := Int.le_of_dvd one_pos hdvd
  omega

/-- **★★★★★`X^ℓ - p` は `ℚ_p` 上既約**(`ℓ` 素数)。 -/
theorem irreducible_X_pow_sub_C_p (p : ℕ) [Fact p.Prime] {ℓ : ℕ} (hℓ : ℓ.Prime) :
    Irreducible ((X : Polynomial ℚ_[p]) ^ ℓ - C (p : ℚ_[p])) :=
  (X_pow_sub_C_irreducible_iff_of_prime hℓ).mpr (fun b => not_isPow_p p hℓ.two_le b)

end ABC3.Found.PGC
