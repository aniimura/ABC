import ABC3.Found.PGC.RamifiedUnramifiedDisjoint
import ABC3.Found.PGC.PrimeToPTorsion

/-!
# `K^ur` は `p` と素な 1 の冪根をすべて含む

古典的な事実:**`m` が `p` と素なら `μ_m ⊆ K^ur`**。

理由は初等的で、`K^ur` の剰余体が `𝔽̄_p` だからだが、ここでは
既に持っている道具だけで組む:

1. `k := φ(m)`(Euler)とすると `m ∣ q^k − 1`(`q = |𝓀_K|` は `p` と素)
2. 次数 `k` の不分岐拡大 `A = K(x)` を取る
   (`exists_isUnramifiedAdjoin`)——その剰余体の元は `q^k` 個
3. `primeToPTorsion A ≃* 𝓀_A^×` は位数 `q^k − 1` の巡回群
   (`Found/PGC/PrimeToPTorsion.lean`、第 962)
4. 巡回群は位数の約数ごとに元を持つので、位数ちょうど `m` の元が取れる
5. `A = K(x) ≤ K^ur`(`adjoin_le_unramifiedClosure_of_isUnramifiedAdjoin`、第 999)

★これは**局所類体論を経由しない**——`q` を群論的に読む「暴分岐を割った
惰性群 `I_K/P_K` に Frobenius が `q` 倍で作用する」という道(Prop 1.2 への
LCFT を通らない経路)の土台。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped Valued

/-! ## 補助(一般) -/

/-- 有限巡回群は位数の約数ごとに、位数ちょうどその値の元を持つ。 -/
theorem exists_orderOf_eq_of_dvd_card {G : Type*} [Group G] [Finite G] [IsCyclic G] {m : ℕ}
    (hdvd : m ∣ Nat.card G) : ∃ g : G, orderOf g = m := by
  obtain ⟨a, ha⟩ := IsCyclic.exists_generator (α := G)
  have hoa : orderOf a = Nat.card G := orderOf_eq_card_of_forall_mem_zpowers ha
  obtain ⟨c, hc⟩ := hdvd
  have hNpos : 0 < Nat.card G := Nat.card_pos
  have hcpos : 0 < c := by
    rcases Nat.eq_zero_or_pos c with h | h
    · rw [h, Nat.mul_zero] at hc; omega
    · exact h
  refine ⟨a ^ c, ?_⟩
  rw [orderOf_pow, hoa, hc]
  have hg : Nat.gcd (m * c) c = c := (Nat.gcd_comm _ _).trans (Nat.gcd_mul_left_right m c)
  rw [hg, Nat.mul_div_cancel _ hcpos]

/-- Euler の定理の言い換え——`q` と `m` が素なら `m ∣ q^{φ(m)} − 1`。 -/
theorem dvd_pow_totient_sub_one {q m : ℕ} (hq : 0 < q) (h : Nat.Coprime q m) :
    m ∣ q ^ Nat.totient m - 1 :=
  (Nat.modEq_iff_dvd' (Nat.one_le_pow _ _ hq)).mp (Nat.ModEq.pow_totient h).symm

variable {p : ℕ} [Fact p.Prime]

/-- `primeToPTorsion A` から `K.closure` への単射準同型
(`A` は `K.closure` の中間体)。 -/
noncomputable def primeToPTorsionToClosure (K : PAdicLocalField p) (x : K.closure) :
    (primeToPTorsion (adjoinField K x)) →* K.closure where
  toFun u := algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure
    (((u : (𝒪[(adjoinField K x).carrier])ˣ) : 𝒪[(adjoinField K x).carrier]) :
      (adjoinField K x).carrier)
  map_one' := rfl
  map_mul' _ _ := rfl

theorem primeToPTorsionToClosure_injective (K : PAdicLocalField p) (x : K.closure) :
    Function.Injective (primeToPTorsionToClosure K x) := by
  intro a b h
  have h' : (((a : (𝒪[(adjoinField K x).carrier])ˣ) : 𝒪[(adjoinField K x).carrier]) :
      (adjoinField K x).carrier)
      = (((b : (𝒪[(adjoinField K x).carrier])ˣ) : 𝒪[(adjoinField K x).carrier]) :
      (adjoinField K x).carrier) :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      K.closure).injective h
  exact Subtype.ext (Units.ext (Subtype.ext h'))

/-- **★★★★★★★★★★★★★★★★★★`p` と素な `m` に対し、位数 `m` の 1 の冪根が
`K^ur` にある**。 -/
theorem exists_isPrimitiveRoot_mem_unramifiedClosure (K : PAdicLocalField p) (m : ℕ)
    (hm : 0 < m) (hmp : ¬ p ∣ m) :
    ∃ ζ : K.closure, ζ ∈ unramifiedClosure K ∧ IsPrimitiveRoot ζ m := by
  obtain ⟨ff, hffpos, hqf⟩ := (residueCardinality p).isPrimePow K
  rw [residueCardinality_card] at hqf
  set q := Nat.card 𝓀[K.carrier] with hqdef
  have hq2 : 2 ≤ q := by
    have h1 : 1 < Fintype.card 𝓀[K.carrier] := Fintype.one_lt_card
    rw [hqdef, Nat.card_eq_fintype_card]; omega
  have hcop : Nat.Coprime q m := by
    rw [hqf]
    exact Nat.Coprime.pow_left _ ((Nat.Prime.coprime_iff_not_dvd Fact.out).mpr hmp)
  set k := Nat.totient m with hkdef
  have hkpos : 0 < k := Nat.totient_pos.mpr hm
  have hdvd : m ∣ q ^ k - 1 := dvd_pow_totient_sub_one (by omega) hcop
  obtain ⟨x, hrank, hu, -, -⟩ := exists_isUnramifiedAdjoin K k hkpos.ne'
  have hAcard : Nat.card 𝓀[(adjoinField K x).carrier] = q ^ k := by
    rw [← residueCardinality_card, residueCardinality_adjoinField,
      residueDegree_eq_residueCard_pow K x, inertiaDegree_eq_finrank_of_isUnramified K x hu,
      hrank]
  have hcardT : Nat.card (primeToPTorsion (adjoinField K x)) = q ^ k - 1 := by
    rw [card_primeToPTorsion, hAcard]
  have hqk2 : 2 ≤ q ^ k := le_trans hq2 (Nat.le_self_pow hkpos.ne' q)
  haveI : Finite (primeToPTorsion (adjoinField K x)) := by
    have hpos : 0 < Nat.card (primeToPTorsion (adjoinField K x)) := by rw [hcardT]; omega
    exact (Nat.card_pos_iff.mp hpos).2
  haveI := isCyclic_primeToPTorsion (adjoinField K x)
  obtain ⟨u, hu'⟩ := exists_orderOf_eq_of_dvd_card (G := primeToPTorsion (adjoinField K x))
    (m := m) (by rw [hcardT]; exact hdvd)
  have hord : orderOf (primeToPTorsionToClosure K x u) = m := by
    rw [orderOf_injective (primeToPTorsionToClosure K x)
      (primeToPTorsionToClosure_injective K x) u]
    exact hu'
  refine ⟨primeToPTorsionToClosure K x u, ?_, ?_, ?_⟩
  · apply adjoin_le_unramifiedClosure_of_isUnramifiedAdjoin K hu
    show algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure
      (((u : (𝒪[(adjoinField K x).carrier])ˣ) : 𝒪[(adjoinField K x).carrier]) :
        (adjoinField K x).carrier) ∈ _
    exact SetLike.coe_mem _
  · rw [← hord]; exact pow_orderOf_eq_one _
  · intro l hl
    rw [← hord]
    exact orderOf_dvd_of_pow_eq_one hl

/-- **★★★★★★★★★★★★★★★★★★★★`μ_m ⊆ K^ur`(`p ∤ m`)**——`m` 乗して `1` に
なる元は**すべて** `K^ur` にある。

原始 `m` 乗根 `ζ₀ ∈ K^ur` を取れば、`X^m − 1` の根は `ζ₀` の冪で尽くされる
(`IsPrimitiveRoot.eq_pow_of_pow_eq_one`)ので、`K^ur` が体であることから従う。 -/
theorem mem_unramifiedClosure_of_pow_eq_one (K : PAdicLocalField p) {m : ℕ}
    (hm : 0 < m) (hmp : ¬ p ∣ m) {ζ : K.closure} (hζ : ζ ^ m = 1) :
    ζ ∈ unramifiedClosure K := by
  haveI : NeZero m := ⟨hm.ne'⟩
  obtain ⟨ζ₀, hζ₀mem, hζ₀⟩ := exists_isPrimitiveRoot_mem_unramifiedClosure K m hm hmp
  obtain ⟨i, -, hi⟩ := hζ₀.eq_pow_of_pow_eq_one hζ
  rw [← hi]
  exact pow_mem hζ₀mem i

end ABC3.Found.PGC
