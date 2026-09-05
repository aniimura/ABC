import ABC3.Found.HenselianSplits
import Mathlib.FieldTheory.Finite.Basic

/-!
# Teichmüller 持ち上げ——剰余体の乗法群は単数群へ分裂する

論文にも我々のモデルにも固有でない、**一般の**結果。mathlib へ出せる形で書く
(`Found/ResidueFieldFinite.lean`・`Found/FiniteFieldIrreducible.lean`・
`Found/HenselianSplits.lean` と同じ位置づけ)。

## 何を言っているか

`A` を Henselian 局所環、剰余体 `𝓀_A` が**有限**(位数 `q`)とする。このとき
剰余写像 `Aˣ → 𝓀_Aˣ` は**群としての切断**を持つ:

```
teichmuller : 𝓀_Aˣ →* Aˣ,   residue ∘ teichmuller = id
```

各 `b ∈ 𝓀_Aˣ` に対し `X^{q-1} - 1` の根で剰余が `b` のものが**ただ一つ**
存在する(Hensel)。これが Teichmüller 代表元 `ω(b)`。

## なぜ要るか

`Found/PGC/UnramifiedExtension.lean` で `𝒪[K.carrier]` が Henselian
局所環で剰余体が有限であることを示したので、p 進局所体に対して
`𝒪_K^× ⊇ μ_{q-1}`(位数 `q-1` の巡回群)が得られる。古典的局所類体論の
`𝒪_K^× ≅ μ_{q-1} × (1+𝔪_K)` の第一因子——[pGC] Proposition 1.2 が
`Γ_K^ab` の群構造から `q` を読み取るときの入口。

## 証明の骨格

* `b ≠ 0` なら `b^{q-1} = 1`(有限体)なので `b` は `X^{q-1}-1` の剰余体での根。
* その微分は `(q-1)X^{q-2}`、`q-1 ≡ -1 (mod p)` なので `b ≠ 0` の所で
  **単根**。よって `exists_root_of_residue_root`(Hensel)が持ち上げる。
* 一意性は `IsLocalRing.eq_of_eval_eq_zero_of_not_isUnit_sub`。
* 乗法性はこの一意性から:`ω(b)ω(c)` も `X^{q-1}-1` の根で剰余が `bc`。
-/

namespace ABC3.Found

open Polynomial

/-- **Teichmüller 代表元の存在**——剰余体の `0` でない元は
`X^{q-1}-1` の根に持ち上がる。 -/
theorem exists_teichmullerRep {A : Type*} [CommRing A] [HenselianLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] (b : IsLocalRing.ResidueField A) (hb : b ≠ 0) :
    ∃ a : A, a ^ (Fintype.card (IsLocalRing.ResidueField A) - 1) = 1
      ∧ IsLocalRing.residue A a = b := by
  set q := Fintype.card (IsLocalRing.ResidueField A) with hq
  have hq2 : 2 ≤ q := Fintype.one_lt_card
  set P : Polynomial A := X ^ (q - 1) - 1 with hP
  have hPm : P.Monic := by
    rw [hP]; simpa using (monic_X_pow_sub_C (1 : A) (n := q - 1) (by omega))
  have hPmap : P.map (IsLocalRing.residue A) = X ^ (q - 1) - 1 := by rw [hP]; simp
  have hroot : Polynomial.eval b (P.map (IsLocalRing.residue A)) = 0 := by
    rw [hPmap]
    simp only [Polynomial.eval_sub, Polynomial.eval_pow, Polynomial.eval_X, Polynomial.eval_one]
    rw [FiniteField.pow_card_sub_one_eq_one b hb, sub_self]
  have hcast : ((q - 1 : ℕ) : IsLocalRing.ResidueField A) = -1 := by
    have h1 : ((q : ℕ) : IsLocalRing.ResidueField A) = 0 :=
      FiniteField.cast_card_eq_zero (IsLocalRing.ResidueField A)
    have h2 : 1 ≤ q := by omega
    push_cast [Nat.cast_sub h2, h1]
    ring
  have hder : Polynomial.eval b (Polynomial.derivative (P.map (IsLocalRing.residue A))) ≠ 0 := by
    rw [hPmap]
    simp only [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_one,
      sub_zero, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow, Polynomial.eval_X]
    rw [hcast]
    exact mul_ne_zero (by simp) (pow_ne_zero _ hb)
  obtain ⟨a, ha, hares⟩ := exists_root_of_residue_root P hPm b hroot hder
  refine ⟨a, ?_, hares⟩
  have hz : Polynomial.eval a P = 0 := ha
  rw [hP] at hz
  simp only [Polynomial.eval_sub, Polynomial.eval_pow, Polynomial.eval_X,
    Polynomial.eval_one] at hz
  linear_combination hz

/-- **Teichmüller 代表元の一意性**——`X^{q-1}-1` の根は剰余で決まる
(`IsLocalRing.eq_of_eval_eq_zero_of_not_isUnit_sub`)。 -/
theorem teichmullerRep_unique {A : Type*} [CommRing A] [IsLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] {a₁ a₂ : A}
    (h1 : a₁ ^ (Fintype.card (IsLocalRing.ResidueField A) - 1) = 1)
    (h2 : a₂ ^ (Fintype.card (IsLocalRing.ResidueField A) - 1) = 1)
    (hres : IsLocalRing.residue A a₁ = IsLocalRing.residue A a₂) : a₁ = a₂ := by
  set q := Fintype.card (IsLocalRing.ResidueField A) with hq
  have hq2 : 2 ≤ q := Fintype.one_lt_card
  set P : Polynomial A := X ^ (q - 1) - 1 with hP
  have hev : ∀ a : A, a ^ (q - 1) = 1 → Polynomial.eval a P = 0 := by
    intro a ha
    rw [hP]
    simp only [Polynomial.eval_sub, Polynomial.eval_pow, Polynomial.eval_X, Polynomial.eval_one]
    rw [ha, sub_self]
  have hbne : IsLocalRing.residue A a₁ ≠ 0 := by
    intro hcon
    have hh := congrArg (IsLocalRing.residue A) h1
    rw [map_pow, hcon, map_one, zero_pow (by omega)] at hh
    exact zero_ne_one hh
  have hcast : ((q - 1 : ℕ) : IsLocalRing.ResidueField A) = -1 := by
    have hc1 : ((q : ℕ) : IsLocalRing.ResidueField A) = 0 :=
      FiniteField.cast_card_eq_zero (IsLocalRing.ResidueField A)
    have hc2 : 1 ≤ q := by omega
    push_cast [Nat.cast_sub hc2, hc1]
    ring
  have hunit : IsUnit (Polynomial.eval a₁ (Polynomial.derivative P)) := by
    refine (IsLocalRing.residue_ne_zero_iff_isUnit _).mp ?_
    have hcomp := Polynomial.hom_eval₂ (Polynomial.derivative P) (RingHom.id A)
      (IsLocalRing.residue A) a₁
    simp only [Polynomial.eval₂_id, RingHom.comp_id] at hcomp
    rw [hcomp, ← Polynomial.eval_map, ← Polynomial.derivative_map]
    have hPmap : P.map (IsLocalRing.residue A) = X ^ (q - 1) - 1 := by rw [hP]; simp
    rw [hPmap]
    simp only [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_one,
      sub_zero, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow, Polynomial.eval_X]
    rw [hcast]
    exact mul_ne_zero (by simp) (pow_ne_zero _ hbne)
  refine IsLocalRing.eq_of_eval_eq_zero_of_not_isUnit_sub (hev a₁ h1) (hev a₂ h2) ?_ hunit
  rw [← mem_nonunits_iff, ← IsLocalRing.mem_maximalIdeal, ← IsLocalRing.residue_eq_zero_iff,
    map_sub, hres, sub_self]

/-- Teichmüller 代表元(`A` の元として)。 -/
noncomputable def teichmullerVal (A : Type*) [CommRing A] [HenselianLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] (b : (IsLocalRing.ResidueField A)ˣ) : A :=
  Classical.choose (exists_teichmullerRep (b : IsLocalRing.ResidueField A) b.ne_zero)

theorem teichmullerVal_pow (A : Type*) [CommRing A] [HenselianLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] (b : (IsLocalRing.ResidueField A)ˣ) :
    (teichmullerVal A b) ^ (Fintype.card (IsLocalRing.ResidueField A) - 1) = 1 :=
  (Classical.choose_spec (exists_teichmullerRep (b : IsLocalRing.ResidueField A) b.ne_zero)).1

theorem teichmullerVal_residue (A : Type*) [CommRing A] [HenselianLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] (b : (IsLocalRing.ResidueField A)ˣ) :
    IsLocalRing.residue A (teichmullerVal A b) = (b : IsLocalRing.ResidueField A) :=
  (Classical.choose_spec (exists_teichmullerRep (b : IsLocalRing.ResidueField A) b.ne_zero)).2

/-- **★Teichmüller 持ち上げ**——剰余写像 `Aˣ → 𝓀_Aˣ` の群としての切断。 -/
noncomputable def teichmuller (A : Type*) [CommRing A] [HenselianLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] :
    (IsLocalRing.ResidueField A)ˣ →* Aˣ where
  toFun b := (IsUnit.of_pow_eq_one (teichmullerVal_pow A b)
    (by have : 2 ≤ Fintype.card (IsLocalRing.ResidueField A) := Fintype.one_lt_card; omega)).unit
  map_one' := by
    apply Units.ext
    rw [IsUnit.unit_spec]
    refine teichmullerRep_unique (teichmullerVal_pow A 1) (by simp) ?_
    simp [teichmullerVal_residue A 1]
  map_mul' b c := by
    apply Units.ext
    rw [Units.val_mul, IsUnit.unit_spec, IsUnit.unit_spec, IsUnit.unit_spec]
    refine teichmullerRep_unique (teichmullerVal_pow A (b * c)) ?_ ?_
    · rw [mul_pow, teichmullerVal_pow A b, teichmullerVal_pow A c, one_mul]
    · rw [map_mul, teichmullerVal_residue A (b * c), teichmullerVal_residue A b,
        teichmullerVal_residue A c, Units.val_mul]

/-- **切断であること**——`residue ∘ teichmuller = id`。 -/
theorem residue_teichmuller (A : Type*) [CommRing A] [HenselianLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] (b : (IsLocalRing.ResidueField A)ˣ) :
    Units.map (IsLocalRing.residue A : A →* IsLocalRing.ResidueField A) (teichmuller A b) = b :=
  Units.ext (teichmullerVal_residue A b)

/-- したがって `teichmuller` は単射——`𝓀_Aˣ` は `Aˣ` の部分群として実現される。 -/
theorem teichmuller_injective (A : Type*) [CommRing A] [HenselianLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] :
    Function.Injective (teichmuller A) := by
  intro a b hab
  have h := congrArg (Units.map (IsLocalRing.residue A : A →* IsLocalRing.ResidueField A)) hab
  rwa [residue_teichmuller, residue_teichmuller] at h


/-! ## 分裂する全射は直積を与える(可換群)

`f : G →* H` が群としての切断 `s` を持てば `G ≃* H × ker f`。
可換群なので半直積でなく**直積**。Teichmüller 持ち上げと組み合わせて
`Aˣ ≃* 𝓀_Aˣ × (1 + 𝔪_A)` を出すのに使う。 -/

/-- 切断を持つ準同型は直積分解を与える(可換群)。 -/
noncomputable def prodKerOfRightInverse {G H : Type*} [CommGroup G] [CommGroup H]
    (f : G →* H) (s : H →* G) (hs : Function.RightInverse s f) : G ≃* H × f.ker where
  toFun g := (f g, ⟨g * (s (f g))⁻¹, by
    rw [MonoidHom.mem_ker, map_mul, map_inv, hs (f g), mul_inv_cancel]⟩)
  invFun q := s q.1 * (q.2 : G)
  left_inv g := by
    simp only
    rw [mul_comm, mul_assoc, inv_mul_cancel, mul_one]
  right_inv q := by
    obtain ⟨h, k, hk⟩ := q
    have hfk : f (k : G) = 1 := hk
    apply Prod.ext
    · simp only
      rw [map_mul, hs h, hfk, mul_one]
    · apply Subtype.ext
      simp only
      rw [map_mul, hs h, hfk, mul_one, mul_comm (s h) (k : G), mul_assoc, mul_inv_cancel, mul_one]
  map_mul' a b := by
    apply Prod.ext
    · simp
    · apply Subtype.ext
      show a * b * (s (f (a * b)))⁻¹ = (a * (s (f a))⁻¹) * (b * (s (f b))⁻¹)
      rw [map_mul, map_mul, mul_inv]
      simp only [mul_assoc, mul_left_comm]

end ABC3.Found
