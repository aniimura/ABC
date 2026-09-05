import ABC3.Found.PGC.QpNonAbelian
import ABC3.Found.PGC.InertiaIdentification

/-!
# [pGC] Corollary 1.3 の**旧**形は偽だった——`SC` が自由だと `I_K` が非正規になる

`Skeleton/PGC/Section1Cor13.lean` の設計メモは、こう書いていた:

> 反証には**非正規な** `inertia` が要り(`Check.PGC.map_conj_of_normal`)、
> それには `Γ_K` の非 Galois な3次以上の拡大を構成する必要がある。

`Found/PGC/QpNonAbelian.lean` でそれを構成した(`ℚ_p(p^{1/ℓ})` は normal で
ない)ので、**反証が届いた**。

## 反例(`cor_1_3_statement_false`、`sorry` 無し)

`H₀ := Gal(K̄/ℚ_p(p^{1/ℓ}))` は開で指数 `ℓ ≥ 2`、そして**正規でない**
(`exists_nonNormal_openSubgroup`)。`RD` は**実物** `residueCardinality p` の
まま、`SC` だけを病的に取る:

```
SC.field K H _ := if (HEq H H₀ ∧ H ≠ ⊤) then A else K
```

ここで `A` は剰余体の元の個数が `q^{[Γ:H₀]}` の p進局所体
(= 次数 `[Γ:H₀]` の不分岐拡大、`exists_residueCard_pow`)。すると判定条件
`q_{L_H} = q^{[Γ_K:H]}` を満たす開部分群はちょうど `{⊤, H₀}` になり
(`inertia_badSC`)、

```
inertia RD SC K = ⊤ ⊓ H₀ = H₀
```

——**正規でない**。ところが `RecoverableFromAbsGal` は任意の
`α : Γ_K ≅ Γ_K`(特に内部自己同型)で `I_K` が保たれることを要求するので、
`I_K` は正規でなければならない。矛盾。

★`SC.field_top` は満たされている(`H = ⊤` の分岐は `H ≠ ⊤` が偽なので
`else` に落ちる)ことに注意——原文が課している条件は全部守ったうえで偽になる。

## 修理

`Found/` は既に**実物**を構成している:
`residueCardinality p`(`ResidueCardinalityConstruction.lean`)と
`subgroupCorrespondence p`(`SubgroupCorrespondenceConstruction.lean`)。
そして実物の下では

```
inertia (residueCardinality p) (subgroupCorrespondence p) K = Gal(K̄/K^ur)
```

であり(`InertiaIdentification.lean::inertia_eq_fixingSubgroup_unramifiedClosure`)、
これは**正規**(`normal_inertia`)。したがって上の反例は塞がる。
`Skeleton/PGC/Section1Cor13.lean::inertia_recoverable` は実物について述べる形に
直した——設計メモの「Track B が本物の `ResidueCardinality` を構成した時点で、
ここに依存する全ての statement が一斉に非空虚性の検査を受ける」がまさに
起きたということ。

★これで「落とした条件は、主張を偽にするか自明にするかのどちらかになる」例は
6 つ目。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC
open scoped Valued Classical

variable {p : ℕ} [Fact p.Prime]

/-- 剰余体の元の個数が `q^n` の p進局所体(= 次数 `n` の不分岐拡大)。 -/
theorem exists_residueCard_pow (K : PAdicLocalField p) (n : ℕ) (hn : n ≠ 0) :
    ∃ A : PAdicLocalField p,
      (residueCardinality p).card A = (Nat.card 𝓀[K.carrier]) ^ n := by
  obtain ⟨z, hrank, hu, -, -⟩ := exists_isUnramifiedAdjoin K n hn
  refine ⟨adjoinField K z, ?_⟩
  rw [residueCardinality_adjoinField, residueDegree_eq_residueCard_pow K z,
    inertiaDegree_eq_finrank_of_isUnramified K z hu, hrank]

/-- **非正規な開部分群が存在する**——`ℚ_p(p^{1/ℓ})` の固定部分群
(`Found/PGC/QpNonAbelian.lean::not_normal_adjoin`)。 -/
theorem exists_nonNormal_openSubgroup (p : ℕ) [Fact p.Prime] :
    ∃ (H : Subgroup (selfField p).absGal) (_ : IsOpen (H : Set (selfField p).absGal)),
      2 ≤ H.index ∧ ¬ H.Normal := by
  haveI := isGalois_closure (selfField p)
  obtain ⟨ℓ, hℓ, hℓgt⟩ := exists_prime_gt p
  obtain ⟨x, hxpow, hrank⟩ := exists_root_pow_eq_p p hℓ
  refine ⟨(IntermediateField.adjoin (selfField p).carrier
    ({x} : Set (selfField p).closure)).fixingSubgroup,
    IntermediateField.fixingSubgroup_isOpen _, ?_, ?_⟩
  · rw [index_fixingSubgroup_adjoin, hrank]
    exact hℓ.two_le
  · intro hnormal
    have hgal := (InfiniteGalois.normal_iff_isGalois _).mp hnormal
    haveI := hgal
    exact not_normal_adjoin p hℓ hℓgt hxpow hrank inferInstance

/-- 病的な「開部分群と中間体の対応」——`H₀` にだけ別の体を返す。
`field_top` は満たしている。 -/
noncomputable def badSC (p : ℕ) [Fact p.Prime]
    (H₀ : Subgroup (selfField p).absGal) (A : PAdicLocalField p) : SubgroupCorrespondence p where
  field K H _ := if (HEq H H₀ ∧ H ≠ ⊤) then A else K
  field_top K h := by simp

theorem badSC_field_of_ne (p : ℕ) [Fact p.Prime]
    (H₀ : Subgroup (selfField p).absGal) (A : PAdicLocalField p)
    {K : PAdicLocalField p} {H : Subgroup K.absGal} (hH : IsOpen (H : Set K.absGal))
    (hne : ¬ (HEq H H₀ ∧ H ≠ ⊤)) :
    (badSC p H₀ A).field K H hH = K := if_neg hne

theorem badSC_field_H₀ (p : ℕ) [Fact p.Prime]
    (H₀ : Subgroup (selfField p).absGal) (A : PAdicLocalField p)
    (hH : IsOpen (H₀ : Set (selfField p).absGal)) (hne : H₀ ≠ ⊤) :
    (badSC p H₀ A).field (selfField p) H₀ hH = A := if_pos ⟨HEq.rfl, hne⟩

/-- **判定条件を満たす開部分群はちょうど `{⊤, H₀}`**——だから `I_K = H₀`。 -/
theorem inertia_badSC (p : ℕ) [Fact p.Prime]
    {H₀ : Subgroup (selfField p).absGal} (hopen : IsOpen (H₀ : Set (selfField p).absGal))
    (hidx : 2 ≤ H₀.index) {A : PAdicLocalField p}
    (hA : (residueCardinality p).card A
      = (Nat.card 𝓀[(selfField p).carrier]) ^ H₀.index) :
    inertia (residueCardinality p) (badSC p H₀ A) (selfField p) = H₀ := by
  have hq2 : 2 ≤ Nat.card 𝓀[(selfField p).carrier] := by
    have h1 : 1 < Fintype.card 𝓀[(selfField p).carrier] := Fintype.one_lt_card
    rw [Nat.card_eq_fintype_card]
    omega
  have hnetop : H₀ ≠ ⊤ := by
    intro h
    rw [h, Subgroup.index_top] at hidx
    omega
  have hmemS : H₀ ∈ {H : Subgroup (selfField p).absGal |
      ∃ hH : IsOpen (H : Set (selfField p).absGal),
        IsUnramifiedAt (residueCardinality p) (badSC p H₀ A) (selfField p) H hH} := by
    refine ⟨hopen, ?_⟩
    show (residueCardinality p).card ((badSC p H₀ A).field (selfField p) H₀ hopen)
      = (residueCardinality p).card (selfField p) ^ H₀.index
    rw [badSC_field_H₀ p H₀ A hopen hnetop, hA, residueCardinality_card]
  refine le_antisymm (sInf_le hmemS) (le_sInf ?_)
  rintro H ⟨hH, hcrit⟩
  by_cases hcase : (HEq H H₀ ∧ H ≠ ⊤)
  · rw [eq_of_heq hcase.1]
  · have hfield : (badSC p H₀ A).field (selfField p) H hH = selfField p :=
      badSC_field_of_ne p H₀ A hH hcase
    have hc : (residueCardinality p).card (selfField p)
        = (residueCardinality p).card (selfField p) ^ H.index := by
      have h := hcrit
      rw [IsUnramifiedAt, hfield] at h
      exact h
    rw [residueCardinality_card] at hc
    have hidx1 : H.index = 1 := by
      rcases Nat.lt_or_ge H.index 2 with h1 | h1
      · by_cases h2 : H.index = 0
        · rw [h2, pow_zero] at hc; omega
        · omega
      · exfalso
        have hlt : (Nat.card 𝓀[(selfField p).carrier]) ^ 1
            < (Nat.card 𝓀[(selfField p).carrier]) ^ H.index :=
          Nat.pow_lt_pow_right (by omega) (by omega)
        rw [pow_one] at hlt
        omega
    rw [Subgroup.index_eq_one.mp hidx1]
    exact le_top

/-- 内部自己同型は連続な群同型。 -/
noncomputable def conjCE (K : PAdicLocalField p) (c : K.absGal) :
    ContinuousMulEquiv K.absGal K.absGal where
  toFun g := c * g * c⁻¹
  invFun g := c⁻¹ * g * c
  left_inv g := by group
  right_inv g := by group
  map_mul' a b := by group
  continuous_toFun := by fun_prop
  continuous_invFun := by fun_prop

/-- **★★★★★★★[pGC] Corollary 1.3 の旧形(自由な `SC`)は偽**。 -/
theorem cor_1_3_statement_false (p : ℕ) [Fact p.Prime] :
    ¬ (∀ (RD : ResidueCardinality p) (SC : SubgroupCorrespondence p),
        (inertiaObject RD SC).RecoverableFromAbsGal) := by
  intro h
  obtain ⟨H₀, hopen, hidx, hnn⟩ := exists_nonNormal_openSubgroup p
  obtain ⟨A, hA⟩ := exists_residueCard_pow (selfField p) H₀.index (by omega)
  have hI := inertia_badSC p hopen hidx hA
  refine hnn ⟨fun n hn g => ?_⟩
  have hrec := h (residueCardinality p) (badSC p H₀ A) (conjCE (selfField p) g)
  have hrec' : Subgroup.map (conjCE (selfField p) g).toMulEquiv.toMonoidHom H₀ = H₀ := by
    rw [← hI]; exact hrec
  have hmem : g * n * g⁻¹ ∈ Subgroup.map (conjCE (selfField p) g).toMulEquiv.toMonoidHom H₀ :=
    ⟨n, hn, rfl⟩
  rwa [hrec'] at hmem

end ABC3.Check.PGC
