import ABC3.Found.PGC.UnramifiedCriterion

/-!
# 構成した惰性群は本物——`I_K = Gal(K̄/K^ur)`

`Skeleton/PGC/Section1Cor13.lean` は、原文が判定条件までしか書いていない
ところを補って、惰性群を

```
inertia RD SC K = sInf { H ⊆ Γ_K : H は開、かつ q_{L_H} = q^{[Γ_K:H]} }
```

として**構成した**(★原文に無い段、暗黙の段 #7)。設計メモは
「`RD` が自由なままだと退化する」と警告していた。

`UnramifiedCriterion.lean` で判定条件の正しさを示したので、本ファイルは
その構成が**本物の惰性群と一致する**ことを示す:

```
inertia (residueCardinality p) (subgroupCorrespondence p) K
  = (unramifiedClosure K).fixingSubgroup   ( = Gal(K̄/K^ur) )
```

## 証明

* `≥`(`fixingSubgroup_unramifiedClosure_le_inertia`、前ファイル):
  判定条件を満たす `H` の固定体は `K^ur` に入るので、Galois 対応で
  `Gal(K̄/K^ur) ≤ H`。すべての `H` について成り立つので `sInf` にも。
* `≤`(本ファイル): `z ∈ K^ur` なら、ある不分岐な `x` について
  `z ∈ K(x)`(`mem_unramifiedClosure_iff`)。その `K(x)` の固定部分群は
  **開**(`IntermediateField.fixingSubgroup_isOpen`、`K(x)/K` は有限次)で
  **判定条件を満たす**(`isUnramifiedAt_iff_isUnramifiedAdjoin`)ので
  `sInf` の材料に入っている。ゆえに `g ∈ inertia` は `K(x)` を固定し、
  特に `g z = z`。

★「閉部分群は自分を含む開部分群すべての共通部分」という副有限群の
一般論は**要らなかった**——`K^ur` が単項不分岐拡大の**有向和**である
(`mem_unramifiedClosure_iff`)ことが、その役割をそのまま果たす。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- 不分岐な `K(x)` の固定部分群は、`inertia` の `sInf` の材料に入っている
——開(`K(x)/K` は有限次)で、判定条件を満たす。 -/
theorem isUnramifiedAt_fixingSubgroup_adjoin (K : PAdicLocalField p) {x : K.closure}
    (hu : IsUnramifiedAdjoin K x) :
    ∃ hH : IsOpen
        (((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup)
          : Set K.absGal),
      IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup hH := by
  haveI := isGalois_closure K
  have hopen : IsOpen
      (((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup)
        : Set K.absGal) :=
    IntermediateField.fixingSubgroup_isOpen _
  refine ⟨hopen, ?_⟩
  by_cases hne : (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup = ⊤
  · have hfield : (subgroupCorrespondence p).field K
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup hopen = K := by
      simp [subgroupCorrespondence, hne]
    rw [IsUnramifiedAt, hfield, hne, Subgroup.index_top, pow_one]
  · refine (isUnramifiedAt_iff_isUnramifiedAdjoin K hopen hne x ?_).mpr hu
    exact (InfiniteGalois.fixedField_fixingSubgroup _).symm

/-- **`I_K ≤ Gal(K̄/K^ur)`**——`K^ur` の元は不分岐な単項拡大の中にあり、
その固定部分群が `sInf` の材料だから。 -/
theorem inertia_le_fixingSubgroup_unramifiedClosure (K : PAdicLocalField p) :
    inertia (residueCardinality p) (subgroupCorrespondence p) K
      ≤ (unramifiedClosure K).fixingSubgroup := by
  intro g hg
  rw [IntermediateField.mem_fixingSubgroup_iff]
  intro z hz
  obtain ⟨x, hux, hzx⟩ := (mem_unramifiedClosure_iff K z).mp hz
  obtain ⟨hH, hcrit⟩ := isUnramifiedAt_fixingSubgroup_adjoin K hux
  have hmem : g ∈ (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup :=
    Subgroup.mem_sInf.mp hg _ ⟨hH, hcrit⟩
  rw [IntermediateField.mem_fixingSubgroup_iff] at hmem
  exact hmem z hzx

/-- **★★★★★★構成した惰性群は本物**——
`inertia (residueCardinality p) (subgroupCorrespondence p) K = Gal(K̄/K^ur)`。

`Skeleton` が原文に無い段として補った構成が、実際に古典的な惰性群と
一致することの確定。設計メモが恐れていた「`RD` の自由度による退化」は、
実物 `residueCardinality p` を入れた時点で**起きていない**。 -/
theorem inertia_eq_fixingSubgroup_unramifiedClosure (K : PAdicLocalField p) :
    inertia (residueCardinality p) (subgroupCorrespondence p) K
      = (unramifiedClosure K).fixingSubgroup :=
  le_antisymm (inertia_le_fixingSubgroup_unramifiedClosure K)
    (fixingSubgroup_unramifiedClosure_le_inertia K)

/-- **`I_K` は正規部分群**——`K^ur/K` が Galois だから
(`InfiniteGalois.normal_iff_isGalois`)。

★これは `Check/PGC/RefutationAttempts.lean` の記録に対する答えでもある:
Corollary 1.3 を共役(内部自己同型)で反証するには**非正規な** `inertia`
が要る(`Check.PGC.map_conj_of_normal`)。構成した実物の下では
`inertia` は正規なので、**その反証経路は閉じている**。 -/
theorem normal_inertia (K : PAdicLocalField p) :
    (inertia (residueCardinality p) (subgroupCorrespondence p) K).Normal := by
  haveI := isGalois_closure K
  rw [inertia_eq_fixingSubgroup_unramifiedClosure K]
  exact (InfiniteGalois.normal_iff_isGalois _).mpr (isGalois_unramifiedClosure K)

/-! ## `Γ_K / I_K ≅ Gal(K^ur/K)` -/

/-- `K^ur/K` は normal——`AlgEquiv.restrictNormalHom` を**式に書くために**
instance が要る(`haveI` を証明の中に置くのでは遅い)。 -/
instance instNormalUnramifiedClosure (K : PAdicLocalField p) :
    Normal K.carrier (unramifiedClosure K) := normal_unramifiedClosure K

/-- `K^ur` への制限射の核はちょうど `I_K`
(`IntermediateField.restrictNormalHom_ker` + 惰性群の同定)。 -/
theorem ker_restrictNormalHom_unramifiedClosure (K : PAdicLocalField p) :
    (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (unramifiedClosure K)).ker
      = inertia (residueCardinality p) (subgroupCorrespondence p) K := by
  haveI := normal_unramifiedClosure K
  rw [inertia_eq_fixingSubgroup_unramifiedClosure K]
  exact IntermediateField.restrictNormalHom_ker _

/-- **★★★★`Γ_K / I_K ≅ Gal(K^ur/K)`**——古典的な「惰性群による商は
不分岐 Galois 群」。核が `I_K` であること(上)と制限射の全射性
(`AlgEquiv.restrictNormalHom_surjective`、`K^ur/K` は normal)から。 -/
noncomputable def absGalQuotKerEquivUnramifiedGal (K : PAdicLocalField p) :
    (K.absGal ⧸ (AlgEquiv.restrictNormalHom
        (F := K.carrier) (K₁ := K.closure) (unramifiedClosure K)).ker)
      ≃* (unramifiedClosure K ≃ₐ[K.carrier] unramifiedClosure K) := by
  haveI := normal_unramifiedClosure K
  haveI := IsAlgClosure.normal K.carrier K.closure
  exact QuotientGroup.quotientKerEquivOfSurjective _
    (AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
      (K₁ := (unramifiedClosure K)) (E := K.closure))

/-- **`Γ_K / I_K` は任意の `n ≥ 1` に対して `ℤ/n` へ全射する**——
`Gal(K^ur/K) ≅ Ẑ` の、`Ẑ` を経由しない具体的な言い換え
(`Ẑ` は mathlib に不在、2026-09-05 実測)。 -/
theorem exists_surjective_quotKer_to_zmod (K : PAdicLocalField p) (n : ℕ) (hn : n ≠ 0) :
    ∃ φ : (K.absGal ⧸ (AlgEquiv.restrictNormalHom
        (F := K.carrier) (K₁ := K.closure) (unramifiedClosure K)).ker)
      →* Multiplicative (ZMod n), Function.Surjective φ := by
  obtain ⟨ψ, hψ⟩ := exists_surjective_unramifiedClosureGal_to_zmod K n hn
  exact ⟨ψ.comp (absGalQuotKerEquivUnramifiedGal K).toMonoidHom,
    hψ.comp (absGalQuotKerEquivUnramifiedGal K).surjective⟩

end ABC3.Found.PGC
