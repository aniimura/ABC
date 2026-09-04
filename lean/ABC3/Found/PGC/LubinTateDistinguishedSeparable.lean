import ABC3.Found.PGC.LubinTatePsiNorm

/-!
# `D_n` は分離的——`Λ_n` は真に `q^n` 個の相異なる元(`sorry` 無し)

古典的な Lubin-Tate 理論の核心的な帰結——`n`-段の全捩れ点を統べる
distinguished 多項式 `D_n`(次数 `q^n`)が**分離的**であり、したがって
`K` の代数閉包 `K.closure` の中に**ちょうど `q^n` 個の相異なる根**を
持つこと——を、実在の p進局所体 `K:PAdicLocalField p` について確立する。

## 証明の筋

`D_n = D_{n-1}・ψ_n`(`iteratedLubinTateDistinguished_eq_mul_psi`)・
`D_0 = X`(`iteratedLubinTateDistinguished_zero`)という漸化式を、次の
不変量を保つ帰納法で押し進める:

> `P(n) := D_n(K.carrier へ写したもの)は分離的、かつ
>   `n < j`(`j≥1`)ならば `D_n` と `ψ_j`(同じく写したもの)は互いに素

- 基底(`n=0`): `D_0=X` は分離的(自明)。`X` と `ψ_j`(`j≥1`)が互いに素
  であることは、`ψ_j` の定数項が非零(`0<‖ψ_j.coeff0‖=‖π‖`、
  `norm_iteratedLubinTatePsi_coeff_zero`・`norm_pi_pos_lt_one`)ことから
  `X∤ψ_j` として出る(`X_isCoprime_iteratedLubinTatePsi`)。
- 帰納段: `D_{n+1}=D_n・ψ_{n+1}` について、
  - 分離性は `D_n` の分離性(IH)・`ψ_{n+1}` の分離性(既出、混標数)・
    `IsCoprime D_n ψ_{n+1}`(IH の後半、`j:=n+1`)から
    `Polynomial.Separable.mul` で出る。
  - `j>n+1` に対する `IsCoprime D_{n+1} ψ_j` は、`IsCoprime D_n ψ_j`
    (IH)と `IsCoprime ψ_{n+1} ψ_j`(`n+1≠j` なので既出の
    `isCoprime_iteratedLubinTatePsi`)から `IsCoprime.mul_left` で出る。

これで `D_n` の分離性が確立され、mathlib の `Polynomial.nodup_roots`・
`IsAlgClosed.card_roots_eq_natDegree` と組み合わせて、`K.closure` の中に
`D_n` の根がちょうど `q^n` 個、しかも**互いに相異なる**ことが従う——
古典的な Lubin-Tate 理論の `|Λ_n|=q^n`(真の集合として)そのもの。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-- `X` と `ψ_j`(`j≥1`)は互いに素——`ψ_j` の定数項が非零
(`0<‖ψ_j.coeff0‖=‖π‖`)ことから `X∤ψ_j`。 -/
theorem X_isCoprime_iteratedLubinTatePsi {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (j : ℕ) (hj : 1 ≤ j) :
    IsCoprime (Polynomial.X : Polynomial K.carrier) (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf j hj)) := by
  haveI := valuationRing_isDVR K
  haveI : UniqueFactorizationMonoid (𝒪[K.carrier]) := uniqueFactorizationMonoid_valuationRing K
  rw [Polynomial.irreducible_X.coprime_iff_not_dvd, Polynomial.X_dvd_iff]
  intro hcon
  have hnormcoeff := norm_iteratedLubinTatePsi_coeff_zero K hq hπmax hπne0 f hf0 hf1 hf j hj
  have hcoeff0 : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf j hj)).coeff 0 =
      ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf j hj).coeff 0 : K.carrier) := by
    rw [Polynomial.coeff_map]; rfl
  rw [hcoeff0] at hcon
  have hcoeffzero : (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf j hj).coeff 0 = 0 := by
    exact_mod_cast hcon
  rw [hcoeffzero] at hnormcoeff
  simp only [norm_zero] at hnormcoeff
  obtain ⟨hπpos, _⟩ := norm_pi_pos_lt_one K hπmax hπne0
  rw [← hnormcoeff] at hπpos
  exact lt_irrefl 0 hπpos

/-- ★★★★★★★★★★★★★**`D_n` は分離的、かつ将来の段の `ψ_j`(`j>n`)とも
互いに素**——`D_n=D_{n-1}・ψ_n`・`D_0=X` の漸化式を、この2つの性質を
同時に保つ帰納法で押し進める。 -/
theorem separable_and_coprime_iteratedLubinTateDistinguished_map {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).Separable ∧
    ∀ (j : ℕ) (hj : 1 ≤ j), n < j → IsCoprime (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
        (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n))
      (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
        (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf j hj)) := by
  induction n with
  | zero =>
    constructor
    · rw [iteratedLubinTateDistinguished_zero]
      simp [Polynomial.separable_X]
    · intro j hj _
      rw [iteratedLubinTateDistinguished_zero]
      simp only [Polynomial.map_X]
      exact X_isCoprime_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf j hj
  | succ n ih =>
    obtain ⟨ihsep, ihcop⟩ := ih
    have hn1 : 1 ≤ n + 1 := by omega
    have hcop_n_n1 := ihcop (n + 1) hn1 (by omega)
    have heqmul : Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
        (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n + 1)) =
        Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
          (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) *
        Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
          (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1) := by
      rw [← Polynomial.map_mul]
      congr 1
      have := iteratedLubinTateDistinguished_eq_mul_psi hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1
      simpa using this
    constructor
    · rw [heqmul]
      exact ihsep.mul
        (separable_iteratedLubinTatePsi_map_carrier K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1) hcop_n_n1
    · intro j hj hlt
      rw [heqmul]
      apply IsCoprime.mul_left
      · exact ihcop j hj (by omega)
      · exact isCoprime_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf hn1 hj (by omega)

/-- ★★★★★★★★★★★★★**`D_n`(`K.carrier` へ写したもの)は分離的**。 -/
theorem separable_iteratedLubinTateDistinguished_map {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).Separable :=
  (separable_and_coprime_iteratedLubinTateDistinguished_map K hq hπmax hπne0 f hf0 hf1 hf n).1

/-- `D_n`(`K.closure` へ写したもの)の根は互いに相異なる。 -/
theorem nodup_roots_iteratedLubinTateDistinguished_map {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).roots.Nodup := by
  apply Polynomial.nodup_roots
  rw [show algebraMap (𝒪[K.carrier]) K.closure =
      (algebraMap K.carrier K.closure).comp (algebraMap (𝒪[K.carrier]) K.carrier) from
    IsScalarTower.algebraMap_eq (𝒪[K.carrier]) K.carrier K.closure, ← Polynomial.map_map]
  exact (separable_iteratedLubinTateDistinguished_map K hq hπmax hπne0 f hf0 hf1 hf n).map

/-- `D_n`(`K.closure` へ写したもの)の根の個数はちょうど `q^n`
——`IsAlgClosed.card_roots_eq_natDegree` と `natDegree_
iteratedLubinTateDistinguished` から。`nodup_roots_iteratedLubinTate
Distinguished_map` と合わせて、`K.closure` の中に `D_n` の根が**ちょうど
`q^n` 個、しかも互いに相異なる**ことが従う——古典的な Lubin-Tate 理論の
`|Λ_n|=q^n`(真の集合として)そのもの。 -/
theorem card_roots_iteratedLubinTateDistinguished_map {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    Multiset.card (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).roots = (pp ^ ff) ^ n := by
  rw [IsAlgClosed.card_roots_eq_natDegree,
    (isDistinguishedAt_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).monic.natDegree_map,
    natDegree_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n]


open scoped Classical in
/-- `n`-段の全捩れ点 `Λ_n`——`D_n`(`K.closure` へ写したもの)の根の
`Finset`。`card_iteratedLubinTateTorsionPoints` によりその濃度はちょうど `q^n`。 -/
noncomputable def iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) : Finset K.closure :=
  (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
    (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).roots.toFinset

open scoped Classical in
/-- ★★★★★★★★★★★★★**`|Λ_n| = q^n`**——`iteratedLubinTateTorsionPoints` の濃度。
`card_roots_iteratedLubinTateDistinguished_map`(根の個数)と
`nodup_roots_iteratedLubinTateDistinguished_map`(重複無し)から
`Multiset.toFinset_card_of_nodup` で従う。 -/
theorem card_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    (iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n).card = (pp ^ ff) ^ n := by
  rw [iteratedLubinTateTorsionPoints, Multiset.toFinset_card_of_nodup
    (nodup_roots_iteratedLubinTateDistinguished_map K hq hπmax hπne0 f hf0 hf1 hf n),
    card_roots_iteratedLubinTateDistinguished_map K hq hπmax hπne0 f hf0 hf1 hf n]

open scoped Classical in
/-- `0 ∈ Λ_n`——`D_n` の定数項が `0` であること
(`iteratedLubinTateDistinguished_coeff_zero`)から。 -/
theorem zero_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) : (0 : K.closure) ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
  haveI := valuationRing_isDVR K
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
  refine ⟨(isDistinguishedAt_iteratedLubinTateDistinguished
    hq hπmax hπne0 f hf0 hf1 hf n).monic.map _ |>.ne_zero, ?_⟩
  show (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).eval 0 = 0
  rw [Polynomial.eval_map, Polynomial.eval₂_at_zero,
    iteratedLubinTateDistinguished_coeff_zero hq hπmax hπne0 f hf0 hf1 hf n, map_zero]

open scoped Classical in
/-- `Λ_n ⊆ Λ_m`(`n≤m`)——`D_n ∣ D_m`(`iteratedLubinTateDistinguished_
dvd_of_le`)から `D_n` の根がすべて `D_m` の根であることが
`Polynomial.roots.le_of_dvd` で従う。 -/
theorem iteratedLubinTateTorsionPoints_subset_of_le {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {n m : ℕ} (hnm : n ≤ m) :
    iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n ⊆ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m := by
  haveI := valuationRing_isDVR K
  intro x hx
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset] at hx ⊢
  have hdvd : Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) ∣
      Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
        (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf m) :=
    Polynomial.map_dvd _ (iteratedLubinTateDistinguished_dvd_of_le hq hπmax hπne0 f hf0 hf1 hf hnm)
  have hmne0 : Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf m) ≠ 0 :=
    (isDistinguishedAt_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf m).monic.map _ |>.ne_zero
  exact Multiset.mem_of_le (Polynomial.roots.le_of_dvd hmne0 hdvd) hx
end ABC3.Found.PGC
