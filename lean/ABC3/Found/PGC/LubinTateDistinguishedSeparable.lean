import ABC3.Found.PGC.LubinTatePsiNorm
import ABC3.Found.PGC.AdjoinPAdicLocalField

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

open scoped Classical in
/-- ★★★★★★★★★★**`Λ_n` の元のスペクトルノルムは `1` 未満**——`𝒪_K` 加群
構造の構成(`[a]_f` による作用)へ向けた前段: `D_n` の根はすべて
「位相的冪零」であることの言い換え。`D_n=D_{n-1}・ψ_n`・`D_0=X` の
漸化式に沿った帰納法: `D_0=X` の根は `0`(`spectralNorm 0=0`)。
`D_{n+1}` の根は `D_n` の根(IH)か `ψ_{n+1}` の根で、後者は
`spectralNorm_root_iteratedLubinTatePsi` により
`‖π‖^(1/(q^{n+1}-q^n))` に等しく、`0<‖π‖<1` かつ指数が正なので `1`
未満(`Real.rpow_lt_one`)。 -/
theorem spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    spectralNorm K.carrier K.closure x < 1 := by
  haveI := valuationRing_isDVR K
  induction n with
  | zero =>
    rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hx
    obtain ⟨_, hxroot⟩ := hx
    rw [iteratedLubinTateDistinguished_zero] at hxroot
    simp only [Polynomial.map_X, Polynomial.IsRoot.def, Polynomial.eval_X] at hxroot
    rw [hxroot, spectralNorm_zero]
    exact zero_lt_one
  | succ n ih =>
    rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hx
    obtain ⟨hne, hxroot⟩ := hx
    have hn1 : 1 ≤ n + 1 := by omega
    have heqmul : Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
        (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n + 1)) =
        Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
          (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) *
        Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
          (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1) := by
      rw [← Polynomial.map_mul]
      congr 1
      have := iteratedLubinTateDistinguished_eq_mul_psi hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1
      simpa using this
    rw [heqmul] at hxroot
    rw [Polynomial.IsRoot.def, Polynomial.eval_mul, mul_eq_zero] at hxroot
    rcases hxroot with h1 | h2
    · apply ih
      rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
      refine ⟨?_, h1⟩
      exact (isDistinguishedAt_iteratedLubinTateDistinguished
        hq hπmax hπne0 f hf0 hf1 hf n).monic.map _ |>.ne_zero
    · have haeval : Polynomial.aeval x (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
          (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1)) = 0 := by
        rw [Polynomial.aeval_def, Polynomial.eval₂_eq_eval_map, Polynomial.map_map,
          ← IsScalarTower.algebraMap_eq (𝒪[K.carrier]) K.carrier K.closure]
        exact h2
      have hspec := spectralNorm_root_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf
        (n + 1) hn1 x haeval
      rw [hspec]
      obtain ⟨hπpos, hπlt1⟩ := norm_pi_pos_lt_one K hπmax hπne0
      have hdegpos : (0 : ℝ) <
          ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1).natDegree : ℝ) := by
        have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
        have hlt : (pp ^ ff) ^ n < (pp ^ ff) ^ (n + 1) := Nat.pow_lt_pow_right h2 (by omega)
        rw [natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1]
        simp only [Nat.add_sub_cancel]
        exact_mod_cast (by omega : 0 < (pp ^ ff) ^ (n + 1) - (pp ^ ff) ^ n)
      exact Real.rpow_lt_one hπpos.le hπlt1 (by positivity)

open scoped Classical in
/-- `Λ_n` の元は `𝒪_K` 上整——`D_n`(`𝒪_K` 係数、モニック)の根なので、
`D_n` 自身が`IsIntegral`の定義に要求されるモニック多項式そのものと
なる。 -/
theorem isIntegral_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    IsIntegral (𝒪[K.carrier]) x := by
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hx
  obtain ⟨hne, hxroot⟩ := hx
  refine ⟨iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n,
    (isDistinguishedAt_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).monic, ?_⟩
  rw [Polynomial.eval₂_eq_eval_map]
  exact hxroot

/-- `Λ_n` の元は `K.carrier`(局所体本体)上でも整——`𝒪_K` 上の整性
(`isIntegral_of_mem_iteratedLubinTateTorsionPoints`)から
`IsIntegral.tower_top` で中間の環 `K.carrier` へ持ち上げる。有限次拡大
`K.carrier⟮x⟯`(`IntermediateField.adjoin.finiteDimensional`)を経由して、
将来 `x` での冪級数評価に必要な完備性(有限次拡大は完備)へ繋がる。 -/
theorem isIntegral_carrier_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    IsIntegral K.carrier x :=
  (isIntegral_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx).tower_top

open scoped Classical in
/-- `Λ_n` の元 `x` を添加した単純拡大 `K.carrier⟮x⟯` は `K.carrier`
上有限次——`isIntegral_carrier_of_mem_iteratedLubinTateTorsionPoints`
(`x` は `K.carrier` 上整)から
`IntermediateField.adjoin.finiteDimensional` で直接従う。`K.closure`
自体は完備でないので、`[a]_f` の実際の評価はこの有限次拡大の中で
行う、というのが見通している次の一歩。 -/
theorem finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  IntermediateField.adjoin.finiteDimensional
    (isIntegral_carrier_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx)

open scoped Classical in
/-- ★★★★★★★★★★**`Λ_n` の元を添加した単純拡大 `K.carrier⟮x⟯` は完備**
——`K.carrier` 上有限次
(`finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints`)・
`K.carrier` 自身が完備(`ℚ_[p]` の有限次拡大)・`K.closure` にスペクトル
ノルムで与えたノルム体構造(`LocalFieldNorm.lean::closureNormedField`
等)を組み合わせて、mathlib の一般論
`FiniteDimensional.complete`(完備な体上の有限次ノルム空間は完備)を
適用する。これで初めて `PowerSeries.aeval` による `[a]_f` の `x` での
評価が可能になる舞台が整った。 -/
theorem completeSpace_adjoin_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    CompleteSpace (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  haveI := finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n x hx
  FiniteDimensional.complete K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))

/-- ★★★★★★★★★★**`Λ_n` の元は位相的冪零(`PowerSeries.HasEval`)**——
`spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints`(ノルムが
`1`未満)と mathlib の
`tendsto_pow_atTop_nhds_zero_of_norm_lt_one`(ノルム`1`未満の元の冪は
`0`へ収束)を組み合わせるだけ。`PowerSeries.HasEval`は定義的に
`IsTopologicallyNilpotent`なので、これで `PowerSeries.aeval` に渡す
「評価点の条件」の**片方**(位相的冪零性)が揃った。

★注記(2026-09-04 の発見、次の一歩の見通し): 残るもう片方の条件
`IsLinearTopology S S`(評価先 `S` 自身の位相がイデアルで生成される
「線形位相」であること)は、体 `K.closure` や `K.carrier⟮x⟯` それ自身
では**成り立たない**(体のイデアルは `{0}` と全体しかないので、非自明な
位相と両立しない)。この条件は mathlib では
`Ideal.isLinearTopology`(付値環の極大イデアルによる adic 位相)から
出るのが自然で、古典的な Lubin-Tate 理論でも `[a]_f` は実際には
「体」でなく「付値環」の間の写像として評価される。見通している次の
一歩: `L:=K.carrier⟮x⟯` は `K.carrier` 上有限次(かつ `ℚ_[p]` 上も
有限次、推移律で)なので、**`L` 自身を新たな `PAdicLocalField p` として
再構成**すれば、`Found/PGC/LocalFieldNorm.lean`・
`Found/PGC/ValuationRingDVR.lean`・`ValuationRingComplete.lean` の
機構(`𝒪[L]`・`Valued L`・その adic 完備性)がそのまま流用でき、
`IsLinearTopology 𝒪[L] 𝒪[L]` が `Ideal.isLinearTopology` から従う
はず。 -/
theorem hasEval_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    PowerSeries.HasEval x :=
  tendsto_pow_atTop_nhds_zero_of_norm_lt_one
    (spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx)

/-- `Λ_n` の元 `x` を添加した単純拡大 `K.carrier⟮x⟯` 自身を新たな
`PAdicLocalField p` として得る(`AdjoinPAdicLocalField.lean::
adjoinPAdicLocalField` の torsion-point 版)——必要な
`FiniteDimensional K.carrier _` は
`finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints` から。
これで `Found/PGC/LocalFieldNorm.lean` 等が `K` に与えていた機構
(`𝒪[·]`・`Valued`・adic 完備性・`Ideal.isLinearTopology` の土台)が
この新しい局所体へそのまま流用できる——ただし
`AdjoinPAdicLocalField.lean` の docstring に記した通り、この新しい
ノルム構造と `K.closure` から部分体として継承されるノルム構造が
definitionally 一致するとは限らないので、両者の一致は今後の課題。 -/
noncomputable def torsionPointAdjoinPAdicLocalField {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    PAdicLocalField p :=
  haveI := finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n x hx
  adjoinPAdicLocalField K x

/-- ★★★★★★★★★★**訂正・前進**: `AdjoinPAdicLocalField.lean` の docstring
で「`K.carrier⟮x⟯` に `LocalFieldNorm.lean` のノルムを適用すると
instance diamond になりうる」と記した懸念は、**単純拡大の体そのもの
のレベルでは杞憂だった**——`IntermediateField.adjoin K.carrier {x}` が
`K.closure` の部分体として mathlib の一般論から自動的に持つ
`NormedField` 構造は、**そのままの定義で**
`‖(⟨x,_⟩ : K.carrier⟮x⟯)‖ = spectralNorm K.carrier K.closure x` を
`rfl` で満たす(部分体の `NormedField` は単に周囲のノルムの制限として
定義されているため)。したがって `spectralNorm_lt_one_of_mem_
iteratedLubinTateTorsionPoints` が直接この体自身のノルムについての
事実になり、`Λ_n` の元を `K.carrier⟮x⟯` の元として見たときも
`PowerSeries.HasEval`(位相的冪零性)が成り立つ。

(`adjoinPAdicLocalField`——`K.carrier⟮x⟯` を新たな `PAdicLocalField p`
として`ℚ_[p]`から再構成する経路——はこの体レベルの事実には不要
だったが、`IsLinearTopology`(評価先が付値**環**であることを要求する
`PowerSeries.aeval` の残る条件)を得るには依然として有用な見通しで
ある——体でなく付値環の話であることに変わりはない。) -/
theorem hasEval_mem_adjoin_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    PowerSeries.HasEval (⟨x, hmem⟩ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  apply tendsto_pow_atTop_nhds_zero_of_norm_lt_one
  show spectralNorm K.carrier K.closure x < 1
  exact spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx

/-! ## `Λ_n = Λ_{n-1} ∪ (ψ_n の根)`——「原始的な」π^n 捩れ点の切り出し

`D_n = D_{n-1}・ψ_n`(`iteratedLubinTateDistinguished_eq_mul_psi`)を
`K.closure` へ写した上で `Polynomial.roots_mul` を適用するだけで、
`Λ_n`(`iteratedLubinTateTorsionPoints`)が `Λ_{n-1}` と「`ψ_n` の根の
なす集合」の(集合としての)和集合に分解されることが分かる——
`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×` へ向けての足がかり:「原始的な」
π^n-捩れ点(`Λ_n` から `Λ_{n-1}` を除いたもの)がちょうど `ψ_n` の根に
対応する、という古典的な事実の Finset レベルでの定式化。 -/

open scoped Classical in
/-- `ψ_n`(`K.closure` へ写したもの)の根のなす `Finset`——
「原始的な」π^n-捩れ点全体。`iteratedLubinTateTorsionPoints`(`Λ_n`)
と同じパターン(`Multiset.toFinset`)で定義する。 -/
noncomputable def iteratedLubinTatePsiTorsionPoints {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) : Finset K.closure :=
  (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)).roots.toFinset

/-- `ψ_n`(`K.closure` へ写したもの)の根は相異なる——`D_n` の分離性の
証明(`separable_iteratedLubinTatePsi_map_carrier`、`LubinTatePsiNorm.lean`
既出)を `card_roots_iteratedLubinTateDistinguished_map` と同じ手筋
(`algebraMap` の分解・`Polynomial.map_map`)で `K.closure` レベルへ運ぶ。 -/
theorem nodup_roots_iteratedLubinTatePsi_map {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)).roots.Nodup := by
  apply Polynomial.nodup_roots
  rw [show algebraMap (𝒪[K.carrier]) K.closure =
      (algebraMap K.carrier K.closure).comp (algebraMap (𝒪[K.carrier]) K.carrier) from
    IsScalarTower.algebraMap_eq (𝒪[K.carrier]) K.carrier K.closure, ← Polynomial.map_map]
  exact (separable_iteratedLubinTatePsi_map_carrier K hq hπmax hπne0 f hf0 hf1 hf n hn).map

/-- `ψ_n`(`K.closure` へ写したもの)の根の個数は `q^n-q^{n-1}`——
代数閉体上の分離多項式の根の個数(`IsAlgClosed.card_roots_eq_
natDegree`)と次数の公式(`natDegree_iteratedLubinTatePsi`、既出)から。 -/
theorem card_roots_iteratedLubinTatePsi_map {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    Multiset.card (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)).roots
      = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  rw [IsAlgClosed.card_roots_eq_natDegree,
    (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic.natDegree_map,
    natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]

open scoped Classical in
/-- ★★★★★★★★★★**`Λ_n = Λ_{n-1} ∪ (ψ_n の根)`**——`D_n=D_{n-1}・ψ_n`
(`iteratedLubinTateDistinguished_eq_mul_psi`)を `K.closure` へ写して
`Polynomial.roots_mul`(積の根はそれぞれの根の合併、multiset の加法)
を適用し、`Multiset.toFinset` の言葉に戻すだけ。`Gal(K(Λ_n)/K)≅
(𝒪_K/π^n)^×` へ向けて、「原始的な」π^n-捩れ点(`ψ_n` の根)を
`Λ_n` の中に位置づける最初の一歩。 -/
theorem iteratedLubinTateTorsionPoints_eq_union {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n =
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n - 1) ∪
        iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn := by
  have hmul : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)) =
      (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
        (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1))) *
      (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
        (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) := by
    rw [← Polynomial.map_mul]
    congr 1
    exact iteratedLubinTateDistinguished_eq_mul_psi hq hπmax hπne0 f hf0 hf1 hf n hn
  ext y
  simp only [iteratedLubinTateTorsionPoints, iteratedLubinTatePsiTorsionPoints,
    Multiset.mem_toFinset, Finset.mem_union]
  rw [hmul, Polynomial.roots_mul, Multiset.mem_add]
  rw [← hmul]
  exact (isDistinguishedAt_iteratedLubinTateDistinguished
    hq hπmax hπne0 f hf0 hf1 hf n).monic.map _ |>.ne_zero

open scoped Classical in
/-- `ψ_n` の根のなす `Finset` の濃度は `q^n-q^{n-1}`——
`nodup_roots_iteratedLubinTatePsi_map`(重複なし)と
`card_roots_iteratedLubinTatePsi_map`(multiset としての濃度)から。
`|Λ_n|=|Λ_{n-1}|+|(ψ_n の根)|`(`card_iteratedLubinTateTorsionPoints`
の `q^n=q^{n-1}+(q^n-q^{n-1})` としての分解)の後半にあたる。 -/
theorem card_iteratedLubinTatePsiTorsionPoints {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn).card =
      (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  rw [iteratedLubinTatePsiTorsionPoints, Multiset.toFinset_card_of_nodup
    (nodup_roots_iteratedLubinTatePsi_map K hq hπmax hπne0 f hf0 hf1 hf n hn),
    card_roots_iteratedLubinTatePsi_map K hq hπmax hπne0 f hf0 hf1 hf n hn]

end ABC3.Found.PGC
