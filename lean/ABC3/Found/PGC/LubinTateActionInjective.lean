import ABC3.Found.PGC.AdjoinIntegers
import ABC3.Found.PGC.LubinTateFormalGroupLawEstimate

/-!
# `𝒪_K` 加群の作用の単射性(`sorry` 無し)

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★ **`a·x=b·x ⟹ π^n∣(a-b)`**
——`𝒪_K` の元 `a,b` の Lubin-Tate 作用が原始的な `π^n`-捩れ点 `x` の
上で一致すれば、`a≡b (mod π^n)`。これまで「`F_f` の形式逆元・
対数のどちらかのフルな構成が要る、既存の存在補題と同規模の新規
構築」と判断されていた単射性の壁を、**`F_f` の引き算・逆元を
一切経由せず**に突破する。

## 証明の鍵

`c:=a-b` とおくと `a=b+c`。既存の加法公式
`lubinTateAction_add`(`(b+c)·x=F_f(b·x,c·x)`)と仮定 `a·x=b·x` を
合わせると `b·x=F_f(b·x,c·x)` が出る。ここへ新しい評価レベルの
不等式 `Found/PGC/LubinTateFormalGroupLawEstimate.lean::
norm_aeval_formalGroupLaw_sub_le`(`‖F_f(z,w)-z-w‖≤‖z‖*‖w‖`)を
`(z,w):=(b·x,c·x)` に適用すると

  `‖c·x‖ = ‖F_f(b·x,c·x)-b·x-c·x‖ ≤ ‖b·x‖*‖c·x‖`

`b·x` は捩れ点なので `‖b·x‖<1`(`spectralNorm_lt_one_of_mem_
iteratedLubinTateTorsionPoints`)。`c·x≠0` と仮定すると
`‖c·x‖≤‖b·x‖*‖c·x‖<‖c·x‖` で矛盾——よって `c·x=0`、すなわち
`(a-b)·x=0`。既存の核の特徴づけ(`lubinTateActionAtTorsionPoint_
eq_zero_iff_dvd_of_mem_iteratedLubinTatePsiTorsionPoints`)から
`π^n∣(a-b)` が出る。

`F_f` の減法・逆元・対数のどれも一切構築しない、加法公式+この
評価不等式だけで閉じる経路——単射性に関する当初の3つの候補
(形式逆元・Y-線形係数・Lubin-Tate対数)のいずれとも異なる、より
軽い第4の経路が実際に完成した形。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`𝒪_K` 加群の作用の
単射性**: `x` が原始的な `π^n`-捩れ点(`ψ_n` の根)のとき、
`a·x=b·x` ならば `π^n∣(a-b)`。 -/
theorem lubinTateActionAtTorsionPoint_injective_of_eq
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a b : 𝒪[K.carrier])
    (heq : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a =
        lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b) :
    π ^ n ∣ (a - b) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  rw [← lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_mem_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hmem hxn]
  set c := a - b with hc
  have hac : a = b + c := by rw [hc]; ring
  have hadd : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a =
      MvPowerSeries.aeval (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b c)
        (formalGroupLaw hq hπmax f hf0 hf1 hf) := by
    rw [hac]; exact lubinTateAction_add K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b c
  have hcombine : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b =
      MvPowerSeries.aeval (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b c)
        (formalGroupLaw hq hπmax f hf0 hf1 hf) :=
    heq.symm.trans hadd
  set bx := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b with hbx
  set cx := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem c with hcx
  show cx = 0
  have hbound : ‖MvPowerSeries.aeval (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b c)
      (formalGroupLaw hq hπmax f hf0 hf1 hf) -
      (fun i : Fin 2 => if i = 0 then bx else cx) 0 - (fun i : Fin 2 => if i = 0 then bx else cx) 1‖ ≤
      ‖(fun i : Fin 2 => if i = 0 then bx else cx) 0‖ * ‖(fun i : Fin 2 => if i = 0 then bx else cx) 1‖ :=
    norm_aeval_formalGroupLaw_sub_le hq hπmax hπne0 f hf0 hf1 hf
      (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b c)
      bx.2 cx.2 (algebraMap_mem_adjoinIntegers K x)
  simp only [if_pos, if_neg (by decide : (1 : Fin 2) ≠ 0)] at hbound
  rw [← hcombine, sub_self, zero_sub, norm_neg] at hbound
  have hbxlt1 : ‖bx‖ < 1 := by
    show spectralNorm K.carrier K.closure
      (↑(↑(bx : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) < 1
    exact spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n _
      (lubinTateActionAtTorsionPoint_mem K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem b)
  by_contra hcxne
  have hcxpos : 0 < ‖cx‖ := norm_pos_iff.mpr hcxne
  have : ‖cx‖ < ‖cx‖ := by
    calc ‖cx‖ ≤ ‖bx‖ * ‖cx‖ := hbound
      _ < 1 * ‖cx‖ := by apply mul_lt_mul_of_pos_right hbxlt1 hcxpos
      _ = ‖cx‖ := one_mul _
  exact absurd this (lt_irrefl _)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`unitActionQuotientLift`
の単射性**——上の`lubinTateActionAtTorsionPoint_injective_of_eq`を
`u,v:(𝒪_K)^×`(単数)に適用し、`π^n∣(u-v)`から`u⁻¹*v∈principalUnits
K π n`(`QuotientGroup.eq`+`mem_principalUnits_iff`+`u⁻¹*u=1`だけの
純代数計算)を導くだけ。既存の濃度の一致
(`card_principalUnitsQuotient`=`card_iteratedLubinTatePsiTorsionPoints`、
`QuotientCardinality.lean`)と合わせれば、この単射写像は有限集合
`(𝒪_K)^×⧸principalUnits K π n`から同じ濃度の集合への単射——
すなわち**全単射**であることまであと僅かの距離になった。 -/
theorem unitActionQuotientLift_injective
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Function.Injective (unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem) := by
  intro U V
  induction U using QuotientGroup.induction_on with
  | H u =>
    induction V using QuotientGroup.induction_on with
    | H v =>
      intro heq
      rw [unitActionQuotientLift_mk, unitActionQuotientLift_mk] at heq
      obtain ⟨d, hd⟩ := lubinTateActionAtTorsionPoint_injective_of_eq K hq hπmax hπne0 f hf0 hf1 hf
        n hn x hxψ hxn hmem (u : 𝒪[K.carrier]) (v : 𝒪[K.carrier]) heq
      rw [QuotientGroup.eq, mem_principalUnits_iff]
      have hinv : (↑u⁻¹ * ↑u : 𝒪[K.carrier]) = 1 := by
        exact_mod_cast u.inv_mul
      refine ⟨-(d * (↑u⁻¹ : 𝒪[K.carrier])), ?_⟩
      show (↑u⁻¹ * ↑v : 𝒪[K.carrier]) = 1 + -(d * ↑u⁻¹) * π ^ n
      linear_combination -(↑u⁻¹ : 𝒪[K.carrier]) * hd + hinv

end ABC3.Found.PGC
