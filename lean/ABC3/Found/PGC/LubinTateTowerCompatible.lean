import ABC3.Found.PGC.AdjoinIntegers

/-!
# 捩れ塔が compatible system をなす(`sorry` 無し)

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)へ向けた最初の一歩。古典的な
Lubin-Tate 理論で `L_n:=K(Λ_n)` の塔が意味を持つのは、`Λ_n` の
「原始的な」生成元(`ψ_n` の根)の列が `[π]_f` で互いに繋がっている
——`x∈ψ_{n+1}` の根なら `π·x∈ψ_n` の根——という compatible system の
事実による。本ファイルはこれを2段階で確立する:

1. `lubinTateActionAtTorsionPoint_pi_mem_pred`: `x∈Λ_{n+1}` ならば
   `π·x∈Λ_n`(1段下がる)。
2. `lubinTateActionAtTorsionPoint_pi_mem_iteratedLubinTatePsiTorsionPoints`:
   `x` が「原始的な」π^{n+1}-捩れ点(`ψ_{n+1}` の根)ならば、`π·x` も
   「原始的な」π^n-捩れ点(`ψ_n` の根)——原始性(生成元であること)も
   1段下がる。

どちらも新しい数学的内容は無く、`AdjoinIntegers.lean` で確立済みの
部品(`lubinTateAction_mul`(乗法性)・`pi_pow_action_eq_zero`/
`eq_zero_of_pi_pow_action_eq_zero`(`π^n·x=0⟺D_n(x)=0`の橋渡し)・
`lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_
iteratedLubinTatePsiTorsionPoints`(原始性の特徴づけ)・
`iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints`
(`Λ_n\Λ_{n-1}=ψ_nの根`))の組み立てのみで閉じる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- ★★★★★★★★★★**`Λ_{n+1}` の元に `π` を作用させると `Λ_n` に落ちる**
(捩れ塔が compatible system をなすことの核心)。`x∈Λ_{n+1}` のとき
`z:=π·x`(`x` 自身の座標系での値)は、`π^n·z=π^{n+1}·x=0`
(`lubinTateAction_mul`の乗法性+`pi_pow_action_eq_zero`)から
`D_n(z)=0`(`eq_zero_of_pi_pow_action_eq_zero`)を満たし、これを
`K.closure` へ押し出すと `z`(の`K.closure`での値)が`Λ_n`の元
であることが従う。`lubinTateActionAtTorsionPoint_mem`(同じ段に
留まる事実)と同じ証明の骨格を、段をひとつ下げる形に流用しただけ。 -/
theorem lubinTateActionAtTorsionPoint_pi_mem_pred {p : ℕ} [Fact p.Prime]
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
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  set z := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π with hz_def
  have hz0 : lubinTateEvalAtPoint K x z
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π)
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)) = 0 := by
    rw [lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem (π ^ n) π, ← pow_succ]
    exact pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem
  have hDn0 : Polynomial.aeval z (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0 :=
    eq_zero_of_pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x z _ hz0
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
  refine ⟨(isDistinguishedAt_iteratedLubinTateDistinguished
    hq hπmax hπne0 f hf0 hf1 hf n).monic.map _ |>.ne_zero, ?_⟩
  set g : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hg_def
  have key := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)
    (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g z
  rw [← Polynomial.aeval_def] at key
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro y; rfl
  rw [hgcomp, hDn0, map_zero] at key
  show (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).eval (g z) = 0
  rw [Polynomial.eval_map]
  exact key.symm

/-- ★★★★★★★★★★★★★★★★**原始性も1段下がる**: `x` が「原始的な」
π^{n+1}-捩れ点(`ψ_{n+1}` の根)ならば、`π·x` も「原始的な」π^n-捩れ点
(`ψ_n` の根)。`lubinTateActionAtTorsionPoint_pi_mem_pred` で `π·x∈Λ_n`
は分かっているので、`Λ_n\Λ_{n-1}=ψ_nの根`
(`iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints`)
により `π·x∉Λ_{n-1}` だけ示せばよい。背理法: `π·x∈Λ_{n-1}` と仮定すると
`D_{n-1}(π·x)=0`(`K.closure` レベル)を `x` 自身の座標系へ引き戻して
`π^{n-1}·(π·x)=0`、すなわち乗法性で `π^n·x=0` が出るが、これは `x` が
`ψ_{n+1}` の根であること(`lubinTateActionAtTorsionPoint_pi_pow_pred_
ne_zero_of_mem_iteratedLubinTatePsiTorsionPoints`、`(n+1)-1=n` の場合)
に矛盾する。 -/
theorem lubinTateActionAtTorsionPoint_pi_mem_iteratedLubinTatePsiTorsionPoints
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
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n + 1) (by omega))
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  have hyn : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n :=
    lubinTateActionAtTorsionPoint_pi_mem_pred K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
  rw [← iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn, Finset.mem_sdiff]
  refine ⟨hyn, ?_⟩
  intro hyn1
  set z := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π with hz_def
  set g : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hg_def
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro y; rfl
  have hyn1' : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1))).eval (g z) = 0 := by
    rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hyn1
    exact hyn1.2
  rw [Polynomial.eval_map] at hyn1'
  have key := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1))
    (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g z
  rw [← Polynomial.aeval_def, hgcomp] at key
  have hDn1z : Polynomial.aeval z (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1)) = 0 := by
    have hginj : Function.Injective g := fun a b h => Subtype.ext (Subtype.ext h)
    apply hginj
    rw [key, hyn1', map_zero]
  have hz : PowerSeries.HasEval z :=
    hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π
  have hπnm1z : lubinTateEvalAtPoint K x z hz (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (n - 1))) = 0 := by
    show PowerSeries.aeval hz (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (n - 1))) = 0
    rw [LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf (n - 1),
      iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf (n - 1),
      map_mul, PowerSeries.aeval_coe, hDn1z, zero_mul]
  have hπnx : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem (π ^ n) = 0 := by
    rw [show π ^ n = π ^ (n - 1) * π by rw [← pow_succ]; congr 1; omega]
    rw [← lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem (π ^ (n - 1)) π]
    exact hπnm1z
  have hne := lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf (n + 1) (by omega) x hxψ hmem hx
  rw [show n + 1 - 1 = n by omega] at hne
  exact hne hπnx

end ABC3.Found.PGC
