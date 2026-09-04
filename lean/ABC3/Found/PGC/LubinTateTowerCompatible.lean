import ABC3.Found.PGC.AdjoinIntegers
import ABC3.Found.PGC.LubinTateActionInclusion
import ABC3.Found.PGC.LubinTateActionEquivariance

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
3. `adjoin_pi_mem_pred_le`: 体の言葉での言い換え——
   `K.carrier⟮π·x⟯≤K.carrier⟮x⟯`。すなわち塔 `L_n:=K.carrier⟮π·x⟯≤
   L_{n+1}:=K.carrier⟮x⟯` が実際に体の包含として成り立つ。
4. `principalUnits_succ_le`: `principalUnits K π (n+1)≤principalUnits
   K π n`——単位群側の「モジュラス」の単調性、純環論的事実。

1〜4はいずれも新しい数学的内容は無く、`AdjoinIntegers.lean` で
確立済みの部品(`lubinTateAction_mul`(乗法性)・`pi_pow_action_eq_zero`/
`eq_zero_of_pi_pow_action_eq_zero`(`π^n·x=0⟺D_n(x)=0`の橋渡し)・
`lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_
iteratedLubinTatePsiTorsionPoints`(原始性の特徴づけ)・
`iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints`
(`Λ_n\Λ_{n-1}=ψ_nの根`))の組み立てのみで閉じる。

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★節目——`reciprocityMap`の`n`跨ぎの可換性(cross-point bridging を突破)

続く2定理は`LubinTateActionInclusion.lean`(体の包含に沿った冪級数
評価のcross-point bridging、新規)を使って、当初「壁」として記録
されていた障害を実際に突破する:

5. `lubinTateActionAtTorsionPoint_pi_mul_eq_pred`: `(a*π)·x`(`x`の
   座標系、level `n+1`)と`a·y`(`y:=π·x`自身の座標系、level `n`)が
   同じ`K.closure`の値を与える——`lubinTateEvalAtPoint_inclusion_comm`
   (cross-point bridging)と`lubinTateAction_mul`(乗法性)の組み合わせ。
6. `reciprocityMap_pred_eq_map_succ`: **`reciprocityMap`の`n`跨ぎの
   可換性**——同じ大域的`σ`について、level `n+1`での`reciprocityMap`
   を`principalUnits K π n`まで落とした(`QuotientGroup.map`)ものが
   level `n`での`reciprocityMap`に一致する。5・Galois同変性
   (`algEquiv_lubinTateActionAtTorsionPoint_comm`)・一意性
   (`existsUnique_unitActionQuotient_eq_algEquiv`)を組み合わせ、
   「`u_σ·y=u_σ·(π·x)=(u_σ*π)·x=π·(u_σ·x)=π·σ(x)=σ(π·x)=σ(y)`」
   という等式の連鎖(乗法性2回+Galois同変性)を実際に完成させた。

これで節目(5)(射影極限`Gal(L_π/K)≅𝒪_K^×`)の核心的な数学的内容
——各段の`reciprocityMap`(ひいては`galoisReciprocityEquiv`)が
射影系をなすこと——が確立された。残るのは`Mathlib.FieldTheory.
Galois.Infinite`を使って実際に逆極限を組み立てる、というAPI調査・
組み立て工程。
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

/-- ★★★★★★★★★★**`L_n:=K(π·x) ≤ L_{n+1}:=K(x)`**——塔の体としての
包含。`π·x`(`x` 自身の座標系での値、`lubinTateActionAtTorsionPoint`)
は定義から常に`adjoinIntegers K x ⊆ K.carrier⟮x⟯`の元なので、
`IntermediateField.adjoin_simple_le_iff`(生成元が属せば adjoin は
包含される)を適用するだけ。`algEquiv_map_adjoin_eq`
(`LubinTateActionEquivariance.lean`)で使った「`SetLike.coe_mem`+
`adjoin_simple_le_iff`」と全く同じ型の議論。 -/
theorem adjoin_pi_mem_pred_le {p : ℕ} [Fact p.Prime]
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
    IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  apply IntermediateField.adjoin_simple_le_iff.mpr
  exact SetLike.coe_mem _

/-- **`principalUnits K π (n+1) ≤ principalUnits K π n`**——単位群の
「モジュラス」`π^n`が`n`について単調に細かくなること。`v-1∈span{π^{n+1}}
⟹ v-1∈span{π^n}`を`mem_principalUnits_iff`で`v=1+cπ^{n+1}=1+(cπ)π^n`
と書き直すだけの、純粋に環論的な事実——`Gal(L_π/K)≅𝒪_K^×`への射影
極限で、`(𝒪_K)^×⧸principalUnits K π (n+1)→(𝒪_K)^×⧸principalUnits K π n`
という自然な射影(`QuotientGroup.map`)を定義するために要る。 -/
theorem principalUnits_succ_le {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (n : ℕ) :
    principalUnits K π (n + 1) ≤ principalUnits K π n := by
  intro v hv
  rw [mem_principalUnits_iff] at hv ⊢
  obtain ⟨c, hc⟩ := hv
  exact ⟨c * π, by rw [hc, pow_succ]; ring⟩

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**cross-point bridging の初めての実戦
投入**: `(a*π)·x`(`x`の座標系、level `n+1`)と`a·y`(`y:=π·x`自身の
座標系、level `n`)は同じ`K.closure`の値。`w:=⟨⟨y,hymem⟩,_⟩:
adjoinIntegers K y`(`y`自身の座標系での基準点)と`z:=π·x:adjoinIntegers
K x`(`x`の座標系での同じ点)が、`adjoinIntegersInclusionAlgHom`を
経由すると一致すること(`coe_adjoinIntegersInclusion`)を`hwz`として
確立し、`lubinTateEvalAtPoint_inclusion_comm`(cross-point bridging)
・`lubinTateEvalAtPoint_congr`(証明項の違いの吸収)・`lubinTateAction_
mul`(乗法性)を繋ぐだけ。 -/
theorem lubinTateActionAtTorsionPoint_pi_mul_eq_pred {p : ℕ} [Fact p.Prime]
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
    (hx1 : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (hyn : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hymem : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)} : Set K.closure))]
    (a : 𝒪[K.carrier]) :
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem (a * π)) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) =
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n
        (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) hyn hymem a) :
        IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)} : Set K.closure)) :
        K.closure) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  set y : K.closure := (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) with hy_def
  haveI := completeSpace_adjoinIntegers K y
  haveI := isLinearTopology_adjoinIntegers K y
  haveI := continuousSMul_adjoinIntegers K y
  set hle := adjoin_pi_mem_pred_le K hq hπmax hπne0 f hf0 hf1 hf n x hx1 hmem with hle_def
  set z : adjoinIntegers K x := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π
    with hz_def
  set w : adjoinIntegers K y := ⟨⟨y, hymem⟩,
      mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n y hyn hymem⟩
    with hw_def
  have hwz : adjoinIntegersInclusionAlgHom K x y hle w = z := by
    apply Subtype.ext; apply Subtype.ext
    show (↑(↑((adjoinIntegersInclusion K x y hle) w) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = _
    rw [coe_adjoinIntegersInclusion]
  have hw : PowerSeries.HasEval w :=
    hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n y hyn hymem
  have hw' : PowerSeries.HasEval (adjoinIntegersInclusionAlgHom K x y hle w) := by
    rw [hwz]; exact hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π
  have hcomm := lubinTateEvalAtPoint_inclusion_comm K x y hle w hw hw'
    (LubinTateAction hq hπmax f hf0 hf1 hf a)
  rw [lubinTateEvalAtPoint_congr K x hw' hwz (LubinTateAction hq hπmax f hf0 hf1 hf a)] at hcomm
  have hcanon : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n y hyn hymem a =
      lubinTateEvalAtPoint K y w hw (LubinTateAction hq hπmax f hf0 hf1 hf a) := rfl
  rw [hcanon, ← hcomm]
  have hmul := lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem a π
  exact congrArg (fun v : adjoinIntegers K x =>
    (↑(↑v : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) hmul.symm

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMap`の`n`跨ぎの可換性**: 同じ大域的`σ`について、
level `n+1`での`reciprocityMap`を`principalUnits K π n`まで落とした
(`QuotientGroup.map`)ものが level `n`での`reciprocityMap`に一致する。
`u_σ`(level `n+1`での`reciprocityMap σ`の代表元)を取り、
`u_σ·y=u_σ·(π·x)=(u_σ*π)·x`(`lubinTateActionAtTorsionPoint_pi_mul_
eq_pred`)`=π·(u_σ·x)`(`lubinTateAction_mul`)`=π·σ(x)`
(`reciprocityMap_spec`)`=σ(π·x)=σ(y)`(Galois同変性
`algEquiv_lubinTateActionAtTorsionPoint_comm`)という等式の連鎖から
`u_σ·y=σ(y)`を導き、一意性(`existsUnique_unitActionQuotient_eq_
algEquiv`)で`reciprocityMap`(level `n`)`=mk(u_σ)=QuotientGroup.map
(mk(u_σ)at level n+1)`を結論する。 -/
theorem reciprocityMap_pred_eq_map_succ {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ1 : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n + 1) (by omega))
    (hx1 : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (hyψ : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hyn : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hymem : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)} : Set K.closure))]
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    QuotientGroup.map (principalUnits K π (n + 1)) (principalUnits K π n) (MonoidHom.id _)
        (principalUnits_succ_le K π n)
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (n + 1) (by omega) x hxψ1 hx1 hmem σ) =
      reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn
        (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) hyψ hyn hymem σ := by
  set y : K.closure := (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) with hy_def
  have hn1 : 1 ≤ n + 1 := by omega
  obtain ⟨uσ, huσ⟩ := QuotientGroup.mk_surjective
    (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1 x hxψ1 hx1 hmem σ)
  have hspecσ : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem
        (uσ : 𝒪[K.carrier])) : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = σ x := by
    have h := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1 x hxψ1 hx1 hmem σ
    rw [← huσ, unitActionQuotientLift_mk] at h
    exact h
  have hw : adjoinIntegersRestrictSelfAlgHom K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1 x hxψ1 hx1 hmem σ
      (⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
        K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem⟩ : adjoinIntegers K x) =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem (uσ : 𝒪[K.carrier]) := by
    show adjoinIntegersRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1 x hxψ1 hx1 hmem σ _ = _
    apply Subtype.ext; apply Subtype.ext
    rw [coe_adjoinIntegersRestrictSelf, hspecσ]
  have hval : lubinTateActionAtAlgEquivPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1 x hxψ1 hx1 hmem σ π =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem
        (π * (uσ : 𝒪[K.carrier])) := by
    show lubinTateEvalAtPoint K x _
      (hasEval_adjoinIntegersRestrictSelfAlgHom_mk K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1 x hxψ1 hx1 hmem σ)
      (LubinTateAction hq hπmax f hf0 hf1 hf π) = _
    rw [lubinTateEvalAtPoint_congr K x
      (hasEval_adjoinIntegersRestrictSelfAlgHom_mk K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1 x hxψ1 hx1 hmem σ)
      hw (LubinTateAction hq hπmax f hf0 hf1 hf π)]
    exact lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf (n + 1) x hx1 hmem π (uσ : 𝒪[K.carrier])
  have hcomm := algEquiv_lubinTateActionAtTorsionPoint_comm K hq hπmax hπne0 f hf0 hf1 hf (n + 1) hn1
    x hxψ1 hx1 hmem σ π
  rw [hval, mul_comm] at hcomm
  have hbridge := lubinTateActionAtTorsionPoint_pi_mul_eq_pred K hq hπmax hπne0 f hf0 hf1 hf n x hx1 hmem
    hyn hymem (uσ : 𝒪[K.carrier])
  rw [hbridge] at hcomm
  have hspecy := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn y hyψ hyn hymem σ
  have hspecmk : (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n y hyn hymem
      (QuotientGroup.mk (s := principalUnits K π n) uσ)) :
      IntermediateField.adjoin K.carrier ({y} : Set K.closure)) : K.closure) = σ y := by
    rw [unitActionQuotientLift_mk]; exact hcomm.symm
  have hfinal := (existsUnique_unitActionQuotient_eq_algEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn
    y hyψ hyn hymem σ).unique hspecy hspecmk
  rw [hfinal, ← huσ, QuotientGroup.map_mk]
  rfl

end ABC3.Found.PGC
