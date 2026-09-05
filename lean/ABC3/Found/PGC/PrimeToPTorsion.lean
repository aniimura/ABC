import ABC3.Found.PGC.AbsGalUnitsSurjective

/-!
# `𝒪_K^×` の「`p` と素な捩れ」は `𝓀^×` そのもの——群から `q` が読める

[pGC] Proposition 1.2 は「剰余体の元の個数 `q` と絶対次数 `[K:ℚ_p]` は
`Γ_K` から群論的に回復できる」と主張する。原典の論拠は相互律
`Γ_K^ab ≅ (K^×)^∧` と、そこからの**計数**:

* **捩れの prime-to-p 部分が `q-1` 個**  ← 本ファイル
* pro-p 商の階数が `[K:ℚ_p]+1`  ← `PrincipalUnitsRank.lean`(階数 `[K:ℚ_p]`)
  と `UnitsSplit.lean`(`K^× ≅ ℤ × 𝒪_K^×` の `ℤ` の分)

本ファイルは前者を確立する。すなわち

```
μ^{(p')}(𝒪_K) := {u ∈ 𝒪_K^× | ある p と素な m について u^m = 1}  ≃*  𝓀^×
```

——**群 `𝒪_K^×` の中だけで定義できる部分群**の位数が、ちょうど `q-1`。

## 証明の二本の柱

* **単射**(`eq_of_unitsResidueHom_eq_of_pow_eq_one`): 主単数群 `1+𝔪_K` に
  `p` と素な捩れは無い。超距離での初等的な議論で足りる——
  `u = 1+x`、`p ∤ m` のとき幾何級数 `S = Σ_{i<m} u^i` は
  `‖S - m‖ ≤ ‖x‖ < 1 = ‖m‖` を満たすので `‖S‖ = 1`、したがって
  `‖u^m - 1‖ = ‖S‖·‖x‖ = ‖x‖`(`norm_pow_sub_one_of_not_dvd`)。
  p進対数も二項展開も要らない。
* **全射**: Teichmüller 持ち上げ(`Found/Teichmuller.lean`)が剰余写像の
  群としての切断を与える。`ω(v)^{q-1} = ω(v^{q-1}) = 1` で `p ∤ q-1`。

## `Γ_K` との接続

`AbsGalUnitsSurjective.lean::exists_surjective_absGal_units` により
`𝒪_K^×` は**無条件に** `Γ_K` の商である。したがって `q-1` は
「`Γ_K` のある商の、`p` と素な捩れ部分群の位数」として現れる
——Proposition 1.2 が要求する「群論的な回復」の、計数の側の中身。
残るのは商の**標準性**(相互律そのもの)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 主単数群に `p` と素な捩れが無いこと -/

/-- `p ∤ m` なら `‖(m : K)‖ = 1`(`Padic.norm_natCast_eq_one_iff`)。 -/
theorem norm_natCast_eq_one_of_not_dvd (K : PAdicLocalField p) {m : ℕ} (hm : ¬ (p ∣ m)) :
    ‖((m : ℕ) : K.carrier)‖ = 1 := by
  have h1 : ((m : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((m : ℕ) : ℚ_[p]) := by
    push_cast; ring
  rw [h1, norm_algebraMap K]
  exact Padic.norm_natCast_eq_one_iff.mpr ((Nat.Prime.coprime_iff_not_dvd Fact.out).mpr hm)

/-- **`‖u-1‖ < 1` かつ `p ∤ m` なら `‖u^m - 1‖ = ‖u - 1‖`**。

`u^m - 1 = S·(u-1)`(`S = Σ_{i<m} u^i`)で、`‖S - m‖ ≤ ‖u-1‖ < 1 = ‖m‖`
だから超距離で `‖S‖ = ‖m‖ = 1`。二項展開を経由しない。 -/
theorem norm_pow_sub_one_of_not_dvd (K : PAdicLocalField p) {u : K.carrier}
    (hu : ‖u - 1‖ < 1) {m : ℕ} (hm : ¬ (p ∣ m)) : ‖u ^ m - 1‖ = ‖u - 1‖ := by
  have hnu : ‖u‖ = 1 := norm_eq_one_of_norm_sub_one_lt K hu
  have hpow : ∀ i : ℕ, ‖u ^ i‖ = 1 := by
    intro i; rw [norm_pow, hnu, one_pow]
  have hgeom : ∀ i : ℕ, ‖∑ j ∈ Finset.range i, u ^ j‖ ≤ 1 :=
    fun i => IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg zero_le_one
      (fun j _ => le_of_eq (hpow j))
  have hui : ∀ i : ℕ, ‖u ^ i - 1‖ ≤ ‖u - 1‖ := by
    intro i
    rw [← geom_sum_mul u i, norm_mul]
    calc ‖∑ j ∈ Finset.range i, u ^ j‖ * ‖u - 1‖ ≤ 1 * ‖u - 1‖ :=
          mul_le_mul_of_nonneg_right (hgeom i) (norm_nonneg _)
      _ = ‖u - 1‖ := one_mul _
  set S := ∑ i ∈ Finset.range m, u ^ i with hS
  have hSm : ‖S - ((m : ℕ) : K.carrier)‖ ≤ ‖u - 1‖ := by
    have hexp : S - ((m : ℕ) : K.carrier) = ∑ i ∈ Finset.range m, (u ^ i - 1) := by
      rw [Finset.sum_sub_distrib, hS]; simp
    rw [hexp]
    exact IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (norm_nonneg _) (fun i _ => hui i)
  have hmnorm : ‖((m : ℕ) : K.carrier)‖ = 1 := norm_natCast_eq_one_of_not_dvd K hm
  have hSlt : ‖S - ((m : ℕ) : K.carrier)‖ < 1 := lt_of_le_of_lt hSm hu
  have hSne : ‖S - ((m : ℕ) : K.carrier)‖ ≠ ‖((m : ℕ) : K.carrier)‖ := by
    rw [hmnorm]; exact ne_of_lt hSlt
  have hSnorm : ‖S‖ = 1 := by
    have hdec : S = (S - ((m : ℕ) : K.carrier)) + ((m : ℕ) : K.carrier) := by ring
    rw [hdec, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier) hSne, hmnorm]
    exact max_eq_right (le_of_lt hSlt)
  rw [← geom_sum_mul u m, ← hS, norm_mul, hSnorm, one_mul]

/-- **主単数群には `p` と素な捩れが無い**。 -/
theorem eq_one_of_pow_eq_one_of_not_dvd (K : PAdicLocalField p) {u : K.carrier}
    (hu : ‖u - 1‖ < 1) {m : ℕ} (hm : ¬ (p ∣ m)) (hpow : u ^ m = 1) : u = 1 := by
  have h := norm_pow_sub_one_of_not_dvd K hu hm
  rw [hpow, sub_self, norm_zero] at h
  have hz : u - 1 = 0 := norm_eq_zero.mp h.symm
  linear_combination hz

/-! ## 単数群から剰余体の乗法群へ -/

/-- 剰余体への還元(単数群の準同型)。 -/
noncomputable abbrev unitsResidueHom (K : PAdicLocalField p) :
    (𝒪[K.carrier])ˣ →* (𝓀[K.carrier])ˣ :=
  Units.map (IsLocalRing.residue (𝒪[K.carrier]) : 𝒪[K.carrier] →* 𝓀[K.carrier])

/-- 剰余が 1 になる単数は主単数——`‖w-1‖ < 1`。 -/
theorem norm_sub_one_lt_of_residue_eq_one (K : PAdicLocalField p) {w : (𝒪[K.carrier])ˣ}
    (hw : unitsResidueHom K w = 1) :
    ‖((w : 𝒪[K.carrier]) : K.carrier) - 1‖ < 1 := by
  have h0 : IsLocalRing.residue (𝒪[K.carrier]) ((w : 𝒪[K.carrier]) - 1) = 0 := by
    rw [map_sub, map_one]
    have hv := congrArg (Units.val) hw
    rw [Units.val_one] at hv
    exact (by rw [show IsLocalRing.residue (𝒪[K.carrier]) (w : 𝒪[K.carrier])
      = ((unitsResidueHom K w : (𝓀[K.carrier])ˣ) : 𝓀[K.carrier]) from rfl, hv, sub_self])
  rw [IsLocalRing.residue_eq_zero_iff, IsLocalRing.mem_maximalIdeal, mem_nonunits_iff,
    Valued.integer.isUnit_iff_norm_eq_one] at h0
  have hle : ‖((w : 𝒪[K.carrier]) - 1 : 𝒪[K.carrier])‖ ≤ 1 := by
    have hm := ((w : 𝒪[K.carrier]) - 1).2
    rw [Valued.integer.mem_iff] at hm
    exact hm
  have heq : ‖((w : 𝒪[K.carrier]) - 1 : 𝒪[K.carrier])‖
      = ‖((w : 𝒪[K.carrier]) : K.carrier) - 1‖ := rfl
  rw [heq] at h0 hle
  exact lt_of_le_of_ne hle h0

/-- **`p` と素な冪根の上で還元は単射**——差の商が主単数で捩れを持つから 1。 -/
theorem eq_of_unitsResidueHom_eq_of_pow_eq_one (K : PAdicLocalField p) {m : ℕ} (hm : ¬ (p ∣ m))
    {u v : (𝒪[K.carrier])ˣ} (hu : u ^ m = 1) (hv : v ^ m = 1)
    (h : unitsResidueHom K u = unitsResidueHom K v) : u = v := by
  set w : (𝒪[K.carrier])ˣ := u * v⁻¹ with hw
  have hres : unitsResidueHom K w = 1 := by
    rw [hw, map_mul, map_inv, h, mul_inv_cancel]
  have hwpow : w ^ m = 1 := by
    rw [hw, mul_pow, inv_pow, hu, hv, inv_one, mul_one]
  have hval : (((w : 𝒪[K.carrier]) : K.carrier)) ^ m = 1 := by
    have h1 : ((w : 𝒪[K.carrier])) ^ m = 1 := by
      have hc := congrArg (Units.val) hwpow
      rwa [Units.val_pow_eq_pow_val, Units.val_one] at hc
    have hc := congrArg (fun z : 𝒪[K.carrier] => (z : K.carrier)) h1
    simpa using hc
  have hone := eq_one_of_pow_eq_one_of_not_dvd K (norm_sub_one_lt_of_residue_eq_one K hres) hm hval
  have hw1 : w = 1 := Units.ext (Subtype.ext hone)
  rw [hw, mul_inv_eq_one] at hw1
  exact hw1

/-! ## `p` と素な捩れ部分群 -/

/-- **`p` と素な位数を持つ捩れ部分群** `μ^{(p')}(𝒪_K)`。
`𝒪_K^×` の**群構造だけ**で定義されている(`p` はパラメータ)。 -/
def primeToPTorsion (K : PAdicLocalField p) : Subgroup (𝒪[K.carrier])ˣ where
  carrier := {u : (𝒪[K.carrier])ˣ | ∃ m : ℕ, ¬ (p ∣ m) ∧ u ^ m = 1}
  one_mem' := ⟨1, by simpa using (Nat.Prime.one_lt (Fact.out : p.Prime)).ne', one_pow 1⟩
  mul_mem' := by
    rintro a b ⟨m, hm, ha⟩ ⟨n, hn, hb⟩
    refine ⟨m * n, ?_, ?_⟩
    · intro hdvd
      rcases (Nat.Prime.dvd_mul (Fact.out : p.Prime)).mp hdvd with h | h
      exacts [hm h, hn h]
    · rw [mul_pow, pow_mul, pow_mul', ha, hb, one_pow, one_pow, mul_one]
  inv_mem' := by
    rintro a ⟨m, hm, ha⟩
    exact ⟨m, hm, by rw [inv_pow, ha, inv_one]⟩

/-- `q - 1` は `p` と素(`q = p^f`、`f ≥ 1`)。 -/
theorem not_dvd_residueCard_sub_one (K : PAdicLocalField p) :
    ¬ (p ∣ (Fintype.card 𝓀[K.carrier] - 1)) := by
  have hcard : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]; exact residueCard_eq_pow K
  have h2 : 1 < Fintype.card 𝓀[K.carrier] := Fintype.one_lt_card
  have hf : absoluteInertiaDegree K ≠ 0 := by
    intro h0; rw [h0, pow_zero] at hcard; omega
  intro hdvd
  have h1 : p ∣ Fintype.card 𝓀[K.carrier] := hcard ▸ dvd_pow_self p hf
  have hone : p ∣ 1 := by
    have h3 := Nat.dvd_sub h1 hdvd
    rwa [Nat.sub_sub_self (le_of_lt h2)] at h3
  exact (Nat.Prime.one_lt (Fact.out : p.Prime)).ne' (Nat.dvd_one.mp hone)

/-- **★★★`p` と素な捩れ部分群は剰余体の乗法群そのもの** `μ^{(p')}(𝒪_K) ≃* 𝓀^×`。

単射は「主単数群に `p` と素な捩れが無い」こと、全射は Teichmüller 持ち上げ。 -/
noncomputable def primeToPTorsionEquiv (K : PAdicLocalField p) :
    primeToPTorsion K ≃* (𝓀[K.carrier])ˣ := by
  haveI := henselianLocalRing_carrierIntegers K
  haveI : DecidableEq 𝓀[K.carrier] := Classical.decEq _
  refine MulEquiv.ofBijective ((unitsResidueHom K).comp (primeToPTorsion K).subtype) ⟨?_, ?_⟩
  · rintro ⟨a, m, hm, ham⟩ ⟨b, n, hn, hbn⟩ h
    have hmn : ¬ (p ∣ m * n) := by
      intro hdvd
      rcases (Nat.Prime.dvd_mul (Fact.out : p.Prime)).mp hdvd with h' | h'
      exacts [hm h', hn h']
    have ha' : a ^ (m * n) = 1 := by rw [pow_mul, ham, one_pow]
    have hb' : b ^ (m * n) = 1 := by rw [pow_mul', hbn, one_pow]
    exact Subtype.ext (eq_of_unitsResidueHom_eq_of_pow_eq_one K hmn ha' hb' h)
  · intro v
    have hv : v ^ (Fintype.card 𝓀[K.carrier] - 1) = 1 := by
      rw [← Fintype.card_units]; exact pow_card_eq_one
    exact ⟨⟨teichmuller (𝒪[K.carrier]) v, ⟨Fintype.card 𝓀[K.carrier] - 1,
      not_dvd_residueCard_sub_one K, by rw [← map_pow, hv, map_one]⟩⟩,
      residue_teichmuller (𝒪[K.carrier]) v⟩

/-- **★★`q - 1` は群 `𝒪_K^×` から読める**——`p` と素な捩れ部分群の位数がちょうど `q-1`。
[pGC] Proposition 1.2 の「`q` の回復」の計数の側の中身。 -/
theorem card_primeToPTorsion (K : PAdicLocalField p) :
    Nat.card (primeToPTorsion K) = Nat.card 𝓀[K.carrier] - 1 := by
  haveI : DecidableEq 𝓀[K.carrier] := Classical.decEq _
  rw [Nat.card_congr (primeToPTorsionEquiv K).toEquiv, Nat.card_eq_fintype_card,
    Fintype.card_units, Nat.card_eq_fintype_card]

/-- `μ^{(p')}(𝒪_K)` は巡回群(`𝓀^×` が巡回だから)。 -/
theorem isCyclic_primeToPTorsion (K : PAdicLocalField p) : IsCyclic (primeToPTorsion K) :=
  (MulEquiv.isCyclic (primeToPTorsionEquiv K).symm).mp inferInstance

/-- **★★★★`q-1` は `Γ_K` のある商の群論的な不変量として現れる**——
`Γ_K ↠ 𝒪_K^×`(`exists_surjective_absGal_units`、無条件)と
`|μ^{(p')}(𝒪_K)| = q-1`(`card_primeToPTorsion`)を並べたもの。

これは Proposition 1.2 そのものではない(商の**標準性**——相互律
`Γ_K^ab ≅ (K^×)^∧`——がまだ無いので、`Γ_K ≅ Γ_{K'}` から
`𝒪_K^× ≅ 𝒪_{K'}^×` は従わない)。だが原典が計数に使う量が
実際にこの商の中に居ることは、これで確定した。 -/
theorem residueCard_sub_one_in_absGal_quotient (K : PAdicLocalField p) :
    (∃ φ : (K.closure ≃ₐ[K.carrier] K.closure) →* (𝒪[K.carrier])ˣ, Function.Surjective φ)
      ∧ Nat.card (primeToPTorsion K) = Nat.card 𝓀[K.carrier] - 1 :=
  ⟨exists_surjective_absGal_units K, card_primeToPTorsion K⟩

end ABC3.Found.PGC
