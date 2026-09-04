import ABC3.Found.PGC.AdjoinIntegers

/-!
# `|𝒪_K/π^n𝒪_K| = q^n`・`|(𝒪_K/π^n𝒪_K)^×| = q^n-q^{n-1}`(`sorry` 無し)

`Found/PGC/AdjoinIntegers.lean::principalUnitsQuotientEquiv`
(`(𝒪_K)^×⧸principalUnits K π n ≃* (𝒪_K/π^n𝒪_K)^×`)の**濃度**を、
純粋に環論的な事実として計算する——`F_f`・`Λ_n`・`ψ_n`固有の議論は
最後の`card_principalUnitsQuotient`まで一切不要。

## 証明の筋

`𝒪_K`が(実在のp進局所体の整数環として)完備離散付値環であることから、
`x↦π^n*x`が定める加群の写像の核・像を追跡する古典的な「次数付き商」
の議論——`𝒪_K/π𝒪_K ≅ π^n𝒪_K/π^{n+1}𝒪_K`(乗法による同型、
`gradedPieceMap`)を`n`について積み上げ、`Submodule.card_quotient_
mul_card_quotient`(第三同型定理の濃度版)で`|𝒪_K/π^n𝒪_K|=q^n`
(`card_quotient_span_pi_pow`)へ帰納的に結論する。単数の濃度は、
局所環`R`が`R=R^×⊔maximalIdeal R`(集合として)に分かれることと、
`maximalIdeal(𝒪_K/π^{m+1}𝒪_K)`が`𝒪_K/π^m𝒪_K`に同型であること
(`maxIdealCardMap`、上と全く同じ手法だが`π^n`ではなく`π`を掛ける
写像)を組み合わせて`|(𝒪_K/π^n𝒪_K)^×|=q^n-q^{n-1}`
(`card_units_quotient_span_pi_pow`)を得る。

最後に`principalUnitsQuotientEquiv`と組み合わせ、
`|(𝒪_K)^×⧸principalUnits K π n|=q^n-q^{n-1}`
(`card_principalUnitsQuotient`)——`card_iteratedLubinTatePsiTorsionPoints`
(`|ψ_nの根|=q^n-q^{n-1}`、`Found/PGC/LubinTateDistinguishedSeparable.lean`)
と**ぴったり一致する濃度**——を確立する。`unitActionQuotientLift`が
単射・全射のどちらか一方を示せば、有限集合・同じ濃度の写像は
「単射⟺全射⟺全単射」という一般論だけでもう一方が自動的に従う、
という数え上げの土台がこれで揃った。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## `𝒪_K/π𝒪_K ≅ π^n𝒪_K/π^{n+1}𝒪_K`(次数付き商、乗法による同型) -/

/-- `x↦π^n*x`が誘導する`𝒪_K/π𝒪_K → 𝒪_K/π^{n+1}𝒪_K`の写像
(`Submodule.mapQ`)——`x≡x' mod π`なら`π^n*x≡π^n*x' mod π^{n+1}`。 -/
noncomputable def gradedPieceMap {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (n : ℕ) :
    (𝒪[K.carrier] ⧸ Ideal.span ({π} : Set (𝒪[K.carrier]))) →ₗ[𝒪[K.carrier]]
      (𝒪[K.carrier] ⧸ Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier]))) :=
  Submodule.mapQ _ _ (LinearMap.lsmul (𝒪[K.carrier]) (𝒪[K.carrier]) (π ^ n))
    (by
      intro x hx
      simp only [Submodule.mem_comap, LinearMap.lsmul_apply]
      rw [Ideal.mem_span_singleton] at hx ⊢
      obtain ⟨c, hc⟩ := hx
      rw [hc]
      exact ⟨c, by ring⟩)

/-- `gradedPieceMap`は単射(核が`⊥`)——`𝒪_K`が整域・`π≠0`から
`π^n*x≡π^n*x' mod π^{n+1}`は`x≡x' mod π`を強制する
(`π^n`で約分できる)。 -/
theorem gradedPieceMap_ker {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (hπne0 : π ≠ 0) (n : ℕ) :
    LinearMap.ker (gradedPieceMap K π n) = ⊥ := by
  rw [LinearMap.ker_eq_bot']
  intro x hx
  induction x using Submodule.Quotient.induction_on with
  | _ x =>
  rw [gradedPieceMap, Submodule.mapQ_apply] at hx
  rw [Submodule.Quotient.mk_eq_zero] at hx ⊢
  simp only [LinearMap.lsmul_apply, smul_eq_mul] at hx
  rw [Ideal.mem_span_singleton] at hx ⊢
  obtain ⟨c, hc⟩ := hx
  have hc' : π ^ n * x = π ^ n * (π * c) := by rw [hc]; ring
  have := mul_left_cancel₀ (pow_ne_zero n hπne0) hc'
  exact ⟨c, this⟩

/-- `gradedPieceMap`の像はちょうど`π^n𝒪_K/π^{n+1}𝒪_K`
(`(span{π^n}).map mkQ`)。 -/
theorem gradedPieceMap_range {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (n : ℕ) :
    LinearMap.range (gradedPieceMap K π n) =
      Submodule.map (Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier]))).mkQ
        (Ideal.span ({π ^ n} : Set (𝒪[K.carrier])) : Submodule (𝒪[K.carrier]) (𝒪[K.carrier])) := by
  unfold gradedPieceMap
  rw [Submodule.range_mapQ]
  congr 1
  rw [LinearMap.range_eq_map]
  ext y
  simp only [Submodule.mem_map, LinearMap.lsmul_apply, Submodule.mem_top, true_and]
  constructor
  · rintro ⟨x, hx⟩
    rw [Ideal.mem_span_singleton]
    exact ⟨x, by rw [← hx]; ring⟩
  · intro hy
    rw [Ideal.mem_span_singleton] at hy
    obtain ⟨c, hc⟩ := hy
    exact ⟨c, by rw [hc]; ring⟩

/-- `|π^n𝒪_K/π^{n+1}𝒪_K| = |𝒪_K/π𝒪_K|`——`gradedPieceMap`の核が`⊥`
であることから(`LinearMap.quotKerEquivRange`)線型同型が出る。 -/
theorem card_gradedPieceMap_range {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (hπne0 : π ≠ 0) (n : ℕ) :
    Nat.card (LinearMap.range (gradedPieceMap K π n)) =
      Nat.card (𝒪[K.carrier] ⧸ Ideal.span ({π} : Set (𝒪[K.carrier]))) := by
  have hequiv := LinearMap.quotKerEquivRange (gradedPieceMap K π n)
  rw [gradedPieceMap_ker K π hπne0 n] at hequiv
  have h1 := Nat.card_congr hequiv.toEquiv
  have h2 := Nat.card_congr (Submodule.quotEquivOfEqBot
    (⊥ : Submodule (𝒪[K.carrier]) (𝒪[K.carrier] ⧸ Ideal.span ({π} : Set (𝒪[K.carrier])))) rfl).toEquiv
  rw [← h1, h2]

/-! ## `𝒪_K/π^n𝒪_K`は有限(非自明性`nontrivial_quotient_span_pi_pow`は
`Found/PGC/AdjoinIntegers.lean`に既出) -/

/-- `𝒪_K/π^n𝒪_K`は有限——mathlibの`Valued.integer.finite_quotient_
maximalIdeal_pow_of_finite_residueField`(剰余体が有限な完備離散付値
環の`maximalIdeal^n`による商は有限)を`hπmax`(`maximalIdeal=span{π}`)
+`Ideal.span_singleton_pow`(`span{π}^n=span{π^n}`)で橋渡しするだけ。 -/
theorem finite_quotient_span_pi_pow {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (n : ℕ) :
    Finite (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) := by
  haveI := valuationRing_isDVR K
  have h1 : Finite (𝒪[K.carrier] ⧸ (IsLocalRing.maximalIdeal (𝒪[K.carrier])) ^ n) :=
    Valued.integer.finite_quotient_maximalIdeal_pow_of_finite_residueField (Finite.of_fintype _) n
  rw [hπmax, Ideal.span_singleton_pow] at h1
  exact h1

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★`|𝒪_K/π^n𝒪_K| = q^n` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`|𝒪_K/π^n𝒪_K| = q^n`**——`n`についての
帰納法。`n=0`: `span{1}=⊤`、商は自明(濃度1)。`n→n+1`:
`Submodule.card_quotient_mul_card_quotient`(第三同型定理の濃度版、
`span{π^{n+1}}≤span{π^n}`に適用)と`card_gradedPieceMap_range`
(`|π^n𝒪_K/π^{n+1}𝒪_K|=q`)を組み合わせる。純粋に環論的な事実
——`f`・Lubin-Tate固有の議論は不要。 -/
theorem card_quotient_span_pi_pow {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0) (n : ℕ) :
    Nat.card (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) = (pp ^ ff) ^ n := by
  induction n with
  | zero =>
    simp only [pow_zero]
    rw [show Ideal.span ({(1 : 𝒪[K.carrier])} : Set (𝒪[K.carrier])) = ⊤ from Ideal.span_singleton_one]
    exact Nat.card_eq_one_iff_unique.mpr ⟨inferInstance, ⟨0⟩⟩
  | succ n ih =>
    have hST : (Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier])) : Submodule (𝒪[K.carrier]) (𝒪[K.carrier])) ≤
        (Ideal.span ({π ^ n} : Set (𝒪[K.carrier])) : Submodule (𝒪[K.carrier]) (𝒪[K.carrier])) := by
      rw [Ideal.span_le]
      intro x hx
      simp only [Set.mem_singleton_iff] at hx
      rw [hx]
      exact Ideal.mem_span_singleton.mpr (pow_dvd_pow π (Nat.le_succ n))
    have hkey := Submodule.card_quotient_mul_card_quotient
      (Ideal.span ({π ^ n} : Set (𝒪[K.carrier])) : Submodule (𝒪[K.carrier]) (𝒪[K.carrier]))
      (Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier])) : Submodule (𝒪[K.carrier]) (𝒪[K.carrier])) hST
    rw [← gradedPieceMap_range K π n, card_gradedPieceMap_range K π hπne0 n] at hkey
    show Nat.card (𝒪[K.carrier] ⧸ Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier]))) = (pp ^ ff) ^ (n + 1)
    rw [← hkey, ih]
    have hcard : Nat.card (𝒪[K.carrier] ⧸ Ideal.span ({π} : Set (𝒪[K.carrier]))) = pp ^ ff := by
      rw [← hπmax]
      rw [show (𝒪[K.carrier] ⧸ IsLocalRing.maximalIdeal (𝒪[K.carrier])) =
        IsLocalRing.ResidueField (𝒪[K.carrier]) from rfl]
      rw [Nat.card_eq_fintype_card, hq]
    rw [hcard]
    ring

/-! ## `π𝒪_K/π^{n+1}𝒪_K ≅ 𝒪_K/π^n𝒪_K`(最大イデアルの次数付き商) -/

/-- `x↦π*x`が誘導する`𝒪_K/π^n𝒪_K → 𝒪_K/π^{n+1}𝒪_K`の写像——
`gradedPieceMap`と同じ構成だが`π^n`ではなく`π`を掛ける。この写像の
**像**が`𝒪_K/π^{n+1}𝒪_K`の**最大イデアル**になる(`card_units_
quotient_span_pi_pow_succ`で使う)。 -/
noncomputable def maxIdealCardMap {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (n : ℕ) :
    (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) →ₗ[𝒪[K.carrier]]
      (𝒪[K.carrier] ⧸ Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier]))) :=
  Submodule.mapQ _ _ (LinearMap.lsmul (𝒪[K.carrier]) (𝒪[K.carrier]) π)
    (by
      intro x hx
      simp only [Submodule.mem_comap, LinearMap.lsmul_apply]
      rw [Ideal.mem_span_singleton] at hx ⊢
      obtain ⟨c, hc⟩ := hx
      rw [hc]
      exact ⟨c, by ring⟩)

/-- `maxIdealCardMap`は単射(核が`⊥`)——`gradedPieceMap_ker`と同じ
議論、`π`で約分するだけ。 -/
theorem maxIdealCardMap_ker {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (hπne0 : π ≠ 0) (n : ℕ) :
    LinearMap.ker (maxIdealCardMap K π n) = ⊥ := by
  rw [LinearMap.ker_eq_bot']
  intro x hx
  induction x using Submodule.Quotient.induction_on with
  | _ x =>
  rw [maxIdealCardMap, Submodule.mapQ_apply] at hx
  rw [Submodule.Quotient.mk_eq_zero] at hx ⊢
  simp only [LinearMap.lsmul_apply, smul_eq_mul] at hx
  rw [Ideal.mem_span_singleton] at hx ⊢
  obtain ⟨c, hc⟩ := hx
  have hc' : π * x = π * (π ^ n * c) := by rw [hc]; ring
  have := mul_left_cancel₀ hπne0 hc'
  exact ⟨c, this⟩

/-- `maxIdealCardMap`の像はちょうど`π𝒪_K/π^{n+1}𝒪_K`
(`(span{π}).map mkQ`)。 -/
theorem maxIdealCardMap_range {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (n : ℕ) :
    LinearMap.range (maxIdealCardMap K π n) =
      Submodule.map (Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier]))).mkQ
        (Ideal.span ({π} : Set (𝒪[K.carrier])) : Submodule (𝒪[K.carrier]) (𝒪[K.carrier])) := by
  unfold maxIdealCardMap
  rw [Submodule.range_mapQ]
  congr 1
  rw [LinearMap.range_eq_map]
  ext y
  simp only [Submodule.mem_map, LinearMap.lsmul_apply, Submodule.mem_top, true_and]
  constructor
  · rintro ⟨x, hx⟩
    rw [Ideal.mem_span_singleton]
    exact ⟨x, by rw [← hx]; ring⟩
  · intro hy
    rw [Ideal.mem_span_singleton] at hy
    obtain ⟨c, hc⟩ := hy
    exact ⟨c, by rw [hc]; ring⟩

/-- `|π𝒪_K/π^{n+1}𝒪_K| = |𝒪_K/π^n𝒪_K|`——`maxIdealCardMap`の核が
`⊥`であることから。 -/
theorem card_maxIdealCardMap_range {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (hπne0 : π ≠ 0) (n : ℕ) :
    Nat.card (LinearMap.range (maxIdealCardMap K π n)) =
      Nat.card (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) := by
  have hequiv := LinearMap.quotKerEquivRange (maxIdealCardMap K π n)
  rw [maxIdealCardMap_ker K π hπne0 n] at hequiv
  have h1 := Nat.card_congr hequiv.toEquiv
  have h2 := Nat.card_congr (Submodule.quotEquivOfEqBot
    (⊥ : Submodule (𝒪[K.carrier]) (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))) rfl).toEquiv
  rw [← h1, h2]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★`|(𝒪_K/π^n𝒪_K)^×| = q^n-q^{n-1}` -/

/-- `|(𝒪_K/π^{m+1}𝒪_K)^×| = q^{m+1}-q^m`——局所環`S:=𝒪_K/π^{m+1}𝒪_K`
は`S=S^×⊔maximalIdeal S`(集合として、`Equiv.sumCompl IsUnit`)に
分かれ、`maximalIdeal S`は元の局所環の最大イデアルの像に一致する
(`IsLocalRing.map_maximalIdeal_of_surjective`)ことと`maxIdealCardMap`
(`|maximalIdeal S|=|𝒪_K/π^m𝒪_K|=q^m`)を組み合わせる。`m+1`の形で
述べるのは、`n-1+1=n`が`n`が変数のときは定義的に成り立たず(自然数
減法)、`LinearMap.range (maxIdealCardMap K π (n-1))`の周辺の型が
`𝒪_K/π^n𝒪_K`と構文的に一致しないことで`isDefEq`が極端に遅くなる
罠を避けるため——`m+1`の形なら最初から型が一致する。 -/
theorem card_units_quotient_span_pi_pow_succ {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0) (m : ℕ) :
    Nat.card (𝒪[K.carrier] ⧸ Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier])))ˣ =
      (pp ^ ff) ^ (m + 1) - (pp ^ ff) ^ m := by
  set S := 𝒪[K.carrier] ⧸ Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier])) with hS_def
  show Nat.card Sˣ = (pp ^ ff) ^ (m + 1) - (pp ^ ff) ^ m
  haveI hfin : Finite S := finite_quotient_span_pi_pow K hπmax (m + 1)
  haveI hNT : Nontrivial S := nontrivial_quotient_span_pi_pow K hπmax (m + 1) (by omega)
  haveI hloc : IsLocalRing S := IsLocalRing.of_surjective' (Ideal.Quotient.mk _) Ideal.Quotient.mk_surjective
  have hpart : Nat.card S = Nat.card {x : S // IsUnit x} + Nat.card {x : S // ¬ IsUnit x} := by
    rw [← Nat.card_sum]
    exact Nat.card_congr (Equiv.sumCompl IsUnit).symm
  have hunits : Nat.card {x : S // IsUnit x} = Nat.card Sˣ :=
    (Nat.card_congr (Submonoid.unitsTypeEquivIsUnitSubmonoid (M := S)).toEquiv).symm
  have hmax : (IsLocalRing.maximalIdeal (𝒪[K.carrier])).map (Ideal.Quotient.mk
      (Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier])))) = IsLocalRing.maximalIdeal S :=
    IsLocalRing.map_maximalIdeal_of_surjective _ Ideal.Quotient.mk_surjective
  rw [hπmax] at hmax
  have hbridge : ∀ x : S, x ∈ Submodule.map (Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier]))).mkQ
        (Ideal.span ({π} : Set (𝒪[K.carrier])) : Submodule (𝒪[K.carrier]) (𝒪[K.carrier])) ↔
      x ∈ Ideal.map (Ideal.Quotient.mk (Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier]))))
        (Ideal.span ({π} : Set (𝒪[K.carrier]))) := by
    intro x
    rw [Submodule.mem_map, Ideal.mem_map_iff_of_surjective _ Ideal.Quotient.mk_surjective]
    rfl
  have hnu_iff : ∀ x : S, ¬ IsUnit x ↔
      x ∈ LinearMap.range (maxIdealCardMap K π m) := by
    intro x
    rw [maxIdealCardMap_range, hbridge, hmax]
    exact (IsLocalRing.mem_maximalIdeal x).symm
  have hnonunits_equiv : {x : S // ¬ IsUnit x} ≃ LinearMap.range (maxIdealCardMap K π m) :=
    Equiv.subtypeEquivRight hnu_iff
  have hcard_nonunit : Nat.card {x : S // ¬ IsUnit x} = (pp ^ ff) ^ m := by
    rw [Nat.card_congr hnonunits_equiv]
    rw [card_maxIdealCardMap_range K π hπne0 m]
    exact card_quotient_span_pi_pow K hq hπmax hπne0 m
  have hcard_S : Nat.card S = (pp ^ ff) ^ (m + 1) := card_quotient_span_pi_pow K hq hπmax hπne0 (m + 1)
  rw [hcard_S, hcard_nonunit, hunits] at hpart
  omega

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`|(𝒪_K/π^n𝒪_K)^×| = q^n-q^{n-1}`**
(`n≥1`)——`card_units_quotient_span_pi_pow_succ`の`n,hn:1≤n`版への
言い換え(`n=m+1`と書き直すだけ)。 -/
theorem card_units_quotient_span_pi_pow {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0) (n : ℕ) (hn : 1 ≤ n) :
    Nat.card (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ =
      (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le hn
  rw [Nat.add_comm 1 m, Nat.add_sub_cancel]
  exact card_units_quotient_span_pi_pow_succ K hq hπmax hπne0 m

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**目標達成——`|(𝒪_K)^×⧸principalUnits
K π n| = q^n-q^{n-1}`**。`principalUnitsQuotientEquiv`
(`(𝒪_K)^×⧸principalUnits K π n≃*(𝒪_K/π^n𝒪_K)^×`)と`card_units_
quotient_span_pi_pow`を組み合わせるだけ。これは
`card_iteratedLubinTatePsiTorsionPoints`(`Found/PGC/
LubinTateDistinguishedSeparable.lean`、`|ψ_nの根|=q^n-q^{n-1}`)と
**ぴったり一致する濃度**——`unitActionQuotientLift`(定義域は
`(𝒪_K)^×⧸principalUnits K π n`、値域は`ψ_nの根`を含む`adjoinIntegers
K x`)が単射・全射のどちらか一方さえ示せば、有限集合・同じ濃度の
写像は「単射⟺全射⟺全単射」という一般論だけでもう一方が自動的に
従う——`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`への数え上げの土台が完成した。 -/
theorem card_principalUnitsQuotient {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0) (n : ℕ) (hn : 1 ≤ n) :
    Nat.card ((𝒪[K.carrier])ˣ ⧸ principalUnits K π n) = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  rw [Nat.card_congr (principalUnitsQuotientEquiv K hπmax n hn).toEquiv]
  exact card_units_quotient_span_pi_pow K hq hπmax hπne0 n hn

end ABC3.Found.PGC
