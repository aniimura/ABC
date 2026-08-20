/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.ProL.LPart
import ABC3.Found.ProL.FinitePrimary

/-!
# `M ≅ ∏_l M[l]` —— 可換副有限群の pro-`l` 分解(チェーン `prol` の葉 `decomposition`)

原文 (FrdI p.52):
> Thus, M decomposes as a direct product of pro-l groups M [l], where l varies over

## ★★★測って選んだ道(2026-08-19)

★**逆向き `Φ : ∏_l M[l] → M` を先に作る**。
`Φ y` は「各有限商で `∏_{p ∣ |M/U|} y_p` を取る」整合族の一意な持ち上げである
(`existsUnique_mk_eq_of_compatible`)。

★★**連続性は圏論を経由しない** —— 副有限群では
「開正規部分群の剰余類」が開基をなす(`exist_openNormalSubgroup_sub_open_nhds_of_one`)ので、
`y ↦ f y U`(**有限**積、離散空間への写像)が連続であることだけで足りる。

★★★全単射を示したら、**コンパクト→Hausdorff の連続全単射は同相**
(`Continuous.homeoOfEquivCompactToT2`)で終わる。

| 段 | 中身 |
|---|---|
| `lComp` | `x` の `l`-成分(各商で `primProj`、整合性は `primProj_naturality`) |
| `Φ` | `∏_l M[l] → M` の構成と群準同型性 |
| 単射 | `primProj` を当てて成分ごとに分離 |
| 全射 | `lComp` で作った族が `x` に戻る(`prod_primProj`) |
| 連続 | 剰余類が開基／`f y U` が有限積 |
-/

namespace ABC3.Found.ProL

open CategoryTheory

universe u

/-- ★`ℕ` から `Nat.Primes` への全域写像(素数でなければ `2`)。 -/
noncomputable def toPrimes (p : ℕ) : Nat.Primes :=
  if h : p.Prime then ⟨p, h⟩ else ⟨2, Nat.prime_two⟩

@[simp] theorem toPrimes_val {p : ℕ} (h : p.Prime) : (toPrimes p).1 = p := by
  simp [toPrimes, h]

section CommProfinite

variable {M : Type u} [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
  [CompactSpace M] [TotallyDisconnectedSpace M]

/-- ★`M` を `ProfiniteGrp` として見る。 -/
abbrev asProfiniteGrp (M : Type u) [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
    [CompactSpace M] [TotallyDisconnectedSpace M] : ProfiniteGrp.{u} := ProfiniteGrp.of M

/-! ## ★1. `l`-成分 -/

variable (M) in
/-- ★`x` の `l`-成分を与える整合族。 -/
noncomputable def lCompFam (l : Nat.Primes) (x : M) (U : OpenNormalSubgroup M) :
    M ⧸ U.toSubgroup :=
  primProj (M ⧸ U.toSubgroup) l.1 (QuotientGroup.mk x)

variable (M) in
theorem lCompFam_compatible (l : Nat.Primes) (x : M) :
    ∀ (U V : OpenNormalSubgroup M) (h : U.toSubgroup ≤ V.toSubgroup),
      QuotientGroup.map U.toSubgroup V.toSubgroup (MonoidHom.id M) h (lCompFam M l x U)
        = lCompFam M l x V := by
  intro U V h
  show QuotientGroup.map _ _ _ h (primProj _ l.1 (QuotientGroup.mk x))
    = primProj _ l.1 (QuotientGroup.mk x)
  rw [primProj_naturality _ l.2]
  congr 1

variable (M) in
/-- ★★**`x` の `l`-成分**。 -/
noncomputable def lComp (l : Nat.Primes) (x : M) : M :=
  (existsUnique_mk_eq_of_compatible (asProfiniteGrp M) (lCompFam M l x)
    (lCompFam_compatible M l x)).choose

variable (M) in
theorem lComp_spec (l : Nat.Primes) (x : M) (U : OpenNormalSubgroup M) :
    (QuotientGroup.mk (lComp M l x) : M ⧸ U.toSubgroup)
      = primProj (M ⧸ U.toSubgroup) l.1 (QuotientGroup.mk x) :=
  (existsUnique_mk_eq_of_compatible (asProfiniteGrp M) (lCompFam M l x)
    (lCompFam_compatible M l x)).choose_spec.1 U

variable (M) in
/-- ★`l`-成分は `M[l]` に入る。 -/
theorem lComp_mem (l : Nat.Primes) (x : M) : lComp M l x ∈ lPart M l.1 := by
  intro U
  have h := lComp_spec M l x U
  obtain ⟨k, hk⟩ := primProj_mem_primaryComponent (A := M ⧸ U.toSubgroup) l.2
    (QuotientGroup.mk x)
  exact ⟨k, by rw [h]; exact hk⟩

/-! ## ★2. `Φ : ∏_l M[l] → M` -/

instance quotFinite (U : OpenNormalSubgroup M) : Finite (M ⧸ U.toSubgroup) :=
  inferInstanceAs (Finite (M ⧸ U.toOpenSubgroup.toSubgroup))

/-- ★`M[p]` の元は `p ∤ |M/U|` なら `M/U` で自明。 -/
theorem mk_eq_one_of_not_mem_primeFactors {p : ℕ} (hp : p.Prime)
    {z : M} (hz : z ∈ lPart M p) (U : OpenNormalSubgroup M)
    (h : p ∉ (Nat.card (M ⧸ U.toSubgroup)).primeFactors) :
    (QuotientGroup.mk z : M ⧸ U.toSubgroup) = 1 := by
  obtain ⟨k, hk⟩ := hz U
  have hord : orderOf (QuotientGroup.mk z : M ⧸ U.toSubgroup) ∣ p ^ k :=
    orderOf_dvd_of_pow_eq_one hk
  have hcard : orderOf (QuotientGroup.mk z : M ⧸ U.toSubgroup)
      ∣ Nat.card (M ⧸ U.toSubgroup) := orderOf_dvd_natCard _
  have hpn : ¬ p ∣ Nat.card (M ⧸ U.toSubgroup) := fun hd =>
    h (Nat.mem_primeFactors.mpr ⟨hp, hd, Nat.card_pos.ne'⟩)
  have hcop : Nat.Coprime (p ^ k) (Nat.card (M ⧸ U.toSubgroup)) :=
    Nat.Coprime.pow_left k ((Nat.Prime.coprime_iff_not_dvd hp).mpr hpn)
  have hd1 : orderOf (QuotientGroup.mk z : M ⧸ U.toSubgroup)
      ∣ Nat.gcd (p ^ k) (Nat.card (M ⧸ U.toSubgroup)) := Nat.dvd_gcd hord hcard
  rw [hcop] at hd1
  exact orderOf_eq_one_iff.mp (Nat.dvd_one.mp hd1)

variable (M) in
/-- ★`Φ` を与える整合族 —— 各有限商での**有限**積。 -/
noncomputable def prodFam (y : (l : Nat.Primes) → lPart M l.1) (U : OpenNormalSubgroup M) :
    M ⧸ U.toSubgroup :=
  ∏ p ∈ (Nat.card (M ⧸ U.toSubgroup)).primeFactors,
    (QuotientGroup.mk ((y (toPrimes p)).1) : M ⧸ U.toSubgroup)

set_option maxHeartbeats 1000000 in
variable (M) in
theorem prodFam_compatible (y : (l : Nat.Primes) → lPart M l.1) :
    ∀ (U V : OpenNormalSubgroup M) (h : U.toSubgroup ≤ V.toSubgroup),
      QuotientGroup.map U.toSubgroup V.toSubgroup (MonoidHom.id M) h (prodFam M y U)
        = prodFam M y V := by
  intro U V h
  have hsurj : Function.Surjective
      (QuotientGroup.map U.toSubgroup V.toSubgroup (MonoidHom.id M) h) := by
    intro z
    obtain ⟨w, rfl⟩ := QuotientGroup.mk_surjective z
    exact ⟨QuotientGroup.mk w, rfl⟩
  have hdvd : Nat.card (M ⧸ V.toSubgroup) ∣ Nat.card (M ⧸ U.toSubgroup) :=
    Subgroup.card_dvd_of_surjective _ hsurj
  have hsub : (Nat.card (M ⧸ V.toSubgroup)).primeFactors
      ⊆ (Nat.card (M ⧸ U.toSubgroup)).primeFactors :=
    Nat.primeFactors_mono hdvd Nat.card_pos.ne'
  show QuotientGroup.map _ _ _ h (∏ p ∈ _, _) = _
  rw [map_prod]
  refine (Finset.prod_subset hsub ?_).symm.trans ?_
  · intro p hpU hpV
    have hp : p.Prime := Nat.prime_of_mem_primeFactors hpU
    have hz : ((y (toPrimes p)).1 : M) ∈ lPart M p := by
      have h3 : lPart M ((toPrimes p).1) = lPart M p := by rw [toPrimes_val hp]
      exact h3 ▸ (y (toPrimes p)).2
    exact mk_eq_one_of_not_mem_primeFactors hp hz V hpV
  · refine Finset.prod_congr rfl (fun p _ => ?_)
    rfl

/-! ## ★3. `Φ` の構成と全単射性 -/

variable (M) in
/-- ★★**`Φ : ∏_l M[l] → M`**。 -/
noncomputable def prodToM (y : (l : Nat.Primes) → lPart M l.1) : M :=
  (existsUnique_mk_eq_of_compatible (asProfiniteGrp M) (prodFam M y)
    (prodFam_compatible M y)).choose

variable (M) in
theorem prodToM_spec (y : (l : Nat.Primes) → lPart M l.1) (U : OpenNormalSubgroup M) :
    (QuotientGroup.mk (prodToM M y) : M ⧸ U.toSubgroup) = prodFam M y U :=
  (existsUnique_mk_eq_of_compatible (asProfiniteGrp M) (prodFam M y)
    (prodFam_compatible M y)).choose_spec.1 U

variable (M) in
/-- ★`Φ` は群準同型。 -/
theorem prodToM_mul (y z : (l : Nat.Primes) → lPart M l.1) :
    prodToM M (y * z) = prodToM M y * prodToM M z := by
  refine eq_of_forall_quotient_eq (asProfiniteGrp M) (fun U => ?_)
  have hm : (QuotientGroup.mk (prodToM M y * prodToM M z) : M ⧸ U.toSubgroup)
      = QuotientGroup.mk (prodToM M y) * QuotientGroup.mk (prodToM M z) := rfl
  rw [prodToM_spec M (y * z) U, hm, prodToM_spec M y U, prodToM_spec M z U]
  show (∏ p ∈ _, (QuotientGroup.mk (((y * z) (toPrimes p)).1) : M ⧸ U.toSubgroup))
    = (∏ p ∈ _, _) * (∏ p ∈ _, _)
  rw [← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl (fun p _ => ?_)
  rfl

variable (M) in
/-- ★★★**成分の取り出し** —— `primProj l` を当てると `l`-成分が出る。 -/
theorem primProj_prodFam (y : (l : Nat.Primes) → lPart M l.1) (U : OpenNormalSubgroup M)
    (l : Nat.Primes) :
    primProj (M ⧸ U.toSubgroup) l.1 (prodFam M y U)
      = (QuotientGroup.mk ((y l).1) : M ⧸ U.toSubgroup) := by
  show primProj _ l.1 (∏ p ∈ _, _) = _
  rw [map_prod]
  by_cases hl : l.1 ∈ (Nat.card (M ⧸ U.toSubgroup)).primeFactors
  · rw [Finset.prod_eq_single l.1]
    · have hmem : (QuotientGroup.mk ((y (toPrimes l.1)).1) : M ⧸ U.toSubgroup)
          ∈ CommMonoid.primaryComponent (M ⧸ U.toSubgroup) l.1 := by
        have hz : ((y l).1 : M) ∈ lPart M l.1 := (y l).2
        obtain ⟨k, hk⟩ := hz U
        have hlp : toPrimes l.1 = l := by
          refine Subtype.ext ?_
          exact toPrimes_val l.2
        rw [hlp]
        exact ⟨k, hk⟩
      have hlp : toPrimes l.1 = l := Subtype.ext (toPrimes_val l.2)
      rw [primProj_eq_self_of_mem l.2 hmem, hlp]
    · intro p hp hne
      have hpp : p.Prime := Nat.prime_of_mem_primeFactors hp
      have hz : ((y (toPrimes p)).1 : M) ∈ lPart M p := by
        have h3 : lPart M ((toPrimes p).1) = lPart M p := by rw [toPrimes_val hpp]
        exact h3 ▸ (y (toPrimes p)).2
      obtain ⟨k, hk⟩ := hz U
      exact primProj_eq_one_of_mem_ne l.2 hpp hne ⟨k, hk⟩
    · intro hnm
      exact absurd hl hnm
  · have hz : ((y l).1 : M) ∈ lPart M l.1 := (y l).2
    rw [mk_eq_one_of_not_mem_primeFactors l.2 hz U hl]
    refine Finset.prod_eq_one (fun p hp => ?_)
    have hpp : p.Prime := Nat.prime_of_mem_primeFactors hp
    have hne : p ≠ l.1 := fun hE => hl (hE ▸ hp)
    have hz2 : ((y (toPrimes p)).1 : M) ∈ lPart M p := by
      have h3 : lPart M ((toPrimes p).1) = lPart M p := by rw [toPrimes_val hpp]
      exact h3 ▸ (y (toPrimes p)).2
    obtain ⟨k, hk⟩ := hz2 U
    exact primProj_eq_one_of_mem_ne l.2 hpp hne ⟨k, hk⟩

variable (M) in
/-- ★★`Φ` は単射。 -/
theorem prodToM_injective : Function.Injective (prodToM M) := by
  intro y z hyz
  funext l
  refine Subtype.ext (eq_of_forall_quotient_eq (asProfiniteGrp M) (fun U => ?_))
  have h := congrArg (fun t : M => (QuotientGroup.mk t : M ⧸ U.toSubgroup)) hyz
  rw [prodToM_spec M y U, prodToM_spec M z U] at h
  have h2 := congrArg (primProj (M ⧸ U.toSubgroup) l.1) h
  rwa [primProj_prodFam, primProj_prodFam] at h2

variable (M) in
/-- ★★`Φ` は全射。 -/
theorem prodToM_surjective : Function.Surjective (prodToM M) := by
  intro x
  refine ⟨fun l => ⟨lComp M l x, lComp_mem M l x⟩, ?_⟩
  refine eq_of_forall_quotient_eq (asProfiniteGrp M) (fun U => ?_)
  rw [prodToM_spec M _ U]
  show (∏ p ∈ (Nat.card (M ⧸ U.toSubgroup)).primeFactors,
    (QuotientGroup.mk (lComp M (toPrimes p) x) : M ⧸ U.toSubgroup)) = QuotientGroup.mk x
  have hterm : ∀ p ∈ (Nat.card (M ⧸ U.toSubgroup)).primeFactors,
      (QuotientGroup.mk (lComp M (toPrimes p) x) : M ⧸ U.toSubgroup)
        = primProj (M ⧸ U.toSubgroup) p (QuotientGroup.mk x) := by
    intro p hp
    have hpp : p.Prime := Nat.prime_of_mem_primeFactors hp
    rw [lComp_spec, toPrimes_val hpp]
  rw [Finset.prod_congr rfl hterm]
  exact prod_primProj _

/-! ## ★4. 連続性と同相 -/

variable (M) in
/-- ★`y ↦ prodFam M y U` は連続(**有限**積、行き先は離散)。 -/
theorem prodFam_continuous (U : OpenNormalSubgroup M) :
    Continuous (fun y : (l : Nat.Primes) → lPart M l.1 => prodFam M y U) := by
  refine continuous_finset_prod _ (fun p _ => ?_)
  exact continuous_quotient_mk'.comp
    ((continuous_subtype_val).comp (continuous_apply (toPrimes p)))

variable (M) in
set_option maxHeartbeats 1000000 in
/-- ★★★★**`Φ` は連続** —— 開正規部分群の剰余類が開基をなすことだけを使う。 -/
theorem prodToM_continuous : Continuous (prodToM M) := by
  rw [continuous_def]
  intro V hV
  rw [isOpen_iff_forall_mem_open]
  intro y₀ hy₀
  -- ★`x₀ := Φ y₀ ∈ V`。`x₀⁻¹ • V` は `1` の開近傍。
  set x₀ := prodToM M y₀ with hx₀
  have hmem : x₀ ∈ V := hy₀
  have hopen : IsOpen ((fun t : M => x₀ * t) ⁻¹' V) :=
    hV.preimage (continuous_const.mul continuous_id)
  have h1 : (1 : M) ∈ (fun t : M => x₀ * t) ⁻¹' V := by
    simpa using hmem
  obtain ⟨U, hU⟩ := ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one
    (G := M) hopen h1
  refine ⟨{y | prodFam M y U = prodFam M y₀ U}, ?_, ?_, ?_⟩
  · -- ★中身は `V` に入る
    intro y hy
    have hq : (QuotientGroup.mk (prodToM M y) : M ⧸ U.toSubgroup)
        = QuotientGroup.mk x₀ := by
      rw [prodToM_spec M y U, hx₀, prodToM_spec M y₀ U]
      exact hy
    have hdiv : x₀⁻¹ * prodToM M y ∈ U.toSubgroup := by
      have h2 : (prodToM M y) / x₀ ∈ U.toSubgroup := QuotientGroup.eq_iff_div_mem.mp hq
      rwa [div_eq_mul_inv, mul_comm] at h2
    have : x₀ * (x₀⁻¹ * prodToM M y) ∈ V := hU hdiv
    simpa using this
  · -- ★開集合
    exact (isOpen_discrete {prodFam M y₀ U}).preimage (prodFam_continuous M U)
  · exact rfl

variable (M) in
/-- ★`Φ` を群準同型として。 -/
noncomputable def prodToMHom : ((l : Nat.Primes) → lPart M l.1) →* M :=
  MonoidHom.mk' (prodToM M) (prodToM_mul M)

variable (M) in
/-- ★`Φ` は群同型。 -/
noncomputable def prodToMEquiv : ((l : Nat.Primes) → lPart M l.1) ≃* M :=
  MulEquiv.ofBijective (prodToMHom M) ⟨prodToM_injective M, prodToM_surjective M⟩

variable (M) in
/-- ★★★★★**`∏_l M[l] ≅ M`**(位相群の同型) —— コンパクト→Hausdorff の連続全単射。 -/
noncomputable def prodToMTopEquiv : ((l : Nat.Primes) → lPart M l.1) ≃ₜ* M :=
  { prodToMEquiv M with
    continuous_toFun := prodToM_continuous M
    continuous_invFun :=
      (Continuous.homeoOfEquivCompactToT2 (f := (prodToMEquiv M).toEquiv)
        (prodToM_continuous M)).continuous_invFun }

variable (M) in
/-- ★★★★★★**[FrdI] Definition 2.8, (ii)** —— `M ≅ ∏_l M[l]`(位相群の同型)。 -/
noncomputable def decompEquiv : M ≃ₜ* ((l : Nat.Primes) → lPart M l.1) :=
  (prodToMTopEquiv M).symm

def decompEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (ii) — M ≅ ∏_l M[l]",
    sectionId := "frdi-def-2-8" }

end CommProfinite

end ABC3.Found.ProL
