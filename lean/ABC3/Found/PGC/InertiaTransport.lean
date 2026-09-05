import ABC3.Found.PGC.InertiaIdentification
import ABC3.Found.PGC.AdjoinFieldClosure

/-!
# Corollary 1.3 は Proposition 1.2 に帰着する

原文 (pGC p.3) の Corollary 1.3 の論拠は

> By applying Proposition 1.2 to L and H, we see that the number q_L of elements in the
> residue field of O[scr]_L can be recovered group-theoretically from H ⊆ Γ_K.

——つまり「Prop 1.2 を各開部分群 `H` とその固定体 `L_H` に適用する」。
`Skeleton/PGC/Section1Cor13.lean::inertia_recoverable.needs` はこれを
`.derivation "Proposition 1.2 を (L, H) に適用して q_L を得る段"` として
記録していた。**本ファイルはその段を形式化する**。

## 主張

`α : Γ_K ≅ Γ_{K'}`(位相群として)に沿って

* `q_K = q_{K'}`(Prop 1.2 の `q` の半分、`K` 自身に適用)
* `q_{L_H} = q_{L_{α(H)}}`(同じく各開部分群に適用)

が成り立てば、**惰性群は `α` で移る**——すなわち Corollary 1.3 が従う。

判定条件 `q_{L_H} = q_K^{[Γ_K:H]}` の三つの量がすべて `α` で保たれる
(`q_K`・`q_{L_H}` は仮定、指数は `Subgroup.index_map_of_bijective`)ので、
`α` は「判定条件を満たす開部分群」の集合を全単射に移し、したがって
その `sInf` である `I_K` も移る。

## ★★★2026-09-05: 第二の仮定は**消えた**——Cor 1.3 ⟸ Prop 1.2(仮定なし)

上に「`Γ_{L_H} ≅ H` を位相群の同型として使うには Krull 位相の比較
(`fixingSubgroupEquiv` が同相であること)がもう一段要る」と書いていた。
その一段が `AdjoinFieldClosure.lean` で付いた
(`continuous_fixingSubgroupEquiv` / `continuous_fixingSubgroupEquiv_symm` /
`absGalFixedFieldCME : Γ_{L_H} ≃ₜ* H`)。

したがって **`htr`(第二の仮定)は `htop`(Prop 1.2 そのもの)から従う**
——`residueCard_transport_of_prop12`。`α` を `H` に制限して
`Γ_{L_H} ≃ₜ* H ≃ₜ* α(H) ≃ₜ* Γ_{L_{α(H)}}` を作り、そこに Prop 1.2 を当てる。
これが原文の "By applying Proposition 1.2 to L and H" そのもの。

結論(`inertia_recoverable_of_prop12`):

> **Prop 1.2(`q` が位相群 `Γ_K` から決まる)が成り立てば、Cor 1.3 が成り立つ。**

残るのは Prop 1.2 のみ。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- 位相群の同型は開部分群を開部分群に移す。 -/
theorem isOpen_map_of_continuousMulEquiv {G H : Type*} [Group G] [Group H]
    [TopologicalSpace G] [TopologicalSpace H]
    (α : ContinuousMulEquiv G H) (S : Subgroup G) (hS : IsOpen (S : Set G)) :
    IsOpen ((S.map α.toMulEquiv.toMonoidHom : Subgroup H) : Set H) := by
  have h : ((S.map α.toMulEquiv.toMonoidHom : Subgroup H) : Set H) = α '' (S : Set G) := by
    ext y
    constructor
    · rintro ⟨x, hx, rfl⟩; exact ⟨x, hx, rfl⟩
    · rintro ⟨x, hx, rfl⟩; exact ⟨x, hx, rfl⟩
  rw [h]
  exact α.toHomeomorph.isOpenMap _ hS

theorem isOpen_comap_of_continuousMulEquiv {G H : Type*} [Group G] [Group H]
    [TopologicalSpace G] [TopologicalSpace H]
    (α : ContinuousMulEquiv G H) (S : Subgroup H) (hS : IsOpen (S : Set H)) :
    IsOpen ((S.comap α.toMulEquiv.toMonoidHom : Subgroup G) : Set G) := by
  have h : ((S.comap α.toMulEquiv.toMonoidHom : Subgroup G) : Set G) = α ⁻¹' (S : Set H) := rfl
  rw [h]
  exact hS.preimage α.continuous_toFun

/-- **★★★★★★★惰性群は `q` の移送さえあれば `α` で移る**
——Corollary 1.3 の、原文が「Prop 1.2 を (L, H) に適用する」と言っている段。 -/
theorem transport_inertia_of_residueCard_transport (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal)
    (htop : Nat.card 𝓀[K.carrier] = Nat.card 𝓀[K'.carrier])
    (htr : ∀ (H : Subgroup K.absGal) (hH : IsOpen (H : Set K.absGal))
      (hH' : IsOpen ((H.map α.toMulEquiv.toMonoidHom : Subgroup K'.absGal) : Set K'.absGal)),
      (residueCardinality p).card ((subgroupCorrespondence p).field K H hH)
        = (residueCardinality p).card ((subgroupCorrespondence p).field K'
            (H.map α.toMulEquiv.toMonoidHom) hH')) :
    Subgroup.map α.toMulEquiv.toMonoidHom
        (inertia (residueCardinality p) (subgroupCorrespondence p) K)
      = inertia (residueCardinality p) (subgroupCorrespondence p) K' := by
  have hcrit2 : ∀ (H : Subgroup K.absGal) (hH : IsOpen (H : Set K.absGal))
      (H' : Subgroup K'.absGal) (hH' : IsOpen (H' : Set K'.absGal))
      (_heq : H.map α.toMulEquiv.toMonoidHom = H'),
      (IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K H hH ↔
        IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K' H' hH') := by
    intro H hH H' hH' heq
    subst heq
    rw [IsUnramifiedAt, IsUnramifiedAt, htr H hH hH']
    simp only [residueCardinality_card]
    rw [htop, Subgroup.index_map_of_bijective α.bijective]
  have hinj : Function.Injective α.toMulEquiv.toMonoidHom := α.injective
  have hsurj : Function.Surjective α.toMulEquiv.toMonoidHom := α.surjective
  refine le_antisymm ?_ ?_
  · refine le_sInf ?_
    rintro H' ⟨hH', hcrit'⟩
    rw [Subgroup.map_le_iff_le_comap]
    refine sInf_le ⟨isOpen_comap_of_continuousMulEquiv α H' hH', ?_⟩
    exact (hcrit2 _ (isOpen_comap_of_continuousMulEquiv α H' hH') H' hH'
      (Subgroup.map_comap_eq_self_of_surjective hsurj H')).mpr hcrit'
  · have hkey : Subgroup.comap α.toMulEquiv.toMonoidHom
        (inertia (residueCardinality p) (subgroupCorrespondence p) K')
        ≤ inertia (residueCardinality p) (subgroupCorrespondence p) K := by
      refine le_sInf ?_
      rintro H ⟨hH, hcritH⟩
      have hmem : (H.map α.toMulEquiv.toMonoidHom) ∈ {J : Subgroup K'.absGal |
          ∃ hJ : IsOpen (J : Set K'.absGal),
            IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K' J hJ} :=
        ⟨isOpen_map_of_continuousMulEquiv α H hH,
          (hcrit2 H hH _ (isOpen_map_of_continuousMulEquiv α H hH) rfl).mp hcritH⟩
      calc Subgroup.comap α.toMulEquiv.toMonoidHom
            (inertia (residueCardinality p) (subgroupCorrespondence p) K')
          ≤ Subgroup.comap α.toMulEquiv.toMonoidHom (H.map α.toMulEquiv.toMonoidHom) :=
            Subgroup.comap_mono (sInf_le hmem)
        _ = H := Subgroup.comap_map_eq_self_of_injective hinj H
    calc inertia (residueCardinality p) (subgroupCorrespondence p) K'
        = Subgroup.map α.toMulEquiv.toMonoidHom (Subgroup.comap α.toMulEquiv.toMonoidHom
            (inertia (residueCardinality p) (subgroupCorrespondence p) K')) :=
          (Subgroup.map_comap_eq_self_of_surjective hsurj _).symm
      _ ≤ Subgroup.map α.toMulEquiv.toMonoidHom
            (inertia (residueCardinality p) (subgroupCorrespondence p) K) :=
          Subgroup.map_mono hkey

/-- **★★★★★★★★Corollary 1.3 は「`q` が `α` で移る」に帰着する**
——`Skeleton/PGC/Section1Cor13.lean::inertia_recoverable.needs` の
`.derivation "Proposition 1.2 を (L, H) に適用して q_L を得る段"` の形式化。 -/
theorem inertia_recoverable_of_residueCard_transport
    (htop : ∀ (K K' : PAdicLocalField p) (_α : ContinuousMulEquiv K.absGal K'.absGal),
      Nat.card 𝓀[K.carrier] = Nat.card 𝓀[K'.carrier])
    (htr : ∀ (K K' : PAdicLocalField p) (α : ContinuousMulEquiv K.absGal K'.absGal)
      (H : Subgroup K.absGal) (hH : IsOpen (H : Set K.absGal))
      (hH' : IsOpen ((H.map α.toMulEquiv.toMonoidHom : Subgroup K'.absGal) : Set K'.absGal)),
      (residueCardinality p).card ((subgroupCorrespondence p).field K H hH)
        = (residueCardinality p).card ((subgroupCorrespondence p).field K'
            (H.map α.toMulEquiv.toMonoidHom) hH')) :
    (inertiaObject (residueCardinality p)
      (subgroupCorrespondence p)).RecoverableFromAbsGal := by
  intro K K' α
  exact transport_inertia_of_residueCard_transport K K' α (htop K K' α) (htr K K' α)

/-- 位相群の同型を部分群に制限したもの。 -/
noncomputable def subgroupMapCME {G G' : Type*} [Group G] [Group G']
    [TopologicalSpace G] [TopologicalSpace G'] (α : ContinuousMulEquiv G G')
    (S : Subgroup G) : ContinuousMulEquiv S (S.map α.toMulEquiv.toMonoidHom) where
  toMulEquiv := Subgroup.equivMapOfInjective S α.toMulEquiv.toMonoidHom α.injective
  continuous_toFun := continuous_induced_rng.2 (α.continuous_toFun.comp continuous_subtype_val)
  continuous_invFun := by
    refine continuous_induced_rng.2 ?_
    refine (α.continuous_invFun.comp continuous_subtype_val).congr ?_
    intro y
    show α.symm (y : G') = _
    apply α.injective
    rw [ContinuousMulEquiv.apply_symm_apply]
    exact (congrArg Subtype.val
      ((Subgroup.equivMapOfInjective S α.toMulEquiv.toMonoidHom
        α.injective).apply_symm_apply y)).symm

/-- **★★★★★★★★★`htr` は Prop 1.2 から従う**——原文の
"By applying Proposition 1.2 to L and H" そのもの。

`H ≠ ⊤` のとき `L_H = fixedFieldLocalField K H`、`L_{α(H)} = fixedFieldLocalField K' α(H)`
で、`absGalFixedFieldCME`(`Γ_{L_H} ≃ₜ* H`)と `α` の制限を繋げば
`Γ_{L_H} ≃ₜ* Γ_{L_{α(H)}}` が得られる。`H = ⊤` のときは `L_⊤ = K` なので
Prop 1.2 を `K`・`K'` に当てるだけ。 -/
theorem residueCard_transport_of_prop12
    (htop : ∀ (K K' : PAdicLocalField p) (_α : ContinuousMulEquiv K.absGal K'.absGal),
      Nat.card 𝓀[K.carrier] = Nat.card 𝓀[K'.carrier])
    (K K' : PAdicLocalField p) (α : ContinuousMulEquiv K.absGal K'.absGal)
    (H : Subgroup K.absGal) (hH : IsOpen (H : Set K.absGal))
    (hH' : IsOpen ((H.map α.toMulEquiv.toMonoidHom : Subgroup K'.absGal) : Set K'.absGal)) :
    (residueCardinality p).card ((subgroupCorrespondence p).field K H hH)
      = (residueCardinality p).card ((subgroupCorrespondence p).field K'
          (H.map α.toMulEquiv.toMonoidHom) hH') := by
  classical
  by_cases hTop : H = ⊤
  · have h2 : H.map α.toMulEquiv.toMonoidHom = ⊤ := by
      rw [hTop]; exact Subgroup.map_top_of_surjective _ α.surjective
    have e1 : (subgroupCorrespondence p).field K H hH = K := if_pos hTop
    have e2 : (subgroupCorrespondence p).field K' (H.map α.toMulEquiv.toMonoidHom) hH' = K' :=
      if_pos h2
    rw [e1, e2]
    simp only [residueCardinality_card]
    exact htop K K' α
  · have h2 : H.map α.toMulEquiv.toMonoidHom ≠ ⊤ := by
      intro hc
      apply hTop
      have hcm := Subgroup.comap_map_eq_self_of_injective
        (f := α.toMulEquiv.toMonoidHom) α.injective H
      rw [hc] at hcm
      rw [← hcm]
      exact Subgroup.comap_top _
    have e1 : (subgroupCorrespondence p).field K H hH = fixedFieldLocalField K H hH := if_neg hTop
    have e2 : (subgroupCorrespondence p).field K' (H.map α.toMulEquiv.toMonoidHom) hH'
        = fixedFieldLocalField K' (H.map α.toMulEquiv.toMonoidHom) hH' := if_neg h2
    rw [e1, e2]
    simp only [residueCardinality_card]
    exact htop _ _ ((absGalFixedFieldCME K H hH).trans
      ((subgroupMapCME α H).trans (absGalFixedFieldCME K' _ hH').symm))

/-- **★★★★★★★★★★★★[pGC] Corollary 1.3 ⟸ Proposition 1.2**(他に仮定なし)。

原文の論拠「Prop 1.2 を (L, H) に適用する」を、そのまま Lean の含意にした形。 -/
theorem inertia_recoverable_of_prop12
    (htop : ∀ (K K' : PAdicLocalField p) (_α : ContinuousMulEquiv K.absGal K'.absGal),
      Nat.card 𝓀[K.carrier] = Nat.card 𝓀[K'.carrier]) :
    (inertiaObject (residueCardinality p)
      (subgroupCorrespondence p)).RecoverableFromAbsGal :=
  inertia_recoverable_of_residueCard_transport htop (residueCard_transport_of_prop12 htop)

end ABC3.Found.PGC
