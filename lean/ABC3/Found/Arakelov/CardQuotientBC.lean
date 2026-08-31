/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.Ideal.Norm.RelNorm
import Mathlib.NumberTheory.NumberField.Basic
import ABC3.Found.Arakelov.TensorIndex
import ABC3.Found.Arakelov.DegArith
import ABC3.Meta.Claim

/-!
# 底変換の**数え上げ**（巡回加群の場合）——`#(𝓞_K/𝔞𝓞_K) = #(𝓞_F/𝔞)^{[K:F]}`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★これは何か

底変換 `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の**有限側**は、`§9-791` の同型

    `Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L)`

と合わせると

    `#(𝓞_K ⊗_{𝓞_F} Q) = #(Q)^{[K:F]}`   （`Q` は有限 `𝓞_F`-加群）

に帰着する。★本ファイルはその**巡回加群の場合**（`Q = 𝓞_F/𝔞`）を入れる。

## ★★★測定の記録（2026-08-28）

mathlib に **`Ideal.absNorm_algebraMap`** がある:

    `absNorm (I.map (algebraMap R S)) = (absNorm I) ^ finrank (FractionRing R) (FractionRing S)`

★`Ideal.absNorm I = Nat.card (R ⧸ I)` は定義そのものなので、そのまま数え上げになる。
★★使うには `attribute [local instance] FractionRing.liftAlgebra` が要る
——mathlib の当該節がそうしている。

## ★残っている段（明示）

★★`finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K)) = [K:F]` の同一視。
★★★一般の有限 `Q` への持ち上げ（巡回商による濾過での dévissage）。
-/

namespace ABC3.Found.Arakelov

open NumberField

attribute [local instance] FractionRing.liftAlgebra in
/-- ★★★★★★**イデアルを延長したときの商の位数**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

    `#(𝓞_K / 𝔞·𝓞_K) = #(𝓞_F / 𝔞) ^ [FractionRing 𝓞_K : FractionRing 𝓞_F]`

★mathlib の `Ideal.absNorm_algebraMap` そのままである
（`Ideal.absNorm I = Nat.card (R ⧸ I)` は定義）。 -/
theorem card_quotient_map_algebraMap (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (a : Ideal (𝓞 F)) :
    Nat.card ((𝓞 K) ⧸ (a.map (algebraMap (𝓞 F) (𝓞 K))))
      = (Nat.card ((𝓞 F) ⧸ a))
        ^ (Module.finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K))) :=
  Ideal.absNorm_algebraMap (R := 𝓞 F) (S := 𝓞 K) a

attribute [local instance] FractionRing.liftAlgebra FractionRing.isScalarTower_liftAlgebra in
/-- ★★★★★**`FractionRing` の次数は体の次数である**。

★`FractionRing (𝓞 F) ≃ₐ F`（`FractionRing.algEquiv`）を 2 つ使って
`Algebra.finrank_eq_of_equiv_equiv` に渡す。
★★互換性は局所化の一意性（`IsLocalization.ringHom_ext`）で出る。 -/
theorem finrank_fractionRing_ringOfIntegers (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] :
    Module.finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K)) = Module.finrank F K := by
  refine Algebra.finrank_eq_of_equiv_equiv
    (FractionRing.algEquiv (𝓞 F) F).toRingEquiv
    (FractionRing.algEquiv (𝓞 K) K).toRingEquiv ?_
  apply IsLocalization.ringHom_ext (nonZeroDivisors (𝓞 F))
  ext x
  simp only [RingHom.coe_comp, Function.comp_apply, RingEquiv.toRingHom_eq_coe,
    RingEquiv.coe_toRingHom, AlgEquiv.coe_ringEquiv, AlgEquiv.commutes]
  rw [← IsScalarTower.algebraMap_apply (𝓞 F) (FractionRing (𝓞 F)) (FractionRing (𝓞 K)),
    IsScalarTower.algebraMap_apply (𝓞 F) (𝓞 K) (FractionRing (𝓞 K)),
    AlgEquiv.commutes,
    ← IsScalarTower.algebraMap_apply (𝓞 F) F K,
    IsScalarTower.algebraMap_apply (𝓞 F) (𝓞 K) K]

attribute [local instance] FractionRing.liftAlgebra in
/-- ★★★★★★★**イデアルを延長したときの商の位数**（`[K:F]` で書いた形）。

    `#(𝓞_K / 𝓞·𝓞_K) = #(𝓞_F / 𝓞) ^ [K:F]` -/
theorem card_quotient_map_algebraMap' (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (a : Ideal (𝓞 F)) :
    Nat.card ((𝓞 K) ⧸ (a.map (algebraMap (𝓞 F) (𝓞 K))))
      = (Nat.card ((𝓞 F) ⧸ a)) ^ (Module.finrank F K) := by
  rw [card_quotient_map_algebraMap F K a, finrank_fractionRing_ringOfIntegers F K]

/-! ## ★★★★★★★★★★一般の有限加群へ（dévissage） -/

open scoped TensorProduct

attribute [local instance] FractionRing.liftAlgebra in
/-- ★★★★★★**巡回加群の場合**——`𝓞_K ⊗ (𝓞_F/𝔞) ≅ 𝓞_K/𝔞𝓞_K`。 -/
theorem card_tensor_cyclic (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (a : Ideal (𝓞 F)) (ha : a ≠ ⊥) :
    Finite (TensorProduct (𝓞 F) (𝓞 K) ((𝓞 F) ⧸ a)) ∧
    Nat.card (TensorProduct (𝓞 F) (𝓞 K) ((𝓞 F) ⧸ a))
      = (Nat.card ((𝓞 F) ⧸ a)) ^ (Module.finrank F K) := by
  have hmapne : (a.map (algebraMap (𝓞 F) (𝓞 K))) ≠ ⊥ := by
    obtain ⟨x, hxa, hx0⟩ := Submodule.exists_mem_ne_zero_of_ne_bot ha
    intro h
    have hx : (algebraMap (𝓞 F) (𝓞 K)) x ∈ a.map (algebraMap (𝓞 F) (𝓞 K)) :=
      Ideal.mem_map_of_mem _ hxa
    rw [h, Ideal.mem_bot] at hx
    exact hx0 ((map_eq_zero_iff _ (FaithfulSMul.algebraMap_injective (𝓞 F) (𝓞 K))).mp hx)
  haveI hfinq : Finite ((𝓞 K) ⧸ (a.map (algebraMap (𝓞 F) (𝓞 K)))) :=
    Ideal.finiteQuotientOfFreeOfNeBot _ hmapne
  have e1 : TensorProduct (𝓞 F) (𝓞 K) ((𝓞 F) ⧸ a)
      ≃ₗ[𝓞 F] TensorProduct (𝓞 F) ((𝓞 F) ⧸ a) (𝓞 K) := TensorProduct.comm _ _ _
  have e2 : TensorProduct (𝓞 F) ((𝓞 F) ⧸ a) (𝓞 K)
      ≃ₗ[𝓞 F] ((𝓞 K) ⧸ (a • (⊤ : Submodule (𝓞 F) (𝓞 K)))) :=
    TensorProduct.quotTensorEquivQuotSMul (𝓞 K) a
  have e3 : ((𝓞 K) ⧸ (a • (⊤ : Submodule (𝓞 F) (𝓞 K))))
      ≃ₗ[𝓞 F] ((𝓞 K) ⧸ ((a.map (algebraMap (𝓞 F) (𝓞 K))).restrictScalars (𝓞 F))) :=
    Submodule.quotEquivOfEq _ _ (Ideal.smul_top_eq_map a)
  have hequiv : Nonempty (TensorProduct (𝓞 F) (𝓞 K) ((𝓞 F) ⧸ a)
      ≃ ((𝓞 K) ⧸ (a.map (algebraMap (𝓞 F) (𝓞 K))))) :=
    ⟨(e1.trans (e2.trans e3)).toEquiv⟩
  constructor
  · exact Finite.of_equiv _ hequiv.some.symm
  · rw [Nat.card_congr hequiv.some, card_quotient_map_algebraMap' F K a]

/-- ★★★★★★★**短完全列に沿った位数の乗法性**（`T` は平坦）。 -/
theorem card_tensor_quot (R T Q : Type) [CommRing R] [CommRing T] [Algebra R T] [Module.Flat R T]
    [AddCommGroup Q] [Module R Q] (N : Submodule R Q)
    (hN : Finite (T ⊗[R] N)) (hQN : Finite (T ⊗[R] (Q ⧸ N))) :
    Finite (T ⊗[R] Q) ∧
      Nat.card (T ⊗[R] Q) = Nat.card (T ⊗[R] N) * Nat.card (T ⊗[R] (Q ⧸ N)) := by
  haveI := hN
  haveI := hQN
  have hinj : Function.Injective (LinearMap.lTensor T N.subtype) :=
    Module.Flat.lTensor_preserves_injective_linearMap _ N.subtype_injective
  have hsurj : Function.Surjective (LinearMap.lTensor T N.mkQ) :=
    LinearMap.lTensor_surjective T (Submodule.mkQ_surjective N)
  have hker := lTensor_mkQ (R := R) (M := Q) T N
  have e1 : (T ⊗[R] N) ≃ₗ[R] (LinearMap.range (LinearMap.lTensor T N.subtype)) :=
    LinearEquiv.ofInjective _ hinj
  have e2 : ((T ⊗[R] Q) ⧸ LinearMap.ker (LinearMap.lTensor T N.mkQ)) ≃ₗ[R] (T ⊗[R] (Q ⧸ N)) :=
    LinearMap.quotKerEquivOfSurjective _ hsurj
  haveI hr : Finite (LinearMap.range (LinearMap.lTensor T N.subtype)) :=
    Finite.of_equiv _ e1.toEquiv
  have e3 : ((T ⊗[R] Q) ⧸ LinearMap.range (LinearMap.lTensor T N.subtype))
      ≃ₗ[R] (T ⊗[R] (Q ⧸ N)) := by
    rw [← hker]; exact e2
  haveI hq : Finite ((T ⊗[R] Q)
      ⧸ (LinearMap.range (LinearMap.lTensor T N.subtype)).toAddSubgroup) :=
    Finite.of_equiv _ e3.toEquiv.symm
  haveI hfin : Finite (T ⊗[R] Q) :=
    Finite.of_addSubgroup_quotient
      (LinearMap.range (LinearMap.lTensor T N.subtype)).toAddSubgroup
  refine ⟨hfin, ?_⟩
  have hcard := AddSubgroup.card_eq_card_quotient_mul_card_addSubgroup
    (LinearMap.range (LinearMap.lTensor T N.subtype)).toAddSubgroup
  have h2 : Nat.card ((T ⊗[R] Q)
      ⧸ (LinearMap.range (LinearMap.lTensor T N.subtype)).toAddSubgroup)
      = Nat.card (T ⊗[R] (Q ⧸ N)) := Nat.card_congr e3.toEquiv
  have h3 : Nat.card ((LinearMap.range (LinearMap.lTensor T N.subtype)).toAddSubgroup : Type)
      = Nat.card (T ⊗[R] N) := (Nat.card_congr e1.toEquiv).symm
  rw [hcard, h2, h3]
  ring

/-- ★★★★★★★★★**dévissage の本体**（`Nat.card Q` についての帰納法）。 -/
theorem card_tensor_of_finite_aux (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] :
    ∀ (m : ℕ) (Q : Type) [AddCommGroup Q] [Module (𝓞 F) Q], Finite Q → Nat.card Q ≤ m →
      Finite ((𝓞 K) ⊗[𝓞 F] Q) ∧
      Nat.card ((𝓞 K) ⊗[𝓞 F] Q) = (Nat.card Q) ^ (Module.finrank F K) := by
  intro m
  induction m with
  | zero =>
      intro Q _ _ hfin hle
      haveI := hfin
      exfalso
      have : 0 < Nat.card Q := Nat.card_pos
      omega
  | succ m ih =>
      intro Q _ _ hfin hle
      haveI := hfin
      rcases subsingleton_or_nontrivial Q with hs | hns
      · haveI := hs
        haveI : Subsingleton ((𝓞 K) ⊗[𝓞 F] Q) := by infer_instance
        refine ⟨inferInstance, ?_⟩
        have h1 : Nat.card ((𝓞 K) ⊗[𝓞 F] Q) = 1 :=
          Nat.card_eq_one_iff_unique.mpr ⟨inferInstance, ⟨0⟩⟩
        have h2 : Nat.card Q = 1 := Nat.card_eq_one_iff_unique.mpr ⟨hs, ⟨0⟩⟩
        rw [h1, h2, one_pow]
      · obtain ⟨q, hq⟩ := exists_ne (0 : Q)
        set f := LinearMap.toSpanSingleton (𝓞 F) Q q with hf
        set N : Submodule (𝓞 F) Q := LinearMap.range f with hNdef
        haveI hQNfin : Finite (Q ⧸ N) := Finite.of_surjective _ (Submodule.mkQ_surjective N)
        have eN : ((𝓞 F) ⧸ LinearMap.ker f) ≃ₗ[𝓞 F] N := LinearMap.quotKerEquivRange f
        haveI hNfin : Finite N := Subtype.finite
        have hker : LinearMap.ker f ≠ ⊥ := by
          intro h
          have e0 : ((𝓞 F) ⧸ LinearMap.ker f) ≃ₗ[𝓞 F] (𝓞 F) := Submodule.quotEquivOfEqBot _ h
          haveI : Finite ((𝓞 F) ⧸ LinearMap.ker f) := Finite.of_equiv _ eN.toEquiv.symm
          haveI : Finite (𝓞 F) := Finite.of_equiv _ e0.toEquiv
          exact (not_finite (𝓞 F))
        have hcyc := card_tensor_cyclic F K (LinearMap.ker f) hker
        haveI := hcyc.1
        have eT : ((𝓞 K) ⊗[𝓞 F] ((𝓞 F) ⧸ LinearMap.ker f)) ≃ₗ[𝓞 F] ((𝓞 K) ⊗[𝓞 F] N) :=
          TensorProduct.congr (LinearEquiv.refl _ _) eN
        haveI hTN : Finite ((𝓞 K) ⊗[𝓞 F] N) := Finite.of_equiv _ eT.toEquiv
        have hTNcard : Nat.card ((𝓞 K) ⊗[𝓞 F] N) = (Nat.card N) ^ (Module.finrank F K) := by
          rw [← Nat.card_congr eT.toEquiv, hcyc.2, Nat.card_congr eN.toEquiv]
        have hsplit : Nat.card Q = Nat.card (Q ⧸ N) * Nat.card N :=
          AddSubgroup.card_eq_card_quotient_mul_card_addSubgroup N.toAddSubgroup
        have hqN : q ∈ N := ⟨1, by simp [hf, LinearMap.toSpanSingleton]⟩
        have hN2 : 2 ≤ Nat.card N := by
          have hnt : Nontrivial N := ⟨⟨⟨q, hqN⟩, 0, by simpa using hq⟩⟩
          have h1 : Nat.card N ≠ 1 := by
            intro h
            rw [Nat.card_eq_one_iff_unique] at h
            exact (not_subsingleton N) h.1
          have h0 : 0 < Nat.card N := Nat.card_pos
          omega
        have hQNlt : Nat.card (Q ⧸ N) ≤ m := by
          have h0 : 0 < Nat.card (Q ⧸ N) := Nat.card_pos
          nlinarith [hsplit, hle, hN2, h0]
        obtain ⟨hfq, hcq⟩ := ih (Q ⧸ N) hQNfin hQNlt
        obtain ⟨hfinal, hcard⟩ := card_tensor_quot (𝓞 F) (𝓞 K) Q N hTN hfq
        refine ⟨hfinal, ?_⟩
        rw [hcard, hTNcard, hcq, hsplit, mul_pow]
        ring

/-- ★★★★★★★★★★**底変換の数え上げ**——`#(𝓞_K ⊗_{𝓞_F} Q) = #(Q)^{[K:F]}`。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これと `§9-791` の同型 `Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L)` を合わせると
底変換の**有限側**が出る。 -/
theorem card_tensor_ringOfIntegers (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (Q : Type) [AddCommGroup Q] [Module (𝓞 F) Q]
    [Finite Q] :
    Nat.card ((𝓞 K) ⊗[𝓞 F] Q) = (Nat.card Q) ^ (Module.finrank F K) :=
  (card_tensor_of_finite_aux F K (Nat.card Q) Q inferInstance le_rfl).2

/-! ### ★出典の紐付け(`.src`) -/

def card_quotient_map_algebraMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(イデアルを延長したときの商の位数——底変換の数え上げの巡回の場合)",
    sectionId := "genell-def-1-1-ii" }

def card_quotient_map_algebraMap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.absNorm_algebraMap(absNorm (I.map) = absNorm I ^ finrank)"
      (.inMathlib "Ideal.absNorm_algebraMap") 4,
    .implicitStep
      ("★残っている段: (a) finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K)) = [K:F] の同一視、" ++
       "(b) 一般の有限 Q への持ち上げ(巡回商による濾過での dévissage)。" ++
       "★★これと §9-791 の同型 Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L) を合わせると" ++
       "底変換の有限側が出る") 4 ]

end ABC3.Found.Arakelov
