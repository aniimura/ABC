/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SplitDensity
import ABC3.Found.NumberField.SplitExponent
import ABC3.Found.NumberField.SplCompositum
import ABC3.Found.NumberField.DirichletDensity
import Mathlib.NumberTheory.Padics.HeightOneSpectrum

/-!
# `ℚ` の高さ 1 素イデアルと有理素数の同一視(鎖 `cheb` の配管)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★何のためか

`Found/NumberField/` の完全分解の結果は、素点を **`Nat.Primes`**(有理素数)で添字づけて
述べてある(`SplQ`)。一方 `Skeleton/NumberField/Chebotarev.lean` は
**`HeightOneSpectrum (𝓞 ℚ)`** で述べている。★本ファイルはその 2 つの語彙を繋ぎ、
`Skeleton` 側が 1 行で委譲できるようにする(判断 D11、2026-09-06)。

## ★★橋の中身は 3 つだけ

| 補題 | 中身 |
|---|---|
| `comap_asIdeal` | `𝔭 ⊆ 𝓞 ℚ` を `ℤ` へ引き戻すと `(p)` になる(`natGenerator`) |
| `absNorm_asIdeal` | `N𝔭 = p`(剰余環を `ℤ ⧸ (p) ≃ ZMod p` へ移す) |
| `primesOver_asIdeal_eq` | `𝔭` の上にある `𝓞 L` の素イデアルと `(p)` の上にあるものは同じ |

★`ℤ ≃+* 𝓞 ℚ` は mathlib の `Rat.IsIntegralClosure.intEquiv`、
素点の対応は `Rat.HeightOneSpectrum.primesEquiv` である
(`Mathlib/NumberTheory/Padics/HeightOneSpectrum.lean`)。
★`ℤ` から任意の環への環準同型は一意(`RingHom.ext_int`)なので、
`intEquiv` と `algebraMap ℤ (𝓞 ℚ)` の同一視は 1 行で済む。

## ★底は `ℚ` だけである

[FrdI] Theorem 6.4, (iv) の底は `ℚ` であり、一般の底 `K` の Chebotarev は原典が要求しない。
★本ファイルも `ℚ` の場合しか作らない。
-/

namespace ABC3.Found.NF

open _root_.NumberField IsDedekindDomain Ideal Filter Topology
open scoped _root_.NumberField

/-! ## ★1. `ℤ` と `𝓞 ℚ` の同一視 -/

/-- ★`algebraMap ℤ (𝓞 ℚ)` は全射(`ℤ ≃+* 𝓞 ℚ` だから)。 -/
theorem algebraMap_int_ringOfIntegersRat_surjective :
    Function.Surjective (algebraMap ℤ (𝓞 ℚ)) := by
  have h : (algebraMap ℤ (𝓞 ℚ)) = ((Rat.IsIntegralClosure.intEquiv (𝓞 ℚ)).symm : ℤ →+* 𝓞 ℚ) :=
    RingHom.ext_int _ _
  rw [h]
  exact (Rat.IsIntegralClosure.intEquiv (𝓞 ℚ)).symm.surjective

/-- ★★**高さ 1 素イデアルを `ℤ` へ引き戻すと `(p)`** —— `p` は `natGenerator`。 -/
theorem comap_asIdeal (𝔭 : HeightOneSpectrum (𝓞 ℚ)) :
    Ideal.comap (algebraMap ℤ (𝓞 ℚ)) 𝔭.asIdeal
      = Ideal.span {((Rat.HeightOneSpectrum.natGenerator 𝔭 : ℕ) : ℤ)} := by
  have h1 : Ideal.comap (algebraMap ℤ (𝓞 ℚ)) 𝔭.asIdeal
      = Ideal.comap ((Rat.IsIntegralClosure.intEquiv (𝓞 ℚ)).symm : ℤ →+* 𝓞 ℚ) 𝔭.asIdeal := by
    congr 1
    exact RingHom.ext_int _ _
  rw [h1, Rat.HeightOneSpectrum.span_natGenerator]
  exact Ideal.comap_symm _

/-- ★★**`N𝔭 = p`** —— 剰余環を `ℤ ⧸ (p) ≃+* ZMod p` へ移して数える。 -/
theorem absNorm_asIdeal (𝔭 : HeightOneSpectrum (𝓞 ℚ)) :
    Ideal.absNorm 𝔭.asIdeal = Rat.HeightOneSpectrum.natGenerator 𝔭 := by
  set p : ℕ := Rat.HeightOneSpectrum.natGenerator 𝔭 with hp
  have hmap : Ideal.span {((p : ℕ) : ℤ)}
      = Ideal.map ((Rat.IsIntegralClosure.intEquiv (𝓞 ℚ)) : 𝓞 ℚ →+* ℤ) 𝔭.asIdeal := by
    rw [Rat.HeightOneSpectrum.span_natGenerator]
    rfl
  have he : (𝓞 ℚ ⧸ 𝔭.asIdeal) ≃+* (ℤ ⧸ Ideal.span {((p : ℕ) : ℤ)}) :=
    Ideal.quotientEquiv _ _ (Rat.IsIntegralClosure.intEquiv (𝓞 ℚ)) hmap
  have h1 : Ideal.absNorm 𝔭.asIdeal = Nat.card (𝓞 ℚ ⧸ 𝔭.asIdeal) := by
    rw [Ideal.absNorm_apply, Submodule.cardQuot_apply]
  have h2 : Nat.card (ℤ ⧸ Ideal.span {((p : ℕ) : ℤ)}) = p := by
    rw [Nat.card_congr (Int.quotientSpanEquivZMod ((p : ℕ) : ℤ)).toEquiv, Nat.card_zmod]
    simp
  rw [h1, Nat.card_congr he.toEquiv, h2]

/-- ★`primesEquiv` の値は `natGenerator` そのもの。 -/
theorem coe_primesEquiv (𝔭 : HeightOneSpectrum (𝓞 ℚ)) :
    ((Rat.HeightOneSpectrum.primesEquiv 𝔭 : Nat.Primes) : ℕ)
      = Rat.HeightOneSpectrum.natGenerator 𝔭 := by
  simp [Rat.HeightOneSpectrum.primesEquiv]

variable (L : Type*) [Field L] [NumberField L]

/-- ★★**`𝔭` の上にある素イデアルと `(p)` の上にあるものは同じ**。

★`ℤ → 𝓞 ℚ → 𝓞 L` の塔(`Ideal.under_under`)と、`ℤ ≃+* 𝓞 ℚ` の全射性だけを使う。 -/
theorem primesOver_asIdeal_eq (𝔭 : HeightOneSpectrum (𝓞 ℚ)) :
    Ideal.primesOver 𝔭.asIdeal (𝓞 L)
      = Ideal.primesOver
          (Ideal.span {((Rat.HeightOneSpectrum.natGenerator 𝔭 : ℕ) : ℤ)}) (𝓞 L) := by
  ext P
  simp only [Ideal.primesOver, Set.mem_setOf_eq]
  constructor
  · rintro ⟨hP, hLO⟩
    refine ⟨hP, ⟨?_⟩⟩
    rw [← comap_asIdeal 𝔭, ← Ideal.under_under (A := ℤ) (B := 𝓞 ℚ) P, Ideal.under_def,
      ← hLO.over]
  · rintro ⟨hP, hLO⟩
    refine ⟨hP, ⟨?_⟩⟩
    have h2 : Ideal.comap (algebraMap ℤ (𝓞 ℚ)) 𝔭.asIdeal
        = Ideal.comap (algebraMap ℤ (𝓞 ℚ)) (Ideal.under (𝓞 ℚ) P) := by
      rw [comap_asIdeal, ← Ideal.under_def, Ideal.under_under]
      exact hLO.over
    exact Ideal.comap_injective_of_surjective _ algebraMap_int_ringOfIntegersRat_surjective h2

/-- ★★**完全分解の 2 つの書き方は一致する** —— `HeightOneSpectrum (𝓞 ℚ)` の側と
`Nat.Primes` の側(`SplQ`)。 -/
theorem splitsCompletely_heightOne_iff_mem_SplQ (𝔭 : HeightOneSpectrum (𝓞 ℚ)) :
    (Ideal.primesOver 𝔭.asIdeal (𝓞 L)).ncard = Module.finrank ℚ L
      ↔ Rat.HeightOneSpectrum.primesEquiv 𝔭 ∈ SplQ L := by
  rw [primesOver_asIdeal_eq L 𝔭]
  show _ ↔ (Ideal.primesOver (Ideal.span
    {((((Rat.HeightOneSpectrum.primesEquiv 𝔭 : Nat.Primes) : ℕ)) : ℤ)}) (𝓞 L)).ncard = _
  rw [coe_primesEquiv]

/-! ## ★2. Chebotarev(完全分解の場合)を `HeightOneSpectrum` の語彙で -/

/-- ★★★★★★**[cheb-split-density] 完全分解する素点の Dirichlet 密度は `1/[L:ℚ]`**
——`HeightOneSpectrum (𝓞 ℚ)` の語彙で述べたもの。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★中身は `tendsto_splQ_div_log`(**類体論を使わない**)の添字の張り替えだけである。 -/
theorem hasDirichletDensity_splitsCompletely [IsGalois ℚ L] :
    HasDirichletDensity (K := ℚ)
        {𝔭 : HeightOneSpectrum (𝓞 ℚ) |
          (Ideal.primesOver 𝔭.asIdeal (𝓞 L)).ncard = Module.finrank ℚ L}
        (1 / (Module.finrank ℚ L : ℝ)) := by
  refine Filter.Tendsto.congr (fun s => ?_) (tendsto_splQ_div_log (L := L))
  congr 1
  rw [← Equiv.tsum_eq (Equiv.subtypeEquiv Rat.HeightOneSpectrum.primesEquiv
      (fun 𝔭 => splitsCompletely_heightOne_iff_mem_SplQ L 𝔭))
      (fun p : SplQ L => (((p : Nat.Primes) : ℕ) : ℝ) ^ (-s))]
  refine tsum_congr fun 𝔭 => ?_
  simp only [Equiv.subtypeEquiv_apply]
  rw [coe_primesEquiv, absNorm_asIdeal]

/-- ★★`SplQ L` は無限集合(★Chebotarev を使わない迂回路
`infinite_splitsCompletely_of_isGalois`)。 -/
theorem infinite_SplQ [IsGalois ℚ L] : (SplQ L).Infinite := by
  refine Set.Infinite.of_image (fun p : Nat.Primes => (p : ℕ)) ?_
  have hset : (fun p : Nat.Primes => (p : ℕ)) '' (SplQ L)
      = {p : ℕ | p.Prime ∧
          (Ideal.primesOver (Ideal.span {(p : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L} := by
    ext n
    constructor
    · rintro ⟨p, hp, rfl⟩
      exact ⟨p.2, hp⟩
    · rintro ⟨hn, h⟩
      exact ⟨⟨n, hn⟩, h, rfl⟩
  rw [hset]
  exact infinite_splitsCompletely_of_isGalois L

/-- ★★★**完全分解する素点は無限個**(`HeightOneSpectrum (𝓞 ℚ)` の語彙)。

★★これは密度を要求しない —— 最小多項式の値の素因子を数える初等的な議論
(`infinite_splitsCompletely_of_isGalois`、**仮定なし**)で出る。 -/
theorem infinite_splitsCompletely_heightOne [IsGalois ℚ L] :
    {𝔭 : HeightOneSpectrum (𝓞 ℚ) |
      (Ideal.primesOver 𝔭.asIdeal (𝓞 L)).ncard = Module.finrank ℚ L}.Infinite := by
  have hset : {𝔭 : HeightOneSpectrum (𝓞 ℚ) |
      (Ideal.primesOver 𝔭.asIdeal (𝓞 L)).ncard = Module.finrank ℚ L}
      = Rat.HeightOneSpectrum.primesEquiv ⁻¹' (SplQ L) := by
    ext 𝔭
    exact splitsCompletely_heightOne_iff_mem_SplQ L 𝔭
  rw [hset]
  exact (infinite_SplQ L).preimage (by rw [Equiv.range_eq_univ]; exact Set.subset_univ _)

/-! ## ★3. 完全分解する素点の包含が体の埋め込みを与える -/

/-- ★★★★★★★**`Spl(L₁) ⊆ Spl(L₂)` なら `L₂ ↪ L₁`**(抽象な 2 つの数体の版)。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★中身は「両方を `ℚ̄` へ埋め込んで合成体の中で `le_of_SplQ_subset` を当てる」だけ
(`nonempty_algEquiv_of_SplQ_eq` と同じ手順の、等号を包含に緩めた版)。 -/
theorem nonempty_algHom_of_SplQ_subset (L₁ L₂ : Type) [Field L₁] [Field L₂]
    [NumberField L₁] [NumberField L₂] [IsGalois ℚ L₁] [IsGalois ℚ L₂]
    (h : SplQ L₁ ⊆ SplQ L₂) : Nonempty (L₂ →ₐ[ℚ] L₁) := by
  classical
  let i₁ : L₁ →ₐ[ℚ] AlgebraicClosure ℚ := IsAlgClosed.lift
  let i₂ : L₂ →ₐ[ℚ] AlgebraicClosure ℚ := IsAlgClosed.lift
  let F₁ : IntermediateField ℚ (AlgebraicClosure ℚ) := i₁.fieldRange
  let F₂ : IntermediateField ℚ (AlgebraicClosure ℚ) := i₂.fieldRange
  let e₁ : L₁ ≃ₐ[ℚ] ↥F₁ := AlgHom.equivFieldRange
  let e₂ : L₂ ≃ₐ[ℚ] ↥F₂ := AlgHom.equivFieldRange
  haveI hf1 : FiniteDimensional ℚ ↥F₁ := LinearEquiv.finiteDimensional e₁.toLinearEquiv
  haveI hf2 : FiniteDimensional ℚ ↥F₂ := LinearEquiv.finiteDimensional e₂.toLinearEquiv
  haveI : NumberField ↥F₁ := ⟨⟩
  haveI : NumberField ↥F₂ := ⟨⟩
  haveI hg1 : IsGalois ℚ ↥F₁ := (AlgEquiv.transfer_galois e₁).mp ‹_›
  haveI hg2 : IsGalois ℚ ↥F₂ := (AlgEquiv.transfer_galois e₂).mp ‹_›
  haveI hfsup : FiniteDimensional ℚ ↥(F₁ ⊔ F₂) :=
    IntermediateField.finiteDimensional_sup F₁ F₂
  haveI : NumberField ↥(F₁ ⊔ F₂) := ⟨⟩
  haveI : IsGalois ℚ ↥(F₁ ⊔ F₂) := isGalois_sup_ac F₁ F₂ hg1.to_normal hg2.to_normal
  let G₁ : IntermediateField ℚ ↥(F₁ ⊔ F₂) :=
    (IntermediateField.inclusion (le_sup_left : F₁ ≤ F₁ ⊔ F₂)).fieldRange
  let G₂ : IntermediateField ℚ ↥(F₁ ⊔ F₂) :=
    (IntermediateField.inclusion (le_sup_right : F₂ ≤ F₁ ⊔ F₂)).fieldRange
  let f₁ : ↥F₁ ≃ₐ[ℚ] ↥G₁ := AlgHom.equivFieldRange
  let f₂ : ↥F₂ ≃ₐ[ℚ] ↥G₂ := AlgHom.equivFieldRange
  haveI : FiniteDimensional ℚ ↥G₁ := LinearEquiv.finiteDimensional f₁.toLinearEquiv
  haveI : FiniteDimensional ℚ ↥G₂ := LinearEquiv.finiteDimensional f₂.toLinearEquiv
  haveI : NumberField ↥G₁ := ⟨⟩
  haveI : NumberField ↥G₂ := ⟨⟩
  haveI hgg1 : IsGalois ℚ ↥G₁ := (AlgEquiv.transfer_galois f₁).mp hg1
  haveI hgg2 : IsGalois ℚ ↥G₂ := (AlgEquiv.transfer_galois f₂).mp hg2
  have hsub : SplQ ↥G₁ ⊆ SplQ ↥G₂ := by
    rw [← SplQ_congr f₁, ← SplQ_congr e₁, ← SplQ_congr f₂, ← SplQ_congr e₂]
    exact h
  have hle : G₂ ≤ G₁ := le_of_SplQ_subset G₂ G₁ hgg2 hgg1 hsub
  exact ⟨((e₁.trans f₁).symm.toAlgHom).comp
    ((IntermediateField.inclusion hle).comp ((e₂.trans f₂).toAlgHom))⟩

/-- ★★★★**完全分解する素点の包含が体の埋め込みを与える**
(`HeightOneSpectrum (𝓞 ℚ)` の語彙)。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★`Theorem 6.4, (iv)` が最後に使うのはこれである。 -/
theorem nonempty_algHom_of_splitsCompletely_subset (L₁ L₂ : Type) [Field L₁] [Field L₂]
    [NumberField L₁] [NumberField L₂] [IsGalois ℚ L₁] [IsGalois ℚ L₂]
    (h : {𝔭 : HeightOneSpectrum (𝓞 ℚ) |
            (Ideal.primesOver 𝔭.asIdeal (𝓞 L₁)).ncard = Module.finrank ℚ L₁}
          ⊆ {𝔭 : HeightOneSpectrum (𝓞 ℚ) |
            (Ideal.primesOver 𝔭.asIdeal (𝓞 L₂)).ncard = Module.finrank ℚ L₂}) :
    Nonempty (L₂ →ₐ[ℚ] L₁) := by
  refine nonempty_algHom_of_SplQ_subset L₁ L₂ (fun p hp => ?_)
  obtain ⟨𝔭, h1⟩ : ∃ 𝔭 : HeightOneSpectrum (𝓞 ℚ),
      Rat.HeightOneSpectrum.primesEquiv 𝔭 = p :=
    ⟨Rat.HeightOneSpectrum.primesEquiv.symm p, Equiv.apply_symm_apply _ _⟩
  have hL : (Ideal.primesOver 𝔭.asIdeal (𝓞 L₁)).ncard = Module.finrank ℚ L₁ :=
    (splitsCompletely_heightOne_iff_mem_SplQ L₁ 𝔭).mpr (by rw [h1]; exact hp)
  have hM := (splitsCompletely_heightOne_iff_mem_SplQ L₂ 𝔭).mp (h hL)
  rwa [h1] at hM

/-! ### ★出典の紐付け -/

def hasDirichletDensity_splitsCompletely.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素点の Dirichlet 密度(高さ 1 素イデアルの語彙)",
    sectionId := "frdi-thm-6-4" }

def hasDirichletDensity_splitsCompletely.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tendsto_splQ_div_log(類体論を使わない Chebotarev、完全分解の場合)"
      (.inProject "ABC3" "ABC3.Found.NF.tendsto_splQ_div_log") 116,
    .citation "[mathlib]" "Rat.HeightOneSpectrum.primesEquiv(高さ 1 素イデアル ≃ 有理素数)"
      (.inMathlib "Rat.HeightOneSpectrum.primesEquiv") 116,
    .derivation "添字の張り替え(N𝔭 = p、𝔭 の上の素イデアル = (p) の上の素イデアル)" 116 ]

def infinite_splitsCompletely_heightOne.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素点は無限個(高さ 1 素イデアルの語彙)",
    sectionId := "frdi-thm-6-4" }

def nonempty_algHom_of_splitsCompletely_subset.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素点の包含が体の埋め込みを与える",
    sectionId := "frdi-thm-6-4" }

def nonempty_algHom_of_splitsCompletely_subset.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "le_of_SplQ_subset(Bauer——Spl の包含から体の包含)"
      (.inProject "ABC3" "ABC3.Found.NF.le_of_SplQ_subset") 116,
    .citation "[ABC3]" "nonempty_algEquiv_of_SplQ_eq(同じ埋め込みの手順の等号版)"
      (.inProject "ABC3" "ABC3.Found.NF.nonempty_algEquiv_of_SplQ_eq") 116,
    .derivation "両方を ℚ̄ へ埋め込み、合成体の中で le_of_SplQ_subset を当てて包含を作る" 116 ]

end ABC3.Found.NF
