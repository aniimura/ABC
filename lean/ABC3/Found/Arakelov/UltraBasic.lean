import ABC3.Found.Arakelov.UltraTransfer

/-!
# Arakelov (C2) 第 85–86 ブロック —— **★★★基本開集合と閉点**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★極限点は chart に丸ごと収まる

持ち上げ `l : Spec O ⟶ X` の **`O` は付値環——すなわち局所環**である。
★局所環の `Spec` では、**閉点を含む開集合は全体しかない**。

★★したがって極限点(閉点の像)が `V` に入れば、
**持ち上げ全体が `V` に入る**——一般点も含めて。

★★★閉点が `st` の核であることは、`st` が体への全射だから
(核は極大イデアル、局所環の極大イデアルは一意)。

## ★★★★基本開集合で挟む

超積の点は `U ∩ V` に入る。★そこを**基本開集合 `D(f)` で挟む**と、
「`f` の値が 0 でない」という**評価だけの条件**になり、
★★超フィルターの言葉(ほとんど至る所 0 でない)に翻訳できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mem_basicOpen_comap_field` | ★基本開集合 ⟺ 値が 0 でない |
| `exists_basicOpen_le_of_mem` | ★開集合の中の点は `D(f)` で挟める |
| `range_comap_subset_basicOpen` | ★★`φ f` が単元 ⟹ 像は `D(f)` |
| `stHom_surjective` / `comap_stHom_eq_closedPoint` | ★★`st` の核は閉点 |
| `range_ultraLift_subset` | ★★★**持ち上げは丸ごと chart に入る** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Filter

/-! ## ★基本開集合 -/

/-- ★体に値を取る点が基本開集合に入ることは、値が 0 でないことと同値。 -/
theorem mem_basicOpen_comap_field {A : Type} [CommRing A] {K : Type} [Field K] (ψ : A →+* K)
    (f : A) (P : PrimeSpectrum K) :
    (PrimeSpectrum.comap ψ P) ∈ PrimeSpectrum.basicOpen f ↔ ψ f ≠ 0 := by
  have hP : P.asIdeal = ⊥ := by
    rcases P.asIdeal.eq_bot_or_top with h1 | h1
    · exact h1
    · exact absurd h1 P.isPrime.ne_top
  rw [PrimeSpectrum.mem_basicOpen]
  simp [PrimeSpectrum.comap_asIdeal, hP]

/-- ★開集合の中の点は基本開集合で挟める。 -/
theorem exists_basicOpen_le_of_mem {A : Type} [CommRing A] (W : Set (PrimeSpectrum A))
    (hW : IsOpen W) (Q : PrimeSpectrum A) (hQ : Q ∈ W) :
    ∃ f : A, Q ∈ PrimeSpectrum.basicOpen f ∧
      (PrimeSpectrum.basicOpen f : Set (PrimeSpectrum A)) ⊆ W := by
  obtain ⟨_, ⟨_, ⟨g, rfl⟩, rfl⟩, h1, h2⟩ :=
    (PrimeSpectrum.isBasis_basic_opens).exists_subset_of_mem_open hQ hW
  exact ⟨g, h1, h2⟩

/-- ★★`φ f` が単元なら、像は基本開集合に入る。 -/
theorem range_comap_subset_basicOpen {A R : Type} [CommRing A] [CommRing R] (φ : A →+* R)
    (f : A) (hf : IsUnit (φ f)) :
    Set.range (PrimeSpectrum.comap φ) ⊆ (PrimeSpectrum.basicOpen f : Set (PrimeSpectrum A)) := by
  rintro _ ⟨P, rfl⟩
  show PrimeSpectrum.comap φ P ∈ PrimeSpectrum.basicOpen f
  rw [PrimeSpectrum.mem_basicOpen]
  intro hmem
  rw [PrimeSpectrum.comap_asIdeal, Ideal.mem_comap] at hmem
  exact P.isPrime.ne_top (Ideal.eq_top_of_isUnit_mem _ hmem hf)

/-! ## ★★閉点 -/

variable {α : Type} (𝒰 : Ultrafilter α)

theorem stHom_surjective : Function.Surjective (stHom 𝒰) := by
  intro c
  refine ⟨⟨(((fun _ => c) : α → ℂ) : Germ (𝒰 : Filter α) ℂ), ⟨‖c‖, ?_⟩⟩, ?_⟩
  · rw [Germ.liftPred_coe]
    filter_upwards with a
    exact le_rfl
  · exact stG_const 𝒰 c

/-- ★★`st` の核は極大——したがって閉点である。 -/
theorem comap_stHom_eq_closedPoint :
    PrimeSpectrum.comap (stHom 𝒰) (default : PrimeSpectrum ℂ)
      = IsLocalRing.closedPoint ↥(finGermSub 𝒰) := by
  have hmax : (RingHom.ker (stHom 𝒰)).IsMaximal :=
    RingHom.ker_isMaximal_of_surjective _ (stHom_surjective 𝒰)
  have hP : (default : PrimeSpectrum ℂ).asIdeal = ⊥ := by
    rcases (default : PrimeSpectrum ℂ).asIdeal.eq_bot_or_top with h1 | h1
    · exact h1
    · exact absurd h1 (default : PrimeSpectrum ℂ).isPrime.ne_top
  have hker : (PrimeSpectrum.comap (stHom 𝒰) (default : PrimeSpectrum ℂ)).asIdeal
      = RingHom.ker (stHom 𝒰) := by
    rw [PrimeSpectrum.comap_asIdeal, hP]
    rfl
  apply PrimeSpectrum.ext
  rw [hker]
  exact (IsLocalRing.eq_maximalIdeal hmax).symm ▸ rfl

/-! ## ★★★持ち上げは chart に収まる -/

/-- ★★局所環の Spec からの射は、閉点の像を含む開集合に丸ごと入る。 -/
theorem range_base_subset_of_closedPoint {X : Scheme.{0}} {R : Type} [CommRing R] [IsLocalRing R]
    (l : Spec (CommRingCat.of R) ⟶ X) (V : X.Opens)
    (h : l.base (IsLocalRing.closedPoint R) ∈ V) :
    Set.range l.base ⊆ (V : Set X) := by
  have hopen : IsOpen (l.base ⁻¹' (V : Set X)) := (V.2).preimage (Scheme.Hom.continuous l)
  have htop : (⟨l.base ⁻¹' (V : Set X), hopen⟩ : TopologicalSpace.Opens (PrimeSpectrum R)) = ⊤ :=
    (IsLocalRing.closedPoint_mem_iff _).1 h
  rintro _ ⟨x, rfl⟩
  have hx : x ∈ (⟨l.base ⁻¹' (V : Set X), hopen⟩ : TopologicalSpace.Opens (PrimeSpectrum R)) := by
    rw [htop]; trivial
  exact hx

/-- ★★★**極限点が `V` に入るなら、持ち上げは丸ごと `V` に入る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
theorem range_ultraLift_subset {X : Scheme.{0}}
    (l : Spec (CommRingCat.of ↥(finGermSub 𝒰)) ⟶ X) (V : X.Opens)
    (hq : (Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l).base
        (default : ↥(Spec (CommRingCat.of ℂ))) ∈ V) :
    Set.range l.base ⊆ (V : Set X) := by
  refine range_base_subset_of_closedPoint l V ?_
  have hd : (default : PrimeSpectrum ℂ) = (default : ↥(Spec (CommRingCat.of ℂ))) :=
    Subsingleton.elim _ _
  rw [← comap_stHom_eq_closedPoint 𝒰, hd]
  exact hq

/-! ## ★出典の紐付け(`.src`) -/

def range_ultraLift_subset.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——局所環の閉点から chart への収まり)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
