import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.InfinitePlace.Ramification
import Mathlib.FieldTheory.PrimitiveElement

/-!
# [GenEll] §1 —— 無限素点の**相対**次数和(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

## ★★これは mathlib に**無い**(2026-08-17 実測)

`NumberField.InfinitePlace` には**絶対版** `sum_mult_eq : Σ_w mult w = finrank ℚ K` がある。
★しかし**相対版**——`F ⊆ K` について

> **`Σ_{w | v} mult w = mult v · [K:F]`**

——は無い(`Ramification.lean` の全定理名を確認。`comap` はあるが次数和は無い)。

## ★なぜ要るか —— `deg_F` の正規化が底変換で不変であること

`degNormalized(a) ≝ deg_F(a)/[F:ℚ]` が底変換で不変であることが、
**`ht` を `X(ℚ̄)` の上で well-defined にする**(`Interface/GenEll/ArithLineBundle.lean` の
`base_change_invariant` が posit しているもの)。

★アルキメデス側でそれを言うのに、ちょうどこの相対次数和が要る。
有限側は mathlib の `Ideal.sum_ramification_inertia`(`Σ_{w|v} e_w f_w = [K:F]`)で足りる。

## ★機構 —— 埋め込みの二重数え上げ

`S ≝ {φ : K →+* ℂ | φ|_F が定める素点が v}` を 2 通りに数える:

- `w` ごとに分けると `Σ_{w|v} #{φ | mk φ = w} = Σ_{w|v} mult w`
  (`NumberField.InfinitePlace.card_filter_mk_eq`)
- `ψ` ごとに分けると `Σ_{ψ : mk ψ = v} #{φ | φ|_F = ψ} = mult v · [K:F]`
  (`AlgHom.card`——代数閉体への延長は `[K:F]` 個)
-/

namespace ABC3.Found.GenEll

open NumberField

variable (F K : Type*) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

/-- ★`ψ : F →+* ℂ` を延長する `φ : K →+* ℂ` はちょうど `[K:F]` 個ある。

★`ℂ` を `ψ` によって `F`-代数と見なすと、延長は `K →ₐ[F] ℂ` と 1 対 1 で、
その個数は `AlgHom.card`(代数閉体なので `finrank`)である。 -/
theorem card_extension_eq (ψ : F →+* ℂ) :
    Nat.card {φ : K →+* ℂ // φ.comp (algebraMap F K) = ψ} = Module.finrank F K := by
  classical
  letI : Algebra F ℂ := ψ.toAlgebra
  haveI : FiniteDimensional F K := Module.Finite.of_restrictScalars_finite ℚ F K
  haveI : Algebra.IsSeparable F K := Algebra.IsSeparable.of_integral F K
  -- 延長の集合は `K →ₐ[F] ℂ` と同型
  let e : {φ : K →+* ℂ // φ.comp (algebraMap F K) = ψ} ≃ (K →ₐ[F] ℂ) :=
    { toFun := fun p =>
        { toRingHom := p.1
          commutes' := fun c => by
            have := congrArg (fun r : F →+* ℂ => r c) p.2
            simpa [RingHom.algebraMap_toAlgebra] using this }
      invFun := fun f => ⟨f.toRingHom, by
        ext c
        simpa [RingHom.algebraMap_toAlgebra] using f.commutes c⟩
      left_inv := fun p => by ext x; rfl
      right_inv := fun f => by ext x; rfl }
  rw [Nat.card_congr e, Nat.card_eq_fintype_card, AlgHom.card F K ℂ]

open scoped Classical in
/-- ★★**無限素点の相対次数和** —— `Σ_{w | v} mult w = mult v · [K:F]`。

★mathlib は**絶対版**(`sum_mult_eq : Σ_w mult w = finrank ℚ K`)しか持たない
(2026-08-17 実測)。ここで相対版を取る。

★証明は `S ≝ {φ : K →+* ℂ | φ|_F が定める素点が v}` の**二重数え上げ**である。 -/
theorem sum_mult_comap_eq (v : InfinitePlace F) :
    ∑ w ∈ Finset.univ.filter (fun w : InfinitePlace K =>
        w.comap (algebraMap F K) = v), InfinitePlace.mult w
      = InfinitePlace.mult v * Module.finrank F K := by
  classical
  set S : Finset (K →+* ℂ) :=
    Finset.univ.filter (fun φ : K →+* ℂ =>
      (InfinitePlace.mk φ).comap (algebraMap F K) = v) with hS
  have hmemS : ∀ φ : K →+* ℂ, φ ∈ S ↔ (InfinitePlace.mk φ).comap (algebraMap F K) = v := by
    intro φ; simp [hS]
  -- ★数え方 1 —— `w` ごとに分ける
  have hmem1 : ∀ φ ∈ S, InfinitePlace.mk φ ∈
      Finset.univ.filter (fun w : InfinitePlace K => w.comap (algebraMap F K) = v) := by
    intro φ hφ
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact (hmemS φ).1 hφ
  have hcount1 : S.card
      = ∑ w ∈ Finset.univ.filter (fun w : InfinitePlace K =>
          w.comap (algebraMap F K) = v), InfinitePlace.mult w := by
    rw [Finset.card_eq_sum_card_fiberwise hmem1]
    refine Finset.sum_congr rfl fun w hw => ?_
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hw
    rw [← InfinitePlace.card_filter_mk_eq w]
    congr 1
    ext φ
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨_, h⟩; exact h
    · intro h
      refine ⟨(hmemS φ).2 ?_, h⟩
      rw [h]; exact hw
  -- ★数え方 2 —— `ψ` ごとに分ける
  have hmem2 : ∀ φ ∈ S, φ.comp (algebraMap F K) ∈
      Finset.univ.filter (fun ψ : F →+* ℂ => InfinitePlace.mk ψ = v) := by
    intro φ hφ
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    have := (hmemS φ).1 hφ
    rwa [InfinitePlace.comap_mk] at this
  have hcount2 : S.card = InfinitePlace.mult v * Module.finrank F K := by
    rw [Finset.card_eq_sum_card_fiberwise hmem2]
    have hfib : ∀ ψ ∈ Finset.univ.filter (fun ψ : F →+* ℂ => InfinitePlace.mk ψ = v),
        (S.filter (fun φ : K →+* ℂ => φ.comp (algebraMap F K) = ψ)).card
          = Module.finrank F K := by
      intro ψ hψ
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hψ
      have hset : S.filter (fun φ : K →+* ℂ => φ.comp (algebraMap F K) = ψ)
          = Finset.univ.filter (fun φ : K →+* ℂ => φ.comp (algebraMap F K) = ψ) := by
        ext φ
        simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        constructor
        · rintro ⟨_, h⟩; exact h
        · intro h
          refine ⟨(hmemS φ).2 ?_, h⟩
          rw [InfinitePlace.comap_mk, h, hψ]
      rw [hset, ← card_extension_eq F K ψ, Nat.card_eq_fintype_card, Fintype.card_subtype]
    rw [Finset.sum_congr rfl hfib, Finset.sum_const, smul_eq_mul,
      InfinitePlace.card_filter_mk_eq v]
  rw [← hcount1, hcount2]

/-! ## ★出典の紐付け(`.src`) -/

def sum_mult_comap_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

end ABC3.Found.GenEll
