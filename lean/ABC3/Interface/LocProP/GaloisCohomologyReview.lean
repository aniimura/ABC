import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs

/-!
# [LocProP] §2 —— Review of Galois Cohomology の posit(`Interface`)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、物理 p.20-21。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

`Lemma 2.1`-`2.3`・`Proposition 2.4` は Faltings の almost étale extension 理論
([Falt1]・[Falt2])・Hodge-Tate 分解・log 構造つきスキームの Galois cohomology を扱う。
★★**mathlib にはこの理論が丸ごと無い**(almost ring theory・perfectoid・Hodge-Tate
分解いずれも 0 件、2026-09-04 実測)ので、named な群とその間の同型/完全列として posit する
(`Definition 0.3`・`Lemma 0.4`・§1 と同じ扱い)。
-/

namespace ABC3.Interface.LocProP

universe u

/-- 加法版の自明群(non-vacuous witness 用)。 -/
instance addTrivialGroup : AddCommGroup Unit where
  add _ _ := ()
  zero := ()
  neg _ := ()
  nsmul _ _ := ()
  zsmul _ _ := ()
  add_assoc := by intros; rfl
  zero_add := by intros; rfl
  add_zero := by intros; rfl
  neg_add_cancel := by intros; rfl
  add_comm := by intros; rfl

/-- ★posit —— §2(Review of Galois Cohomology)の骨組み。 -/
structure GaloisCohomologyReviewSetup where
  /-- **[LocProP] Lemma 2.1, (1)**: `H⁰(Δ_R^prf, R̂_K̄(j)) = R̂_K̄(j)`。 -/
  H0DeltaRprf : ℤ → Type u
  H0DeltaRprfGrp : ∀ j, AddCommGroup (H0DeltaRprf j)
  RhatKbar : ℤ → Type u
  RhatKbarGrp : ∀ j, AddCommGroup (RhatKbar j)
  lemma_2_1_i : ∀ j, H0DeltaRprf j ≃+ RhatKbar j
  /-- **[LocProP] Lemma 2.1, (2)**: `H¹(Δ_R^prf, R̂_K̄(j)) = Ω_{S^log/O_K} ⊗_R R̂_K̄(j-1)`。 -/
  H1DeltaRprf : ℤ → Type u
  H1DeltaRprfGrp : ∀ j, AddCommGroup (H1DeltaRprf j)
  OmegaTensorRhatKbar : ℤ → Type u
  OmegaTensorRhatKbarGrp : ∀ j, AddCommGroup (OmegaTensorRhatKbar j)
  lemma_2_1_ii : ∀ j, H1DeltaRprf j ≃+ OmegaTensorRhatKbar j
  /-- **[LocProP] Lemma 2.2, (i)(iii)**: `H⁰(Γ_K, R̂_K̄) = R̂_K = H¹(Γ_K, R̂_K̄)`。 -/
  H0GammaK : Type u
  H0GammaKGrp : AddCommGroup H0GammaK
  RhatK : Type u
  RhatKGrp : AddCommGroup RhatK
  H1GammaK : Type u
  H1GammaKGrp : AddCommGroup H1GammaK
  lemma_2_2_i : H0GammaK ≃+ RhatK
  lemma_2_2_iii : H1GammaK ≃+ RhatK
  /-- **[LocProP] Lemma 2.2, (ii)**: `Hⁿ(Γ_K, R̂_K̄(j)) = 0`(`j ≠ 0`、`n = 0, 1`)。 -/
  HnGammaKTwisted : ℤ → Fin 2 → Type u
  HnGammaKTwistedGrp : ∀ j n, AddCommGroup (HnGammaKTwisted j n)
  lemma_2_2_ii : ∀ j ≠ (0 : ℤ), ∀ n, Subsingleton (HnGammaKTwisted j n)
  /-- **[LocProP] Lemma 2.3, (i)**: `Hⁿ(Γ_R, R̂̂_K) = R̂̂_K`(`n = 0, 1`)。 -/
  HnGammaR : Fin 2 → Type u
  HnGammaRGrp : ∀ n, AddCommGroup (HnGammaR n)
  RhatKHat : Type u
  RhatKHatGrp : AddCommGroup RhatKHat
  lemma_2_3_i : ∀ n, HnGammaR n ≃+ RhatKHat
  /-- **[LocProP] Lemma 2.3, (ii)**: `Hⁿ(Γ_R, R̂̂_K(-1)) = 0`。 -/
  HnGammaRNeg1 : Fin 2 → Type u
  HnGammaRNeg1Grp : ∀ n, AddCommGroup (HnGammaRNeg1 n)
  lemma_2_3_ii : ∀ n, Subsingleton (HnGammaRNeg1 n)
  /-- **[LocProP] Proposition 2.4**: `0 → ℋ₁ ⊗ R̂_K(1) → H ⊗ R̂_K → ℋ₀ ⊗ R̂_K → 0` 完全。 -/
  H1Otimes : Type u
  H1OtimesGrp : AddCommGroup H1Otimes
  HOtimesRhatK : Type u
  HOtimesRhatKGrp : AddCommGroup HOtimesRhatK
  H0Otimes : Type u
  H0OtimesGrp : AddCommGroup H0Otimes
  prop24Left : H1Otimes →+ HOtimesRhatK
  prop24Right : HOtimesRhatK →+ H0Otimes
  prop24_left_injective : Function.Injective prop24Left
  prop24_right_surjective : Function.Surjective prop24Right
  prop24_exact : ∀ y, prop24Right y = 0 ↔ ∃ x, prop24Left x = y

/-- ★★**非退化性の witness**(具体項)—— 同型を要求する箇所は `ℤ ≃+ ℤ`(自明でない)、
`= 0` を要求する箇所は `Unit`(自明群)に潰す。 -/
@[reducible] def GaloisCohomologyReviewSetup.example : GaloisCohomologyReviewSetup.{0} where
  H0DeltaRprf := fun _ => ℤ; H0DeltaRprfGrp := fun _ => inferInstance
  RhatKbar := fun _ => ℤ; RhatKbarGrp := fun _ => inferInstance
  lemma_2_1_i := fun _ => AddEquiv.refl ℤ
  H1DeltaRprf := fun _ => ℤ; H1DeltaRprfGrp := fun _ => inferInstance
  OmegaTensorRhatKbar := fun _ => ℤ; OmegaTensorRhatKbarGrp := fun _ => inferInstance
  lemma_2_1_ii := fun _ => AddEquiv.refl ℤ
  H0GammaK := ℤ; H0GammaKGrp := inferInstance
  RhatK := ℤ; RhatKGrp := inferInstance
  H1GammaK := ℤ; H1GammaKGrp := inferInstance
  lemma_2_2_i := AddEquiv.refl ℤ
  lemma_2_2_iii := AddEquiv.refl ℤ
  HnGammaKTwisted := fun _ _ => Unit; HnGammaKTwistedGrp := fun _ _ => inferInstance
  lemma_2_2_ii := fun _ _ _ => inferInstance
  HnGammaR := fun _ => ℤ; HnGammaRGrp := fun _ => inferInstance
  RhatKHat := ℤ; RhatKHatGrp := inferInstance
  lemma_2_3_i := fun _ => AddEquiv.refl ℤ
  HnGammaRNeg1 := fun _ => Unit; HnGammaRNeg1Grp := fun _ => inferInstance
  lemma_2_3_ii := fun _ => inferInstance
  H1Otimes := ℤ; H1OtimesGrp := inferInstance
  HOtimesRhatK := ℤ; HOtimesRhatKGrp := inferInstance
  H0Otimes := Unit; H0OtimesGrp := inferInstance
  prop24Left := AddMonoidHom.id ℤ
  prop24Right := 0
  prop24_left_injective := fun _ _ h => h
  prop24_right_surjective := fun _ => ⟨0, Subsingleton.elim _ _⟩
  prop24_exact := fun y => ⟨fun _ => ⟨y, rfl⟩, fun _ => Subsingleton.elim _ _⟩

@[reducible] def GaloisCohomologyReviewSetup.nonvacuous : Nonempty (GaloisCohomologyReviewSetup.{0}) :=
  ⟨GaloisCohomologyReviewSetup.example⟩

end ABC3.Interface.LocProP
