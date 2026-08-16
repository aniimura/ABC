import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.RingTheory.Ideal.GoingUp
import Mathlib.NumberTheory.RamificationInertia.Basic

/-!
# [GenEll] §1 —— 有限素点の引き戻し(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

## ★★`HeightOneSpectrum.comap` は mathlib に**無い**(2026-08-17 実測)

`IsDedekindDomain.HeightOneSpectrum` に**素点の引き戻し**は無い
(`DedekindDomain/` 配下を `comap` で検索し、`HeightOneSpectrum` の comap は 0 件)。

★無限素点の側には `NumberField.InfinitePlace.comap` が**ある**のに、
有限素点の側には無い——**非対称である**。

## ★何に要るか

`degNormalized` の**底変換不変性**の有限側。
`w : 𝕍(K)^non` の係数を `e(w|v)·n_v`(`v = w` の引き戻し)と定めるとき、
まず `w ↦ v` の写像が要る。

★アルキメデス側は `Found/GenEll/InfinitePlaceRel.lean` で済んでいる
(`sum_arc_base_change`)。

## ★構成は 3 行で済む

- `asIdeal` は `Ideal.comap`
- 素であることは `Ideal.comap_isPrime`
- ★**非零であることが唯一の中身**で、`Ideal.IsIntegral.comap_ne_bot`
  (整拡大では素イデアルの縮約が非零)が mathlib に**ある**
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable (F K : Type*) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

/-- ★**有限素点の引き戻し** `𝕍(K)^non → 𝕍(F)^non`。

★mathlib に無いので作る(上の docstring)。
非零性は `Ideal.IsIntegral.comap_ne_bot`——**整拡大では素イデアルの縮約が非零**。 -/
def hosComap (w : HeightOneSpectrum (𝓞 K)) : HeightOneSpectrum (𝓞 F) where
  asIdeal := w.asIdeal.comap (algebraMap (𝓞 F) (𝓞 K))
  isPrime := Ideal.comap_isPrime _ _
  ne_bot := Ideal.IsIntegral.comap_ne_bot (𝓞 F) w.ne_bot

@[simp] theorem hosComap_asIdeal (w : HeightOneSpectrum (𝓞 K)) :
    (hosComap F K w).asIdeal = w.asIdeal.comap (algebraMap (𝓞 F) (𝓞 K)) := rfl

/-- ★`w` は自身の引き戻しの上にある(`v ≤ w ∩ 𝓞_F` の形)。 -/
theorem hosComap_le_comap (w : HeightOneSpectrum (𝓞 K)) :
    (hosComap F K w).asIdeal ≤ w.asIdeal.comap (algebraMap (𝓞 F) (𝓞 K)) := le_rfl

/-- ★引き戻した素点は `w` の下にある——`Ideal.LiesOver` の形。 -/
instance hosComap_liesOver (w : HeightOneSpectrum (𝓞 K)) :
    w.asIdeal.LiesOver (hosComap F K w).asIdeal :=
  ⟨rfl⟩

/-! ## ★★引き戻しの fiber —— `w | v` の集合

★mathlib の `Ideal.sum_ramification_inertia` は
**`IsDedekindDomain.primesOverFinset` の上の和**として述べられており、
`HeightOneSpectrum` ではない。ここでその 2 つを繋ぐ。 -/

/-- ★`v` の上にある素点は**有限個**である。

★`primesOverFinset`(mathlib、有限)への単射で押さえる。
`HeightOneSpectrum` は `asIdeal` で決まる(残る 2 つのフィールドは Prop)から単射。 -/
theorem finite_hosComap_fiber (v : HeightOneSpectrum (𝓞 F)) :
    {w : HeightOneSpectrum (𝓞 K) | hosComap F K w = v}.Finite := by
  haveI : v.asIdeal.IsMaximal := v.isPrime.isMaximal v.ne_bot
  refine Set.Finite.of_finite_image (f := HeightOneSpectrum.asIdeal) ?_ ?_
  · refine Set.Finite.subset
      (IsDedekindDomain.primesOverFinset v.asIdeal (𝓞 K)).finite_toSet ?_
    rintro _ ⟨w, hw, rfl⟩
    simp only [Set.mem_setOf_eq] at hw
    rw [Finset.mem_coe, IsDedekindDomain.mem_primesOverFinset_iff v.ne_bot]
    exact ⟨w.isPrime, ⟨(congrArg HeightOneSpectrum.asIdeal hw).symm⟩⟩
  · intro a _ b _ h
    exact HeightOneSpectrum.ext h

open scoped Classical in
/-- ★`v` の上にある素点全体の `Finset`。 -/
noncomputable def hosFiber (v : HeightOneSpectrum (𝓞 F)) :
    Finset (HeightOneSpectrum (𝓞 K)) :=
  (finite_hosComap_fiber F K v).toFinset

open scoped Classical in
@[simp] theorem mem_hosFiber {v : HeightOneSpectrum (𝓞 F)}
    {w : HeightOneSpectrum (𝓞 K)} :
    w ∈ hosFiber F K v ↔ hosComap F K w = v := by
  simp [hosFiber]

open scoped Classical in
/-- ★★**基本等式の `HeightOneSpectrum` 版** —— `Σ_{w | v} e(w|v)·f(w|v) = [K:F]`。

★mathlib の `Ideal.sum_ramification_inertia` は
**`IsDedekindDomain.primesOverFinset` の上の和**であり、
`HeightOneSpectrum` ではない(2026-08-17 実測)。
★本定理がその 2 つを繋ぐ——`asIdeal` による**全単射**である。 -/
theorem sum_ramification_inertia_hos (v : HeightOneSpectrum (𝓞 F)) :
    ∑ w ∈ hosFiber F K v,
        v.asIdeal.ramificationIdx w.asIdeal * v.asIdeal.inertiaDeg w.asIdeal
      = Module.finrank F K := by
  classical
  haveI : v.asIdeal.IsMaximal := v.isPrime.isMaximal v.ne_bot
  rw [← Ideal.sum_ramification_inertia (𝓞 K) F K v.ne_bot]
  refine Finset.sum_bij (i := fun w _ => w.asIdeal) ?_ ?_ ?_ ?_
  · intro w hw
    rw [mem_hosFiber] at hw
    rw [IsDedekindDomain.mem_primesOverFinset_iff v.ne_bot]
    exact ⟨w.isPrime, ⟨(congrArg HeightOneSpectrum.asIdeal hw).symm⟩⟩
  · intro a _ b _ h
    exact HeightOneSpectrum.ext h
  · intro P hP
    rw [IsDedekindDomain.mem_primesOverFinset_iff v.ne_bot] at hP
    have hPne : P ≠ ⊥ := Ideal.ne_bot_of_mem_primesOver v.ne_bot hP
    refine ⟨⟨P, hP.1, hPne⟩, ?_, rfl⟩
    rw [mem_hosFiber]
    exact HeightOneSpectrum.ext hP.2.over.symm
  · intro a _
    rfl

/-! ## ★出典の紐付け(`.src`) -/

def hosComap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

end ABC3.Found.GenEll
