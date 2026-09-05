import ABC3.Found.Falt1.Section1
import ABC3.Found.Falt1.Section4

/-!
# [Falt1] Chapter I §1 —— `Theorem 1.2` の塔の**非空虚性**(2026-09-05)

`thm_1_2_of_ring_tower`(および `thm_1_2_of_ring_tower_monogenic`)は
原文の設定そのままの形で `δₙ → 0` を主張するが、その**仮定の束が
同時に満たせること**(＝空虚でないこと)は別に示す必要がある。

仮定のうち `Ψ`(塔の仮定「`Ω_{V_{n+1}/Vₙ}` が `(V_{n+1}/pV_{n+1})^{d+1}`
を商に持つ」)は `p` が真の非単元であることを要求するので、
**本当に分岐した無限塔**でなければ満たせない——退化した(定数の)塔では
`Ω = 0` になり `Ψ` が全射になれない。

そこで標数 2 の切断多項式環の塔

    Vₙ := 𝔽₂[X]/(X^{2^{n+1}}),   遷移写像 Vₙ → V_{n+1} は `X ↦ 0`

を使う。標数 2 なので `d(X^{2^{n+2}}) = 2^{n+2}·X^{…}·dX = 0` となり、
**`Ω_{V_{n+1}/Vₙ}` は階数 1 の自由加群**(≠ 0)になる。

## 内容

* `Vt`——塔そのもの。`truncPoly_*`(`Section1.lean`)で局所・ネーター・
  アルティン・極大イデアルがすべて揃う。
* `VtHom`/`VtAlgebra`/`VtTower`——遷移写像と代数構造。
* `VtOmegaEquiv`——**`Ω[V_{n+1}⁄Vₙ] ≅ V_{n+1}`**。
  `kaehler_map_bijective_of_D_eq_zero`(底の取り換え)と
  `truncPoly_omegaEquiv`(`Ω[k[X]/(Xᴺ)⁄k] ≅ B`)の合成。
-/

namespace ABC3.Found.Falt1

open Polynomial

/-- **標数 2 の切断多項式環の塔**——`Vt n = 𝔽₂[X]/(X^{2^{n+1}})`。 -/
abbrev Vt (n : ℕ) : Type :=
  Polynomial (ZMod 2) ⧸ Ideal.span {(X : Polynomial (ZMod 2)) ^ (2 ^ (n + 1))}

/-- 遷移写像 `Vt n → Vt (n+1)`——`X ↦ 0`(つまり `𝔽₂` を経由する)。 -/
noncomputable def VtHom (n : ℕ) : Vt n →ₐ[ZMod 2] Vt (n + 1) := by
  refine Ideal.Quotient.liftₐ _ (Polynomial.aeval (0 : Vt (n + 1))) ?_
  intro a ha
  refine Submodule.span_induction ?_ ?_ ?_ ?_ ha
  · rintro _ rfl
    rw [map_pow, Polynomial.aeval_X, zero_pow (by positivity)]
  · simp
  · intro x y _ _ hx hy; rw [map_add, hx, hy, add_zero]
  · intro c x _ hx; rw [smul_eq_mul, map_mul, hx, mul_zero]

noncomputable instance VtAlgebra (n : ℕ) : Algebra (Vt n) (Vt (n + 1)) :=
  (VtHom n).toRingHom.toAlgebra

instance VtTower (n : ℕ) : IsScalarTower (ZMod 2) (Vt n) (Vt (n + 1)) :=
  IsScalarTower.of_algebraMap_eq (fun x => by
    show algebraMap (ZMod 2) (Vt (n + 1)) x = VtHom n (algebraMap (ZMod 2) (Vt n) x)
    rw [AlgHom.commutes])

/-- **遷移写像の像の微分は消える**——`X ↦ 0` なので像は `𝔽₂` の像に入る。 -/
theorem VtHom_D_eq_zero (n : ℕ) (a : Vt n) :
    KaehlerDifferential.D (ZMod 2) (Vt (n + 1)) (algebraMap (Vt n) (Vt (n + 1)) a) = 0 := by
  obtain ⟨f, rfl⟩ := Ideal.Quotient.mk_surjective
    (I := Ideal.span {(X : Polynomial (ZMod 2)) ^ (2 ^ (n + 1))}) a
  show KaehlerDifferential.D (ZMod 2) (Vt (n + 1))
    (VtHom n (Ideal.Quotient.mk _ f)) = 0
  rw [VtHom, Ideal.Quotient.liftₐ_apply, Ideal.Quotient.lift_mk]
  have key : ∀ g : Polynomial (ZMod 2),
      KaehlerDifferential.D (ZMod 2) (Vt (n + 1))
        (Polynomial.aeval (0 : Vt (n + 1)) g) = 0 := by
    intro g
    rw [Polynomial.aeval_def, Polynomial.eval₂_at_zero]
    exact (KaehlerDifferential.D (ZMod 2) (Vt (n + 1))).map_algebraMap _
  exact key f

/-- **`Ω[V_{n+1}⁄Vₙ] ≅ Ω[V_{n+1}⁄𝔽₂]`**——底の取り換え。 -/
theorem VtKaehlerMap_bijective (n : ℕ) :
    Function.Bijective
      (KaehlerDifferential.map (ZMod 2) (Vt n) (Vt (n + 1)) (Vt (n + 1))) :=
  kaehler_map_bijective_of_D_eq_zero (ZMod 2) (Vt n) (Vt (n + 1)) (VtHom_D_eq_zero n)

/-- **`Ω[V_{n+1}⁄Vₙ] ≅ V_{n+1}`**——階数 1 の自由加群。
標数 2 なので `d(X^{2^{n+2}}) = 0` になるのが効いている。 -/
noncomputable def VtOmegaEquiv (n : ℕ) :
    Ω[(Vt (n + 1))⁄(Vt n)] ≃ₗ[Vt (n + 1)] Vt (n + 1) :=
  (LinearEquiv.ofBijective
      (KaehlerDifferential.map (ZMod 2) (Vt n) (Vt (n + 1)) (Vt (n + 1)))
      (VtKaehlerMap_bijective n)).symm.trans
    (truncPoly_omegaEquiv (k := ZMod 2) (2 ^ (n + 2)) (by
      have : ((2 : ℕ) : ZMod 2) = 0 := by decide
      rw [Nat.cast_pow, this, zero_pow (by positivity)]))

/-! ## 局所環としてのインスタンスと生成元 -/

instance VtLocal (n : ℕ) : IsLocalRing (Vt n) :=
  truncPoly_isLocalRing _ (by positivity)

instance VtArtinianRing (n : ℕ) : IsArtinianRing (Vt n) :=
  truncPoly_isArtinianRing _

/-- 極大イデアルの生成元 `X̄`。 -/
noncomputable def VtGen (n : ℕ) : Vt n :=
  Ideal.Quotient.mk _ (X : Polynomial (ZMod 2))

theorem VtGen_span (n : ℕ) :
    Ideal.span ({VtGen n} : Set (Vt n)) = IsLocalRing.maximalIdeal (Vt n) :=
  (truncPoly_maximalIdeal _ (by positivity)).symm

/-- `Ω[V_{n+1}⁄Vₙ]` は `V_{n+1}` と同型なのでアルティン。 -/
instance VtOmegaArtinian (n : ℕ) : IsArtinian (Vt (n + 1)) (Ω[(Vt (n + 1))⁄(Vt n)]) :=
  isArtinian_of_linearEquiv (VtOmegaEquiv n).symm

instance VtOmegaNoetherian (n : ℕ) : IsNoetherian (Vt (n + 1)) (Ω[(Vt (n + 1))⁄(Vt n)]) :=
  isNoetherian_of_linearEquiv (VtOmegaEquiv n).symm

/-! ## 塔の仮定 `Ψ`——`Ω[V_{n+1}⁄Vₙ] ⊗ V_{n+1}` は `V_{n+1}/𝔪` を商に持つ -/

/-- **塔の仮定の実現**——`Ω[V_{n+1}⁄Vₙ] ⊗ V_{n+1} ↠ (V_{n+1}/𝔪)^1`。 -/
noncomputable def VtPsi (n : ℕ) :
    TensorProduct (Vt (n + 1)) (Vt (n + 1)) (Ω[(Vt (n + 1))⁄(Vt n)]) →ₗ[Vt (n + 1)]
      (Fin 1 → Vt (n + 1) ⧸ Ideal.span ({VtGen (n + 1)} : Set (Vt (n + 1)))) :=
  LinearMap.pi (fun _ : Fin 1 =>
    (Ideal.span ({VtGen (n + 1)} : Set (Vt (n + 1)))).mkQ ∘ₗ
      (VtOmegaEquiv n).toLinearMap ∘ₗ
      (TensorProduct.lid (Vt (n + 1)) (Ω[(Vt (n + 1))⁄(Vt n)])).toLinearMap)

set_option maxHeartbeats 4000000 in
/-- 合成が全射であること(`VtPsi` の各成分)。 -/
theorem VtPsiComp_surjective (n : ℕ) : Function.Surjective
    ((Ideal.span ({VtGen (n + 1)} : Set (Vt (n + 1)))).mkQ ∘ₗ
      (VtOmegaEquiv n).toLinearMap ∘ₗ
      (TensorProduct.lid (Vt (n + 1)) (Ω[(Vt (n + 1))⁄(Vt n)])).toLinearMap) := by
  intro z
  obtain ⟨v, hv⟩ := Submodule.mkQ_surjective
    (Ideal.span ({VtGen (n + 1)} : Set (Vt (n + 1)))) z
  refine ⟨(TensorProduct.lid (Vt (n + 1)) (Ω[(Vt (n + 1))⁄(Vt n)])).symm
    ((VtOmegaEquiv n).symm v), ?_⟩
  rw [LinearMap.coe_comp, Function.comp_apply, LinearMap.coe_comp, Function.comp_apply,
    LinearEquiv.coe_coe, LinearEquiv.coe_coe, LinearEquiv.apply_symm_apply,
    LinearEquiv.apply_symm_apply]
  exact hv

set_option maxHeartbeats 4000000 in
theorem VtPsi_surjective (n : ℕ) : Function.Surjective (VtPsi n) := by
  intro y
  obtain ⟨x, hx⟩ := VtPsiComp_surjective n (y 0)
  refine ⟨x, funext fun i => ?_⟩
  have hi : i = 0 := Subsingleton.elim i 0
  subst hi
  exact hx

end ABC3.Found.Falt1
