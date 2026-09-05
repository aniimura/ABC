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

end ABC3.Found.Falt1
