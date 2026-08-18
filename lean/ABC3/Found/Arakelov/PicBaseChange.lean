import ABC3.Meta.Claim
import Mathlib.RingTheory.PicardGroup
import Mathlib.LinearAlgebra.TensorProduct.Tower

/-!
# Arakelov (B2) 第 197 ブロック —— **係数環を同型で取り替えても可逆**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`R` と `Γ(Spec R, ⊤)` の橋

`Spec` を使うと必ず出てくる摩擦がある——`R` と `Γ(Spec R, ⊤)` は**同型だが等しくない**。
`invertible_gammaCarrier`(第 132 系)は `R` 上の可逆性を与えるが、
Interface が要求するのは `Γ(Spec R, ⊤)` 上の可逆性である。

## ★★テンソルで渡すのが一番短い

`Module.Invertible A (A ⊗[R] M)`(mathlib の instance)を使い、

    algebraMap R S が全射  ⟹  S ⊗[R] M ≃ₗ[S] M

を示せば `Module.Invertible.congr` で移る。★★全射性だけで足りる
(同型である必要は無い)——`s ⊗ m = 1 ⊗ (s • m)` が言えるからである。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `bcAct` / `bcAct_tmul` | ★係数拡大の写像(値は `rfl`) |
| `bcActEquiv` | ★★★`algebraMap` が全射なら係数拡大は何もしない |
| `invertible_of_surjective_algebraMap` | ★★★★**係数環を取り替えても可逆** |
-/

universe u

namespace ABC3.Found.Arakelov

open TensorProduct

variable {R S M : Type u} [CommRing R] [CommRing S] [Algebra R S]
  [AddCommGroup M] [Module R M] [Module S M] [IsScalarTower R S M]

noncomputable def bcAct : (S ⊗[R] M) →ₗ[S] M :=
  TensorProduct.AlgebraTensorModule.lift
    { toFun := fun s => s • (LinearMap.id : M →ₗ[R] M)
      map_add' := fun a b => by ext m; simp [add_smul]
      map_smul' := fun a b => by ext m; simp [mul_smul] }

@[simp] theorem bcAct_tmul (s : S) (m : M) : bcAct (s ⊗ₜ[R] m) = s • m := rfl

/-- ★★`algebraMap` が全射なら係数拡大は何もしない。 -/
noncomputable def bcActEquiv (hs : Function.Surjective (algebraMap R S)) :
    (S ⊗[R] M) ≃ₗ[S] M := by
  refine LinearEquiv.ofBijective bcAct ⟨?_, fun m => ⟨(1 : S) ⊗ₜ[R] m, by simp⟩⟩
  intro a b hab
  have key : ∀ t : S ⊗[R] M, t = (1 : S) ⊗ₜ[R] (bcAct t) := by
    intro t
    induction t using TensorProduct.induction_on with
    | zero => simp
    | tmul s m =>
        obtain ⟨r, rfl⟩ := hs s
        rw [bcAct_tmul, algebraMap_smul, ← TensorProduct.smul_tmul]
        congr 1
        rw [Algebra.smul_def, mul_one]
    | add x y hx hy =>
        rw [map_add, TensorProduct.tmul_add]
        exact congrArg₂ (· + ·) hx hy
  rw [key a, key b, hab]


/-- ★★★★**係数環を全射(同型)で取り替えても可逆性は保たれる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `R` と `Γ(Spec R, ⊤)` の橋である。 -/
theorem invertible_of_surjective_algebraMap
    (hs : Function.Surjective (algebraMap R S)) [Module.Invertible R M] :
    Module.Invertible S M :=
  Module.Invertible.congr (bcActEquiv hs)

/-! ## ★出典の紐付け(`.src`) -/

def invertible_of_surjective_algebraMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——係数環を同型で取り替えても可逆)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
