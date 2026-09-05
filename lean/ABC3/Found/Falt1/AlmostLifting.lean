import ABC3.Found.Falt1.AlmostDifferentials
import ABC3.Found.Falt1.AlmostProjective
import ABC3.Found.Falt1.HochschildLowDegree

/-!
# [Falt1] Theorem 2.2 の一意性——完全に証明(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、物理 p.7
(印字 p.260)。

内容 (Falt1 p.7、260dpi 目視): **2.2. Theorem.** Suppose `B = A + mB` is an
almost étale covering of `A`, `C` an `A`-algebra, `I ⊂ C` a nilpotent ideal,
and `φ: B → C/I` an `A`-algebra morphism. Then `φ` **lifts uniquely** to
`B → C`.

## 本ファイルの到達点——**一意性は完全に閉じた**

Faltings の証明の最後は *"Uniqueness has already been established"* で
締められており、その根拠は *"Such a lifting is unique up to `H¹(B/A,I)`,
hence up to `p`-torsion. As `C` has no such torsion..."* である。

これを完全に形式化した(`thm_2_2_uniqueness`・
`thm_2_2_uniqueness_of_isAlmostEtale`)。議論は3段:

1. 2つの持ち上げ `ψ, ψ' : B →ₐ[A] C` の差 `d := ψ' - ψ` は、
   `I² = 0` から **`A`-導分** `B → C`(`C` は `ψ` 経由で `B`-代数)になる:
   `d(xy) = ψ(x)·d(y) + ψ(y)·d(x)`。
2. 導分は `Ω[B⁄A]` の普遍性を経由する(`Derivation.liftKaehlerDifferential`)
   ので、`p^n` が `Ω[B⁄A]` を零化すること
   (`AlmostDifferentials.kaehler_almost_zero`——`Definition 2.1` 条件(iii)の
   witness から直接示した)から `p^n·d = 0`。
2'. これが `H¹(B/A,I)` が `m` で零化されるという Faltings の主張の中身
   である(可換な `B` と対称両側加群では `HH¹ = Der`)。
3. `C` が `p^n` 捩れ無し(Faltings の標準仮定)なので `d = 0`、すなわち
   `ψ = ψ'`。

**存在**の側(`φ` が実際に持ち上がること)は
`HochschildLowDegree.lean` の末尾に記した通り、`ε` 族の極限を取る操作を
要し、そこには `p`-可除な底(`m = ∪p^{1/p^k}`)を持つ almost mathematics の
層が要る。一意性はその層を要しないので、ここで閉じられる。
-/

namespace ABC3.Found.Falt1

universe u

/-- **`Ω` の almost 消滅から導分の almost 消滅へ**。導分は `Ω[B⁄A]` の
普遍性を経由する(`Derivation.liftKaehlerDifferential`)ので、`p^n` が
`Ω[B⁄A]` を零化すれば任意の `A`-導分も `p^n` で零化される。
(`B` 線形写像に `A`-スカラーを通すので `map_smul_of_tower` を使う、
`tools/lean-idioms.md` #50(c)。) -/
theorem derivation_almost_zero {A B M : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [AddCommGroup M] [Module B M] [Module A M] [IsScalarTower A B M]
    (c : A) (hΩ : ∀ x : Ω[B⁄A], c • x = 0)
    (D : Derivation A B M) (b : B) : c • D b = 0 := by
  have h := D.liftKaehlerDifferential.map_smul_of_tower c (KaehlerDifferential.D A B b)
  rw [hΩ (KaehlerDifferential.D A B b), map_zero] at h
  rw [← Derivation.liftKaehlerDifferential_comp_D D b, ← h]

/-- **`Theorem 2.2` の一意性、一般形**。`p^n` が `Ω[B⁄A]` を零化し、
`C` が `p^n` 捩れ無しなら、差が二乗零(`hsq`——`I²=0` から出る)である
ような2つの `A`-代数写像 `ψ, ψ' : B →ₐ[A] C` は一致する。 -/
theorem thm_2_2_uniqueness {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (c : A) (hΩ : ∀ x : Ω[B⁄A], c • x = 0)
    (htors : ∀ x : C, (algebraMap A C) c * x = 0 → x = 0)
    (ψ ψ' : B →ₐ[A] C) (hsq : ∀ x y : B, (ψ' x - ψ x) * (ψ' y - ψ y) = 0) :
    ψ = ψ' := by
  letI : Algebra B C := ψ.toRingHom.toAlgebra
  haveI : IsScalarTower A B C := by
    apply IsScalarTower.of_algebraMap_eq
    intro x
    show algebraMap A C x = ψ (algebraMap A B x)
    rw [AlgHom.commutes]
  let D : Derivation A B C :=
    { toLinearMap := ψ'.toLinearMap - ψ.toLinearMap
      map_one_eq_zero' := by simp
      leibniz' := by
        intro x y
        have h := hsq x y
        show ψ' (x*y) - ψ (x*y) = ψ x * (ψ' y - ψ y) + ψ y * (ψ' x - ψ x)
        rw [map_mul, map_mul]
        linear_combination h }
  have hzero : ∀ b : B, D b = 0 := by
    intro b
    apply htors
    have h := derivation_almost_zero c hΩ D b
    rwa [Algebra.smul_def] at h
  ext b
  have hb := hzero b
  show ψ b = ψ' b
  have h2 : ψ' b - ψ b = 0 := hb
  linear_combination -h2

/-- **`Theorem 2.2` の一意性、原典の仮定の形**。`B` が `A` の almost étale
covering(`Definition 2.1`)であり、`π : C →ₐ[A] D` の核 `I` が二乗零、
`C` が `p^n` 捩れ無しなら、`π` に沿った持ち上げは高々1つしかない。 -/
theorem thm_2_2_uniqueness_of_isAlmostEtale {A B C D : Type u} [CommRing A] [CommRing B]
    [CommRing C] [CommRing D] [Algebra A B] [Algebra A C] [Algebra A D] [Module.Free A B]
    (p : A) (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p)))
    (π : C →ₐ[A] D)
    (hIsq : ∀ u v : C, π u = 0 → π v = 0 → u * v = 0)
    (htors : ∀ x : C, (algebraMap A C) (p ^ n) * x = 0 → x = 0)
    (ψ ψ' : B →ₐ[A] C) (hlift : ∀ b : B, π (ψ b) = π (ψ' b)) :
    ψ = ψ' := by
  refine thm_2_2_uniqueness (p ^ n)
    (fun x => kaehler_almost_zero_of_isAlmostEtale p hAE hf0inj n w hw x)
    htors ψ ψ' (fun x y => ?_)
  refine hIsq _ _ ?_ ?_ <;>
  · rw [map_sub, hlift]
    ring

end ABC3.Found.Falt1
