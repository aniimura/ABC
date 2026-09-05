import ABC3.Found.Falt1.AlmostEtale

/-!
# [Falt1] Hochschild `H²` の almost 消滅、コチェインの明示形(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、
remark 2.1(v)(物理 p.6=印字 p.259)。

`AlmostEtale.lean` の `hochschild_ext_almost_zero` は remark 2.1(v) の
almost 版を **`Ext^{k+1}_{B⊗_AB}(B,M)` の言葉で全次数について**証明して
いる。本ファイルはその `H²` の場合を、**コチェインの明示形**——
「2-コサイクル `c` に対し `p^n·c` はコバウンダリであり、その homotopy を
witness `w` から陽に作れる」——で与える。

## なぜ明示形が要るか

`Theorem 2.2`(冪零イデアルに沿った lifting)の証明の第2段は、
持ち上げた `A`-加群写像 `φ_ε` の積の非可換性の障害

    c(b₁,b₂) := p^ε·φ_ε(b₁b₂) - φ_ε(b₁)·φ_ε(b₂) ∈ I

が Hochschild 2-コサイクルであることを使って `φ_ε` を補正する、という
議論である。ここで必要なのは「`Ext²` が消える」という抽象的な事実では
なく、**具体的なコサイクル `c` に対する具体的なコバウンダリ表示**である。

## 構成(remark(v)の null-homotopy そのもの)

`w = Σᵢ xᵢ⊗yᵢ` を `Definition 2.1` 条件(iii)の witness とすると、

    h(b) := Σᵢ xᵢ·c(yᵢ, b)

が求める homotopy になる。実際、コサイクル条件を `(yᵢ, b₁, b₂)` で使うと

    b₁·h(b₂) - h(b₁b₂) + b₂·h(b₁)
      = Σᵢ b₁xᵢ·c(yᵢ,b₂) + (Σᵢ xᵢyᵢ)·c(b₁,b₂) - Σᵢ xᵢ·c(b₁yᵢ,b₂)

となり、第1項と第3項は `(b₁⊗1)·w = (1⊗b₁)·w`(annihilation、
`almost_swap_annihilate`)で打ち消し合い、第2項は `Σᵢ xᵢyᵢ = p^n`
(augmentation、`almost_swap_augment`)で `p^n·c(b₁,b₂)` になる。
Faltings が remark(v) で書いている縮約ホモトピー
`b₀⊗b₁⊗⋯ ↦ Σxᵢ⊗yᵢb₀⊗b₁⊗⋯` の、コチェイン側での姿である。

`Σᵢ` は Lean では `TensorProduct.lift` で扱う(有限和に分解せずに済む)。

## `H¹` について

`B` が可換で `M` が対称両側加群のときコバウンダリは 0 なので
`HH¹(B,M) = Der_A(B,M) = Hom_B(Ω[B⁄A], M)` であり、その almost 消滅は
`AlmostDifferentials.lean` の `kaehler_almost_zero`(`p^n` が `Ω[B⁄A]` を
零化する)から直ちに従う。したがって本ファイルは `H²` だけを扱う。
-/

namespace ABC3.Found.Falt1

universe u

/-- homotopy `h(b) := Σᵢ xᵢ·c(yᵢ,b)` の台になる双線形写像
(`u ⊗ v ↦ (b ↦ u·c(v,b))`)。`TensorProduct.lift` に渡して使う。 -/
noncomputable def hochH2Bil {A B M : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [AddCommGroup M] [Module B M] [Module A M] [IsScalarTower A B M]
    (c : B →ₗ[A] B →ₗ[A] M) : B →ₗ[A] B →ₗ[A] (B →ₗ[A] M) :=
  LinearMap.mk₂ A (fun u v => (LinearMap.lsmul B M u).restrictScalars A ∘ₗ (c v))
    (by intro u u' v; ext b; simp)
    (by intro a u v; ext b; simp)
    (by intro u v v'; ext b; simp)
    (by intro a u v; ext b; simp)

@[simp]
theorem hochH2Bil_apply {A B M : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [AddCommGroup M] [Module B M] [Module A M] [IsScalarTower A B M]
    (c : B →ₗ[A] B →ₗ[A] M) (u v b : B) :
    (TensorProduct.lift (hochH2Bil c)) (u ⊗ₜ[A] v) b = u • c v b := rfl

/-- **鍵の恒等式**(任意の `z : B⊗_AB` について成立)。
`z` について `A`-線形なので純粋テンソルで確かめればよく、そこは
コサイクル条件を `(v,b₁,b₂)` で使って `u` 倍するだけ。 -/
theorem hochH2_identity {A B M : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [AddCommGroup M] [Module B M] [Module A M] [IsScalarTower A B M]
    (c : B →ₗ[A] B →ₗ[A] M)
    (hc : ∀ v b₁ b₂ : B, v • c b₁ b₂ - c (v*b₁) b₂ + c v (b₁*b₂) - b₂ • c v b₁ = 0)
    (b₁ b₂ : B) (z : TensorProduct A B B) :
    b₁ • (TensorProduct.lift (hochH2Bil c) z) b₂
      - (TensorProduct.lift (hochH2Bil c) z) (b₁*b₂)
      + b₂ • (TensorProduct.lift (hochH2Bil c) z) b₁
    = (Algebra.TensorProduct.lmul' A z) • c b₁ b₂
      + (TensorProduct.lift (hochH2Bil c) ((b₁ ⊗ₜ[A] (1:B)) * z)) b₂
      - (TensorProduct.lift (hochH2Bil c) (((1:B) ⊗ₜ[A] b₁) * z)) b₂ := by
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul u v =>
    rw [Algebra.TensorProduct.tmul_mul_tmul, Algebra.TensorProduct.tmul_mul_tmul,
      Algebra.TensorProduct.lmul'_apply_tmul]
    show b₁ • (u • c v b₂) - u • c v (b₁*b₂) + b₂ • (u • c v b₁)
      = (u*v) • c b₁ b₂ + (b₁*u) • c (1*v) b₂ - (1*u) • c (b₁*v) b₂
    have h := hc v b₁ b₂
    have h2 := congrArg (fun m : M => u • m) h
    simp only [smul_sub, smul_add, smul_zero, smul_smul] at h2
    rw [one_mul, one_mul, mul_comm b₁ v]
    rw [smul_smul, smul_smul, mul_comm b₁ u, mul_comm b₂ u]
    linear_combination (norm := abel) -h2
  | add z z' hz hz' =>
    simp only [map_add, mul_add, LinearMap.add_apply, smul_add, add_smul]
    linear_combination (norm := abel) hz + hz'

/-- **`H²` の almost 消滅、明示形**。`Definition 2.1` 条件(iii)の witness
`w`(annihilation `hw_ann` と augmentation `hw_aug` を満たすもの)があれば、
任意の Hochschild 2-コサイクル `c` について `p^n·c` はコバウンダリであり、
その homotopy は `h := TensorProduct.lift (hochH2Bil c) w` で**陽に**与えられる。 -/
theorem hochschild_H2_almost_coboundary {A B M : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [AddCommGroup M] [Module B M] [Module A M] [IsScalarTower A B M]
    (p : A) (n : ℕ) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = p ^ n • (1 : B))
    (c : B →ₗ[A] B →ₗ[A] M)
    (hc : ∀ v b₁ b₂ : B, v • c b₁ b₂ - c (v*b₁) b₂ + c v (b₁*b₂) - b₂ • c v b₁ = 0) :
    ∃ h : B →ₗ[A] M, ∀ b₁ b₂ : B,
      (p ^ n) • c b₁ b₂ = b₁ • h b₂ - h (b₁*b₂) + b₂ • h b₁ := by
  refine ⟨TensorProduct.lift (hochH2Bil c) w, fun b₁ b₂ => ?_⟩
  have hid := hochH2_identity c hc b₁ b₂ w
  have hcancel : ((b₁ ⊗ₜ[A] (1:B)) * w) = (((1:B) ⊗ₜ[A] b₁) * w) := by
    have h0 := hw_ann b₁
    rw [sub_mul, sub_eq_zero] at h0
    exact h0.symm
  rw [hcancel, hw_aug] at hid
  rw [hid, Algebra.smul_def, mul_one, ← algebraMap_smul B (p ^ n) (c b₁ b₂)]
  abel

/-- `Definition 2.1` の仮定から直接述べた形。 -/
theorem hochschild_H2_almost_coboundary_of_isAlmostEtale {A B M : Type u}
    [CommRing A] [CommRing B] [Algebra A B] [Module.Free A B]
    [AddCommGroup M] [Module B M] [Module A M] [IsScalarTower A B M]
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
    (c : B →ₗ[A] B →ₗ[A] M)
    (hc : ∀ v b₁ b₂ : B, v • c b₁ b₂ - c (v*b₁) b₂ + c v (b₁*b₂) - b₂ • c v b₁ = 0) :
    ∃ h : B →ₗ[A] M, ∀ b₁ b₂ : B,
      (p ^ n) • c b₁ b₂ = b₁ • h b₂ - h (b₁*b₂) + b₂ • h b₁ :=
  hochschild_H2_almost_coboundary p n w
    (fun q => almost_swap_annihilate p hAE hf0inj n w hw q)
    (almost_swap_augment p hAE hf0inj n w hw) c hc

end ABC3.Found.Falt1
