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
    (t : A) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = t • (1 : B))
    (c : B →ₗ[A] B →ₗ[A] M)
    (hc : ∀ v b₁ b₂ : B, v • c b₁ b₂ - c (v*b₁) b₂ + c v (b₁*b₂) - b₂ • c v b₁ = 0) :
    ∃ h : B →ₗ[A] M, ∀ b₁ b₂ : B,
      t • c b₁ b₂ = b₁ • h b₂ - h (b₁*b₂) + b₂ • h b₁ := by
  refine ⟨TensorProduct.lift (hochH2Bil c) w, fun b₁ b₂ => ?_⟩
  have hid := hochH2_identity c hc b₁ b₂ w
  have hcancel : ((b₁ ⊗ₜ[A] (1:B)) * w) = (((1:B) ⊗ₜ[A] b₁) * w) := by
    have h0 := hw_ann b₁
    rw [sub_mul, sub_eq_zero] at h0
    exact h0.symm
  rw [hcancel, hw_aug] at hid
  rw [hid, Algebra.smul_def, mul_one, ← algebraMap_smul B t (c b₁ b₂)]
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
  hochschild_H2_almost_coboundary (p ^ n) w
    (fun q => almost_swap_annihilate p hAE hf0inj n w hw q)
    (almost_swap_augment p hAE hf0inj n w hw) c hc

/-! ## `Theorem 2.2` の障害類——2-コサイクルであることの代数的な核

`Theorem 2.2` の第2段は、第1段で得た `A`-加群写像 `ψ`
(`π∘ψ = p^n·φ`、`almost_lift_of_isAlmostEtale`)の積の非可換性

    c(b₁,b₂) := p^n·ψ(b₁b₂) - ψ(b₁)·ψ(b₂)

を Hochschild 2-コサイクルとして扱う。以下の2本がその代数的な核である。 -/

/-- **障害はイデアル `I = ker π` に入る**。`π(c(b₁,b₂)) = p^{2n}φ(b₁)φ(b₂)
- p^{2n}φ(b₁)φ(b₂) = 0`。 -/
theorem obstruction_mem_ker {A B C D : Type u} [CommRing A] [CommRing B] [CommRing C] [CommRing D]
    [Algebra A B] [Algebra A C] [Algebra A D]
    (c : A) (π : C →ₐ[A] D) (φ : B →ₐ[A] D) (ψ : B →ₗ[A] C)
    (hψ : ∀ b : B, π (ψ b) = c • φ b) (b₁ b₂ : B) :
    π (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂) = 0 := by
  rw [map_sub, map_smul, map_mul π (ψ b₁) (ψ b₂), hψ, hψ, hψ, map_mul φ]
  simp only [Algebra.smul_def]
  ring

/-- **障害の `p^n` 倍されたコサイクル恒等式**。任意の `A`-線形 `ψ` について
恒等式として成り立つ(仮定は一切要らない——`C` の可換性だけで `ring` が
閉じる)。`x ∈ I := ker π` に対して `ψ(b)·x = p^n·(b•x)`(`b•x` は
`φ` 経由の `B`-作用)なので、これは**コサイクル条件を `p^n` 倍したもの**
である。したがって `C` が `p` 捩れ無し(Faltings の標準仮定)であれば
`c` は honest な 2-コサイクルであり、`hochschild_H2_almost_coboundary` を
適用できる。 -/
theorem obstruction_identity {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (c : A) (ψ : B →ₗ[A] C) (b₁ b₂ b₃ : B) :
    ψ b₁ * (c • ψ (b₂ * b₃) - ψ b₂ * ψ b₃)
      - c • (c • ψ (b₁ * b₂ * b₃) - ψ (b₁ * b₂) * ψ b₃)
      + c • (c • ψ (b₁ * (b₂ * b₃)) - ψ b₁ * ψ (b₂ * b₃))
      - ψ b₃ * (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂) = 0 := by
  rw [mul_assoc b₁ b₂ b₃]
  simp only [Algebra.smul_def]
  ring

/-! ## Faltings の「`ε` を倍にする」段(物理p.7=印字p.260)

原文: *"Doubling `ε` and then enlarging it a little we may assume that this
class vanishes, and then we can modify `φ_ε` so that it becomes
multiplicative, i.e. such that `p^ε φ_ε(xy) = φ_ε(x)φ_ε(y)`."*

この2文が指す操作は、以下の2本にちょうど分解される。 -/

/-- **「`ε` を少し大きくする」段**。`ψ` を `p^d` 倍すると、水準は `n` から
`n+d` に上がり、障害は `p^{2d}` 倍になる。したがって
`H²` の almost 消滅が与える `p^d·c = d h` と合わせると、
**`p^d ψ` の障害 `p^{2d}c = d(p^d h)` はちょうどコバウンダリになる**
(= class が消える)。 -/
theorem rescale_obstruction {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (c t : A) (ψ : B →ₗ[A] C) (b₁ b₂ : B) :
    (c * t) • ((t • ψ) (b₁ * b₂)) - ((t • ψ) b₁) * ((t • ψ) b₂)
      = (t * t) • (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂) := by
  simp only [LinearMap.smul_apply, Algebra.smul_def, map_mul]
  ring

/-- **「`ε` を倍にする」段——Faltings の `Theorem 2.2` の核心**。
障害がコバウンダリ(`hcob`——`ψ` 側にスケールした形で書いてある。
`x ∈ I` に対し `ψ(b)·x = p^n·(b•x)` なので、これは抽象的な
コバウンダリ条件 `c = dh` を `p^n` 倍したものである)であり、
`h` が `I`(`I²=0`)に値を取る(`hIsq`)なら、
**`ψ'' := p^n·ψ + h` は水準 `2n` で厳密に乗法的**になる:

    p^{2n}·ψ''(b₁b₂) = ψ''(b₁)·ψ''(b₂).

証明は可換環の恒等式1本(`linear_combination q^n * hcob - hIsq`、
`q := algebraMap A C p`)——`h(b₁)h(b₂)=0` と コバウンダリ関係で
すべての項が打ち消し合う。 -/
theorem doubling_multiplicative {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (c : A) (ψ h : B →ₗ[A] C)
    (hIsq : ∀ b₁ b₂ : B, h b₁ * h b₂ = 0)
    (hcob : ∀ b₁ b₂ : B, c • (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂)
        = ψ b₁ * h b₂ - c • h (b₁ * b₂) + ψ b₂ * h b₁)
    (b₁ b₂ : B) :
    (c * c) • ((c • ψ + h) (b₁ * b₂))
      = ((c • ψ + h) b₁) * ((c • ψ + h) b₂) := by
  have hc := hcob b₁ b₂
  have hsq := hIsq b₁ b₂
  simp only [LinearMap.add_apply, LinearMap.smul_apply, Algebra.smul_def, map_mul] at hc ⊢
  linear_combination ((algebraMap A C) c) * hc - hsq

/-- **一意性の代数的な核**(Faltings の *"Such a lifting is unique up to
`H¹(B/A,I)`"*)。水準 `n` の乗法的持ち上げが2つあれば、その差
`d := ψ' - ψ` は `p^n·d(xy) = ψ(x)·d(y) + ψ(y)·d(x)` を満たす——
`x ∈ I` に対し `ψ(b)·x = p^n·(b•x)` なので、これは `d` が(`p^n` 倍
された意味で)**導分**であることを言っている。`C` が `p` 捩れ無しなら
honest な導分になり、`Ω[B⁄A]` の almost 消滅
(`AlmostDifferentials.kaehler_almost_zero`)から `d` は `p` 冪で消える。 -/
theorem uniqueness_derivation_eq {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (c : A) (ψ ψ' : B →ₗ[A] C)
    (hmul : ∀ x y : B, c • ψ (x * y) = ψ x * ψ y)
    (hmul' : ∀ x y : B, c • ψ' (x * y) = ψ' x * ψ' y)
    (hsq : ∀ x y : B, (ψ' x - ψ x) * (ψ' y - ψ y) = 0)
    (x y : B) :
    c • ((ψ' - ψ) (x * y)) = ψ x * ((ψ' - ψ) y) + ψ y * ((ψ' - ψ) x) := by
  have h1 := hmul x y
  have h2 := hmul' x y
  have h3 := hsq x y
  simp only [LinearMap.sub_apply, smul_sub] at *
  linear_combination h2 - h1 + h3

/-- **貼り合わせの整合性の核**(Faltings の *"the different `φ_ε` glue
together"*)。水準 `n` の乗法的持ち上げ `ψ` を `p^k` 倍したものは、
水準 `n+k` の乗法的持ち上げになる。したがって一意性と合わせると
`ψ_{n+k} = p^k·ψ_n` が従い、族が整合する。 -/
theorem rescale_multiplicative {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (c t : A) (ψ : B →ₗ[A] C)
    (hmul : ∀ x y : B, c • ψ (x * y) = ψ x * ψ y) (x y : B) :
    (c * t) • ((t • ψ) (x * y)) = ((t • ψ) x) * ((t • ψ) y) := by
  have h := hmul x y
  simp only [LinearMap.smul_apply, Algebra.smul_def, map_mul] at h ⊢
  linear_combination ((algebraMap A C) t * (algebraMap A C) t) * h

/-! ## `Theorem 2.2` の残り(正直な記録)

上の `doubling_multiplicative` により、`Theorem 2.2` の証明の
「`p^ε φ_ε(xy) = φ_ε(x)φ_ε(y)` を満たす `φ_ε` を作る」という部分は
**完全に閉じた**(材料は `almost_lift_of_isAlmostEtale`(第1段)・
`obstruction_mem_ker`・`obstruction_identity`(障害がコサイクル)・
`hochschild_H2_almost_coboundary`(コバウンダリ表示)・
`rescale_obstruction`・`doubling_multiplicative`)。

残るのは Faltings の証明の最後の2文である:

1. *"Such a lifting is unique up to `H¹(B/A,I)`, hence up to `p`-torsion.
   As `C` has no such torsion, the different `φ_ε` glue together to give a
   multiplicative `A`-linear map `φ₀ : mB → C`."*
   ——一意性の核は `uniqueness_derivation_eq`、整合性の核は
   `rescale_multiplicative` で閉じている。残るのは族の極限を取る操作。
2. *"We can extend to `B = A + mB`, because `φ₀` maps `p^ε` to an element
   `x = p^ε + y` (`y ∈ I`) satisfying `x² = p^ε x`, hence `p^ε y = 0`."*

### ★重要な構造上の発見(2026-09-05)——なぜ `2.4` は閉じて `2.2` は閉じないか

本プロジェクトの `Definition 2.1` の条件(iii)は指数を **`n : ℕ`** で
取っている(`∀ n : ℕ, 0 < n → ∃ e, diagonalCompare p e = p^n • elem`)。
原典は `ε > 0` を**`p` 冪分母を持つ正の有理数**で走らせる——つまり
`A` が `p^{1/p^k}` を含む(perfectoid 的な)底を前提とし、
`m = ∪_ε p^ε A` は `(p)` より真に小さい。

この違いが `2.2`/`2.3` と `2.4` を分ける:

- **`Theorem 2.4` は almost な結論**(「`m` が零化する」)なので、
  整数指数の枠組みでもそのまま意味を持ち、実際に閉じた
  (`GaloisTransfer.thm_2_4_ii`・`AlmostDifferentials.thm_2_4_i_cokernel`)。
- **`Theorem 2.2`/`2.3` は honest な結論**(「`φ` は**一意に**持ち上がる」)
  であり、`p^ε φ_ε` の族から `ε → 0` の極限として honest な `φ₀` を得る
  操作が要る。整数指数では `ε ≥ 1` までしか下がらないので、
  **`p` で割る操作が原理的に構成できない**。すなわち `2.2`/`2.3` を
  閉じるには、まず**`p`-可除な底(`m = ∪ p^{1/p^k}`)を持つ層**を
  新たに設計する必要がある——これは `ResearchPaper/mathlib-gap.json` の
  `falt1-almost-mathematics` そのものである。

したがって次に着手すべきは「`Theorem 2.2` の続き」ではなく、
**`p`-可除な底の上の almost mathematics の層**の設計である。
本ファイルまでに揃えた材料(第1段・障害のコサイクル性・コバウンダリ
表示・倍化・一意性・整合性)は、その層の上でそのまま再利用できる。
-/

end ABC3.Found.Falt1
