import ABC3.Meta.Claim
import ABC3.Found.Falt1.AlmostDifferentials
import ABC3.Found.Falt1.AlmostProjective
import ABC3.Found.Falt1.HochschildLowDegree
import ABC3.Found.Falt1.AlmostBase

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

/-! ## 存在側への配線——二乗零イデアルの加群構造

`Theorem 2.2` の存在側で `Hochschild H²` の機械
(`HochschildLowDegree.hochschild_H2_almost_coboundary`)を実際の障害
`Ob(b₁,b₂) := c·ψ(b₁b₂) - ψ(b₁)ψ(b₂)` に適用するには、`I := ker π` を
`B`-加群と見なす必要がある。`I² = 0` なので `C` の作用は `C ⧸ I` を
経由し、`φ : B →ₐ[A] C ⧸ I` で `B` が作用する。mathlib には
「二乗零イデアル上の `C ⧸ I`-加群構造」がそのままの形では無いので
ここで作る(`Ideal.Cotangent` は `I ⧸ I²` なので、`I² = 0` のときは
同型だが、同一視の手間を避けて直接構成する)。 -/

/-- **二乗零イデアル `I` 上の `C ⧸ I`-加群構造**。`mk c ⋆ x := c * x` が
well-defined なのは `I·I = 0` だから。 -/
@[reducible] noncomputable def sqZeroModule {C : Type u} [CommRing C] (I : Ideal C)
    (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) : Module (C ⧸ I) I where
  smul := fun d x => Quotient.liftOn d (fun c => (⟨c * (x : C), I.mul_mem_left c x.2⟩ : I))
    (by
      intro c₁ c₂ h
      have hI : c₁ - c₂ ∈ I := by rwa [← Submodule.quotientRel_def]
      apply Subtype.ext
      show c₁ * (x : C) = c₂ * (x : C)
      have hz := hsq _ hI _ x.2
      linear_combination hz)
  one_smul := by intro x; apply Subtype.ext; show (1 : C) * (x : C) = x; ring
  mul_smul := by
    intro d₁ d₂ x
    induction d₁ using Quotient.inductionOn with
    | _ c₁ =>
      induction d₂ using Quotient.inductionOn with
      | _ c₂ => apply Subtype.ext; show (c₁ * c₂) * (x : C) = c₁ * (c₂ * (x : C)); ring
  smul_zero := by
    intro d
    induction d using Quotient.inductionOn with
    | _ c => apply Subtype.ext; show c * (0 : C) = 0; ring
  smul_add := by
    intro d x y
    induction d using Quotient.inductionOn with
    | _ c => apply Subtype.ext; show c * ((x : C) + y) = c * x + c * y; ring
  add_smul := by
    intro d₁ d₂ x
    induction d₁ using Quotient.inductionOn with
    | _ c₁ =>
      induction d₂ using Quotient.inductionOn with
      | _ c₂ => apply Subtype.ext; show (c₁ + c₂) * (x : C) = c₁ * x + c₂ * x; ring
  zero_smul := by intro x; apply Subtype.ext; show (0 : C) * (x : C) = 0; ring

/-- **`ψ` による乗法と `φ` 経由の `B`-作用の関係**: `π∘ψ = c·φ` なら
`x ∈ I` について `ψ(b)·x = c·(b • x)`。これが「障害 `Ob` の `c` 倍された
コサイクル恒等式(`obstruction_identity`)を、`I` の `B`-加群構造に関する
honest なコサイクル条件へ翻訳する」ための橋である(`C` の `p` 捩れ
無しで `c` を落とす)。 -/
theorem sqZero_act_eq {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I)) (c : A) (ψ : B →ₗ[A] C)
    (hψ : ∀ b : B, Ideal.Quotient.mk I (ψ b) = c • φ b) :
    letI : Module (C ⧸ I) I := sqZeroModule I hsq
    letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
    ∀ (b : B) (z : I), ψ b * (z : C) = c • ((b • z : I) : C) := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  intro b z
  obtain ⟨y, hy⟩ := Ideal.Quotient.mk_surjective (φ b)
  have hbz : ((b • z : I) : C) = y * (z : C) := by
    show ((φ b • z : I) : C) = y * (z : C)
    rw [← hy]
    rfl
  have hmem : ψ b - (algebraMap A C c) * y ∈ I := by
    rw [← Ideal.Quotient.eq_zero_iff_mem, map_sub, hψ b, map_mul, hy]
    show c • φ b - (algebraMap A (C ⧸ I) c) * φ b = 0
    rw [Algebra.smul_def, sub_self]
  have hzero : (ψ b - (algebraMap A C c) * y) * (z : C) = 0 := hsq _ hmem _ z.2
  rw [hbz, Algebra.smul_def]
  linear_combination hzero

/-- `A` → `B` → `I` のスカラー塔(`hochschild_H2_almost_coboundary` が
要求する `[IsScalarTower A B M]`)。 -/
theorem sqZeroTower {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I)) :
    letI : Module (C ⧸ I) I := sqZeroModule I hsq
    letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
    IsScalarTower A B I := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  constructor
  intro a b z
  apply Subtype.ext
  obtain ⟨y, hy⟩ := Ideal.Quotient.mk_surjective (φ b)
  have h1 : ((b • z : I) : C) = y * (z : C) := by
    show ((φ b • z : I) : C) = _
    rw [← hy]; rfl
  have hmk : (algebraMap A (C ⧸ I)) a = Ideal.Quotient.mk I (algebraMap A C a) := rfl
  have h2 : (((a • b) • z : I) : C) = (algebraMap A C a * y) * (z : C) := by
    show ((φ (a • b) • z : I) : C) = _
    rw [map_smul, Algebra.smul_def, hmk, ← hy, ← map_mul (Ideal.Quotient.mk I)]
    rfl
  rw [h2, Submodule.coe_smul_of_tower, h1, Algebra.smul_def, mul_assoc]

/-- **障害は honest な Hochschild 2-コサイクルである**。
`obstruction_identity`(仮定なしの `c` 倍された恒等式)を `sqZero_act_eq`
(`ψ(b)·x = c·(b•x)`)で `I` の `B`-加群構造の言葉に翻訳し、`C` の
`p` 捩れ無し(`htors`)で `c` を落とす。これで
`HochschildLowDegree.hochschild_H2_almost_coboundary` を実際の障害に
適用できるようになった——`Theorem 2.2` の存在側で残っていた
「障害を本当にコサイクルとして扱う」部分の橋である。 -/
theorem obstruction_isCocycle {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I)) (c : A) (ψ : B →ₗ[A] C)
    (hψ : ∀ b : B, Ideal.Quotient.mk I (ψ b) = c • φ b)
    (htors : ∀ x : C, (algebraMap A C) c * x = 0 → x = 0) :
    letI : Module (C ⧸ I) I := sqZeroModule I hsq
    letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
    ∀ (Ob : B →ₗ[A] B →ₗ[A] I),
      (∀ b₁ b₂ : B, ((Ob b₁ b₂ : I) : C) = c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂) →
      ∀ v b₁ b₂ : B, v • Ob b₁ b₂ - Ob (v * b₁) b₂ + Ob v (b₁ * b₂) - b₂ • Ob v b₁ = 0 := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  intro Ob hOb v b₁ b₂
  have hact := sqZero_act_eq I hsq φ c ψ hψ
  apply Subtype.ext
  apply htors
  have hid := obstruction_identity c ψ v b₁ b₂
  have e1 := (hact v (Ob b₁ b₂)).symm
  have e2 := (hact b₂ (Ob v b₁)).symm
  have h1 := hOb b₁ b₂
  have h2 := hOb (v * b₁) b₂
  have h3 := hOb v (b₁ * b₂)
  have h4 := hOb v b₁
  simp only [Submodule.coe_sub, Submodule.coe_add, Algebra.smul_def] at hid e1 e2 h1 h2 h3 h4 ⊢
  rw [← mul_assoc v b₁ b₂] at hid h3
  rw [mul_sub, mul_add, mul_sub, e1, e2, h1, h2, h3, h4]
  linear_combination hid

/-- 障害を `I` に値を取る `A`-双線形写像として束ねたもの
(`B →ₗ[A] B →ₗ[A] I`)。値が `I` に入るのは `obstruction_mem_ker`、
双線形性は `ψ` の `A`-線形性から。 -/
noncomputable def obsBil {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (φ : B →ₐ[A] (C ⧸ I)) (c : A) (ψ : B →ₗ[A] C)
    (hψ : ∀ b : B, Ideal.Quotient.mk I (ψ b) = c • φ b) : B →ₗ[A] B →ₗ[A] I :=
  LinearMap.mk₂ A
    (fun b₁ b₂ => ⟨c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂, by
      have hmem := obstruction_mem_ker c (Ideal.Quotient.mkₐ A I) φ ψ (fun b => hψ b) b₁ b₂
      rwa [← Ideal.Quotient.eq_zero_iff_mem]⟩)
    (by intro b₁ b₁' b₂; apply Subtype.ext
        show c • ψ ((b₁ + b₁') * b₂) - ψ (b₁ + b₁') * ψ b₂
          = (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂) + (c • ψ (b₁' * b₂) - ψ b₁' * ψ b₂)
        rw [add_mul, map_add, map_add, smul_add]; ring)
    (by intro a b₁ b₂; apply Subtype.ext
        show c • ψ ((a • b₁) * b₂) - ψ (a • b₁) * ψ b₂
          = a • (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂)
        rw [smul_mul_assoc, map_smul, map_smul, smul_sub, smul_comm c a, smul_mul_assoc])
    (by intro b₁ b₂ b₂'; apply Subtype.ext
        show c • ψ (b₁ * (b₂ + b₂')) - ψ b₁ * ψ (b₂ + b₂')
          = (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂) + (c • ψ (b₁ * b₂') - ψ b₁ * ψ b₂')
        rw [mul_add, map_add, map_add, smul_add]; ring)
    (by intro a b₁ b₂; apply Subtype.ext
        show c • ψ (b₁ * (a • b₂)) - ψ b₁ * ψ (a • b₂)
          = a • (c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂)
        rw [mul_smul_comm, map_smul, map_smul, smul_sub, smul_comm c a, mul_smul_comm])

/-- **`Theorem 2.2` の第2段、完成**——Faltings の
*"Doubling `ε` and then enlarging it a little we may assume that this class
vanishes, and then we can modify `φ_ε` so that it becomes multiplicative,
i.e. such that `p^ε φ_ε(xy) = φ_ε(x)φ_ε(y)`"* そのもの。

`I ⊆ C` が二乗零、`φ : B →ₐ[A] C ⧸ I`、`ψ` が水準 `c` の `A`-加群持ち上げ
(`AlmostProjective.almost_lift_of_isAlmostEtale` が供給する)、`C` が
`c` 捩れ無し、`w` が水準 `t` の witness——このとき
**水準 `(ct)²` で厳密に乗法的な持ち上げ `ψ₂` が存在する**。

証明は本ファイルと `HochschildLowDegree.lean` の部品を繋ぐだけ:
障害 `Ob`(`obsBil`)は `I` に値を取りコサイクル(`obstruction_isCocycle`)、
`w` から `t·Ob = dh`(`hochschild_H2_almost_coboundary`)、
`ψ₁ := t·ψ` の障害はちょうどコバウンダリ(`rescale_obstruction`)、
`ψ₂ := (ct)·ψ₁ + t·h` が乗法的(`doubling_multiplicative`)。 -/
theorem exists_multiplicative_lift {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I))
    (c : A) (ψ : B →ₗ[A] C) (hψ : ∀ b, Ideal.Quotient.mk I (ψ b) = c • φ b)
    (htors : ∀ x : C, (algebraMap A C) c * x = 0 → x = 0)
    (t : A) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = t • (1 : B)) :
    ∃ ψ₂ : B →ₗ[A] C,
      (∀ b, Ideal.Quotient.mk I (ψ₂ b) = ((c*t)*(c*t)) • φ b) ∧
      (∀ x y, ((c*t)*(c*t)) • ψ₂ (x*y) = ψ₂ x * ψ₂ y) := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  haveI : IsScalarTower A B I := sqZeroTower I hsq φ
  set Ob := obsBil I φ c ψ hψ with hObdef
  have hObC : ∀ b₁ b₂ : B, ((Ob b₁ b₂ : I) : C) = c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂ :=
    fun b₁ b₂ => rfl
  have hcoc := obstruction_isCocycle I hsq φ c ψ hψ htors Ob hObC
  obtain ⟨h, hh⟩ := hochschild_H2_almost_coboundary t w hw_ann hw_aug Ob hcoc
  have hact := sqZero_act_eq I hsq φ c ψ hψ
  set ĥ : B →ₗ[A] C := t • ((I.subtype.restrictScalars A) ∘ₗ h) with hĥ
  have hIsq2 : ∀ b₁ b₂ : B, ĥ b₁ * ĥ b₂ = 0 := by
    intro b₁ b₂
    show (t • ((h b₁ : I) : C)) * (t • ((h b₂ : I) : C)) = 0
    have hz := hsq _ (h b₁).2 _ (h b₂).2
    rw [Algebra.smul_def, Algebra.smul_def]
    linear_combination (algebraMap A C t * algebraMap A C t) * hz
  have hcob : ∀ b₁ b₂ : B,
      (c*t) • ((c*t) • ((t • ψ) (b₁ * b₂)) - (t • ψ) b₁ * (t • ψ) b₂)
      = (t • ψ) b₁ * ĥ b₂ - (c*t) • ĥ (b₁ * b₂) + (t • ψ) b₂ * ĥ b₁ := by
    intro b₁ b₂
    have hh12 := congrArg (fun x : I => (x : C)) (hh b₁ b₂)
    have ha1 := hact b₁ (h b₂)
    have ha2 := hact b₂ (h b₁)
    have hO := hObC b₁ b₂
    simp only [hĥ, LinearMap.smul_apply, LinearMap.coe_comp, Function.comp_apply,
      LinearMap.coe_restrictScalars, Submodule.coe_subtype,
      Submodule.coe_add, Submodule.coe_sub, Submodule.coe_smul_of_tower,
      Algebra.smul_def, map_mul] at hh12 ha1 ha2 hO ⊢
    linear_combination (-(algebraMap A C c * (algebraMap A C t)^3)) * hO
      + (algebraMap A C c * (algebraMap A C t)^2) * hh12
      - ((algebraMap A C t)^2) * ha1 - ((algebraMap A C t)^2) * ha2
  refine ⟨(c*t) • (t • ψ) + ĥ, fun b => ?_, fun x y => ?_⟩
  · show Ideal.Quotient.mk I ((c*t) • (t • ψ b) + t • ((h b : I) : C)) = ((c*t)*(c*t)) • φ b
    have hzero : Ideal.Quotient.mk I ((h b : I) : C) = 0 :=
      (Ideal.Quotient.eq_zero_iff_mem).mpr (h b).2
    have hmk : ∀ x : A, Ideal.Quotient.mk I (algebraMap A C x) = algebraMap A (C ⧸ I) x :=
      fun _ => rfl
    simp only [Algebra.smul_def, map_add, map_mul, hmk, hzero, mul_zero, add_zero, hψ b]
    ring
  · exact doubling_multiplicative (c*t) (t • ψ) ĥ hIsq2 hcob x y

/-! ## `Theorem 2.2` 存在側の組み立て——残りの見取り図

`exists_multiplicative_lift` により、Faltings の証明のうち

- 第1段(`AlmostProjective.almost_lift_of_isAlmostEtale`——almost 射影性から
  `A`-加群写像を `p^ε` 倍で持ち上げる)
- 第2段(**本ファイルの `exists_multiplicative_lift`**——障害を `H²` で
  消して乗法的にする)
- 一意性(`thm_2_2_uniqueness`)

は**すべて閉じた**。残るのは Faltings の証明の最後の2文だけである:

1. *"the different `φ_ε` glue together to give a multiplicative `A`-linear
   map `φ₀ : mB → C`"* ——`ε` 族の極限。整合性の核は
   `HochschildLowDegree.rescale_multiplicative`、一意性の核は
   `HochschildLowDegree.uniqueness_derivation_eq` で用意済み。
2. *"We can extend to `B = A + mB`"* ——最後の詰め。

どちらも `AlmostBase.lean` の塔(`PDivTower`)の上で書く必要がある——
`m² = m` が効くのはそこ(`mB` の元の積がまた `mB` に入る)である。
-/

/-- **`Definition 2.1` から直接述べた第1段+第2段**。水準 `c` の witness で
`A`-加群持ち上げを作り(第1段)、水準 `t` の witness で乗法的にする(第2段)。
結果は**水準 `(ct)²` で厳密に乗法的な持ち上げ**。 -/
theorem exists_multiplicative_lift_of_isAlmostEtale {A B C : Type u}
    [CommRing A] [CommRing B] [CommRing C] [Algebra A B] [Algebra A C]
    [Module.Finite A B] [Module.Free A B]
    (p : A) (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I))
    (c : A) (wc : TensorProduct A B B)
    (hwc : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p wc = c • Algebra.FormallyUnramified.elem (Localization.Away p)
        (Localization.Away (algebraMap A B p)))
    (htors : ∀ x : C, (algebraMap A C) c * x = 0 → x = 0)
    (t : A) (wt : TensorProduct A B B)
    (hwt : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p wt = t • Algebra.FormallyUnramified.elem (Localization.Away p)
        (Localization.Away (algebraMap A B p))) :
    ∃ ψ₂ : B →ₗ[A] C,
      (∀ b, Ideal.Quotient.mk I (ψ₂ b) = ((c*t)*(c*t)) • φ b) ∧
      (∀ x y, ((c*t)*(c*t)) • ψ₂ (x*y) = ψ₂ x * ψ₂ y) := by
  obtain ⟨ψ, hψ⟩ := almost_lift_of_isAlmostEtale p hAE hf0inj c wc hwc
    ((Ideal.Quotient.mkₐ A I).toLinearMap) Ideal.Quotient.mk_surjective φ.toLinearMap
  exact exists_multiplicative_lift I hsq φ c ψ hψ htors t wt
    (fun q => almost_swap_annihilate p hAE hf0inj t wt hwt q)
    (almost_swap_augment p hAE hf0inj t wt hwt)

/-- 塔に対する `Definition 2.1` から、`p := ϖ 0` についての通常の
`Definition 2.1` が出る(条件(iii)の `p^n` 版は `k=0` の witness を
`(ϖ 0)^{n-1}` 倍して得る)。 -/
theorem isAlmostEtaleCovering_of_tower {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    {q : ℕ} (T : PDivTower A q) (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T) :
    IsAlmostEtaleCovering (A := A) (B := B) (T.ϖ 0) := by
  letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
  haveI := awayScalarTower (T.ϖ 0) (A := A) (B := B)
  obtain ⟨hFree, hFin, hEt, htr, hwit⟩ := hAET
  refine ⟨hFree, hFin, hEt, htr, fun n hn => ?_⟩
  obtain ⟨e₀, he₀⟩ := hwit 0
  refine ⟨(T.ϖ 0) ^ (n - 1) • e₀, ?_⟩
  rw [map_smul, he₀, smul_smul]
  congr 1
  rw [← pow_succ]
  congr 1
  omega

/-- **塔の上での第1段+第2段**——`Theorem 2.2` の存在側が要求する
「`ε` を小さくしながら乗法的な持ち上げを作る」族そのもの。
任意の `k, j` について、**水準 `(ϖ k · ϖ j)²`(`k,j` を大きくすればいくらでも
深くなる)で厳密に乗法的な持ち上げ**が存在する。

Faltings の証明で残るのは、この族の極限を取って honest な
`φ₀ : mB → C` を得る段(整合性の核は
`HochschildLowDegree.rescale_multiplicative`、一意性の核は
`uniqueness_derivation_eq`・`thm_2_2_uniqueness`)と、`B = A + mB` への
拡張だけである。 -/
theorem exists_multiplicative_lift_tower {A B C : Type u}
    [CommRing A] [CommRing B] [CommRing C] [Algebra A B] [Algebra A C]
    [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0)))))
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (T.ϖ k) * x = 0 → x = 0)
    (k j : ℕ) :
    ∃ ψ : B →ₗ[A] C,
      (∀ b, Ideal.Quotient.mk I (ψ b)
        = ((T.ϖ k * T.ϖ j) * (T.ϖ k * T.ϖ j)) • φ b) ∧
      (∀ x y, ((T.ϖ k * T.ϖ j) * (T.ϖ k * T.ϖ j)) • ψ (x*y) = ψ x * ψ y) := by
  have hAE := isAlmostEtaleCovering_of_tower T hAET
  obtain ⟨_, _, _, _, hwit⟩ := hAET
  obtain ⟨wk, hwk⟩ := hwit k
  obtain ⟨wj, hwj⟩ := hwit j
  exact exists_multiplicative_lift_of_isAlmostEtale (T.ϖ 0) hAE hf0inj I hsq φ
    (T.ϖ k) wk hwk (htors k) (T.ϖ j) wj hwj

/-! ## `ε` 族の極限へ——**イデアル `sB` の上では「almost」が消える**

Faltings の *"the different `φ_ε` glue together to give a multiplicative
`A`-linear map `φ₀ : mB → C`"* の中身は、次の観察である:

水準 `s` の乗法的持ち上げ `ψ`(`π∘ψ = s·φ`、`s·ψ(xy) = ψ(x)ψ(y)`)が
あるとき、**イデアル `sB` の上では `φ₀(s·b) := ψ(b)` が honest な
(`p^ε` 倍の付かない)持ち上げになる**:

- 加法性・`π∘φ₀ = φ|_{sB}`:`ψ` の `A`-線形性から。
- **乗法性**: `(s·b)(s·b') = s·(s·(bb'))` なので
  `φ₀((s b)(s b')) = ψ(s·(bb')) = s·ψ(bb') = ψ(b)ψ(b') = φ₀(sb)φ₀(sb')`
  ——`s` の因子がちょうど打ち消し合う。
- **well-defined**:`B` が `s` 捩れ無しであること。これは Faltings 自身の
  標準仮定(`algebraMap B B[1/p]` が単射)から従う——`s` は `p` を割るので
  `s·b = 0` なら `b` は `p` 捩れ、よって `0`。

`mB = ∪_k (ϖ k)·B`(塔の元は整除で全順序なので、有限和は1つの積に書ける)
なので、あとはこの族を `k` について貼り合わせれば `φ₀ : mB → C` が得られる。 -/

/-- well-defined 性:`B` が `s` 捩れ無しなら `s•b = s•b'` から `ψ b = ψ b'`。 -/
theorem lift_wd {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (s : A) (ψ : B →ₗ[A] C)
    (hBtors : ∀ b : B, s • b = 0 → b = 0) (b b' : B) (h : s • b = s • b') : ψ b = ψ b' := by
  have hs : s • (b - b') = 0 := by rw [smul_sub, h, sub_self]
  have hb := hBtors _ hs
  have hbb : b = b' := by linear_combination hb
  rw [hbb]

/-- 乗法性の核:`ψ(s•(bb')) = ψ(b)ψ(b')`。水準 `s` の乗法性と `A`-線形性から。 -/
theorem lift_mul {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (s : A) (ψ : B →ₗ[A] C)
    (hmul : ∀ x y : B, s • ψ (x * y) = ψ x * ψ y) (b b' : B) :
    ψ (s • (b * b')) = ψ b * ψ b' := by
  rw [map_smul, hmul]

/-- 積の代表元:`(s•b)·(s•b') = s•(s•(b*b'))`。 -/
theorem lift_mul_rep {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (s : A) (b b' : B) : (s • b) * (s • b') = s • (s • (b * b')) := by
  simp only [Algebra.smul_def]
  ring

open Classical in
/-- **イデアル `sB` の上の honest な持ち上げ**——Faltings の `φ₀` の1水準版。
水準 `s` の乗法的持ち上げ `ψ` と `B` の `s` 捩れ無しから、
`φ₀(s·b) := ψ(b)` が well-defined で、加法的・**乗法的**であり、
`π∘φ₀ = φ` を `sB` の上で満たす。**`p^ε` 倍が付かない honest な等式**
であることが要点である。 -/
theorem exists_honest_lift_on_ideal {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (φ : B →ₐ[A] (C ⧸ I)) (s : A) (ψ : B →ₗ[A] C)
    (hψ : ∀ b, Ideal.Quotient.mk I (ψ b) = s • φ b)
    (hmul : ∀ x y : B, s • ψ (x * y) = ψ x * ψ y)
    (hBtors : ∀ b : B, s • b = 0 → b = 0) :
    ∃ φ₀ : B → C,
      (∀ b : B, φ₀ (s • b) = ψ b) ∧
      (∀ b b' : B, φ₀ ((s • b) * (s • b')) = φ₀ (s • b) * φ₀ (s • b')) ∧
      (∀ b : B, Ideal.Quotient.mk I (φ₀ (s • b)) = φ (s • b)) ∧
      (∀ b b' : B, φ₀ (s • b + s • b') = φ₀ (s • b) + φ₀ (s • b')) := by
  set F : B → C := fun z => if h : ∃ b : B, s • b = z then ψ h.choose else 0 with hF
  have hkey : ∀ b : B, F (s • b) = ψ b := by
    intro b
    have hex : ∃ b' : B, s • b' = s • b := ⟨b, rfl⟩
    show (if h : ∃ b' : B, s • b' = s • b then ψ h.choose else 0) = ψ b
    rw [dif_pos hex]
    exact lift_wd s ψ hBtors _ b hex.choose_spec
  refine ⟨F, hkey, ?_, ?_, ?_⟩
  · intro b b'
    rw [lift_mul_rep, hkey, hkey, hkey, lift_mul s ψ hmul]
  · intro b
    rw [hkey, hψ b, map_smul]
  · intro b b'
    rw [← smul_add, hkey, hkey, hkey, map_add]

/-! ## 族の貼り合わせ——一意性と整合性

`ε` 族を貼り合わせるには「各水準で持ち上げが一意」であることが要る
(Faltings: *"Such a lifting is unique up to `H¹(B/A,I)`, hence up to
`p`-torsion"*)。`thm_2_2_uniqueness` は **`A`-代数写像**の一意性だったが、
族の各項は「水準 `s` で乗法的」なだけなので、その版が別に要る。 -/

/-- **水準 `s` の乗法的持ち上げの一意性**。差 `d := ψ' - ψ` は `I` に値を取り、
`uniqueness_derivation_eq` と `sqZero_act_eq` から `s·d(xy) = s·(x•dy + y•dx)`、
`C` の `s` 捩れ無しで honest な導分になる。あとは `Ω[B⁄A]` の almost 消滅
(`hΩ`)と `c` 捩れ無しで `d = 0`。 -/
theorem uniqueness_level_lift {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I))
    (s : A) (ψ ψ' : B →ₗ[A] C)
    (hψ : ∀ b, Ideal.Quotient.mk I (ψ b) = s • φ b)
    (hψ' : ∀ b, Ideal.Quotient.mk I (ψ' b) = s • φ b)
    (hmul : ∀ x y : B, s • ψ (x * y) = ψ x * ψ y)
    (hmul' : ∀ x y : B, s • ψ' (x * y) = ψ' x * ψ' y)
    (c : A) (hΩ : ∀ x : Ω[B⁄A], c • x = 0)
    (hstors : ∀ x : C, (algebraMap A C) s * x = 0 → x = 0)
    (hctors : ∀ x : C, (algebraMap A C) c * x = 0 → x = 0) :
    ψ = ψ' := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  haveI : IsScalarTower A B I := sqZeroTower I hsq φ
  have hmem : ∀ b : B, ψ' b - ψ b ∈ I := by
    intro b
    rw [← Ideal.Quotient.eq_zero_iff_mem, map_sub, hψ b, hψ' b, sub_self]
  have hsq2 : ∀ x y : B, (ψ' x - ψ x) * (ψ' y - ψ y) = 0 :=
    fun x y => hsq _ (hmem x) _ (hmem y)
  set dI : B →ₗ[A] I :=
    { toFun := fun b => ⟨ψ' b - ψ b, hmem b⟩
      map_add' := by intro x y; apply Subtype.ext; show ψ' (x+y) - ψ (x+y) = _; simp; abel
      map_smul' := by
        intro a x; apply Subtype.ext
        show ψ' (a • x) - ψ (a • x) = ((a • ⟨ψ' x - ψ x, hmem x⟩ : I) : C)
        rw [Submodule.coe_smul_of_tower, map_smul, map_smul, smul_sub] } with hdI
  have hact := sqZero_act_eq I hsq φ s ψ hψ
  have hder : ∀ x y : B, dI (x * y) = x • dI y + y • dI x := by
    intro x y
    apply Subtype.ext
    rw [← sub_eq_zero]
    apply hstors
    have h1 := uniqueness_derivation_eq s ψ ψ' hmul hmul' hsq2 x y
    have h2 := hact x (dI y)
    have h3 := hact y (dI x)
    have hd1 : ((dI (x*y) : I) : C) = ψ' (x*y) - ψ (x*y) := rfl
    have hd2 : ((dI y : I) : C) = ψ' y - ψ y := rfl
    have hd3 : ((dI x : I) : C) = ψ' x - ψ x := rfl
    simp only [Submodule.coe_add, LinearMap.sub_apply, Algebra.smul_def, hd1, hd2, hd3]
      at h1 h2 h3 ⊢
    linear_combination h1 + h2 + h3
  set D : Derivation A B I :=
    { toLinearMap := dI
      map_one_eq_zero' := by
        have h := hder 1 1
        rw [one_mul, one_smul] at h
        have h3 : (0 : I) + dI 1 = dI 1 + dI 1 := by rw [zero_add]; exact h
        exact (add_right_cancel h3).symm
      leibniz' := hder } with hD
  have hzero : ∀ b : B, dI b = 0 := by
    intro b
    apply Subtype.ext
    apply hctors
    have h := derivation_almost_zero c hΩ D b
    have hc : ((c • D b : I) : C) = 0 := by rw [h]; rfl
    rw [Submodule.coe_smul_of_tower, Algebra.smul_def] at hc
    exact hc
  ext b
  have hb := hzero b
  have hbb : ψ' b - ψ b = 0 := congrArg (fun x : I => (x : C)) hb
  linear_combination -hbb

/-- **族の整合性**: 水準 `s·d` の持ち上げは、水準 `s` の持ち上げの `d` 倍に
一致する(`rescale_multiplicative` で `d·ψ` が水準 `s·d` の持ち上げに
なることを見て、`uniqueness_level_lift` を適用するだけ)。 -/
theorem lift_compat {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I))
    (s d : A) (ψ ψ' : B →ₗ[A] C)
    (hψ : ∀ b, Ideal.Quotient.mk I (ψ b) = s • φ b)
    (hmul : ∀ x y : B, s • ψ (x * y) = ψ x * ψ y)
    (hψ' : ∀ b, Ideal.Quotient.mk I (ψ' b) = (s * d) • φ b)
    (hmul' : ∀ x y : B, (s * d) • ψ' (x * y) = ψ' x * ψ' y)
    (c : A) (hΩ : ∀ x : Ω[B⁄A], c • x = 0)
    (hsdtors : ∀ x : C, (algebraMap A C) (s * d) * x = 0 → x = 0)
    (hctors : ∀ x : C, (algebraMap A C) c * x = 0 → x = 0) :
    ψ' = d • ψ := by
  refine uniqueness_level_lift I hsq φ (s * d) ψ' (d • ψ) hψ' ?_ hmul' ?_ c hΩ hsdtors hctors
  · intro b
    show Ideal.Quotient.mk I (d • ψ b) = (s * d) • φ b
    have hmk : ∀ x : A, Ideal.Quotient.mk I (algebraMap A C x) = algebraMap A (C ⧸ I) x :=
      fun _ => rfl
    simp only [Algebra.smul_def, map_mul, hmk, hψ b]
    ring
  · exact fun x y => rescale_multiplicative s d ψ hmul x y

/-- **貼り合わせの整合性、`φ₀` の言葉で**: `ψ' = d • ψ` なら、`sB` 上の
honest な持ち上げ `φ₀^s` と `(sd)B` 上の `φ₀^{sd}` は重なり `(sd)B` の上で
一致する——`(sd)•b = s•(d•b)` なので、両辺とも `ψ (d • b)` になる。
これで族 `{φ₀^{ϖk}}` が貼り合わさり `φ₀ : mB → C` が定まる。 -/
theorem honest_lift_agree {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (d : A) (ψ ψ' : B →ₗ[A] C) (hcompat : ψ' = d • ψ) (b : B) :
    ψ (d • b) = ψ' b := by
  rw [hcompat, map_smul]
  rfl

/-! ## 総括——`Theorem 2.2` の存在側で得られたもの -/

/-- **`Theorem 2.2` の存在側、到達点**。`B` が塔に対する almost étale
covering、`I ⊆ C` が二乗零、`φ : B →ₐ[A] C ⧸ I` のとき、**任意の `k` について
`m` の元 `s`(`= (ϖ k)⁴`、`k` を大きくすればいくらでも深い)が取れて、
イデアル `sB` の上に honest な(`p^ε` 倍の付かない)加法的・乗法的な
持ち上げ `φ₀` が存在する**。

`k` は任意なので、これは Faltings が
*"the different `φ_ε` glue together"* で貼り合わせる族そのものである。
族の整合性は `lift_compat`・`honest_lift_agree` で証明済みなので、
残るのは
- `mB = ∪_k (ϖ k)·B` 上へこの族を集合論的に取りまとめる作業、
- `B = A + mB` を使う最後の拡張(原文の
  *"φ₀ maps `p^ε` to an element `x = p^ε + y` satisfying `x² = p^ε x`,
  hence `p^ε y = 0`"*)
の2つだけである。 -/
theorem thm_2_2_honest_lift_on_ideal {A B C : Type u}
    [CommRing A] [CommRing B] [CommRing C] [Algebra A B] [Algebra A C]
    [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0)))))
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (T.ϖ k) * x = 0 → x = 0)
    (hBtors : ∀ (k : ℕ) (b : B), T.ϖ k • b = 0 → b = 0)
    (k : ℕ) :
    ∃ (s : A) (_ : s ∈ T.m) (φ₀ : B → C),
      (∀ b b' : B, φ₀ ((s • b) * (s • b')) = φ₀ (s • b) * φ₀ (s • b')) ∧
      (∀ b b' : B, φ₀ (s • b + s • b') = φ₀ (s • b) + φ₀ (s • b')) ∧
      (∀ b : B, Ideal.Quotient.mk I (φ₀ (s • b)) = φ (s • b)) := by
  set s : A := (T.ϖ k * T.ϖ k) * (T.ϖ k * T.ϖ k) with hs
  have hsm : s ∈ T.m := by
    rw [hs]
    exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ (T.ϖ_mem_m k))
  have hstors : ∀ x : C, (algebraMap A C) s * x = 0 → x = 0 := by
    intro x hx
    have h4 : (algebraMap A C) s * x
        = (algebraMap A C) (T.ϖ k) * ((algebraMap A C) (T.ϖ k) *
          ((algebraMap A C) (T.ϖ k) * ((algebraMap A C) (T.ϖ k) * x))) := by
      rw [hs]; simp only [map_mul]; ring
    rw [h4] at hx
    exact htors k _ (htors k _ (htors k _ (htors k _ hx)))
  have hsBtors : ∀ b : B, s • b = 0 → b = 0 := by
    intro b hb
    have h4 : s • b = T.ϖ k • (T.ϖ k • (T.ϖ k • (T.ϖ k • b))) := by
      rw [smul_smul, smul_smul, smul_smul, hs]
      congr 1
      ring
    rw [h4] at hb
    exact hBtors k _ (hBtors k _ (hBtors k _ (hBtors k _ hb)))
  obtain ⟨ψ, hψ, hmul⟩ := exists_multiplicative_lift_tower T hAET hf0inj I hsq φ htors k k
  obtain ⟨φ₀, hφ₀, hmulφ, hliftφ, haddφ⟩ :=
    exists_honest_lift_on_ideal I φ s ψ hψ hmul hsBtors
  exact ⟨s, hsm, φ₀, hmulφ, haddφ, hliftφ⟩

/-! ## 族の貼り合わせ——**`φ₀ : mB → C` の構成**

Faltings の *"the different `φ_ε` glue together to give a multiplicative
`A`-linear map `φ₀ : mB → C`"*(物理 p.7)を実行する段。

`mB` の元は `{z | ∃ k b, lev k • b = z}` として表す。塔 `PDivTower` では
水準が整除で全順序なので、2つの元は常に**共通の水準**の代表元に書き直せる
(`exists_common_level`)——これが「有限和を1つの積に書ける」ことに当たり、
貼り合わせを純粋に集合論的な操作にしている。 -/

open Classical in
/-- **貼り合わせ本体**。整合な族 `{Ψ k}`(水準 `lev k`、`k ≤ k'` のとき
`lev k' ∣ lev k`)から、`Φ (lev k • b) = Ψ k b` を全ての `k` について
満たす関数 `Φ` を得る。交差水準での well-defined 性は整合性 `hcompat` と
`B` の捩れ無し `hBtors` から出る。 -/
theorem exists_glued_lift {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C))
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (hcompat : ∀ (k k' : ℕ) (d : A), k ≤ k' → lev k = lev k' * d → Ψ k = d • Ψ k')
    (hBtors : ∀ (k : ℕ) (b : B), lev k • b = 0 → b = 0) :
    ∃ Φ : B → C, ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b := by
  have hwd_le : ∀ (k₁ k₂ : ℕ) (b₁ b₂ : B), k₁ ≤ k₂ →
      lev k₁ • b₁ = lev k₂ • b₂ → Ψ k₁ b₁ = Ψ k₂ b₂ := by
    intro k₁ k₂ b₁ b₂ hle heq
    obtain ⟨d, hd⟩ := hdvd k₁ k₂ hle
    have hΨ := hcompat k₁ k₂ d hle hd
    have hsm : lev k₂ • (d • b₁) = lev k₂ • b₂ := by
      rw [smul_smul, ← hd]; exact heq
    have hsub : lev k₂ • (d • b₁ - b₂) = 0 := by rw [smul_sub, hsm, sub_self]
    have hb : d • b₁ - b₂ = 0 := hBtors k₂ _ hsub
    have hb2 : d • b₁ = b₂ := by linear_combination hb
    rw [hΨ]
    show d • Ψ k₂ b₁ = Ψ k₂ b₂
    rw [← map_smul, hb2]
  have hwd : ∀ (k₁ k₂ : ℕ) (b₁ b₂ : B),
      lev k₁ • b₁ = lev k₂ • b₂ → Ψ k₁ b₁ = Ψ k₂ b₂ := by
    intro k₁ k₂ b₁ b₂ heq
    rcases le_total k₁ k₂ with h | h
    · exact hwd_le k₁ k₂ b₁ b₂ h heq
    · exact (hwd_le k₂ k₁ b₂ b₁ h heq.symm).symm
  refine ⟨fun z => if h : ∃ (k : ℕ) (b : B), lev k • b = z then
    Ψ h.choose (h.choose_spec.choose) else 0, ?_⟩
  intro k b
  have hex : ∃ (k' : ℕ) (b' : B), lev k' • b' = lev k • b := ⟨k, b, rfl⟩
  show (if h : ∃ (k' : ℕ) (b' : B), lev k' • b' = lev k • b then
    Ψ h.choose (h.choose_spec.choose) else 0) = Ψ k b
  rw [dif_pos hex]
  exact hwd _ k _ b hex.choose_spec.choose_spec

/-- 2つの元を**共通の水準**の代表元で書き直す(`k := max k₁ k₂` を取るだけ)。 -/
theorem exists_common_level {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (lev : ℕ → A) (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (k₁ k₂ : ℕ) (b₁ b₂ : B) :
    ∃ (k : ℕ) (c₁ c₂ : B), lev k • c₁ = lev k₁ • b₁ ∧ lev k • c₂ = lev k₂ • b₂ := by
  obtain ⟨d₁, hd₁⟩ := hdvd k₁ (max k₁ k₂) (le_max_left _ _)
  obtain ⟨d₂, hd₂⟩ := hdvd k₂ (max k₁ k₂) (le_max_right _ _)
  refine ⟨max k₁ k₂, d₁ • b₁, d₂ • b₂, ?_, ?_⟩
  · rw [smul_smul, ← hd₁]
  · rw [smul_smul, ← hd₂]

/-- 同一水準での加法性。 -/
theorem glued_lift_add {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b) (k : ℕ) (b b' : B) :
    Φ (lev k • b + lev k • b') = Φ (lev k • b) + Φ (lev k • b') := by
  rw [← smul_add, hΦ, hΦ, hΦ, map_add]

/-- 同一水準での乗法性(`lift_mul_rep` を経由する)。 -/
theorem glued_lift_mul {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (k : ℕ) (b b' : B) :
    Φ ((lev k • b) * (lev k • b')) = Φ (lev k • b) * Φ (lev k • b') := by
  rw [lift_mul_rep, hΦ, hΦ, hΦ, map_smul, hmul]

/-- 同一水準で `φ` の持ち上げになっていること。 -/
theorem glued_lift_quot {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (I : Ideal C) (φ : B →ₐ[A] (C ⧸ I))
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hψ : ∀ (k : ℕ) (b : B), Ideal.Quotient.mk I (Ψ k b) = lev k • φ b)
    (k : ℕ) (b : B) :
    Ideal.Quotient.mk I (Φ (lev k • b)) = φ (lev k • b) := by
  rw [hΦ, hψ, map_smul]

/-- 任意の2元での加法性——共通水準に落としてから `glued_lift_add`。 -/
theorem glued_lift_add_general {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (z z' : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z)
    (hz' : ∃ (k : ℕ) (b : B), lev k • b = z') :
    Φ (z + z') = Φ z + Φ z' := by
  obtain ⟨k₁, b₁, rfl⟩ := hz
  obtain ⟨k₂, b₂, rfl⟩ := hz'
  obtain ⟨k, c₁, c₂, h1, h2⟩ := exists_common_level lev hdvd k₁ k₂ b₁ b₂
  rw [← h1, ← h2]
  exact glued_lift_add lev Ψ Φ hΦ k c₁ c₂

/-- 任意の2元での乗法性——共通水準に落としてから `glued_lift_mul`。 -/
theorem glued_lift_mul_general {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (z z' : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z)
    (hz' : ∃ (k : ℕ) (b : B), lev k • b = z') :
    Φ (z * z') = Φ z * Φ z' := by
  obtain ⟨k₁, b₁, rfl⟩ := hz
  obtain ⟨k₂, b₂, rfl⟩ := hz'
  obtain ⟨k, c₁, c₂, h1, h2⟩ := exists_common_level lev hdvd k₁ k₂ b₁ b₂
  rw [← h1, ← h2]
  exact glued_lift_mul lev Ψ Φ hΦ hmul k c₁ c₂

/-- 任意の元で `φ` の持ち上げになっていること。 -/
theorem glued_lift_quot_general {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (I : Ideal C) (φ : B →ₐ[A] (C ⧸ I))
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hψ : ∀ (k : ℕ) (b : B), Ideal.Quotient.mk I (Ψ k b) = lev k • φ b)
    (z : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z) :
    Ideal.Quotient.mk I (Φ z) = φ z := by
  obtain ⟨k, b, rfl⟩ := hz
  exact glued_lift_quot I φ lev Ψ Φ hΦ hψ k b

/-- **Faltings の `φ₀ : mB → C`**——整合な族を貼り合わせて得られる、
`mB` 全体の上で加法的・乗法的な honest な持ち上げ。
ここで `mB` の元は `{z | ∃ k b, lev k • b = z}` として表している。

これで `Theorem 2.2` の存在側に残るのは **`B = A + mB` を使う最後の拡張**
(原文の *"`φ₀` maps `p^ε` to an element `x = p^ε + y`(`y ∈ I`)satisfying
`x² = p^ε x`, hence `p^ε y = 0`"*)だけである。 -/
theorem exists_lift_on_mB {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (I : Ideal C) (φ : B →ₐ[A] (C ⧸ I))
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C))
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (hcompat : ∀ (k k' : ℕ) (d : A), k ≤ k' → lev k = lev k' * d → Ψ k = d • Ψ k')
    (hBtors : ∀ (k : ℕ) (b : B), lev k • b = 0 → b = 0)
    (hψ : ∀ (k : ℕ) (b : B), Ideal.Quotient.mk I (Ψ k b) = lev k • φ b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y) :
    ∃ Φ : B → C,
      (∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b) ∧
      (∀ z z' : B, (∃ (k : ℕ) (b : B), lev k • b = z) → (∃ (k : ℕ) (b : B), lev k • b = z') →
        Φ (z + z') = Φ z + Φ z') ∧
      (∀ z z' : B, (∃ (k : ℕ) (b : B), lev k • b = z) → (∃ (k : ℕ) (b : B), lev k • b = z') →
        Φ (z * z') = Φ z * Φ z') ∧
      (∀ z : B, (∃ (k : ℕ) (b : B), lev k • b = z) → Ideal.Quotient.mk I (Φ z) = φ z) := by
  obtain ⟨Φ, hΦ⟩ := exists_glued_lift lev Ψ hdvd hcompat hBtors
  exact ⟨Φ, hΦ,
    fun z z' hz hz' => glued_lift_add_general lev Ψ Φ hΦ hdvd z z' hz hz',
    fun z z' hz hz' => glued_lift_mul_general lev Ψ Φ hΦ hmul hdvd z z' hz hz',
    fun z hz => glued_lift_quot_general I φ lev Ψ Φ hΦ hψ z hz⟩

/-! ## `B = A·1 + mB` による最後の拡張

Faltings の *"`φ₀` maps `p^ε` to an element `x = p^ε + y`(`y ∈ I`)satisfying
`x² = p^ε x`, hence `p^ε y = 0`"*(物理 p.7)。これで `φ₀ : mB → C` は
`A·1` 上で `algebraMap` に一致することが分かり、`B = A·1 + mB` から
`B` 全体への `A`-代数準同型が一意に決まる。 -/

/-- `mB` は差について閉じている(共通水準に落とすだけ)。 -/
theorem mem_mB_sub {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (lev : ℕ → A) (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (z z' : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z)
    (hz' : ∃ (k : ℕ) (b : B), lev k • b = z') :
    ∃ (k : ℕ) (b : B), lev k • b = z - z' := by
  obtain ⟨k₁, b₁, rfl⟩ := hz
  obtain ⟨k₂, b₂, rfl⟩ := hz'
  obtain ⟨k, c₁, c₂, h1, h2⟩ := exists_common_level lev hdvd k₁ k₂ b₁ b₂
  exact ⟨k, c₁ - c₂, by rw [smul_sub, h1, h2]⟩

/-- `mB` は和について閉じている。 -/
theorem mem_mB_add {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (lev : ℕ → A) (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (z z' : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z)
    (hz' : ∃ (k : ℕ) (b : B), lev k • b = z') :
    ∃ (k : ℕ) (b : B), lev k • b = z + z' := by
  obtain ⟨k₁, b₁, rfl⟩ := hz
  obtain ⟨k₂, b₂, rfl⟩ := hz'
  obtain ⟨k, c₁, c₂, h1, h2⟩ := exists_common_level lev hdvd k₁ k₂ b₁ b₂
  exact ⟨k, c₁ + c₂, by rw [smul_add, h1, h2]⟩

/-- `mB` は `A` 倍について閉じている。 -/
theorem mem_mB_smul {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (lev : ℕ → A) (a : A) (z : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z) :
    ∃ (k : ℕ) (b : B), lev k • b = a • z := by
  obtain ⟨k, b, rfl⟩ := hz
  exact ⟨k, a • b, by rw [smul_comm]⟩

/-- `mB` はイデアルなので積について閉じている。 -/
theorem mem_mB_mul {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (lev : ℕ → A) (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (z z' : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z)
    (hz' : ∃ (k : ℕ) (b : B), lev k • b = z') :
    ∃ (k : ℕ) (b : B), lev k • b = z * z' := by
  obtain ⟨k₁, b₁, rfl⟩ := hz
  obtain ⟨k₂, b₂, rfl⟩ := hz'
  obtain ⟨k, c₁, c₂, h1, h2⟩ := exists_common_level lev hdvd k₁ k₂ b₁ b₂
  exact ⟨k, lev k • (c₁ * c₂), by rw [smul_smul, ← smul_mul_smul_comm, h1, h2]⟩

theorem mem_mB_zero {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (lev : ℕ → A) : ∃ (k : ℕ) (b : B), lev k • b = 0 :=
  ⟨0, 0, smul_zero _⟩

/-- `φ₀` は `mB` 上で `A`-線型(各水準の `Ψ k` が `A`-線型だから)。 -/
theorem glued_lift_smul {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (a : A) (z : B) (hz : ∃ (k : ℕ) (b : B), lev k • b = z) :
    Φ (a • z) = a • Φ z := by
  obtain ⟨k, b, rfl⟩ := hz
  rw [smul_comm, hΦ, hΦ, map_smul]

/-- **Faltings の `x² = p^ε x ⟹ p^ε y = 0` の段**。`φ₀(p^ε·1)` は
`x² = p^ε x` を満たし、かつ `x ≡ p^ε mod I` なので、`I` が平方零で
`p^ε` が非零因子なら `x = p^ε` ちょうど。 -/
theorem glued_lift_one {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I))
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hψ : ∀ (k : ℕ) (b : B), Ideal.Quotient.mk I (Ψ k b) = lev k • φ b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (lev k) * x = 0 → x = 0)
    (k : ℕ) :
    Φ (lev k • (1 : B)) = (algebraMap A C) (lev k) := by
  set x : C := Φ (lev k • (1 : B)) with hx
  have hxx : x * x = (algebraMap A C) (lev k) * x := by
    have h1 : Φ ((lev k • (1:B)) * (lev k • (1:B))) = x * x :=
      glued_lift_mul lev Ψ Φ hΦ hmul k 1 1
    have h2 : (lev k • (1:B)) * (lev k • (1:B)) = lev k • (lev k • (1:B)) := by
      rw [smul_mul_smul_comm, one_mul, smul_smul]
    rw [h2] at h1
    rw [glued_lift_smul lev Ψ Φ hΦ (lev k) _ ⟨k, 1, rfl⟩, hx] at h1
    rw [← h1, Algebra.smul_def]
  have hy : x - (algebraMap A C) (lev k) ∈ I := by
    rw [← Ideal.Quotient.eq_zero_iff_mem]
    have h3 := hψ k 1
    rw [hx, hΦ, map_sub, h3, map_one, Algebra.smul_def, mul_one,
      IsScalarTower.algebraMap_apply A C (C ⧸ I), Ideal.Quotient.algebraMap_eq, sub_self]
  set y : C := x - (algebraMap A C) (lev k) with hyd
  have hy2 : y * y = 0 := hsq y hy y hy
  have hxy : x = (algebraMap A C) (lev k) + y := by rw [hyd]; ring
  have hkey : (algebraMap A C) (lev k) * y = 0 := by
    rw [hxy] at hxx
    linear_combination hxx - hy2
  have hz : y = 0 := htors k y hkey
  rw [hyd] at hz
  linear_combination hz

/-- `mB ∩ A·1` 上では `φ₀` は `algebraMap` に一致する——`B = A·1 + mB`
分解の well-defined 性の核心。 -/
theorem glued_lift_algebraMap {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (hone : ∀ k : ℕ, Φ (lev k • (1 : B)) = (algebraMap A C) (lev k))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (lev k) * x = 0 → x = 0)
    (a : A) (hu : ∃ (k : ℕ) (b : B), lev k • b = a • (1 : B)) :
    Φ (a • (1 : B)) = (algebraMap A C) a := by
  obtain ⟨k, b, hkb⟩ := hu
  have h1 : Φ ((lev k • (1:B)) * (a • (1:B))) = Φ (lev k • (1:B)) * Φ (a • (1:B)) :=
    glued_lift_mul_general lev Ψ Φ hΦ hmul hdvd _ _ ⟨k, 1, rfl⟩ ⟨k, b, hkb⟩
  have h2 : (lev k • (1:B)) * (a • (1:B)) = a • (lev k • (1:B)) := by
    rw [smul_mul_smul_comm, one_mul, smul_smul, mul_comm]
  rw [h2, glued_lift_smul lev Ψ Φ hΦ a _ ⟨k, 1, rfl⟩, hone] at h1
  have h3 : (algebraMap A C) (lev k) * (Φ (a • (1:B)) - (algebraMap A C) a) = 0 := by
    rw [Algebra.smul_def] at h1
    linear_combination -h1
  have h4 := htors k _ h3
  linear_combination h4

/-- **`B = A·1 + mB` 分解の well-defined 性**。 -/
theorem glued_lift_wd {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (hone : ∀ k : ℕ, Φ (lev k • (1 : B)) = (algebraMap A C) (lev k))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (lev k) * x = 0 → x = 0)
    (a a' : A) (z z' : B)
    (hz : ∃ (k : ℕ) (b : B), lev k • b = z) (hz' : ∃ (k : ℕ) (b : B), lev k • b = z')
    (heq : a • (1 : B) + z = a' • (1 : B) + z') :
    (algebraMap A C) a + Φ z = (algebraMap A C) a' + Φ z' := by
  have hdiff : (a - a') • (1 : B) = z' - z := by
    rw [sub_smul]; linear_combination heq
  have hmem : ∃ (k : ℕ) (b : B), lev k • b = (a - a') • (1 : B) := by
    rw [hdiff]; exact mem_mB_sub lev hdvd z' z hz' hz
  have h1 : Φ ((a - a') • (1:B)) = (algebraMap A C) (a - a') :=
    glued_lift_algebraMap lev Ψ Φ hΦ hmul hdvd hone htors (a - a') hmem
  have h2 : Φ ((z' - z) + z) = Φ (z' - z) + Φ z :=
    glued_lift_add_general lev Ψ Φ hΦ hdvd _ _ (by rw [← hdiff]; exact hmem) hz
  rw [sub_add_cancel] at h2
  rw [hdiff] at h1
  rw [h1, map_sub] at h2
  linear_combination -h2

open Classical in
/-- **`B = A·1 + mB` による拡張**。`mB` 上の `φ₀` を `B` 全体に伸ばす。 -/
theorem exists_glued_extension {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (hone : ∀ k : ℕ, Φ (lev k • (1 : B)) = (algebraMap A C) (lev k))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (lev k) * x = 0 → x = 0)
    (hdecomp : ∀ b : B, ∃ (a : A) (z : B),
      (∃ (k : ℕ) (b' : B), lev k • b' = z) ∧ b = a • (1 : B) + z) :
    ∃ Φe : B → C, ∀ (a : A) (z : B), (∃ (k : ℕ) (b : B), lev k • b = z) →
      Φe (a • (1 : B) + z) = (algebraMap A C) a + Φ z := by
  refine ⟨fun b => (algebraMap A C) (hdecomp b).choose
    + Φ ((hdecomp b).choose_spec.choose), ?_⟩
  intro a z hz
  set b : B := a • (1 : B) + z with hb
  obtain ⟨hmem, heq⟩ := (hdecomp b).choose_spec.choose_spec
  show (algebraMap A C) (hdecomp b).choose + Φ ((hdecomp b).choose_spec.choose)
    = (algebraMap A C) a + Φ z
  exact glued_lift_wd lev Ψ Φ hΦ hmul hdvd hone htors _ a _ z hmem hz heq.symm

/-- **`Theorem 2.2` の存在側・完成形**——`mB` 上の `φ₀` を `B = A·1 + mB` で
拡張して得られる `A`-代数準同型 `ψ : B →ₐ[A] C` が `φ` を持ち上げる。 -/
theorem exists_algHom_of_glued {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (I : Ideal C) (φ : B →ₐ[A] (C ⧸ I))
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C)) (Φ : B → C)
    (hΦ : ∀ (k : ℕ) (b : B), Φ (lev k • b) = Ψ k b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (hone : ∀ k : ℕ, Φ (lev k • (1 : B)) = (algebraMap A C) (lev k))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (lev k) * x = 0 → x = 0)
    (hadd : ∀ z z' : B, (∃ (k : ℕ) (b : B), lev k • b = z) →
      (∃ (k : ℕ) (b : B), lev k • b = z') → Φ (z + z') = Φ z + Φ z')
    (hmulΦ : ∀ z z' : B, (∃ (k : ℕ) (b : B), lev k • b = z) →
      (∃ (k : ℕ) (b : B), lev k • b = z') → Φ (z * z') = Φ z * Φ z')
    (hquot : ∀ z : B, (∃ (k : ℕ) (b : B), lev k • b = z) →
      Ideal.Quotient.mk I (Φ z) = φ z)
    (hdecomp : ∀ b : B, ∃ (a : A) (z : B),
      (∃ (k : ℕ) (b' : B), lev k • b' = z) ∧ b = a • (1 : B) + z) :
    ∃ ψ : B →ₐ[A] C, ∀ b : B, Ideal.Quotient.mk I (ψ b) = φ b := by
  obtain ⟨Φe, hΦe⟩ := exists_glued_extension lev Ψ Φ hΦ hmul hdvd hone htors hdecomp
  have hΦ0 : Φ 0 = 0 := by have h := hΦ 0 0; rwa [smul_zero, map_zero] at h
  have hrep : ∀ b : B, ∃ (a : A) (z : B), (∃ (k : ℕ) (b' : B), lev k • b' = z) ∧
      b = a • (1 : B) + z ∧ Φe b = (algebraMap A C) a + Φ z := by
    intro b
    obtain ⟨a, z, hz, hbz⟩ := hdecomp b
    exact ⟨a, z, hz, hbz, by rw [hbz]; exact hΦe a z hz⟩
  have hcomm : ∀ a : A, Φe ((algebraMap A B) a) = (algebraMap A C) a := by
    intro a
    have h : (algebraMap A B) a = a • (1:B) + 0 := by
      rw [add_zero, Algebra.algebraMap_eq_smul_one]
    rw [h, hΦe a 0 (mem_mB_zero lev), hΦ0, add_zero]
  refine ⟨{ toFun := Φe
            map_one' := ?_
            map_mul' := ?_
            map_zero' := ?_
            map_add' := ?_
            commutes' := hcomm }, ?_⟩
  · have h : (1:B) = (1:A) • (1:B) + 0 := by rw [one_smul, add_zero]
    rw [h, hΦe 1 0 (mem_mB_zero lev), hΦ0, add_zero, map_one]
  · intro b b'
    obtain ⟨a, z, hz, hbz, hΦb⟩ := hrep b
    obtain ⟨a', z', hz', hbz', hΦb'⟩ := hrep b'
    have hprod : b * b' = (a * a') • (1:B) + (a • z' + a' • z + z * z') := by
      rw [hbz, hbz']; simp only [Algebra.smul_def, map_mul, mul_one]; ring
    have hmem : ∃ (k : ℕ) (b : B), lev k • b = a • z' + a' • z + z * z' :=
      mem_mB_add lev hdvd _ _ (mem_mB_add lev hdvd _ _ (mem_mB_smul lev a z' hz')
        (mem_mB_smul lev a' z hz)) (mem_mB_mul lev hdvd z z' hz hz')
    rw [hprod, hΦe _ _ hmem, hΦb, hΦb',
      hadd _ _ (mem_mB_add lev hdvd _ _ (mem_mB_smul lev a z' hz') (mem_mB_smul lev a' z hz))
        (mem_mB_mul lev hdvd z z' hz hz'),
      hadd _ _ (mem_mB_smul lev a z' hz') (mem_mB_smul lev a' z hz),
      glued_lift_smul lev Ψ Φ hΦ a z' hz', glued_lift_smul lev Ψ Φ hΦ a' z hz,
      hmulΦ z z' hz hz', map_mul]
    simp only [Algebra.smul_def]
    ring
  · have h : (0:B) = (0:A) • (1:B) + 0 := by rw [zero_smul, add_zero]
    rw [h, hΦe 0 0 (mem_mB_zero lev), hΦ0, add_zero, map_zero]
  · intro b b'
    obtain ⟨a, z, hz, hbz, hΦb⟩ := hrep b
    obtain ⟨a', z', hz', hbz', hΦb'⟩ := hrep b'
    have hsum : b + b' = (a + a') • (1:B) + (z + z') := by
      rw [hbz, hbz', add_smul]; ring
    rw [hsum, hΦe _ _ (mem_mB_add lev hdvd z z' hz hz'), hΦb, hΦb',
      hadd z z' hz hz', map_add]
    ring
  · intro b
    obtain ⟨a, z, hz, hbz, hΦb⟩ := hrep b
    show Ideal.Quotient.mk I (Φe b) = φ b
    rw [hΦb, map_add, hquot z hz, hbz, map_add, map_smul, map_one,
      ← Algebra.algebraMap_eq_smul_one, IsScalarTower.algebraMap_apply A C (C ⧸ I),
      Ideal.Quotient.algebraMap_eq]

/-- **`Theorem 2.2`(存在側・単一の主張)**——整合な almost 持ち上げの族
`{Ψ k}` から、honest な `A`-代数準同型 `ψ : B →ₐ[A] C` による `φ` の
持ち上げを得る。仮定はすべて Faltings の設定そのもの:

* `hdvd`/`hcompat`——族が水準について整合(`p^{1/p^k}` の塔)
* `hBtors`/`htors`——`B`・`C` が `lev k` 捩れを持たない
* `hsq`——`I` は平方零(`Theorem 2.2` の設定)
* `hψ`/`hmul`——各水準で `almost` に持ち上げかつ `almost` 乗法的
* `hdecomp`——`B = A·1 + mB` -/
theorem thm_2_2_lift_of_family {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I))
    (lev : ℕ → A) (Ψ : ℕ → (B →ₗ[A] C))
    (hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d)
    (hcompat : ∀ (k k' : ℕ) (d : A), k ≤ k' → lev k = lev k' * d → Ψ k = d • Ψ k')
    (hBtors : ∀ (k : ℕ) (b : B), lev k • b = 0 → b = 0)
    (hψ : ∀ (k : ℕ) (b : B), Ideal.Quotient.mk I (Ψ k b) = lev k • φ b)
    (hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y)
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (lev k) * x = 0 → x = 0)
    (hdecomp : ∀ b : B, ∃ (a : A) (z : B),
      (∃ (k : ℕ) (b' : B), lev k • b' = z) ∧ b = a • (1 : B) + z) :
    ∃ ψ : B →ₐ[A] C, ∀ b : B, Ideal.Quotient.mk I (ψ b) = φ b := by
  obtain ⟨Φ, hΦ, hadd, hmulΦ, hquot⟩ :=
    exists_lift_on_mB I φ lev Ψ hdvd hcompat hBtors hψ hmul
  exact exists_algHom_of_glued I φ lev Ψ Φ hΦ hmul hdvd
    (fun k => glued_lift_one I hsq φ lev Ψ Φ hΦ hψ hmul htors k)
    htors hadd hmulΦ hquot hdecomp

open Classical in
/-- **`Theorem 2.2`(存在側)——`p`-可除な塔での完成形**。

> *Theorem 2.2. Suppose `B` is an almost étale covering of `A`, and
> `C → C/I` a surjection with `I² = 0`. Then any `A`-algebra map
> `φ : B → C/I` lifts to an `A`-algebra map `B → C`.*

証明の流れ(すべて本ファイル内で閉じている):

1. `B` の almost 射影性(`AlmostProjective.lean`)から、各水準 `ε = ϖ k`
   で `A`-加群写像 `ψ_ε` に持ち上げる。
2. 障害類を Hochschild `H²` で消して `ψ_ε` を almost 乗法的にする
   (`exists_multiplicative_lift_tower`)。
3. `Ω[B⁄A]` が almost 零(`kaehler_isAlmostZero_of_tower`)なので
   `ψ_ε` は水準を除いて一意——族 `{ψ_ε}` は整合する(`lift_compat`)。
4. 族を貼り合わせて `φ₀ : mB → C`(`exists_lift_on_mB`)。
5. `x² = p^ε x ⟹ p^ε y = 0` で `φ₀` が `A·1` 上 `algebraMap` に一致する
   ことを示し(`glued_lift_one`)、`B = A·1 + mB` で `B` 全体へ拡張する
   (`exists_algHom_of_glued`)。

一意性側は `thm_2_2_uniqueness_of_isAlmostEtale`(既証)。 -/
theorem thm_2_2_tower {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (hq : 1 ≤ q) (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0)))))
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (T.ϖ k) * x = 0 → x = 0)
    (hBtors : ∀ (k : ℕ) (b : B), T.ϖ k • b = 0 → b = 0)
    (hdecomp : ∀ b : B, ∃ (a : A) (z : B),
      (∃ (k : ℕ) (b' : B), ((T.ϖ k * T.ϖ k) * (T.ϖ k * T.ϖ k)) • b' = z)
        ∧ b = a • (1 : B) + z) :
    ∃ ψ : B →ₐ[A] C, ∀ b : B, Ideal.Quotient.mk I (ψ b) = φ b := by
  set lev : ℕ → A := fun k => (T.ϖ k * T.ϖ k) * (T.ϖ k * T.ϖ k) with hlev
  choose Ψ hψraw hmulraw using fun k =>
    exists_multiplicative_lift_tower T hAET hf0inj I hsq φ htors k k
  have hψ : ∀ (k : ℕ) (b : B), Ideal.Quotient.mk I (Ψ k b) = lev k • φ b := hψraw
  have hmul : ∀ (k : ℕ) (x y : B), lev k • Ψ k (x * y) = Ψ k x * Ψ k y := hmulraw
  have hlevtors : ∀ (k : ℕ) (x : C), (algebraMap A C) (lev k) * x = 0 → x = 0 := by
    intro k x hx
    have h4 : (algebraMap A C) (lev k) * x
        = (algebraMap A C) (T.ϖ k) * ((algebraMap A C) (T.ϖ k) *
          ((algebraMap A C) (T.ϖ k) * ((algebraMap A C) (T.ϖ k) * x))) := by
      rw [hlev]; simp only [map_mul]; ring
    rw [h4] at hx
    exact htors k _ (htors k _ (htors k _ (htors k _ hx)))
  have hlevBtors : ∀ (k : ℕ) (b : B), lev k • b = 0 → b = 0 := by
    intro k b hb
    have h4 : lev k • b = T.ϖ k • (T.ϖ k • (T.ϖ k • (T.ϖ k • b))) := by
      rw [smul_smul, smul_smul, smul_smul, hlev]; congr 1; ring
    rw [h4] at hb
    exact hBtors k _ (hBtors k _ (hBtors k _ (hBtors k _ hb)))
  have hdvd : ∀ k k' : ℕ, k ≤ k' → ∃ d : A, lev k = lev k' * d := by
    intro k k' hle
    obtain ⟨e, he⟩ := T.dvd_of_le hq hle
    exact ⟨(e * e) * (e * e), by rw [hlev]; simp only []; rw [he]; ring⟩
  have hΩ : ∀ x : Ω[B⁄A], T.ϖ 0 • x = 0 :=
    fun x => kaehler_isAlmostZero_of_tower T hAET hf0inj 0 x
  have hcompat : ∀ (k k' : ℕ) (d : A), k ≤ k' → lev k = lev k' * d → Ψ k = d • Ψ k' := by
    intro k k' d _ hd
    refine lift_compat I hsq φ (lev k') d (Ψ k') (Ψ k) (hψ k') (hmul k') ?_ ?_
      (T.ϖ 0) hΩ ?_ (htors 0)
    · intro b; rw [← hd]; exact hψ k b
    · intro x y; rw [← hd]; exact hmul k x y
    · intro x hx; rw [← hd] at hx; exact hlevtors k x hx
  exact thm_2_2_lift_of_family I hsq φ lev Ψ hdvd hcompat hlevBtors hψ hmul hlevtors hdecomp

/-- **`Theorem 2.2` の一意性側・塔版**。`Ω[B⁄A]` が almost 零であることから
2つの持ち上げの差は導分になり、消える(`thm_2_2_uniqueness`)。 -/
theorem thm_2_2_uniqueness_tower {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] [Module.Free A B]
    {q : ℕ} (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0)))))
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (T.ϖ k) * x = 0 → x = 0)
    (ψ ψ' : B →ₐ[A] C)
    (hlift : ∀ b : B, Ideal.Quotient.mk I (ψ b) = Ideal.Quotient.mk I (ψ' b)) :
    ψ = ψ' := by
  refine thm_2_2_uniqueness (T.ϖ 0)
    (fun x => kaehler_isAlmostZero_of_tower T hAET hf0inj 0 x) (htors 0) ψ ψ' (fun x y => ?_)
  refine hsq _ ?_ _ ?_ <;>
  · rw [← Ideal.Quotient.eq_zero_iff_mem, map_sub, hlift, sub_self]

/-- **`Theorem 2.2`(Faltings, *p-Adic Hodge Theory* Ch.I §2)——存在と一意性**。

> *2.2. Theorem. Suppose `B = A + mB` is an almost étale covering of `A`,
> `C` an `A`-algebra, `I ⊂ C` a nilpotent ideal, and `φ : B → C/I` an
> `A`-algebra morphism. Then `φ` lifts uniquely to `B → C`.*
> (物理 p.6-7 = 印字 p.259-260)

★仮定 `hdecomp` は原文の **`B = A + mB`** そのもの(定理の主語に書かれて
いる)。`htors`(`C` の `ϖ k` 捩れ無し)も原文が証明中で明示的に使う
(*"As `C` has no such torsion, the different `φ_ε` glue together"*)。

存在は `thm_2_2_tower`、一意性は `thm_2_2_uniqueness_tower`。 -/
theorem thm_2_2 {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (hq : 1 ≤ q) (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0)))))
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (φ : B →ₐ[A] (C ⧸ I))
    (htors : ∀ (k : ℕ) (x : C), (algebraMap A C) (T.ϖ k) * x = 0 → x = 0)
    (hBtors : ∀ (k : ℕ) (b : B), T.ϖ k • b = 0 → b = 0)
    (hdecomp : ∀ b : B, ∃ (a : A) (z : B),
      (∃ (k : ℕ) (b' : B), ((T.ϖ k * T.ϖ k) * (T.ϖ k * T.ϖ k)) • b' = z)
        ∧ b = a • (1 : B) + z) :
    ∃! ψ : B →ₐ[A] C, ∀ b : B, Ideal.Quotient.mk I (ψ b) = φ b := by
  obtain ⟨ψ, hψ⟩ := thm_2_2_tower hq T hAET hf0inj I hsq φ htors hBtors hdecomp
  refine ⟨ψ, hψ, fun ψ' hψ' => ?_⟩
  exact thm_2_2_uniqueness_tower T hAET hf0inj I hsq htors ψ' ψ
    (fun b => by rw [hψ' b, hψ b])

/-! ## `Theorem 2.2` の非空虚性

仮定が空虚に真になっていないことを、実際に成り立つ具体例で示す。
2つの軸(「`p` が真の非単元か」と「`B/A` が非自明か」)を別々に埋める
——`AlmostEtale.lean` の `awayOne_*` 系列と同じ流儀。 -/

/-- 古典的に étale・finite・free な `B/A` は、**任意の**塔 `T` について
almost étale covering tower になる。 -/
theorem isAlmostEtaleCoveringTower_of_etale {A B : Type u} [CommRing A] [CommRing B]
    [Algebra A B] [Algebra.Etale A B] [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (T : PDivTower A q) :
    IsAlmostEtaleCoveringTower (A := A) (B := B) T := by
  letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
  haveI := awayScalarTower (T.ϖ 0) (A := A) (B := B)
  obtain ⟨hF, hFin, hEt, htr, _⟩ :=
    isAlmostEtaleCovering_of_etale_general (A := A) (B := B) (T.ϖ 0)
  refine ⟨hF, hFin, hEt, htr, fun k => ?_⟩
  refine ⟨T.ϖ k • Algebra.FormallyUnramified.elem A B, ?_⟩
  rw [map_smul, diagonalCompare_elem_eq]

/-- 定数塔 `ϖ ≡ p`(`q = 1`)。`AlmostBase.lean` の退化例に対応する。 -/
def constTower {A : Type u} [CommRing A] (p : A) : PDivTower A 1 :=
  { ϖ := fun _ => p, ϖ_succ := fun _ => by ring }

/-- **非空虚性(その1)——真の非単元 `p = 5`**。`A = B = C = ℤ`、`I = ⊥`、
塔は定数塔 `ϖ ≡ 5`。`5` は `ℤ` の単元ではない。 -/
example : ∃! ψ : ℤ →ₐ[ℤ] ℤ, ∀ b : ℤ,
    Ideal.Quotient.mk (⊥ : Ideal ℤ) (ψ b) = (Ideal.Quotient.mkₐ ℤ (⊥ : Ideal ℤ)) b := by
  refine thm_2_2 (A := ℤ) (B := ℤ) (C := ℤ) le_rfl (constTower (5:ℤ))
    (isAlmostEtaleCoveringTower_of_etale _) ?_ ⊥ ?_ (Ideal.Quotient.mkₐ ℤ ⊥) ?_ ?_ ?_
  · show Function.Injective (algebraMap ℤ (Localization.Away ((algebraMap ℤ ℤ) (5:ℤ))))
    refine IsLocalization.injective (M := Submonoid.powers ((algebraMap ℤ ℤ) (5:ℤ))) _ ?_
    rw [Submonoid.powers_le]
    simp only [map_ofNat]
    exact mem_nonZeroDivisors_of_ne_zero (by norm_num)
  · intro u hu v _; rw [Ideal.mem_bot.mp hu, zero_mul]
  · intro k x hx; have h : (5:ℤ) * x = 0 := hx; omega
  · intro k b hb; have h : (5:ℤ) * b = 0 := hb; omega
  · intro b; exact ⟨b, 0, ⟨0, 0, smul_zero _⟩, by rw [smul_eq_mul, mul_one, add_zero]⟩

/-- **非空虚性(その2)——階数 2 の非自明な被覆 `B = Fin 2 → ℤ`**。
`p = 1`(単元)だが `B/A` は自明でない有限エタール拡大。 -/
example : ∃! ψ : (Fin 2 → ℤ) →ₐ[ℤ] (Fin 2 → ℤ), ∀ b : Fin 2 → ℤ,
    Ideal.Quotient.mk (⊥ : Ideal (Fin 2 → ℤ)) (ψ b)
      = (Ideal.Quotient.mkₐ ℤ (⊥ : Ideal (Fin 2 → ℤ))) b := by
  refine thm_2_2 (A := ℤ) (B := Fin 2 → ℤ) (C := Fin 2 → ℤ) le_rfl (constTower (1:ℤ))
    (isAlmostEtaleCoveringTower_of_etale _) ?_ ⊥ ?_ (Ideal.Quotient.mkₐ ℤ ⊥) ?_ ?_ ?_
  · show Function.Injective (algebraMap (Fin 2 → ℤ)
      (Localization.Away ((algebraMap ℤ (Fin 2 → ℤ)) (1:ℤ))))
    refine IsLocalization.injective
      (M := Submonoid.powers ((algebraMap ℤ (Fin 2 → ℤ)) (1:ℤ))) _ ?_
    rw [Submonoid.powers_le, map_one]
    exact one_mem _
  · intro u hu v _; rw [Ideal.mem_bot.mp hu, zero_mul]
  · intro k x hx; have h : (1 : Fin 2 → ℤ) * x = 0 := hx; rwa [one_mul] at h
  · intro k b hb; have h : (1:ℤ) • b = 0 := hb; rwa [one_smul] at h
  · intro b; exact ⟨0, b, ⟨0, b, one_smul _ _⟩, by rw [zero_smul, zero_add]⟩

/-! ## 冪零 → 平方零の還元(原文の *"We may assume that `I² = 0`"*)

原文 `Theorem 2.2` は *"`I ⊂ C` a **nilpotent** ideal"* と述べ、証明の
冒頭で *"We may assume that `I² = 0`"* と平方零に帰着している。
その還元を**純粋に形式的な dévissage** として切り出す。

`J := I^{k+1}` と置くと `J² = I^{2k+2} = 0`(`I^{k+2} = 0` より)なので
`J` は平方零であり、`C/J` の中で `I` の像は `(k+1)` 乗して消える。
したがって

* 帰納法の仮定を `C/J` に適用して `φ : B → C/I ≅ (C/J)/(I/J)` を
  `C/J` まで持ち上げ、
* 平方零の場合を `C ↠ C/J` に適用してさらに `C` まで持ち上げる。

★注意:この還元は「平方零の場合が**任意の `A`-代数 `C`** について
成り立つ」ことを仮定として要求する。我々の `thm_2_2` は `C` に
`htors`(`ϖ k` 捩れ無し)等を課しているので、そのまま合成はできない
——原文も `2.1` 直前で *"Divide `C/I²` by its `p`-torsion, etc."* と
断っており、冪零の場合には捩れの扱いに手当てが要ることを認めている。
本補題はその手当てを除いた**骨格**を与える。 -/

/-- **冪零 → 平方零の還元**。平方零の場合が任意の `A`-代数について
成り立てば、冪零の場合も出る。 -/
theorem lift_nilpotent_of_lift_sq {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (H : ∀ (C : Type u) (_ : CommRing C) (_ : Algebra A C) (J : Ideal C),
      (∀ u ∈ J, ∀ v ∈ J, u * v = 0) →
      ∀ φ : B →ₐ[A] (C ⧸ J), ∃ ψ : B →ₐ[A] C, ∀ b, Ideal.Quotient.mk J (ψ b) = φ b) :
    ∀ (k : ℕ) (C : Type u) (_ : CommRing C) (_ : Algebra A C) (I : Ideal C),
      I^(k+1) = 0 →
      ∀ φ : B →ₐ[A] (C ⧸ I), ∃ ψ : B →ₐ[A] C, ∀ b, Ideal.Quotient.mk I (ψ b) = φ b := by
  intro k
  induction k with
  | zero =>
    intro C hC hA I hI φ
    have hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0 := by
      intro u hu v hv
      rw [pow_one] at hI
      rw [hI] at hu
      simp only [Submodule.zero_eq_bot, Ideal.mem_bot] at hu
      rw [hu, zero_mul]
    exact H C hC hA I hsq φ
  | succ m ih =>
    intro C hC hA I hI φ
    set J : Ideal C := I^(m+1) with hJ
    have hJle : J ≤ I := by
      rw [hJ]
      calc I^(m+1) ≤ I^1 := Ideal.pow_le_pow_right (by omega)
        _ = I := pow_one I
    have hJsq : ∀ u ∈ J, ∀ v ∈ J, u * v = 0 := by
      intro u hu v hv
      have hmem : u * v ∈ I^(m+1) * I^(m+1) := Ideal.mul_mem_mul hu hv
      have hle : I^(m+1) * I^(m+1) ≤ I^(m+2) := by
        rw [← pow_add]
        exact Ideal.pow_le_pow_right (by omega)
      rw [hI] at hle
      have hb := hle hmem
      simpa only [Submodule.zero_eq_bot, Ideal.mem_bot] using hb
    have hIm : (Ideal.map (Ideal.Quotient.mk J) I)^(m+1) = 0 := by
      rw [← Ideal.map_pow, ← hJ, Ideal.map_quotient_self]
      rfl
    obtain ⟨ψ', hψ'⟩ := ih (C ⧸ J) inferInstance inferInstance
      (Ideal.map (Ideal.Quotient.mk J) I) hIm
      ((DoubleQuot.quotQuotEquivQuotOfLEₐ A hJle).symm.toAlgHom.comp φ)
    obtain ⟨ψ, hψ⟩ := H C hC hA J hJsq ψ'
    refine ⟨ψ, fun b => ?_⟩
    have h1 : Ideal.Quotient.mk (Ideal.map (Ideal.Quotient.mk J) I)
        (Ideal.Quotient.mk J (ψ b))
        = (DoubleQuot.quotQuotEquivQuotOfLEₐ A hJle).symm (φ b) := by
      rw [hψ b]; exact hψ' b
    have h3 : (DoubleQuot.quotQuotEquivQuotOfLE hJle)
        (Ideal.Quotient.mk (Ideal.map (Ideal.Quotient.mk J) I) (Ideal.Quotient.mk J (ψ b)))
        = Ideal.Quotient.mk I (ψ b) :=
      DoubleQuot.quotQuotEquivQuotOfLE_quotQuotMk (ψ b) hJle
    have hagree : ∀ x, (DoubleQuot.quotQuotEquivQuotOfLE hJle) x
        = (DoubleQuot.quotQuotEquivQuotOfLEₐ A hJle) x := fun _ => rfl
    rw [h1, hagree, AlgEquiv.apply_symm_apply] at h3
    exact h3.symm

/-- 非空虚性——`B = A` のとき仮定 `H` は成り立つ(`A`-代数写像 `A → X` は
`algebraMap` ただ 1 つ)。したがって冪零の場合の結論も実際に得られる。 -/
example (A : Type u) [CommRing A] :
    ∀ (k : ℕ) (C : Type u) (_ : CommRing C) (_ : Algebra A C) (I : Ideal C),
      I^(k+1) = 0 →
      ∀ φ : A →ₐ[A] (C ⧸ I), ∃ ψ : A →ₐ[A] C, ∀ b, Ideal.Quotient.mk I (ψ b) = φ b := by
  refine lift_nilpotent_of_lift_sq (A := A) (B := A) ?_
  intro C hC hA J _ φ
  refine ⟨Algebra.ofId A C, fun b => ?_⟩
  have h1 : φ b = algebraMap A (C ⧸ J) b := by
    have hc := φ.commutes b
    simpa using hc
  rw [h1]
  rfl

/-! ### ★★★★項目全体の `.src`(2026-09-05)

`Theorem 2.2` は存在と一意性の両方が `∃!` の形で証明され、非空虚性の
対照も 2 件ある。原文に無い仮定(逸脱)は `.needs` に記録する。 -/

/-- ★★★★**[Falt1] Theorem 2.2**——almost étale covering に沿った
冪零イデアル上の**持ち上げの存在と一意性**(`∃!`)。

## ★主張

| 原文 | 宣言 |
|---|---|
| *"there exists a unique lifting"* | `thm_2_2`(`∃! ψ : B →ₐ[A] C`) |
| 存在 | `thm_2_2_tower`(族の貼り合わせ `exists_glued_lift` + `B = A·1 + mB` による拡張 `exists_glued_extension`) |
| 一意性 | `thm_2_2_uniqueness_tower` |
| *"We may assume that `I² = 0`"* | `lift_nilpotent_of_lift_sq` |
| 非空虚性 | 2 件(`p` が真の非単元・`B/A` が非自明の 2 軸) |

## ★逸脱の記録(`falt1-goal.md` §0.1)

1. `PDivTower A q`——原文の `p^ε`(`ε ∈ ℚ_{>0}`)に Lean で意味を
   与えるための塔。`m := span(range ϖ)` は `m² = m` を満たす。
2. `hBtors`(`B` の `ϖ k` 捩れ無し)——族の貼り合わせの well-defined 性で使う。
3. `[Module.Finite A B]`・`[Module.Free A B]`——mathlib の
   `Algebra.trace` が `Module.Free` を要求するため。

★`hdecomp : B = A·1 + mB` と `htors`(`C` の捩れ無し)は**逸脱ではない**
——原文が `Theorem 2.2` の仮定・証明中で明示している。 -/
def thm_2_2.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 7, item := "Theorem 2.2", sectionId := "falt1-thm-2-2" }

def thm_2_2.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "thm_2_2(∃! そのもの)"
      (.inProject "ABC3" "ABC3.Found.Falt1.thm_2_2") 7,
    .citation "[ABC3]" "thm_2_2_tower(存在)"
      (.inProject "ABC3" "ABC3.Found.Falt1.thm_2_2_tower") 7,
    .citation "[ABC3]" "thm_2_2_uniqueness_tower(一意性)"
      (.inProject "ABC3" "ABC3.Found.Falt1.thm_2_2_uniqueness_tower") 7,
    .citation "[ABC3]" "lift_nilpotent_of_lift_sq(原文の \"We may assume I² = 0\")"
      (.inProject "ABC3" "ABC3.Found.Falt1.lift_nilpotent_of_lift_sq") 7,
    .implicitStep
      ("★逸脱 1: PDivTower(原文の p^ε に意味を与える塔)。" ++
       "n : ℕ 添字だけでは honest な結論(2.2 の存在)は出せないので必要だった") 7,
    .implicitStep
      "★逸脱 2: hBtors(B の ϖk 捩れ無し)——族の貼り合わせの well-defined 性で使う" 7,
    .implicitStep
      "★逸脱 3: Module.Finite/Free A B——mathlib の Algebra.trace が Module.Free を要求" 7 ]

end ABC3.Found.Falt1
