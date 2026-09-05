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

end ABC3.Found.Falt1
