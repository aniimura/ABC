import ABC3.Found.Falt1.AlmostEtale

/-!
# [Falt1] almost projectivity と almost lifting——`Theorem 2.2` の第1段(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、物理 p.6-7
(印字 p.259-260)。

`Theorem 2.2`(冪零イデアルに沿った lifting)の証明の第1段は
「`B` の almost 射影性から、任意の `ε>0` で `p^ε φ` を `A` 加群写像として
持ち上げる」である(物理p.7、260dpi 目視)。本ファイルはその段を閉じる。

## 内容

**almost 射影性は remark(iii) から出る**。`Definition 2.1` 条件(iii)の
witness `w = Σᵢ xᵢ⊗yᵢ` について remark(iii)
(`remark_iii_trace_identity`、既証)は

    p^n·b = Σᵢ tr_{B/A}(b·xᵢ)·yᵢ   (任意の b∈B)

を与える——これはまさに **`{xᵢ}` と `{tr(·xᵢ)}` が almost dual basis を
なす**という主張であり、`p^n·id_B` が有限自由 `A`-加群を経由して分解する
(`almost_projective_factor`)ことに他ならない。

そこから **almost lifting**(`almost_lift_of_surjective`)は純粋な加群論
1行で出る:`f : B → A^ι`・`g : A^ι → B` で `g∘f = p^n·id` なら、
自由加群 `A^ι` は射影的なので `φ∘g` は持ち上がり(`Module.
projective_lifting_property`)、その持ち上げに `f` を合成すれば
`π∘ψ = p^n·φ` が得られる。

これは今セッションで確立した **「almost split ⟹ almost な性質」**
(`ext_smul_eq_zero_of_almost_split`・`transfer_groupCohomology_smul_eq_zero`)
と同じパターンの3度目の適用である(`tools/lean-idioms.md` #54)。

## `Theorem 2.2` の残り

第1段が閉じたので、残るのは第2段以降——持ち上げた `A`-加群写像 `φ_ε` の
積の非可換性の障害 `b₀(p^ε φ_ε(b₁b₂) - φ_ε(b₁)φ_ε(b₂))b₃` を
Hochschild cohomology `H²` の類として**実際に取り出す**変形理論である。
その入力(`p^n` が `Ext^{k+1}_{B⊗_AB}(B,M)` を零化する)は
`AlmostEtale.lean` の `hochschild_ext_almost_zero` で完成しているので、
残るのは「障害類を `Ext²` の元として構成する」枠組みの設計だけである。
-/

namespace ABC3.Found.Falt1

universe u

/-- **almost dual basis**(remark(iii)の言い換え)。`Definition 2.1`
条件(iii)の witness `w` を有限和 `Σᵢ xᵢ⊗yᵢ` に書けば、
`p^n·b = Σᵢ tr_{B/A}(b·xᵢ)·yᵢ` が任意の `b` について成り立つ。 -/
theorem almost_dual_basis {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Finite A B] [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p))) :
    ∃ S : Finset (B × B), ∀ b : B,
      (p ^ n) • b = ∑ i ∈ S, (Algebra.trace A B (b * i.1)) • i.2 := by
  obtain ⟨S, hS⟩ := TensorProduct.exists_finset w
  refine ⟨S, fun b => ?_⟩
  have h := remark_iii_trace_identity p hAE hf0inj n w b hw
  rw [← Algebra.smul_def] at h
  rw [h, hS, Finset.mul_sum, map_sum]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Algebra.TensorProduct.tmul_mul_tmul, one_mul, Tr1map_tmul]

/-- **`B` は `A` 上 almost 射影的**: `p^n·id_B` が有限自由 `A`-加群
`A^ι` を経由して分解する。`f b := (tr(b·xᵢ))ᵢ`・`g a := Σᵢ aᵢ·yᵢ` と
置けばよい(`almost_dual_basis`)。 -/
theorem almost_projective_factor {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Finite A B] [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p))) :
    ∃ (ι : Type u) (_ : Fintype ι) (f : B →ₗ[A] (ι → A)) (g : (ι → A) →ₗ[A] B),
      ∀ b : B, g (f b) = (p ^ n) • b := by
  obtain ⟨S, hS⟩ := almost_dual_basis p hAE hf0inj n w hw
  refine ⟨{x // x ∈ S}, inferInstance,
    { toFun := fun b i => Algebra.trace A B (b * (i : B × B).1)
      map_add' := by intro b b'; funext i; simp [add_mul]
      map_smul' := by
        intro a b; funext i
        simp only [RingHom.id_apply, Pi.smul_apply, smul_eq_mul]
        rw [smul_mul_assoc, map_smul, smul_eq_mul] },
    { toFun := fun a => ∑ i : {x // x ∈ S}, a i • (i : B × B).2
      map_add' := by
        intro a a'
        simp only [Pi.add_apply, add_smul]
        rw [Finset.sum_add_distrib]
      map_smul' := by
        intro c a
        simp only [Pi.smul_apply, RingHom.id_apply, smul_eq_mul, Finset.smul_sum, mul_smul] },
    fun b => ?_⟩
  show (∑ i : {x // x ∈ S}, (Algebra.trace A B (b * (i : B × B).1)) • (i : B × B).2) = (p ^ n) • b
  rw [Finset.sum_coe_sort S (fun i => (Algebra.trace A B (b * i.1)) • i.2)]
  exact (hS b).symm

/-- **almost lifting**(`Theorem 2.2` の第1段、純粋な加群論)。
`c·id_B` が有限自由加群を経由して分解するなら、任意の全射 `π : C ↠ D` と
任意の `φ : B → D` について、`π∘ψ = c·φ` を満たす `ψ : B → C` が取れる。
自由加群の射影性(`Module.projective_lifting_property`)一発。 -/
theorem almost_lift_of_surjective {A B C D : Type u} [CommRing A] [AddCommGroup B] [Module A B]
    [AddCommGroup C] [Module A C] [AddCommGroup D] [Module A D]
    (c : A)
    (hfact : ∃ (ι : Type u) (_ : Fintype ι) (f : B →ₗ[A] (ι → A)) (g : (ι → A) →ₗ[A] B),
      ∀ b : B, g (f b) = c • b)
    (π : C →ₗ[A] D) (hπ : Function.Surjective π) (φ : B →ₗ[A] D) :
    ∃ ψ : B →ₗ[A] C, ∀ b : B, π (ψ b) = c • φ b := by
  obtain ⟨ι, _, f, g, hgf⟩ := hfact
  obtain ⟨χ, hχ⟩ := Module.projective_lifting_property π (φ ∘ₗ g) hπ
  refine ⟨χ ∘ₗ f, fun b => ?_⟩
  have h := LinearMap.congr_fun hχ (f b)
  simp only [LinearMap.coe_comp, Function.comp_apply] at h ⊢
  rw [h, hgf b, map_smul]

/-- **`Theorem 2.2` の第1段そのもの**: `B` が `A` の almost étale covering
なら、任意の全射 `π : C ↠ D`(`A`-加群)と任意の `A`-加群写像 `φ : B → D`
について、`π∘ψ = p^n·φ` を満たす `ψ` が取れる。 -/
theorem almost_lift_of_isAlmostEtale {A B C D : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Finite A B] [Module.Free A B]
    [AddCommGroup C] [Module A C] [AddCommGroup D] [Module A D]
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
    (π : C →ₗ[A] D) (hπ : Function.Surjective π) (φ : B →ₗ[A] D) :
    ∃ ψ : B →ₗ[A] C, ∀ b : B, π (ψ b) = (p ^ n) • φ b :=
  almost_lift_of_surjective (p ^ n) (almost_projective_factor p hAE hf0inj n w hw) π hπ φ

end ABC3.Found.Falt1
