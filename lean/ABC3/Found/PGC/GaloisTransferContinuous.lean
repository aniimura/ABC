import ABC3.Found.PGC.GaloisTransfer
import ABC3.Found.PGC.UnramifiedSubextension

/-!
# 体の同型から誘導される `Γ_K ≅ Γ_{K'}` は**連続**である

`Found/PGC/GaloisTransfer.lean` は `α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier` から
**裸の**群同型 `galMulEquivOf α ᾱ hfwd : K.absGal ≃* K'.absGal` を作り、
延長 `ᾱ` の取り方によらない(内部自己同型で繋がる)ことまで示していたが、
同ファイルの docstring は

> continuity(`ContinuousMulEquiv` にするための Krull 位相の連続性)も未確認

と記録していた。**本ファイルはその穴を埋める。**

## 証明

群準同型なので `1` での連続性でよい(`continuous_of_continuousAt_one`)。
`s ∈ 𝓝 1` を取ると `krullTopology_mem_nhds_one_iff` から有限次中間体 `E'`
で `E'.fixingSubgroup ⊆ s` なるものが取れる。標数 0 の原始元定理
(`exists_adjoin_eq_of_finiteDimensional`)で `E' = K'⟮y⟯` と書き、
`x := ᾱ⁻¹ y` として `K⟮x⟯.fixingSubgroup` を取ればよい:

`g` が `x` を固定すれば `(ᾱ∘g∘ᾱ⁻¹) y = ᾱ (g x) = ᾱ x = y` なので、
`ᾱ∘g∘ᾱ⁻¹` は `K'⟮y⟯` を固定する(生成元を固定する自己同型は生成体を固定する
——`mem_fixingSubgroup_adjoin_simple`)。

★逆写像も同じ補題を `α.symm`・`ᾱ.symm` に適用するだけで済む
(`galMulEquivOf` の `invFun` がちょうどその形に定義されているため、
証明項は**そのまま**通る)。

## これで何が言えるか

`Skeleton/PGC/Section4.lean::theorem_4_2` の「自然な射」を構成するのに
必要な部品は、これで

* 延長 `extendToClosure`(mathlib `IsAlgClosure.equivOfEquiv`)
* 共役 `galMulEquivOf`
* 選択非依存 `galMulEquivOf_indep`(外部同型としては一意)
* **連続性(本ファイル)**

まで揃った。残るのは `map_Gv`——`Interface.PGC.RamificationFiltration` に
「体の同型から誘導される共役が `Gv` を保つ」という**自然性の公理が無い**
ことに帰着する(`memory/pgc-ramification-naturality-gap.md`)。
本物の高次分岐群なら成り立つ性質だが、現在の `Interface` では課されていない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-- 生成元を固定する自己同型は、その生成体を(点ごとに)固定する。 -/
theorem mem_fixingSubgroup_adjoin_simple {F L : Type*} [Field F] [Field L] [Algebra F L]
    {y : L} {σ : L ≃ₐ[F] L} (h : σ y = y) :
    σ ∈ (IntermediateField.adjoin F ({y} : Set L)).fixingSubgroup := by
  have hle : Subgroup.closure ({σ} : Set (L ≃ₐ[F] L))
      ≤ (IntermediateField.adjoin F ({y} : Set L)).fixingSubgroup := by
    rw [← IntermediateField.le_iff_le, IntermediateField.adjoin_simple_le_iff,
      IntermediateField.mem_fixedField_iff]
    intro f hf
    induction hf using Subgroup.closure_induction with
    | mem g hg => rw [Set.mem_singleton_iff] at hg; rw [hg]; exact h
    | one => simp
    | mul a b _ _ ha hb => rw [AlgEquiv.mul_apply, hb, ha]
    | inv a _ ha =>
        have h1 : (a⁻¹ * a) y = a⁻¹ y := by rw [AlgEquiv.mul_apply, ha]
        rw [inv_mul_cancel, AlgEquiv.one_apply] at h1
        exact h1.symm
  exact hle (Subgroup.subset_closure rfl)

variable {p : ℕ} [Fact p.Prime]

/-- **★★★★★体の同型から誘導される `Γ_K ≅ Γ_{K'}` は連続**。 -/
theorem continuous_galMulEquivOf {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    Continuous (galMulEquivOf α ᾱ hfwd) := by
  refine continuous_of_continuousAt_one (galMulEquivOf α ᾱ hfwd).toMonoidHom ?_
  rw [ContinuousAt, map_one, Filter.tendsto_def]
  intro s hs
  obtain ⟨E', hE'fin, hE'sub⟩ := (krullTopology_mem_nhds_one_iff K'.carrier K'.closure s).mp hs
  haveI := hE'fin
  obtain ⟨y, hy⟩ := exists_adjoin_eq_of_finiteDimensional K' E'
  rw [krullTopology_mem_nhds_one_iff]
  refine ⟨IntermediateField.adjoin K.carrier ({ᾱ.symm y} : Set K.closure), inferInstance, ?_⟩
  intro g hg
  apply hE'sub
  rw [← hy]
  refine mem_fixingSubgroup_adjoin_simple ?_
  have hgx : g (ᾱ.symm y) = ᾱ.symm y :=
    (IntermediateField.mem_fixingSubgroup_iff _ g).mp hg _
      (IntermediateField.mem_adjoin_simple_self K.carrier (ᾱ.symm y))
  show (conjGalOfEquiv α ᾱ hfwd g) y = y
  rw [conjGalOfEquiv_apply]
  show ᾱ (g (ᾱ.symm y)) = y
  rw [hgx, RingEquiv.apply_symm_apply]

/-- 逆写像も連続——`galMulEquivOf` の `invFun` は
`conjGalOfEquiv α.symm ᾱ.symm _` そのものなので、上の補題を
`α.symm`・`ᾱ.symm` に適用した証明項がそのまま通る。 -/
theorem continuous_galMulEquivOf_symm {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    Continuous (galMulEquivOf α ᾱ hfwd).symm :=
  continuous_galMulEquivOf α.symm ᾱ.symm (hfwd_symm_of_hfwd α ᾱ hfwd)

/-- **★★★★★★体の同型から誘導される絶対 Galois 群の `ContinuousMulEquiv`**
(互換な延長 `ᾱ` を指定する版)。 -/
noncomputable def galContinuousMulEquivOf {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    ContinuousMulEquiv K.absGal K'.absGal where
  toMulEquiv := galMulEquivOf α ᾱ hfwd
  continuous_toFun := continuous_galMulEquivOf α ᾱ hfwd
  continuous_invFun := continuous_galMulEquivOf_symm α ᾱ hfwd

/-- 既定の延長 `extendToClosure α` を使った版。
`galMulEquiv_conj_indep` により、他の延長を選んでも内部自己同型のずれしか
生まない——すなわち**外部同型としては一意**。 -/
noncomputable def galContinuousMulEquiv {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) : ContinuousMulEquiv K.absGal K'.absGal :=
  galContinuousMulEquivOf α (extendToClosure α) (extendToClosure_algebraMap α)

@[simp] theorem galContinuousMulEquiv_apply {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (g : K.absGal) :
    galContinuousMulEquiv α g = galMulEquiv α g := rfl

end ABC3.Found.PGC
