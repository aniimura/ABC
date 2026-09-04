import ABC3.Found.PGC.LubinTateReciprocitySurjective
import Mathlib.FieldTheory.Galois.Profinite

/-!
# `Gal(K.closure/K.carrier)` はコンパクト、`{σ:σx=y}` は閉(`sorry` 無し)

節目(5)のさらに先(古典的Lubin-Tate理論が実際に主張する`Gal(K_π/K)
≅𝒪_K^×`、`K_π:=K(Λ_∞)`)へ向けた、**より単純な**布石。

`reciprocityMapLimit`の全射性を示す当初の計画(`K_π`を構成し、
両立する自己同型の無限列を`K_π`上へ「貼り合わせる」)は、より単純な
**コンパクト性論法**に置き換えられる見込みが分かった: `Gal(K.closure/
K.carrier)`は副有限群(profinite)ゆえコンパクトであり、各レベル`n`
での条件`C_n:={σ:σ(x_n)=u·x_n}`は(`reciprocityMap_surjective`から)
空でなく、(`fixingSubgroup`の開性・閉性から)**閉**であり、n跨ぎで
`C_{n+1}⊆C_n`となる(入れ子)——`IsCompact.nonempty_iInter_of_
sequence_nonempty_isCompact_isClosed`(mathlib)で交わりが空でないこと
が従い、その交わりの元`σ`が**すべての`n`で同時に**`σ(x_n)=u·x_n`を
満たす。`K_π`という新しい対象を構成する必要も、「両立する自己同型の
無限列の貼り合わせ」という難所も回避できる。

## 本ファイルで確立する部品

- `Gal(K.closure/K.carrier)`が`CompactSpace`であること
  (`InfiniteGalois.continuousMulEquivToLimit`——副有限群の標準的な
  射影極限との同相を経由)。
- `{σ:σ(x)=y}`が(適当な`σ₀`で`σ₀(x)=y`が成り立つとき)閉集合である
  こと——`(K.carrier⟮x⟯).fixingSubgroup`が開部分群(`Intermediate
  Field.fixingSubgroup_isOpen`、`FiniteDimensional`から)かつ
  (`Subgroup.isClosed_of_isOpen`から)閉部分群であることと、
  `{σ:σx=y}`がこの部分群の`σ₀`によるコセット(`Homeomorph.mulLeft
  σ₀`による像)に一致することを組み合わせる。
- `g∈(K.carrier⟮x⟯).fixingSubgroup ↔ g(x)=x`の言い換え(`x`だけを
  固定すれば`K.carrier⟮x⟯`全体を固定することが、`x`が生成元である
  ことと`IntermediateField.algHom_ext_of_eq_adjoin`から従う)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- `K.carrier`は`ℚ_[p]`の有限次拡大なので標数`0`。 -/
theorem charZero_carrier {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) : CharZero K.carrier :=
  charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective

/-- `K.closure/K.carrier`は(標数0なので分離的、既存の`IsAlgClosure.
normal`と合わせて)Galois拡大。 -/
theorem isGalois_carrier_closure {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    IsGalois K.carrier K.closure := by
  haveI := IsAlgClosure.normal K.carrier K.closure
  haveI := charZero_carrier K
  haveI : PerfectField K.carrier := PerfectField.ofCharZero
  haveI : Algebra.IsSeparable K.carrier K.closure := Algebra.IsAlgebraic.isSeparable_of_perfectField
  exact ⟨⟩

/-- **`Gal(K.closure/K.carrier)`はコンパクト**——副有限群としての
標準的な射影極限との同相(`InfiniteGalois.continuousMulEquivToLimit`)
を経由する。 -/
theorem compactSpace_algEquiv {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    CompactSpace (K.closure ≃ₐ[K.carrier] K.closure) := by
  haveI := isGalois_carrier_closure K
  exact (InfiniteGalois.continuousMulEquivToLimit K.carrier K.closure).toHomeomorph.symm.compactSpace

/-- **`x`を固定する自己同型は`K.carrier⟮x⟯`全体を固定する**——`x`が
`K.carrier⟮x⟯`の生成元であることと`IntermediateField.algHom_ext_of_
eq_adjoin`(生成集合上で一致する2つの代数準同型は全体で一致する)から。 -/
theorem mem_fixingSubgroup_adjoin_of_eq_self
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (g : K.closure ≃ₐ[K.carrier] K.closure) (hg : g x = x) :
    g ∈ (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup := by
  rw [IntermediateField.mem_fixingSubgroup_iff]
  have hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) :=
    IntermediateField.mem_adjoin_simple_self K.carrier x
  set φ₁ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) →ₐ[K.carrier] K.closure :=
    g.toAlgHom.comp (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val with hφ1
  set φ₂ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) →ₐ[K.carrier] K.closure :=
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val with hφ2
  have heq : φ₁ = φ₂ := by
    apply IntermediateField.algHom_ext_of_eq_adjoin K.carrier rfl
    intro w hw
    simp only [Set.mem_singleton_iff] at hw
    show φ₁ ⟨w, hw ▸ hmem⟩ = φ₂ ⟨w, hw ▸ hmem⟩
    show g w = w
    rw [hw]; exact hg
  intro z hz
  have h := congrFun (congrArg DFunLike.coe heq) ⟨z, hz⟩
  show φ₁ ⟨z, hz⟩ = z
  rw [h]
  rfl

/-- **`{σ:σ(x)=y}`は閉集合**(`σ₀(x)=y`となる`σ₀`が1つ存在すれば)
——この集合は`(K.carrier⟮x⟯).fixingSubgroup`の`σ₀`による左コセット
(`Homeomorph.mulLeft σ₀`の像)に一致し、`fixingSubgroup`は開かつ閉
(`IntermediateField.fixingSubgroup_isOpen`+`Subgroup.isClosed_of_
isOpen`)、コセットは同相写像の像なので閉性を保つ。 -/
theorem isClosed_algEquiv_eq
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) [IsGalois K.carrier K.closure]
    (x y : K.closure) [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ₀ : K.closure ≃ₐ[K.carrier] K.closure) (hσ₀ : σ₀ x = y) :
    IsClosed {σ : K.closure ≃ₐ[K.carrier] K.closure | σ x = y} := by
  have heq : {σ : K.closure ≃ₐ[K.carrier] K.closure | σ x = y} =
      (Homeomorph.mulLeft σ₀) '' (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup := by
    ext τ
    simp only [Set.mem_setOf_eq, Set.mem_image, SetLike.mem_coe]
    constructor
    · intro hτ
      refine ⟨σ₀⁻¹ * τ, mem_fixingSubgroup_adjoin_of_eq_self K x (σ₀⁻¹ * τ) ?_, ?_⟩
      · show σ₀⁻¹ (τ x) = x
        rw [hτ, ← hσ₀]
        exact AlgEquiv.symm_apply_apply σ₀ x
      · show σ₀ * (σ₀⁻¹ * τ) = τ
        rw [← mul_assoc, mul_inv_cancel, one_mul]
    · rintro ⟨g, hg, rfl⟩
      rw [IntermediateField.mem_fixingSubgroup_iff] at hg
      show σ₀ (g x) = y
      rw [hg x (IntermediateField.mem_adjoin_simple_self K.carrier x), hσ₀]
  rw [heq]
  have hopen : IsOpen ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup :
      Set (K.closure ≃ₐ[K.carrier] K.closure)) :=
    IntermediateField.fixingSubgroup_isOpen (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  have hclosed : IsClosed ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup :
      Set (K.closure ≃ₐ[K.carrier] K.closure)) :=
    Subgroup.isClosed_of_isOpen _ hopen
  exact (Homeomorph.mulLeft σ₀).isClosed_image.mpr hclosed

end ABC3.Found.PGC
