import ABC3.Skeleton.PGC.Section4

/-!
# 体の同型から絶対 Galois 群の同型へ(`sorry` 無し、局所類体論を経由しない)

Theorem 4.2(`Skeleton/PGC/Section4.lean`)の `implicitStep`:

> Q_p-代数同型 α : K ≃ K′ を代数閉包への同型 ᾱ : K̄ ≃ K̄′ へ延長し、共役で誘導される
> Γ_K ≅ Γ_K′ が ᾱ の取り方によらない外部同型類を定めることの構成

を実際に検討した(2026-09-04)。**「延長+共役+選択非依存」の部分は局所類体論を
一切経由しない**——一般の Galois 理論だけで閉じることを、本ファイルで実際に
sorry 無しで確認した。

## できたこと

1. `extendToClosure` —— 体の同型 `α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier` を
   代数閉包への同型 `K.closure ≃+* K'.closure` へ延長する(`IsAlgClosure.equivOfEquiv`、
   mathlib 既存)。
2. `galMulEquivOf` —— 延長 `ᾱ`(基礎体との両立性 `hfwd` 付き)から、共役
   `g ↦ ᾱ∘g∘ᾱ⁻¹` によって誘導される **絶対 Galois 群の(裸の)群同型**
   `K.absGal ≃* K'.absGal` を構成する。
3. **`galMulEquivOf_indep`**(★★本ファイルの核心)—— 2つの異なる延長 `ᾱ₁, ᾱ₂`
   (どちらも同じ `α` を延長する)から誘導される2つの群同型は、**内部自己同型で
   ちょうど繋がる**: ある `c : K'.absGal` があって
   `galMulEquivOf α ᾱ₁ g = c * galMulEquivOf α ᾱ₂ g * c⁻¹`。
   これは `FilteredGroup.Iso.setoid`(`Found/PGC/FilteredGroup.lean`)の同値関係
   ——「後合成による内部自己同型」——とちょうど同じ形である。

## まだ無いもの(★★別の穴、`memory/pgc-ramification-naturality-gap.md` 参照)

上の構成が実際に返すのは Γ_K の**裸の**群同型であって、`FilteredGroup.Iso`
(高次分岐群のフィルトレーション `Gv` を保つ同型)ではない。`Interface.PGC.
RamificationFiltration` の現在の定義は K ごとに完全に独立な `Gv` を許す抽象データ
であり、「K の同型から誘導される共役が `Gv` を保つ」という自然性の公理を持たない
——ここで作った `galMulEquivOf` を `FilteredGroup.Iso` へ持ち上げるには、
`RamificationFiltration` にこの自然性を追加で課すか、具体的な(本物の)
`RamificationFiltration` を構成する必要がある。continuity(`ContinuousMulEquiv`
にするための Krull 位相の連続性)も未確認。したがって本ファイルは
Theorem 4.2 の `implicitStep` を**完全には**解消しない——「局所類体論を経由しない
部分」だけを切り出して確立したものである。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ### 部品1: 体の同型を代数閉包へ延長する -/

/-- `α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier` を `K.closure ≃+* K'.closure` へ延長する。
`K.closure = AlgebraicClosure K.carrier` は自動的に `IsAlgClosure K.carrier K.closure`
インスタンスを持つ(`infer_instance` で実測)ので、mathlib の一般定理
`IsAlgClosure.equivOfEquiv` がそのまま使える。 -/
noncomputable def extendToClosure {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) : K.closure ≃+* K'.closure :=
  IsAlgClosure.equivOfEquiv K.closure K'.closure α.toRingEquiv

/-- 延長は基礎体上で `α` と一致する。 -/
theorem extendToClosure_algebraMap {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (x : K.carrier) :
    extendToClosure α (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x) :=
  IsAlgClosure.equivOfEquiv_algebraMap K.closure K'.closure α.toRingEquiv x

/-! ### 部品2: 任意の互換な延長からの共役(一般形)

`galMulEquivOf` は特定の `extendToClosure α` にではなく、`α` と両立する**任意の**
`ᾱ : K.closure ≃+* K'.closure` にパラメータ化する——`galMulEquivOf_indep`
(部品3)で「異なる `ᾱ` を選んでも内部自己同型のずれしか生まない」ことを示すため、
この一般形が必要になる。 -/

/-- `hfwd`(`ᾱ` が基礎体上で `α` と一致する)から、`ᾱ` が `K'.carrier` を固定する
逆方向の一致も導ける——`ᾱ.symm` の側から見ても同じ形の両立性になる、という事実。 -/
theorem hfwd_symm_of_hfwd {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    ∀ y : K'.carrier, ᾱ.symm (algebraMap K'.carrier K'.closure y)
      = algebraMap K.carrier K.closure (α.symm y) := by
  intro y
  have h := hfwd (α.symm y)
  rw [AlgEquiv.apply_symm_apply] at h
  rw [← h, RingEquiv.symm_apply_apply]

/-- `g : K.absGal` を `ᾱ` で共役して得られる `K'.closure` 上の環同型
`ᾱ ∘ g ∘ ᾱ⁻¹`。 -/
noncomputable def conjGalOf {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure) (g : K.absGal) : K'.closure ≃+* K'.closure :=
  ᾱ.symm.trans (g.toRingEquiv.trans ᾱ)

theorem conjGalOf_apply {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure) (g : K.absGal) (x : K'.closure) :
    conjGalOf α ᾱ g x = ᾱ (g.toRingEquiv (ᾱ.symm x)) := rfl

/-- `hfwd` があれば `conjGalOf` は実際に `K'.carrier` を固定する
——つまり `K'.absGal` の元に昇格できる。 -/
theorem conjGalOf_fixes {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g : K.absGal) (y : K'.carrier) :
    conjGalOf α ᾱ g (algebraMap K'.carrier K'.closure y) = algebraMap K'.carrier K'.closure y := by
  rw [conjGalOf_apply, hfwd_symm_of_hfwd α ᾱ hfwd]
  have hg : g.toRingEquiv (algebraMap K.carrier K.closure (α.symm y))
      = algebraMap K.carrier K.closure (α.symm y) := g.commutes (α.symm y)
  rw [hg, hfwd, AlgEquiv.apply_symm_apply]

/-- `conjGalOf` を `K'.absGal`(`K'.carrier` を固定する `AlgEquiv`)へ昇格したもの。 -/
noncomputable def conjGalOfEquiv {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g : K.absGal) : K'.absGal :=
  AlgEquiv.ofRingEquiv (f := conjGalOf α ᾱ g) (conjGalOf_fixes α ᾱ hfwd g)

theorem conjGalOfEquiv_apply {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g : K.absGal) (x : K'.closure) :
    conjGalOfEquiv α ᾱ hfwd g x = ᾱ (g.toRingEquiv (ᾱ.symm x)) := rfl

theorem conjGalOfEquiv_map_mul {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g1 g2 : K.absGal) :
    conjGalOfEquiv α ᾱ hfwd (g1 * g2) = conjGalOfEquiv α ᾱ hfwd g1 * conjGalOfEquiv α ᾱ hfwd g2 := by
  apply AlgEquiv.ext
  intro x
  show conjGalOfEquiv α ᾱ hfwd (g1 * g2) x
      = conjGalOfEquiv α ᾱ hfwd g1 (conjGalOfEquiv α ᾱ hfwd g2 x)
  rw [conjGalOfEquiv_apply, conjGalOfEquiv_apply, conjGalOfEquiv_apply]
  show ᾱ ((g1 * g2).toRingEquiv (ᾱ.symm x)) = _
  have heq : (g1 * g2).toRingEquiv (ᾱ.symm x) = g1.toRingEquiv (g2.toRingEquiv (ᾱ.symm x)) := rfl
  rw [heq]
  congr 2
  rw [RingEquiv.symm_apply_apply]

theorem conjGalOfEquiv_symm_conjGalOfEquiv {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g : K.absGal) :
    conjGalOfEquiv α.symm ᾱ.symm (hfwd_symm_of_hfwd α ᾱ hfwd) (conjGalOfEquiv α ᾱ hfwd g) = g := by
  apply AlgEquiv.ext
  intro x
  show ᾱ.symm ((conjGalOfEquiv α ᾱ hfwd g).toRingEquiv (ᾱ.symm.symm x)) = g x
  show ᾱ.symm (ᾱ (g.toRingEquiv (ᾱ.symm (ᾱ.symm.symm x)))) = g x
  rw [RingEquiv.symm_symm, RingEquiv.symm_apply_apply, RingEquiv.symm_apply_apply]
  rfl

theorem conjGalOfEquiv_conjGalOfEquiv_symm {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g : K'.absGal) :
    conjGalOfEquiv α ᾱ hfwd (conjGalOfEquiv α.symm ᾱ.symm (hfwd_symm_of_hfwd α ᾱ hfwd) g) = g := by
  apply AlgEquiv.ext
  intro x
  show ᾱ ((conjGalOfEquiv α.symm ᾱ.symm (hfwd_symm_of_hfwd α ᾱ hfwd) g).toRingEquiv (ᾱ.symm x)) = g x
  show ᾱ (ᾱ.symm (g.toRingEquiv (ᾱ.symm.symm (ᾱ.symm x)))) = g x
  rw [RingEquiv.symm_symm, RingEquiv.apply_symm_apply, RingEquiv.apply_symm_apply]
  rfl

theorem conjGalOfEquiv_map_one {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    conjGalOfEquiv α ᾱ hfwd (1 : K.absGal) = 1 := by
  apply AlgEquiv.ext
  intro x
  rw [conjGalOfEquiv_apply]
  show ᾱ (ᾱ.symm x) = x
  simp

/-- ★**互換な延長 `ᾱ` から誘導される絶対 Galois 群の(裸の)群同型**。 -/
noncomputable def galMulEquivOf {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    K.absGal ≃* K'.absGal where
  toFun := conjGalOfEquiv α ᾱ hfwd
  invFun := conjGalOfEquiv α.symm ᾱ.symm (hfwd_symm_of_hfwd α ᾱ hfwd)
  left_inv := conjGalOfEquiv_symm_conjGalOfEquiv α ᾱ hfwd
  right_inv := conjGalOfEquiv_conjGalOfEquiv_symm α ᾱ hfwd
  map_mul' := conjGalOfEquiv_map_mul α ᾱ hfwd

/-! ### 部品3: ★★★選択非依存(内部自己同型で繋がる) -/

theorem galMulEquivOf_indep_fix {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ1 ᾱ2 : K.closure ≃+* K'.closure)
    (hfwd1 : ∀ x : K.carrier, ᾱ1 (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (hfwd2 : ∀ x : K.carrier, ᾱ2 (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    ∀ x : K.carrier, (ᾱ1.trans ᾱ2.symm) (algebraMap K.carrier K.closure x)
      = algebraMap K.carrier K.closure x := by
  intro x
  show ᾱ2.symm (ᾱ1 (algebraMap K.carrier K.closure x)) = algebraMap K.carrier K.closure x
  rw [hfwd1]
  have := hfwd_symm_of_hfwd α ᾱ2 hfwd2 (α x)
  rwa [AlgEquiv.symm_apply_apply] at this

/-- ★★★**選択非依存**: 2つの異なる延長 `ᾱ₁, ᾱ₂`(どちらも同じ `α` を延長する)から
誘導される群同型は、`K'.absGal` の内部自己同型でちょうど繋がる。

`c := ᾱ₁∘ᾱ₂⁻¹`(`K.closure` 上、`K.carrier` を固定するので `K.absGal` の元)を
取り、`galMulEquivOf α ᾱ₂ hfwd2 c` がその内部自己同型を与える元になる。

これは Theorem 4.2 の原文が要求する「Φ は延長 ᾱ の取り方によらない外部同型類を
定める」という主張の核心を、局所類体論を一切経由せず捉えたものである
(`FilteredGroup.Iso.setoid` の同値関係——後合成による内部自己同型——と
ちょうど同じ形)。 -/
theorem galMulEquivOf_indep {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ1 ᾱ2 : K.closure ≃+* K'.closure)
    (hfwd1 : ∀ x : K.carrier, ᾱ1 (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (hfwd2 : ∀ x : K.carrier, ᾱ2 (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g : K.absGal) :
    ∃ c : K'.absGal, galMulEquivOf α ᾱ1 hfwd1 g = c * galMulEquivOf α ᾱ2 hfwd2 g * c⁻¹ := by
  set c0 : K.closure ≃+* K.closure := ᾱ1.trans ᾱ2.symm with hc0
  set c : K.absGal := AlgEquiv.ofRingEquiv (f := c0) (galMulEquivOf_indep_fix α ᾱ1 ᾱ2 hfwd1 hfwd2)
    with hc
  have hc_apply : ∀ z : K.closure, c z = ᾱ2.symm (ᾱ1 z) := fun z => rfl
  have hcinv_apply : ∀ z : K.closure, c⁻¹ z = ᾱ1.symm (ᾱ2 z) := by
    intro z
    have hcc : c (c⁻¹ z) = z := by
      show (c * c⁻¹) z = z
      rw [mul_inv_cancel]; rfl
    rw [hc_apply] at hcc
    have h2 : ᾱ1 (c⁻¹ z) = ᾱ2 z := by
      have := congrArg ᾱ2 hcc
      rwa [RingEquiv.apply_symm_apply] at this
    have := congrArg ᾱ1.symm h2
    rwa [RingEquiv.symm_apply_apply] at this
  refine ⟨galMulEquivOf α ᾱ2 hfwd2 c, ?_⟩
  apply AlgEquiv.ext
  intro x
  show ᾱ1 (g.toRingEquiv (ᾱ1.symm x))
      = (galMulEquivOf α ᾱ2 hfwd2 c * galMulEquivOf α ᾱ2 hfwd2 g
          * (galMulEquivOf α ᾱ2 hfwd2 c)⁻¹) x
  rw [← map_inv, ← map_mul, ← map_mul]
  show ᾱ1 (g.toRingEquiv (ᾱ1.symm x)) = ᾱ2 ((c * g * c⁻¹).toRingEquiv (ᾱ2.symm x))
  have hmul : (c * g * c⁻¹).toRingEquiv (ᾱ2.symm x) = c (g (c⁻¹ (ᾱ2.symm x))) := rfl
  rw [hmul, hcinv_apply, RingEquiv.apply_symm_apply, hc_apply, RingEquiv.apply_symm_apply]
  rfl

/-! ### まとめ: `extendToClosure` を既定の選択に取った版 -/

/-- **体の同型から誘導される絶対 Galois 群の同型**(既定の延長 `extendToClosure α` を使う)。 -/
noncomputable def galMulEquiv {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    K.absGal ≃* K'.absGal :=
  galMulEquivOf α (extendToClosure α) (extendToClosure_algebraMap α)

/-- 上の選択非依存性を `galMulEquiv`(既定の選択)を基準にして書き直したもの:
**任意の**互換な延長 `ᾱ` から誘導される同型は、`galMulEquiv α` と内部自己同型で
ちょうど繋がる。 -/
theorem galMulEquiv_conj_indep {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (g : K.absGal) :
    ∃ c : K'.absGal, galMulEquivOf α ᾱ hfwd g = c * galMulEquiv α g * c⁻¹ :=
  galMulEquivOf_indep α ᾱ (extendToClosure α) hfwd (extendToClosure_algebraMap α) g

end ABC3.Found.PGC
