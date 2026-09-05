import ABC3.Skeleton.PGC.Section2
import ABC3.Found.PGC.QpNonAbelian

/-!
# [pGC] Proposition 2.2 の**旧**形は偽だった——`SMul` は「加群」ではない

原文 (pGC p.5):

> Suppose that we are given the following group-theoretic data: the topological group Γ_K,
> together with the indexed filtration Γ_K^v for all v > 0. Then the Γ_K-modules O[scr]_K[bar], and
> K[bar]∧ can be recovered group-theoretically from this group-theoretic data.

原文は「Γ_K-**加群**」と言っている。ところが旧形は `O_K̄`・`K̄^` を
**自由な型族**として抽象化したうえで、その作用を **`SMul`**
——`one_smul` も `mul_smul` も `smul_add` も要求しない、公理ゼロのクラス——
で受け取っていた。これは「加群」ではない。

## 反例(`prop_2_2_statement_false`、`sorry` 無し)

`Γ_{ℚ_p}` は非可換なので `c g₀ c⁻¹ ≠ g₀` なる `c`・`g₀` が取れる
(`Found/PGC/QpNonAbelian.lean::exists_not_commute_absGal`——`ℚ_p(p^{1/ℓ})` が
normal でないことから)。そこで

* `IntKbar K := ℤ`
* `g • n := if HEq g g₀ then 0 else n`(★`SMul` は公理が無いので合法)
* `α := c による共役`(位相群なので連続、`ContinuousMulEquiv`)

とすると、`φ : ℤ ≃+ ℤ` が `φ (g₀ • 1) = (c g₀ c⁻¹) • φ 1` を満たすことは不可能:
左辺 `= φ 0 = 0`、右辺 `= φ 1`(`c g₀ c⁻¹ ≠ g₀` なので潰れない)。
`φ 1 = 0` は `φ` が同型であることに反する。

★この反例は**非可換性を本質的に使う**——`Γ_K` が可換だと `c g₀ c⁻¹ = g₀` で
右辺も潰れてしまう。`Check/PGC/RefutationAttempts.lean` が
「共役による反証には非正規な部分群が要る」と記録していたのはこのこと。

## 修理

`SMul` を **`DistribMulAction`**(= 加法的自己同型として作用する、原文の
「Γ_K-加群」)に強めた。上の病的な作用は `MulAction` ですらない
(`g₀ • (g₀⁻¹ • n) = 0 ≠ n = (g₀ * g₀⁻¹) • n`)ので反例は塞がる。
`Prop 2.1` が使う `K.closure` への自然な作用は `DistribMulAction` を満たす
(`closureDistribMulAction`)。

★これで「落とした条件は、主張を偽にするか自明にするかのどちらかになる」例は
5 つ目(`InertiaDegeneracy`・`Theorem42Degenerate`・`Def32Degenerate`・
`Cor33Degenerate`・本件)。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC
open scoped Classical

variable {p : ℕ} [Fact p.Prime]

/-- **旧** `RecoverableAsAddModule`——作用が公理ゼロの `SMul` だった形。 -/
def RecoverableAsAddModuleOld (Obj : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (Obj K)] [∀ K, SMul K.absGal (Obj K)] : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : ContinuousMulEquiv K.absGal K'.absGal),
    ∃ φ : Obj K ≃+ Obj K', ∀ (g : K.absGal) (x : Obj K),
      φ (g • x) = (α.toMulEquiv g) • (φ x)

/-- 内部自己同型は連続な群同型。 -/
noncomputable def conjContinuousEquiv (K : PAdicLocalField p) (c : K.absGal) :
    ContinuousMulEquiv K.absGal K.absGal where
  toFun g := c * g * c⁻¹
  invFun g := c⁻¹ * g * c
  left_inv g := by group
  right_inv g := by group
  map_mul' a b := by group
  continuous_toFun := by fun_prop
  continuous_invFun := by fun_prop

@[simp] theorem conjContinuousEquiv_apply (K : PAdicLocalField p) (c g : K.absGal) :
    conjContinuousEquiv K c g = c * g * c⁻¹ := rfl

/-- 病的な「作用」——`g₀` だけがすべてを潰す。`SMul` には公理が無いので合法。 -/
@[implicit_reducible] noncomputable def badSMul (p : ℕ) [Fact p.Prime] (g₀ : (selfField p).absGal)
    (K : PAdicLocalField p) : SMul K.absGal ℤ :=
  ⟨fun g n => if HEq g g₀ then 0 else n⟩

noncomputable def topFiltration' (p : ℕ) [Fact p.Prime] : RamificationFiltration p where
  Gv _ _ := ⊤
  isClosed _ _ := by simp
  isNormal _ _ := inferInstance
  antitone _ := fun _ _ _ => le_rfl

/-- **★★★★★★[pGC] Proposition 2.2 の旧形は偽**——作用が `SMul`(公理ゼロ)
だったので、非可換性を使った病的な作用が反例になる。 -/
theorem prop_2_2_statement_false (p : ℕ) [Fact p.Prime] :
    ¬ (∀ (_RF : RamificationFiltration p)
        (IntKbar CompKbar : PAdicLocalField p → Type)
        (i1 : ∀ K, AddCommGroup (IntKbar K)) (i2 : ∀ K, SMul K.absGal (IntKbar K))
        (i3 : ∀ K, AddCommGroup (CompKbar K)) (i4 : ∀ K, SMul K.absGal (CompKbar K)),
        @RecoverableAsAddModuleOld p _ IntKbar i1 i2
          ∧ @RecoverableAsAddModuleOld p _ CompKbar i3 i4) := by
  intro h
  obtain ⟨c, g₀, hne⟩ := exists_not_commute_absGal p
  letI : SMul (selfField p).absGal ℤ := badSMul p g₀ (selfField p)
  obtain ⟨key, -⟩ := h (topFiltration' p) (fun _ => ℤ) (fun _ => ℤ)
      (fun _ => inferInstance) (badSMul p g₀) (fun _ => inferInstance) (badSMul p g₀)
  obtain ⟨φ, hφ⟩ := key (K := selfField p) (K' := selfField p) (conjContinuousEquiv (selfField p) c)
  have h1 := hφ g₀ 1
  have hL : (g₀ : (selfField p).absGal) • (1 : ℤ) = 0 := by
    show (if HEq g₀ g₀ then (0:ℤ) else 1) = 0
    simp
  have hR : ((conjContinuousEquiv (selfField p) c).toMulEquiv g₀) • (φ 1) = φ 1 := by
    show (if HEq (c * g₀ * c⁻¹) g₀ then (0:ℤ) else φ 1) = φ 1
    rw [if_neg]
    intro hh
    exact hne (eq_of_heq hh)
  rw [hL, hR, map_zero] at h1
  exact one_ne_zero (φ.injective (by rw [map_zero]; exact h1.symm) : (1:ℤ) = 0)

end ABC3.Check.PGC
