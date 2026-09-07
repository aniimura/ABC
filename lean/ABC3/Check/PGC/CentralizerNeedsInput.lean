import ABC3.Found.PGC.Theorem42Bijectivity
import Mathlib.LinearAlgebra.Complex.Module

/-!
# [pGC] `CentralizerActsTriviallyOnBase` は体の一般論からは出ない —— 反例つき

`Found/PGC/Theorem42Bijectivity.lean::CentralizerActsTriviallyOnBase` は
pGC Theorem 4.2 の単射性と**同値**な 1 文である（`injective_naturalOuterIso_iff`）。
本ファイルはその文について「p 進の入力なしに証明できるか」を測る。

原文 (pGC 物理 p.7, Theorem 4.2 の証明の括弧書き):

> (Alternatively, the more group-theoretically
> oriented reader may prefer to regard the injectivity of this morphism as a consequence of
> the fact that the centralizer of ΓK in ΓQp is trivial.)

## 測った結果

| 何を測ったか | 結果 |
|---|---|
| 抽象版が現行の定義と同じ形か | `centralizerActsTriviallyOnBase_eq_abstract`（`rfl` で一致） |
| 抽象版は一般の体の 3 つ組で真か | ★**偽**。`(k, K, L) = (ℝ, ℂ, ℂ)` が反例（`not_centralizerActsTriviallyOnBaseAbstract_real`） |
| 純群論版（`C_G(H) ⊆ H`）は一般に真か | ★**偽**。可換群と `H = ⊥` が反例（`not_forall_centralizer_le`） |

★したがって `CentralizerActsTriviallyOnBase K` の証明は
**体の一般論・群の一般論のどちらからも出ない**。必ず p 進局所体に固有の入力
（少なくとも「`K` が代数閉体でない」より強い何か）を通る。
★これは `Found/PGC/CentralizerReduction.lean` が「十分条件」しか与えていないことの
**言い訳ではなく理由**である。

## 反例の中身

`k = ℝ`, `K = ℂ`, `L = ℂ` を取る。`Γ_K = Aut(ℂ/ℂ)` は 1 点なので、
複素共役 `σ = conj` は `Γ_K` を（空虚に）中心化する。`σ` は `ℝ` 上の代数同型で
`K` の上では `ρ = conj ≠ id` として働く。ゆえに結論 `ρ = id` は破れる。

★この反例が効くのは「`K` が代数閉じている」ためであり、p 進局所体では起こらない。
★逆に言えば、`CentralizerActsTriviallyOnBase` の証明は
**`Γ_K` が十分大きいこと**を必ず使わなければならない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Found.PGC

/-- `CentralizerActsTriviallyOnBase` の**抽象版**——基礎体 `k`、その上の体 `K`、
`K` の上の体 `L` の 3 つ組で書いたもの。p 進性も有限次性も代数閉性も要求しない。 -/
def CentralizerActsTriviallyOnBaseAbstract (k K L : Type) [Field k] [Field K] [Field L]
    [Algebra k K] [Algebra K L] : Prop :=
  ∀ (ρ : K ≃ₐ[k] K) (σ : L ≃+* L),
    (∀ x : K, σ (algebraMap K L x) = algebraMap K L (ρ x)) →
    (∀ (g : L ≃ₐ[K] L) (x : L), σ (g x) = g (σ x)) →
    ρ = AlgEquiv.refl

/-- ★抽象版は現行の定義**そのもの**である（`(k, K, L) = (ℚ_p, K, K̄)` を代入すると
定義的に一致する）。★これで下の反例が「同じ文についての反例」であることが担保される。 -/
theorem centralizerActsTriviallyOnBase_eq_abstract {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) :
    CentralizerActsTriviallyOnBase K
      = CentralizerActsTriviallyOnBaseAbstract ℚ_[p] K.carrier K.closure :=
  rfl

/-- `E` の `E`-代数自己同型は恒等だけ（`Theorem42Bijectivity.lean::subsingleton_selfAlgEquiv`
の一般体版）。 -/
theorem subsingleton_selfAlgEquiv_field (E : Type) [Field E] : Subsingleton (E ≃ₐ[E] E) :=
  ⟨fun e f => AlgEquiv.ext fun x => by
    have he : e x = x := by simpa using e.commutes x
    have hf : f x = x := by simpa using f.commutes x
    rw [he, hf]⟩

/-- ★★★**抽象版は偽**。`(k, K, L) = (ℝ, ℂ, ℂ)` が反例。

`Aut(ℂ/ℂ)` が 1 点なので複素共役は空虚に中心化条件を満たすが、`ℂ` の上で恒等ではない。 -/
theorem not_centralizerActsTriviallyOnBaseAbstract_real :
    ¬ CentralizerActsTriviallyOnBaseAbstract ℝ ℂ ℂ := by
  intro h
  have hconj : (Complex.conjAe : ℂ ≃ₐ[ℝ] ℂ) = AlgEquiv.refl := by
    refine h Complex.conjAe Complex.conjAe.toRingEquiv (fun x => by simp) ?_
    intro g x
    haveI := subsingleton_selfAlgEquiv_field ℂ
    have hg : g = AlgEquiv.refl := Subsingleton.elim _ _
    rw [hg]
    rfl
  have hI := congrArg (fun e : ℂ ≃ₐ[ℝ] ℂ => e Complex.I) hconj
  simp only [Complex.conjAe_coe, Complex.conj_I, AlgEquiv.coe_refl, id_eq] at hI
  exact Complex.I_ne_zero (by linear_combination (-1 / 2 : ℂ) * hI)

/-- ★★**純群論版も偽**。「部分群の中心化群はその部分群に入る」は一般には成り立たない
（可換群で `H = ⊥` を取ればよい）。

★原典の言い方（`C_{Γ_Qp}(Γ_K)` が自明）は `(G, H) = (Γ_ℚp, Γ_K)` という**特定の組**の
性質であって、群論の一般論ではない。 -/
theorem not_forall_centralizer_le :
    ¬ ∀ (G : Type) (_ : Group G) (H : Subgroup G), Subgroup.centralizer (H : Set G) ≤ H := by
  intro h
  have hle := h (Multiplicative (ZMod 2)) inferInstance ⊥
  have hmem : (Multiplicative.ofAdd (1 : ZMod 2))
      ∈ Subgroup.centralizer ((⊥ : Subgroup (Multiplicative (ZMod 2))) : Set _) := by
    intro y _
    exact mul_comm _ _
  have hbot := hle hmem
  simp [Subgroup.mem_bot] at hbot

#print axioms centralizerActsTriviallyOnBase_eq_abstract
#print axioms subsingleton_selfAlgEquiv_field
#print axioms not_centralizerActsTriviallyOnBaseAbstract_real
#print axioms not_forall_centralizer_le

end ABC3.Check.PGC
