import ABC3.Found.PGC.SubgroupCorrespondenceConstruction
import ABC3.Found.PGC.AbsClosureModules
import ABC3.Found.PGC.AxTowerDecay
import ABC3.Found.PGC.WildDepthDescent

/-!
# [pGC] `AxWildDescent K (fun _ => p)` —— wild 深さ 1 段の降下を体の層で閉じる

`Found/PGC/WildDepthDescent.lean` が抽象 Galois の段で作った
`exists_natDegree_minpoly_descent_div` を `PAdicLocalField` に代入し、
★★**`AxWildDescent K (fun _ => (p : ℝ))` を無条件に証明する**。

## 位置づけ(何がどれだけ良くなったか)

| 出典 | 定数 `c k` | `∏_{k ∈ Icc 1 n} c k` |
|---|---|---|
| `AxTowerDecay.axWildDescent_pow`(前の波) | `p^k` | `p^{n(n+1)/2}` |
| ★本ファイル `axWildDescent_prime` | ★**`p`**(`k` に依らない) | `p^n` |
| `AxEpsilonDecay.axDecay`(目標) | `p^{(1/(p−1))p^{1−k}}` | `p^{p/(p−1)²}`(有界) |

★★**まだ `AxLemma` は出ない**(`p^n` は非有界)。
★残っているのは「1 段の損失を `p` から `p^{(1/(p−1))p^{1−k}}` に絞る」ことだけで、
それには**分岐**(異なるイデアル・跳び)が要る ——
`Found/PGC/CyclicJumpNorm.lean` の `k = 1` の段と同じ話である。
★★**逆に、深さの降り方・平均化・指数の勘定はもう詰まっていない。**

## ★★持ち場の目標が偽だったこと

本ファイルの上流 `WildDepthDescent.lean` の冒頭 docstring を読むこと。要点だけ:

* 配られた目標「`wildDepth K x = k ≥ 2` なら深さ `k−1` の**中間体** `M` が在る」は
  ★**偽**である(`not_forall_exists_relIndex_padicValNat_eq`、反例 `A₄ ⊃ C₃`)。
* 正しい降下は `x'` を `K(x)` の**外**に出す。本ファイルの `x'` は
  `x` の Galois 閉包 `M` の中の元であって、`K(x)` の中には無い。
* ★★`k ≥ 2` と `k = 1` を分ける必要は無かった。1 本の定理が両方を扱う。

## 配管(★測った。前の波の「中間体 2 層で越えられない」は**ここでは当たらない**)

```
haveI := isGalois_closure K            -- Found/PGC/SubgroupCorrespondenceConstruction.lean:50
FiniteGaloisIntermediateField.adjoin K.carrier ({x} : Set K.closure)   -- lean-idioms.md #153
  → FiniteDimensional / Normal / Algebra.IsSeparable が全部 inferInstance(最後だけ
    IntermediateField.isSeparable_tower_bot を 1 行)
  → NormedField ↥M / IsUltrametricDist ↥M / DistribMulAction (↥M ≃ₐ[K.carrier] ↥M) ↥M
    も全部 inferInstance
‖(a : ↥M)‖ = ‖(a : K.closure)‖                          --  rfl
((a - b : ↥M) : K.closure) = (a : K.closure) - (b : K.closure)  --  rfl
AlgEquiv.restrictNormalHom_surjective / AlgEquiv.restrictNormal_commutes  -- 作用の往復
IntermediateField.minpoly_eq                            -- minpoly の往復
```

★★**`IntermediateField` を 1 つしか作っていない**(`K(x)` を中間体として作らない)。
`wildDepth` は `MulAction.stabilizer` の指数として測るので、
`lean-idioms.md` #59 の「中間体 2 層の `rfl` が kernel を止める」に**一度も触らない**。
★これが前の波と違う唯一の点である。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 原典 (Ax 1970) は定数の勘定を地の文で畳む。本ファイルの `c k ≡ p` という中間結果は
   原典に対応する文が無い。`.src` は `AxTowerDecay.lean` / `AxEpsilonDecay.lean` と同じ項目
   (pGC 物理 p.6 Corollary 3.1)を指す。
2. 持ち場が名指しした「中間体で 1 段下がる」は偽なので捨てた(上記)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**wild 深さ 1 段の降下が体の層で閉じた**(損失 `‖(p : K^al)‖⁻¹`)。

`x` の Galois 閉包 `M` の中で `WildDepthDescent.exists_natDegree_minpoly_descent_div`
を使う。★`x'` は `M` の元であって `K(x)` の元ではない ——
「中間体で降りる」形が偽であること(`not_forall_exists_relIndex_padicValNat_eq`)と両立する。 -/
theorem axWildDescent_normInv (K : PAdicLocalField p) :
    AxWildDescent K (fun _ => ‖((p : ℕ) : K.closure)‖⁻¹) := by
  intro ε hε x hxε hdvd
  haveI := isGalois_closure K
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  set M : IntermediateField K.carrier K.closure :=
    (FiniteGaloisIntermediateField.adjoin K.carrier ({x} : Set K.closure)).toIntermediateField
      with hMdef
  have hxM : x ∈ M :=
    FiniteGaloisIntermediateField.subset_adjoin K.carrier ({x} : Set K.closure) rfl
  haveI : FiniteDimensional K.carrier M := inferInstance
  haveI : Normal K.carrier M := inferInstance
  haveI : Algebra.IsSeparable K.carrier M := IntermediateField.isSeparable_tower_bot K.carrier M
  haveI : CharZero (M : Type _) :=
    charZero_of_injective_algebraMap (algebraMap K.carrier (M : Type _)).injective
  set xM : (M : Type _) := ⟨x, hxM⟩ with hxMdef
  have hmin : minpoly K.carrier xM = minpoly K.carrier x := IntermediateField.minpoly_eq xM
  have hdeg0 : (minpoly K.carrier x).natDegree ≠ 0 :=
    (minpoly.natDegree_pos (Algebra.IsIntegral.isIntegral x)).ne'
  have hk1 : 1 ≤ wildDepth K x := by
    rw [wildDepth]
    exact (padicValNat_dvd_iff_le hdeg0).mp (by simpa using hdvd)
  have hk : padicValNat p (minpoly K.carrier xM).natDegree = (wildDepth K x - 1) + 1 := by
    rw [hmin, ← wildDepth]
    omega
  have hxMbound : ∀ g : (M : Type _) ≃ₐ[K.carrier] (M : Type _), ‖g • xM - xM‖ ≤ ε := by
    intro g
    obtain ⟨σ, hσ⟩ :=
      AlgEquiv.restrictNormalHom_surjective (F := K.carrier) (K₁ := (M : Type _)) K.closure g
    have hcoe : ((g • xM : (M : Type _)) : K.closure) = σ • x := by
      rw [← hσ]
      exact AlgEquiv.restrictNormal_commutes σ (M : Type _) xM
    show ‖((g • xM - xM : (M : Type _)) : K.closure)‖ ≤ ε
    rw [show ((g • xM - xM : (M : Type _)) : K.closure)
        = ((g • xM : (M : Type _)) : K.closure) - x from rfl, hcoe]
    exact hxε σ
  have hpne : ((p : ℕ) : (M : Type _)) ≠ 0 := Nat.cast_ne_zero.mpr (Fact.out : p.Prime).pos.ne'
  obtain ⟨y, hy1, hy2, hy3⟩ :=
    exists_natDegree_minpoly_descent_div (p := p) hε hpne hxMbound (wildDepth K x - 1) hk
  have hminy : minpoly K.carrier y = minpoly K.carrier (y : K.closure) :=
    IntermediateField.minpoly_eq y
  have hnormp : ‖((p : ℕ) : (M : Type _))‖ = ‖((p : ℕ) : K.closure)‖ := rfl
  refine ⟨(y : K.closure), ?_, ?_, ?_⟩
  · rw [wildDepth, ← hminy]
    omega
  · show ‖x - (y : K.closure)‖ ≤ ‖((p : ℕ) : K.closure)‖⁻¹ * ε
    rw [← hnormp]
    exact hy2
  · intro σ
    have hcoe : σ • (y : K.closure)
        = ((AlgEquiv.restrictNormalHom (M : Type _) σ • y : (M : Type _)) : K.closure) :=
      (AlgEquiv.restrictNormal_commutes σ (M : Type _) y).symm
    show ‖σ • (y : K.closure) - (y : K.closure)‖ ≤ ‖((p : ℕ) : K.closure)‖⁻¹ * ε
    rw [hcoe, ← hnormp]
    exact hy3 _

def axWildDescent_normInv.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**`AxWildDescent K (fun _ => p)`** —— `‖(p : K^al)‖ = 1/p` を入れただけ。

★`AxTowerDecay.axWildDescent_pow`(`c k = p^k`)より真に良い。
★★それでも `∏_{k ∈ Icc 1 n} p = p^n` は非有界なので `AxLemma` は出ない。
残っているのは**分岐による絞り込み**だけである(モジュール docstring の表)。 -/
theorem axWildDescent_prime (K : PAdicLocalField p) :
    AxWildDescent K (fun _ => (p : ℝ)) := by
  have h := axWildDescent_normInv K
  rw [show (fun _ : ℕ => (p : ℝ)) = (fun _ : ℕ => ‖((p : ℕ) : K.closure)‖⁻¹) by
    funext k
    rw [norm_natCast_p_closure, inv_inv]]
  exact h

def axWildDescent_prime.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end ABC3.Found.PGC
