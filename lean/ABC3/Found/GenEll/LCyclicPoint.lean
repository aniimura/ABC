/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.PadicRedVec
import ABC3.Found.GaloisRep.TateKerLevel
import ABC3.Meta.Claim

/-!
# 第 1205 ブロック —— **`HasLCyclicJ` は `Gal`-安定な位数 `l` の点を与える**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——節点 2 の橋

`HasLCyclicJ`（`T_l E` の `mod l` 還元の中の `Gal`-安定な直線）から
**位数ちょうど `l` の点 `Q` で、その直線 `⟨Q⟩` が `Gal`-安定なもの**を作る。

| 段 | 材料 | 第 |
|---|---|---|
| `v` を `ℤ_l²` へ持ち上げる | `exists_redVec_eq` | 1204 |
| `Q ≔ tateProj₁ (e⁻¹ w)` が `≠ 0` | `tateProj_one_eq_zero_iff` ＋ `redVec_eq_zero_iff` | 1203・1204 |
| 位数がちょうど `l` | `addOrderOf_tateProj_one` | 1201 |
| `σ` の作用が `Q` の直線に入る | `redVec_mulVec_glRed` ＋ 核の記述 | 1204・1203 |

★★★これで `Skeleton/GenEll/LCyclicReading.lean` の**節点 2 の橋**が架かる
——`HasLCyclicJ` から `Lemma 3.5` が要る `addOrderOf Q = l` が出る。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve Matrix ABC3.Meta
open scoped MatrixGroups Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`HasLCyclicJ` は `Gal`-安定な位数 `l` の点を与える**——★**無条件**（第 1205）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`mod l` の安定直線の生成元 `v` を `ℤ_l²` に持ち上げて `e⁻¹` で `T_l E` に戻し、
第 1 層に落とすと**位数ちょうど `l` の点 `Q`** が得られる。
`Gal` の作用は `mod l` で `c` 倍だから、`T_l E` では `c` 倍と `l · T_l E` の差、
すなわち第 1 層では**ちょうど `c` 倍**になる。

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の**節点 2 の橋**である。 -/
theorem exists_stablePoint_of_hasLCyclicJ (E : SSCurve) (l : ℕ) (hl : l.Prime)
    (h : HasLCyclicJ E l) :
    ∃ Q : (E.W.baseChange E.alg).toAffine.Point, addOrderOf Q = l ∧
      ∀ σ : E.alg ≃ₐ[E.fld] E.alg, ∃ k : ℕ, galPoint E.W σ Q = k • Q := by
  haveI hfact : Fact l.Prime := ⟨hl⟩
  haveI : NeZero l := ⟨hl.ne_zero⟩
  obtain ⟨e, v, hv0, hstab⟩ := h hfact
  obtain ⟨w, hw⟩ := exists_redVec_eq l v
  set t : E.tate l := e.symm w with ht
  have het : e t = w := e.apply_symm_apply w
  set Q : (E.W.baseChange E.alg).toAffine.Point :=
    ((tateProj (E.W.baseChange E.alg) l 1 t : (E.W.baseChange E.alg).toAffine.Point)) with hQ
  -- ☆`Q ≠ 0`——`mod l` で `v ≠ 0` だから
  have hQne : Q ≠ 0 := by
    intro h0
    obtain ⟨g, hg⟩ := (tateProj_one_eq_zero_iff (E.W.baseChange E.alg) l t).1 h0
    have hwl : w = l • e g := by rw [← het, ← hg, map_nsmul]
    have hz : redVec l w = 0 := (redVec_eq_zero_iff l w).2 ⟨e g, hwl⟩
    rw [hw] at hz
    exact hv0 hz
  refine ⟨Q, addOrderOf_tateProj_one E.W hl t hQne, ?_⟩
  intro σ
  -- ☆`σ` の行列の `mod l` 還元は `v` を `c` 倍する
  have hmem : (glRedPadic l (galRep E.W l e σ))
      ∈ ((galRep E.W l e).range).map (glRedPadic l) :=
    Subgroup.mem_map_of_mem _ (MonoidHom.mem_range.2 ⟨σ, rfl⟩)
  obtain ⟨c, hc⟩ := hstab _ hmem
  set k : ℕ := c.val with hk
  have hcv : ((k : ℕ) : ZMod l) = c := by simp [hk]
  -- ☆`ℤ_l` の側で `galMat σ ·ᵥ w − k · w` は `l` で割れる
  have hgal : e (galTate E.W l σ t) = (galMat E.W l e σ).mulVec w := by
    rw [galMat_apply, het]
  have hmat : galMat E.W l e σ
      = ((galRep E.W l e σ : GL (Fin 2) ℤ_[l]) : Matrix (Fin 2) (Fin 2) ℤ_[l]) :=
    (galRep_coe E.W l e σ).symm
  have hred : redVec l ((galMat E.W l e σ).mulVec w) = c • v := by
    rw [hmat, redVec_mulVec_glRed, hw, hc]
  have hkv : (k : ℕ) • v = c • v := by
    rw [← Nat.cast_smul_eq_nsmul (ZMod l) k v, hcv]
  have hdiff : redVec l ((galMat E.W l e σ).mulVec w - k • w) = 0 := by
    rw [redVec_sub, hred, redVec_nsmul, hw, hkv, sub_self]
  obtain ⟨u, hu⟩ := (redVec_eq_zero_iff l _).1 hdiff
  -- ☆`T_l E` の側へ戻す
  have hsub : galTate E.W l σ t - k • t = l • e.symm u := by
    apply e.injective
    rw [map_sub, map_nsmul, map_nsmul, hgal, e.apply_symm_apply u, het]
    exact hu
  have hzero : ((tateProj (E.W.baseChange E.alg) l 1
      (galTate E.W l σ t - k • t) : (E.W.baseChange E.alg).toAffine.Point)) = 0 := by
    rw [hsub]
    have hs := tateProj_smul_pow_eq_zero (E.W.baseChange E.alg) l 1 (e.symm u)
    simpa using hs
  refine ⟨k, ?_⟩
  rw [map_sub, map_nsmul] at hzero
  have hcoe : ((tateProj (E.W.baseChange E.alg) l 1 (galTate E.W l σ t) : _) -
      k • (tateProj (E.W.baseChange E.alg) l 1 t : _) : (E.W.baseChange E.alg).toAffine.Point)
      = 0 := by simpa using hzero
  have hfin := sub_eq_zero.mp hcoe
  rw [← tateProj_galTate E.W l 1 σ t]
  exact hfin

/-! ## ★出典の紐付け(`.src`) -/

def exists_stablePoint_of_hasLCyclicJ.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(HasLCyclicJ は Gal-安定な位数 l の点を与える。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_stablePoint_of_hasLCyclicJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateProj_one_eq_zero_iff(第 1203、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateProj_one_eq_zero_iff") 1,
    .citation "[ABC3]" "redVec_eq_zero_iff(第 1204、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.redVec_eq_zero_iff") 1,
    .citation "[ABC3]" "addOrderOf_tateProj_one(第 1201、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.addOrderOf_tateProj_one") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1205）**——`Skeleton/GenEll/LCyclicReading.lean` の" ++
       "**節点 2 の橋**である。☆`mod l` の安定直線の生成元 `v` を `ℤ_l²` に持ち上げ、" ++
       "`e⁻¹` で `T_l E` に戻して第 1 層に落とすと位数ちょうど `l` の点 `Q` になる。" ++
       "★`Gal` の作用は `mod l` で `c` 倍だから、`T_l E` では `c` 倍と `l · T_l E` の差、" ++
       "すなわち第 1 層では**ちょうど `c` 倍**である。") 3 ]

end ABC3.Found.GenEll
