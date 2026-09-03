/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.FullImageWitness
import ABC3.Found.GenEll.PadicRedVec
import ABC3.Meta.Claim

/-!
# 第 1295 ブロック —— **`ζ_l` が基礎体にあれば `det` は `mod l` で 1**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1292 の入力（`det = 1`）

`det ρ(σ)` は円分指標である（`det_cyclotomic_full`、在庫）。
★`σ` が原始 `l` 乗根 `ζ` を固定するなら、`ζ = ζ^k`（`k = (det mod l).val`）から
`k = 1`、すなわち **`det ρ(σ) ≡ 1 (mod l)`** である。

☆`ZMod (l^1)` と `ZMod l` の型の食い違いは、**イデアルの言葉**を経由して避ける
（`ker (toZModPow 1) = span {l^1} = span {l} = maximalIdeal = ker toZMod`）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

variable {K L : Type} [Field K] [DecidableEq K] [CharZero K]
  [Field L] [DecidableEq L] [Algebra K L] [IsAlgClosed L]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★
**`σ` が原始 `l` 乗根を固定すれば `det ρ(σ) ≡ 1 (mod l)`**——★**無条件**（第 1295）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`det_cyclotomic_full`（在庫）に `n = 1` を当てるだけである。 -/
theorem toZMod_det_eq_one_of_fixed_root
    (W : WeierstrassCurve K) [WeierstrassCurve.IsElliptic (W.baseChange L).toAffine]
    (l : ℕ) [Fact l.Prime] (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (ζ : L) (hζ : IsPrimitiveRoot ζ l) (hfix : σ ζ = ζ) :
    PadicInt.toZMod ((galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det) = 1 := by
  have hl : l.Prime := Fact.out
  haveI : NeZero (l ^ 1) := ⟨pow_ne_zero _ hl.ne_zero⟩
  have hζl : ζ ^ (l ^ 1) = 1 := by
    rw [pow_one]
    exact hζ.pow_eq_one
  have h := det_cyclotomic_full W l e σ 1 ζ hζl
  rw [hfix] at h
  -- h : ζ = ζ ^ (toZModPow 1 D).val
  set D := (galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det with hD
  set k := (PadicInt.toZModPow 1 D).val with hk
  have h1lt : 1 < l := hl.one_lt
  have hpow : l ^ 1 = l := pow_one l
  have hklt : k < l := by
    have h0 := ZMod.val_lt (PadicInt.toZModPow 1 D)
    rw [← hk] at h0
    omega
  have heq : ζ ^ (1 : ℕ) = ζ ^ k := by
    rw [pow_one]
    exact h
  have hk1 : (1 : ℕ) = k := hζ.pow_inj h1lt hklt heq
  haveI : Fact (1 < l ^ 1) := ⟨by rw [hpow]; exact h1lt⟩
  have hval : (PadicInt.toZModPow 1 D) = 1 := by
    refine ZMod.val_injective _ ?_
    rw [ZMod.val_one, ← hk]
    exact hk1.symm
  -- イデアルの言葉で `toZMod` へ移す
  have hker : D - 1 ∈ (Ideal.span {(l : ℤ_[l]) ^ 1} : Ideal ℤ_[l]) := by
    rw [← PadicInt.ker_toZModPow]
    show PadicInt.toZModPow 1 (D - 1) = 0
    rw [map_sub, hval, map_one, sub_self]
  have hmax : D - 1 ∈ (Ideal.span {(l : ℤ_[l])} : Ideal ℤ_[l]) := by
    rwa [pow_one] at hker
  have hz : PadicInt.toZMod (D - 1) = 0 := by
    rw [← RingHom.mem_ker, PadicInt.ker_toZMod, PadicInt.maximalIdeal_eq_span_p]
    exact hmax
  rw [map_sub, map_one, sub_eq_zero] at hz
  exact hz

/-- ★★★★★★★★★★★★★★★★
**`mod l` に落とした行列の `det` は 1**——★**無条件**（第 1295）。

☆`glRedPadic` は成分ごとの `toZMod` なので、`det` も `toZMod` で移る。 -/
theorem det_glRed_eq_one_of_fixed_root
    (W : WeierstrassCurve K) [WeierstrassCurve.IsElliptic (W.baseChange L).toAffine]
    (l : ℕ) [Fact l.Prime] (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (ζ : L) (hζ : IsPrimitiveRoot ζ l) (hfix : σ ζ = ζ) :
    (((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
      Matrix (Fin 2) (Fin 2) (ZMod l))).det = 1 := by
  rw [coe_glRedPadic]
  rw [show ((galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).map PadicInt.toZMod)
      = PadicInt.toZMod.mapMatrix (galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]) from rfl]
  rw [← RingHom.map_det]
  exact toZMod_det_eq_one_of_fixed_root W l e σ ζ hζ hfix

/-! ## ★出典の紐付け(`.src`) -/

def toZMod_det_eq_one_of_fixed_root.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が原始 l 乗根を固定すれば det ρ(σ) ≡ 1 (mod l)。★無条件)",
    sectionId := "genell-thm-3-8" }

def det_glRed_eq_one_of_fixed_root.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(mod l に落とした行列の det は 1。★無条件)",
    sectionId := "genell-thm-3-8" }

def toZMod_det_eq_one_of_fixed_root.needs : List ProofObligation :=
  [ .citation "[ABC3]" "det_cyclotomic_full(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.det_cyclotomic_full") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1295）**——第 1292 の入力（`det = 1`）である。" ++
       "☆`ZMod (l^1)` と `ZMod l` の型の食い違いは**イデアルの言葉**で避けた" ++
       "（`ker (toZModPow 1) = span {l^1} = span {l} = maximalIdeal = ker toZMod`）。") 2 ]

end ABC3.Found.GenEll
