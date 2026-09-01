/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicEval
import ABC3.Meta.Claim

/-!
# 第 1106 ブロック —— **係数環つきの adic 評価**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か

`AdicEval.lean` の `evalAdic` は **`PowerSeries ℤ` 専用**である（第 1105 で実測）。
`Lemma 3.5` の分母払い（第 1092-1105）は万有な環 `PowerSeries (ℤ[ζ_l])` を経由するので、
**係数環 `S` と環準同型 `φ : S →+* R` を受ける版**が要る。

☆中身は `evalAdic` と同じ——部分和が `I` 進 Cauchy であることを使い、
`IsPrecomplete` で極限を取る。★係数を `φ (coeff n g)` に差し替えるだけである。

## ★★★★これで何が繋がるか

第 1104 で測ったとおり mathlib は
`instance : IsAdicComplete (.span {X}) (PowerSeries R)` を持つ。したがって

    `A₀ ≔ PowerSeries (ℤ[ζ_l])`  →（`φ : ℤ[ζ_l] →+* R`、`ζ_l ↦ ζ`）→  `R`

という特殊化が本ファイルの `evalAdicMap` で書ける。
★`A₁ ≔ PowerSeries (ℚ(ζ_l))` では `l` が可逆なので `1 − ζ^i` も可逆であり、
既存の `hu` つきの補題がそのまま使える。☆`A₀ → A₁` の単射性で降ろせば `hu` が消える。
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

variable {R S : Type} [CommRing R] [CommRing S] {I : Ideal R}

/-- ★係数環つきの部分和 `∑_{n<N} φ(c_n)·q^n`。 -/
noncomputable def partialEvalMap (φ : S →+* R) (g : PowerSeries S) (q : R) (N : ℕ) : R :=
  ∑ n ∈ Finset.range N, φ (PowerSeries.coeff n g) * q ^ n

/-- ★★**係数環つきの部分和も `I` 進 Cauchy である**。 -/
theorem partialEvalMap_cauchy (φ : S →+* R) (g : PowerSeries S) {q : R} (hq : q ∈ I)
    {m n : ℕ} (hmn : m ≤ n) :
    partialEvalMap φ g q m ≡ partialEvalMap φ g q n [SMOD (I ^ m • ⊤ : Submodule R R)] := by
  rw [SModEq.sub_mem]
  have hsub : partialEvalMap φ g q m - partialEvalMap φ g q n
      = -∑ k ∈ Finset.Ico m n, φ (PowerSeries.coeff k g) * q ^ k := by
    rw [partialEvalMap, partialEvalMap, Finset.range_eq_Ico, Finset.range_eq_Ico,
      ← Finset.sum_Ico_consecutive (fun k => φ (PowerSeries.coeff k g) * q ^ k)
        (Nat.zero_le m) hmn]
    ring
  rw [hsub]
  have hmem : ∀ k ∈ Finset.Ico m n, φ (PowerSeries.coeff k g) * q ^ k ∈ I ^ m := by
    intro k hk
    have hkm : m ≤ k := (Finset.mem_Ico.1 hk).1
    have hpow : q ^ k = q ^ m * q ^ (k - m) := by
      rw [← pow_add]
      congr 1
      omega
    rw [hpow, ← mul_assoc]
    exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq m))
  have hs : (∑ k ∈ Finset.Ico m n, φ (PowerSeries.coeff k g) * q ^ k) ∈ I ^ m :=
    Submodule.sum_mem _ hmem
  simpa using neg_mem hs

/-- ★★★**係数環つきの adic 評価**——部分和の `I` 進極限。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def evalAdicMap [IsPrecomplete I R] (φ : S →+* R) (g : PowerSeries S)
    (q : R) (hq : q ∈ I) : R :=
  Classical.choose (IsPrecomplete.prec' (partialEvalMap φ g q) (partialEvalMap_cauchy φ g hq))

theorem evalAdicMap_spec [IsPrecomplete I R] (φ : S →+* R) (g : PowerSeries S)
    (q : R) (hq : q ∈ I) (n : ℕ) :
    partialEvalMap φ g q n ≡ evalAdicMap φ g q hq [SMOD (I ^ n • ⊤ : Submodule R R)] :=
  Classical.choose_spec
    (IsPrecomplete.prec' (partialEvalMap φ g q) (partialEvalMap_cauchy φ g hq)) n

/-- ★★★★**一意性**——`IsHausdorff` の下で極限は 1 つ。 -/
theorem evalAdicMap_unique [IsAdicComplete I R] (φ : S →+* R) (g : PowerSeries S)
    (q : R) (hq : q ∈ I) (L : R)
    (hL : ∀ n, partialEvalMap φ g q n ≡ L [SMOD (I ^ n • ⊤ : Submodule R R)]) :
    evalAdicMap φ g q hq = L :=
  (IsHausdorff.eq_iff_smodEq (I := I)).2
    (fun n => (evalAdicMap_spec φ g q hq n).symm.trans (hL n))

/-- ★★★★★**`S = ℤ` なら元の `evalAdic` に一致する**。 -/
theorem evalAdicMap_int [IsAdicComplete I R] (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) :
    evalAdicMap (Int.castRingHom R) f q hq = evalAdic f q hq := by
  refine evalAdicMap_unique (Int.castRingHom R) f q hq _ (fun n => ?_)
  have h : partialEvalMap (Int.castRingHom R) f q n = partialEval f q n := by
    rw [partialEvalMap, partialEval]
    exact Finset.sum_congr rfl (fun k _ => by simp)
  rw [h]
  exact evalAdic_spec f q hq n

/-! ## ★出典の紐付け(`.src`) -/

def evalAdicMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線——係数環つきの完備環への特殊化)",
    sectionId := "genell-def-3-3" }

def evalAdicMap_int.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(係数環つきの評価は S = ℤ で元の evalAdic に一致する)",
    sectionId := "genell-def-3-3" }

def evalAdicMap.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★**2026-09-01（第 1106）**——これが `Lemma 3.5` の分母払い（第 1092-1105）の" ++
       "**特殊化の段**である。☆`A₀ ≔ PowerSeries (ℤ[ζ_l]) → R` を `evalAdicMap` で書き、" ++
       "`A₁ ≔ PowerSeries (ℚ(ζ_l))` で `hu` つきの補題を使ってから" ++
       "`A₀ → A₁` の単射性で降ろす。★これで `d + 1 < l` が外れる。") 15 ]

end ABC3.Found.GaloisRep
