/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38KummerExists
import ABC3.Found.GenEll.Thm38ZetaPi
import ABC3.Meta.Claim

/-!
# 第 1212 ブロック —— **Kummer 拡大の上で `α` の座標の式が出る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——`TateSetup` の底変換は要らなかった

第 1211 で `K = AdjoinRoot (Xˡ − C q)` の上に
`π`（`πˡ = q`）と `σ`（`σ ζ = ζ`・`σ π = ζ·π`）を取った。

☆`α` の座標の式は `sigma_coord_alpha`（第 1174 の下地）が与えるが、
**それは任意の可換群で成り立つ**——`TateSetup` を要求していない。
★要るのは「`ζ` が原始 `l` 乗根」「`πˡ = Q`」「`Q` が無限位数」だけである。

★★★したがって **`TateSetup` を `K` へ底変換する段は要らない**——
基礎体 `K₀` の側の仮説（`l ∤ v(q)`、`q` が無限位数、`ζ₀` が原始 `l` 乗根）だけで
`K` の上の座標の式が出る。

## ★★★★到達点

    σ (ζᵃ πᵇ) = ζ^{a'} π^{b'} Qⁿ  ⟹  l ∣ (a + b − a')  かつ  l ∣ (b − b')

☆すなわち `σ` は `mod l` の座標で `α = (1 1 / 0 1)` として作用する。
-/

namespace ABC3.Found.GenEll

open Polynomial ABC3.Meta

/-! ## ★★★★★★★★★★★★★★★★Kummer 拡大の上の `α` -/

open AdjoinRoot in
set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**Kummer 拡大の上で `σ` は `α` として作用する**——★**無条件**（第 1212）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆仮説は基礎体 `K₀` の側だけ:
`l ∤ v(q)`・`q` が無限位数・`ζ₀` が原始 `l` 乗根。

★★★`TateSetup` を `K` へ底変換する段は**要らなかった**
——`sigma_coord_alpha` は任意の可換群で成り立つからである。 -/
theorem kummer_sigma_coord_alpha {K₀ : Type} [Field K₀] {l : ℕ} (hl : Nat.Prime l)
    (v : K₀ˣ →* Multiplicative ℤ) (q : K₀ˣ)
    (hnd : ¬ ((l : ℤ) ∣ Multiplicative.toAdd (v q)))
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0)
    (ζ₀ : K₀ˣ) (hζ₀ : IsPrimitiveRoot ((ζ₀ : K₀)) l) :
    ∃ (π : (AdjoinRoot (X ^ l - C (q : K₀)))ˣ)
      (σ : (AdjoinRoot (X ^ l - C (q : K₀)))ˣ →* (AdjoinRoot (X ^ l - C (q : K₀)))ˣ),
      π ^ l = Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) q ∧
      ∀ a b a' b' : ℤ,
        (∃ n : ℤ,
          σ ((Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) ζ₀) ^ a
              * π ^ b)
            = (Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) ζ₀) ^ a'
              * π ^ b'
              * (Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) q) ^ n) →
        ((l : ℤ) ∣ (a + b - a')) ∧ ((l : ℤ) ∣ (b - b')) := by
  haveI : Fact (Irreducible (X ^ l - C (q : K₀))) :=
    fact_irreducible_X_pow_sub_C_of_not_pow hl (not_lth_power_of_val hl v q hnd)
  obtain ⟨π, σ, hπl, hσζ, hσπ⟩ := exists_units_sigma_kummer hl v q hnd ζ₀ hζ₀
  refine ⟨π, σ, hπl, ?_⟩
  have hinj : Function.Injective
      (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀)))) := (algebraMap K₀ _).injective
  have hminj : Function.Injective
      (Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀)))
        : K₀ →* AdjoinRoot (X ^ l - C (q : K₀)))) := by
    intro x y hxy
    refine Units.ext (hinj ?_)
    exact congrArg (fun u : (AdjoinRoot (X ^ l - C (q : K₀)))ˣ =>
      (u : AdjoinRoot (X ^ l - C (q : K₀)))) hxy
  intro a b a' b' h
  refine sigma_coord_alpha hl.pos ?_ ?_ hπl ?_ σ hσζ hσπ a b a' b' h
  · rw [← map_pow, units_pow_eq_one_of_isPrimitiveRoot hζ₀, map_one]
  · intro n hn hnl hcon
    refine units_pow_ne_one_of_isPrimitiveRoot hζ₀ n hn hnl ?_
    refine hminj ?_
    rw [map_pow, hcon, map_one]
  · intro j hj
    refine hqinf j ?_
    refine hminj ?_
    rw [map_zpow, hj, map_one]

/-! ## ★出典の紐付け(`.src`) -/

def kummer_sigma_coord_alpha.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Kummer 拡大の上で σ は α として作用する。★無条件)",
    sectionId := "genell-thm-3-8" }

def kummer_sigma_coord_alpha.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_units_sigma_kummer(第 1211、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_units_sigma_kummer") 1,
    .citation "[ABC3]" "sigma_coord_alpha(第 1174 の下地、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sigma_coord_alpha") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1212）**——**`TateSetup` を `K` へ底変換する段は" ++
       "要らなかった**。☆`sigma_coord_alpha` は**任意の可換群**で成り立ち、" ++
       "要るのは「`ζ` が原始 `l` 乗根」「`πˡ = Q`」「`Q` が無限位数」だけだからである。" ++
       "★仮説はすべて基礎体 `K₀` の側にある。") 3 ]

end ABC3.Found.GenEll
