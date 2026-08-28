/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArithCartierComap
import ABC3.Found.GenEll.HeightAdditive
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★`ht_{D^n} = n·ht_D`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★これは何か —— 一般の `X` へ渡る道の 2 本目

`§9-879`（段 C2f）で `ht_{ψ^*E}(x) = ht_E(ψ ∘ x)` が入った。
★一般の `X`（`L_ℚ` が豊富）への還元にはもう 1 本、

    `ht_{L^n} = n·ht_L`

が要る——原文が「[some positive tensor power of] the ample line bundle」と書く
**その冪の分**である。★★本ファイルはそれを取る。

## ★★★機構 —— 既存の加法性の帰納

★`htArith_tensor_unconditional`（`§9-HeightAdditive`）を `n` 回繰り返すだけである。
★★各段で「引き戻しイデアルが `0` でない」ことが要るが、
`pullbackMul_all`（引き戻しは積を保つ、**無条件**）と
`𝓞_F` が整域であることから帰納で出る（`pullbackIdeal_npow_ne_zero`）。

## ★★★★これで一般の `X` への還元は 1 本になる

| 本 | 内容 | 状態 |
|---|---|---|
| (1) | `ht_{ψ^*E}(x) = ht_E(ψ ∘ x)` | ★済（`§9-879`） |
| (2) | `ht_{L^n} = n·ht_L` | ★★**本ファイル** |
| (3) | `L^n` が**非常に豊富**（`ψ^*O(1) ≅ L^n`） | 段 E3 |
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★算術因子の冪 -/

/-- ★**算術因子の `n` 乗** —— 因子は積、Green 関数は和の `n` 倍。 -/
def ArithCartier.npow {X : Scheme.{0}} (D : ArithCartier X) : ℕ → ArithCartier X
  | 0 => ArithCartier.one X
  | n + 1 => (D.npow n).tensor D

@[simp] theorem ArithCartier.npow_zero {X : Scheme.{0}} (D : ArithCartier X) :
    D.npow 0 = ArithCartier.one X := rfl

@[simp] theorem ArithCartier.npow_succ {X : Scheme.{0}} (D : ArithCartier X) (n : ℕ) :
    D.npow (n + 1) = (D.npow n).tensor D := rfl

/-! ## ★★引き戻しイデアルは `0` にならない -/

/-- ★★**冪の引き戻しイデアルも `0` でない**。

★`pullbackMul_all`（引き戻しは積を保つ、**無条件**）と `𝓞_F` が整域であることによる。 -/
theorem pullbackIdeal_npow_ne_zero (F : Type) [Field F] [NumberField F] {X : Scheme.{0}}
    (D : ArithCartier X) (xF : specRingOfIntegers F ⟶ X)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (n : ℕ) :
    pullbackIdeal F (D.npow n).divisor xF ≠ 0 := by
  induction n with
  | zero =>
      show pullbackIdeal F (ArithCartier.one X).divisor xF ≠ 0
      rw [show (ArithCartier.one X).divisor = (⊤ : X.IdealSheafData) from rfl,
        pullbackIdeal_top]
      exact top_ne_bot
  | succ n ih =>
      show pullbackIdeal F ((D.npow n).tensor D).divisor xF ≠ 0
      rw [show ((D.npow n).tensor D).divisor = (D.npow n).divisor * D.divisor from rfl,
        pullbackMul_all F xF]
      exact mul_ne_zero ih hD

/-! ## ★★★★★★★★★★★`ht_{D^n} = n·ht_D` -/

/-- ★★★★★★★★★★★**`ht_{D^n} = n·ht_D`**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★原文が「[some positive tensor power of] the ample line bundle」と書く**その冪の分**である。
★★機構は `htArith_tensor_unconditional` の帰納だけ。 -/
theorem htArith_npow (F : Type) [Field F] [NumberField F] {X : Scheme.{0}}
    (D : ArithCartier X) (xF : specRingOfIntegers F ⟶ X)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (n : ℕ) :
    htArith F (D.npow n) xF = (n : ℝ) * htArith F D xF := by
  induction n with
  | zero =>
      show htArith F (ArithCartier.one X) xF = ((0 : ℕ) : ℝ) * htArith F D xF
      rw [htArith_one, Nat.cast_zero, zero_mul]
  | succ n ih =>
      show htArith F ((D.npow n).tensor D) xF = _
      rw [htArith_tensor_unconditional F (D.npow n) D xF
        (pullbackIdeal_npow_ne_zero F D xF hD n) hD, ih]
      push_cast
      ring

/-- ★★**高さで有界な点の集合は冪でスケールする**。 -/
theorem setOf_htArith_npow (F : Type) [Field F] [NumberField F] {X : Scheme.{0}}
    (D : ArithCartier X) (n : ℕ) (C : ℝ)
    (hD : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0) :
    {xF : specRingOfIntegers F ⟶ X | htArith F (D.npow n) xF ≤ C}
      = {xF | (n : ℝ) * htArith F D xF ≤ C} := by
  ext xF
  rw [Set.mem_setOf_eq, Set.mem_setOf_eq, htArith_npow F D xF (hD xF) n]

/-! ## ★出典の紐付け(`.src`) -/

def ArithCartier.npow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(算術因子の n 乗)",
    sectionId := "genell-prop-1-4" }

def pullbackIdeal_npow_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(冪の引き戻しイデアルも 0 でない)",
    sectionId := "genell-prop-1-4" }

def htArith_npow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(ht_{D^n} = n·ht_D)",
    sectionId := "genell-prop-1-4" }

def htArith_npow.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_tensor_unconditional(高さの加法性、無条件)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_tensor_unconditional") 2,
    .citation "[ABC3]" "pullbackMul_all(引き戻しは積を保つ、無条件)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackMul_all") 2,
    .implicitStep
      ("★原文が「[some positive tensor power of] the ample line bundle」と書く" ++
       "**その冪の分**である") 2,
    .implicitStep
      ("★★これで一般の X への還元は 1 本になった: " ++
       "(1) ht_{ψ^*E}(x) = ht_E(ψ ∘ x)——§9-879 済、" ++
       "(2) ht_{L^n} = n·ht_L——本ファイル、" ++
       "(3) L^n が**非常に豊富**(ψ^*O(1) ≅ L^n)——段 E3 が残る") 4 ]

end ABC3.Found.GenEll
