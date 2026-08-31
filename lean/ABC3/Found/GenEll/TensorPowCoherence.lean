/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SecPow
import ABC3.Meta.Claim

/-!
# ★★★★★★テンソル冪のコヒーレンス `M^{⊗(n·k)} ≅ (M^{⊗n})^{⊗k}`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★これは段 E2 の**最後の道具**である

`§9-817` の有限被覆は切断が別々の `M^{⊗n_j}` に属する。
`§9-818` が `X_{s^{⊗k}} = X_s`（非消失軌跡は冪で変わらない）を与えたが、
★`secPow (M^{⊗n}) s k` の行き先は **`(M^{⊗n})^{⊗k}`** であって `M^{⊗(n·k)}` ではない。

★★射 `X ⟶ ℙᴺ_R` は**ひとつの**束の切断を要求するので、この 2 つを繋ぐ必要がある。
それが本ファイルである。

## ★★★★★機構 —— モノイダル構造だけ

| 補題 | 中身 |
|---|---|
| `presheafTensorPowAdd` | `M^{⊗(a+b)} ≅ M^{⊗a} ⊗ M^{⊗b}`（`a` についての帰納、結合子と左単位子） |
| `presheafTensorPowMul` | `M^{⊗(n·k)} ≅ (M^{⊗n})^{⊗k}`（`k` についての帰納、`Add` 版を使う） |

★`presheafTensorPow M (n+1) = M ⊗ presheafTensorPow M n` は `rfl` なので、
指数の書き換えは `eqToIso` で吸収する。

★★`Mul` の側では `n·(k+1) = n + n·k` と**左に `n` を出す**のが要点である
——そうすれば `M^{⊗n} ⊗ M^{⊗(n·k)}` になり、
`whiskerLeftIso` で `(M^{⊗n})^{⊗(k+1)}` に届く。★★★**組み紐（braiding）は要らない**。

## ★★★これで段 E2 の道具は揃った

    §9-817  有限被覆
    §9-818  X_{s^{⊗k}} = X_s
    ★本ファイル  M^{⊗(n·k)} ≅ (M^{⊗n})^{⊗k}
    既存      nonVanishing_of_iso（同型で X_s は移る、段 D1）

★残るのは **`lcm` の帳簿**——`L = lcm(n_j)`（または `∏ n_j`）を取り、
各 `j` で `k_j = L / n_j` として `secPow` と本ファイルの同型を合成する段である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-- ★★★★★**`M^{⊗(a+b)} ≅ M^{⊗a} ⊗ M^{⊗b}`**。

★`a` についての帰納。`a = 0` は左単位子、`a+1` は結合子である。 -/
noncomputable def presheafTensorPowAdd (M : X.PresheafOfModules) :
    ∀ a b : ℕ, presheafTensorPow M (a + b) ≅ presheafTensorPow M a ⊗ presheafTensorPow M b
  | 0, b =>
      (eqToIso (by rw [Nat.zero_add])).trans (λ_ (presheafTensorPow M b)).symm
  | a + 1, b =>
      (eqToIso (show presheafTensorPow M (a + 1 + b) = M ⊗ presheafTensorPow M (a + b) by
          rw [show a + 1 + b = (a + b) + 1 from by omega]; rfl)).trans
        ((whiskerLeftIso M (presheafTensorPowAdd M a b)).trans
          ((α_ M (presheafTensorPow M a) (presheafTensorPow M b)).symm))

/-- ★★★★★★**`M^{⊗(n·k)} ≅ (M^{⊗n})^{⊗k}`** —— 次数揃えのコヒーレンス。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`n·(k+1) = n + n·k` と**左に `n` を出す**のが要点である
——`whiskerLeftIso` だけで届き、組み紐は要らない。 -/
noncomputable def presheafTensorPowMul (M : X.PresheafOfModules) (n : ℕ) :
    ∀ k : ℕ, presheafTensorPow M (n * k) ≅ presheafTensorPow (presheafTensorPow M n) k
  | 0 => eqToIso (by rw [Nat.mul_zero]; rfl)
  | k + 1 =>
      (eqToIso (show presheafTensorPow M (n * (k + 1)) = presheafTensorPow M (n + n * k) by
          rw [show n * (k + 1) = n + n * k from by ring])).trans
        ((presheafTensorPowAdd M n (n * k)).trans
          (whiskerLeftIso (presheafTensorPow M n) (presheafTensorPowMul M n k)))

/-! ## ★出典の紐付け(`.src`) -/

def presheafTensorPowAdd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(M^{⊗(a+b)} ≅ M^{⊗a} ⊗ M^{⊗b})",
    sectionId := "genell-prop-1-4" }

def presheafTensorPowMul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(M^{⊗(n·k)} ≅ (M^{⊗n})^{⊗k}——次数揃えのコヒーレンス)",
    sectionId := "genell-prop-1-4" }

def presheafTensorPowMul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "PresheafOfModules.monoidalCategory(結合子・単位子)"
      (.inMathlib "PresheafOfModules.monoidalCategory") 6,
    .citation "[ABC3]" "nonVanishing_secPow(X_{s^{⊗k}} = X_s、§9-818)"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_secPow") 6,
    .citation "[ABC3]" "nonVanishing_of_iso(同型で X_s は移る、段 D1)"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_of_iso") 6,
    .implicitStep
      ("★残るのは **lcm の帳簿**——L = lcm(n_j)(または ∏ n_j)を取り、" ++
       "各 j で k_j = L / n_j として secPow と本ファイルの同型を合成する段である。" ++
       "★★これが済めば段 E2 が閉じ、残るのは E1d(貼り合わせ)と E3(immersion 性)になる") 6 ]

end ABC3.Found.GenEll
