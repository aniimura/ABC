/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentTowerInduction
import ABC3.Meta.Claim

/-!
# ★★★★★★★★中間体へ降ろす —— `ζ_p` の還元の前半（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

## ★★★★★★★★これは何か —— 「reduce immediately to the case where ζ ∈ K」の前半

原文は elementary claim の証明で

> By separating the extension L[ζ]/K … we reduce immediately to the case of
> wildly ramified L/K such that K contains a primitive p-th root of unity ζ

と書く。★この「還元」を型にすると 2 段になる:

| 段 | 内容 | 状態 |
|---|---|---|
| **前半** | `𝔡_{A,C} ⊆ (𝔡_{A,B})·C`（大きい塔の different は中間の different に含まれる） | ★★**本ファイル** |
| 後半 | `x ∈ (𝔡_{A,B})·C` かつ `x ∈ B` ⟹ `x ∈ 𝔡_{A,B}`（降下） | ★開 |

★★前半は乗法性 `𝔡_{A,C} = 𝔡_{B,C}·(𝔡_{A,B}·C)` から
`Ideal.mul_le_left` だけで出る。

## ★★★これで何が言えるか

`A = O_K`、`B = O_L`、`C = O_{L(ζ)}` と取ると

    `p^n ∈ 𝔡_{K, L(ζ)}` ⟹ `p^n ∈ (𝔡_{K,L})·O_{L(ζ)}`

★したがって「`L(ζ)/K` について claim が言えれば、`L/K` については
**イデアルの拡大を降ろす段だけ**が残る」という形になる。

## ★★★★★★後日の追記（2026-08-28、`§9-904`）—— 後半は**要らなかった**

★降下（`(I·C) ∩ B = I`）を `EC8` として立て、mathlib に無いことまで測った（`§9-903`）。
★★★★**しかし消費側は降下を要求しない**——`Proposition 1.7` が `𝔡_{K,L}` を使うのは
`log-diff` を通してだけであり、`log-diff` は単調（`logDiffOfField_le`）なので

    `log-diff(L) − log-diff(K) ≤ log-diff(M) − log-diff(K) ≤ n·log p`

で `M = L(ζ)` での主張がそのまま `L` での主張を含む。
★詳細は `Found/GenEll/LogDiffDescent.lean`。
★★本ファイルの前半はイデアルの段としては正しく、参照として残す。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★大きい塔の different は中間の different に含まれる -/

/-- ★★★★★★★★**`𝔡_{A,C} ⊆ (𝔡_{A,B})·C`**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★乗法性 `𝔡_{A,C} = 𝔡_{B,C}·(𝔡_{A,B}·C)` から `Ideal.mul_le_left` だけで出る。 -/
theorem differentIdeal_le_map (A B C : Type*) [CommRing A] [CommRing B] [Algebra A B]
    [CommRing C] [Algebra B C] [Algebra A C] [IsScalarTower A B C]
    [IsDomain A] [IsDedekindDomain B] [IsDedekindDomain C]
    [Module.IsTorsionFree A B] [Module.IsTorsionFree B C] [Module.IsTorsionFree A C]
    (hmul : differentIdeal A C
      = differentIdeal B C * Ideal.map (algebraMap B C) (differentIdeal A B)) :
    differentIdeal A C ≤ Ideal.map (algebraMap B C) (differentIdeal A B) := by
  rw [hmul]
  exact Ideal.mul_le_left

/-- ★★**元の形** —— `x ∈ 𝔡_{A,C}` なら `x ∈ (𝔡_{A,B})·C`。 -/
theorem mem_map_of_mem_differentIdeal (A B C : Type*) [CommRing A] [CommRing B] [Algebra A B]
    [CommRing C] [Algebra B C] [Algebra A C] [IsScalarTower A B C]
    [IsDomain A] [IsDedekindDomain B] [IsDedekindDomain C]
    [Module.IsTorsionFree A B] [Module.IsTorsionFree B C] [Module.IsTorsionFree A C]
    (hmul : differentIdeal A C
      = differentIdeal B C * Ideal.map (algebraMap B C) (differentIdeal A B))
    {x : C} (hx : x ∈ differentIdeal A C) :
    x ∈ Ideal.map (algebraMap B C) (differentIdeal A B) :=
  differentIdeal_le_map A B C hmul hx

/-! ## ★出典の紐付け(`.src`) -/

def differentIdeal_le_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——大きい塔の different は中間の different に含まれる)",
    sectionId := "genell-prop-1-7" }

def mem_map_of_mem_differentIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——ζ_p の還元の前半)",
    sectionId := "genell-prop-1-7" }

def differentIdeal_le_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "differentIdeal_eq_differentIdeal_mul_differentIdeal(different の乗法性)"
      (.inMathlib "differentIdeal_eq_differentIdeal_mul_differentIdeal") 1,
    .implicitStep
      ("★原文の「we reduce immediately to the case of wildly ramified L/K such that " ++
       "K contains a primitive p-th root of unity」を型にすると 2 段になる: " ++
       "前半(本ファイル)は 𝔡_{A,C} ⊆ (𝔡_{A,B})·C、" ++
       "後半は降下((I·C) ∩ B = I)である") 2,
    .implicitStep
      ("★★★★降下は B → C が忠実平坦であることから出る" ++
       "——数体の整数環の間では Module.Flat と Ideal.comap_map の組合せになる。" ++
       "mathlib での正確な形はまだ測っていない(2026-08-28)") 4 ]

end ABC3.Found.GenEll
