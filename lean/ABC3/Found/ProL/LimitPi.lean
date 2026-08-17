/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Algebra.Group.Subgroup.Basic
import Mathlib.Algebra.Group.Pi.Lemmas

/-!
# 極限と積の交換 —— チェーン `prol` の葉 `limit-pi-comm`

★`ResearchPaper/frdi-decomposition.json` の `prol` チェーンの葉。
最終目標は [FrdI] `Definition 2.8, (ii)`。

## ★★設計を変えた記録(2026-08-18)

★台帳では `CategoryTheory.Limits.limitFlipCompLimIsoLimitCompLim`
(`CategoryTheory/Limits/Fubini.lean:545`)を使うと書いた。**実装で変えた。**

理由: 圏論版を `ProfiniteGrp` に当てるには
「`OpenNormalSubgroup M × Discrete Nat.Primes ⥤ ProfiniteGrp`」の双関手を組む必要があり、
そこから**使える形**(整合族の言葉)へ戻す手間が、直接書くより大きい。
★下流(`M[l]` の構成)が要るのは**整合族の入れ替え**そのものなので、
それを直接書いた。★台帳の `★限界` が「実装すると変わりうる」と書いていた通りである。

## ★中身

`J` を前順序の添字、`G : J → ι → Type` を逆系の族、`T` を遷移射とすると

  **積の逆系の整合族 ≃* 各成分の整合族の族**

これが「極限と積は交換する」の中身であり、証明は**引数の入れ替えだけ**である
(`left_inv` / `right_inv` はいずれも `rfl`)。

★★**向きに注意**: `OpenNormalSubgroup M` を添字にすると、`U ≤ V` に対して
遷移射は `M/U → M/V` と**大きい方へ**行く。ここではその向きで書いてある。
-/

namespace ABC3.Found.ProL

universe uJ ui uG

variable {J : Type uJ} {ι : Type ui} [Preorder J]

/-- ★逆系の**整合族**がなす部分群。 -/
def compatSubgroup {G : J → Type uG} [∀ j, Group (G j)]
    (T : ∀ {j k : J}, j ≤ k → G j →* G k) : Subgroup (∀ j, G j) where
  carrier := {x | ∀ (j k : J) (h : j ≤ k), T h (x j) = x k}
  one_mem' := by intro j k h; simp
  mul_mem' := by
    intro a b ha hb j k h
    show T h (a j * b j) = a k * b k
    rw [map_mul, ha j k h, hb j k h]
  inv_mem' := by
    intro a ha j k h
    show T h (a j)⁻¹ = (a k)⁻¹
    rw [map_inv, ha j k h]

@[simp] theorem mem_compatSubgroup {G : J → Type uG} [∀ j, Group (G j)]
    (T : ∀ {j k : J}, j ≤ k → G j →* G k) (x : ∀ j, G j) :
    x ∈ compatSubgroup T ↔ ∀ (j k : J) (h : j ≤ k), T h (x j) = x k := Iff.rfl

/-- 積の逆系の遷移射(成分ごとに `T` を当てる)。 -/
def piTrans {G : J → ι → Type uG} [∀ j i, Group (G j i)]
    (T : ∀ {j k : J}, j ≤ k → ∀ i, G j i →* G k i) {j k : J} (h : j ≤ k) :
    (∀ i, G j i) →* (∀ i, G k i) where
  toFun x := fun i => T h i (x i)
  map_one' := funext fun _ => map_one _
  map_mul' _ _ := funext fun _ => map_mul _ _ _

@[simp] theorem piTrans_apply {G : J → ι → Type uG} [∀ j i, Group (G j i)]
    (T : ∀ {j k : J}, j ≤ k → ∀ i, G j i →* G k i) {j k : J} (h : j ≤ k)
    (x : ∀ i, G j i) (i : ι) : piTrans T h x i = T h i (x i) := rfl

/-- ★★★**極限と積の交換** —— 積の逆系の整合族 ≃* 各成分の整合族の族。

★証明は引数の入れ替えだけである(`left_inv` / `right_inv` は `rfl`)。
**それが「極限と積が交換する」ことの中身**である。 -/
def limitPiMulEquiv {G : J → ι → Type uG} [∀ j i, Group (G j i)]
    (T : ∀ {j k : J}, j ≤ k → ∀ i, G j i →* G k i) :
    compatSubgroup (G := fun j => ∀ i, G j i) (fun h => piTrans T h)
      ≃* (∀ i : ι, compatSubgroup (G := fun j => G j i) (fun h => T h i)) where
  toFun x := fun i => ⟨fun j => x.1 j i, fun j k h => congrFun (x.2 j k h) i⟩
  invFun y := ⟨fun j i => (y i).1 j, fun j k h => funext fun i => (y i).2 j k h⟩
  left_inv _ := rfl
  right_inv _ := rfl
  map_mul' _ _ := rfl

@[simp] theorem limitPiMulEquiv_apply_coe {G : J → ι → Type uG} [∀ j i, Group (G j i)]
    (T : ∀ {j k : J}, j ≤ k → ∀ i, G j i →* G k i)
    (x : compatSubgroup (G := fun j => ∀ i, G j i) (fun h => piTrans T h)) (i : ι) (j : J) :
    ((limitPiMulEquiv T x) i : ∀ j, G j i) j = (x : ∀ j, ∀ i, G j i) j i := rfl

end ABC3.Found.ProL
