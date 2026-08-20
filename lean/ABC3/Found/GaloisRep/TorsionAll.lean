import ABC3.Found.GaloisRep.MulOrder

/-!
# Galois (G1) 第 54 ブロック —— **★★★★★★★捩れの有限性(一般の `n`)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★`E[n]` は有限である

    (n : F) ≠ 0  ⟹  { P | n•P = 0 } は有限

★これで `Interface/GaloisRep/Torsion.lean` の `torsion_finite` が埋まる。

## ★★使う多項式

第 53 ブロックより、`n•P = 0` なら **`n` の約数** `k` で `ΨSq_k(x_P) = 0` となる。
★したがって

    p := ∏_{k ∣ n} ΨSq_k

の根に `x_P` が入る。★★`k ∣ n` と `(n : F) ≠ 0` から `(k : F) ≠ 0` なので、
mathlib の `ΨSq_ne_zero` により **各因子が `0` でない**——積も `0` でない ✅

## ★★★標数の扱い

★`(n : F) ≠ 0` だけを仮定する。`k ∣ n` に限ったことが効いている——
`k ≤ n` のままだと標数 `p ≤ n` のとき `ΨSq_p = 0` の可能性を排除できない。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `finite_torsion_of_root` | ★★一般の多項式版の有限性 |
| `finite_torsion` | ★★★★★★★**`E[n]` は有限** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★一般の多項式版——`p ≠ 0` で `m • P = 0 → p(x_P) = 0` なら `E[m]` は有限。 -/
theorem finite_torsion_of_root {m : ℕ} (p : F[X]) (hp : p ≠ 0)
    (hform : ∀ (x y : F) (h : W.toAffine.Nonsingular x y),
      m • (Point.some x y h) = 0 → p.IsRoot x) :
    {P : W.toAffine.Point | m • P = 0}.Finite := by
  classical
  have hroots : {x : F | p.IsRoot x}.Finite := Polynomial.finite_setOf_isRoot hp
  have hpairs : {xy : F × F | p.IsRoot xy.1 ∧ W.toAffine.Equation xy.1 xy.2}.Finite := by
    refine Set.Finite.subset (Set.Finite.biUnion hroots
      (fun x _ => Set.Finite.image (fun y => (x, y)) (finite_y_of_x W x))) ?_
    rintro ⟨x, y⟩ ⟨hx, hy⟩
    exact Set.mem_biUnion hx ⟨y, hy, rfl⟩
  refine Set.Finite.of_finite_image (f := xyOf W) ?_ ?_
  · refine Set.Finite.subset (Set.Finite.insert none
      (Set.Finite.image (fun xy => (some xy : Option (F × F))) hpairs)) ?_
    rintro _ ⟨P, hP, rfl⟩
    rcases P with _ | ⟨x, y, h⟩
    · exact Set.mem_insert _ _
    · exact Set.mem_insert_of_mem _ ⟨(x, y), ⟨hform x y h hP, h.left⟩, rfl⟩
  · rintro (_ | ⟨x, y, h⟩) hP (_ | ⟨x', y', h'⟩) hP' hxx
    · rfl
    · exact absurd hxx (by simp [xyOf])
    · exact absurd hxx (by simp [xyOf])
    · simp only [xyOf, Option.some.injEq, Prod.mk.injEq] at hxx
      obtain ⟨rfl, rfl⟩ := hxx
      rfl

/-- ★★★★★★★**捩れは有限**——`(n : F) ≠ 0` なら `E[n]` は有限集合。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが `Interface/GaloisRep/Torsion.lean` の `TorsionStructureData.torsion_finite` である。 -/
theorem finite_torsion (n : ℕ) (hn : 1 ≤ n) (hchar : (n : F) ≠ 0) :
    {P : W.toAffine.Point | n • P = 0}.Finite := by
  classical
  have hdvd : ∀ k : ℕ, k ∣ n → (k : F) ≠ 0 := by
    intro k hk hk0
    obtain ⟨t, rfl⟩ := hk
    apply hchar
    push_cast
    rw [hk0, zero_mul]
  refine finite_torsion_of_root W (∏ k ∈ n.divisors, W.ΨSq (k : ℤ)) ?_ ?_
  · refine Finset.prod_ne_zero_iff.2 (fun k hk => ?_)
    refine W.ΨSq_ne_zero ?_
    have hkn : k ∣ n := (Nat.mem_divisors.mp hk).1
    push_cast
    exact hdvd k hkn
  · intro x y h hP
    obtain ⟨k, hk1, hkn, hkroot⟩ := exists_divisor_root W h n hn hP
    have hmem : k ∈ n.divisors := Nat.mem_divisors.mpr ⟨hkn, by omega⟩
    simp only [Polynomial.IsRoot, eval_prod]
    exact Finset.prod_eq_zero hmem hkroot

/-! ## ★出典の紐付け(`.src`) -/

def finite_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(n-捩れ部分群の有限性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
