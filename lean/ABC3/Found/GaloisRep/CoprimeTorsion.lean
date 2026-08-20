import ABC3.Found.GaloisRep.ThreeTorsion

/-!
# Galois (G1) 第 37 ブロック —— **★★★★★互いに素な `n` への還元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★これは**群論だけ**で出る

    gcd(m, n) = 1,  E[m] 有限,  E[n] 有限  ⟹  E[mn] 有限

★機構: `P ↦ (m•P, n•P)` は `E[mn] → E[n] × E[m]` の準同型で、
**核が自明**(`m•P = n•P = 0` なら位数が `gcd(m,n) = 1` を割る)。

★★**分点多項式も乗法公式も要らない**——`addOrderOf` の性質だけである。

## ★★これで `E[6]` が出る

第 33(`E[2]`)と第 36(`E[3]`)を合わせるだけ。

## ★★★一般の `n` への含意

★これで **`n` を素冪に分解できる**——残るのは `E[p^k]` の有限性である。
★★`E[p]` から `E[p^k]` へは `E[p^k]/E[p^{k-1}]` の議論が要り、
そこはやはり乗法公式(または同種写像の次数)に落ちる。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- 互いに素な `m, n` について `E[mn]` は `E[n] × E[m]` に単射で入る。 -/
theorem finite_torsion_mul_of_coprime {m n : ℕ} (hmn : Nat.Coprime m n)
    (hm : {P : W.toAffine.Point | m • P = 0}.Finite)
    (hn : {P : W.toAffine.Point | n • P = 0}.Finite) :
    {P : W.toAffine.Point | (m * n) • P = 0}.Finite := by
  classical
  haveI := hm.to_subtype
  haveI := hn.to_subtype
  refine Set.Finite.of_finite_image
    (f := fun P : W.toAffine.Point => (m • P, n • P)) ?_ ?_
  · refine Set.Finite.subset (Set.Finite.prod hn hm) ?_
    rintro _ ⟨P, hP, rfl⟩
    refine ⟨?_, ?_⟩
    · show n • (m • P) = 0
      rw [smul_smul, Nat.mul_comm]; exact hP
    · show m • (n • P) = 0
      rw [smul_smul]; exact hP
  · rintro P hP Q hQ hPQ
    simp only [Prod.mk.injEq] at hPQ
    have hd : m • (P - Q) = 0 := by rw [smul_sub, hPQ.1, sub_self]
    have hd' : n • (P - Q) = 0 := by rw [smul_sub, hPQ.2, sub_self]
    have h1 : addOrderOf (P - Q) ∣ Nat.gcd m n :=
      Nat.dvd_gcd (addOrderOf_dvd_of_nsmul_eq_zero hd) (addOrderOf_dvd_of_nsmul_eq_zero hd')
    rw [hmn] at h1
    have h2 : (1 : ℕ) • (P - Q) = 0 := addOrderOf_dvd_iff_nsmul_eq_zero.1 h1
    have : P - Q = 0 := by simpa using h2
    exact sub_eq_zero.1 this


/-- ★★★★**`E[6]` は有限**——第 33 と第 36 を合わせるだけ。 -/
theorem finite_six_torsion (h2 : (2 : F) ≠ 0) (h3 : (3 : F) ≠ 0) :
    {P : W.toAffine.Point | (6 : ℕ) • P = 0}.Finite := by
  have h : (6 : ℕ) = 2 * 3 := by norm_num
  rw [h]
  exact finite_torsion_mul_of_coprime W (by norm_num)
    (finite_two_torsion W h2) (finite_three_torsion W h2 h3)

/-! ## ★出典の紐付け(`.src`) -/

def finite_torsion_mul_of_coprime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(互いに素な n への還元)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
