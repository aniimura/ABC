import ABC3.Found.GaloisRep.PrimeGenerators
import Mathlib.RingTheory.Localization.AtPrime.Basic

/-!
# Galois (G5) 第 136 ブロック —— **★★★★★★局所化の極大イデアルは単項**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★2 つの場合を局所化へ移す

第 135 で `P = (x − c, y − y₀)` が出た。★本ブロックは、これを `P` での局所化へ送り、
**極大イデアルが 1 元で生成される**ことを示す。

| 場合 | 一意化元 | 使う分解 |
|---|---|---|
| `Ψ₂Sq(c) = 0`(分岐、`s = 0`) | `z` | 第 132 `z² = (x−c)·e(x)`、`e(c) ≠ 0` |
| `Ψ₂Sq(c) ≠ 0`(不分岐、`s ≠ 0`) | `x − c` | 第 133 `(z−s)(z+s) = (x−c)·g(x)` |

★★機構はどちらも同じである——**分解の片側が局所化で単元になる**:

* 分岐: `e(x) ≡ e(c) ≠ 0 (mod P)` だから `e(x) ∉ P`。よって `x − c = z²·e(x)⁻¹ ∈ (z)`。
* 不分岐: `z + s ∉ P`(入れば `2s ∈ P` だが `2s` は単元)。よって `z − s = (x−c)·g(x)·(z+s)⁻¹ ∈ (x−c)`。

★★★もう一方の生成元 `y − y₀` は `2(y − y₀) = (z − s) − a₁(x − c)` から従う(`2` が単元)。

## ★★★★抽象的な局所化で述べる——実測で効いた

当初 `Localization.AtPrime P` の**具体型**で述べたところ、
**文の詰め込みだけで 200000 heartbeats を超えた**(2026-08-20 実測、8.3 秒でタイムアウト)。
★`(S : Type) [IsLocalization P.primeCompl S]` と抽象化したら **0.24 秒**になった。
★★具体型のインスタンス探索が重い場合は、抽象化が単なる一般化以上の効果を持つ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mem_of_isUnit_mul` | ★単元倍はイデアル所属を保つ |
| `isUnit_loc_of_not_mem` | ★`P` の外の元は局所化で単元 |
| `eval₂_genX_eq_algebraMap` | ★`eval₂(x)` は `F[X]` からの `algebraMap` |
| `algebraMap_poly_not_mem` | ★★`q(c) ≠ 0` なら `q(x) ∉ P` |
| `map_pair_eq_span_single` | ★★2 元生成の像が単項になる十分条件 |
| `map_eq_span_genZ` | ★★★★★**分岐——極大イデアル `= (z)`** |
| `map_eq_span_genX` | ★★★★★**不分岐——極大イデアル `= (x − c)`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

/-! ## ★汎用の小道具 -/

/-- ★単元倍はイデアル所属を保つ。 -/
theorem mem_of_isUnit_mul {R : Type} [CommRing R] {I : Ideal R} {a x : R} (ha : IsUnit a)
    (h : a * x ∈ I) : x ∈ I := by
  have := Ideal.mul_mem_left I (↑ha.unit⁻¹) h
  rwa [← mul_assoc, IsUnit.val_inv_mul, one_mul] at this

/-- ★`P` の外の元は `P` での局所化で単元になる。 -/
theorem isUnit_loc_of_not_mem {R : Type} [CommRing R] (P : Ideal R) [P.IsPrime]
    (S : Type) [CommRing S] [Algebra R S] [IsLocalization P.primeCompl S]
    {x : R} (hx : x ∉ P) : IsUnit (algebraMap R S x) :=
  IsLocalization.map_units _ (⟨x, hx⟩ : P.primeCompl)

/-- ★2 元生成イデアルの像が単項になるための十分条件。 -/
theorem map_pair_eq_span_single {R S : Type} [CommRing R] [CommRing S] [Algebra R S]
    {u v π : R}
    (hu : algebraMap R S u ∈ Ideal.span {algebraMap R S π})
    (hv : algebraMap R S v ∈ Ideal.span {algebraMap R S π})
    (hπ : algebraMap R S π ∈ Ideal.map (algebraMap R S) (Ideal.span ({u, v} : Set R))) :
    Ideal.map (algebraMap R S) (Ideal.span ({u, v} : Set R))
      = Ideal.span {algebraMap R S π} := by
  refine le_antisymm ?_ (by rw [Ideal.span_le]; rintro z (rfl : z = _); exact hπ)
  rw [Ideal.map_span, Ideal.span_le]
  rintro z ⟨w, (rfl | rfl), rfl⟩
  · exact hu
  · exact hv

variable {F : Type} [Field F]

/-- ★`eval₂` による `x` への代入は `F[X]` からの `algebraMap` と一致する。 -/
theorem eval₂_genX_eq_algebraMap (W : WeierstrassCurve.Affine F) (p : Polynomial F) :
    Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W) p
      = algebraMap (Polynomial F) W.CoordinateRing p := by
  rw [eval₂_genX, algebraMap_polynomial_coordinateRing]

/-- ★★`q(c) ≠ 0` なら `q(x) ∉ P`——第 135 の「`mod P` で定数」から。 -/
theorem algebraMap_poly_not_mem (W : WeierstrassCurve.Affine F) (P : Ideal W.CoordinateRing)
    [hPp : P.IsPrime] {c : F} (hx : genX W - algebraMap F W.CoordinateRing c ∈ P)
    {q : Polynomial F} (hq : q.eval c ≠ 0) :
    algebraMap (Polynomial F) W.CoordinateRing q ∉ P := by
  intro hmem
  have h := quotient_algebraMap_poly W P hx q
  rw [Ideal.Quotient.eq_zero_iff_mem.2 hmem, eq_comm, Ideal.Quotient.eq_zero_iff_mem] at h
  exact hPp.ne_top (Ideal.eq_top_of_isUnit_mem P h
    ((isUnit_iff_ne_zero.2 hq).map (algebraMap F W.CoordinateRing)))

variable [DecidableEq F]

set_option maxHeartbeats 1000000 in
/-- ★★★★★**分岐(2 等分点)の場合——局所化の極大イデアルは `z` で生成される**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 132 の `z² = (x − c)·e(x)`(`e(c) ≠ 0`)で `e(x)` が局所化の単元になるので
`x − c = z²·e(x)⁻¹ ∈ (z)`。★★`y − y₀` は `2(y − y₀) = z − a₁(x − c)` から。 -/
theorem map_eq_span_genZ [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    (P : Ideal W.CoordinateRing) [hPp : P.IsPrime]
    (S : Type) [CommRing S] [Algebra W.CoordinateRing S]
    [IsLocalization P.primeCompl S] {c : F}
    (hc0 : W.Ψ₂Sq.eval c = 0)
    (hx : genX W - algebraMap F W.CoordinateRing c ∈ P)
    (hz : genZ W ∈ P)
    (hspan : P = Ideal.span ({genX W - algebraMap F W.CoordinateRing c,
        genY W - algebraMap F W.CoordinateRing ((0 - W.a₁ * c - W.a₃) / 2)} : Set _)) :
    Ideal.map (algebraMap W.CoordinateRing S) P
      = Ideal.span {algebraMap W.CoordinateRing S (genZ W)} := by
  set A := algebraMap F W.CoordinateRing with hA
  set φ := algebraMap W.CoordinateRing S with hφ
  set y₀ : F := (0 - W.a₁ * c - W.a₃) / 2 with hy₀
  obtain ⟨e, hec, hze⟩ := genZ_sq_factor W h2 hc0
  rw [eval₂_genX_eq_algebraMap] at hze
  have hEu : IsUnit (φ (algebraMap (Polynomial F) W.CoordinateRing e)) :=
    isUnit_loc_of_not_mem P S (algebraMap_poly_not_mem W P hx hec)
  have hxmem : φ (genX W - A c) ∈ Ideal.span {φ (genZ W)} := by
    rw [Ideal.mem_span_singleton']
    refine ⟨φ (genZ W) * ↑hEu.unit⁻¹, ?_⟩
    have h1 : φ (genZ W) * ↑hEu.unit⁻¹ * φ (genZ W)
        = (φ (genZ W) * φ (genZ W)) * ↑hEu.unit⁻¹ := by ring
    rw [h1, ← map_mul, ← sq, hze, map_mul, mul_assoc, IsUnit.mul_val_inv, mul_one]
  have h2' : (2 : F) ≠ 0 := h2.ne_zero
  have hy2 : (2 : F) * y₀ = 0 - W.a₁ * c - W.a₃ := by rw [hy₀]; field_simp
  have hidCR : (2 : W.CoordinateRing) * (genY W - A y₀) = genZ W - A W.a₁ * (genX W - A c) := by
    have hkey : (2 : W.CoordinateRing) * A y₀ = A 0 - A W.a₁ * A c - A W.a₃ := by
      rw [show (2 : W.CoordinateRing) = A 2 from (map_ofNat _ _).symm, ← map_mul, hy2,
        map_sub, map_sub, map_mul]
    rw [mul_sub, hkey, genZ, map_zero]; ring
  have h2CR : IsUnit (2 : W.CoordinateRing) := by
    rw [show (2 : W.CoordinateRing) = A 2 from (map_ofNat _ _).symm]; exact h2.map A
  have h2S : IsUnit (2 : S) := by
    rw [show (2 : S) = φ 2 from (map_ofNat _ _).symm]; exact h2CR.map φ
  have hymem : φ (genY W - A y₀) ∈ Ideal.span {φ (genZ W)} := by
    refine mem_of_isUnit_mul h2S ?_
    have hrw : (2 : S) * φ (genY W - A y₀)
        = φ ((2 : W.CoordinateRing) * (genY W - A y₀)) := by rw [map_mul, map_ofNat]
    have hsplit : φ (genZ W - A W.a₁ * (genX W - A c))
        = φ (genZ W) - φ (A W.a₁) * φ (genX W - A c) := by rw [map_sub, map_mul]
    rw [hrw, hidCR, hsplit]
    exact Ideal.sub_mem _ (Ideal.mem_span_singleton_self _)
      (Ideal.mul_mem_left _ _ hxmem)
  rw [hspan]
  exact map_pair_eq_span_single hxmem hymem (Ideal.mem_map_of_mem _ (hspan ▸ hz))

set_option maxHeartbeats 1000000 in
/-- ★★★★★**不分岐の場合——局所化の極大イデアルは `x − c` で生成される**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 133 の `(z − s)(z + s) = (x − c)·g(x)` で `z + s` が局所化の単元になるので
`z − s = (x − c)·g(x)·(z + s)⁻¹ ∈ (x − c)`。 -/
theorem map_eq_span_genX (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F)
    (P : Ideal W.CoordinateRing) [hPp : P.IsPrime]
    (S : Type) [CommRing S] [Algebra W.CoordinateRing S]
    [IsLocalization P.primeCompl S] {c s : F}
    (hs : s ^ 2 = W.Ψ₂Sq.eval c) (hsne : s ≠ 0)
    (hx : genX W - algebraMap F W.CoordinateRing c ∈ P)
    (hz : genZ W - algebraMap F W.CoordinateRing s ∈ P)
    (hspan : P = Ideal.span ({genX W - algebraMap F W.CoordinateRing c,
        genY W - algebraMap F W.CoordinateRing ((s - W.a₁ * c - W.a₃) / 2)} : Set _)) :
    Ideal.map (algebraMap W.CoordinateRing S) P
      = Ideal.span {algebraMap W.CoordinateRing S
          (genX W - algebraMap F W.CoordinateRing c)} := by
  set A := algebraMap F W.CoordinateRing with hA
  set φ := algebraMap W.CoordinateRing S with hφ
  set y₀ : F := (s - W.a₁ * c - W.a₃) / 2 with hy₀
  obtain ⟨g, hg⟩ := genZ_sub_factor W hs
  rw [eval₂_genX_eq_algebraMap] at hg
  have hznot : genZ W + A s ∉ P := by
    intro hmem
    have hd : (genZ W + A s) - (genZ W - A s) = A (2 * s) := by
      rw [map_mul, show A 2 = (2 : W.CoordinateRing) from map_ofNat _ _]; ring
    have hmem2 : A (2 * s) ∈ P := by rw [← hd]; exact P.sub_mem hmem hz
    exact hPp.ne_top (Ideal.eq_top_of_isUnit_mem P hmem2
      ((isUnit_iff_ne_zero.2 (mul_ne_zero h2.ne_zero hsne)).map A))
  have hZu : IsUnit (φ (genZ W + A s)) := isUnit_loc_of_not_mem P S hznot
  have hzmem : φ (genZ W - A s) ∈ Ideal.span {φ (genX W - A c)} := by
    rw [Ideal.mem_span_singleton']
    refine ⟨φ (algebraMap (Polynomial F) W.CoordinateRing g) * ↑hZu.unit⁻¹, ?_⟩
    have h1 : φ (algebraMap (Polynomial F) W.CoordinateRing g) * ↑hZu.unit⁻¹
          * φ (genX W - A c)
        = (φ (genX W - A c) * φ (algebraMap (Polynomial F) W.CoordinateRing g))
          * ↑hZu.unit⁻¹ := by ring
    rw [h1, ← map_mul, ← hg, map_mul, mul_assoc, IsUnit.mul_val_inv, mul_one]
  have h2' : (2 : F) ≠ 0 := h2.ne_zero
  have hy2 : (2 : F) * y₀ = s - W.a₁ * c - W.a₃ := by rw [hy₀]; field_simp
  have hidCR : (2 : W.CoordinateRing) * (genY W - A y₀)
      = (genZ W - A s) - A W.a₁ * (genX W - A c) := by
    have hkey : (2 : W.CoordinateRing) * A y₀ = A s - A W.a₁ * A c - A W.a₃ := by
      rw [show (2 : W.CoordinateRing) = A 2 from (map_ofNat _ _).symm, ← map_mul, hy2,
        map_sub, map_sub, map_mul]
    rw [mul_sub, hkey, genZ]; ring
  have h2CR : IsUnit (2 : W.CoordinateRing) := by
    rw [show (2 : W.CoordinateRing) = A 2 from (map_ofNat _ _).symm]; exact h2.map A
  have h2S : IsUnit (2 : S) := by
    rw [show (2 : S) = φ 2 from (map_ofNat _ _).symm]; exact h2CR.map φ
  have hymem : φ (genY W - A y₀) ∈ Ideal.span {φ (genX W - A c)} := by
    refine mem_of_isUnit_mul h2S ?_
    have hrw : (2 : S) * φ (genY W - A y₀)
        = φ ((2 : W.CoordinateRing) * (genY W - A y₀)) := by rw [map_mul, map_ofNat]
    have hsplit : φ ((genZ W - A s) - A W.a₁ * (genX W - A c))
        = φ (genZ W - A s) - φ (A W.a₁) * φ (genX W - A c) := by rw [map_sub, map_mul]
    rw [hrw, hidCR, hsplit]
    exact Ideal.sub_mem _ hzmem (Ideal.mul_mem_left _ _ (Ideal.mem_span_singleton_self _))
  rw [hspan]
  exact map_pair_eq_span_single (Ideal.mem_span_singleton_self _) hymem
    (Ideal.mem_map_of_mem _ (hspan ▸ hx))

/-! ## ★出典の紐付け(`.src`) -/

def map_eq_span_genZ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——局所化の極大イデアルが単項であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
