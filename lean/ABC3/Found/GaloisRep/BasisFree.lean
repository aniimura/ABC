import ABC3.Found.GaloisRep.TateSurj
import Mathlib.RingTheory.RootsOfUnity.PrimitiveRoots

/-!
# Galois (G5) 第 210 ブロック —— **★★★★★★★★★★`hgen` は非退化性だけから出る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★基底性が完全に閉じた

第 209 で `P, Q` が `E[l^n]` を**生成する**ことが取れた。★本ブロックは残りの半分——
**`P` の位数がちょうど `l^n` であること**——を、`ℤ_l` 線型性を一切使わずに
**数え上げだけ**で出す。

### ★★★★★★★数え上げで自由性が出る

`e : T_l E ≃+ ℤ_l²` は**加法同型でしかない**(`ℤ_l` 加群の構造を入れていない)。
★したがって「`single₀` は `l·ℤ_l²` に入らない」という線型代数は使えない。
★★しかし次の数え上げで足りる:

| 段 | 内容 |
|---|---|
| 1 | `ord P = l^a`・`ord Q = l^b` と置く(`a, b ≤ n`) |
| 2 | 生成性から `Fin l^a × Fin l^b → E[l^n]` が**全射** |
| 3 | `#E[l^n] = l^{2n}` なので `l^{2n} ≤ l^{a+b}` |
| 4 | `a, b ≤ n` と合わせて **`a = b = n`** |

★★★`a < n` を仮定すると `l^{n-1}·P = 0` となり、上限が `l^{2n-1} < l^{2n}` に落ちて矛盾する。

### ★★★★★★★★非退化性から原始根性へ

`ζ := e_{l^n}(P, Q)` とする。★`ζ^m = 1` なら、任意の `R = aP + bQ` に対し

    e(m·P, R) = e(P, P)^{am} · e(P, Q)^{bm} = 1 · (ζ^m)^b = 1

なので**非退化性**から `m·P = 0`、すなわち `l^n ∣ m`(位数がちょうど `l^n` だから)。
★★これは `IsPrimitiveRoot ζ (l^n)` そのものであり、
`IsPrimitiveRoot.eq_pow_of_pow_eq_one` から第 207 の `hgen` が出る。

## ★★★★★★★★★★残件 (ii) は消えた

第 207 の `.needs` に書いた**残件 (ii)「Lean 上の配線」は本ブロックで完了**した。
★`det_galRep_eq_cyclotomic` に残っているのは**非退化性ただ 1 つ**であり、
それは `normEDS` が楕円列であること(mathlib 自身の TODO、残件 (i))に帰着している。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `card_torsion_le_of_two_gen` | ★★★★★2 元で生成される捩れ群の位数の上限 |
| `order_exact_of_gen` | ★★★★★★★**生成すれば位数はちょうど `l^n`** |
| `dvd_of_nsmul_eq_zero` | 位数がちょうど `l^n` なら `m·P = 0 ⟹ l^n ∣ m` |
| `nsmul_eq_zero_of_pow_eq_one` | ★★★★★★非退化性で `ζ^m = 1 ⟹ m·P = 0` |
| `isPrimitiveRoot_weilPairingVal` | ★★★★★★★★**Weil 対の値は原始 `l^n` 乗根** |
| `exists_pow_weilPairingVal_eq` | ★★★★★★★★★**`hgen` そのもの** |
| `hgen_of_nondeg` | ★★★★★★★★★**`hgen` は非退化性だけから出る** |
| `det_cyclotomic_of_nondeg` | ★★★★★★★★★★**`det ρ = 円分指標`(非退化性だけ)** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

universe u

/-! ## ★★★★★数え上げの段 -/

/-- ★★★★★2 元で生成される捩れ群の位数の上限。 -/
theorem card_torsion_le_of_two_gen {A : Type u} [AddCommGroup A] (N : ℕ) (P Q : A) (s t : ℕ)
    (hs : 0 < s) (ht : 0 < t) (hP : s • P = 0) (hQ : t • Q = 0)
    (hPN : N • P = 0) (hQN : N • Q = 0)
    (hgen : ∀ R : A, N • R = 0 → ∃ a b : ℕ, R = a • P + b • Q) :
    Nat.card {R : A // N • R = 0} ≤ s * t := by
  have hmem : ∀ a b : ℕ, N • ((a • P + b • Q : A)) = 0 := by
    intro a b
    rw [smul_add, smul_comm, hPN, smul_comm, hQN, smul_zero, smul_zero, add_zero]
  set f : Fin s × Fin t → {R : A // N • R = 0} :=
    fun p => ⟨((p.1 : ℕ) • P + (p.2 : ℕ) • Q : A), hmem _ _⟩ with hf
  have hsurj : Function.Surjective f := by
    rintro ⟨R, hR⟩
    obtain ⟨a, b, hab⟩ := hgen R hR
    refine ⟨(⟨a % s, Nat.mod_lt _ hs⟩, ⟨b % t, Nat.mod_lt _ ht⟩), ?_⟩
    refine Subtype.ext ?_
    show (a % s) • P + (b % t) • Q = R
    rw [nsmul_mod s hP, nsmul_mod t hQ, hab]
  have := Nat.card_le_card_of_surjective f hsurj
  simpa using this

/-- ★★★★★★★**生成すれば位数はちょうど `l^n`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ℤ_l` 線型性は要らない——`#E[l^n] = l^{2n}` の数え上げだけで自由性が出る。 -/
theorem order_exact_of_gen {A : Type u} [AddCommGroup A] (l n : ℕ) (hl : 2 ≤ l) (hn : 1 ≤ n)
    (P Q : A)
    (hcard : Nat.card {R : A // (l ^ n) • R = 0} = (l ^ n) ^ 2)
    (hPN : (l ^ n) • P = 0) (hQN : (l ^ n) • Q = 0)
    (hgen : ∀ R : A, (l ^ n) • R = 0 → ∃ a b : ℕ, R = a • P + b • Q) :
    (l ^ (n - 1)) • P ≠ 0 := by
  intro h
  have hpos : 0 < l ^ (n - 1) := pow_pos (by omega) _
  have hpos' : 0 < l ^ n := pow_pos (by omega) _
  have hle := card_torsion_le_of_two_gen (l ^ n) P Q (l ^ (n - 1)) (l ^ n)
    hpos hpos' h hQN hPN hQN hgen
  rw [hcard] at hle
  have he1 : (l ^ n) ^ 2 = l ^ (2 * n) := by rw [← pow_mul]; ring_nf
  have he2 : l ^ (n - 1) * l ^ n = l ^ (2 * n - 1) := by
    rw [← pow_add, show n - 1 + n = 2 * n - 1 by omega]
  rw [he1, he2] at hle
  have hlt : l ^ (2 * n - 1) < l ^ (2 * n) := Nat.pow_lt_pow_right (by omega) (by omega)
  omega

/-- 位数がちょうど `l^n` なら `m·P = 0 ⟹ l^n ∣ m`。 -/
theorem dvd_of_nsmul_eq_zero {A : Type u} [AddCommGroup A] (l n : ℕ) (hl : l.Prime) (hn : 1 ≤ n)
    (P : A) (hPN : (l ^ n) • P = 0) (hPne : (l ^ (n - 1)) • P ≠ 0) {m : ℕ} (hm : m • P = 0) :
    l ^ n ∣ m := by
  have hdvd : addOrderOf P ∣ l ^ n := addOrderOf_dvd_iff_nsmul_eq_zero.2 hPN
  obtain ⟨k, hk, hke⟩ := (Nat.dvd_prime_pow hl).1 hdvd
  have hkn : k = n := by
    by_contra hne
    exact hPne (addOrderOf_dvd_iff_nsmul_eq_zero.1
      (by rw [hke]; exact pow_dvd_pow l (by omega)))
  have hord : addOrderOf P = l ^ n := by rw [hke, hkn]
  rw [← hord]
  exact addOrderOf_dvd_iff_nsmul_eq_zero.2 hm

/-! ## ★★★★★★★★非退化性から原始根性へ -/

section Pairing

variable {F : Type} [Field F] [DecidableEq F] [CharZero F] [IsAlgClosed F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-- ★★★★★★**非退化性で `ζ^m = 1 ⟹ m·P = 0`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem nsmul_eq_zero_of_pow_eq_one (N : ℕ) (hN : 1 ≤ N) (P Q : W.Point)
    (hPN : N • P = 0) (hQN : N • Q = 0)
    (hgen : ∀ R : W.Point, N • R = 0 → ∃ a b : ℕ, R = a • P + b • Q)
    (hnd : ∀ S : W.Point, N • S = 0 →
      (∀ R : W.Point, N • R = 0 → weilPairingVal W N S R = 1) → S = 0)
    {m : ℕ} (hm : (weilPairingVal W N P Q) ^ m = 1) :
    m • P = 0 := by
  refine hnd (m • P) (by rw [smul_comm, hPN, smul_zero]) ?_
  intro R hR
  obtain ⟨a, b, rfl⟩ := hgen R hR
  have haP : N • (a • P) = 0 := by rw [smul_comm, hPN, smul_zero]
  have hbQ : N • (b • Q) = 0 := by rw [smul_comm, hQN, smul_zero]
  rw [weilPairingVal_nsmul_left W N hN P _ hPN hR m,
    weilPairingVal_add_right W N hN P (a • P) (b • Q) hPN haP hbQ,
    weilPairingVal_nsmul_right W N hN P P hPN hPN a,
    weilPairingVal_nsmul_right W N hN P Q hPN hQN b,
    weilPairingVal_alt W N hN P hPN, one_pow, one_mul, ← pow_mul, mul_comm b m, pow_mul, hm,
    one_pow]

/-- ★★★★★★★★**Weil 対の値は原始 `l^n` 乗根である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem isPrimitiveRoot_weilPairingVal (l n : ℕ) (hl : l.Prime) (hn : 1 ≤ n) (P Q : W.Point)
    (hPN : (l ^ n) • P = 0) (hQN : (l ^ n) • Q = 0) (hPne : (l ^ (n - 1)) • P ≠ 0)
    (hgen : ∀ R : W.Point, (l ^ n) • R = 0 → ∃ a b : ℕ, R = a • P + b • Q)
    (hnd : ∀ S : W.Point, (l ^ n) • S = 0 →
      (∀ R : W.Point, (l ^ n) • R = 0 → weilPairingVal W (l ^ n) S R = 1) → S = 0) :
    IsPrimitiveRoot (weilPairingVal W (l ^ n) P Q) (l ^ n) := by
  have hN : 1 ≤ l ^ n := Nat.one_le_pow _ _ hl.pos
  refine IsPrimitiveRoot.mk (weilPairingVal_pow_eq_one W P Q hQN) ?_
  intro m hm
  exact dvd_of_nsmul_eq_zero l n hl hn P hPN hPne
    (nsmul_eq_zero_of_pow_eq_one W (l ^ n) hN P Q hPN hQN hgen hnd hm)

/-- ★★★★★★★★★**`hgen` そのもの**——1 の `l^n` 乗根はすべて Weil 対の値の冪である。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem exists_pow_weilPairingVal_eq (l n : ℕ) (hl : l.Prime) (hn : 1 ≤ n) (P Q : W.Point)
    (hPN : (l ^ n) • P = 0) (hQN : (l ^ n) • Q = 0) (hPne : (l ^ (n - 1)) • P ≠ 0)
    (hgen : ∀ R : W.Point, (l ^ n) • R = 0 → ∃ a b : ℕ, R = a • P + b • Q)
    (hnd : ∀ S : W.Point, (l ^ n) • S = 0 →
      (∀ R : W.Point, (l ^ n) • R = 0 → weilPairingVal W (l ^ n) S R = 1) → S = 0)
    (z : F) (hz : z ^ (l ^ n) = 1) :
    ∃ k : ℕ, (weilPairingVal W (l ^ n) P Q) ^ k = z := by
  haveI : NeZero (l ^ n) := ⟨by have := Nat.one_le_pow n l hl.pos; omega⟩
  obtain ⟨i, _, hi⟩ := (isPrimitiveRoot_weilPairingVal W l n hl hn P Q hPN hQN hPne hgen
    hnd).eq_pow_of_pow_eq_one hz
  exact ⟨i, hi⟩

end Pairing

/-! ## ★★★★★★★★★★配線を閉じる -/

section Wiring

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★★★★★★★**`hgen` は非退化性だけから出る**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★生成性は第 209、位数の自由性は本ブロックの数え上げ、原始根性は非退化性から。 -/
theorem hgen_of_nondeg [IsAlgClosed L] [CharZero L]
    (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (n : ℕ) (hn : 1 ≤ n)
    (hnd : ∀ S : (W.baseChange L).toAffine.Point, (l ^ n) • S = 0 →
      (∀ R : (W.baseChange L).toAffine.Point, (l ^ n) • R = 0 →
        weilPairingVal (W.baseChange L).toAffine (l ^ n) S R = 1) → S = 0)
    (z : L) (hz : z ^ (l ^ n) = 1) :
    ∃ k : ℕ, (weilPairingVal (W.baseChange L).toAffine (l ^ n)
      (tateVec W l e n (Pi.single 0 1)) (tateVec W l e n (Pi.single 1 1))) ^ k = z := by
  have hl : l.Prime := Fact.out
  have hΔ : (W.baseChange L).Δ ≠ 0 := (W.baseChange L).isUnit_Δ.ne_zero
  have hchar : ∀ k : ℕ, 1 ≤ k → ((k : L) ≠ 0) := fun k hk => Nat.cast_ne_zero.2 (by omega)
  have hN : 1 ≤ l ^ n := Nat.one_le_pow _ _ hl.pos
  have hstep : ∀ (m : ℕ) (T : (W.baseChange L).toAffine.Point), (l ^ m) • T = 0 →
      ∃ S : (W.baseChange L).toAffine.Point, (l ^ (m + 1)) • S = 0 ∧ l • S = T :=
    fun m T hT => exists_smul_step (W.baseChange L) hΔ l m hl.one_lt.le
      (fun k hk1 _ => hchar k hk1) hT
  have hgenR : ∀ R : (W.baseChange L).toAffine.Point, (l ^ n) • R = 0 →
      ∃ a b : ℕ, R = a • tateVec W l e n (Pi.single 0 1)
        + b • tateVec W l e n (Pi.single 1 1) :=
    fun R hR => exists_smul_repr W l e n hstep hR
  have hcard : Nat.card {R : (W.baseChange L).toAffine.Point // (l ^ n) • R = 0} = (l ^ n) ^ 2 :=
    torsion_card (W.baseChange L) hΔ (l ^ n) hN (fun k hk1 _ => hchar k hk1)
  have hPne := order_exact_of_gen l n hl.two_le hn
    (tateVec W l e n (Pi.single 0 1)) (tateVec W l e n (Pi.single 1 1)) hcard
    (tateVec_torsion W l e n _) (tateVec_torsion W l e n _) hgenR
  exact exists_pow_weilPairingVal_eq (W.baseChange L).toAffine l n hl hn _ _
    (tateVec_torsion W l e n _) (tateVec_torsion W l e n _) hPne hgenR hnd z hz

/-- ★★★★★★★★★★**`det ρ = 円分指標`(非退化性だけを仮定した形)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 207 の `hgen` を本ブロックの `hgen_of_nondeg` で消した形である。
★★`n = 0` は `ζ = 1` なので自明に処理でき、`1 ≤ n` の仮定も要らない。 -/
theorem det_cyclotomic_of_nondeg [IsAlgClosed L] [CharZero L]
    (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (n : ℕ) (ζ : L) (hζ : ζ ^ (l ^ n) = 1)
    (hnd : ∀ S : (W.baseChange L).toAffine.Point, (l ^ n) • S = 0 →
      (∀ R : (W.baseChange L).toAffine.Point, (l ^ n) • R = 0 →
        weilPairingVal (W.baseChange L).toAffine (l ^ n) S R = 1) → S = 0) :
    σ ζ = ζ ^ ((PadicInt.toZModPow n
      ((galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det)).val) := by
  rcases Nat.eq_zero_or_pos n with hn0 | hn
  · subst hn0
    rw [pow_zero, pow_one] at hζ
    rw [hζ, map_one, one_pow]
  · exact det_cyclotomic_of_gen W l e σ n ζ hζ
      (fun z hz => hgen_of_nondeg W l e n hn hnd z hz)

end Wiring

/-! ## ★出典の紐付け(`.src`) -/

def order_exact_of_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——基底の位数がちょうど l^n であること)",
    sectionId := "genell-thm-3-8" }

def isPrimitiveRoot_weilPairingVal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の値が原始 l^n 乗根であること)",
    sectionId := "genell-thm-3-8" }

def hgen_of_nondeg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——基底性が非退化性に帰着すること)",
    sectionId := "genell-thm-3-8" }

def det_cyclotomic_of_nondeg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式が円分指標であること——非退化性のみを仮定した形)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
