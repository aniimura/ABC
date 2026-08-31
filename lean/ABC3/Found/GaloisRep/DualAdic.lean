/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicSeries
import Mathlib.Algebra.TrivSqZeroExt.Basic

/-!
# Galois (G6) 第 850 ブロック —— **★★★★★★★★★★`R` の上の双対数は `I` 進完備**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★なぜ双対数か（`tate_ode` への道）

残る葉 `tate_ode`（`DY = 3X² − Y + a₄`）は `tate_equation` を「微分」して出る。
★しかし `R` は任意の完備環であり、微分は存在しない。

★★**双対数を使えば `R` の中で微分できる**:

    `a′ := a + εa`,  `w′ := w − εw`   ⟹  `a′w′ = aw(1 − ε²) = q`  ★★q は動かない

`1 − a′ = (1−a) − εa` は `1 − a` が単元なら単元（単元＋冪零）。よって `tate_equation`
（証明済み）を `R[ε]` の中で `(a′, w′, q)` に適用でき、その `ε` 成分が

    `2(2Y+X)(2DY+DX) = (12X² + 2X + 4a₄)·DX`

を与える。★`DX = 2Y+X`（第 846、証明済み）を入れて `2DX` で割れば `tate_ode` である。

☆実際 `tateXterm(a′) = tateXterm(a) + ε·tateDXterm(a)` が成り立つ:
`(1−a′)⁻¹ = (1−a)⁻¹ + ε·a(1−a)⁻²` なので
`a′(1−a′)⁻² = a(1−a)⁻² + ε[a(1−a)⁻² + 2a²(1−a)⁻³] = f(a) + ε·Df(a)` ✓
（`Df(t) = t(1+t)/(1−t)³` は第 846 の `tateDXterm`）。

## ★★★★★本ブロックの内容——`R[ε]` が `I[ε]` 進完備であること

`tate_equation` は `[IsAdicComplete I R]` を要求するので、まず

    `IsAdicComplete (dualIdeal I) (DualNum R)`

を作る。鍵は **`(dualIdeal I)^n = dualIdeal (I^n)`**（両成分が `I^n` に入る）である。

| 定義・定理 | 内容 |
|---|---|
| `DualNum` | `R[ε]`（`TrivSqZeroExt` を `def` で包む——第 845 の配管） |
| `DualNum.re` / `eps` / `mk` | 成分 |
| `dualIdeal` | `I[ε] = {x ∣ re x ∈ I ∧ eps x ∈ I}` |
| `dualIdeal_mul` | `I[ε]·J[ε] = (IJ)[ε]` |
| `dualIdeal_pow` | ★★★**`I[ε]^n = (I^n)[ε]`** |
| `DualNum.instIsAdicComplete` | ★★★★★★**`R[ε]` は完備** |
-/

namespace ABC3.Found.GaloisRep

/-- ★★★★**双対数 `R[ε]`**（`ε² = 0`）。

★インスタンスの菱形を避けるため `TrivSqZeroExt` を **`def`** で包む（第 845 の配管）。 -/
def DualNum (R : Type) [CommRing R] : Type := TrivSqZeroExt R R

namespace DualNum

variable {R : Type} [CommRing R]

instance : CommRing (DualNum R) := inferInstanceAs (CommRing (TrivSqZeroExt R R))

/-- ★実部。 -/
def re (x : DualNum R) : R := TrivSqZeroExt.fst (show TrivSqZeroExt R R from x)

/-- ★`ε` 部。 -/
def eps (x : DualNum R) : R := TrivSqZeroExt.snd (show TrivSqZeroExt R R from x)

/-- ★`a + εb`。 -/
def mk (a b : R) : DualNum R :=
  show TrivSqZeroExt R R from TrivSqZeroExt.inl a + TrivSqZeroExt.inr b

@[ext] theorem ext {x y : DualNum R} (hre : x.re = y.re) (heps : x.eps = y.eps) : x = y :=
  TrivSqZeroExt.ext hre heps

@[simp] theorem re_mk (a b : R) : (mk a b).re = a := by simp [re, mk]
@[simp] theorem eps_mk (a b : R) : (mk a b).eps = b := by simp [eps, mk]

@[simp] theorem re_zero : (0 : DualNum R).re = 0 := rfl
@[simp] theorem eps_zero : (0 : DualNum R).eps = 0 := rfl
@[simp] theorem re_one : (1 : DualNum R).re = 1 := rfl
@[simp] theorem eps_one : (1 : DualNum R).eps = 0 := rfl

@[simp] theorem re_add (x y : DualNum R) : (x + y).re = x.re + y.re := rfl
@[simp] theorem eps_add (x y : DualNum R) : (x + y).eps = x.eps + y.eps := rfl
@[simp] theorem re_neg (x : DualNum R) : (-x).re = -x.re := rfl
@[simp] theorem eps_neg (x : DualNum R) : (-x).eps = -x.eps := rfl
@[simp] theorem re_sub (x y : DualNum R) : (x - y).re = x.re - y.re := rfl
@[simp] theorem eps_sub (x y : DualNum R) : (x - y).eps = x.eps - y.eps := rfl

theorem mul_eq (x y : DualNum R) :
    (show TrivSqZeroExt R R from x * y)
      = (show TrivSqZeroExt R R from x) * (show TrivSqZeroExt R R from y) := rfl

@[simp] theorem re_mul (x y : DualNum R) : (x * y).re = x.re * y.re := by
  simp only [re, mul_eq, TrivSqZeroExt.fst_mul]

@[simp] theorem eps_mul (x y : DualNum R) : (x * y).eps = x.re * y.eps + y.re * x.eps := by
  simp only [eps, re, mul_eq, TrivSqZeroExt.snd_mul, TrivSqZeroExt.fst_mul]
  simp [smul_eq_mul, mul_comm]

theorem mk_re_eps (x : DualNum R) : mk x.re x.eps = x := by ext <;> simp

end DualNum

/-! ## ★★★★★イデアル `I[ε]` -/

variable {R : Type} [CommRing R]

/-- ★★**`I[ε] = {x ∣ re x ∈ I ∧ eps x ∈ I}`**。 -/
def dualIdeal (I : Ideal R) : Ideal (DualNum R) where
  carrier := {x : DualNum R | x.re ∈ I ∧ x.eps ∈ I}
  add_mem' := by
    rintro x y ⟨hx1, hx2⟩ ⟨hy1, hy2⟩
    exact ⟨Ideal.add_mem I hx1 hy1, Ideal.add_mem I hx2 hy2⟩
  zero_mem' := by
    exact ⟨Ideal.zero_mem I, Ideal.zero_mem I⟩
  smul_mem' := by
    rintro c x ⟨hx1, hx2⟩
    refine ⟨?_, ?_⟩
    · rw [smul_eq_mul, DualNum.re_mul]
      exact Ideal.mul_mem_left I _ hx1
    · rw [smul_eq_mul, DualNum.eps_mul]
      exact Ideal.add_mem I (Ideal.mul_mem_left I _ hx2) (Ideal.mul_mem_right _ I hx1)

@[simp] theorem mem_dualIdeal {I : Ideal R} {x : DualNum R} :
    x ∈ dualIdeal I ↔ x.re ∈ I ∧ x.eps ∈ I := Iff.rfl

theorem dualIdeal_top : dualIdeal (⊤ : Ideal R) = ⊤ := by
  ext x
  simp

/-- ★★★**`I[ε]·J[ε] = (IJ)[ε]`**。 -/
theorem dualIdeal_mul (I J : Ideal R) :
    dualIdeal I * dualIdeal J = dualIdeal (I * J) := by
  refine le_antisymm ?_ ?_
  · refine Ideal.mul_le.2 (fun x hx y hy => ?_)
    rw [mem_dualIdeal] at hx hy ⊢
    refine ⟨Ideal.mul_mem_mul hx.1 hy.1, ?_⟩
    rw [DualNum.eps_mul]
    refine Ideal.add_mem _ (Ideal.mul_mem_mul hx.1 hy.2) ?_
    have hcm : y.re * x.eps = x.eps * y.re := mul_comm _ _
    rw [hcm]
    exact Ideal.mul_mem_mul hx.2 hy.1
  · intro x hx
    rw [mem_dualIdeal] at hx
    have hsplit : x = DualNum.mk x.re 0 + DualNum.mk 0 x.eps := by
      ext <;> simp
    rw [hsplit]
    refine Ideal.add_mem _ ?_ ?_
    · refine Submodule.mul_induction_on hx.1 ?_ ?_
      · intro i hi j hj
        have : DualNum.mk (i * j) 0 = DualNum.mk i 0 * DualNum.mk j 0 := by
          ext <;> simp
        rw [this]
        exact Ideal.mul_mem_mul ⟨by simpa using hi, by simp⟩ ⟨by simpa using hj, by simp⟩
      · intro y z hy hz
        have : DualNum.mk (y + z) 0 = DualNum.mk y 0 + DualNum.mk z 0 := by
          ext <;> simp
        rw [this]
        exact Ideal.add_mem _ hy hz
    · refine Submodule.mul_induction_on hx.2 ?_ ?_
      · intro i hi j hj
        have : DualNum.mk 0 (i * j) = DualNum.mk i 0 * DualNum.mk 0 j := by
          ext <;> simp
        rw [this]
        exact Ideal.mul_mem_mul ⟨by simpa using hi, by simp⟩ ⟨by simp, by simpa using hj⟩
      · intro y z hy hz
        have : DualNum.mk 0 (y + z) = DualNum.mk 0 y + DualNum.mk 0 z := by
          ext <;> simp
        rw [this]
        exact Ideal.add_mem _ hy hz

/-- ★★★★★★**`I[ε]^n = (I^n)[ε]`**。 -/
theorem dualIdeal_pow (I : Ideal R) (n : ℕ) :
    (dualIdeal I) ^ n = dualIdeal (I ^ n) := by
  induction n with
  | zero => simpa using (dualIdeal_top (R := R)).symm
  | succ m ih => rw [pow_succ, ih, dualIdeal_mul, pow_succ]

/-! ## ★★★★★★★★`R[ε]` は完備 -/

theorem mem_smul_top_self {J : Ideal R} {x : R} :
    x ∈ (J • (⊤ : Submodule R R)) ↔ x ∈ J := by simp

theorem mem_dual_smul_top {I : Ideal R} {n : ℕ} {x : DualNum R} :
    x ∈ ((dualIdeal I) ^ n • (⊤ : Submodule (DualNum R) (DualNum R)))
      ↔ x.re ∈ I ^ n ∧ x.eps ∈ I ^ n := by
  have h : ((dualIdeal I) ^ n • (⊤ : Submodule (DualNum R) (DualNum R)))
      = ((dualIdeal I) ^ n : Ideal (DualNum R)) := by simp
  rw [h, dualIdeal_pow]
  exact mem_dualIdeal

instance instIsHausdorffDual (I : Ideal R) [IsHausdorff I R] :
    IsHausdorff (dualIdeal I) (DualNum R) where
  haus' := by
    intro x hx
    have hcomp : ∀ n : ℕ, x.re ∈ I ^ n ∧ x.eps ∈ I ^ n := by
      intro n
      have h := hx n
      rw [SModEq.sub_mem, sub_zero, mem_dual_smul_top] at h
      exact h
    have h1 : x.re = 0 := by
      refine IsHausdorff.haus' (I := I) _ (fun n => ?_)
      rw [SModEq.sub_mem, sub_zero, mem_smul_top_self]
      exact (hcomp n).1
    have h2 : x.eps = 0 := by
      refine IsHausdorff.haus' (I := I) _ (fun n => ?_)
      rw [SModEq.sub_mem, sub_zero, mem_smul_top_self]
      exact (hcomp n).2
    ext
    · simpa using h1
    · simpa using h2

instance instIsPrecompleteDual (I : Ideal R) [IsPrecomplete I R] :
    IsPrecomplete (dualIdeal I) (DualNum R) where
  prec' := by
    intro f hf
    have hcomp : ∀ {m n : ℕ}, m ≤ n →
        ((f m).re - (f n).re ∈ I ^ m ∧ (f m).eps - (f n).eps ∈ I ^ m) := by
      intro m n hmn
      have h := hf hmn
      rw [SModEq.sub_mem, mem_dual_smul_top] at h
      simpa using h
    have hre : ∀ {m n : ℕ}, m ≤ n →
        (f m).re ≡ (f n).re [SMOD (I ^ m • ⊤ : Submodule R R)] := by
      intro m n hmn
      rw [SModEq.sub_mem, mem_smul_top_self]
      exact (hcomp hmn).1
    have heps : ∀ {m n : ℕ}, m ≤ n →
        (f m).eps ≡ (f n).eps [SMOD (I ^ m • ⊤ : Submodule R R)] := by
      intro m n hmn
      rw [SModEq.sub_mem, mem_smul_top_self]
      exact (hcomp hmn).2
    obtain ⟨L1, hL1⟩ := IsPrecomplete.prec' (I := I) (fun n => (f n).re) hre
    obtain ⟨L2, hL2⟩ := IsPrecomplete.prec' (I := I) (fun n => (f n).eps) heps
    refine ⟨DualNum.mk L1 L2, fun n => ?_⟩
    have h1 := hL1 n
    have h2 := hL2 n
    rw [SModEq.sub_mem, mem_smul_top_self] at h1 h2
    rw [SModEq.sub_mem, mem_dual_smul_top]
    constructor
    · simpa using h1
    · simpa using h2

instance instIsAdicCompleteDual (I : Ideal R) [IsAdicComplete I R] :
    IsAdicComplete (dualIdeal I) (DualNum R) where

/-! ## ★★★★★★`a′ = a + εa`、`w′ = w − εw` -/

theorem DualNum.isUnit_iff {x : DualNum R} : IsUnit x ↔ IsUnit x.re :=
  TrivSqZeroExt.isUnit_iff_isUnit_fst (x := show TrivSqZeroExt R R from x)

/-- ★★★★★★★★**`a′w′ = q`**——`ε` 成分は `a(−w) + wa = 0`。
★★**`q` は動かない**ことがこの変形の心臓である。 -/
theorem dualA_mul_dualW (a w q : R) (haw : a * w = q) :
    DualNum.mk a a * DualNum.mk w (-w) = DualNum.mk q 0 := by
  ext
  · simpa using haw
  · simp
    ring

theorem mk_mem_dualIdeal {I : Ideal R} {q : R} (hq : q ∈ I) :
    DualNum.mk q 0 ∈ dualIdeal I := ⟨by simpa using hq, by simp⟩

/-- ★★★`1 − a′` は `1 − a` が単元なら単元（単元＋冒零）。 -/
theorem isUnit_one_sub_dualA {a : R} (ha : IsUnit (1 - a)) :
    IsUnit (1 - DualNum.mk a a) := by
  rw [DualNum.isUnit_iff]
  simpa using ha

theorem isUnit_one_sub_dualW {w : R} (hw : IsUnit (1 - w)) :
    IsUnit (1 - DualNum.mk w (-w)) := by
  rw [DualNum.isUnit_iff]
  simpa using hw

/-! ## ★出典の紐付け(`.src`) -/

def DualNum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(双対数 R[ε]。★無条件)",
    sectionId := "genell-lemma-3-2" }

def dualIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(イデアル I[ε]。★無条件)",
    sectionId := "genell-lemma-3-2" }

def dualIdeal_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(I[ε]^n = (I^n)[ε]。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
