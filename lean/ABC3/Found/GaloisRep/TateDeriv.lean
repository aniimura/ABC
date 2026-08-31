/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateUniversal
import Mathlib.Algebra.MvPolynomial.PDeriv
import Mathlib.Algebra.TrivSqZeroExt.Basic

/-!
# Galois (G6) 第 845 ブロック —— **★★★★★★★★万有な環の上の微分 `D`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★なぜ微分が要るのか（第 842 の道）

葉 1（Vélu の商が `E_{q^l}` であること）の残りは `∑_ζ X(ζ)²` の値であり、
係数を直に比べると **Besge の畳み込み恒等式**に落ちる。
★★微分 `D = u∂_u`（`q` を固定）を使えば、畳み込みを**経由せずに**済む:

    `(DX)² = 4X³ + X² + 4a₄X + 4a₆`   （= `tate_equation`、証明済み）
    ↓ `D` を適用して `2DX` で割る
    `D²X = 6X² + X + 2a₄`

## ★★★★★なぜ万有な環の上で作るのか

`tateXpair a w q hq : R` の `R` は**任意の** `I` 進完備環であり、
そこに `u` など無い——`D` は存在しない。
★しかし**万有な環** `TateBase = ℤ[A, W]`（`q = AW`）の上には

    `D₀ = A∂_A − W∂_W`

があり、★★**`D₀ q = A·(−W) + W·A = 0`**——つまり `q` は `D` の定数である。
（`A = u`、`W = u⁻¹` と思えば `D₀ = u∂_u` そのもの。）

## ★★★★★★★局所化への延長は**双対数**で行う

mathlib には「局所化への微分の延長」は無い
（`Mathlib/RingTheory/Derivation/` を `Localization` で grep して 0 件、2026-08-31）。
★そこで**双対数のトリック**を使う:

    `Φ₀ : TateBase →+* TateUniv[ε]`,  `f ↦ ι f + ε·ι(D₀ f)`

は Leibniz 則により**環準同型**であり、★★分母 `s` の行先は
`fst = ι s` が単元なので**自動的に単元**である（`isUnit_iff_isUnit_fst`）。
よって `IsLocalization.lift` が `Φ : TateUniv →+* TateUniv[ε]` を与え、

    `D x := (Φ x).snd`

が求める微分である。★★★**商の規則 `D(a/s) = (sDa − aDs)/s²` を手で書かずに済む**。

## ★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `univD0` | ★★`A∂_A − W∂_W`（`TateBase` の微分） |
| `univD0_A` / `univD0_W` / `univD0_Q` | `DA = A`、`DW = −W`、★★**`Dq = 0`** |
| `univDual` | 双対数への環準同型 |
| `univDualLift` | ★局所化への持ち上げ |
| `univD` | ★★★★★★**`TateUniv` の上の微分** |
| `univD_mul` | ★★★Leibniz 則 |
| `univD_inv` | ★★商の規則（逆元の微分） |
-/

namespace ABC3.Found.GaloisRep

open MvPolynomial

/-! ## ★★★★★`TateBase` の上の微分 -/

/-- ★★★★★★**万有な微分 `D₀ = A∂_A − W∂_W`**。

`A = u`、`W = u⁻¹` と思えば `u∂_u` そのものであり、`q = AW` を死に履く。 -/
noncomputable def univD0 : Derivation ℤ TateBase TateBase :=
  univA • (pderiv (0 : Fin 2)) - univW • (pderiv (1 : Fin 2))

theorem univD0_apply (f : TateBase) :
    univD0 f = univA * pderiv 0 f - univW * pderiv 1 f := by
  simp only [univD0, Derivation.sub_apply]
  rfl

@[simp] theorem univD0_A : univD0 univA = univA := by
  have h0 : pderiv (0 : Fin 2) (univA : TateBase) = 1 := by
    rw [univA]; exact pderiv_X_self 0
  have h1 : pderiv (1 : Fin 2) (univA : TateBase) = 0 := by
    rw [univA]; exact pderiv_X_of_ne (by decide)
  rw [univD0_apply, h0, h1]
  ring

@[simp] theorem univD0_W : univD0 univW = -univW := by
  have h0 : pderiv (0 : Fin 2) (univW : TateBase) = 0 := by
    rw [univW]; exact pderiv_X_of_ne (by decide)
  have h1 : pderiv (1 : Fin 2) (univW : TateBase) = 1 := by
    rw [univW]; exact pderiv_X_self 1
  rw [univD0_apply, h0, h1]
  ring

/-- ★★★★★★★★**`q` は `D₀` の定数**——これがこの微分を使う理由である。 -/
@[simp] theorem univD0_Q : univD0 univQ = 0 := by
  rw [univQ, Derivation.leibniz, univD0_A, univD0_W, smul_eq_mul, smul_eq_mul]
  ring

/-! ## ★★★★★★★双対数 -/

/-- ★★★★**双対数 `TateUniv[ε]`**（`ε² = 0`）。

★インスタンスの菱形を避けるため `TrivSqZeroExt` を**新しい型**（`def`）で包む。 -/
def TateDual : Type := TrivSqZeroExt TateUniv TateUniv

noncomputable instance : CommRing TateDual :=
  inferInstanceAs (CommRing (TrivSqZeroExt TateUniv TateUniv))

namespace TateDual

/-- ★実部。 -/
noncomputable def re (x : TateDual) : TateUniv :=
  TrivSqZeroExt.fst (show TrivSqZeroExt TateUniv TateUniv from x)

/-- ★`ε` 部。 -/
noncomputable def eps (x : TateDual) : TateUniv :=
  TrivSqZeroExt.snd (show TrivSqZeroExt TateUniv TateUniv from x)

/-- ★`a + εb`。 -/
noncomputable def mk (a b : TateUniv) : TateDual :=
  show TrivSqZeroExt TateUniv TateUniv from
    TrivSqZeroExt.inl a + TrivSqZeroExt.inr b

@[simp] theorem re_mk (a b : TateUniv) : (mk a b).re = a := by
  simp [re, mk]

@[simp] theorem eps_mk (a b : TateUniv) : (mk a b).eps = b := by
  simp [eps, mk]

@[simp] theorem re_zero : (0 : TateDual).re = 0 := rfl
@[simp] theorem eps_zero : (0 : TateDual).eps = 0 := rfl
@[simp] theorem re_one : (1 : TateDual).re = 1 := rfl
@[simp] theorem eps_one : (1 : TateDual).eps = 0 := rfl

@[simp] theorem re_add (x y : TateDual) : (x + y).re = x.re + y.re := rfl
@[simp] theorem eps_add (x y : TateDual) : (x + y).eps = x.eps + y.eps := rfl
@[simp] theorem re_neg (x : TateDual) : (-x).re = -x.re := rfl
@[simp] theorem eps_neg (x : TateDual) : (-x).eps = -x.eps := rfl
@[simp] theorem re_sub (x y : TateDual) : (x - y).re = x.re - y.re := rfl
@[simp] theorem eps_sub (x y : TateDual) : (x - y).eps = x.eps - y.eps := rfl

theorem mul_eq (x y : TateDual) :
    (show TrivSqZeroExt TateUniv TateUniv from x * y)
      = (show TrivSqZeroExt TateUniv TateUniv from x) *
        (show TrivSqZeroExt TateUniv TateUniv from y) := rfl

@[simp] theorem re_mul (x y : TateDual) : (x * y).re = x.re * y.re := by
  simp only [re, mul_eq, TrivSqZeroExt.fst_mul]

@[simp] theorem eps_mul (x y : TateDual) : (x * y).eps = x.re * y.eps + y.re * x.eps := by
  simp only [eps, re, mul_eq, TrivSqZeroExt.snd_mul, TrivSqZeroExt.fst_mul]
  simp [smul_eq_mul, mul_comm]

@[simp] theorem re_intCast (a : ℤ) : ((a : TateDual)).re = (a : TateUniv) := by
  induction a using Int.induction_on with
  | zero => rfl
  | succ k ih => push_cast; push_cast at ih; simp [ih]
  | pred k ih => push_cast; push_cast at ih; simp [ih]

@[simp] theorem eps_intCast (a : ℤ) : ((a : TateDual)).eps = 0 := by
  induction a using Int.induction_on with
  | zero => rfl
  | succ k ih => push_cast; push_cast at ih; simp [ih]
  | pred k ih => push_cast; push_cast at ih; simp [ih]

end TateDual

/-! ## ★★★★★★★双対数への持ち上げ -/

/-- ★`TateBase → TateUniv` の構造写像。 -/
noncomputable abbrev univIota : TateBase →+* TateUniv := algebraMap TateBase TateUniv

/-- ★`A` の行先 `A + εA`。 -/
noncomputable def dualA : TateDual := TateDual.mk (univIota univA) (univIota univA)

/-- ★`W` の行先 `W − εW`。 -/
noncomputable def dualW : TateDual := TateDual.mk (univIota univW) (-(univIota univW))

/-- ★★★★★**双対数への環準同型** `f ↦ ι f + ε·ι(D₀ f)`。

★`eval₂Hom` で作るのは、戻り値の型のインスタンスを `IsLocalization.lift` と
揃えるためである（本ファイル冒頭の「配管」を参照）。 -/
noncomputable def univDual :=
  MvPolynomial.eval₂Hom (Int.castRingHom TateDual) ![dualA, dualW]

@[simp] theorem univDual_A : univDual univA = dualA := by
  rw [univDual, univA, MvPolynomial.eval₂Hom_X']
  rfl

@[simp] theorem univDual_W : univDual univW = dualW := by
  rw [univDual, univW, MvPolynomial.eval₂Hom_X']
  rfl

@[simp] theorem univDual_C (a : ℤ) : univDual (MvPolynomial.C a) = (a : TateDual) := by
  rw [univDual, MvPolynomial.eval₂Hom_C]
  simp

theorem univIota_C (a : ℤ) : univIota (MvPolynomial.C a) = (a : TateUniv) := by
  have h : (MvPolynomial.C a : TateBase) = (a : TateBase) := by simp
  rw [h, map_intCast]

/-- ★生成元の上での実部。 -/
theorem re_univDual_X : ∀ i : Fin 2,
    (univDual (MvPolynomial.X i)).re = univIota (MvPolynomial.X i)
  | 0 => by
      rw [show (MvPolynomial.X (0 : Fin 2) : TateBase) = univA from rfl, univDual_A]
      simp [dualA]
  | 1 => by
      rw [show (MvPolynomial.X (1 : Fin 2) : TateBase) = univW from rfl, univDual_W]
      simp [dualW]

/-- ★生成元の上での `ε` 部。 -/
theorem eps_univDual_X : ∀ i : Fin 2,
    (univDual (MvPolynomial.X i)).eps = univIota (univD0 (MvPolynomial.X i))
  | 0 => by
      rw [show (MvPolynomial.X (0 : Fin 2) : TateBase) = univA from rfl, univDual_A,
        univD0_A]
      simp [dualA]
  | 1 => by
      rw [show (MvPolynomial.X (1 : Fin 2) : TateBase) = univW from rfl, univDual_W,
        univD0_W]
      simp [dualW]

/-- ★★★★**実部は `ι` である**。 -/
theorem re_univDual (f : TateBase) : (univDual f).re = univIota f := by
  induction f using MvPolynomial.induction_on with
  | C a => rw [univDual_C, TateDual.re_intCast, univIota_C]
  | add p q hp hq => rw [map_add, TateDual.re_add, hp, hq, map_add]
  | mul_X p i hp => rw [map_mul, TateDual.re_mul, hp, re_univDual_X, map_mul]

/-- ★★★★★★**`ε` 部は `ι ∘ D₀` である**。 -/
theorem eps_univDual (f : TateBase) : (univDual f).eps = univIota (univD0 f) := by
  induction f using MvPolynomial.induction_on with
  | C a =>
      rw [univDual_C, TateDual.eps_intCast]
      have hC : (MvPolynomial.C a : TateBase) = algebraMap ℤ TateBase a := rfl
      rw [hC, Derivation.map_algebraMap, map_zero]
  | add p q hp hq => rw [map_add, TateDual.eps_add, hp, hq, map_add, map_add]
  | mul_X p i hp =>
      rw [map_mul, TateDual.eps_mul, re_univDual, re_univDual, hp, eps_univDual_X,
        Derivation.leibniz, map_add]
      simp [smul_eq_mul, map_mul]

/-- ★★★★**分母の行先は自動的に単元**である——`re` が単元だから。 -/
theorem univDual_isUnit (y : tateDenoms) : IsUnit (univDual (y : TateBase)) := by
  have h : IsUnit ((univDual (y : TateBase)).re) := by
    rw [re_univDual]; exact IsLocalization.map_units TateUniv y
  exact (TrivSqZeroExt.isUnit_iff_isUnit_fst
    (x := show TrivSqZeroExt TateUniv TateUniv from univDual (y : TateBase))).2 h

/-- ★★★★★**局所化への持ち上げ**。 -/
noncomputable def univDualLift :=
  IsLocalization.lift (M := tateDenoms) (S := TateUniv) univDual_isUnit

@[simp] theorem univDualLift_algebraMap (f : TateBase) :
    univDualLift (univIota f) = univDual f :=
  IsLocalization.lift_eq univDual_isUnit f

/-! ## ★★★★★★★★`TateUniv` の上の微分 -/

/-- ★★★★★★**`re` 成分は恒等写像**。 -/
@[simp] theorem re_univDualLift (x : TateUniv) : (univDualLift x).re = x := by
  obtain ⟨⟨a, s⟩, hz⟩ := IsLocalization.surj (M := tateDenoms) x
  have h1 : univDualLift x * univDual (s : TateBase) = univDual a := by
    rw [← univDualLift_algebraMap (s : TateBase), ← map_mul, hz, univDualLift_algebraMap]
  have h2 := congrArg TateDual.re h1
  rw [TateDual.re_mul, re_univDual, re_univDual] at h2
  refine (IsLocalization.map_units TateUniv s).mul_left_injective ?_
  show (univDualLift x).re * univIota (s : TateBase) = x * univIota (s : TateBase)
  rw [h2, hz]

/-- ★★★★★★★★★★**万有な環の上の微分 `D`**。 -/
noncomputable def univD (x : TateUniv) : TateUniv := (univDualLift x).eps

@[simp] theorem univD_algebraMap (f : TateBase) : univD (univIota f) = univIota (univD0 f) := by
  rw [univD, univDualLift_algebraMap, eps_univDual]

@[simp] theorem univD_zero : univD 0 = 0 := by
  rw [univD, map_zero]; rfl

@[simp] theorem univD_one : univD 1 = 0 := by
  rw [univD, map_one]; rfl

theorem univD_add (x y : TateUniv) : univD (x + y) = univD x + univD y := by
  rw [univD, map_add, TateDual.eps_add]; rfl

theorem univD_neg (x : TateUniv) : univD (-x) = - univD x := by
  rw [univD, map_neg, TateDual.eps_neg]; rfl

theorem univD_sub (x y : TateUniv) : univD (x - y) = univD x - univD y := by
  rw [sub_eq_add_neg, univD_add, univD_neg, sub_eq_add_neg]

/-- ★★★★★★★★**Leibniz 則**。 -/
theorem univD_mul (x y : TateUniv) : univD (x * y) = x * univD y + y * univD x := by
  have h : univDualLift (x * y) = univDualLift x * univDualLift y := map_mul _ _ _
  have h2 := congrArg TateDual.eps h
  rw [TateDual.eps_mul, re_univDualLift, re_univDualLift] at h2
  exact h2

theorem univD_pow (x : TateUniv) (n : ℕ) :
    univD (x ^ (n + 1)) = ((n : TateUniv) + 1) * x ^ n * univD x := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [pow_succ, univD_mul, ih]
      push_cast
      rw [pow_succ]
      ring

/-- ★★★★★★**商の規則**——逆元の微分。 -/
theorem univD_inv (x y : TateUniv) (h : x * y = 1) :
    univD y * x ^ 2 = - univD x := by
  have h0 : x * univD y + y * univD x = 0 := by
    have h1 := univD_mul x y
    rw [h, univD_one] at h1
    exact h1.symm
  have h1 : x * univD y = - (y * univD x) := by linear_combination h0
  calc univD y * x ^ 2 = x * (x * univD y) := by ring
    _ = x * (- (y * univD x)) := by rw [h1]
    _ = - ((x * y) * univD x) := by ring
    _ = - univD x := by rw [h]; ring

/-- ★★★★**`q` は `D` の定数**。 -/
@[simp] theorem univD_Q : univD (univIota univQ) = 0 := by
  rw [univD_algebraMap, univD0_Q, map_zero]

/-! ## ★出典の紐付け(`.src`) -/

def univD0.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(万有な微分 D₀ = A∂_A − W∂_W。★無条件)",
    sectionId := "genell-lemma-3-2" }

def univD0_Q.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(q = AW は D₀ の定数。★無条件)",
    sectionId := "genell-lemma-3-2" }

def univDual.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(双対数への環準同型。★無条件)",
    sectionId := "genell-lemma-3-2" }

def univD.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(TateUniv の上の微分 D。★無条件)",
    sectionId := "genell-lemma-3-2" }

def univD_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D の Leibniz 則。★無条件)",
    sectionId := "genell-lemma-3-2" }

def univD_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D の商の規則。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
