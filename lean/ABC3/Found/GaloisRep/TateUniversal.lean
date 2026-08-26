import ABC3.Found.GaloisRep.TateTruncate
import Mathlib.RingTheory.Localization.Away.Basic
import Mathlib.Algebra.MvPolynomial.CommRing

/-!
# Galois (G6) 第 223 ブロック —— **★★★★★★★★万有な環と特殊化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★万有な環を作る

第 222 で葉 (b) は「切り詰め `tateDefectTrunc n a w q` が `I^n` に入る」に落ちた。
切り詰めは**有限の式**(`Ring.inverse` を除けば多項式)なので、**万有な環**を作って
そこで一度示せばよい。

    TateBase  := ℤ[A, W]                          (`MvPolynomial (Fin 2) ℤ`)
    TateUniv  := TateBase を次の元で局所化したもの
                 1 − A、1 − W、1 − (AW)^{m+1}A、1 − (AW)^{m+1}W  (m ≥ 0)

★★**完備化は要らない**——切り詰めが有限式だからである。局所化だけで足りる。

## ★★★★★特殊化はどんな完備環にも届く

`R` が `I` 進完備で `a·w = q ∈ I`、`1 − a` と `1 − w` が単元なら、
分母のすべてが `R` で単元になる:

| 分母 | `R` での単元性 |
|---|---|
| `1 − A ↦ 1 − a` | ★仮定 |
| `1 − W ↦ 1 − w` | ★仮定(`w ∈ I` なら自動) |
| `1 − (AW)^{m+1}A ↦ 1 − q^{m+1}a` | ★★`q^{m+1}a ∈ I` なので自動(`isUnit_one_sub`) |
| `1 − (AW)^{m+1}W ↦ 1 − q^{m+1}w` | ★★同上 |

よって `IsLocalization.lift` で `TateUniv →+* R` が得られる(`tateSpecialize`)。

## ★★★★★★★★到達点

> **`TateUniv` の中で `(AW)^n ∣ tateDefectTrunc n A W (AW)` を示せば、
> 任意の完備環で Weierstrass 方程式が成り立つ。**

これが `tateDefect_eq_zero_of_univ` である。★葉 (b) は**万有な環での整除性ひとつ**に
帰着した。あとは `ℤ[A,W]` が UFD であることを使い、`A^n ∣` と `W^n ∣` に分けて、
解析側の恒等式(第 220)から係数の消滅を出す。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `TateBase` / `univA` / `univW` / `univQ` | 万有な多項式環と生成元 |
| `tateDenomSet` / `tateDenoms` / `TateUniv` | ★★★★★**万有な環** |
| `tateEval` | `A ↦ a`、`W ↦ w` の評価 |
| `tateEval_isUnit_denoms` | ★★★★分母はすべて `R` で単元になる |
| `tateSpecialize` | ★★★★★★**特殊化 `TateUniv →+* R`** |
| `isUnit_one_sub_uA` 他 | 万有な環の側での単元性 |
| `tateDefect_eq_zero_of_univ` | ★★★★★★★★**葉 (b) は整除性ひとつに帰着** |
-/

namespace ABC3.Found.GaloisRep

open MvPolynomial

/-! ## ★★★★★万有な環 -/

/-- ★万有な多項式環 `ℤ[A, W]`。 -/
abbrev TateBase : Type := MvPolynomial (Fin 2) ℤ

/-- ★万有な `A`。 -/
noncomputable def univA : TateBase := MvPolynomial.X 0
/-- ★万有な `W`。 -/
noncomputable def univW : TateBase := MvPolynomial.X 1
/-- ★万有な `q = AW`。 -/
noncomputable def univQ : TateBase := univA * univW

/-- ★分母となる元の集合——これらが単元になる局所化を取る。 -/
noncomputable def tateDenomSet : Set TateBase :=
  {1 - univA, 1 - univW} ∪ (Set.range fun m : ℕ => 1 - univQ ^ (m + 1) * univA)
    ∪ (Set.range fun m : ℕ => 1 - univQ ^ (m + 1) * univW)

noncomputable def tateDenoms : Submonoid TateBase := Submonoid.closure tateDenomSet

/-- ★★★★★**Tate 級数の万有な環**——`ℤ[A,W]` を分母で局所化したもの。

★完備化は要らない(切り詰めが有限式だから)。 -/
abbrev TateUniv : Type := Localization tateDenoms

/-! ## ★★★★特殊化 -/

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★万有な環からの評価 `A ↦ a`、`W ↦ w`。 -/
noncomputable def tateEval (a w : R) : TateBase →+* R :=
  MvPolynomial.eval₂Hom (Int.castRingHom R) ![a, w]

@[simp] theorem tateEval_A (a w : R) : tateEval a w univA = a := by
  rw [tateEval, univA, MvPolynomial.eval₂Hom_X']
  rfl

@[simp] theorem tateEval_W (a w : R) : tateEval a w univW = w := by
  rw [tateEval, univW, MvPolynomial.eval₂Hom_X']
  rfl

@[simp] theorem tateEval_Q (a w : R) : tateEval a w univQ = a * w := by
  rw [univQ, map_mul, tateEval_A, tateEval_W]

/-- ★★★★**分母の生成元はすべて `R` で単元になる**。 -/
theorem tateEval_isUnit_denomSet [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    {x : TateBase} (hx : x ∈ tateDenomSet) : IsUnit (tateEval a w x) := by
  rcases hx with (hx | ⟨m, rfl⟩) | ⟨m, rfl⟩
  · rcases hx with rfl | rfl
    · simpa using ha
    · simpa using hw
  · have hval : tateEval a w (1 - univQ ^ (m + 1) * univA) = 1 - q ^ (m + 1) * a := by
      simp [haw]
    rw [hval]
    exact isUnit_one_sub (I := I) (pow_succ_mul_mem hq m)
  · have hval : tateEval a w (1 - univQ ^ (m + 1) * univW) = 1 - q ^ (m + 1) * w := by
      simp [haw]
    rw [hval]
    exact isUnit_one_sub (I := I) (pow_succ_mul_mem hq m)

theorem tateEval_isUnit_denoms [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (y : tateDenoms) : IsUnit (tateEval a w (y : TateBase)) := by
  obtain ⟨x, hx⟩ := y
  refine Submonoid.closure_induction (s := tateDenomSet)
    (motive := fun z _ => IsUnit (tateEval a w z)) ?_ ?_ ?_ hx
  · intro z hz
    exact tateEval_isUnit_denomSet a w q hq haw ha hw hz
  · simp
  · intro z₁ z₂ _ _ h₁ h₂
    rw [map_mul]
    exact h₁.mul h₂

/-- ★★★★★★**万有な環から任意の完備環への特殊化**。 -/
noncomputable def tateSpecialize [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) : TateUniv →+* R :=
  IsLocalization.lift (tateEval_isUnit_denoms (I := I) a w q hq haw ha hw)

theorem tateSpecialize_base [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) (x : TateBase) :
    tateSpecialize a w q hq haw ha hw (algebraMap TateBase TateUniv x) = tateEval a w x :=
  IsLocalization.lift_eq _ x

/-! ## ★万有な環の側の元と単元性 -/

/-- ★万有な環の中の `A`。 -/
noncomputable def uA : TateUniv := algebraMap TateBase TateUniv univA
/-- ★万有な環の中の `W`。 -/
noncomputable def uW : TateUniv := algebraMap TateBase TateUniv univW

theorem uA_mul_uW : uA * uW = algebraMap TateBase TateUniv univQ := by
  rw [uA, uW, ← map_mul, univQ]

theorem mem_tateDenoms_of_mem (x : TateBase) (hx : x ∈ tateDenomSet) : x ∈ tateDenoms :=
  Submonoid.subset_closure hx

theorem isUnit_one_sub_uA : IsUnit (1 - uA) := by
  have h : (1 : TateUniv) - uA = algebraMap TateBase TateUniv (1 - univA) := by
    rw [map_sub, map_one, uA]
  rw [h]
  exact IsLocalization.map_units TateUniv
    ⟨1 - univA, mem_tateDenoms_of_mem _ (by left; left; left; rfl)⟩

theorem isUnit_one_sub_uW : IsUnit (1 - uW) := by
  have h : (1 : TateUniv) - uW = algebraMap TateBase TateUniv (1 - univW) := by
    rw [map_sub, map_one, uW]
  rw [h]
  exact IsLocalization.map_units TateUniv
    ⟨1 - univW, mem_tateDenoms_of_mem _ (by left; left; right; rfl)⟩

theorem isUnit_one_sub_uQA (m : ℕ) : IsUnit (1 - (uA * uW) ^ (m + 1) * uA) := by
  have h : (1 : TateUniv) - (uA * uW) ^ (m + 1) * uA
      = algebraMap TateBase TateUniv (1 - univQ ^ (m + 1) * univA) := by
    rw [map_sub, map_one, map_mul, map_pow, ← uA_mul_uW, uA]
  rw [h]
  exact IsLocalization.map_units TateUniv
    ⟨_, mem_tateDenoms_of_mem _ (by left; right; exact ⟨m, rfl⟩)⟩

theorem isUnit_one_sub_uQW (m : ℕ) : IsUnit (1 - (uA * uW) ^ (m + 1) * uW) := by
  have h : (1 : TateUniv) - (uA * uW) ^ (m + 1) * uW
      = algebraMap TateBase TateUniv (1 - univQ ^ (m + 1) * univW) := by
    rw [map_sub, map_one, map_mul, map_pow, ← uA_mul_uW, uW]
  rw [h]
  exact IsLocalization.map_units TateUniv
    ⟨_, mem_tateDenoms_of_mem _ (by right; exact ⟨m, rfl⟩)⟩

/-! ## ★★★★★★★★葉 (b) の帰着 -/

theorem tateSpecialize_uA [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateSpecialize a w q hq haw ha hw uA = a := by
  rw [uA, tateSpecialize_base, tateEval_A]

theorem tateSpecialize_uW [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateSpecialize a w q hq haw ha hw uW = w := by
  rw [uW, tateSpecialize_base, tateEval_W]

/-- ★★★★★★★★**葉 (b) は万有な環での整除性ひとつに帰着した**。

`ℤ[A,W]` を `1−A`・`1−W`・`1−(AW)^{m+1}A`・`1−(AW)^{m+1}W` で局所化した環の中で

    (AW)^n ∣ tateDefectTrunc n A W (AW)

を示せば、**任意の完備環で Weierstrass 方程式が成り立つ**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateDefect_eq_zero_of_univ [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (H : ∀ n : ℕ, ((uA * uW) ^ n) ∣ tateDefectTrunc n uA uW (uA * uW)) :
    tateDefect a w q hq = 0 :=
  tateDefect_eq_zero_of_specialize a w q hq haw uA uW (tateSpecialize a w q hq haw ha hw)
    (tateSpecialize_uA a w q hq haw ha hw) (tateSpecialize_uW a w q hq haw ha hw)
    isUnit_one_sub_uA isUnit_one_sub_uW isUnit_one_sub_uQA isUnit_one_sub_uQW H

/-! ## ★出典の紐付け(`.src`) -/

def TateUniv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——万有な環)",
    sectionId := "genell-def-3-3" }

def tateSpecialize.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——万有な環からの特殊化)",
    sectionId := "genell-def-3-3" }

def tateDefect_eq_zero_of_univ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——葉 (b) の万有な環への帰着)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
