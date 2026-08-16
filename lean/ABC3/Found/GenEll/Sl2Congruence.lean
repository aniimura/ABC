import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.SpecialLinearGroup
import Mathlib.Tactic.LinearCombination
import Mathlib.Data.Nat.Choose.Dvd
import Mathlib.Data.ZMod.Basic
import Mathlib.LinearAlgebra.Matrix.Trace
import ABC3.Meta.Claim

/-!
# 合同核の代数 —— [GenEll] Lemma 3.1, (iv) の部品(段 1)

★**これは原典の項目の実装ではない。** `Lemma 3.1, (iv)` の道具であり、`.src` は付けない。
位置づけは `ResearchPaper/genell-goal.md` の段 1・柱 2 の「段 1」。

## ★何を示すか

(iv) の標準的な証明は合同フィルトレーション
`K_n = ker(SL₂(ℤ_l) → SL₂(ℤ/l^n))` に沿う帰納法で、その各段で使うのが

> **`K_n / K_{n+1}` は `𝔰𝔩₂(𝔽_l)`(跡 0 の行列のなす**加法**群)である**

という事実である。その**代数の核**は「`p² = 0` なる `p` に対して
`1 + p·A` の形の行列が**加法的に**掛かる」ことであり、それは環の計算だけで出る。
本ファイルはそこを取る。

★段の全体は `Sl2Adjoint.lean` の docstring の表を参照。本ファイルは段 1 の**代数部分**で、
`ℤ_l` の位相(逆極限としての `K_n` の定義)は含まない。
-/

namespace ABC3.Found.GenEll

open Matrix

section Congruence

variable {R : Type*} [CommRing R] {p : R}

/-- ★**合同核が加法群になることの核**。

`p² = 0` なら `(1 + p·A)(1 + p·B) = 1 + p·(A + B)`。
★これが「`K_n / K_{n+1}` が**加法**群である」ことの中身である——
積が和に化けるのは `p²` の項が消えるからで、それ以外の理由は無い。 -/
theorem mul_one_add_smul (hp : p * p = 0) (A B : Matrix (Fin 2) (Fin 2) R) :
    (1 + p • A) * (1 + p • B) = 1 + p • (A + B) := by
  have hAB : (p • A) * (p • B) = (p * p) • (A * B) := by
    rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  rw [add_mul, mul_add, mul_add, one_mul, mul_one, hAB, hp, zero_smul, add_zero,
    smul_add, add_assoc, one_mul, add_comm (p • B) (p • A)]

/-- `1 + p·A` の逆元は `1 + p·(−A)`。 -/
theorem one_add_smul_mul_neg (hp : p * p = 0) (A : Matrix (Fin 2) (Fin 2) R) :
    (1 + p • A) * (1 + p • (-A)) = 1 := by
  rw [mul_one_add_smul hp, add_neg_cancel, smul_zero, add_zero]

/-- ★**行列式が跡を測る**。

2×2 では `det(1 + p·A) = 1 + p·tr(A)`(`p² = 0` のとき)。
★したがって **`1 + p·A ∈ SL₂` ⟺ `p·tr(A) = 0`** であり、
`p` が `𝔽_l` 上の単元倍として効く状況では **`tr(A) = 0`**——
すなわち合同核の元は**跡 0** の `A` で書ける。これが `𝔰𝔩₂` が現れる理由である。 -/
theorem det_one_add_smul (hp : p * p = 0) (A : Matrix (Fin 2) (Fin 2) R) :
    (1 + p • A).det = 1 + p * (A 0 0 + A 1 1) := by
  have hform : (1 : Matrix (Fin 2) (Fin 2) R) + p • A
      = !![1 + p * A 0 0, p * A 0 1; p * A 1 0, 1 + p * A 1 1] := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp
  rw [hform, Matrix.det_fin_two_of]
  linear_combination (A 0 0 * A 1 1 - A 0 1 * A 1 0) * hp

/-- ★`p·tr(A) = 0` なら `1 + p·A` は `SL(2, R)` の元。 -/
def slOfSmul (hp : p * p = 0) (A : Matrix (Fin 2) (Fin 2) R)
    (htr : p * (A 0 0 + A 1 1) = 0) : Matrix.SpecialLinearGroup (Fin 2) R :=
  ⟨1 + p • A, by rw [det_one_add_smul hp, htr, add_zero]⟩

@[simp] theorem coe_slOfSmul (hp : p * p = 0) (A : Matrix (Fin 2) (Fin 2) R)
    (htr : p * (A 0 0 + A 1 1) = 0) :
    ((slOfSmul hp A htr : Matrix.SpecialLinearGroup (Fin 2) R)
      : Matrix (Fin 2) (Fin 2) R) = 1 + p • A := rfl

/-- ★★**合同核は加法群である**——`SL` の積が `𝔰𝔩₂` の和に対応する。

これが段 1 の主張「`K_n / K_{n+1} ≅ 𝔰𝔩₂`」の**群構造の部分**である。 -/
theorem slOfSmul_mul (hp : p * p = 0) (A B : Matrix (Fin 2) (Fin 2) R)
    (hA : p * (A 0 0 + A 1 1) = 0) (hB : p * (B 0 0 + B 1 1) = 0)
    (hAB : p * ((A + B) 0 0 + (A + B) 1 1) = 0) :
    slOfSmul hp A hA * slOfSmul hp B hB = slOfSmul hp (A + B) hAB := by
  apply Subtype.ext
  rw [Matrix.SpecialLinearGroup.coe_mul, coe_slOfSmul, coe_slOfSmul, coe_slOfSmul]
  exact mul_one_add_smul hp A B

end Congruence

section WhereFiveEnters

/-! ## ★★`l ≥ 5` が効く正確な場所(段 4 の算術核)

(iv) の証明で `x` を `upper 1` の**持ち上げ**とすると、`x^l` の
`K_1 / K_2` における類は

```
E₁₂ + (Σ_{i<l} Ad(u^i)) · A     (u = upper 1 mod l、A は持ち上げ方の差)
```

になる。`Ad(u) = 1 + D` で `D` は `𝔰𝔩₂`(3 次元)上冪零、`D³ = 0`。
**hockey-stick 恒等式** `Σ_{i<l} C(i,j) = C(l, j+1)` により

```
Σ_{i<l} (1 + D)^i = C(l,1) + C(l,2)·D + C(l,3)·D²
```

★したがって **`C(l,1) = C(l,2) = C(l,3) = 0 (mod l)`** なら第 2 項が消え、
**`x^l` の類は持ち上げ方に依らず `E₁₂`(≠ 0)** になる。それが段 4 の要である。

★**そしてそこに `l ≥ 5` が要る**——`C(l,3)` を落とすには `3 < l` が必要で、
`l = 3` では `C(3,3) = 1 ≢ 0` である(下の負の対照)。 -/

/-- ★`l ≥ 5` が素数なら `C(l,1) = C(l,2) = C(l,3) = 0` in `𝔽_l`。 -/
theorem choose_one_two_three_eq_zero (l : ℕ) (hl : Nat.Prime l) (h5 : 5 ≤ l) :
    ((l.choose 1 : ℕ) : ZMod l) = 0 ∧ ((l.choose 2 : ℕ) : ZMod l) = 0
      ∧ ((l.choose 3 : ℕ) : ZMod l) = 0 := by
  refine ⟨(ZMod.natCast_eq_zero_iff _ _).2 ?_, (ZMod.natCast_eq_zero_iff _ _).2 ?_,
    (ZMod.natCast_eq_zero_iff _ _).2 ?_⟩
  · exact Nat.Prime.dvd_choose_self hl one_ne_zero (by omega)
  · exact Nat.Prime.dvd_choose_self hl two_ne_zero (by omega)
  · exact Nat.Prime.dvd_choose_self hl (by norm_num) (by omega)

/-- ★★**負の対照**: `l = 3` では `C(3,3) = 1 ≠ 0` なので、段 4 の消去が**破れる**。

すなわち `l ≥ 5` は飾りではない。 -/
theorem choose_three_three_ne_zero_mod_three :
    ((Nat.choose 3 3 : ℕ) : ZMod 3) ≠ 0 := by decide

/-- **hockey-stick 恒等式**の `range` 版。

`∑_{i<n} C(i, j) = C(n, j+1)`。mathlib の `Nat.sum_Icc_choose` は `Icc` 版なので、
`range` 版を帰納法で直に取る(`Icc` への変換より短い)。 -/
theorem sum_range_choose_eq (n j : ℕ) :
    ∑ i ∈ Finset.range n, i.choose j = n.choose (j + 1) := by
  induction n with
  | zero => simp
  | succ k ih =>
    rw [Finset.sum_range_succ, ih, Nat.choose_succ_succ' k j, Nat.add_comm]

/-- ★冪単元の冪展開。`D³ = 0` なら `(1 + D)^i = 1 + i·D + C(i,2)·D²`。 -/
theorem one_add_pow_of_cube_eq_zero {R : Type*} [CommRing R] {D : R} (hD : D ^ 3 = 0) (i : ℕ) :
    (1 + D) ^ i = 1 + (i : R) * D + (i.choose 2 : R) * D ^ 2 := by
  induction i with
  | zero => simp
  | succ k ih =>
    rw [pow_succ, ih, Nat.choose_succ_succ' k 1, Nat.choose_one_right]
    push_cast
    have h3 : D ^ 3 = 0 := hD
    linear_combination ((k.choose 2 : R)) * h3

/-- ★★**段 4 の代数的な核**。

`l` が `5` 以上の素数で、環の標数が `l`、`D³ = 0` なら

```
∑_{i<l} (1 + D)^i = 0
```

★これが「`x^l` の `K_1/K_2` における類が**持ち上げ方に依らない**」ことの中身である。
3 つの係数 `C(l,1)`・`C(l,2)`・`C(l,3)` がすべて `𝔽_l` で消えることを使う——
**`C(l,3)` を消すのに `l ≥ 5` が要る**(`choose_three_three_ne_zero_mod_three` を参照)。 -/
theorem sum_range_one_add_pow_eq_zero {R : Type*} [CommRing R] (l : ℕ)
    (hl : Nat.Prime l) (h5 : 5 ≤ l) (hchar : (l : R) = 0) {D : R} (hD : D ^ 3 = 0) :
    ∑ i ∈ Finset.range l, (1 + D) ^ i = 0 := by
  have hexp : ∀ i ∈ Finset.range l,
      (1 + D) ^ i = 1 + (i : R) * D + (i.choose 2 : R) * D ^ 2 :=
    fun i _ => one_add_pow_of_cube_eq_zero hD i
  rw [Finset.sum_congr rfl hexp]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.sum_mul]
  -- 3 つの係数を `C(l,1)`・`C(l,2)`・`C(l,3)` に直す
  have c1 : ∑ _i ∈ Finset.range l, (1 : R) = (l.choose 1 : ℕ) := by
    simp [Nat.choose_one_right]
  have c2 : ∑ i ∈ Finset.range l, (i : R) = ((l.choose 2 : ℕ) : R) := by
    have := sum_range_choose_eq l 1
    simp only [Nat.choose_one_right] at this
    calc ∑ i ∈ Finset.range l, (i : R)
        = ((∑ i ∈ Finset.range l, i : ℕ) : R) := by push_cast; ring
      _ = ((l.choose 2 : ℕ) : R) := by rw [this]
  have c3 : ∑ i ∈ Finset.range l, ((i.choose 2 : ℕ) : R) = ((l.choose 3 : ℕ) : R) := by
    calc ∑ i ∈ Finset.range l, ((i.choose 2 : ℕ) : R)
        = ((∑ i ∈ Finset.range l, i.choose 2 : ℕ) : R) := by push_cast; ring
      _ = ((l.choose 3 : ℕ) : R) := by rw [sum_range_choose_eq l 2]
  rw [c1, c2, c3]
  -- 3 つとも `𝔽_l` で消える(= `l ∣ C(l,k)`)ので、標数 `l` の環でも消える
  have d1 : ((l.choose 1 : ℕ) : R) = 0 := by
    obtain ⟨m, hm⟩ := Nat.Prime.dvd_choose_self hl one_ne_zero (by omega)
    rw [hm]; push_cast; rw [hchar, zero_mul]
  have d2 : ((l.choose 2 : ℕ) : R) = 0 := by
    obtain ⟨m, hm⟩ := Nat.Prime.dvd_choose_self hl two_ne_zero (by omega)
    rw [hm]; push_cast; rw [hchar, zero_mul]
  have d3 : ((l.choose 3 : ℕ) : R) = 0 := by
    obtain ⟨m, hm⟩ := Nat.Prime.dvd_choose_self (k := 3) hl (by norm_num) (by omega)
    rw [hm]; push_cast; rw [hchar, zero_mul]
  rw [d1, d2, d3, zero_mul, zero_mul, add_zero, add_zero]

end WhereFiveEnters

section FiniteLevel

/-! ## ★合同フィルトレーションを**有限段**に降ろす(位相を使わない部分)

★(iv) の証明は `SL₂(ℤ_l)` の閉部分群についての主張だが、**位相が要るのは最後の 1 歩だけ**である:

| | 内容 | 位相 |
|---|---|---|
| (A) | `H ≤ SL₂(ℤ/l^n)` が `SL₂(𝔽_l)` へ全射なら `H = ⊤` | **要らない**(有限群論) |
| (B) | 閉部分群 `J ⊆ SL₂(ℤ_l)` は各 `mod l^n` の像で決まる | 要る(逆極限) |

(A) の帰納法の各段で使うのが `mul_one_add_smul` などであり、そこで要る仮定 `p² = 0` は
**`p = l^n` を `ℤ/l^{n+1}` で見たときに実際に成り立つ**。それを取る。 -/

/-- ★`ℤ/l^{n+1}` の中で `l^n` は**平方が消える**(`n ≥ 1`)。

これが `mul_one_add_smul` / `det_one_add_smul` の仮定 `p * p = 0` を、
実際の合同フィルトレーション `K_n ⊆ SL₂(ℤ/l^{n+1})` について満たす。 -/
theorem pow_self_mul_self_eq_zero (l n : ℕ) (hn : 1 ≤ n) :
    ((l : ZMod (l ^ (n + 1))) ^ n) * ((l : ZMod (l ^ (n + 1))) ^ n) = 0 := by
  rw [← pow_add]
  have hcast : (l : ZMod (l ^ (n + 1))) ^ (n + n) = ((l ^ (n + n) : ℕ) : ZMod (l ^ (n + 1))) := by
    push_cast
    ring
  rw [hcast, ZMod.natCast_eq_zero_iff]
  exact pow_dvd_pow l (by omega)

/-- ★上と同じことを、`mul_one_add_smul` に直結する形で述べたもの。

`ℤ/l^{n+1}` 上の 2×2 行列について、`1 + l^n·A` の形の元は**加法的に**掛かる。 -/
theorem congruence_mul (l n : ℕ) (hn : 1 ≤ n)
    (A B : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) :
    (1 + ((l : ZMod (l ^ (n + 1))) ^ n) • A) * (1 + ((l : ZMod (l ^ (n + 1))) ^ n) • B)
      = 1 + ((l : ZMod (l ^ (n + 1))) ^ n) • (A + B) :=
  mul_one_add_smul (pow_self_mul_self_eq_zero l n hn) A B

/-- ★同じく、行列式が跡を測ることの合同フィルトレーション版。 -/
theorem congruence_det (l n : ℕ) (hn : 1 ≤ n)
    (A : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) :
    (1 + ((l : ZMod (l ^ (n + 1))) ^ n) • A).det
      = 1 + ((l : ZMod (l ^ (n + 1))) ^ n) * (A 0 0 + A 1 1) :=
  det_one_add_smul (pow_self_mul_self_eq_zero l n hn) A

end FiniteLevel

section ConjugationIsAdjoint

variable {R : Type*} [CommRing R] {p : R}

/-- ★★**共役作用は随伴作用に落ちる**——段 1b の中身。

`g · (1 + p·A) · g⁻¹ = 1 + p·(g A g⁻¹)`。

★これが「合同核 `K_n/K_{n+1}` への共役作用が `𝔰𝔩₂` の**随伴表現**である」ことの正体である。
★**`p² = 0` すら要らない**——効いているのは「スカラー `p` が行列と可換」ということだけ。
位相も逆極限も要らない。 -/
theorem conj_one_add_smul (g : Matrix.SpecialLinearGroup (Fin 2) R)
    (A : Matrix (Fin 2) (Fin 2) R) :
    (g : Matrix (Fin 2) (Fin 2) R) * (1 + p • A) * ((g⁻¹ : Matrix.SpecialLinearGroup (Fin 2) R)
        : Matrix (Fin 2) (Fin 2) R)
      = 1 + p • ((g : Matrix (Fin 2) (Fin 2) R) * A
          * ((g⁻¹ : Matrix.SpecialLinearGroup (Fin 2) R) : Matrix (Fin 2) (Fin 2) R)) := by
  have hgg : (g : Matrix (Fin 2) (Fin 2) R)
      * ((g⁻¹ : Matrix.SpecialLinearGroup (Fin 2) R) : Matrix (Fin 2) (Fin 2) R) = 1 := by
    rw [← Matrix.SpecialLinearGroup.coe_mul, mul_inv_cancel]
    rfl
  rw [mul_add, mul_one, add_mul, hgg, Matrix.mul_smul, Matrix.smul_mul, mul_assoc]

/-- ★随伴作用は**跡を保つ**——したがって `𝔰𝔩₂`(跡 0)は不変部分空間である。

これが「共役作用が `𝔰𝔩₂` の上の表現として閉じる」ことの中身。 -/
theorem trace_conj (g : Matrix.SpecialLinearGroup (Fin 2) R)
    (A : Matrix (Fin 2) (Fin 2) R) :
    Matrix.trace ((g : Matrix (Fin 2) (Fin 2) R) * A
        * ((g⁻¹ : Matrix.SpecialLinearGroup (Fin 2) R) : Matrix (Fin 2) (Fin 2) R))
      = Matrix.trace A := by
  have hgg : ((g⁻¹ : Matrix.SpecialLinearGroup (Fin 2) R) : Matrix (Fin 2) (Fin 2) R)
      * (g : Matrix (Fin 2) (Fin 2) R) = 1 := by
    rw [← Matrix.SpecialLinearGroup.coe_mul, inv_mul_cancel]
    rfl
  rw [Matrix.trace_mul_comm, ← mul_assoc, hgg, one_mul]

end ConjugationIsAdjoint

section KernelShape

/-! ## ★合同核の元が `1 + l^n·A` の形に**書けること**

段 (A) の帰納法では「`K_n` の元は `1 + l^n·A` と書ける」を使う。
その根は `ZMod` の側の事実——**`ℤ/l^{n+1} → ℤ/l^n` の核は `(l^n)`** ——である。
★これも位相を要しない。 -/

/-- ★`ℤ/l^{n+1} → ℤ/l^n` で `0` に行く元は `l^n` の倍数である。 -/
theorem exists_mul_of_castHom_eq_zero (l n : ℕ) (hl : 0 < l) (x : ZMod (l ^ (n + 1)))
    (hx : ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n)) x = 0) :
    ∃ y : ZMod (l ^ (n + 1)), x = ((l : ZMod (l ^ (n + 1))) ^ n) * y := by
  haveI : NeZero (l ^ (n + 1)) := ⟨pow_ne_zero _ (by omega)⟩
  -- `castHom` を `x.val` の像として読む
  have hcast : ((x.val : ℕ) : ZMod (l ^ n)) = 0 := by
    rw [ZMod.natCast_val, ← ZMod.castHom_apply (h := pow_dvd_pow l (Nat.le_succ n))]
    exact hx
  obtain ⟨m, hm⟩ := (ZMod.natCast_eq_zero_iff _ _).1 hcast
  refine ⟨(m : ZMod (l ^ (n + 1))), ?_⟩
  have hx' : ((x.val : ℕ) : ZMod (l ^ (n + 1))) = x := by
    rw [ZMod.natCast_val, ZMod.cast_id]
  rw [← hx', hm]
  push_cast
  ring

/-- ★★**合同核の元は `1 + l^n·A` の形に書ける**(`ℤ/l^{n+1}` の中で)。

`ℤ/l^{n+1} → ℤ/l^n` の核に入る行列は、成分ごとに `l^n` の倍数だからである。
★これと `congruence_mul` / `congruence_det` / `conj_one_add_smul` が揃って、
段 (A) の帰納法の 1 段ぶんの材料になる。 -/
theorem exists_matrix_of_castHom_eq_zero (l n : ℕ) (hl : 0 < l)
    (M : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))))
    (hM : ∀ i j, ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n))
      (M i j - (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j) = 0) :
    ∃ A : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))),
      M = 1 + ((l : ZMod (l ^ (n + 1))) ^ n) • A := by
  choose A hA using fun i j => exists_mul_of_castHom_eq_zero l n hl
    (M i j - (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j) (hM i j)
  refine ⟨Matrix.of A, ?_⟩
  ext i j
  have := hA i j
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, Matrix.of_apply]
  linear_combination this

end KernelShape

end ABC3.Found.GenEll
