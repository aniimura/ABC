import ABC3.Found.GenEll.Sl2Congruence
import ABC3.Found.GenEll.Sl2Adjoint
import ABC3.Found.GenEll.Lemma31
import Mathlib.Data.Nat.Multiplicity
import Mathlib.Algebra.Module.ZMod

/-!
# [GenEll] Lemma 3.1, (iv) の**有限段**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.14。**目視確認済み**(`ResearchPaper/1_Structured/…/lemma-3-1.html`)。

原文 (GenEll p.14):
> (iv) Let J ⊆ GL2(Zl) be a closed subgroup whose image HJ in GL2(Fl) contains

## ★★位相が要るのは最後の 1 歩だけ —— その手前を全部取る

`Sl2Congruence.lean` の切り分け(段 (A)/(B))のうち、本ファイルは **(A) 全体**である:

| | 内容 | 位相 |
|---|---|---|
| **(A)** | `H ≤ SL₂(ℤ/l^n)` が `SL₂(𝔽_l)` へ全射なら `H = ⊤` | ★**要らない**(本ファイル) |
| (B) | 閉部分群 `J ⊆ SL₂(ℤ_l)` は各 `mod l^n` の像で決まる | 要る(逆極限) |

★**(A) は有限群論だけで閉じる。** 原文は [Serre] Chapter IV, §3.4, Lemma 3 を引くが、
その [Serre] は `0_Source` に無いので**自分で証明する**——それが本ファイルである。

## ★証明の骨格(`n → n+1` の帰納法)

`K ≝ ker(SL₂(ℤ/l^{n+1}) → SL₂(ℤ/l^n))` と置く。帰納法の仮定から `H·K = ⊤`。

`H ∩ K` の像 `M ⊆ 𝔰𝔩₂(𝔽_l)` は `𝔽_l`-部分空間で、`H` が `SL₂(𝔽_l)` へ全射なので
**随伴作用で不変**。ゆえに `sl2_adjoint_irreducible`(`2 ≠ 0` だけで足りる)により

- `M = 𝔰𝔩₂` なら `K ≤ H`、`H·K = ⊤` と併せて `H = ⊤`。
- `M = 0` なら `H ∩ K = 1`。

★★**場合分けは要らなかった。** 当初は `M = 0` を背理法で潰すつもりだったが、
**`M ≠ 0` を直接示せる**——`slE ∈ M` が always 成り立つ:

`H·K = ⊤` から `π(h) = π(u)`(`u ≝ upper 1`)なる `h ∈ H` を取り、`h = k·u`(`k ∈ K`)と書くと

```
h^{l^n} = (1 + l^n·Σ_{i<l^n} Ad(u^i)(B)) · u^{l^n}
```

で、**`Σ_{i<l^n} Ad(u^i)` は `𝔽_l` 上ゼロ写像**(`sum_adU_eq_zero`)。
したがって **`h^{l^n} = u^{l^n} = upper(l^n) = 1 + l^n·e`** であり、`e ∈ M`。
★これで `M ≠ 0` が出るので、あとは既約性から `M = 𝔰𝔩₂` が直ちに従う。

★**そしてそこに `l ≥ 5` が要る**——`C(l^n,3)` を落とすのに `3 < l^n` が要り、
`l = 3, n = 1` では `C(3,3) = 1 ≠ 0` である
(`Sl2Congruence.lean` の `choose_three_three_ne_zero_mod_three`)。

## ★★段 (A) は完成した(2026-08-17)

本体は **`subgroup_eq_top_of_redPow_surj`**:

> `l ≥ 5` 素数、`H ≤ SL₂(ℤ/l^{k+1})` が `SL₂(𝔽_l)` へ全射なら `H = ⊤`。

★**`sorry` 0・公理 0。** mathlib は `SL₂` の合同部分群の理論を持たないので、
`ZMod` の水準から組み上げた。

★残るのは **段 (B)** ——`ℤ_l` の位相で「閉部分群は各 `mod l^n` の像で決まる」。
そこだけが位相を要する。

★**`.src` は付けない。** `Lemma 3.1, (iv)` はまだ完成していない
——`.src` の 2 値規則(条なし = 命題全体が済んだ)を守る。
-/

namespace ABC3.Found.GenEll

open Matrix
open scoped MatrixGroups

/-! ## ★段 1 —— 冪和の `l^n` 版

`Sl2Congruence.lean` の `sum_range_one_add_pow_eq_zero` は項数 `l` の版だった。
帰納法の第 `n` 段では**項数 `l^n`** が要る。 -/

/-- ★**項数が `l^n` でも冪和は消える**(`l ≥ 5` 素数、`n ≥ 1`、標数 `l`、`D³ = 0`)。

`Σ_{i<m} (1+D)^i = C(m,1) + C(m,2)·D + C(m,3)·D²`(hockey-stick)であり、
`m = l^n` のとき 3 つの二項係数がすべて `l` で割れる
(`Nat.Prime.dvd_choose_pow`——`k ≠ 0` かつ `k ≠ l^n` が条件で、
`k ≤ 3 < 5 ≤ l ≤ l^n` から後者が出る)。

★これが「`x^{l^n}` の類が**持ち上げ方に依らない**」ことの中身である。 -/
theorem sum_range_one_add_pow_pow_eq_zero {R : Type*} [CommRing R] (l n : ℕ)
    (hl : Nat.Prime l) (h5 : 5 ≤ l) (hn : 1 ≤ n) (hchar : (l : R) = 0)
    {D : R} (hD : D ^ 3 = 0) :
    ∑ i ∈ Finset.range (l ^ n), (1 + D) ^ i = 0 := by
  have hlm : l ≤ l ^ n := Nat.le_self_pow (by omega) l
  -- 3 つの二項係数が `l` で割れる
  have hdvd : ∀ k : ℕ, k ≠ 0 → k ≤ 3 → l ∣ (l ^ n).choose k := by
    intro k hk0 hk3
    exact hl.dvd_choose_pow hk0 (by omega)
  have hzero : ∀ k : ℕ, k ≠ 0 → k ≤ 3 → (((l ^ n).choose k : ℕ) : R) = 0 := by
    intro k hk0 hk3
    obtain ⟨c, hc⟩ := hdvd k hk0 hk3
    rw [hc]; push_cast; rw [hchar, zero_mul]
  -- hockey-stick で 3 項に落とす
  have hexp : ∀ i ∈ Finset.range (l ^ n),
      (1 + D) ^ i = 1 + (i : R) * D + (i.choose 2 : R) * D ^ 2 :=
    fun i _ => one_add_pow_of_cube_eq_zero hD i
  rw [Finset.sum_congr rfl hexp, Finset.sum_add_distrib, Finset.sum_add_distrib,
    ← Finset.sum_mul, ← Finset.sum_mul]
  have c1 : ∑ _i ∈ Finset.range (l ^ n), (1 : R) = (((l ^ n).choose 1 : ℕ) : R) := by
    simp [Nat.choose_one_right]
  have c2 : ∑ i ∈ Finset.range (l ^ n), (i : R) = (((l ^ n).choose 2 : ℕ) : R) := by
    have hs := sum_range_choose_eq (l ^ n) 1
    simp only [Nat.choose_one_right] at hs
    calc ∑ i ∈ Finset.range (l ^ n), (i : R)
        = ((∑ i ∈ Finset.range (l ^ n), i : ℕ) : R) := by push_cast; ring
      _ = (((l ^ n).choose 2 : ℕ) : R) := by rw [hs]
  have c3 : ∑ i ∈ Finset.range (l ^ n), ((i.choose 2 : ℕ) : R)
      = (((l ^ n).choose 3 : ℕ) : R) := by
    calc ∑ i ∈ Finset.range (l ^ n), ((i.choose 2 : ℕ) : R)
        = ((∑ i ∈ Finset.range (l ^ n), i.choose 2 : ℕ) : R) := by push_cast; ring
      _ = (((l ^ n).choose 3 : ℕ) : R) := by rw [sum_range_choose_eq (l ^ n) 2]
  rw [c1, c2, c3, hzero 1 one_ne_zero (by omega), hzero 2 two_ne_zero (by omega),
    hzero 3 (by norm_num) (by omega), zero_mul, zero_mul, add_zero, add_zero]

/-! ## ★段 2 —— 段の間の還元射 -/

variable (l : ℕ) [Fact (Nat.Prime l)]

/-- `SL₂(ℤ/l^b) → SL₂(ℤ/l^a)`(`a ≤ b`)。 -/
def redLevel {a b : ℕ} (h : a ≤ b) :
    SL(2, ZMod (l ^ b)) →* SL(2, ZMod (l ^ a)) :=
  Matrix.SpecialLinearGroup.map (ZMod.castHom (pow_dvd_pow l h) (ZMod (l ^ a)))

@[simp] theorem redLevel_coe {a b : ℕ} (h : a ≤ b) (g : SL(2, ZMod (l ^ b))) (i j : Fin 2) :
    ((redLevel l h g : SL(2, ZMod (l ^ a))) : Matrix (Fin 2) (Fin 2) (ZMod (l ^ a))) i j
      = ZMod.castHom (pow_dvd_pow l h) (ZMod (l ^ a)) ((g : Matrix (Fin 2) (Fin 2) _) i j) :=
  rfl

/-- ★還元射は合成で閉じる。 -/
theorem redPow_comp_redLevel {n : ℕ} (hn : n ≠ 0) :
    (redPow l n hn).comp (redLevel l (Nat.le_succ n)) = redPow l (n + 1) (Nat.succ_ne_zero n) := by
  ext g i j
  show ZMod.castHom (dvd_pow_self l hn) (ZMod l)
      (ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n))
        ((g : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j))
    = ZMod.castHom (dvd_pow_self l (Nat.succ_ne_zero n)) (ZMod l) _
  rw [← RingHom.comp_apply, ZMod.castHom_comp]

/-- ★`ZMod (l^1) → ZMod l` は**単射**である。

★`l^1 = l` は自然数の等式なので、`x.val < l^1` を `x.val < l` に読み替えるだけでよい
——**型の依存性に触れずに済む**。 -/
theorem castHom_pow_one_injective (l : ℕ) [NeZero l] :
    Function.Injective (ZMod.castHom (dvd_pow_self l one_ne_zero) (ZMod l)) := by
  haveI : NeZero (l ^ 1) := ⟨by simpa using (NeZero.ne l)⟩
  refine (injective_iff_map_eq_zero _).2 fun x hx => ?_
  have hval : ((x.val : ℕ) : ZMod l) = 0 := by
    rw [ZMod.natCast_val, ← ZMod.castHom_apply (h := dvd_pow_self l one_ne_zero)]
    exact hx
  obtain ⟨m, hm⟩ := (ZMod.natCast_eq_zero_iff _ _).1 hval
  have hlt : x.val < l := by
    have := ZMod.val_lt x
    simpa using this
  have : x.val = 0 := by
    rcases Nat.eq_zero_or_pos m with rfl | hm0
    · simpa using hm
    · exfalso; nlinarith [hm, hm0, hlt]
  exact (ZMod.val_eq_zero x).1 this

/-! ## ★段 3 —— `mod l` の核と `l^n` 倍の関係

★合同核 `K` の元 `1 + l^n·A` について、**`A` は `mod l` でしか効かない**。
その 2 方向を先に取る。 -/

section ModPrime

variable {l : ℕ}

/-- ★`ZMod (l^m) → ZMod l` の核は `l` の倍数からなる。 -/
theorem exists_mul_of_castPrime_eq_zero {m : ℕ} (hm : m ≠ 0) (hl : 0 < l)
    (x : ZMod (l ^ m))
    (hx : ZMod.castHom (dvd_pow_self l hm) (ZMod l) x = 0) :
    ∃ y : ZMod (l ^ m), x = (l : ZMod (l ^ m)) * y := by
  haveI : NeZero (l ^ m) := ⟨pow_ne_zero _ (by omega)⟩
  have hcast : ((x.val : ℕ) : ZMod l) = 0 := by
    rw [ZMod.natCast_val, ← ZMod.castHom_apply (h := dvd_pow_self l hm)]
    exact hx
  obtain ⟨c, hc⟩ := (ZMod.natCast_eq_zero_iff _ _).1 hcast
  refine ⟨(c : ZMod (l ^ m)), ?_⟩
  have hx' : ((x.val : ℕ) : ZMod (l ^ m)) = x := by rw [ZMod.natCast_val, ZMod.cast_id]
  rw [← hx', hc]
  push_cast
  ring

/-- ★逆向き: `l^n·x = 0`(`ℤ/l^{n+1}` の中)なら `x ≡ 0 (mod l)`。 -/
theorem castPrime_eq_zero_of_pow_mul_eq_zero (n : ℕ) (hl : 0 < l) (hn : 1 ≤ n)
    (x : ZMod (l ^ (n + 1)))
    (hx : ((l : ZMod (l ^ (n + 1)))) ^ n * x = 0) :
    ZMod.castHom (dvd_pow_self l (Nat.succ_ne_zero n)) (ZMod l) x = 0 := by
  haveI : NeZero (l ^ (n + 1)) := ⟨pow_ne_zero _ (by omega)⟩
  have hx' : ((x.val : ℕ) : ZMod (l ^ (n + 1))) = x := by rw [ZMod.natCast_val, ZMod.cast_id]
  have hnat : (((l ^ n * x.val : ℕ)) : ZMod (l ^ (n + 1))) = 0 := by
    push_cast
    rw [hx']
    exact hx
  obtain ⟨c, hc⟩ := (ZMod.natCast_eq_zero_iff _ _).1 hnat
  -- `l^{n+1} ∣ l^n · x.val` から `l ∣ x.val`
  have hdvd : l ∣ x.val := by
    have h1 : l ^ n * x.val = l ^ n * (l * c) := by rw [hc]; ring
    have h2 : x.val = l * c := Nat.eq_of_mul_eq_mul_left (pow_pos hl n) h1
    exact ⟨c, h2⟩
  rw [← hx', map_natCast]
  exact (ZMod.natCast_eq_zero_iff _ _).2 hdvd

/-- ★★**`l^n·(−)` は `mod l` の情報しか見ない。**

これが「`K ≅ 𝔰𝔩₂(𝔽_l)`」の well-defined 性の中身である。 -/
theorem pow_mul_congr_of_castPrime_eq (n : ℕ) (hl : 0 < l) (x y : ZMod (l ^ (n + 1)))
    (h : ZMod.castHom (dvd_pow_self l (Nat.succ_ne_zero n)) (ZMod l) x
       = ZMod.castHom (dvd_pow_self l (Nat.succ_ne_zero n)) (ZMod l) y) :
    ((l : ZMod (l ^ (n + 1)))) ^ n * x = ((l : ZMod (l ^ (n + 1)))) ^ n * y := by
  have hsub : ZMod.castHom (dvd_pow_self l (Nat.succ_ne_zero n)) (ZMod l) (x - y) = 0 := by
    rw [map_sub, h, sub_self]
  obtain ⟨z, hz⟩ := exists_mul_of_castPrime_eq_zero (Nat.succ_ne_zero n) hl (x - y) hsub
  have hzero : ((l : ZMod (l ^ (n + 1)))) ^ (n + 1) = 0 := by
    have : ((l : ZMod (l ^ (n + 1)))) ^ (n + 1) = ((l ^ (n + 1) : ℕ) : ZMod (l ^ (n + 1))) := by
      push_cast; ring
    rw [this, ZMod.natCast_self]
  have : ((l : ZMod (l ^ (n + 1)))) ^ n * (x - y) = 0 := by
    rw [hz, ← mul_assoc, ← pow_succ, hzero, zero_mul]
  linear_combination this

end ModPrime

/-! ## ★段 4 —— 合同核の `𝔰𝔩₂(𝔽_l)` における像

★ここが「`K ≅ 𝔰𝔩₂(𝔽_l)`」を**加法群として**取る場所である。
`Sl2Congruence.lean` の `congruence_mul`(積が和になる)がそのまま効く。 -/

section KernelModule

open Matrix

variable (l : ℕ) [Fact (Nat.Prime l)] (n : ℕ)

/-- `ZMod (l^{n+1}) → ZMod l`。 -/
def toPrime : ZMod (l ^ (n + 1)) →+* ZMod l :=
  ZMod.castHom (dvd_pow_self l (Nat.succ_ne_zero n)) (ZMod l)

/-- 行列の `mod l` 還元(環準同型)。 -/
def redMatHom : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))) →+*
    Matrix (Fin 2) (Fin 2) (ZMod l) :=
  (toPrime l n).mapMatrix

@[simp] theorem redMatHom_apply (B : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) (i j : Fin 2) :
    redMatHom l n B i j = toPrime l n (B i j) := rfl

/-- 原文の `1 + l^n·A`(合同核の元の形)。 -/
def congElt (A : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) :
    Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))) :=
  1 + ((l : ZMod (l ^ (n + 1))) ^ n) • A

theorem congElt_mul (hn : 1 ≤ n) (A B : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) :
    congElt l n A * congElt l n B = congElt l n (A + B) :=
  congruence_mul l n hn A B

@[simp] theorem congElt_zero : congElt l n 0 = 1 := by simp [congElt]

/-- ★**`H` の合同核部分が `𝔰𝔩₂(𝔽_l)` に落とす像**。加法部分群になる。 -/
def kerImage (H : Subgroup SL(2, ZMod (l ^ (n + 1)))) (hn : 1 ≤ n) :
    AddSubgroup (Matrix (Fin 2) (Fin 2) (ZMod l)) where
  carrier := {A | ∃ h ∈ H, ∃ B, (h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))))
      = congElt l n B ∧ redMatHom l n B = A}
  zero_mem' := ⟨1, H.one_mem, 0, by simp, by simp⟩
  add_mem' := by
    rintro A₁ A₂ ⟨h₁, hh₁, B₁, hB₁, rfl⟩ ⟨h₂, hh₂, B₂, hB₂, rfl⟩
    refine ⟨h₁ * h₂, H.mul_mem hh₁ hh₂, B₁ + B₂, ?_, by rw [map_add]⟩
    rw [Matrix.SpecialLinearGroup.coe_mul, hB₁, hB₂, congElt_mul l n hn]
  neg_mem' := by
    rintro A ⟨h, hh, B, hB, rfl⟩
    refine ⟨h⁻¹, H.inv_mem hh, -B, ?_, by rw [map_neg]⟩
    have hone : (h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) * congElt l n (-B) = 1 := by
      rw [hB, congElt_mul l n hn]
      simp
    have hinv : (h⁻¹ : SL(2, ZMod (l ^ (n + 1)))).1
        * (h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) = 1 := by
      rw [← Matrix.SpecialLinearGroup.coe_mul, inv_mul_cancel]
      rfl
    calc (h⁻¹ : SL(2, ZMod (l ^ (n + 1)))).1
        = (h⁻¹ : SL(2, ZMod (l ^ (n + 1)))).1
            * ((h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) * congElt l n (-B)) := by
          rw [hone, mul_one]
      _ = congElt l n (-B) := by rw [← mul_assoc, hinv, one_mul]

theorem mem_kerImage {H : Subgroup SL(2, ZMod (l ^ (n + 1)))} {hn : 1 ≤ n}
    {A : Matrix (Fin 2) (Fin 2) (ZMod l)} :
    A ∈ kerImage l n H hn ↔ ∃ h ∈ H, ∃ B, (h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))))
      = congElt l n B ∧ redMatHom l n B = A := Iff.rfl

/-- ★**像は跡 0 に入る**。`det = 1` が `l^n·tr(B) = 0` を強制するからである。 -/
theorem kerImage_trace_eq_zero (H : Subgroup SL(2, ZMod (l ^ (n + 1)))) (hn : 1 ≤ n)
    (A : Matrix (Fin 2) (Fin 2) (ZMod l)) (hA : A ∈ kerImage l n H hn) :
    A 0 0 + A 1 1 = 0 := by
  obtain ⟨h, _, B, hB, rfl⟩ := hA
  have hdet : (h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))).det = 1 := h.2
  rw [hB] at hdet
  have := congruence_det l n hn B
  rw [congElt] at hdet
  rw [this] at hdet
  have hzero : ((l : ZMod (l ^ (n + 1))) ^ n) * (B 0 0 + B 1 1) = 0 := by
    linear_combination hdet
  have hl0 : 0 < l := (Fact.out : Nat.Prime l).pos
  have := castPrime_eq_zero_of_pow_mul_eq_zero n hl0 hn _ hzero
  simpa [redMatHom, toPrime, map_add] using this

end KernelModule

/-! ## ★段 5 —— 随伴作用の `l^n` 個の和が消える

★★**これが `M ≠ 0` を出す機構である。**

原文の証明は `x^l` の類が持ち上げ方に依らないことを使う。その一般段は
`Σ_{i<l^n} Ad(u^i)` が `𝔰𝔩₂(𝔽_l)` の上で**ゼロ写像**であること。

★**演算子環(非可換)を経由せず、成分で計算する。**
`adU_traceZero` により 3 つのべき和 `Σ 1`・`Σ i`・`Σ i²` に落ち、
それぞれ `C(m,1)`・`C(m,2)`・`2C(m,3)+C(m,2)`(`m = l^n`)である。 -/

section AdjointSum

/-- `i² = 2·C(i,2) + i`(自然数の恒等式)。 -/
theorem sq_eq_two_mul_choose_two_add (i : ℕ) : i ^ 2 = 2 * i.choose 2 + i := by
  induction i with
  | zero => simp
  | succ k ih =>
    rw [Nat.choose_succ_succ' k 1, Nat.choose_one_right]
    ring_nf
    ring_nf at ih
    omega

/-- `Σ_{i<m} i² = 2·C(m,3) + C(m,2)`。 -/
theorem sum_range_sq (m : ℕ) :
    ∑ i ∈ Finset.range m, i ^ 2 = 2 * m.choose 3 + m.choose 2 := by
  have h : ∀ i ∈ Finset.range m, i ^ 2 = 2 * i.choose 2 + i :=
    fun i _ => sq_eq_two_mul_choose_two_add i
  rw [Finset.sum_congr rfl h, Finset.sum_add_distrib, ← Finset.mul_sum,
    sum_range_choose_eq m 2]
  have := sum_range_choose_eq m 1
  simp only [Nat.choose_one_right] at this
  rw [this]

variable (l n : ℕ) [Fact (Nat.Prime l)]

/-- ★`m = l^n` のとき `C(m,1)`・`C(m,2)`・`C(m,3)` はすべて `𝔽_l` で消える。 -/
theorem choose_pow_eq_zero (hl : Nat.Prime l) (h5 : 5 ≤ l) (hn : 1 ≤ n) (k : ℕ)
    (hk0 : k ≠ 0) (hk3 : k ≤ 3) : (((l ^ n).choose k : ℕ) : ZMod l) = 0 := by
  have hlm : l ≤ l ^ n := Nat.le_self_pow (by omega) l
  exact (ZMod.natCast_eq_zero_iff _ _).2 (hl.dvd_choose_pow hk0 (by omega))

/-- ★★**`Σ_{i<l^n} adU(i)` は跡 0 の行列の上でゼロ写像**。 -/
theorem sum_adU_eq_zero (hl : Nat.Prime l) (h5 : 5 ≤ l) (hn : 1 ≤ n)
    (X : Matrix (Fin 2) (Fin 2) (ZMod l)) (hX : X 0 0 + X 1 1 = 0) :
    ∑ i ∈ Finset.range (l ^ n), adU ((i : ℕ) : ZMod l) X = 0 := by
  have hform : X = !![X 0 0, X 0 1; X 1 0, -(X 0 0)] := eq_traceZero_form X hX
  set a := X 0 0
  set b := X 0 1
  set c := X 1 0
  -- 3 つのべき和
  have hs0 : ∑ _i ∈ Finset.range (l ^ n), (1 : ZMod l) = 0 := by
    rw [Finset.sum_const, nsmul_eq_mul, mul_one, Finset.card_range]
    have := choose_pow_eq_zero l n hl h5 hn 1 one_ne_zero (by omega)
    simpa [Nat.choose_one_right] using this
  have hs1 : ∑ i ∈ Finset.range (l ^ n), ((i : ℕ) : ZMod l) = 0 := by
    have hnat : ∑ i ∈ Finset.range (l ^ n), i = (l ^ n).choose 2 := by
      have := sum_range_choose_eq (l ^ n) 1
      simpa [Nat.choose_one_right] using this
    calc ∑ i ∈ Finset.range (l ^ n), ((i : ℕ) : ZMod l)
        = ((∑ i ∈ Finset.range (l ^ n), i : ℕ) : ZMod l) := by push_cast; ring
      _ = 0 := by rw [hnat]; exact choose_pow_eq_zero l n hl h5 hn 2 two_ne_zero (by omega)
  have hs2 : ∑ i ∈ Finset.range (l ^ n), ((i : ℕ) : ZMod l) ^ 2 = 0 := by
    calc ∑ i ∈ Finset.range (l ^ n), ((i : ℕ) : ZMod l) ^ 2
        = ((∑ i ∈ Finset.range (l ^ n), i ^ 2 : ℕ) : ZMod l) := by push_cast; ring
      _ = ((2 * (l ^ n).choose 3 + (l ^ n).choose 2 : ℕ) : ZMod l) := by rw [sum_range_sq]
      _ = 0 := by
          push_cast
          rw [choose_pow_eq_zero l n hl h5 hn 3 (by norm_num) (by omega),
            choose_pow_eq_zero l n hl h5 hn 2 two_ne_zero (by omega)]
          ring
  -- 成分ごとに計算する
  have hval : ∀ i : ℕ, adU ((i : ℕ) : ZMod l) X
      = !![a + (i : ZMod l) * c, b - 2 * a * (i : ZMod l) - c * (i : ZMod l) ^ 2;
           c, -(a + (i : ZMod l) * c)] := by
    intro i
    rw [hform]
    exact adU_traceZero _ _ _ _
  have hm0 : ((l ^ n : ℕ) : ZMod l) = 0 := by
    have := choose_pow_eq_zero l n hl h5 hn 1 one_ne_zero (by omega)
    simpa [Nat.choose_one_right] using this
  have key0 : ∑ x ∈ Finset.range (l ^ n), (a + (x : ZMod l) * c) = 0 := by
    rw [Finset.sum_add_distrib, ← Finset.sum_mul, hs1, zero_mul, add_zero,
      Finset.sum_const, Finset.card_range, nsmul_eq_mul, hm0, zero_mul]
  rw [Finset.sum_congr rfl (fun i _ => hval i)]
  ext j k
  fin_cases j <;> fin_cases k <;> rw [Matrix.sum_apply] <;> simp
  · exact key0
  · rw [zero_pow (by omega : n ≠ 0), zero_mul, ← Finset.mul_sum, hs1, mul_zero,
      ← Finset.mul_sum, hs2, mul_zero, sub_zero, sub_zero]
  · exact Or.inl (by omega)
  · calc ∑ x ∈ Finset.range (l ^ n), (-((x : ZMod l) * c) + -a)
        = ∑ x ∈ Finset.range (l ^ n), -(a + (x : ZMod l) * c) :=
          Finset.sum_congr rfl fun x _ => by ring
      _ = -∑ x ∈ Finset.range (l ^ n), (a + (x : ZMod l) * c) := by
          rw [Finset.sum_neg_distrib]
      _ = 0 := by rw [key0, neg_zero]

end AdjointSum


/-! ## ★段 6 —— 可換環上の随伴作用と、`(k·u)^m` の積公式

★`Sl2Adjoint.lean` の `adU` は**体上でしか定義されていない**(`[Field K]`)。
`ℤ/l^{n+1}` は体ではないので、同じ式を可換環上で取り直す。 -/

section AdjointRing

variable {R : Type*} [CommRing R]

/-- 上三角基本行列 `!![1,t;0,1]`(`SL` を経由しない生の行列)。 -/
def upperM (t : R) : Matrix (Fin 2) (Fin 2) R := !![1, t; 0, 1]

@[simp] theorem upperM_zero : upperM (0 : R) = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [upperM, Matrix.one_fin_two]

theorem upperM_mul (s t : R) : upperM s * upperM t = upperM (s + t) := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [upperM] <;> ring

@[simp] theorem upperM_mul_neg (t : R) : upperM t * upperM (-t) = 1 := by
  rw [upperM_mul]; simp

@[simp] theorem upperM_neg_mul (t : R) : upperM (-t) * upperM t = 1 := by
  rw [upperM_mul]; simp

theorem upperM_pow (t : R) (m : ℕ) : upperM t ^ m = upperM ((m : R) * t) := by
  induction m with
  | zero => simp
  | succ k ih => rw [pow_succ, ih, upperM_mul]; push_cast; ring_nf

/-- 可換環上の随伴作用。★体上では `Sl2Adjoint.lean` の `adU` と**同じ式**である。 -/
def adUR (t : R) (X : Matrix (Fin 2) (Fin 2) R) : Matrix (Fin 2) (Fin 2) R :=
  upperM t * X * upperM (-t)

theorem adUR_eq_adU {K : Type*} [Field K] (t : K) (X : Matrix (Fin 2) (Fin 2) K) :
    adUR t X = adU t X := rfl

end AdjointRing

section ProductFormula

variable (l : ℕ) [Fact (Nat.Prime l)] (n : ℕ)

/-- ★`upperM t` を左から通すと随伴作用が出る。 -/
theorem upperM_mul_congElt (t : ZMod (l ^ (n + 1)))
    (B : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) :
    upperM t * congElt l n B = congElt l n (adUR t B) * upperM t := by
  rw [congElt, congElt, adUR, mul_add, mul_one, add_mul, one_mul,
    Matrix.mul_smul, Matrix.smul_mul, mul_assoc, mul_assoc, upperM_neg_mul, mul_one]

/-- ★★**積公式** —— `(k·u)^m = (1 + l^n·Σ_{i<m} Ad(u^i)B) · u^m`。

★これが「`x^{l^n}` の類が持ち上げ方に依らない」ことを計算に落とす形である。 -/
theorem pow_congElt_mul_upperM (hn : 1 ≤ n)
    (B : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) (m : ℕ) :
    (congElt l n B * upperM (1 : ZMod (l ^ (n + 1)))) ^ m
      = congElt l n (∑ i ∈ Finset.range m, adUR ((i : ZMod (l ^ (n + 1)))) B)
          * upperM ((m : ZMod (l ^ (n + 1)))) := by
  induction m with
  | zero => simp
  | succ k ih =>
    rw [pow_succ (congElt l n B * upperM (1 : ZMod (l ^ (n + 1)))) k, ih]
    calc congElt l n (∑ i ∈ Finset.range k, adUR ((i : ZMod (l ^ (n + 1)))) B)
            * upperM ((k : ZMod (l ^ (n + 1))))
            * (congElt l n B * upperM (1 : ZMod (l ^ (n + 1))))
        = congElt l n (∑ i ∈ Finset.range k, adUR ((i : ZMod (l ^ (n + 1)))) B)
            * (upperM ((k : ZMod (l ^ (n + 1)))) * congElt l n B)
            * upperM (1 : ZMod (l ^ (n + 1))) := by simp only [mul_assoc]
      _ = congElt l n (∑ i ∈ Finset.range k, adUR ((i : ZMod (l ^ (n + 1)))) B)
            * (congElt l n (adUR ((k : ZMod (l ^ (n + 1)))) B)
                * upperM ((k : ZMod (l ^ (n + 1)))))
            * upperM (1 : ZMod (l ^ (n + 1))) := by rw [upperM_mul_congElt]
      _ = congElt l n (∑ i ∈ Finset.range k, adUR ((i : ZMod (l ^ (n + 1)))) B)
            * congElt l n (adUR ((k : ZMod (l ^ (n + 1)))) B)
            * (upperM ((k : ZMod (l ^ (n + 1)))) * upperM (1 : ZMod (l ^ (n + 1)))) := by
          simp only [mul_assoc]
      _ = congElt l n (∑ i ∈ Finset.range (k + 1), adUR ((i : ZMod (l ^ (n + 1)))) B)
            * upperM (((k + 1 : ℕ) : ZMod (l ^ (n + 1)))) := by
          rw [congElt_mul l n hn, upperM_mul, Finset.sum_range_succ]
          push_cast
          ring_nf

end ProductFormula


/-! ## ★段 7 —— `kerImage` の随伴不変性

★`H` が `SL₂(𝔽_l)` へ全射であることが、ここで初めて効く。
`conj_one_add_smul`(共役作用が随伴作用に落ちる)をそのまま使う。 -/

section Stability

variable (l : ℕ) [Fact (Nat.Prime l)] (n : ℕ)

@[simp] theorem redMatHom_coe (x : SL(2, ZMod (l ^ (n + 1)))) :
    redMatHom l n (x : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))))
      = ((redPow l (n + 1) (Nat.succ_ne_zero n) x : SL(2, ZMod l))
          : Matrix (Fin 2) (Fin 2) (ZMod l)) := rfl

/-- ★`kerImage` は `adU t` で不変。 -/
theorem kerImage_adU_mem (H : Subgroup SL(2, ZMod (l ^ (n + 1)))) (hn : 1 ≤ n)
    (hsurj : ∀ g : SL(2, ZMod l), ∃ h ∈ H, redPow l (n + 1) (Nat.succ_ne_zero n) h = g)
    (t : ZMod l) {A : Matrix (Fin 2) (Fin 2) (ZMod l)} (hA : A ∈ kerImage l n H hn) :
    adU t A ∈ kerImage l n H hn := by
  obtain ⟨h, hh, B, hB, rfl⟩ := hA
  obtain ⟨g, hg, hgu⟩ := hsurj (upper t)
  refine ⟨g * h * g⁻¹, H.mul_mem (H.mul_mem hg hh) (H.inv_mem hg),
    (g : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) * B
      * ((g⁻¹ : SL(2, ZMod (l ^ (n + 1)))) : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))),
    ?_, ?_⟩
  · rw [Matrix.SpecialLinearGroup.coe_mul, Matrix.SpecialLinearGroup.coe_mul, hB]
    exact conj_one_add_smul g B
  · rw [map_mul, map_mul, redMatHom_coe, redMatHom_coe, map_inv, hgu, upper_inv]
    simp only [coe_upper]
    rfl

/-- ★`kerImage` は `adL t` で不変。 -/
theorem kerImage_adL_mem (H : Subgroup SL(2, ZMod (l ^ (n + 1)))) (hn : 1 ≤ n)
    (hsurj : ∀ g : SL(2, ZMod l), ∃ h ∈ H, redPow l (n + 1) (Nat.succ_ne_zero n) h = g)
    (t : ZMod l) {A : Matrix (Fin 2) (Fin 2) (ZMod l)} (hA : A ∈ kerImage l n H hn) :
    adL t A ∈ kerImage l n H hn := by
  obtain ⟨h, hh, B, hB, rfl⟩ := hA
  obtain ⟨g, hg, hgu⟩ := hsurj (lower t)
  refine ⟨g * h * g⁻¹, H.mul_mem (H.mul_mem hg hh) (H.inv_mem hg),
    (g : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) * B
      * ((g⁻¹ : SL(2, ZMod (l ^ (n + 1)))) : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))),
    ?_, ?_⟩
  · rw [Matrix.SpecialLinearGroup.coe_mul, Matrix.SpecialLinearGroup.coe_mul, hB]
    exact conj_one_add_smul g B
  · rw [map_mul, map_mul, redMatHom_coe, redMatHom_coe, map_inv, hgu, lower_inv]
    simp only [coe_lower]
    rfl

end Stability

/-! ## ★段 8 —— `slE` が像に入る(`M ≠ 0`)

★★**ここが原文の「`x^l` の類は持ち上げ方に依らない」に対応する段である。**

`redLevel` が全射なので `π(h) = π(u)`(`u = upper 1`)なる `h ∈ H` が取れる。
`h = k·u`(`k ∈ K`)と書いて `h^{l^n}` を積公式で展開すると、
`Σ_{i<l^n} Ad(u^i)` が `𝔽_l` 上ゼロなので **`h^{l^n} = u^{l^n} = 1 + l^n·e`** になる。
★ゆえに `e ∈ M` であり、とくに **`M ≠ 0`**。 -/

section Nonvanishing

variable (l : ℕ) [Fact (Nat.Prime l)] (n : ℕ)

theorem coe_pow_sl {R : Type*} [CommRing R] (g : SL(2, R)) (m : ℕ) :
    ((g ^ m : SL(2, R)) : Matrix (Fin 2) (Fin 2) R)
      = (g : Matrix (Fin 2) (Fin 2) R) ^ m := by
  induction m with
  | zero => simp
  | succ k ih => rw [pow_succ, pow_succ, Matrix.SpecialLinearGroup.coe_mul, ih]

@[simp] theorem redMatHom_upperM (t : ZMod (l ^ (n + 1))) :
    redMatHom l n (upperM t) = upperM (toPrime l n t) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [upperM, redMatHom, RingHom.mapMatrix]

theorem redMatHom_adUR (t : ZMod (l ^ (n + 1)))
    (B : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) :
    redMatHom l n (adUR t B) = adUR (toPrime l n t) (redMatHom l n B) := by
  rw [adUR, adUR, map_mul, map_mul, redMatHom_upperM, ← map_neg, redMatHom_upperM]

/-- ★`mod l` で消える `S` については `1 + l^n·S = 1`。 -/
theorem congElt_eq_one_of_redMatHom_eq_zero (hn : 1 ≤ n)
    (S : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) (hS : redMatHom l n S = 0) :
    congElt l n S = 1 := by
  have hl0 : 0 < l := (Fact.out : Nat.Prime l).pos
  have hzero : ((l : ZMod (l ^ (n + 1)))) ^ n • S = 0 := by
    ext i j
    have h0 : toPrime l n (S i j) = 0 := by
      have := congrFun (congrFun hS i) j
      simpa [redMatHom, RingHom.mapMatrix] using this
    have := pow_mul_congr_of_castPrime_eq n hl0 (S i j) 0 (by simpa [toPrime] using h0)
    simpa using this
  rw [congElt, hzero, add_zero]

/-- ★`(x : SL)` が `congElt B` の形なら `redMatHom B` は跡 0。 -/
theorem trace_redMatHom_eq_zero (hn : 1 ≤ n) (x : SL(2, ZMod (l ^ (n + 1))))
    (B : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))))
    (hB : (x : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) = congElt l n B) :
    (redMatHom l n B) 0 0 + (redMatHom l n B) 1 1 = 0 := by
  have hdet : (x : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))).det = 1 := x.2
  rw [hB, congElt, congruence_det l n hn B] at hdet
  have hzero : ((l : ZMod (l ^ (n + 1))) ^ n) * (B 0 0 + B 1 1) = 0 := by
    linear_combination hdet
  have hl0 : 0 < l := (Fact.out : Nat.Prime l).pos
  have := castPrime_eq_zero_of_pow_mul_eq_zero n hl0 hn _ hzero
  simpa [redMatHom, RingHom.mapMatrix, toPrime, map_add] using this

/-- 可換環上の `e = !![0,1;0,0]`(`Sl2Adjoint.lean` の `slE` は体上のみ)。 -/
def slER {R : Type*} [CommRing R] : Matrix (Fin 2) (Fin 2) R := !![0, 1; 0, 0]

theorem upperM_eq_congElt (l n : ℕ) [Fact (Nat.Prime l)] :
    upperM (((l ^ n : ℕ) : ZMod (l ^ (n + 1)))) = congElt l n slER := by
  have hcast : ((l ^ n : ℕ) : ZMod (l ^ (n + 1))) = (l : ZMod (l ^ (n + 1))) ^ n := by
    push_cast; ring
  rw [congElt, hcast]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [upperM, slER, Matrix.one_fin_two]

theorem redMatHom_slER (l n : ℕ) [Fact (Nat.Prime l)] :
    redMatHom l n slER = (slE : Matrix (Fin 2) (Fin 2) (ZMod l)) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [redMatHom, RingHom.mapMatrix, slER, slE, toPrime]

/-- ★★**`slE` は像に入る** —— `M ≠ 0` の中身。 -/
theorem slE_mem_kerImage (H : Subgroup SL(2, ZMod (l ^ (n + 1)))) (hn : 1 ≤ n) (h5 : 5 ≤ l)
    (hfull : ∀ y : SL(2, ZMod (l ^ n)), ∃ h ∈ H, redLevel l (Nat.le_succ n) h = y) :
    (slE : Matrix (Fin 2) (Fin 2) (ZMod l)) ∈ kerImage l n H hn := by
  have hl : Nat.Prime l := Fact.out
  have hl0 : 0 < l := hl.pos
  obtain ⟨h, hh, hhu⟩ := hfull
    (redLevel l (Nat.le_succ n) (upper (1 : ZMod (l ^ (n + 1)))))
  have hker : redLevel l (Nat.le_succ n) (h * (upper (1 : ZMod (l ^ (n + 1))))⁻¹) = 1 := by
    rw [map_mul, map_inv, hhu, mul_inv_cancel]
  have hentry : ∀ i j, ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n))
      (((h * (upper (1 : ZMod (l ^ (n + 1))))⁻¹ : SL(2, ZMod (l ^ (n + 1))))
          : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j
        - (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j) = 0 := by
    intro i j
    have h1 := congrFun (congrFun (congrArg
      (fun x : SL(2, ZMod (l ^ n)) => (x : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n)))) hker) i) j
    rw [redLevel_coe] at h1
    have h2 : ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n))
        ((1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j)
        = (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n))) i j := by
      by_cases hij : i = j
      · subst hij; rw [Matrix.one_apply_eq, Matrix.one_apply_eq, map_one]
      · simp [Matrix.one_apply_ne hij]
    rw [map_sub, h1, h2]
    show ((1 : SL(2, ZMod (l ^ n))) : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n))) i j
      - (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n))) i j = 0
    simp
  obtain ⟨B, hB⟩ := exists_matrix_of_castHom_eq_zero l n hl0 _ hentry
  have hcoe : (h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))))
      = congElt l n B * upperM (1 : ZMod (l ^ (n + 1))) := by
    have hhk : h = (h * (upper (1 : ZMod (l ^ (n + 1))))⁻¹) * upper (1 : ZMod (l ^ (n + 1))) := by
      rw [inv_mul_cancel_right]
    rw [hhk, Matrix.SpecialLinearGroup.coe_mul, hB]
    rfl
  have htr : (redMatHom l n B) 0 0 + (redMatHom l n B) 1 1 = 0 :=
    trace_redMatHom_eq_zero l n hn _ B hB
  have hSzero : redMatHom l n
      (∑ i ∈ Finset.range (l ^ n), adUR ((i : ZMod (l ^ (n + 1)))) B) = 0 := by
    rw [map_sum]
    have hcongr : ∀ i ∈ Finset.range (l ^ n),
        redMatHom l n (adUR ((i : ZMod (l ^ (n + 1)))) B)
          = adU ((i : ℕ) : ZMod l) (redMatHom l n B) := by
      intro i _
      rw [redMatHom_adUR]
      congr 1
      simp [toPrime]
    rw [Finset.sum_congr rfl hcongr]
    exact sum_adU_eq_zero l n hl h5 hn _ htr
  refine ⟨h ^ (l ^ n), H.pow_mem hh _, slER, ?_, redMatHom_slER l n⟩
  rw [coe_pow_sl, hcoe, pow_congElt_mul_upperM l n hn B (l ^ n),
    congElt_eq_one_of_redMatHom_eq_zero l n hn _ hSzero, one_mul, upperM_eq_congElt]

end Nonvanishing

/-! ## ★段 9 —— 有限段の本体

★`slE ∈ M` と随伴不変性から `sl2_adjoint_irreducible` で `M = 𝔰𝔩₂` が出る。
そこから合同核が `H` に入り、帰納法の仮定と併せて `H = ⊤`。 -/

section FiniteLevelMain

variable (l : ℕ) [Fact (Nat.Prime l)]

theorem two_ne_zero_of_five_le (h5 : 5 ≤ l) : (2 : ZMod l) ≠ 0 := by
  have hl : Nat.Prime l := Fact.out
  intro hcon
  have : ((2 : ℕ) : ZMod l) = 0 := by exact_mod_cast hcon
  have hdvd : l ∣ 2 := (ZMod.natCast_eq_zero_iff _ _).1 this
  have := Nat.le_of_dvd (by norm_num) hdvd
  omega

theorem slE_ne_zero (h5 : 5 ≤ l) : (slE : Matrix (Fin 2) (Fin 2) (ZMod l)) ≠ 0 := by
  intro hcon
  have h1 : (slE : Matrix (Fin 2) (Fin 2) (ZMod l)) 0 1 = 0 := by rw [hcon]; rfl
  have h2 : (slE : Matrix (Fin 2) (Fin 2) (ZMod l)) 0 1 = 1 := rfl
  have h3 : (1 : ZMod l) = 0 := by rw [← h2, h1]
  haveI : Fact (1 < l) := ⟨by omega⟩
  exact one_ne_zero h3

variable (n : ℕ)

/-- ★★**合同核は `H` に含まれる** —— 既約性の帰結。 -/
theorem ker_le_of_slE_mem (H : Subgroup SL(2, ZMod (l ^ (n + 1)))) (hn : 1 ≤ n) (h5 : 5 ≤ l)
    (hsurj : ∀ g : SL(2, ZMod l), ∃ h ∈ H, redPow l (n + 1) (Nat.succ_ne_zero n) h = g)
    (hslE : (slE : Matrix (Fin 2) (Fin 2) (ZMod l)) ∈ kerImage l n H hn)
    (x : SL(2, ZMod (l ^ (n + 1))))
    (hx : redLevel l (Nat.le_succ n) x = 1) : x ∈ H := by
  have hl : Nat.Prime l := Fact.out
  have hl0 : 0 < l := hl.pos
  -- `M` を `𝔽_l`-部分加群として見る
  set V : Submodule (ZMod l) (Matrix (Fin 2) (Fin 2) (ZMod l)) :=
    AddSubgroup.toZModSubmodule l (kerImage l n H hn) with hV
  have hmemV : ∀ A : Matrix (Fin 2) (Fin 2) (ZMod l), A ∈ V ↔ A ∈ kerImage l n H hn :=
    fun A => Iff.rfl
  have hU : ∀ (t : ZMod l) {X : Matrix (Fin 2) (Fin 2) (ZMod l)}, X ∈ V → adU t X ∈ V :=
    fun t _ hX => (hmemV _).2 (kerImage_adU_mem l n H hn hsurj t ((hmemV _).1 hX))
  have hL : ∀ (t : ZMod l) {X : Matrix (Fin 2) (Fin 2) (ZMod l)}, X ∈ V → adL t X ∈ V :=
    fun t _ hX => (hmemV _).2 (kerImage_adL_mem l n H hn hsurj t ((hmemV _).1 hX))
  have htr : ∀ X ∈ V, X 0 0 + X 1 1 = 0 :=
    fun X hX => kerImage_trace_eq_zero l n H hn X ((hmemV _).1 hX)
  obtain ⟨heV, hfV, hhV⟩ := sl2_adjoint_irreducible V hU hL htr (two_ne_zero_of_five_le l h5)
    ((hmemV _).2 hslE) (slE_ne_zero l h5)
  -- 跡 0 の行列はすべて `V` に入る
  have hall : ∀ A : Matrix (Fin 2) (Fin 2) (ZMod l), A 0 0 + A 1 1 = 0 → A ∈ V := by
    intro A hA
    have hform : A = A 0 0 • (slH : Matrix (Fin 2) (Fin 2) (ZMod l))
        + A 0 1 • (slE : Matrix (Fin 2) (Fin 2) (ZMod l))
        + A 1 0 • (slF : Matrix (Fin 2) (Fin 2) (ZMod l)) := by
      ext i j
      fin_cases i <;> fin_cases j <;> simp [slE, slF, slH] <;> linear_combination hA
    rw [hform]
    exact V.add_mem (V.add_mem (V.smul_mem _ hhV) (V.smul_mem _ heV)) (V.smul_mem _ hfV)
  -- `x` を `congElt` の形に書く
  have hentry : ∀ i j, ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n))
      ((x : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j
        - (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j) = 0 := by
    intro i j
    have h1 := congrFun (congrFun (congrArg
      (fun z : SL(2, ZMod (l ^ n)) => (z : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n)))) hx) i) j
    rw [redLevel_coe] at h1
    have h2 : ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n))
        ((1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1)))) i j)
        = (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n))) i j := by
      by_cases hij : i = j
      · subst hij; rw [Matrix.one_apply_eq, Matrix.one_apply_eq, map_one]
      · simp [Matrix.one_apply_ne hij]
    rw [map_sub, h1, h2]
    show ((1 : SL(2, ZMod (l ^ n))) : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n))) i j
      - (1 : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n))) i j = 0
    simp
  obtain ⟨B, hB⟩ := exists_matrix_of_castHom_eq_zero l n hl0 _ hentry
  have htrB : (redMatHom l n B) 0 0 + (redMatHom l n B) 1 1 = 0 :=
    trace_redMatHom_eq_zero l n hn x B hB
  obtain ⟨h, hh, B', hB', hBB'⟩ := (hmemV _).1 (hall _ htrB)
  have hB'' : (h : Matrix (Fin 2) (Fin 2) (ZMod (l ^ (n + 1))))
      = 1 + ((l : ZMod (l ^ (n + 1)))) ^ n • B' := hB'
  -- ★`l^n·B = l^n·B'` —— `mod l` で一致していれば十分である
  have hsm : ((l : ZMod (l ^ (n + 1)))) ^ n • B = ((l : ZMod (l ^ (n + 1)))) ^ n • B' := by
    ext i j
    have hpt : toPrime l n (B i j) = toPrime l n (B' i j) := by
      have := congrFun (congrFun hBB' i) j
      simpa [redMatHom, RingHom.mapMatrix] using this.symm
    simpa using
      pow_mul_congr_of_castPrime_eq n hl0 (B i j) (B' i j) (by simpa [toPrime] using hpt)
  have hxh : x = h := by
    apply Subtype.ext
    rw [hB, hB'', hsm]
  rw [hxh]
  exact hh

/-- ★★**段 (A) の本体**(`n ≥ 1`)。

`SL₂(ℤ/l^n)` の部分群が `SL₂(𝔽_l)` へ全射なら、それは全体である。 -/
theorem subgroup_eq_top_of_redPow_surj (h5 : 5 ≤ l) :
    ∀ (k : ℕ) (H : Subgroup SL(2, ZMod (l ^ (k + 1)))),
      (∀ g : SL(2, ZMod l), ∃ h ∈ H, redPow l (k + 1) (Nat.succ_ne_zero k) h = g) → H = ⊤ := by
  have hl : Nat.Prime l := Fact.out
  haveI : NeZero l := ⟨hl.ne_zero⟩
  intro k
  induction k with
  | zero =>
    -- `l^1` の段: `redPow` は単射
    intro H hsurj
    have hinj : Function.Injective (redPow l 1 (Nat.succ_ne_zero 0)) := by
      intro x y hxy
      apply Subtype.ext
      ext i j
      have := congrFun (congrFun (congrArg
        (fun z : SL(2, ZMod l) => (z : Matrix (Fin 2) (Fin 2) (ZMod l))) hxy) i) j
      exact castHom_pow_one_injective l this
    refine (Subgroup.eq_top_iff' H).2 fun g => ?_
    obtain ⟨h, hh, hhg⟩ := hsurj (redPow l 1 (Nat.succ_ne_zero 0) g)
    rwa [hinj hhg] at hh
  | succ n ih =>
    intro H hsurj
    have hn : 1 ≤ n + 1 := Nat.one_le_iff_ne_zero.2 (Nat.succ_ne_zero n)
    -- 帰納法の仮定で 1 段下は全体
    have hlower : H.map (redLevel l (Nat.le_succ (n + 1))) = ⊤ := by
      refine ih _ ?_
      intro g
      obtain ⟨h, hh, hhg⟩ := hsurj g
      refine ⟨redLevel l (Nat.le_succ (n + 1)) h, ⟨h, hh, rfl⟩, ?_⟩
      have := congrFun (congrArg (fun f : SL(2, ZMod (l ^ (n + 2))) →* SL(2, ZMod l) => f.toFun)
        (redPow_comp_redLevel l (Nat.succ_ne_zero n))) h
      simpa [hhg] using this
    have hfull : ∀ y : SL(2, ZMod (l ^ (n + 1))),
        ∃ h ∈ H, redLevel l (Nat.le_succ (n + 1)) h = y := by
      intro y
      have : y ∈ H.map (redLevel l (Nat.le_succ (n + 1))) := by rw [hlower]; trivial
      exact this
    -- `slE ∈ M` から合同核が `H` に入る
    have hslE := slE_mem_kerImage l (n + 1) H hn h5 hfull
    have hker := ker_le_of_slE_mem l (n + 1) H hn h5 hsurj hslE
    refine (Subgroup.eq_top_iff' H).2 fun g => ?_
    obtain ⟨h, hh, hhg⟩ := hfull (redLevel l (Nat.le_succ (n + 1)) g)
    have hmem : g * h⁻¹ ∈ H := by
      refine hker _ ?_
      rw [map_mul, map_inv, hhg, mul_inv_cancel]
    have : g = (g * h⁻¹) * h := by rw [inv_mul_cancel_right]
    rw [this]
    exact H.mul_mem hmem hh

end FiniteLevelMain

end ABC3.Found.GenEll
