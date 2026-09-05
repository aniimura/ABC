import ABC3.Found.Falt1.AlmostDifferentials
import ABC3.Found.Falt1.HochschildLowDegree

/-!
# [Falt1] `p`-可除な底と almost zero(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2 冒頭
(物理 p.5-6=印字 p.258-259)——`A` は `p` のすべての `p` 冪根を含み、
`m = ∪_ε p^ε A`(`ε` は `p` 冪分母の正有理数)を「almost zero」の
判定に使うイデアルとする、という設定。

## なぜこの層が要るか(2026-09-05 の構造上の発見)

本プロジェクトの `Definition 2.1` 条件(iii)は指数を `n : ℕ` で取っており、
そこから得られる主張は「`p^n` が零化する」という形になる。これで
**`Theorem 2.4`(almost な結論)は閉じた**が、**`Theorem 2.2`/`2.3`
(honest な結論)は閉じない**——`p^ε φ_ε` の族から `ε→0` の極限として
honest な `φ₀` を得る操作が要り、整数指数では `p` で割れないからである。

本ファイルはその不足を埋める最小の層を用意する:
`ϖ : ℕ → A` で `(ϖ (k+1))^q = ϖ k`(`ϖ k` が「`p^{1/q^k}`」)という
**塔**を持ち、`m := (ϖ k : k)` を Faltings の `m` とする。

## 内容

- `PDivTower`:塔そのもの。`dvd_of_le` で「添字が大きいほど深く割れる」。
- `m` と **`m² = m`**(`m_sq_eq`)——almost ring theory の
  *basic setup* が要求する冪等性。`ϖ k = ϖ(k+1)·(ϖ(k+1))^{q-1}` から出る。
- `IsAlmostZero`:`m` が零化すること。
- 既存の結果(`kaehler_almost_zero`・`hochschild_H2_almost_coboundary`)を
  この言葉に翻訳したもの。前コミットで下流補題を一般スカラー化したので、
  翻訳は各 `ϖ k` に適用するだけで済む。
-/

namespace ABC3.Found.Falt1

universe u v

/-- **`p`-可除な塔**。`ϖ k` は `p := ϖ 0` の `q^k` 乗根に相当し、
`(ϖ (k+1))^q = ϖ k` を満たす。Faltings の `m = ∪_ε p^ε` を
有限的(可算な塔)に表現する最小の道具。 -/
structure PDivTower (A : Type u) [CommRing A] (q : ℕ) where
  /-- `ϖ k` は「`p^{1/q^k}`」に相当する元。 -/
  ϖ : ℕ → A
  /-- `ϖ (k+1)` の `q` 乗が `ϖ k`。 -/
  ϖ_succ : ∀ k, (ϖ (k + 1)) ^ q = ϖ k

namespace PDivTower

variable {A : Type u} [CommRing A] {q : ℕ} (T : PDivTower A q)

theorem dvd_succ (hq : 1 ≤ q) (k : ℕ) : T.ϖ (k + 1) ∣ T.ϖ k := by
  rw [← T.ϖ_succ k]
  exact dvd_pow_self _ (by omega)

/-- 添字が大きいほど「小さい」(= より深く割れる)。`ε` が単調に減る
という Faltings の直観に対応する。 -/
theorem dvd_of_le (hq : 1 ≤ q) {j k : ℕ} (h : j ≤ k) : T.ϖ k ∣ T.ϖ j := by
  induction k with
  | zero => simp_all
  | succ k ih =>
    rcases Nat.lt_or_ge j (k + 1) with hj | hj
    · exact dvd_trans (T.dvd_succ hq k) (ih (by omega))
    · have hjk : j = k + 1 := by omega
      subst hjk
      exact dvd_rfl

/-- **Faltings の `m`** —— 塔の元すべてが生成するイデアル。 -/
def m : Ideal A := Ideal.span (Set.range T.ϖ)

theorem ϖ_mem_m (k : ℕ) : T.ϖ k ∈ T.m := Ideal.subset_span ⟨k, rfl⟩

/-- **`m` は冪等**(`m² = m`)——almost ring theory の *basic setup* が
要求する条件。`ϖ k = ϖ(k+1)·(ϖ(k+1))^{q-1}` と `q ≥ 2` から、
各生成元が `m²` に入る。 -/
theorem m_sq_eq (hq : 2 ≤ q) : T.m ^ 2 = T.m := by
  refine le_antisymm (by rw [sq]; exact Ideal.mul_le_right) ?_
  rw [PDivTower.m, Ideal.span_le]
  rintro _ ⟨k, rfl⟩
  have h : T.ϖ k = T.ϖ (k + 1) * (T.ϖ (k + 1)) ^ (q - 1) := by
    rw [← pow_succ', show q - 1 + 1 = q by omega]
    exact (T.ϖ_succ k).symm
  rw [h, sq]
  exact Ideal.mul_mem_mul (T.ϖ_mem_m (k + 1))
    (Ideal.pow_mem_of_mem _ (T.ϖ_mem_m (k + 1)) _ (by omega))

/-- **almost zero**: `m` が `M` を零化する(= 各 `ϖ k` が零化する)。
Faltings の「`m` annihilates …」はこの述語である。 -/
def IsAlmostZero (M : Type v) [AddCommGroup M] [Module A M] : Prop :=
  ∀ (k : ℕ) (x : M), T.ϖ k • x = 0

theorem isAlmostZero_iff_m_smul_eq_bot (M : Type v) [AddCommGroup M] [Module A M] :
    T.IsAlmostZero M ↔ ∀ (x : M) (a : A), a ∈ T.m → a • x = 0 := by
  constructor
  · intro h x a ha
    refine Submodule.span_induction ?_ ?_ ?_ ?_ ha
    · rintro _ ⟨k, rfl⟩; exact h k x
    · simp
    · intro y z _ _ hy hz; rw [add_smul, hy, hz, add_zero]
    · intro y z _ hz; rw [smul_eq_mul, mul_smul, hz, smul_zero]
  · intro h k x
    exact h x _ (T.ϖ_mem_m k)

/-- **塔に沿った witness の族**——`Definition 2.1` 条件(iii)を各 `ϖ k`
について要求したもの(`w` の annihilation と augmentation だけを抜き出した
形にしてある。`almost_swap_annihilate`・`almost_swap_augment` がこの形の
witness を供給する)。 -/
def HasAlmostWitnesses (B : Type u) [CommRing B] [Algebra A B] : Prop :=
  ∀ k : ℕ, ∃ w : TensorProduct A B B,
    (∀ b : B, (1 ⊗ₜ[A] b - b ⊗ₜ[A] 1) * w = 0) ∧
    Algebra.TensorProduct.lmul' A w = T.ϖ k • (1 : B)

/-- **`Ω[B⁄A]` は almost zero**(`Theorem 2.4(i)` の余核側を、
almost mathematics の言葉で述べたもの)。 -/
theorem kaehler_isAlmostZero {B : Type u} [CommRing B] [Algebra A B]
    (hwit : T.HasAlmostWitnesses B) : T.IsAlmostZero (Ω[B⁄A]) := by
  intro k x
  obtain ⟨w, hann, haug⟩ := hwit k
  exact kaehler_almost_zero (T.ϖ k) w hann haug x

/-- **Hochschild `H²` は almost zero**(`Theorem 2.2` の第2段の入力を、
almost mathematics の言葉で述べたもの):任意の 2-コサイクル `c` と
任意の `k` について、`ϖ k · c` はコバウンダリである。 -/
theorem hochschild_H2_isAlmostCoboundary {B M : Type u} [CommRing B] [Algebra A B]
    [AddCommGroup M] [Module B M] [Module A M] [IsScalarTower A B M]
    (hwit : T.HasAlmostWitnesses B)
    (c : B →ₗ[A] B →ₗ[A] M)
    (hc : ∀ v b₁ b₂ : B, v • c b₁ b₂ - c (v * b₁) b₂ + c v (b₁ * b₂) - b₂ • c v b₁ = 0)
    (k : ℕ) :
    ∃ h : B →ₗ[A] M, ∀ b₁ b₂ : B,
      T.ϖ k • c b₁ b₂ = b₁ • h b₂ - h (b₁ * b₂) + b₂ • h b₁ := by
  obtain ⟨w, hann, haug⟩ := hwit k
  exact hochschild_H2_almost_coboundary (T.ϖ k) w hann haug c hc

end PDivTower

/-! ## 非空虚性の対照

塔が実際に存在すること、しかも **`m` が真のイデアル**(`m ≠ ⊤`)に
なるような非退化な例があることを示す。 -/

/-- `q = 1` の塔(`ϖ` は定数)——**現行の枠組み(単一の `p` で進めてきた
これまでの形式化)がちょうどこの退化例に当たる**。`m = (5)` は真の
非自明イデアルだが、`q = 1` なので `m² = m` は成り立たない
(`m_sq_eq` は `q ≥ 2` を要求する)。 -/
example : PDivTower ℤ 1 := { ϖ := fun _ => 5, ϖ_succ := fun _ => by ring }

/-- `X_{k+1}^q - X_k` で割った「完全化」——`ϖ k := X_k` が塔をなす。
`q ≥ 2` ならこれが almost ring theory の *basic setup*(`m² = m` かつ
`m ≠ ⊤`)の非退化な実例になる。 -/
noncomputable def perfRel (q : ℕ) : Ideal (MvPolynomial ℕ ℤ) :=
  Ideal.span (Set.range fun k : ℕ => (MvPolynomial.X (k+1) : MvPolynomial ℕ ℤ) ^ q
    - MvPolynomial.X k)

noncomputable def perfTower (q : ℕ) : PDivTower (MvPolynomial ℕ ℤ ⧸ perfRel q) q where
  ϖ k := Ideal.Quotient.mk _ (MvPolynomial.X k)
  ϖ_succ k := by
    rw [← map_pow, Ideal.Quotient.eq]
    exact Ideal.subset_span ⟨k, rfl⟩

open MvPolynomial in
/-- **完全化の `m` は真のイデアル**。すべての `X_k` を `0` に送る評価写像が
関係式を殺すので商を経由し、`m` を `0` に送る。したがって `1 ∉ m`。 -/
theorem perfTower_m_ne_top (q : ℕ) (hq : q ≠ 0) : (perfTower q).m ≠ ⊤ := by
  have hker : ∀ x ∈ perfRel q, (eval (fun _ : ℕ => (0:ℤ))) x = 0 := by
    intro x hx
    refine Submodule.span_induction ?_ ?_ ?_ ?_ hx
    · rintro _ ⟨k, rfl⟩
      rw [map_sub, map_pow, eval_X, eval_X, zero_pow hq, sub_zero]
    · rw [map_zero]
    · intro y z _ _ hy hz; rw [map_add, hy, hz, add_zero]
    · intro a y _ hy; rw [smul_eq_mul, map_mul, hy, mul_zero]
  set g : (MvPolynomial ℕ ℤ ⧸ perfRel q) →+* ℤ :=
    Ideal.Quotient.lift _ (eval (fun _ : ℕ => (0:ℤ))) hker with hg
  intro htop
  have h1 : (1 : MvPolynomial ℕ ℤ ⧸ perfRel q) ∈ (perfTower q).m := htop ▸ Submodule.mem_top
  have hzero : ∀ y ∈ (perfTower q).m, g y = 0 := by
    intro y hy
    refine Submodule.span_induction ?_ ?_ ?_ ?_ hy
    · rintro _ ⟨k, rfl⟩
      show g (Ideal.Quotient.mk _ (X k)) = 0
      rw [hg, Ideal.Quotient.lift_mk, eval_X]
    · rw [map_zero]
    · intro u v _ _ hu hv; rw [map_add, hu, hv, add_zero]
    · intro a u _ hu; rw [smul_eq_mul, map_mul, hu, mul_zero]
  have hcontra := hzero 1 h1
  rw [map_one] at hcontra
  exact one_ne_zero hcontra

/-- **`q ≥ 2` の完全化は almost ring theory の *basic setup* の非退化な実例**:
`m² = m` かつ `m ≠ ⊤`。 -/
theorem perfTower_basicSetup (q : ℕ) (hq : 2 ≤ q) :
    (perfTower q).m ^ 2 = (perfTower q).m ∧ (perfTower q).m ≠ ⊤ :=
  ⟨PDivTower.m_sq_eq _ hq, perfTower_m_ne_top q (by omega)⟩

end ABC3.Found.Falt1
