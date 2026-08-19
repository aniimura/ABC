import ABC3.Found.GaloisRep.TateWitness

/-!
# Galois (G3) 第 78 ブロック —— **★★★★加法写像は自動的に `ℤ_l` 線型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★行列が出るための鍵

Galois 群は `T_l E` に**加法的に**作用する。★原文が言う `GL₂(ℤ_l)` に値を取るには
その作用が **`ℤ_l` 線型**でなければならない——普通はここで連続性を要求する。

★★しかし `ℤ_l` 加群では**線型性は自動である**:

    f(p^n y) = p^n f(y)   (加法写像だから)

★★★これが「`f` は `p` 進的に連続」を意味し、`ℤ` が `ℤ_p` で稠密だから線型になる。
★★★★書き下すと `α = α_n + p^n γ`(`α_n ∈ ℤ`)と分けて

    f(αx) − α f(x) = p^n (f(γx) − γ f(x)) ∈ p^n N

——★これが全ての `n` で成り立つので、分離性から `0` である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `padic_dvd_of_toZModPow` | ★`toZModPow n x = 0` ⟹ `p^n ∣ x` |
| `addMonoidHom_padic_smul` | ★★★★**加法写像は `ℤ_p` 線型** |
| `padic_pi_sep` | ★`ℤ_p^k` の分離性 |
| `exists_matrix_of_addHom` | ★★★★**加法自己準同型は行列** |
| `matOf` ほか | ★★行列の割り当てと、その合成・恒等 |
| `prodEquivPiTwo` | ★`ℤ_l × ℤ_l ≃+ (Fin 2 → ℤ_l)` |
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★`toZModPow n x = 0` なら `p^n` で割れる。 -/
theorem padic_dvd_of_toZModPow (p : ℕ) [Fact p.Prime] (n : ℕ) (x : ℤ_[p])
    (h : PadicInt.toZModPow n x = 0) : ∃ γ : ℤ_[p], x = (p : ℤ_[p]) ^ n * γ := by
  have hmem : x ∈ RingHom.ker (PadicInt.toZModPow (p := p) n) := h
  rw [PadicInt.ker_toZModPow, Ideal.mem_span_singleton] at hmem
  obtain ⟨γ, hγ⟩ := hmem
  exact ⟨γ, hγ⟩

/-- ★★★★**`ℤ_p` 加群のあいだの加法写像は自動的に `ℤ_p` 線型である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★連続性を仮定する必要がない——`f(p^n y) = p^n f(y)` がそれを与える。 -/
theorem addMonoidHom_padic_smul {p : ℕ} [Fact p.Prime] {M N : Type u}
    [AddCommGroup M] [AddCommGroup N] [Module ℤ_[p] M] [Module ℤ_[p] N]
    (hsep : ∀ y : N, (∀ n : ℕ, ∃ z : N, y = ((p : ℤ_[p]) ^ n) • z) → y = 0)
    (f : M →+ N) (α : ℤ_[p]) (x : M) : f (α • x) = α • f x := by
  have hnat : ∀ (n : ℕ) (y : M), f (((p : ℤ_[p]) ^ n) • y) = ((p : ℤ_[p]) ^ n) • f y := by
    intro n y
    have h1 : ((p : ℤ_[p]) ^ n) = (((p ^ n : ℕ) : ℤ_[p])) := by push_cast; ring
    rw [h1, Nat.cast_smul_eq_nsmul, map_nsmul, Nat.cast_smul_eq_nsmul]
  have hint : ∀ (c : ℤ) (y : M), f ((c : ℤ_[p]) • y) = (c : ℤ_[p]) • f y := by
    intro c y
    rw [Int.cast_smul_eq_zsmul, map_zsmul, Int.cast_smul_eq_zsmul]
  have hd : f (α • x) - α • f x = 0 := by
    refine hsep _ (fun n => ?_)
    have hz : PadicInt.toZModPow n (α - ((padicRep p α n : ℤ) : ℤ_[p])) = 0 := by
      rw [map_sub, map_intCast, padicRep_cast]
      ring
    obtain ⟨γ, hγ⟩ := padic_dvd_of_toZModPow p n _ hz
    have hsplit : α = ((padicRep p α n : ℤ) : ℤ_[p]) + (p : ℤ_[p]) ^ n * γ := by
      rw [← hγ]; ring
    refine ⟨f (γ • x) - γ • f x, ?_⟩
    rw [hsplit]
    simp only [add_smul, map_add, mul_smul, hint, hnat, smul_sub]
    abel
  exact sub_eq_zero.mp hd

/-- ★`ℤ_p^k` は `p` 進的に分離している。 -/
theorem padic_pi_sep {p : ℕ} [Fact p.Prime] {k : ℕ} (y : Fin k → ℤ_[p])
    (h : ∀ n : ℕ, ∃ z : Fin k → ℤ_[p], y = ((p : ℤ_[p]) ^ n) • z) : y = 0 := by
  funext i
  refine PadicInt.ext_of_toZModPow.1 (fun n => ?_)
  obtain ⟨z, hz⟩ := h n
  have hyi : y i = (p : ℤ_[p]) ^ n * z i := by rw [hz]; simp
  have hp0 : PadicInt.toZModPow (p := p) n ((p : ℤ_[p]) ^ n) = 0 := by
    rw [map_pow, map_natCast, ← Nat.cast_pow, ZMod.natCast_self]
  rw [hyi, map_mul, hp0, zero_mul]
  simp

/-- ★★★★**`ℤ_p^k` の加法自己準同型は行列である**。 -/
theorem exists_matrix_of_addHom {p : ℕ} [Fact p.Prime] {k : ℕ}
    (f : (Fin k → ℤ_[p]) →+ (Fin k → ℤ_[p])) :
    ∃ M : Matrix (Fin k) (Fin k) ℤ_[p], ∀ x, f x = Matrix.mulVec M x := by
  have hlin : ∀ (α : ℤ_[p]) (x : Fin k → ℤ_[p]), f (α • x) = α • f x :=
    fun α x => addMonoidHom_padic_smul padic_pi_sep f α x
  refine ⟨Matrix.of (fun i j => f (Pi.single j 1) i), fun x => ?_⟩
  have hx : x = ∑ j : Fin k, x j • (Pi.single j 1 : Fin k → ℤ_[p]) := by
    funext i
    simp [Finset.sum_apply, Pi.single_apply]
  funext i
  calc f x i = f (∑ j : Fin k, x j • (Pi.single j 1 : Fin k → ℤ_[p])) i := by rw [← hx]
    _ = (∑ j : Fin k, f (x j • (Pi.single j 1 : Fin k → ℤ_[p]))) i := by rw [map_sum]
    _ = ∑ j : Fin k, (x j • f (Pi.single j 1 : Fin k → ℤ_[p])) i := by
        simp only [hlin, Finset.sum_apply]
    _ = Matrix.mulVec (Matrix.of (fun i j => f (Pi.single j 1) i)) x i := by
        simp [Matrix.mulVec, dotProduct, mul_comm]

/-! ## ★★行列の割り当て -/

/-- ★ベクトルへの作用が同じ行列は等しい。 -/
theorem matrix_ext_of_mulVec {p : ℕ} [Fact p.Prime] {k : ℕ}
    {M N : Matrix (Fin k) (Fin k) ℤ_[p]} (h : ∀ x, Matrix.mulVec M x = Matrix.mulVec N x) :
    M = N := by
  ext i j
  have h1 := congrFun (h (Pi.single j 1)) i
  simpa [Matrix.mulVec, dotProduct, Pi.single_apply] using h1

/-- ★★`T_l` の加法自己準同型に対応する行列。 -/
noncomputable def matOf {A : Type u} [AddCommGroup A] {l : ℕ} [Fact l.Prime]
    (e : limTors A l ≃+ (Fin 2 → ℤ_[l])) (T : limTors A l →+ limTors A l) :
    Matrix (Fin 2) (Fin 2) ℤ_[l] :=
  Classical.choose (exists_matrix_of_addHom
    ((e.toAddMonoidHom.comp T).comp e.symm.toAddMonoidHom))

theorem matOf_apply {A : Type u} [AddCommGroup A] {l : ℕ} [Fact l.Prime]
    (e : limTors A l ≃+ (Fin 2 → ℤ_[l])) (T : limTors A l →+ limTors A l)
    (y : limTors A l) : e (T y) = Matrix.mulVec (matOf e T) (e y) := by
  have h := Classical.choose_spec (exists_matrix_of_addHom
    ((e.toAddMonoidHom.comp T).comp e.symm.toAddMonoidHom)) (e y)
  simpa only [matOf, AddMonoidHom.coe_comp, Function.comp_apply, AddEquiv.coe_toAddMonoidHom,
    AddEquiv.symm_apply_apply] using h

theorem matOf_comp {A : Type u} [AddCommGroup A] {l : ℕ} [Fact l.Prime]
    (e : limTors A l ≃+ (Fin 2 → ℤ_[l])) (T S : limTors A l →+ limTors A l) :
    matOf e (T.comp S) = matOf e T * matOf e S := by
  refine matrix_ext_of_mulVec (fun x => ?_)
  have h1 := matOf_apply e (T.comp S) (e.symm x)
  have h2 := matOf_apply e T (S (e.symm x))
  have h3 := matOf_apply e S (e.symm x)
  simp only [AddMonoidHom.coe_comp, Function.comp_apply, AddEquiv.apply_symm_apply] at h1 h2 h3
  rw [← Matrix.mulVec_mulVec, ← h3, ← h2, ← h1]

theorem matOf_id {A : Type u} [AddCommGroup A] {l : ℕ} [Fact l.Prime]
    (e : limTors A l ≃+ (Fin 2 → ℤ_[l])) :
    matOf e (AddMonoidHom.id (limTors A l)) = 1 := by
  refine matrix_ext_of_mulVec (fun x => ?_)
  have h1 := matOf_apply e (AddMonoidHom.id (limTors A l)) (e.symm x)
  simp only [AddMonoidHom.id_apply, AddEquiv.apply_symm_apply] at h1
  rw [← h1, Matrix.one_mulVec]

/-- ★`ℤ_l × ℤ_l ≃+ (Fin 2 → ℤ_l)`。 -/
def prodEquivPiTwo (l : ℕ) [Fact l.Prime] : (ℤ_[l] × ℤ_[l]) ≃+ (Fin 2 → ℤ_[l]) where
  toFun := fun p => ![p.1, p.2]
  invFun := fun f => (f 0, f 1)
  left_inv := by intro p; simp
  right_inv := by
    intro f
    funext i
    fin_cases i <;> simp
  map_add' := by
    intro p q
    funext i
    fin_cases i <;> simp

/-! ## ★出典の紐付け(`.src`) -/

def addMonoidHom_padic_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現——加法写像の Z_l 線型性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
