import ABC3.Found.GaloisRep.TorsionStructure

/-!
# Galois (G2) 第 73 ブロック —— **★★★捩れの塔は全射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★Tate 加群へ向かう最初の一歩

`T_l E = lim_n E[l^n]` を作るには、塔

    E[l] ← E[l²] ← E[l³] ← ⋯   (射は `l` 倍)

が**全射**であることが要る。★これも数え上げで出る:

    #range(l·) · #ker(l·) = #A[l^{n+1}]
      ⟹ #range = l^{2(n+1)} / l² = l^{2n} = #A[l^n]
    range(l·) ⊆ A[l^n]  ⟹  **range(l·) = A[l^n]**

★★第 68 ブロック(`range_eq_ker`)の**鏡像**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `range_nsmul_eq_ker` | ★★`l` 倍の像がちょうど `A[l^n]` |
| `exists_smul_lift` | ★★★**塔の全射性**(持ち上げの存在) |
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★★指数 `l^(n+1)`、位数 `l^(2n+2)` なら、`l` 倍の像はちょうど `A[l^n]`。 -/
theorem range_nsmul_eq_ker {A : Type u} [AddCommGroup A] [Finite A] (l n : ℕ)
    (hcardA : Nat.card A = (l ^ (n + 1)) ^ 2)
    (hker1 : Nat.card (nsmulHom A l).ker = l ^ 2)
    (hkern : Nat.card (nsmulHom A (l ^ n)).ker = (l ^ n) ^ 2)
    (hl : 1 < l)
    (hexp : ∀ x : A, (l ^ (n + 1)) • x = 0) :
    (nsmulHom A l).range = (nsmulHom A (l ^ n)).ker := by
  have hle : (nsmulHom A l).range ≤ (nsmulHom A (l ^ n)).ker := by
    rintro _ ⟨y, rfl⟩
    rw [AddMonoidHom.mem_ker, nsmulHom_apply, nsmulHom_apply, smul_smul,
      show l ^ n * l = l ^ (n + 1) by ring]
    exact hexp y
  have hrng := card_range_mul_card_ker (A := A) l
  rw [hker1, hcardA] at hrng
  have hcr : Nat.card (nsmulHom A l).range = (l ^ n) ^ 2 := by
    have hsplit : ((l ^ (n + 1)) ^ 2 : ℕ) = (l ^ n) ^ 2 * l ^ 2 := by ring
    rw [hsplit] at hrng
    exact Nat.eq_of_mul_eq_mul_right (by positivity) hrng
  exact AddSubgroup.eq_of_le_of_card_ge hle (by rw [hcr, hkern])

/-- ★★★**捩れの塔は全射である**——`A[l^n]` の点は `A[l^{n+1}]` から `l` 倍で来る。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem exists_smul_lift {A : Type u} [AddCommGroup A]
    (hfin : ∀ m : ℕ, 1 ≤ m → Finite (nsmulHom A m).ker)
    (hcard : ∀ m : ℕ, 1 ≤ m → Nat.card (nsmulHom A m).ker = m ^ 2)
    (l n : ℕ) (hl : 1 < l) (x : A) (hx : (l ^ n) • x = 0) :
    ∃ y : A, l • y = x ∧ (l ^ (n + 1)) • y = 0 := by
  haveI := hfin (l ^ (n + 1)) (Nat.one_le_pow _ _ (by omega))
  have hxB : x ∈ (nsmulHom A (l ^ (n + 1))).ker := by
    show (l ^ (n + 1)) • x = 0
    rw [pow_succ, mul_comm, mul_smul, hx, smul_zero]
  have hcardA : Nat.card (nsmulHom A (l ^ (n + 1))).ker = (l ^ (n + 1)) ^ 2 :=
    hcard _ (Nat.one_le_pow _ _ (by omega))
  have hker1 : Nat.card (nsmulHom ((nsmulHom A (l ^ (n + 1))).ker) l).ker = l ^ 2 := by
    rw [Nat.card_congr (kerSubEquiv (l ^ (n + 1)) l (dvd_pow_self l (by omega))).toEquiv]
    exact hcard l (by omega)
  have hkern : Nat.card (nsmulHom ((nsmulHom A (l ^ (n + 1))).ker) (l ^ n)).ker = (l ^ n) ^ 2 := by
    rw [Nat.card_congr (kerSubEquiv (l ^ (n + 1)) (l ^ n) (pow_dvd_pow l (by omega))).toEquiv]
    exact hcard _ (Nat.one_le_pow _ _ (by omega))
  have hexp : ∀ z : (nsmulHom A (l ^ (n + 1))).ker, (l ^ (n + 1)) • z = 0 := fun z => by
    ext; exact z.2
  have hrange := range_nsmul_eq_ker (A := (nsmulHom A (l ^ (n + 1))).ker) l n hcardA hker1 hkern
    hl hexp
  have hmem : (⟨x, hxB⟩ : (nsmulHom A (l ^ (n + 1))).ker) ∈ (nsmulHom _ (l ^ n)).ker := by
    show (l ^ n) • (⟨x, hxB⟩ : (nsmulHom A (l ^ (n + 1))).ker) = 0
    ext; exact hx
  rw [← hrange] at hmem
  obtain ⟨y, hy⟩ := hmem
  refine ⟨(y : A), ?_, y.2⟩
  have hy' := congrArg (Subtype.val) hy
  simpa [nsmulHom_apply] using hy'

/-! ## ★出典の紐付け(`.src`) -/

def exists_smul_lift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 加群——捩れの塔の全射性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
