import ABC3.Found.GaloisRep.DetCyclotomic

/-!
# Galois (G5) 第 208 ブロック —— **★★★★★★塔の各段の全射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★基底性に要るもの

第 207 の `hgen` のうち**基底性**は、`tateProj n : T_l E → E[l^n]` が**全射**であることに帰着する
(`e.symm` は全単射なので `tateVec n` の全射性と同値、そこから
`R = a·P + b·Q` が第 206 で出る)。

★逆極限の全射性は、塔の各段 **`[l] : E[l^{m+1}] → E[l^m]` が全射**であることから出る。

### ★★★★★★第 186 の一般化で出る

第 186 は `#A[n] = n²`・`#A[n²] = n⁴` から `A[n]` が `n` で割れることを出した。
★同じ数え上げを `[n] : A[nm] → A[m]` に当てると:

| 段 | 内容 |
|---|---|
| 1 | `φ : A[nm] →+ A`、`Q ↦ n·Q`。像は `A[m]` に入る |
| 2 | `ker φ = A[n]`(`A[n] ≤ A[nm]` だから) |
| 3 | Lagrange と第一同型で `(nm)² = #image · n²`、すなわち `#image = m²` |
| 4 | `image ≤ A[m]` で位数が等しい ⟹ `image = A[m]` |

★★`n := l`・`m := l^k` と取れば **`[l] : E[l^{k+1}] → E[l^k]` の全射**になる。

## ★★残るのは逆極限の構成だけ

`R ∈ E[l^n]` から両立列 `f`(`l·f(m+1) = f(m)`)を作る段。
★`m ≤ n` では `f m := l^{n−m}·R`、`m > n` では本ブロックの全射性で**再帰的に選ぶ**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_nsmul_eq_of_card_gen` | ★★★★★★★**位数から `[n] : A[nm] → A[m]` の全射** |
| `exists_smul_step` | ★★★★★★**`[l] : E[l^{m+1}] → E[l^m]` は全射** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

universe u

/-- ★★★★★★★**位数から `[n] : A[nm] → A[m]` の全射**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 186 の `exists_nsmul_eq_of_card` の一般化(そこは `m = n` の場合)。 -/
theorem exists_nsmul_eq_of_card_gen {A : Type u} [AddCommGroup A] (n m : ℕ)
    (hn : 1 ≤ n) (hm : 1 ≤ m)
    (hc1 : Nat.card {P : A // m • P = 0} = m ^ 2)
    (hc2 : Nat.card {P : A // (n * m) • P = 0} = (n * m) ^ 2)
    (hck : Nat.card {P : A // n • P = 0} = n ^ 2)
    {P : A} (hP : m • P = 0) : ∃ Q : A, (n * m) • Q = 0 ∧ n • Q = P := by
  have hn2 : (0 : ℕ) < n ^ 2 := pow_pos (by omega) 2
  set K1 := (nsmulEndo (A := A) m).ker with hK1
  set K2 := (nsmulEndo (A := A) (n * m)).ker with hK2
  set Kn := (nsmulEndo (A := A) n).ker with hKn
  have hle : Kn ≤ K2 := by
    intro x hx
    rw [hKn, mem_ker_nsmulEndo] at hx
    rw [hK2, mem_ker_nsmulEndo, mul_comm, mul_smul, hx, smul_zero]
  set φ : K2 →+ A := (nsmulEndo (A := A) n).comp K2.subtype with hφ
  have hkerφ : φ.ker = Kn.addSubgroupOf K2 := rfl
  have hckn : Nat.card Kn = n ^ 2 := (card_ker_nsmulEndo n).trans hck
  have hcker : Nat.card φ.ker = n ^ 2 := by
    rw [hkerφ, Nat.card_congr (AddSubgroup.addSubgroupOfEquivOfLe hle).toEquiv]
    exact hckn
  have hck2 : Nat.card K2 = (n * m) ^ 2 := (card_ker_nsmulEndo (n * m)).trans hc2
  have hck1 : Nat.card K1 = m ^ 2 := (card_ker_nsmulEndo m).trans hc1
  have hlag := AddSubgroup.card_eq_card_quotient_mul_card_addSubgroup (φ.ker)
  rw [Nat.card_congr (QuotientAddGroup.quotientKerEquivRange φ).toEquiv, hcker, hck2] at hlag
  have hcrange : Nat.card φ.range = m ^ 2 := by
    have hsq : (n * m) ^ 2 = m ^ 2 * n ^ 2 := by ring
    rw [hsq] at hlag
    exact Nat.eq_of_mul_eq_mul_right hn2 hlag.symm
  have hrle : φ.range ≤ K1 := by
    rintro y ⟨⟨x, hx⟩, rfl⟩
    rw [hK1, mem_ker_nsmulEndo]
    rw [hK2, mem_ker_nsmulEndo] at hx
    show m • (n • x) = 0
    rw [← mul_smul, mul_comm, hx]
  haveI : Finite K1 := Nat.finite_of_card_ne_zero (by rw [hck1]; positivity)
  have htop : φ.range.addSubgroupOf K1 = ⊤ := by
    refine AddSubgroup.eq_top_of_card_eq _ ?_
    rw [Nat.card_congr (AddSubgroup.addSubgroupOfEquivOfLe hrle).toEquiv, hcrange, hck1]
  have hle2 : K1 ≤ φ.range := AddSubgroup.addSubgroupOf_eq_top.mp htop
  have hPK1 : P ∈ K1 := by rw [hK1, mem_ker_nsmulEndo]; exact hP
  obtain ⟨⟨Q, hQ⟩, hQP⟩ := hle2 hPK1
  refine ⟨Q, ?_, hQP⟩
  rw [hK2, mem_ker_nsmulEndo] at hQ
  exact hQ

/-- ★★★★★★**`[l] : E[l^{m+1}] → E[l^m]` は全射**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★逆極限 `T_l E → E[l^n]` の全射性の各段である。 -/
theorem exists_smul_step {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F]
    (W : WeierstrassCurve F) (hΔ : W.Δ ≠ 0) (l m : ℕ) (hl : 1 ≤ l)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ l ^ (m + 1) → (k : F) ≠ 0)
    {T : W.toAffine.Point} (hT : (l ^ m) • T = 0) :
    ∃ S : W.toAffine.Point, (l ^ (m + 1)) • S = 0 ∧ l • S = T := by
  have hlm : 1 ≤ l ^ m := Nat.one_le_pow _ _ (by omega)
  have hlm1 : 1 ≤ l ^ (m + 1) := Nat.one_le_pow _ _ (by omega)
  have hmul : l * l ^ m = l ^ (m + 1) := by rw [← pow_succ']
  have hcharle : ∀ j : ℕ, j ≤ l ^ (m + 1) → ∀ k : ℕ, 1 ≤ k → k ≤ j → (k : F) ≠ 0 :=
    fun j hj k hk1 hk2 => hchar k hk1 (le_trans hk2 hj)
  have hc1 : Nat.card {P : W.toAffine.Point // (l ^ m) • P = 0} = (l ^ m) ^ 2 :=
    torsion_card W hΔ (l ^ m) hlm (hcharle _ (Nat.pow_le_pow_right (by omega) (by omega)))
  have hc2 : Nat.card {P : W.toAffine.Point // (l * l ^ m) • P = 0} = (l * l ^ m) ^ 2 := by
    rw [hmul]
    exact torsion_card W hΔ (l ^ (m + 1)) hlm1 (hcharle _ le_rfl)
  have hck : Nat.card {P : W.toAffine.Point // l • P = 0} = l ^ 2 :=
    torsion_card W hΔ l hl (hcharle _ (by
      calc l = l ^ 1 := (pow_one l).symm
      _ ≤ l ^ (m + 1) := Nat.pow_le_pow_right (by omega) (by omega)))
  obtain ⟨Q, hQ1, hQ2⟩ := exists_nsmul_eq_of_card_gen l (l ^ m) hl hlm hc1 hc2 hck hT
  exact ⟨Q, by rw [← hmul]; exact hQ1, hQ2⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_nsmul_eq_of_card_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——位数から [n] : A[nm] → A[m] の全射)",
    sectionId := "genell-thm-3-8" }

def exists_smul_step.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——塔の各段の全射性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
