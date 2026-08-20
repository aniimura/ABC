import ABC3.Found.GaloisRep.CountFP

/-!
# Galois (G5) 第 186 ブロック —— **★★★★★★★捩れ点は `n` で割れる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★交代性が要求する唯一のもの

交代性 `e_n(P,P) = 1` の古典的な証明は、`P` の **`n` 等分点** `P'`(`n·P' = P`)を取る。
★これは `[n] : E(F̄) → E(F̄)` の全射性であり、普通は**分割多項式**
`Φ_n(x) − a·Ψ_n(x)² = 0` の根の存在(次数 `n²`、`F` は代数閉)で示す。
★★しかし我々の `Skeleton/GaloisRep/WeilFunctionField.lean` の `x([n]P) = Φ_n/ΨSq_n` は
**臨界路の外**として残してある。

### ★★★★★★★位数だけで足りる

我々には第 (TorsionCard) の

    #E[n] = n²      (F 代数閉、標数が n! を割らない)

がある。★★これだけで `[n] : E[n²] → E[n]` の全射が**数え上げで**出る:

| 段 | 内容 |
|---|---|
| 1 | `φ : E[n²] →+ E`、`Q ↦ n·Q`。像は `E[n]` に入る |
| 2 | `ker φ = E[n]`(`E[n] ≤ E[n²]` だから) |
| 3 | Lagrange: `#E[n²] = #(E[n²]/ker φ) · #ker φ`、第一同型で `#image · n²` |
| 4 | `n⁴ = #image · n²` すなわち `#image = n² = #E[n]` |
| 5 | `image ≤ E[n]` で位数が等しい ⟹ `image = E[n]` |

★★★**分割多項式は要らない**。★これは「臨界路の外」と印を付けた項目を
**迂回できた**という記録である。

### ★逸脱(記録)

`hchar` の範囲が `k ≤ n` から **`k ≤ n²`** に広がる(`#E[n²] = n⁴` を使うため)。
★最終消費者 `det_cyclotomic` は `[CharZero K]` の下で述べるので無害である。
★★第 175/179 で入れた `hchar : ∀ k, 1 ≤ k → k ≤ n → (k:F) ≠ 0` と同じ理由の逸脱。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `nsmulEndo` | `n • ·` を自己準同型として見る |
| `exists_nsmul_eq_of_card` | ★★★★★★★**位数から割り算**(抽象群) |
| `exists_nsmul_eq_point` | ★★★★★★★**捩れ点は `n` で割れる**(楕円曲線) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

universe u

/-! ## ★抽象群の段 -/

/-- `n • ·` を加法群の自己準同型として見る。 -/
def nsmulEndo {A : Type u} [AddCommGroup A] (n : ℕ) : A →+ A :=
  AddMonoidHom.mk' (fun P : A => n • P) (fun a b => smul_add n a b)

@[simp] theorem nsmulEndo_apply {A : Type u} [AddCommGroup A] (n : ℕ) (P : A) :
    nsmulEndo n P = n • P := rfl

theorem mem_ker_nsmulEndo {A : Type u} [AddCommGroup A] (n : ℕ) (P : A) :
    P ∈ (nsmulEndo (A := A) n).ker ↔ n • P = 0 := AddMonoidHom.mem_ker

theorem card_ker_nsmulEndo {A : Type u} [AddCommGroup A] (n : ℕ) :
    Nat.card ((nsmulEndo (A := A) n).ker) = Nat.card {P : A // n • P = 0} :=
  Nat.card_congr (Equiv.subtypeEquivRight (fun P => mem_ker_nsmulEndo n P))

/-- `A[n] ≤ A[n²]`。 -/
theorem ker_nsmulEndo_le_sq {A : Type u} [AddCommGroup A] (n : ℕ) :
    (nsmulEndo (A := A) n).ker ≤ (nsmulEndo (A := A) (n ^ 2)).ker := by
  intro P hP
  rw [mem_ker_nsmulEndo] at hP ⊢
  rw [pow_two, mul_smul, hP, smul_zero]

/-- ★★★★★★★**位数から割り算が出る**——`#A[n] = n²`・`#A[n²] = n⁴` なら
`A[n]` の元はすべて `n` で割れる。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★分割多項式を使わずに `[n]` の全射性(捩れへの制限)を得る道。 -/
theorem exists_nsmul_eq_of_card {A : Type u} [AddCommGroup A] (n : ℕ) (hn : 1 ≤ n)
    (hc1 : Nat.card {P : A // n • P = 0} = n ^ 2)
    (hc2 : Nat.card {P : A // (n ^ 2) • P = 0} = (n ^ 2) ^ 2)
    {P : A} (hP : n • P = 0) : ∃ Q : A, n • Q = P := by
  have hn2 : (0 : ℕ) < n ^ 2 := pow_pos (by omega) 2
  set K1 := (nsmulEndo (A := A) n).ker with hK1
  set K2 := (nsmulEndo (A := A) (n ^ 2)).ker with hK2
  set φ : K2 →+ A := (nsmulEndo (A := A) n).comp K2.subtype with hφ
  have hkerφ : φ.ker = K1.addSubgroupOf K2 := rfl
  have hck1 : Nat.card K1 = n ^ 2 := (card_ker_nsmulEndo n).trans hc1
  have hcker : Nat.card φ.ker = n ^ 2 := by
    rw [hkerφ, Nat.card_congr (AddSubgroup.addSubgroupOfEquivOfLe (ker_nsmulEndo_le_sq n)).toEquiv]
    exact hck1
  have hck2 : Nat.card K2 = (n ^ 2) ^ 2 := (card_ker_nsmulEndo (n ^ 2)).trans hc2
  have hlag := AddSubgroup.card_eq_card_quotient_mul_card_addSubgroup (φ.ker)
  rw [Nat.card_congr (QuotientAddGroup.quotientKerEquivRange φ).toEquiv, hcker, hck2] at hlag
  have hcrange : Nat.card φ.range = n ^ 2 := by
    have hsq : (n ^ 2) ^ 2 = n ^ 2 * n ^ 2 := by ring
    rw [hsq] at hlag
    exact Nat.eq_of_mul_eq_mul_right hn2 hlag.symm
  have hrle : φ.range ≤ K1 := by
    rintro y ⟨⟨x, hx⟩, rfl⟩
    rw [hK1, mem_ker_nsmulEndo]
    rw [hK2, mem_ker_nsmulEndo] at hx
    show n • (n • x) = 0
    rw [← mul_smul, ← pow_two, hx]
  haveI : Finite K1 := Nat.finite_of_card_ne_zero (by rw [hck1]; omega)
  have htop : φ.range.addSubgroupOf K1 = ⊤ := by
    refine AddSubgroup.eq_top_of_card_eq _ ?_
    rw [Nat.card_congr (AddSubgroup.addSubgroupOfEquivOfLe hrle).toEquiv, hcrange, hck1]
  have hle2 : K1 ≤ φ.range := AddSubgroup.addSubgroupOf_eq_top.mp htop
  have hPK1 : P ∈ K1 := by rw [hK1, mem_ker_nsmulEndo]; exact hP
  obtain ⟨⟨Q, _⟩, hQP⟩ := hle2 hPK1
  exact ⟨Q, hQP⟩

/-! ## ★楕円曲線の段 -/

/-- ★★★★★★★**捩れ点は `n` で割れる**——`n·P = 0` なら `n·Q = P` なる `Q` がある。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`#E[n] = n²`(`torsion_card`)を `n` と `n²` の 2 回使うだけ。
★★交代性 `e_n(P,P) = 1` が要求する唯一の新しい入力である。 -/
theorem exists_nsmul_eq_point {F : Type u} [Field F] [DecidableEq F] [IsAlgClosed F]
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n ^ 2 → (k : F) ≠ 0)
    {P : W.Point} (hP : n • P = 0) : ∃ Q : W.Point, n • Q = P := by
  have hΔ : W.Δ ≠ 0 := W.isElliptic_iff.mp inferInstance |>.ne_zero
  refine exists_nsmul_eq_of_card n hn ?_ ?_ hP
  · exact torsion_card W hΔ n hn
      (fun k h1 h2 => hchar k h1 (le_trans h2 (Nat.le_self_pow two_ne_zero n)))
  · exact torsion_card W hΔ (n ^ 2) (Nat.one_le_pow 2 n (by omega)) hchar

/-! ## ★出典の紐付け(`.src`) -/

def exists_nsmul_eq_of_card.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——位数からの割り算)",
    sectionId := "genell-thm-3-8" }

def exists_nsmul_eq_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——捩れ点の n 等分)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
