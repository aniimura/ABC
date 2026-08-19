import ABC3.Found.GaloisRep.CoprimeTorsion

/-!
# Galois (G1) 第 38 ブロック —— **★★★★★★`E[mn]` は互いに素でなくても出る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★第 37 を一般化した

第 37 は `gcd(m,n) = 1` を要求していたが、★**要らない**:

    E[m] 有限,  E[n] 有限  ⟹  E[mn] 有限

★★機構: `P ↦ m•P` は `E[mn] → E[n]` で、**核が `E[m]`**。
核と像が有限なら、もとの集合も有限(剰余類が有限個の有限集合)。

## ★★★これで素冪が出る

    E[m] 有限  ⟹  E[m^k] 有限   (k に関する帰納)

★★★★したがって **`E[p]`(素数)さえ示せば `E[n]` はすべて出る**
——`n` を素因数分解して第 38 を繰り返すだけ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `finite_of_ker_image` | ★★★核と像が有限なら有限(一般の可換群) |
| `finite_torsion_mul` | ★★★★**`E[mn]`(互いに素不要)** |
| `finite_torsion_pow` | ★★★**`E[m^k]`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

universe u v

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- 核と像が有限なら、もとの集合も有限。 -/
theorem finite_of_ker_image {A : Type u} {B : Type v} [AddCommGroup A] [AddCommGroup B]
    (f : A →+ B) {s : Set A} (hker : {a : A | f a = 0}.Finite) (himg : (f '' s).Finite) :
    s.Finite := by
  classical
  have hsub : s ⊆ ⋃ b ∈ f '' s, {a | a ∈ s ∧ f a = b} := by
    intro a ha
    exact Set.mem_biUnion ⟨a, ha, rfl⟩ ⟨ha, rfl⟩
  refine Set.Finite.subset (Set.Finite.biUnion himg ?_) hsub
  intro b _
  rcases Set.eq_empty_or_nonempty {a | a ∈ s ∧ f a = b} with h | ⟨a₀, ha₀⟩
  · rw [h]; exact Set.finite_empty
  · refine Set.Finite.of_finite_image (f := fun a => a - a₀) ?_
      (Set.injOn_of_injective (fun u v huv => by simpa using huv))
    refine hker.subset ?_
    rintro _ ⟨a, ha, rfl⟩
    show f (a - a₀) = 0
    rw [map_sub, ha.2, ha₀.2, sub_self]

/-- ★互いに素でなくてよい版。 -/
theorem finite_torsion_mul (m n : ℕ)
    (hm : {P : W.toAffine.Point | m • P = 0}.Finite)
    (hn : {P : W.toAffine.Point | n • P = 0}.Finite) :
    {P : W.toAffine.Point | (m * n) • P = 0}.Finite := by
  refine finite_of_ker_image
    (AddMonoidHom.mk' (fun P : W.toAffine.Point => m • P) (fun a b => smul_add m a b)) hm ?_
  refine hn.subset ?_
  rintro _ ⟨P, hP, rfl⟩
  show n • (m • P) = 0
  rw [smul_smul, Nat.mul_comm]
  exact hP


/-- ★★★**`E[m^k]` は `E[m]` から出る**。 -/
theorem finite_torsion_pow (m : ℕ) (hm : {P : W.toAffine.Point | m • P = 0}.Finite) :
    ∀ k : ℕ, {P : W.toAffine.Point | (m ^ k) • P = 0}.Finite
  | 0 => by
    simp only [pow_zero, one_smul]
    exact Set.Finite.subset (Set.finite_singleton (0 : W.toAffine.Point))
      (fun P hP => by simpa using hP)
  | (k + 1) => by
    rw [pow_succ]
    exact finite_torsion_mul W (m ^ k) m (finite_torsion_pow m hm k) hm

/-! ## ★出典の紐付け(`.src`) -/

def finite_torsion_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(E[mn] は互いに素でなくても出ること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
