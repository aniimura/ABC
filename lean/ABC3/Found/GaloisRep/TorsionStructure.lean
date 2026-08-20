import ABC3.Found.GaloisRep.ZmodSqEquiv
import ABC3.Interface.GaloisRep.Torsion

/-!
# Galois (G1) 第 72 ブロック —— **★★★★★★★★G1 達成 `E[n] ≃+ (ℤ/n)²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★鎖が繋がった

| 段 | ブロック |
|---|---|
| Somos-4 と (Λ) | 47, 56 |
| **Ward の定理** | 58 |
| 乗法公式 | 52 |
| 零集合 = 位数の倍数 | 60 |
| `Φ_n` と `ΨSq_n` の互いに素性 | 61, 62 |
| Wronskian ≠ 0 → 分離性 | 55, 63 |
| **`#E[n] = n²`** | 65, 66 |
| 有限アーベル群の構造定理 | 67–71 |
| **`E[n] ≃+ (ℤ/n)²`** | ★本ブロック |

## ★★★Interface を強めた(§9-400)

元の `structure_eq` は `≃`(型の全単射)で、`Nat.card` だけで埋まってしまった。
★★原典 Theorem 3.8 は**群同型**を主張するので **`≃+`** に強め、
標数の仮定 `∀ k ≤ n, (k : K) ≠ 0` を足した(§9-398 の記録)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `torsion_addEquiv` | ★★★★★★★★**`(ℤ/n)² ≃+ E[n]`** |
| `torsionStructureDataWitness` | ★`TorsionStructureData` の実装 |
| `TorsionStructureData.nonvacuous` | ★★★★★★★★**G1 達成** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine ABC3.Interface.GaloisRep

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★★★★**`(ℤ/n)² ≃+ E[n]`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem torsion_addEquiv [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    Nonempty ((ZMod n × ZMod n) ≃+ (nsmulHom W.toAffine.Point n).ker) := by
  haveI : Finite ((nsmulHom W.toAffine.Point n).ker) :=
    (finite_torsion W n hn (by exact_mod_cast hchar n hn le_rfl)).to_subtype
  refine addEquiv_zmod_sq n hn (fun x => by ext; exact x.2) ?_
  intro d hd
  have hd1 : 1 ≤ d := Nat.pos_of_dvd_of_pos hd (by omega)
  rw [Nat.card_congr (kerSubEquiv n d hd).toEquiv]
  exact torsion_card_dvd W hΔ n hn hchar d hd1 hd

/-- ★`Interface` の `torsionPoints` は `n` 倍写像の核そのものである。 -/
theorem torsionPoints_eq_ker {K : Type} [Field K] [DecidableEq K] (W : WeierstrassCurve K)
    (n : ℕ) : torsionPoints W n = (nsmulHom W.toAffine.Point n).ker := rfl

/-- ★★★★★★★★**`TorsionStructureData` の実装**。 -/
noncomputable def torsionStructureDataWitness : TorsionStructureData where
  structure_eq := by
    intro K _ _ _ W hell n hn hchar
    haveI := hell
    have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
    obtain ⟨e⟩ := torsion_addEquiv W hΔ n hn hchar
    rw [torsionPoints_eq_ker]
    exact ⟨e.symm⟩
  torsion_finite := by
    intro K _ _ _ W hell n hn hchar
    haveI := hell
    exact (finite_torsion W n hn hchar).to_subtype

/-- ★★★★★★★★**`TorsionStructureData` は非空虚である**——G1 達成。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが Galois 表現論の 8 件のうち **G1** である。 -/
theorem TorsionStructureData.nonvacuous : Nonempty TorsionStructureData :=
  ⟨torsionStructureDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def torsion_addEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(n-捩れ部分群の構造定理)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
