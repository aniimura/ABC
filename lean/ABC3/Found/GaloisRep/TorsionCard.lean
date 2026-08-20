import ABC3.Found.GaloisRep.TorsionCount

/-!
# Galois (G1) 第 66 ブロック —— **★★★★★★★`Nat.card E[n] = n²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★`Interface` の形に合わせる

`Interface/GaloisRep/Torsion.lean` の `torsionPoints W n` は
`{P : W.toAffine.Point // n • P = 0}` である。★第 65 ブロックの `Set.ncard` の結果を
`Nat.card` に直し、`n = 1` の場合も込めた形にする。

## ★★約数についても同じ

`d ∣ n` なら `k ≤ d ≤ n` なので標数の仮定がそのまま効き、**`#E[d] = d²`** が出る。
★これが `E[n] ≅ (ℤ/n)²` を絞り込むのに要る。

## ★★★残る段(§9-399)

`A := E[n]` は

* 位数 `n²` の有限アーベル群
* 指数が `n` を割る
* `#A[d] = d²` (`d ∣ n`)

である。★★これから `A ≅ (ℤ/n)²` を出すのが `structure_eq` の最後の段である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `torsion_card` | ★★★★★★★**`Nat.card E[n] = n²`** |
| `torsion_card_dvd` | ★★約数版 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★★★**`Nat.card E[n] = n²`**(`n = 1` も含む)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem torsion_card [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    Nat.card {P : W.toAffine.Point // n • P = 0} = n ^ 2 := by
  rcases Nat.lt_or_ge n 2 with h1 | h2
  · have hn1 : n = 1 := by omega
    subst hn1
    haveI : Unique {P : W.toAffine.Point // (1 : ℕ) • P = 0} :=
      { default := ⟨0, by simp⟩
        uniq := by rintro ⟨P, hP⟩; simp only [one_smul] at hP; exact Subtype.ext hP }
    rw [Nat.card_unique]
    simp
  · show Nat.card ↑({P : W.toAffine.Point | n • P = 0} : Set _) = n ^ 2
    rw [Nat.card_coe_set_eq]
    exact torsion_ncard W hΔ n h2 hchar

/-- ★★**約数版**——`d ∣ n` なら `#E[d] = d²`。 -/
theorem torsion_card_dvd [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) (d : ℕ) (hd1 : 1 ≤ d) (hd : d ∣ n) :
    Nat.card {P : W.toAffine.Point // d • P = 0} = d ^ 2 :=
  torsion_card W hΔ d hd1 (fun k hk1 hk2 => hchar k hk1 (le_trans hk2 (Nat.le_of_dvd (by omega) hd)))

/-! ## ★出典の紐付け(`.src`) -/

def torsion_card.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(n-捩れ部分群の位数——Interface の形)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
