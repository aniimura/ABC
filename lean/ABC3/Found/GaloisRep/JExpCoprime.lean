/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim

/-!
# 第 1208 ブロック —— **`l` と局所高さの素性は分岐指数が `l` と素なら底変換で保たれる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1199 の最後の仮説

第 1199（`lemma_3_5_height_ineq_over_extension`）の仮説のうち

    ∀ P, jExp P (E ⊗ L'') < 0 → ¬ (l ∣ jExp P (E ⊗ L''))

だけが、`L` の側の `PrimeToLocalHeights` から**自動では出ない**。

☆`jExp_baseChange`（在庫）が `v_P(j) = e(P|p) · v_p(j)` を与えるので、
`l` が素数なら `l ∤ e(P|p)` と `l ∤ v_p(j)` から出る。

## ★★★★測ったこと——**新しい義務が 1 つ見えた**

したがって残るのは **`l ∤ e(P|p)`** である。原文の道では
`L'' = L(H)` の次数が `l−1` を割るのでこれは自動だが、
本形式化ではまだその段（`Gal(L(H)/L) ↪ 𝔽_l^×`）を取っていない。

☆`M` は `{k • Q}` の座標で生成するので `Gal`-安定（第 1205 の直線の安定性）であり、
標数 0 だから `M/L` は Galois、よって `e ∣ [M:L]` である。
★そこから `[M:L] ∣ l−1` を出せば `l ∤ e` が従う——**これが次の義務**である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

variable (L L' : Type) [Field L] [NumberField L] [Field L'] [NumberField L']
  [Algebra L L'] [IsScalarTower ℚ L L']

/-- ★★★★★★★★★★★★
**分岐指数が `l` と素なら `l ∤ v_P(j)` は底変換で保たれる**——★**無条件**（第 1208）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`jExp_baseChange`（在庫）が `v_P(j) = e(P|p) · v_p(j)` を与えるので、
`l` が素数なら積の一方を割ることになり、両方の仮定に反する。

★★★これが第 1199 の最後の仮説
（`∀ P, jExp P < 0 → ¬ (l ∣ jExp P)`）を `L` の側へ降ろす段である。 -/
theorem not_dvd_jExp_baseChange (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) [E.IsElliptic]
    {l : ℕ} (hl : l.Prime)
    (he : ¬ ((l : ℤ) ∣ (p.asIdeal.ramificationIdx P.asIdeal : ℤ)))
    (hj : ¬ ((l : ℤ) ∣ jExp p E)) :
    ¬ ((l : ℤ) ∣ jExp P (E.baseChange L')) := by
  rw [jExp_baseChange L L' p P E]
  intro hdvd
  have hp : Prime (l : ℤ) := Nat.prime_iff_prime_int.mp hl
  rcases hp.dvd_mul.1 hdvd with h | h
  · exact he h
  · exact hj h

/-! ## ★出典の紐付け(`.src`) -/

def not_dvd_jExp_baseChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(分岐指数が l と素なら l ∤ v_P(j) は底変換で保たれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def not_dvd_jExp_baseChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "jExp_baseChange(証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.jExp_baseChange") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1208）**——第 1199 の最後の仮説を `L` の側へ降ろす段である。" ++
       "☆残るのは **`l ∤ e(P|p)`**。原文の道では `L'' = L(H)` の次数が `l−1` を" ++
       "割るので自動だが、本形式化ではまだ `Gal(L(H)/L) ↪ 𝔽_l^×` を取っていない。" ++
       "★`M` は `{k • Q}` の座標で生成するので `Gal`-安定（第 1205）であり、" ++
       "標数 0 だから `M/L` は Galois、よって `e ∣ [M:L]` である——" ++
       "そこから `[M:L] ∣ l−1` を出せばよい。") 3 ]

end ABC3.Found.GaloisRep
