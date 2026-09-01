/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TorsionCard
import ABC3.Found.GaloisRep.TorsionAll
import ABC3.Interface.GaloisRep.Representation
import ABC3.Meta.Claim

/-!
# 第 1271 ブロック —— **代数閉体の埋め込みで `E[n]` は移り、条件も移る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——II 側の配管の第 2 段

第 1270 で `T_l E` の 2 条件は `E[l]` の条件に落ちた。
★次は**局所で作った `σ` の条件を大域の `E[l]` へ運ぶ**段である。

☆道具は 2 つだけ:

| # | 内容 | 由来 |
|---|---|---|
| 1 | `Point.map` は `galPoint` と可換（埋め込みが `σ` と可換なら） | `Point.map_map` |
| 2 | 代数閉体の間では `E[n]` の写像は**全単射** | `torsion_card`（両側で `n²`）＋ `Point.map_injective` |

★★★2 は**次数の勘定だけ**である——単射で、しかも両側の個数が `n²` で等しい。
☆したがって分体多項式も Néron モデルも要らない。

★★これで「局所の `σ` が `E[l]` に `α` として作用する」を確かめれば、
`restrictLocalHom`（第 1167）で大域へ運べる。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

variable {K₀ F K : Type} [Field K₀] [DecidableEq K₀] [Field F] [DecidableEq F]
  [Field K] [DecidableEq K] [Algebra K₀ F] [Algebra K₀ K]

/-- ★★★★★★★★★★★★
**埋め込みは `galPoint` と可換**——★**無条件**（第 1271）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`galPoint W σ` は `Point.map σ.toAlgHom` なので、
`Point.map_map` で合成に直せば `f ∘ σF = σK ∘ f` そのものになる。 -/
theorem point_map_galPoint (W : WeierstrassCurve K₀) (f : F →ₐ[K₀] K)
    (σF : F ≃ₐ[K₀] F) (σK : K ≃ₐ[K₀] K) (hcomm : ∀ x : F, f (σF x) = σK (f x))
    (P : (W.baseChange F).toAffine.Point) :
    Point.map f (galPoint W σF P) = galPoint W σK (Point.map f P) := by
  show Point.map f (Point.map σF.toAlgHom P) = Point.map σK.toAlgHom (Point.map f P)
  rw [Point.map_map, Point.map_map]
  have hfe : f.comp σF.toAlgHom = σK.toAlgHom.comp f := AlgHom.ext fun x => hcomm x
  rw [hfe]

/-- ★★★★★★★★**`n`-捩れの上の写像**（第 1271）。 -/
noncomputable def torsionMap (W : WeierstrassCurve K₀) (f : F →ₐ[K₀] K) (n : ℕ) :
    {P : (W.baseChange F).toAffine.Point // n • P = 0} →
      {P : (W.baseChange K).toAffine.Point // n • P = 0} :=
  fun P => ⟨Point.map f (P : (W.baseChange F).toAffine.Point), by
    rw [← map_nsmul, P.2, map_zero]⟩

/-- ★★★★★★★★★★★★★★★★
**代数閉体の埋め込みは `E[n]` の間の全単射を与える**——★**無条件**（第 1271）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆単射性は `Point.map_injective`、個数はどちらも `n²`（`torsion_card`）。
★★**有限集合の間の単射で個数が等しいから全射**——それだけである。 -/
omit [DecidableEq K₀] in
theorem torsionMap_bijective [IsAlgClosed F] [IsAlgClosed K] (W : WeierstrassCurve K₀)
    (hΔF : (W.baseChange F).Δ ≠ 0) (hΔK : (W.baseChange K).Δ ≠ 0)
    (f : F →ₐ[K₀] K) (n : ℕ) (hn : 1 ≤ n)
    (hcF : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (hcK : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : K) ≠ 0) :
    Function.Bijective (torsionMap W f n) := by
  haveI hfF : Finite {P : (W.baseChange F).toAffine.Point // n • P = 0} :=
    (finite_torsion (W.baseChange F) n hn (hcF n hn le_rfl)).to_subtype
  haveI hfK : Finite {P : (W.baseChange K).toAffine.Point // n • P = 0} :=
    (finite_torsion (W.baseChange K) n hn (hcK n hn le_rfl)).to_subtype
  refine (Nat.bijective_iff_injective_and_card _).2 ⟨?_, ?_⟩
  · intro a b hab
    exact Subtype.ext (Point.map_injective (f := f) (Subtype.ext_iff.1 hab))
  · rw [torsion_card (W.baseChange F) hΔF n hn hcF, torsion_card (W.baseChange K) hΔK n hn hcK]

/-- ★★★★★★★★★★★★★★★★
**幂単性は埋め込みで降りる**——★**無条件**（第 1271）。

☆`Point.map f` は単射なので、像で成り立つ等式は元でも成り立つ。 -/
theorem galPoint_unipotent_of_map (W : WeierstrassCurve K₀) (f : F →ₐ[K₀] K)
    (σF : F ≃ₐ[K₀] F) (σK : K ≃ₐ[K₀] K) (hcomm : ∀ x : F, f (σF x) = σK (f x)) (n : ℕ)
    (h : ∀ Q : (W.baseChange K).toAffine.Point, n • Q = 0 →
      galPoint W σK (galPoint W σK Q) + Q = galPoint W σK Q + galPoint W σK Q) :
    ∀ P : (W.baseChange F).toAffine.Point, n • P = 0 →
      galPoint W σF (galPoint W σF P) + P = galPoint W σF P + galPoint W σF P := by
  intro P hP
  refine Point.map_injective (f := f) ?_
  have hQ : n • (Point.map f P) = 0 := by rw [← map_nsmul, hP, map_zero]
  have hkey := h (Point.map f P) hQ
  simp only [map_add, point_map_galPoint W f σF σK hcomm]
  exact hkey

/-- ★★★★★★★★★★★★★★★★
**非自明性も埋め込みで降りる**——★**無条件**（第 1271）。

☆こちらは**全射性**を使う——像の側の証人 `Q` を引き戻す。 -/
theorem exists_galPoint_ne_of_map [IsAlgClosed F] [IsAlgClosed K] (W : WeierstrassCurve K₀)
    (hΔF : (W.baseChange F).Δ ≠ 0) (hΔK : (W.baseChange K).Δ ≠ 0)
    (f : F →ₐ[K₀] K) (σF : F ≃ₐ[K₀] F) (σK : K ≃ₐ[K₀] K)
    (hcomm : ∀ x : F, f (σF x) = σK (f x)) (n : ℕ) (hn : 1 ≤ n)
    (hcF : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (hcK : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : K) ≠ 0)
    (h : ∃ Q : (W.baseChange K).toAffine.Point, n • Q = 0 ∧ galPoint W σK Q ≠ Q) :
    ∃ P : (W.baseChange F).toAffine.Point, n • P = 0 ∧ galPoint W σF P ≠ P := by
  obtain ⟨Q, hQ, hne⟩ := h
  obtain ⟨P, hP⟩ := (torsionMap_bijective W hΔF hΔK f n hn hcF hcK).2 ⟨Q, hQ⟩
  refine ⟨(P : (W.baseChange F).toAffine.Point), P.2, ?_⟩
  intro hcon
  refine hne ?_
  have hmap := congrArg (fun z => Point.map f z) hcon
  rw [point_map_galPoint W f σF σK hcomm] at hmap
  have hPQ : Point.map f ((P : (W.baseChange F).toAffine.Point)) = Q :=
    congrArg Subtype.val hP
  rwa [hPQ] at hmap

/-! ## ★出典の紐付け(`.src`) -/

def point_map_galPoint.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(埋め込みは galPoint と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def torsionMap.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(n-捩れの上の写像)",
    sectionId := "genell-thm-3-8" }

def torsionMap_bijective.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(代数閉体の埋め込みは E[n] の間の全単射。★無条件)",
    sectionId := "genell-thm-3-8" }

def galPoint_unipotent_of_map.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(幂単性は埋め込みで降りる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galPoint_ne_of_map.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(非自明性は埋め込みで降りる。★無条件)",
    sectionId := "genell-thm-3-8" }

def torsionMap_bijective.needs : List ProofObligation :=
  [ .citation "[ABC3]" "torsion_card(証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.torsion_card") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1271）**——分体多項式も Néron モデルも要らない。" ++
       "☆`Point.map` は単射（mathlib）、`E[n]` の個数はどちらの代数閉体でも `n²`" ++
       "（`torsion_card`、証明済み）なので、**次数の勘定だけで全単射**になる。" ++
       "★★これで局所で作った `σ` の条件を大域の `E[l]` へ運べる。") 2 ]

end ABC3.Found.GaloisRep
